/**
 * @file Kernels.cu
 * @brief CUDA kernels and device functions for VCF file parsing
 * @author Your Name
 * @date 2025-07-16
 *
 * @details This file implements the GPU-accelerated parsing of VCF files through:
 *  - Line indexing: Identifies newline positions in parallel
 *  - Field parsing: Extracts and processes VCF fields (INFO, FORMAT, etc.)
 *  - Memory optimization: Uses shared and device memory efficiently
 *  - Batch processing: Handles large files in configurable chunks
 *
 * @note Uses __ldg for coalesced memory access to improve performance
 * @warning Requires compute capability 3.0 or higher for __ldg operations
 */

#ifndef KERNELS_CU
#define KERNELS_CU

#include "CUDAUtils.cuh"
#include "Utils.h"
#include "DataStructures.h"
#include "DataFrames.h"

#include <cuda_runtime.h>     
#include <cuda_fp16.h>  
#include <thrust/device_ptr.h> 
#include <thrust/sort.h>

#include <chrono>
#include <fstream>
#include <filesystem>
#include <sys/wait.h>
#include <unistd.h>
#include <map>
#include <omp.h> 

using namespace std;

/**
 * @brief CUDA kernel to find newline indices in an input buffer
 *
 * @param input [in] Input character buffer in device memory
 * @param len [in] Length of input buffer
 * @param output [out] Array to store newline indices
 * @param len_output [in] Size of output array
 * @param global_count [out] Atomic counter for output array
 *
 * @details Each thread:
 *  1. Processes one character using coalesced reads
 *  2. If newline found, atomically adds index to output
 *  3. Last thread writes sentinel value
 */
__global__ void cu_find_new_lines_index(const char* input, unsigned int len, unsigned int* output, 
                unsigned int len_output, unsigned int* global_count){
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;  // Indice globale del thread

    if (idx < len && __ldg(&input[idx]) == '\n') { //coaleasced read only
        unsigned int pos = atomicAdd(global_count, 1); //primo spazio libero dove salvare
        output[pos] = idx;  // Salva la posizione trovata nell'array finale == idx
    }else if(idx == len){
        //invece che mettere a 0 il primo valore dell'array che darebbe problemi per synch 
        //o per offset lo metto alla fine tanto poi va ordinato
        output[len_output-1] = 0; 
    }
}

/**
 * @brief Resets a temporary string buffer.
 *
 * Sets the temporary index to zero and initializes the buffer with a null terminator.
 *
 * @param tmp Pointer to the temporary character array.
 * @param tmp_idx Reference to the current index in the temporary buffer.
 */
__device__ void reset_tmp(char* tmp, int& tmp_idx) {
    tmp_idx = 0;
    tmp[0] = '\0';
};

/**
 * @brief Appends a character to a temporary string buffer.
 *
 * If there is space in the buffer (less than MAX_TOKEN_LEN-1 characters),
 * the character is appended and the buffer is null-terminated.
 *
 * @param tmp Pointer to the temporary character array.
 * @param tmp_idx Reference to the current index in the temporary buffer.
 * @param c The character to append.
 */
__device__ void append_tmp(char* tmp, int& tmp_idx, char c) {
    if (tmp_idx < MAX_TMP_LEN-1) tmp[tmp_idx++] = c;
    tmp[tmp_idx] = '\0';
};

/**
 * @brief Device function to parse a VCF line
 *
 * @param params [in,out] Kernel parameters containing all data structures
 * @param my_mem [in] Temporary buffer for string operations
 * @param currBatch [in] Current batch index
 * @param batch_size [in] Number of lines per batch
 * @param hasSamp [in] Whether sample data is present
 *
 * @details Processes VCF fields in sequence:
 *  1. Basic fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER)
 *  2. INFO field key-value pairs
 *  3. FORMAT field specifications (if samples present)
 *  4. Sample data according to FORMAT (if samples present)
 *
 * @note Uses temporary buffers in my_mem for string operations
 * @warning Assumes my_mem size >= thID*MAX_TOKEN_LEN*MAX_TOKENS*3
 */
__device__ void get_vcf_line(KernelParams* params, char* my_mem, int currBatch, int batch_size, bool hasSamp){
    long thID =  threadIdx.x + blockIdx.x * blockDim.x;
    bool find1 = false;
    long iter=0;
    int tmp_idx;
    int num_sample_tokens;
    char tmp[MAX_TMP_LEN]; 
    char* tmp_split = (char*) &my_mem[thID*MAX_TOKEN_LEN*MAX_TOKENS*3];
    char* tmp_values = (char*) &my_mem[thID*MAX_TOKEN_LEN*MAX_TOKENS*3 + MAX_TOKEN_LEN*MAX_TOKENS];
    char* sub_split = (char*) &my_mem[thID*MAX_TOKEN_LEN*MAX_TOKENS*3 + 2*MAX_TOKEN_LEN*MAX_TOKENS];
    char key[MAX_TOKEN_LEN], value[MAX_TOKEN_LEN];

    thID = thID + currBatch*batch_size;

    if(thID>=params->numLines){
        return;
    }
    long start = params->new_lines_index[thID];

    
    //Var Number
    params->var_number[thID] = thID;

    //Chromosome - CPU
    while (__ldg(&params->line[start + iter]) != '\t' && __ldg(&params->line[start + iter]) != ' ') {
        iter++;
    }
    iter++;

    //Position
    reset_tmp(tmp, tmp_idx);
    while (__ldg(&params->line[start + iter]) != '\t' && __ldg(&params->line[start + iter]) != ' ') {
        append_tmp(tmp, tmp_idx, __ldg(&params->line[start + iter]));
        iter++;
    }
    iter++;
    params->pos[thID] = cuda_atol(tmp);

    //ID - CPU
    reset_tmp(tmp, tmp_idx);
    find1=false;
    while (__ldg(&params->line[start + iter]) != '\t' && __ldg(&params->line[start + iter]) != ' ') {
        iter++;
    }
    iter++;

    //Reference - CPU
    reset_tmp(tmp, tmp_idx);
    find1=false;
    while (__ldg(&params->line[start + iter]) != '\t' && __ldg(&params->line[start + iter]) != ' ') {
        iter++;
    }
    iter++;

    //Alternative - CPU
    reset_tmp(tmp, tmp_idx);
    while (__ldg(&params->line[start + iter]) != '\t' && __ldg(&params->line[start + iter]) != ' ') {
        iter++;
    }
    iter++;

    //Quality
    reset_tmp(tmp, tmp_idx);
    while (__ldg(&params->line[start + iter]) != '\t' && __ldg(&params->line[start + iter]) != ' ' && __ldg(&params->line[start + iter]) != '\n') {
        append_tmp(tmp, tmp_idx, __ldg(&params->line[start + iter]));
        ++iter;
    }
    ++iter;
    params->qual[thID] = (cuda_strncmp(tmp, ".", MAX_TMP_LEN) == 0) ? __float2half(-1.0f) : safeStof(tmp);
    
    //Filter
    reset_tmp(tmp, tmp_idx);
    while (__ldg(&params->line[start + iter]) != '\t' && __ldg(&params->line[start + iter]) != ' ') {
        iter++;
    }
    iter++;

    // Info field (semicolon-separated key-value pairs)
    reset_tmp(tmp, tmp_idx);
    while (__ldg(&params->line[start + iter]) != '\t' && __ldg(&params->line[start + iter]) != '\n') {
        append_tmp(tmp, tmp_idx, __ldg(&params->line[start + iter]));
        ++iter;
    }
    ++iter;

    int num_info_tokens = split(tmp, ';', tmp_split);
    for (int i = 0; i < num_info_tokens; ++i) {
        char *key_value = &tmp_split[MAX_TOKEN_LEN*i];
        int j = 0;
        while (key_value[j] != '=' && key_value[j] != '\0') {
            key[j] = key_value[j];
            ++j;
        }
        key[j] = '\0';
        if (key_value[j] == '=') {
            cuda_strncpy(value, key_value + j + 1, MAX_TOKEN_LEN);
        } else {
            value[0] = '\0';
        }

        int type = getValueFromKeyMap1(key);

        if (type == INT) {
            // Process INT type
            int el = 0;
            while (cuda_strncmp(&(params->int_name[el*16]), key, MAX_TOKEN_LEN) != 0 && el < NUM_KEYS_MAP1) ++el;
            if (el < NUM_KEYS_MAP1) {
                if(cuda_strncmp(&(params->int_name[el*16]), "TSA", MAX_TOKEN_LEN)==0){
                    if(!cuda_strncmp(value, "SNV", MAX_TOKEN_LEN)){
                        params->in_int[el*params->numLines+thID] = 0;
                    }else if(!cuda_strncmp(value, "INS", MAX_TOKEN_LEN) || !cuda_strncmp(value, "insertion", MAX_TOKEN_LEN)){
                        params->in_int[el*params->numLines+thID] = 1;
                    }else if(!cuda_strncmp(value, "DEL", MAX_TOKEN_LEN) || !cuda_strncmp(value, "deletion", MAX_TOKEN_LEN)){
                        params->in_int[el*params->numLines+thID] = 2;
                    }else if(!cuda_strncmp(value, "INV", MAX_TOKEN_LEN) || !cuda_strncmp(value, "inversion", MAX_TOKEN_LEN)){
                        params->in_int[el*params->numLines+thID] = 3;
                    }else{
                        params->in_int[el*params->numLines+thID] = 4;
                    }
                }else{
                    params->in_int[el*params->numLines+thID] = cuda_atoi(value);                   
                }  
            }
        } else if (type == FLOAT) {
            // Process FLOAT type
            int el = 0;
            while (cuda_strncmp(&(params->float_name[el*16]), key, MAX_TOKEN_LEN) != 0 && el < NUM_KEYS_MAP1) ++el;
            if (el < NUM_KEYS_MAP1) {
                params->in_float[el*params->numLines+thID] = safeStof(value);
            }
        } else if (type == FLAG) {
            // Process FLAG type
            int el = 0;
            while (cuda_strncmp(&(params->flag_name[el*16]), key, 16) != 0 && el < NUM_KEYS_MAP1) ++el;
            //cuda_strncmp
            if (el < NUM_KEYS_MAP1) {
                params->in_flag[el*params->numLines+thID] = 1;
            }
        }
    }

    if(!hasSamp){
        return;
    }

    //Getting the format fields
    reset_tmp(tmp, tmp_idx);
    while (__ldg(&params->line[start + iter]) != '\t' && __ldg(&params->line[start + iter]) != '\n') {
        append_tmp(tmp, tmp_idx, __ldg(&params->line[start + iter]));
        ++iter;
    }
    ++iter;
    
    //Qui ho GT:AD
    num_sample_tokens = split(tmp, ':', tmp_values);

    // Process each sample
    for (int samp = 0; samp < params->numSample; samp++) {
        reset_tmp(tmp, tmp_idx);
        find1 = false;
        while (!find1) {
            if (__ldg(&params->line[start + iter]) == '\t' || __ldg(&params->line[start + iter]) == ' ' || __ldg(&params->line[start + iter]) == '\n') {
                find1 = true;
                iter++;

                if(cuda_strncmp(tmp, "./.", MAX_TMP_LEN)==0 || cuda_strncmp(tmp, ".|.", MAX_TMP_LEN)==0){
                    params->samp_var_id[thID * params->numSample + samp] = thID;
                    params->samp_id[thID * params->numSample + samp] = static_cast<unsigned short>(samp);
                    params->sample_GT[thID * params->numSample + samp] = getValueFromKeyGT(tmp);
                    continue;
                }

                split(tmp, ':', tmp_split);
                for (int j = 0; j < num_sample_tokens; j++) {
                    bool find_type = false;
                    bool find_elem = false;
                    while (!find_type) {
                        if (cuda_strncmp(&tmp_values[MAX_TOKEN_LEN*j], "GT", MAX_TOKEN_LEN) == 0) {
                            // Process GT (Genotype)
                            params->samp_var_id[thID * params->numSample + samp] = thID;
                            params->samp_id[thID * params->numSample + samp] = static_cast<unsigned short>(samp);
                            if (params->numGT > 1) {
                                int num_gt_tokens = split(&tmp_split[MAX_TOKEN_LEN*j], ',', sub_split);

                                for (int k = 0; k < params->numGT; k++) {
                                    params->sample_GT[(k*params->numLines)+(thID*params->numSample)+samp] = getValueFromKeyGT(&sub_split[MAX_TOKEN_LEN*k]);
                                }
                            } else if (params->numGT == 1) {
                                params->sample_GT[thID * params->numSample + samp] = getValueFromKeyGT(&tmp_split[MAX_TOKEN_LEN*j]);
                            }
                            find_type = true;

                        } else if (getValueFromKeyMap1(&tmp_values[MAX_TOKEN_LEN*j]) == INT_FORMAT) {
                            // Process INT format
                            params->samp_var_id[thID * params->numSample + samp] = thID;
                            params->samp_id[thID * params->numSample + samp] = samp;
                            int el = 0;                     
                            while (!find_elem) {
                                if (cuda_strncmp_custom(&tmp_values[MAX_TOKEN_LEN*j], &(params->samp_int_name[el*16]), MAX_TOKEN_LEN) == 0) {    
                                    //printf("el = %d\n", el);                                
                                    if (params->samp_int_numb[el] == 1) {
                                        params->samp_int[(el*params->numLines)+(thID*params->numSample)+samp] = cuda_atoi(&tmp_split[MAX_TOKEN_LEN*j]);
                                    } else {
                                        int num_int_tokens = split(&tmp_split[MAX_TOKEN_LEN*j], ',', sub_split);
                                        for (int i = 0; i < params->samp_int_numb[el]; i++) {
                                            params->samp_int[((el + i) * params->numLines * params->numSample) + (thID * params->numSample + samp)]= cuda_atoi(&sub_split[MAX_TOKEN_LEN*i]);
                                        }
                                    }
                                    find_elem = true;
                                }
                                el++;
                            }
                            find_type = true;
                        } else if (getValueFromKeyMap1(&tmp_values[MAX_TOKEN_LEN*j]) == FLOAT_FORMAT) {
                            // Process FLOAT format
                            params->samp_var_id[thID * params->numSample + samp] = thID;
                            params->samp_id[thID * params->numSample + samp] = samp;
                            int el = 0;
                            while (!find_elem) {
                                if (cuda_strncmp_custom(&tmp_values[MAX_TOKEN_LEN*j], &(params->samp_float_name[el*16]), MAX_TOKEN_LEN) == 0) {
                                    if (params->samp_float_numb[el] == 1) {
                                        params->samp_float[(el*params->numLines)+(thID*params->numSample)+samp] = safeStof(&tmp_split[MAX_TOKEN_LEN*j]);
                                    } else {
                                        int num_float_tokens = split(&tmp_split[MAX_TOKEN_LEN*j], ',', sub_split);
                                        for (int i = 0; i < params->samp_float_numb[el]; i++) {
                                            params->samp_float[((el+i)*params->numLines) + ((thID*params->numSample)+samp)] = safeStof(&sub_split[MAX_TOKEN_LEN*i]);
                                        }
                                    }
                                    find_elem = true;
                                }
                                el++;
                            }
                            find_type = true;
                        } else {
                            find_type = true;
                        }
                    }
                }
            } else {
                append_tmp(tmp, tmp_idx, __ldg(&params->line[start + iter]));
                iter++;
            }
        }
    }

}

/**
 * @brief Main parsing kernel that processes VCF lines in batches
 *
 * @param params [in,out] Structure containing all parsing parameters and data
 * @param my_mem [in] Device memory for temporary string operations
 * @param batch_size [in] Number of lines to process per batch
 * @param hasSamp [in] Flag indicating presence of sample data
 *
 * @details 
 *  1. Calculates required iterations based on total lines and batch size
 *  2. Processes each batch using get_vcf_line
 *  3. Handles memory management for temporary operations
 *
 * @note Uses batch processing to handle large files efficiently
 * @warning Ensure my_mem is large enough for all concurrent threads
 */
__global__ void kernel (KernelParams* params, char* my_mem, int batch_size, bool hasSamp)
{
    int num_iteration = (params->numLines + batch_size - 1)/batch_size;
    for(int i=0; i<num_iteration; i++){
        get_vcf_line(params, my_mem, i, batch_size, hasSamp);
    }
    
}


#endif