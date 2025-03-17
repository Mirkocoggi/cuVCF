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
/*
__global__ void get_vcf_line_kernel(char *line, unsigned int *var_number, unsigned int *pos, __half *qual,
            __half *in_float, bool *in_flag, int *in_int, char* float_name, char* flag_name, char* int_name, unsigned int *new_lines_index, unsigned int numLines)
{
    long thID =  threadIdx.x + blockIdx.x * blockDim.x;
    if(thID>=numLines){
        return;
    }

    bool find1 = false;
    long iter=0;
    //string tmp="\0";
    char tmp[MAX_TOKEN_LEN];
    char tmp_split[MAX_TOKENS][MAX_TOKEN_LEN];
    //vector<string> tmp_split;
    //vector<string> tmp_format_split;
    //vector<string> tmp_subSplit;
    long start = new_lines_index[thID];
    int tmp_idx;

    // Reset temp string
    auto reset_tmp = [&]() {
        tmp_idx = 0;
        tmp[0] = '\0';
    };

    // Append character to temp string - TODO - Check corretta
    auto append_tmp = [&](char c) {
        if (tmp_idx < MAX_TOKEN_LEN) tmp[tmp_idx++] = c;
        tmp[tmp_idx] = '\0';
    };
    
    //Var Number
    var_number[thID] = thID;

    //Chromosome - CPU
    while (__ldg(&line[start + iter]) != '\t' && __ldg(&line[start + iter]) != ' ') {
        iter++;
    }
    iter++;

    //Position
    reset_tmp();
    while (__ldg(&line[start + iter]) != '\t' && __ldg(&line[start + iter]) != ' ') {
        append_tmp(__ldg(&line[start + iter]));
        iter++;
    }
    iter++;
    pos[thID] = cuda_atol(tmp);
    
    reset_tmp();
    find1=false;
    while(!find1){
        if(__ldg(&line[start+iter])=='\t'||__ldg(&line[start+iter])==' '){
            find1 = true;
            iter++;
            pos[thID] = cuda_stoul(tmp);
        }else{
            append_tmp(__ldg(&line[start+iter]));
            iter++;
        }
    }

    //ID - CPU
    reset_tmp();
    find1=false;
    while (__ldg(&line[start + iter]) != '\t' && __ldg(&line[start + iter]) != ' ') {
        iter++;
    }
    iter++;

    //Reference - CPU
    reset_tmp();
    find1=false;
    while (__ldg(&line[start + iter]) != '\t' && __ldg(&line[start + iter]) != ' ') {
        iter++;
    }
    iter++;

    //Alternative - CPU
    reset_tmp();
    find1=false;
    while (__ldg(&line[start + iter]) != '\t' && __ldg(&line[start + iter]) != ' ') {
        iter++;
    }
    iter++;

    //Quality
    reset_tmp();
    while (__ldg(&line[start + iter]) != '\t' && __ldg(&line[start + iter]) != ' ' && __ldg(&line[start + iter]) != '\n') {
        append_tmp(__ldg(&line[start + iter]));
        ++iter;
    }
    ++iter;
    qual[thID] = (cuda_strcmp(tmp, ".") == 0) ? __float2half(0.0f) : safeStof(tmp);
    
    //Filter
    reset_tmp();
    find1=false;
    while (__ldg(&line[start + iter]) != '\t' && __ldg(&line[start + iter]) != ' ') {
        iter++;
    }
    iter++;

    // Info field (semicolon-separated key-value pairs)
    reset_tmp();
    while (__ldg(&line[start + iter]) != '\t' && __ldg(&line[start + iter]) != '\n') {
        append_tmp(__ldg(&line[start + iter]));
        ++iter;
    }
    ++iter;
    
    int num_info_tokens = split(tmp, ';', tmp_split);
    for (int i = 0; i < num_info_tokens; ++i) {
        char *key_value = tmp_split[MAX_TOKEN_LEN*i];
        char key[MAX_TOKEN_LEN], value[MAX_TOKEN_LEN];
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
            while (cuda_strcmp(&int_name[el*16], key) != 0 && el < NUM_KEYS_MAP1) ++el;
            if (el < NUM_KEYS_MAP1) {
                if(cuda_strcmp(&int_name[el*16], "TSA")){ //TODO tmp implementation per TSA
                    if(!cuda_strcmp(value, "SNV")){
                        in_int[el*numLines+thID] = 0;
                    }else if(!cuda_strcmp(value, "INS") || !cuda_strcmp(value, "insertion")){
                        in_int[el*numLines+thID] = 1;
                    }else if(!cuda_strcmp(value, "DEL") || !cuda_strcmp(value, "deletion")){
                        in_int[el*numLines+thID] = 2;
                    }else if(!cuda_strcmp(value, "INV") || !cuda_strcmp(value, "inversion")){
                        in_int[el*numLines+thID] = 3;
                    }else{
                        in_int[el*numLines+thID] = 4;
                    }
                }else{
                    in_int[el*numLines+thID] = cuda_atoi(value);
                }                
            }
        } else if (type == FLOAT) {
            // Process FLOAT type
            int el = 0;
            while (cuda_strcmp(&float_name[el*16], key) != 0 && el < NUM_KEYS_MAP1) ++el;
            if (el < NUM_KEYS_MAP1) {
                in_float[el*numLines+thID] = safeStof(value);
            }
        } else if (type == FLAG) {
            // Process FLAG type
            int el = 0;
            while (cuda_strcmp(&flag_name[el*16], key) != 0 && el < NUM_KEYS_MAP1) ++el;
            if (el < NUM_KEYS_MAP1) {
                in_flag[el*numLines+thID] = 1;
            }
        }
    }
}
*/

// Reset temp string
__device__ void reset_tmp(char* tmp, int& tmp_idx) {
    tmp_idx = 0;
    tmp[0] = '\0';
};

// Append character to temp string
__device__ void append_tmp(char* tmp, int& tmp_idx, char c) {
    if (tmp_idx < MAX_TOKEN_LEN-1) tmp[tmp_idx++] = c;
    tmp[tmp_idx] = '\0';
};

__global__ void get_vcf_line_format_kernel(KernelParams* params, char* my_mem)
{
    long thID =  threadIdx.x + blockIdx.x * blockDim.x;
    if(thID>=params->numLines){
        return;
    }


    bool find1 = false;
    long iter=0;
    //string tmp="\0"; TODO - sistemare tutte le dimensioni necessarie
    char tmp[MAX_TOKEN_LEN]; 
    //char tmp_split[MAX_TOKENS][MAX_TOKEN_LEN];
    //char tmp_values[MAX_TOKENS][MAX_TOKEN_LEN];
    //char sub_split[MAX_TOKENS][MAX_TOKEN_LEN];

    char* tmp_split = (char*) &my_mem[thID*MAX_TOKEN_LEN*MAX_TOKENS*3];
    char* tmp_values = (char*) &my_mem[thID*MAX_TOKEN_LEN*MAX_TOKENS*3 + MAX_TOKEN_LEN*MAX_TOKENS];
    char* sub_split = (char*) &my_mem[thID*MAX_TOKEN_LEN*MAX_TOKENS*3 + 2*MAX_TOKEN_LEN*MAX_TOKENS];

    long start = params->new_lines_index[thID];
    int tmp_idx;
    int num_sample_tokens;
    
    //Var Number
    params->var_number[thID] = thID;

    //Chromosome - CPU
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        iter++;
    }
    iter++;

    //Position
    reset_tmp(tmp, tmp_idx);
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        append_tmp(tmp, tmp_idx, params->line[start + iter]);
        iter++;
    }
    iter++;
    params->pos[thID] = cuda_atol(tmp);

    reset_tmp(tmp, tmp_idx);
    find1=false;
    while(!find1){
        if(params->line[start+iter]=='\t'||params->line[start+iter]==' '){
            find1 = true;
            iter++;
            params->pos[thID] = cuda_stoul(tmp);
        }else{
            append_tmp(tmp, tmp_idx, params->line[start+iter]);
            iter++;
        }
    }

    //ID - CPU
    reset_tmp(tmp, tmp_idx);
    find1=false;
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        iter++;
    }
    iter++;

    //Reference - CPU
    reset_tmp(tmp, tmp_idx);
    find1=false;
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        iter++;
    }
    iter++;

    //Alternative - CPU
    reset_tmp(tmp, tmp_idx);
    find1=false;
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        iter++;
    }
    iter++;

    //Quality
    reset_tmp(tmp, tmp_idx);
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ' && params->line[start + iter] != '\n') {
        append_tmp(tmp, tmp_idx, params->line[start + iter]);
        ++iter;
    }
    ++iter;
    params->qual[thID] = (cuda_strcmp(tmp, ".") == 0) ? __float2half(0.0f) : safeStof(tmp);
    
    //Filter
    reset_tmp(tmp, tmp_idx);
    find1=false;
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        iter++;
    }
    iter++;

    // Info field (semicolon-separated key-value pairs)
    reset_tmp(tmp, tmp_idx);
    while (params->line[start + iter] != '\t' && params->line[start + iter] != '\n') {
        append_tmp(tmp, tmp_idx, params->line[start + iter]);
        ++iter;
    }
    ++iter;
    
    int num_info_tokens = split(tmp, ';', tmp_split);
    char key[MAX_TOKEN_LEN], value[MAX_TOKEN_LEN];
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
            while (cuda_strcmp(&(params->int_name[el*16]), key) != 0 && el < NUM_KEYS_MAP1) ++el;
            if (el < NUM_KEYS_MAP1) {
                params->in_int[el*params->numLines+thID] = cuda_atoi(value);
            }
        } else if (type == FLOAT) {
            // Process FLOAT type
            int el = 0;
            while (cuda_strcmp(&(params->float_name[el*16]), key) != 0 && el < NUM_KEYS_MAP1) ++el;
            if (el < NUM_KEYS_MAP1) {
                params->in_float[el*params->numLines+thID] = safeStof(value);
            }
        } else if (type == FLAG) {
            // Process FLAG type
            int el = 0;
            while (cuda_strcmp(&(params->flag_name[el*16]), key) != 0 && el < NUM_KEYS_MAP1) ++el;
            if (el < NUM_KEYS_MAP1) {
                params->in_flag[el*params->numLines+thID] = 1;
            }
        }
    }

    //Qui ho GT:AD
    num_sample_tokens = split(tmp, ':', tmp_values);

    // Process each sample
    for (int samp = 0; samp < params->numSample; samp++) {
        reset_tmp(tmp, tmp_idx);
        find1 = false;
        while (!find1) {
            if (params->line[start + iter] == '\t' || params->line[start + iter] == ' ' || params->line[start + iter] == '\n') {
                find1 = true;
                iter++;

                if(cuda_strcmp(tmp, "./.")==0 || cuda_strcmp(tmp, ".|.")==0){
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
                        if (cuda_strcmp(&tmp_values[MAX_TOKEN_LEN*j], "GT") == 0) {
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
                append_tmp(tmp, tmp_idx, params->line[start + iter]);
                iter++;
            }
        }
    }

}

#endif