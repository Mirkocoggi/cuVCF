#ifndef KERNELS_H
#define KERNELS_H

#include "CUDAUtils.cuh"
#include "Utils.h"

#include <cuda_runtime.h>     
#include <cuda_fp16.h>  
#include <thrust/device_ptr.h> 
#include <thrust/sort.h>
#include <boost/algorithm/string.hpp>
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

__global__ static void get_vcf_line_kernel(char *line, unsigned int *var_number, unsigned int *pos, __half *qual,
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

    // Append character to temp string
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
        char *key_value = tmp_split[i];
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

__global__ static void get_vcf_line_format_kernel(KernelParams* params)
{
    long thID =  threadIdx.x + blockIdx.x * blockDim.x;
    if(thID>=params->numLines){
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
    long start = params->new_lines_index[thID];
    int tmp_idx;

    // Reset temp string
    auto reset_tmp = [&]() {
        tmp_idx = 0;
        tmp[0] = '\0';
    };

    // Append character to temp string
    auto append_tmp = [&](char c) {
        if (tmp_idx < MAX_TOKEN_LEN) tmp[tmp_idx++] = c;
        tmp[tmp_idx] = '\0';
    };
    
    //Var Number
    params->var_number[thID] = thID;

    //Chromosome - CPU
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        iter++;
    }
    iter++;

    //Position
    reset_tmp();
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        append_tmp(params->line[start + iter]);
        iter++;
    }
    iter++;
    params->pos[thID] = cuda_atol(tmp);

    reset_tmp();
    find1=false;
    while(!find1){
        if(params->line[start+iter]=='\t'||params->line[start+iter]==' '){
            find1 = true;
            iter++;
            params->pos[thID] = cuda_stoul(tmp);
        }else{
            append_tmp(params->line[start+iter]);
            iter++;
        }
    }

    //ID - CPU
    reset_tmp();
    find1=false;
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        iter++;
    }
    iter++;

    //Reference - CPU
    reset_tmp();
    find1=false;
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        iter++;
    }
    iter++;

    //Alternative - CPU
    reset_tmp();
    find1=false;
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        iter++;
    }
    iter++;

    //Quality
    reset_tmp();
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ' && params->line[start + iter] != '\n') {
        append_tmp(params->line[start + iter]);
        ++iter;
    }
    ++iter;
    params->qual[thID] = (cuda_strcmp(tmp, ".") == 0) ? __float2half(0.0f) : safeStof(tmp);
    
    //Filter
    reset_tmp();
    find1=false;
    while (params->line[start + iter] != '\t' && params->line[start + iter] != ' ') {
        iter++;
    }
    iter++;

    // Info field (semicolon-separated key-value pairs)
    reset_tmp();
    while (params->line[start + iter] != '\t' && params->line[start + iter] != '\n') {
        append_tmp(params->line[start + iter]);
        ++iter;
    }
    ++iter;
    
    int num_info_tokens = split(tmp, ';', tmp_split);
    for (int i = 0; i < num_info_tokens; ++i) {
        char *key_value = tmp_split[i];
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

    /*Info
    reset_tmp();
    find1=false;
    while(!find1){
        if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
            find1 = true;
            iter++;
            vector<string> tmp_el;
            boost::split(tmp_el, tmp, boost::is_any_of(";")); //Info arguments separation
            vector<string> tmp_elems;
            for(int ii=0; ii<tmp_el.size(); ii++){
                boost::split(tmp_elems, tmp_el[ii], boost::is_any_of("=")); //info_id separation from contents
                bool find_info_type = false;
                bool find_info_elem = false;
                if(tmp_elems.size()==2){
                    while(!find_info_type){
                        if(getValueFromKeyMap1(tmp_elems[0].c_str())==INT){
                            //Int
                            bool isAlt = false;
                            int el=0;
                            while(!find_info_elem){
                                if(in_int[el].name == tmp_elems[0]){
                                    in_int[el].i_int[thID] = cuda_atoi(tmp_elems[1]);
                                    find_info_elem = true;
                                }
                                el++; 
                            }
                            find_info_type = true;
                        }else if(getValueFromKeyMap1(tmp_elems[0].c_str())==FLOAT){
                            //Float
                            int el=0;
                            while(!find_info_elem){
                                if(in_float[el].name == tmp_elems[0]){
                                    in_float[el].i_float[thID] = safeStof(tmp_elems[1].c_str());
                                    find_info_elem = true;
                                } 
                                el++;
                            }
                            find_info_type = true;
                        }else{
                            find_info_type = true;
                        }
                    }
                }else{
                    if((getValueFromKeyMap1(tmp_elems[0].c_str())==FLAG) && cuda_strcmp(&tmp_elems[0][0],"")){
                        //Flag
                        int el=0;
                        while(!find_info_elem){
                            if(in_flag[el].name == tmp_elems[0]){
                                in_flag[el].i_flag[thID] = 1;
                                find_info_elem = true;
                            } 
                            el++;
                        }
                            
                    }
                    tmp_elems.clear();
                }
            }
        }else{
            tmp += line[start+iter];
            iter++;
        }
    }
    */

   // Format decomposition
    reset_tmp();
    find1 = false;

    // Parse Format's template
    while (!find1) {
        if (params->line[start + iter] == '\t' || params->line[start + iter] == ' ') {
            find1 = true;
            iter++;
            // Split the `tmp` buffer into format tokens using ':' as delimiter
            split(tmp, ':', tmp_split); // `split` populates `tmp_split` and updates `token_count`
        } else {
            append_tmp(params->line[start + iter]);
            iter++;
        }
    }

    // Process each sample
    for (int samp = 0; samp < params->numSample; samp++) {
        reset_tmp();
        find1 = false;

        while (!find1) {
            if (params->line[start + iter] == '\t' || params->line[start + iter] == ' ' || params->line[start + iter] == '\n') {
                find1 = true;
                iter++;
                // Split the tmp buffer into sample tokens
                int num_sample_tokens = split(tmp, ':', tmp_split);

                for (int j = 0; j < num_sample_tokens; j++) {
                    bool find_type = false;
                    bool find_elem = false;

                    while (!find_type) {
                        if (cuda_strcmp(tmp_split[j], "GT") == 0) {
                            // Process GT (Genotype)
                            params->samp_var_id[thID * params->numSample + samp] = params->var_number[thID];
                            params->samp_id[thID * params->numSample + samp] = static_cast<unsigned short>(samp);

                            if (params->numGT > 1) {
                                char sub_split[MAX_TOKENS][MAX_TOKEN_LEN];
                                int num_gt_tokens = split(tmp_split[j], ',', sub_split);

                                for (int k = 0; k < params->numGT; k++) {
                                    params->sample_GT[(k*params->numLines)+(thID*params->numSample)+samp] = getValueFromKeyGT(sub_split[k]);
                                }
                            } else if (params->numGT == 1) {
                                params->sample_GT[thID * params->numSample + samp] = getValueFromKeyGT(tmp_split[j]);
                            }
                            find_type = true;
                        } else if (getValueFromKeyMap1(tmp_split[j]) == INT_FORMAT) {
                            // Process INT format
                            params->samp_var_id[thID * params->numSample + samp] = params->var_number[thID];
                            params->samp_id[thID * params->numSample + samp] = samp;

                            int el = 0;
                            while (!find_elem) {
                                if (cuda_strncmp(&(params->samp_int_name[el*16]), tmp_split[j], MAX_TOKEN_LEN) == 0) {
                                    if (params->samp_int_numb[el] == 1) {
                                        params->samp_int[(el*params->numLines)+(thID*params->numSample)+samp] = cuda_atoi(tmp_split[j]);
                                    } else {
                                        char sub_split[MAX_TOKENS][MAX_TOKEN_LEN];
                                        int num_int_tokens = split(tmp_split[j], ',', sub_split);

                                        for (int i = 0; i < params->samp_int_numb[el]; i++) {
                                            params->samp_int[((el+i)*params->numLines) + ((thID*params->numSample)+samp)]= cuda_atoi(sub_split[i]);
                                        }
                                    }
                                    find_elem = true;
                                }
                                el++;
                            }
                            find_type = true;
                        } else if (getValueFromKeyMap1(tmp_split[j]) == FLOAT_FORMAT) {
                            // Process FLOAT format
                            params->samp_var_id[thID * params->numSample + samp] = params->var_number[thID];
                            params->samp_id[thID * params->numSample + samp] = samp;

                            int el = 0;
                            while (!find_elem) {
                                if (cuda_strncmp(&(params->samp_float_name[el*16]), tmp_split[j], MAX_TOKEN_LEN) == 0) {
                                    if (params->samp_float_numb[el] == 1) {
                                        params->samp_float[(el*params->numLines)+(thID*params->numSample)+samp] = safeStof(tmp_split[j]);
                                    } else {
                                        char sub_split[MAX_TOKENS][MAX_TOKEN_LEN];
                                        int num_float_tokens = split(tmp_split[j], ',', sub_split);

                                        for (int i = 0; i < params->samp_float_numb[el]; i++) {
                                            params->samp_float[((el+i)*params->numLines) + ((thID*params->numSample)+samp)] = safeStof(sub_split[i]);
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
                append_tmp(params->line[start + iter]);
                iter++;
            }
        }
    }
}

#endif