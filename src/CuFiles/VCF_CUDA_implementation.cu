#ifndef VCF_CUDA_H
#define VCF_CUDA_H

//#include "VCFparser_mt_col_struct_CU.h"
#include "VCF_var_columns_df_CU.h"

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

//GTMap
#define NUM_KEYS_GT 242
#define MAX_KEY_LENGTH_GT 5
#define MAX_TOKENS 10 // Maximum number of tokens to split
#define MAX_TOKEN_LEN 512 // Maximum length of each token

__constant__ char d_keys_gt[NUM_KEYS_GT][MAX_KEY_LENGTH_GT];   // Chiavi nella constant memory
__constant__ char d_values_gt[NUM_KEYS_GT];                    // Valori nella constant memory


__device__ int cuda_strncmp(const char *s1, const char *s2, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (s1[i] != s2[i] || s1[i] == '\0') {
            return (unsigned char)s1[i] - (unsigned char)s2[i];
        }
    }
    return 0;
}

__device__ int cuda_strcmp(const char *s1, const char *s2) {
    while (*s1 && (*s1 == *s2)) {
        s1++;
        s2++;
    }
    return *(unsigned char *)s1 - *(unsigned char *)s2;
}

__device__ int cuda_atoi(const char *str) {
    int num = 0;
    int sign = 1;
    if (*str == '-') {
        sign = -1;
        str++;
    }
    while (*str >= '0' && *str <= '9') {
        num = num * 10 + (*str - '0');
        str++;
    }
    return sign * num;
}

__device__ float cuda_atof(const char *str) {
    float num = 0.0f;
    float fraction = 0.1f;
    int sign = 1;
    bool is_fraction = false;

    if (*str == '-') {
        sign = -1;
        str++;
    }

    while (*str) {
        if (*str == '.') {
            is_fraction = true;
            str++;
            continue;
        }

        if (*str >= '0' && *str <= '9') {
            if (is_fraction) {
                num += (*str - '0') * fraction;
                fraction *= 0.1f;
            } else {
                num = num * 10 + (*str - '0');
            }
        } else {
            break; // Invalid character
        }
        str++;
    }

    return sign * num;
}

__device__ long cuda_atol(const char *str) {
    long num = 0;
    int sign = 1;
    if (*str == '-') {
        sign = -1;
        str++;
    }
    while (*str >= '0' && *str <= '9') {
        num = num * 10 + (*str - '0');
        str++;
    }
    return sign * num;
}

__device__ void cuda_strncpy(char *dest, const char *src, size_t n) {
    size_t i;
    for (i = 0; i < n && src[i] != '\0'; i++) {
        dest[i] = src[i];
    }
    for (; i < n; i++) {
        dest[i] = '\0';
    }
}

__device__ unsigned long cuda_stoul(const char *str) {
    unsigned long num = 0;
    while (*str >= '0' && *str <= '9') { // Ensure only digits are processed
        num = num * 10 + (*str - '0');   // Convert ASCII character to numeric value
        str++;
    }
    return num;
}

__device__ int getValueFromKeyGT(const char* key) {
    for (int i = 0; i < NUM_KEYS_GT; ++i) {
        if (cuda_strncmp(key, d_keys_gt[i], MAX_KEY_LENGTH_GT) == 0) {
            return d_values_gt[i];
        }
    }
    return -1;
}

#define NUM_KEYS_MAP1 128
#define MAX_KEY_LENGTH_MAP1 32

__device__ char d_keys_map1[NUM_KEYS_MAP1][MAX_KEY_LENGTH_MAP1];   // Chiavi nella constant memory
__device__ int d_values_map1[NUM_KEYS_MAP1];                       // Valori nella constant memory

__device__ int getValueFromKeyMap1(const char* key) {
    for (int i = 0; i < NUM_KEYS_MAP1; ++i) {
        if (cuda_strncmp(key, d_keys_map1[i], MAX_KEY_LENGTH_MAP1) == 0) {
            return d_values_map1[i];
        }
    }
    return -1;
}

// Helper function to split a string using a delimiter
__device__ int split(const char *str, char delimiter, char split_array[MAX_TOKENS][MAX_TOKEN_LEN]) {
    int token_count = 0;
    int token_idx = 0;

    for (int i = 0; str[i] != '\0'; ++i) {
        if (str[i] == delimiter) {
            split_array[token_count][token_idx] = '\0'; // Null-terminate the token
            ++token_count;
            token_idx = 0;
            if (token_count >= MAX_TOKENS) break; // Avoid overflow
        } else {
            if (token_idx < MAX_TOKEN_LEN - 1) {
                split_array[token_count][token_idx++] = str[i];
            }
        }
    }
    split_array[token_count][token_idx] = '\0'; // Null-terminate the last token
    return token_count + 1;
};

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

__device__ half safeStof(const char* tmp) {
    float value = 0.0f;
    bool is_valid = true;

    for (int i = 0; tmp[i] != '\0'; i++) {
        if (!((tmp[i] >= '0' && tmp[i] <= '9') || tmp[i] == '.' || tmp[i] == '-')) {
            is_valid = false;
            break;
        }
    }

    if (is_valid) {
        value = cuda_atof(tmp); //like stof but for device
    } else {
        value = 0.0f; // Valore predefinito in caso di errore
    }

    return __float2half(value);
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
    printf("AAAAAAAAAAAAAAA 1\n");
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