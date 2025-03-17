#ifndef CUDA_UTILS_CUH
#define CUDA_UTILS_CUH

#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <stddef.h>

using namespace std;

//GTMap
#define NUM_KEYS_GT 242
#define MAX_KEY_LENGTH_GT 5
#define MAX_TOKENS 4 // Maximum number of tokens to split
#define MAX_TOKEN_LEN 16 // Maximum length of each token

__constant__ char d_keys_gt[NUM_KEYS_GT][MAX_KEY_LENGTH_GT];   // Chiavi nella constant memory
__constant__ char d_values_gt[NUM_KEYS_GT];                    // Valori nella constant memory

#define NUM_KEYS_MAP1 128
#define MAX_KEY_LENGTH_MAP1 32

__device__ char d_keys_map1[NUM_KEYS_MAP1][MAX_KEY_LENGTH_MAP1];   
__device__ int d_values_map1[NUM_KEYS_MAP1];                       

//TODO - queste funzioni sono simili controlla utilità

//Usata per GT e funzia
__device__ int cuda_strncmp(const char *s1, const char *s2, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (s1[i] != s2[i] || s1[i] == 0) {
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

__device__ int cuda_strncmp_custom(const char *s1, const char *s2, size_t n) {
    size_t i = 0;
    //printf("s1= %s\n", s1);
    //printf("s2= %s\n", s2);
    // Confronta fino a quando s1 non termina
    for (; i < n; ++i) {
        if (s1[i] == '\0') {
            // s1 è terminata, esci dal ciclo
            break;
        }
        if (s2[i] == '\0') {
            // s2 è troppo corto rispetto a s1: non match
            return (unsigned char)s1[i];
        }
        if (s1[i] != s2[i]) {
            return (unsigned char)s1[i] - (unsigned char)s2[i];
        }
    }
    // s1 è terminata: controlla che eventuali caratteri in più in s2 siano cifre
    for (; i < n && s2[i] != '\0'; ++i) {
        char c = s2[i];
        if (c < '0' || c > '9') {
            // Se il carattere non è una cifra, le stringhe non sono considerate uguali
            return (unsigned char)c; 
        }
    }
    return 0;
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

__device__ int getValueFromKeyMap1(const char* key) {
    for (int i = 0; i < NUM_KEYS_MAP1; ++i) {
        if (cuda_strncmp_custom(key, d_keys_map1[i], MAX_KEY_LENGTH_MAP1) == 0) {
            return d_values_map1[i];
        }
    }
    return -1;
}

// Helper function to split a string using a delimiter
__device__ int split(const char *str, char delimiter, char* split_array) {
    int token_count = 0;
    int token_idx = 0;

    for (int i = 0; str[i] != '\0'; ++i) {
        if (str[i] == delimiter) {
            split_array[token_count * MAX_TOKEN_LEN + token_idx] = '\0'; // Null-terminate the token
            ++token_count;
            token_idx = 0;
            if (token_count >= MAX_TOKENS) break; // Avoid overflow
        } else {
            if (token_idx < MAX_TOKEN_LEN - 1) {
                split_array[token_count * MAX_TOKEN_LEN + (token_idx++)] = str[i];
            }
        }
    }
    split_array[token_count * MAX_TOKEN_LEN + token_idx] = '\0'; // Null-terminate the last token
    return token_count + 1;
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

#endif