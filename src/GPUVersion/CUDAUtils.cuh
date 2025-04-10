/**
 * @file CUDAUtils.cuh
 * @brief Provides device-side string manipulation and numeric conversion utilities, as well as
 *        constant/device memory mappings for specific keys and values.
 */

#ifndef CUDA_UTILS_CUH
#define CUDA_UTILS_CUH

#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <stddef.h>

using namespace std;

/// Maximum number of keys in the GT map.
#define NUM_KEYS_GT 244
/// Maximum length for each key in the GT map.
#define MAX_KEY_LENGTH_GT 5
/// Maximum number of tokens to split.
#define MAX_TOKENS 4
/// Maximum length for each token.
#define MAX_TOKEN_LEN 32
/// MAximum length for the temporary string
#define MAX_TMP_LEN 128

/**
 * @brief Constant memory holding the GT keys.
 * Each key is up to MAX_KEY_LENGTH_GT characters long.
 */
__constant__ char d_keys_gt[NUM_KEYS_GT][MAX_KEY_LENGTH_GT];

/**
  * @brief Constant memory holding the values corresponding to the GT keys.
  * Each value is a single char for each key.
  */
__constant__ char d_values_gt[NUM_KEYS_GT];
 
/// Maximum number of keys in Map1.
#define NUM_KEYS_MAP1 128
/// Maximum length for each key in Map1.
#define MAX_KEY_LENGTH_MAP1 32
 
/**
  * @brief Device memory holding the keys for Map1.
  * Each key can be up to MAX_KEY_LENGTH_MAP1 characters long.
  */
__device__ char d_keys_map1[NUM_KEYS_MAP1][MAX_KEY_LENGTH_MAP1];
 
/**
  * @brief Device memory holding the integer values corresponding to Map1 keys.
  */
__device__ int d_values_map1[NUM_KEYS_MAP1];

/**
 * @brief Compares two strings (s1 and s2) up to n characters.
 * @param s1 Pointer to the first string.
 * @param s2 Pointer to the second string.
 * @param n Maximum number of characters to compare.
 * @return 0 if the first n characters match, or a negative/positive value indicating s1 < s2 or s1 > s2.
 */
__device__ int cuda_strncmp(const char *s1, const char *s2, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (s1[i] != s2[i] || s1[i] == 0) {
            return (unsigned char)s1[i] - (unsigned char)s2[i];
        }
    }
    return 0;
}

/**
 * @brief Custom comparison of two strings (s1 and s2) up to n characters, allowing trailing digits in s2.
 *
 * Similar to cuda_strncmp, but if s1 ends before s2, it checks whether the remaining
 * characters in s2 are digits ('0'-'9') to consider the strings as matching.
 * @param s1 Pointer to the first string.
 * @param s2 Pointer to the second string.
 * @param n Maximum number of characters to compare.
 * @return 0 if the strings are considered equal, otherwise the difference between non-matching characters.
 */
__device__ int cuda_strncmp_custom(const char *s1, const char *s2, size_t n) {
    size_t i = 0;
    for (; i < n; ++i) {
        if (s1[i] == '\0') {
            break;
        }
        if (s2[i] == '\0') {
            return (unsigned char)s1[i];
        }
        if (s1[i] != s2[i]) {
            return (unsigned char)s1[i] - (unsigned char)s2[i];
        }
    }
    // If s1 ended, check if remaining characters in s2 are digits
    for (; i < n && s2[i] != '\0'; ++i) {
        char c = s2[i];
        if (c < '0' || c > '9') {
            return (unsigned char)c; 
        }
    }
    return 0;
}

/**
 * @brief Converts a C-string to an integer (similar to atoi).
 * @param str Pointer to the input string.
 * @return The integer value represented by the string.
 */
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

/**
 * @brief Converts a C-string to a float (similar to atof).
 * @param str Pointer to the input string.
 * @return The float value represented by the string.
 */
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

/**
 * @brief Converts a C-string to a long (similar to atol).
 * @param str Pointer to the input string.
 * @return The long value represented by the string.
 */
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

/**
 * @brief Copies up to n characters from src to dest, padding with '\0' if needed.
 * @param dest Pointer to the destination buffer.
 * @param src Pointer to the source string.
 * @param n Maximum number of characters to copy.
 */
__device__ void cuda_strncpy(char *dest, const char *src, size_t n) {
    size_t i;
    for (i = 0; i < n && src[i] != '\0'; i++) {
        dest[i] = src[i];
    }
    for (; i < n; i++) {
        dest[i] = '\0';
    }
}

/**
 * @brief Converts a C-string to an unsigned long.
 * @param str Pointer to the input string.
 * @return The unsigned long value represented by the string.
 */
__device__ unsigned long cuda_stoul(const char *str) {
    unsigned long num = 0;
    while (*str >= '0' && *str <= '9') { // Ensure only digits are processed
        num = num * 10 + (*str - '0');   // Convert ASCII character to numeric value
        str++;
    }
    return num;
}

/**
 * @brief Retrieves the value corresponding to a GT key from d_keys_gt and d_values_gt.
 * @param key Pointer to the key string.
 * @return The associated value if found, or -1 otherwise.
 */
__device__ int getValueFromKeyGT(const char* key) {
    for (int i = 0; i < NUM_KEYS_GT; ++i) {
        if (cuda_strncmp(key, d_keys_gt[i], MAX_KEY_LENGTH_GT) == 0) {
            return d_values_gt[i];
        }
    }
    return -1;
}

/**
 * @brief Retrieves the value corresponding to a key in d_keys_map1 and d_values_map1.
 * @param key Pointer to the key string.
 * @return The associated value if found, or -1 otherwise.
 */
__device__ int getValueFromKeyMap1(const char* key) {
    for (int i = 0; i < NUM_KEYS_MAP1; ++i) {
        if (cuda_strncmp_custom(key, d_keys_map1[i], MAX_KEY_LENGTH_MAP1) == 0) {
            return d_values_map1[i];
        }
    }
    return -1;
}

/**
 * @brief Splits a string by a given delimiter, storing tokens in split_array.
 *
 * Each token is stored in a block of size MAX_TOKEN_LEN. The function stops splitting
 * when either the string ends or the maximum number of tokens (MAX_TOKENS) is reached.
 * @param str Pointer to the input string to split.
 * @param delimiter Delimiter character.
 * @param split_array Pointer to the buffer where tokens will be stored.
 * @return The number of tokens extracted.
 */
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

/**
 * @brief Converts a string to half precision (16-bit) after checking if the string is valid.
 *
 * If the string contains invalid characters (non-digit, non-decimal point, or non-minus),
 * it returns 0.0f as a default value.
 * @param tmp Pointer to the input string.
 * @return The corresponding half precision value, or 0.0f if invalid.
 */
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

//TODO: da controllare:
/// Numero di chiavi per la mappa PolyPhen.
#define NUM_KEYS_POLYPHEN 4
/// Lunghezza massima di ogni chiave per la mappa PolyPhen (includendo il terminatore null).
#define MAX_KEY_LENGTH_POLYPHEN 20

/// Numero di chiavi per la mappa CSQ.
#define NUM_KEYS_CSQ 32
/// Lunghezza massima di ogni chiave per la mappa CSQ (includendo il terminatore null).
#define MAX_KEY_LENGTH_CSQ 32

/**
 * @brief Memoria costante contenente le chiavi per la mappa PolyPhen.
 * Ogni chiave è una stringa lunga fino a MAX_KEY_LENGTH_POLYPHEN caratteri.
 */
__constant__ char d_polyphen_keys[NUM_KEYS_POLYPHEN][MAX_KEY_LENGTH_POLYPHEN] = {
    "",
    "benign",
    "possibly_damaging",
    "probably_damaging"
};

/**
 * @brief Memoria costante contenente i valori corrispondenti alle chiavi PolyPhen.
 * Ogni valore è rappresentato da un singolo char.
 */
__constant__ char d_polyphen_values[NUM_KEYS_POLYPHEN] = {0, 1, 2, 3};

/**
 * @brief Memoria costante contenente le chiavi per la mappa CSQ.
 * Ogni chiave è una stringa lunga fino a MAX_KEY_LENGTH_CSQ caratteri.
 */
__constant__ char d_csq_keys[NUM_KEYS_CSQ][MAX_KEY_LENGTH_CSQ] = {
    "synonymous_variant",
    "missense_variant",
    "nonsense_variant",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "splice_region_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "regulatory_region_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "coding_sequence_variant",
    "non_coding_transcript_variant",
    "mature_miRNA_variant",
    "ncRNA_exon_variant",
    "ncRNA_intronic_variant",
    "protein_altering_variant",
    "transcript_ablation",
    "transcript_amplification",
    "intron_variant",
    "intergenic_variant",
    "enhancer_variant",
    "regulatory_region_variant",
    "non_coding_exon_variant",
    "pseudogene_variant"
};

/**
 * @brief Memoria costante contenente i valori corrispondenti alle chiavi CSQ.
 * Ogni valore è rappresentato da un singolo char.
 */
__constant__ char d_csq_values[NUM_KEYS_CSQ] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27
};




#endif