/**
 * @file DataStructures.h
 * @brief Contains fundamental data structures used for VCF parsing and GPU processing.
 *
 * This file defines several structures used to pass parameters to CUDA kernels,
 * as well as to store VCF header information and parsed data for both host and device.
 */

#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#include <chrono>
#include <boost/algorithm/string.hpp>
#include <cuda_runtime.h>     
#include <cuda_fp16.h>  
#include <iostream>
#include <map>
#include <vector>

using namespace std;

/**
 * @struct KernelParams
 * @brief Structure containing parameters and pointers for CUDA kernel execution.
 *
 * This structure holds pointers to arrays representing a VCF line,
 * variant numbers, positions, quality scores, and various info fields.
 * It also contains pointers to sample-specific data and parameters.
 */
 struct KernelParams {
    char *line;                   ///< Pointer to the VCF line.
    unsigned int *var_number;     ///< Pointer to the array of variant numbers.
    unsigned int *pos;            ///< Pointer to the array of variant positions.
    __half *qual;                 ///< Pointer to the array of quality scores (half precision).
    __half *in_float;             ///< Pointer to the array of float info values.
    bool *in_flag;                ///< Pointer to the array of flag info values.
    int *in_int;                  ///< Pointer to the array of integer info values.
    char* float_name;             ///< Pointer to the float field names.
    char* flag_name;              ///< Pointer to the flag field names.
    char* int_name;               ///< Pointer to the integer field names.
    unsigned int *new_lines_index;///< Pointer to the array of new-line indices.
    unsigned int *samp_var_id;    ///< Pointer to the array of sample variant IDs.
    unsigned short *samp_id;      ///< Pointer to the array of sample IDs.
    __half *samp_float;           ///< Pointer to the array of sample float values.
    bool *samp_flag;              ///< Pointer to the array of sample flag values.
    int *samp_int;                ///< Pointer to the array of sample integer values.
    char* samp_float_name;        ///< Pointer to the sample float field names.
    char* samp_flag_name;         ///< Pointer to the sample flag field names.
    char* samp_int_name;          ///< Pointer to the sample integer field names.
    int* samp_float_numb;         ///< Pointer to the number of sample float values per entry.
    int* samp_flag_numb;          ///< Pointer to the number of sample flag values per entry.
    int* samp_int_numb;           ///< Pointer to the number of sample integer values per entry.
    char *sample_GT;              ///< Pointer to the sample genotype array.
    int numSample;                ///< Number of samples.
    unsigned int numLines;        ///< Total number of VCF lines (variants).
    int numGT;                    ///< Number of genotype entries.
};

/**
 * @struct info_flag
 * @brief Structure to store flag information from VCF INFO fields.
 *
 * Contains a vector of flag values (as uint8_t) and the name of the flag field.
 */
struct info_flag {
    vector<uint8_t> i_flag;  ///< Vector storing flag values.
    string name;             ///< Name of the flag field.
};

/**
 * @struct info_flag_d
 * @brief Device-side version of info_flag.
 *
 * Contains a pointer to a boolean array for flag values and a pointer to the flag field name.
 */
struct info_flag_d {
    bool* i_flag;  ///< Pointer to device flag array.
    char* name;    ///< Pointer to the flag field name.
};

/**
 * @struct info_string
 * @brief Structure to store string information from VCF INFO fields.
 *
 * Contains a vector of strings and the corresponding field name.
 */
struct info_string {
    vector<string> i_string;  ///< Vector storing string values.
    string name;              ///< Name of the string field.
};

/**
 * @struct info_float
 * @brief Structure to store float information from VCF INFO fields.
 *
 * Contains a vector of half-precision floats and the corresponding field name.
 */
struct info_float {
    vector<__half> i_float;   ///< Vector storing float values (half precision).
    string name;              ///< Name of the float field.
};

/**
 * @struct info_float_d
 * @brief Device-side version of info_float.
 *
 * Contains a pointer to an array of half-precision floats and a pointer to the field name.
 */
struct info_float_d {
    __half *i_float;  ///< Pointer to device float array.
    char *name;       ///< Pointer to the float field name.
};

/**
 * @struct info_int
 * @brief Structure to store integer information from VCF INFO fields.
 *
 * Contains a vector of integers and the corresponding field name.
 */
struct info_int {
    vector<int> i_int;  ///< Vector storing integer values.
    string name;        ///< Name of the integer field.
};

/**
 * @struct info_int_d
 * @brief Device-side version of info_int.
 *
 * Contains a pointer to an array of integers and a pointer to the field name.
 */
struct info_int_d {
    int *i_int;   ///< Pointer to device integer array.
    char *name;   ///< Pointer to the integer field name.
};

/**
 * @struct samp_Flag
 * @brief Structure to store flag information for sample fields.
 *
 * Contains a vector of flag values, the field name, and the number of entries.
 */
struct samp_Flag {
    vector<uint8_t> i_flag;  ///< Vector storing sample flag values.
    string name;             ///< Name of the sample flag field.
    int numb;                ///< Number of entries per sample.
};

/**
 * @struct samp_Flag_d
 * @brief Device-side version of samp_Flag.
 *
 * Contains pointers for the total count, flag values, field name, and number of entries.
 */
struct samp_Flag_d {   
    int *tot;    ///< Pointer to total count.
    bool *i_flag;///< Pointer to device flag array.
    char *name;  ///< Pointer to the sample flag field name.
    int *numb;   ///< Pointer to the number of entries per sample.
};

/**
 * @struct samp_String
 * @brief Structure to store string information for sample fields.
 *
 * Contains a vector of sample strings, the field name, and the number of entries.
 */
struct samp_String {
    vector<string> i_string;  ///< Vector storing sample string values.
    string name;              ///< Name of the sample string field.
    int numb;                 ///< Number of entries per sample.
};

/**
 * @struct samp_Float
 * @brief Structure to store float information for sample fields.
 *
 * Contains a vector of half-precision float values, the field name, and the number of entries.
 */
struct samp_Float {
    vector<__half> i_float;  ///< Vector storing sample float values (half precision).
    string name;             ///< Name of the sample float field.
    int numb;                ///< Number of entries per sample.
};

/**
 * @struct samp_Float_d
 * @brief Device-side version of samp_Float.
 *
 * Contains pointers for total count, float values, field name, and number of entries.
 */
struct samp_Float_d {   
    int *tot;         ///< Pointer to total count.
    __half *i_float;  ///< Pointer to device float array.
    char *name;       ///< Pointer to the sample float field name.
    int *numb;        ///< Pointer to the number of entries per sample.
};

/**
 * @struct samp_Int
 * @brief Structure to store integer information for sample fields.
 *
 * Contains a vector of sample integer values, the field name, and the number of entries.
 */
struct samp_Int {
    vector<int> i_int;  ///< Vector storing sample integer values.
    string name;        ///< Name of the sample integer field.
    int numb;           ///< Number of entries per sample.
};

/**
 * @struct samp_Int_d
 * @brief Device-side version of samp_Int.
 *
 * Contains pointers for total count, integer values, field name, and number of entries.
 */
struct samp_Int_d {
    int *tot;   ///< Pointer to total count.
    int *i_int; ///< Pointer to device integer array.
    char *name; ///< Pointer to the sample integer field name.
    int *numb;  ///< Pointer to the number of entries per sample.
};

/**
 * @struct samp_GT
 * @brief Structure to store genotype information for a sample.
 *
 * Contains a vector of genotype characters and the number of genotype entries.
 */
struct samp_GT {
    vector<char> GT;  ///< Vector storing genotype values.
    int numb;         ///< Number of genotype entries.
};

/**
 * @struct samp_GT_d
 * @brief Device-side version of samp_GT.
 *
 * Contains a pointer to an array of genotype characters and a pointer to the number of entries.
 */
struct samp_GT_d {
    char *GT;   ///< Pointer to device genotype array.
    int *numb;  ///< Pointer to the number of genotype entries.
};

/**
 * @struct header_element
 * @brief Structure representing a VCF header element.
 *
 * Contains vectors for the ID, Number, and Type fields extracted from the header,
 * as well as counters for various types of values.
 */
struct header_element {
    vector<string> ID;       ///< Vector of header IDs.
    vector<string> Number;   ///< Vector of header Number values.
    vector<string> Type;     ///< Vector of header Type values.
    int total_values = 0;    ///< Total number of header values.
    int alt_values = 0;      ///< Number of alternative allele values.
    int no_alt_values = 0;   ///< Number of non-alternative values.
    int ints_alt = 0;        ///< Count of alternative integer values.
    int floats_alt = 0;      ///< Count of alternative float values.
    int strings_alt = 0;     ///< Count of alternative string values.
    int flags_alt = 0;       ///< Count of alternative flag values.
    int ints = 0;            ///< Count of integer values.
    int floats = 0;          ///< Count of float values.
    int strings = 0;         ///< Count of string values.
    int flags = 0;           ///< Count of flag values.
    bool hasGT = false;      ///< Flag indicating whether genotype (GT) data is present.
    char numGT = 0;          ///< Number of genotype entries.
};

#endif