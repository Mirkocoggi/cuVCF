#ifndef GPU_STRUCTS_H
#define GPU_STRUCTS_H

#include <cuda_fp16.h>
#include <cstdint>

/**
 * @brief Structure containing device pointers for encoding dictionaries.
 * 
 * It holds null-terminated concatenated strings and their relative offset arrays 
 * to quickly map raw bytes into text strings within the CUDA kernel.
 */
struct DeviceMaps {
    char* chrom_strings = nullptr;           ///< Buffer for chromosome strings
    unsigned int* chrom_offsets = nullptr;   ///< Starting offsets for each chromosome string
    char* filter_strings = nullptr;          ///< Buffer for filter strings
    unsigned int* filter_offsets = nullptr;  ///< Starting offsets for each filter string
};

/**
 * @brief Shadow structure for var_columns_df (DF1) allocated on the device.
 * 
 * Contains core variant data (CHROM, POS, ID, REF, QUAL, FILTER)
 * and general INFO fields at the single-variant level.
 */
struct DeviceVarColumns {
    unsigned int* var_number = nullptr; ///< Sequential variant ID
    char*         chrom = nullptr;      ///< Raw chromosome ID
    unsigned int* pos = nullptr;        ///< Position (POS)
    char*         id_data = nullptr;    ///< Flat buffer for ID strings
    unsigned int* id_offsets = nullptr; ///< Offset array for id_data
    char*         ref_data = nullptr;   ///< Flat buffer for REF strings
    unsigned int* ref_offsets = nullptr;///< Offset array for ref_data
    __half*       qual = nullptr;       ///< Quality score in float16 (fp16)
    char*         filter = nullptr;     ///< Raw filter ID
    
    // INFO fields
    int*          in_int = nullptr;     ///< Buffer for integer INFO fields
    __half*       in_float = nullptr;   ///< Buffer for float INFO fields (fp16)
    uint8_t*      in_flag = nullptr;    ///< Buffer for flag (boolean) INFO fields
    char*         int_names = nullptr;  ///< Labels for integer fields
    char*         float_names = nullptr;///< Labels for float fields
    char*         flag_names = nullptr; ///< Labels for flag fields
    
    // Dimensions
    int num_int_fields = 0;             ///< Total number of Integer INFO columns
    int num_float_fields = 0;           ///< Total number of Float INFO columns
    int num_flag_fields = 0;            ///< Total number of Flag INFO columns
};

/**
 * @brief Shadow structure for alt_columns_df (DF2) allocated on the device.
 * 
 * Manages alternative alleles and their associated INFO fields. 
 * Supports the handling of multi-allelic variants.
 */
struct DeviceAltColumns {
    unsigned int* var_id = nullptr;     ///< Variant ID this allele belongs to
    char*         alt_data = nullptr;   ///< Flat buffer for ALT strings
    unsigned int* alt_offsets = nullptr;///< Offset array for alt_data
    unsigned int* alt_start = nullptr;  ///< Starting index in alt_data for the i-th variant
    unsigned int* alt_count = nullptr;  ///< Number of alternative alleles for the i-th variant

    int*          alt_int = nullptr;    ///< Buffer for integer INFO-ALT fields
    __half*       alt_float = nullptr;  ///< Buffer for float INFO-ALT fields (fp16)
    char*         alt_int_names = nullptr;   ///< Labels for integer INFO-ALT fields
    char*         alt_float_names = nullptr; ///< Labels for float INFO-ALT fields
    
    // Dimensions
    int num_entries = 0;                ///< Total number of DF2 entries for the current chunk
    int num_alt_int_fields = 0;         ///< Number of allele-specific integer fields
    int num_alt_float_fields = 0;       ///< Number of allele-specific float fields
};

/**
 * @brief Shadow structure for sample data allocated on the device.
 * 
 * Combines data from DF3 and DF4 (pre-formatted on the host side) into a single 
 * flat packed buffer to maximize memory throughput during GPU writes.
 */
struct DeviceSampleColumns {
    int num_samples = 0;                      ///< Total number of samples

    char*         sample_data    = nullptr;   ///< Flat packed buffer with pre-formatted samples
    unsigned int* sample_offsets = nullptr;   ///< Offset array: size [chunk_size * num_samples]

    char* format_str    = nullptr;            ///< Constant FORMAT string (e.g., "GT:AD:DP")
    int   format_str_len = 0;                 ///< Length of the FORMAT string
};

#endif