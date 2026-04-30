#ifndef GPU_STRUCTS_H
#define GPU_STRUCTS_H

#include <cuda_fp16.h>
#include <cstdint>

struct DeviceMaps {
    char* chrom_strings = nullptr;
    unsigned int* chrom_offsets = nullptr;
    char* filter_strings = nullptr;
    unsigned int* filter_offsets = nullptr;
    // Qui aggiungeremo poi tsa, polyphen, csq
};

// Struttura ombra per var_columns_df (DF1)
struct DeviceVarColumns {
    unsigned int* var_number = nullptr;
    char*         chrom = nullptr;
    unsigned int* pos = nullptr;
    char*         id_data = nullptr;
    unsigned int* id_offsets = nullptr;
    char*         ref_data = nullptr;
    unsigned int* ref_offsets = nullptr;
    __half*       qual = nullptr;
    char*         filter = nullptr;
    
    // INFO fields
    int*          in_int = nullptr;
    __half*       in_float = nullptr;
    uint8_t*      in_flag = nullptr;
    char*         int_names = nullptr;
    char*         float_names = nullptr;
    char*         flag_names = nullptr;
    
    // Dimensioni
    int num_int_fields = 0;
    int num_float_fields = 0;
    int num_flag_fields = 0;
};

// Struttura ombra per alt_columns_df (DF2)
struct DeviceAltColumns {
    unsigned int* var_id = nullptr;
    char*         alt_data = nullptr;
    unsigned int* alt_offsets = nullptr;
    unsigned int* alt_start = nullptr;
    unsigned int* alt_count = nullptr;

    int*          alt_int = nullptr;
    __half*       alt_float = nullptr;
    char*         alt_int_names = nullptr;
    char*         alt_float_names = nullptr;
    
    // Dimensioni
    int num_entries = 0;
    int num_alt_int_fields = 0;
    int num_alt_float_fields = 0;
};

#endif