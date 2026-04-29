#ifndef VCF_RECONSTRUCTOR_GPU_H
#define VCF_RECONSTRUCTOR_GPU_H

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include "VCFDataFrames.h"
#include <cuda_runtime.h>
#include <cuda_fp16.h>

#define MAX_LINE_LEN 1024
#define MAX_NAME_LEN 16
#define CHUNK_SIZE   100000

class VCFReconstructorGPU {
public:
    VCFReconstructorGPU(const std::string& output_vcf_path,
                        const std::string& header_text);
    ~VCFReconstructorGPU();

    void run(const var_columns_df& df1,
             const alt_columns_df& df2);

private:
    std::string output_vcf_path;
    std::string header_text;

    // Inverse maps (host side)
    std::map<char, std::string> inv_chrom_map;
    std::map<char, std::string> inv_filter_map;
    std::map<int,  std::string> inv_tsa_map;
    std::map<char, std::string> inv_polyphen_map;
    std::map<char, std::string> inv_csq_map;

    // Device output buffer
    char*         d_output;       // [CHUNK_SIZE * MAX_LINE_LEN]
    unsigned int* d_line_lens;    // lunghezza di ogni riga

    // Device DF1
    unsigned int* d_var_number;
    char*         d_chrom;
    unsigned int* d_pos;
    char*         d_id_data;
    unsigned int* d_id_offsets;
    char*         d_ref_data;
    unsigned int* d_ref_offsets;
    __half*       d_qual;
    char*         d_filter;
    int*          d_in_int;       // [num_int_fields * chunk_size]
    __half*       d_in_float;     // [num_float_fields * chunk_size]
    uint8_t*      d_in_flag;      // [num_flag_fields * chunk_size]
    char*         d_int_names;    // [num_int_fields * MAX_NAME_LEN]
    char*         d_float_names;
    char*         d_flag_names;

    // Device DF2
    unsigned int* d_alt_var_id;
    char*         d_alt_data;
    unsigned int* d_alt_offsets;
    int*          d_alt_int;
    __half*       d_alt_float;
    char*         d_alt_int_names;
    char*         d_alt_float_names;

    // Device mappe inverse
    char*         d_chrom_strings;
    unsigned int* d_chrom_offsets;
    char*         d_filter_strings;
    unsigned int* d_filter_offsets;

    // Dimensioni (note dopo il primo chunk)
    int num_int_fields;
    int num_float_fields;
    int num_flag_fields;
    int num_alt_int_fields;
    int num_alt_float_fields;

    void buildInverseMaps(const var_columns_df& df1);

    void allocateDevice(const var_columns_df& df1,
                        const alt_columns_df& df2,
                        int chunk_size,
                        int chunk_start,
                        int chunk_end,
                        int df2_start);

    void freeDevice();

    void hostToDevice(const var_columns_df& df1,
                      const alt_columns_df& df2,
                      int chunk_start,
                      int chunk_end,
                      int df2_start,
                      int df2_end);

    void writeChunk(int num_variants, std::ofstream& out);
    
};

#endif