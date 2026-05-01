#ifndef VCF_RECONSTRUCTOR_GPU_H
#define VCF_RECONSTRUCTOR_GPU_H

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include "VCFDataFrames.h"
#include "GPUStructs.h"
#include <cuda_runtime.h>
#include <cuda_fp16.h>

#define MAX_LINE_LEN 1024
#define MAX_NAME_LEN 16
#define MAX_SAMPLE_STRING_LEN 256
#define CHUNK_SIZE   10000

struct HostBuffers {
    // DF1
    std::vector<char>          id_buffer;
    std::vector<unsigned int>  id_offsets;
    std::vector<char>          ref_buffer;
    std::vector<unsigned int>  ref_offsets;
    std::vector<int>           in_int_buffer;
    std::vector<__half>        in_float_buffer;
    std::vector<uint8_t>       in_flag_buffer;

    // DF2
    int                        df2_count = 0;
    std::vector<char>          alt_data_buffer;
    std::vector<unsigned int>  alt_data_offsets;
    std::vector<unsigned int>  alt_start_buf;
    std::vector<unsigned int>  alt_count_buf;
    std::vector<int>           alt_int_buffer;
    std::vector<__half>        alt_float_buffer;

    // DF3
    std::vector<char>          sample_buffer;
    std::vector<unsigned int>  sample_offsets;

    // Metadati
    int chunk_size  = 0;
    int chunk_start = 0;
    int chunk_end   = 0;
    int df2_start   = 0;
};

class VCFReconstructorGPU {
public:
    VCFReconstructorGPU(const std::string& output_vcf_path,
                        const std::string& header_text);
    ~VCFReconstructorGPU();

    void run(const var_columns_df& df1,
             const alt_columns_df& df2,
             const sample_columns_df& df3,
             const alt_format_df& df4);

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

    DeviceVarColumns d_df1;
    DeviceAltColumns d_df2;
    DeviceSampleColumns d_df3;
    DeviceMaps d_maps;

    bool has_gt;
    bool gt_in_df3;
    bool gt_in_df4;
    std::string format_str;
    std::vector<std::string> ordered_samp_names;

    std::map<std::string, std::string> parseFormatNumbers(const std::string& header_text);
    std::map<std::string, std::string> format_numbers;

    HostBuffers host_buffers;

    void buildInverseMaps(const var_columns_df& df1);

    void allocateDevice(const var_columns_df& df1,
                        const alt_columns_df& df2,
                        const sample_columns_df& df3,
                        int chunk_size,
                        int chunk_start,
                        int chunk_end,
                        int df2_start);

    void freeDevice();

    void prepareHostBuffers(const var_columns_df& df1,
                            const alt_columns_df& df2,
                            const sample_columns_df& df3,
                            const alt_format_df& df4,
                            HostBuffers& buffers);

    void uploadToDevice(const var_columns_df& df1,
                        const alt_columns_df& df2,
                        const HostBuffers& buffers);

    void writeChunk(int num_variants, std::ofstream& out);

    void buildSampleNames(const sample_columns_df& df3);

    void buildSampleStrings(const var_columns_df& df1,
                        const alt_columns_df& df2,
                        const sample_columns_df& df3,
                        const alt_format_df& df4,
                        int chunk_size,
                        int chunk_start,
                        int chunk_end,
                        int df2_start,
                        std::vector<char>& buffer,
                        std::vector<unsigned int>& offsets);
    
};

#endif