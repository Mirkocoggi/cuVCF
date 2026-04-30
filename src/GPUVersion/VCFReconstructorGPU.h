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

    DeviceVarColumns d_df1;
    DeviceAltColumns d_df2;
    DeviceMaps d_maps;

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
                      int chunk_size,
                      int chunk_start,
                      int chunk_end,
                      int df2_start);

    void writeChunk(int num_variants, std::ofstream& out);
    
};

#endif