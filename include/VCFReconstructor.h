#ifndef VCF_RECONSTRUCTOR_H
#define VCF_RECONSTRUCTOR_H

#include <string>
#include "VCFDataFrames.h"
#include <cstdint>
#include <cstring>

// Class to rebuild a VCF file from the 4 normalized DataFrames
class VCFReconstructor {
public:
    // Constructor: sets the output file path
    VCFReconstructor(const std::string& output_vcf_path);

    // Main execution method: processes all rows and writes to file
    void run(const var_columns_df& df1, 
             const alt_columns_df& df2, 
             const sample_columns_df& df3, 
             const alt_format_df& df4);

private:
    std::string output_vcf_path; // Path for the reconstructed VCF file

    float decode_fp16(uint16_t h);

    // Formats a single VCF row by joining data from the DataFrames
    std::string formatVariant(int index, 
                              const var_columns_df& df1, 
                              const alt_columns_df& df2, 
                              const sample_columns_df& df3, 
                              const alt_format_df& df4);
};

#endif