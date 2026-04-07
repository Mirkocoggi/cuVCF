#ifndef VCF_RECONSTRUCTOR_H
#define VCF_RECONSTRUCTOR_H

#include <string>
#include "VCFDataFrames.h"
#include <cstdint>
#include <cstring>
#include <map>
#include <vector>

// Class to rebuild a VCF file from the 4 normalized DataFrames
class VCFReconstructor {
public:
    // Constructor: sets the output file path
    VCFReconstructor(const std::string& output_vcf_path,
                     const std::string& header_text,
                     const std::map<std::string, std::string>& field_types,
                     const std::vector<std::string>& samp_names);

    // Main execution method: processes all rows and writes to file
    void run(const var_columns_df& df1, 
             const alt_columns_df& df2, 
             const sample_columns_df& df3, 
             const alt_format_df& df4);

private:
    std::string output_vcf_path; // Path for the reconstructed VCF file
    std::string header_text;
    std::map<std::string, std::string> field_types;
    std::vector<std::string> samp_names;
    std::string format_str;

    size_t df2_cursor;
    size_t df4_cursor;
    bool gt_in_df4;
    bool has_gt;
    size_t df4_start;   

    // Formats a single VCF row by joining data from the DataFrames
    std::string formatVariant(int index, 
                              const var_columns_df& df1, 
                              const alt_columns_df& df2, 
                              const sample_columns_df& df3, 
                              const alt_format_df& df4);
};

#endif