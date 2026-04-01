#include "VCFReconstructor.h"
#include <fstream>
#include <sstream>
#include <cstring>

// Constructor
VCFReconstructor::VCFReconstructor(const std::string& output_vcf_path,
                                   const std::string& header_text,
                                   const std::map<std::string, std::string>& field_types,
                                   const std::vector<std::string>& samp_names
                                  ): 
                                   output_vcf_path(output_vcf_path),
                                   header_text(header_text),
                                   field_types(field_types),
                                   samp_names(samp_names) {}
   
// Decodes a half-precision float (fp16) to a standard single-precision float (fp32)
float VCFReconstructor::decode_fp16(uint16_t h) {
    // Extract sign, exponent, and mantissa bits
    uint16_t sign = (h >> 15) & 0x0001;
    uint16_t exponent = (h >> 10) & 0x001F;
    uint16_t mantissa = h & 0x03FF;

    // Special case: Zero
    if (exponent == 0) return 0.0f;

    // Special case: Infinity or NaN
    if (exponent == 31) return 0.0f;

    // Reconstruct fp32 bit representation
    uint32_t f_bits = (sign << 31) | ((exponent + 127 - 15) << 23) | (mantissa << 13);
    
    float f;
    std::memcpy(&f, &f_bits, sizeof(float));
    return f;
}

// Main execution method: iterates over variants and writes them to the output file
void VCFReconstructor::run(const var_columns_df& df1, const alt_columns_df& df2, const sample_columns_df& df3, const alt_format_df& df4) {
    std::ofstream vcf_file(output_vcf_path);
    if (!vcf_file.is_open()) {
        throw std::runtime_error("Error opening output VCF file: " + output_vcf_path);
    }

    vcf_file << header_text;
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for(size_t i = 0; i < samp_names.size(); i++) {
        vcf_file << "\t" << samp_names[i];
    }
    vcf_file << "\n";

    size_t fsize = df1.var_number.size();

    bool gt_in_df3 = false;
    // Check if GT is present in df3
    for(size_t i = 0; i < df3.samp_string.size(); i++) {
        if((df3.samp_string[i].name == "GT" || df3.samp_string[i].name == "GT0")) {
            gt_in_df3 = true;
            break;
        }
    }

    // Check if GT is present in df4
    bool gt_in_df4 = false;
    if(!gt_in_df3){
        for(size_t i = 0; i < df4.samp_string.size(); i++) {
            if(df4.samp_string[i].name == "GT" || df4.samp_string[i].name == "GT0") {
                gt_in_df4 = true;
                break;
            }
        }
    }

    bool has_gt = gt_in_df3 || gt_in_df4;

    // Iterate through all variants and format each line
    for(size_t i = 0; i < fsize; i++) {
        std::string line = formatVariant(i, df1, df2, df3, df4, gt_in_df4, has_gt);
        vcf_file << line << "\n";
    }
    
    vcf_file.close();
}

// Formats a single variant record into a valid VCF string
std::string VCFReconstructor::formatVariant(int index, const var_columns_df& df1, const alt_columns_df& df2, const sample_columns_df& df3, const alt_format_df& df4, bool gt_in_df4, bool has_gt) {
    std::stringstream ss;
    unsigned int current_var_id = df1.var_number[index];

    // --- 1. CHROM, POS, ID, REF ---
    ss << df1.chrom[index] << "\t";
    ss << df1.pos[index] << "\t";
    ss << df1.id[index] << "\t";
    ss << df1.ref[index] << "\t";

    // --- 2. ALT ---
    // Retrieve all alternative alleles associated with the current variant
    std::string alt_string = "";
    bool first_alt = true;
    std::vector<size_t> alt_indices; // Track DF2 indices for later use in INFO fields

    for(size_t i = 0; i < df2.var_id.size(); i++) {
        if(df2.var_id[i] == current_var_id){
            if(!first_alt) alt_string += ",";
            alt_string += df2.alt[i];
            first_alt = false;
            alt_indices.push_back(i);
        } else if(df2.var_id[i] > current_var_id){
            break; // Stop searching since var_id is sorted
        }
    }

    if(alt_string.empty()) alt_string = ".";
    ss << alt_string << "\t";

    // --- 3. QUAL ---
    if(df1.qual[index] == -1.0f){
        ss << ".\t";
    } else {
        //uint16_t raw_fp16 = static_cast<uint16_t>(df1.qual[index]);
        //ss << decode_fp16(raw_fp16) << "\t";

        ss << df1.qual[index] << "\t";
    }

    // --- 4. FILTER ---
    if(df1.filter[index].empty()){
        ss << ".\t";
    } else {
        ss << df1.filter[index] << "\t";
    }

    // --- 5. INFO ---
    bool first_info = true;

    // Append variant-level INFO fields from DF1
    for(size_t col = 0; col < df1.in_string.size(); col++){
        std::string val = df1.in_string[col].i_string[index];

        if(!val.empty() && val != "."){

            std::string field_name = df1.in_string[col].name;

            // Look up field type in header map, default to String if not found
            std::string type = field_types.count(field_name) ? field_types[field_name] : "String";
            std::string contribution = "";
            if(type == "Flag"){
                if(val == "1") contribution = field_name;
            } else if(type == "Float"){
                try {
                    //uint16_t raw_fp16 = static_cast<uint16_t>(std::stoul(val));
                    //contribution = field_name + "=" + std::to_string(decode_fp16(raw_fp16));

                    contribution = field_name + "=" + val;
                } catch(const std::exception&) {
                    contribution = field_name + "=" + val;
                }
            } else {
                contribution = field_name + "=" + val;
            }

            if(!contribution.empty()){
                if(!first_info) ss << ";";
                ss << contribution;
                first_info = false;
            }
        }
    }

    // Append allele-level INFO fields from DF2 (comma-separated for multiple alleles)
    for(size_t col = 0; col < df2.alt_string.size(); col++){
        std::string val_str = "";
        
        std::string field_name = df2.alt_string[col].name;

        std::string type = field_types.count(field_name) ? field_types[field_name] : "String";

        bool has_data = false;

        for(size_t idx = 0; idx < alt_indices.size(); idx++){

            std::string alt_val = df2.alt_string[col].i_string[alt_indices[idx]];
            if(alt_val.empty()) alt_val = ".";

            if(type == "Float" && alt_val != "."){
                try {
                    //alt_val = std::to_string(decode_fp16(static_cast<uint16_t>(std::stoul(alt_val))));
                } catch(const std::exception&) {}
            } 

           if(idx > 0) val_str += ",";
           val_str += alt_val;
           
           if(alt_val != ".") has_data = true;
        }



        if(has_data){
            if(!first_info) ss << ";";
            ss << field_name << "=" << val_str;
            first_info = false;
        }
    }
    
    if(first_info) ss << "."; // Empty INFO block indicator
    ss << "\t";

    // --- 6. FORMAT ---
    // Construct the FORMAT template string (e.g., GT:AD:DP)
    std::string format_str = "";
    if(has_gt) format_str = "GT";  

    // Append DF3 format keys
    for(size_t col = 0; col < df3.samp_string.size(); col++){
        if(df3.samp_string[col].name == "GT" || df3.samp_string[col].name == "GT0") continue;
        if(!format_str.empty()) format_str += ":";
        format_str += df3.samp_string[col].name;
    }    

    // Append DF4 format keys
    for(size_t col = 0; col < df4.samp_string.size(); col++){
        if(df4.samp_string[col].name == "GT" || df4.samp_string[col].name == "GT0") continue;
        if(!format_str.empty()) format_str += ":";
        format_str += df4.samp_string[col].name;   
    }

    ss << format_str << "\t";

    // --- 7. SAMPLES ---
    // Find starting index for current variant in DF4
    int start_df4 = -1;
    for(size_t i = 0; i < df4.var_id.size(); i++){
        if(df4.var_id[i] == current_var_id){
            start_df4 = i;
            break;
        }
    }

    // Iterate through all samples to build their specific data strings
    for(int i = 0; i < df3.numSample; i++){
        size_t df3_idx = (static_cast<size_t>(index) * df3.numSample) + i;

        // Initialize GT from DF3 or empty if GT lives in DF4
        std::string sample_data = "";
        if(has_gt){
            if(gt_in_df4){
                sample_data = "";
            } else {
                for(size_t col = 0; col < df3.samp_string.size(); col++){
                    if(df3.samp_string[col].name == "GT" || df3.samp_string[col].name == "GT0"){
                        sample_data = df3.samp_string[col].i_string[df3_idx];
                        break;
                    }
                }
            }
        }

        // Gather DF4 data in a single pass: GT and all other per-allele fields
        std::vector<std::string> df4_fields(df4.samp_string.size(), "");

        if(start_df4 != -1){
            size_t j = static_cast<size_t>(start_df4);
            while(j < df4.var_id.size() && df4.var_id[j] == current_var_id){
                if(df4.samp_id[j] == static_cast<unsigned short>(i)){
                    for(size_t col = 0; col < df4.samp_string.size(); col++){
                        if(df4.samp_string[col].name == "GT" || df4.samp_string[col].name == "GT0"){
                            // GT is the same for all alleles, take first occurrence
                            if(gt_in_df4 && sample_data.empty())
                                sample_data = df4.samp_string[col].i_string[j];
                        } else {
                            if(!df4_fields[col].empty()) df4_fields[col] += ",";
                            df4_fields[col] += df4.samp_string[col].i_string[j];
                        }
                    }
                }
                j++;
            }
        }

        // Append sample-level data from DF3
        for(size_t col = 0; col < df3.samp_string.size(); col++){
            if(df3.samp_string[col].name == "GT" || df3.samp_string[col].name == "GT0") continue;
            if(!sample_data.empty()) sample_data += ":";
            sample_data += df3.samp_string[col].i_string[df3_idx];
        }

        // Append per-allele data from DF4
        for(size_t col = 0; col < df4_fields.size(); col++){
            if(df4.samp_string[col].name == "GT" || df4.samp_string[col].name == "GT0") continue;
            sample_data += ":";
            sample_data += (df4_fields[col].empty() ? "." : df4_fields[col]);
        }

        ss << sample_data;
        if(i < df3.numSample - 1) ss << "\t";
    }

    return ss.str();
}