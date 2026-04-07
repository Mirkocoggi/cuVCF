#include "VCFReconstructor.h"
#include "Utils.h"
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
   

// Main execution method: iterates over variants and writes them to the output file
void VCFReconstructor::run(const var_columns_df& df1, const alt_columns_df& df2, const sample_columns_df& df3, const alt_format_df& df4) {
    std::ofstream vcf_file(output_vcf_path);
    if (!vcf_file.is_open()) {
        throw std::runtime_error("Error opening output VCF file: " + output_vcf_path);
    }

    vcf_file << header_text;
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    if(df3.numSample > 0){
        vcf_file << "\tFORMAT";
        for(size_t i = 0; i < samp_names.size(); i++) {
            vcf_file << "\t" << samp_names[i];
        }
    } 
    vcf_file << "\n";

    size_t fsize = df1.var_number.size();

    gt_in_df4 = false;
    bool gt_in_df3 = false;
    if(df3.numSample > 0){
        gt_in_df3 = !df3.samp_string.empty() && isGTField(df3.samp_string[0].name);
        if(!gt_in_df3){
            gt_in_df4 = !df4.samp_string.empty() && isGTField(df4.samp_string[0].name);
        }
    }

    has_gt = gt_in_df3 || gt_in_df4;
    df4_start = (!df4.samp_string.empty() && isGTField(df4.samp_string[0].name)) ? 1 : 0;

    // Build FORMAT string once
    format_str = "";
    if(df3.numSample > 0){
        if(gt_in_df4 && !df4.samp_string.empty() && isGTField(df4.samp_string[0].name)){
            format_str += df4.samp_string[0].name;
        }
        for(size_t col = 0; col < df3.samp_string.size(); col++){
            if(!format_str.empty()) format_str += ":";
            format_str += df3.samp_string[col].name;
        }
        for(size_t col = df4_start; col < df4.samp_string.size(); col++){
            if(!format_str.empty()) format_str += ":";
            format_str += df4.samp_string[col].name;
        }
    }
    
    // Iterate through all variants and format each line
    df2_cursor = 0;
    df4_cursor = 0;
    for(size_t i = 0; i < fsize; i++) {
        std::string line = formatVariant(i, df1, df2, df3, df4);
        vcf_file << line << "\n";
    }
    
    vcf_file.close();
}

// Formats a single variant record into a valid VCF string
std::string VCFReconstructor::formatVariant(int index, const var_columns_df& df1, const alt_columns_df& df2, const sample_columns_df& df3, const alt_format_df& df4) {
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

    while(df2_cursor < df2.var_id.size() && df2.var_id[df2_cursor] < current_var_id){
        df2_cursor++;
    }

    size_t j = df2_cursor ;
    while(j < df2.var_id.size() && df2.var_id[j] == current_var_id){
        if(!first_alt) alt_string += ",";
        alt_string += df2.alt[j];
        first_alt = false;
        alt_indices.push_back(j);
        j++;
    }

    if(alt_string.empty()) alt_string = ".";
    ss << alt_string << "\t";

    // --- 3. QUAL ---
    if(df1.qual[index] == -1.0f){
        ss << ".\t";
    } else {
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

    if(df3.numSample > 0){
        // --- 6. FORMAT ---
        ss << format_str << "\t";

        // --- 7. SAMPLES ---
        // Advance DF4 cursor to current variant
        while(df4_cursor < df4.var_id.size() && df4.var_id[df4_cursor] < current_var_id){
            df4_cursor++;
        }
        int start_df4 = (df4_cursor < df4.var_id.size() && df4.var_id[df4_cursor] < current_var_id) 
                        ? static_cast<int>(df4_cursor) : -1;
        
            

        // Iterate through all samples to build their specific data strings
        size_t df3_col_start = (has_gt && !gt_in_df4) ? 1 : 0;
        for(int i = 0; i < df3.numSample; i++){

            if(i>0) ss<<"\t";

            size_t df3_idx = (static_cast<size_t>(index) * df3.numSample) + i;
            std::string sample_data = "";

            // GT from DF4 (first, if applicable)
            if(gt_in_df4 && start_df4 != -1){
                size_t j = static_cast<size_t>(start_df4);
                while(j < df4.var_id.size() && df4.var_id[j] == current_var_id){
                    if(df4.samp_id[j] == static_cast<unsigned short>(i)){
                        sample_data = df4.samp_string[0].i_string[j];
                        break;
                    }
                    j++;
                }
            }

            // GT from DF3 (position 0)
            if(!gt_in_df4 && has_gt){
                sample_data = df3.samp_string[0].i_string[df3_idx];
            }

            if(sample_data == ".|." || sample_data == "./."){
                ss << sample_data;
                continue;
            }

            // Remaining DF3 fields (skip position 0 if it was GT)
            for(size_t col = df3_col_start; col < df3.samp_string.size(); col++){
                if(!sample_data.empty()) sample_data += ":";
                sample_data += df3.samp_string[col].i_string[df3_idx];
            }

            // DF4 per-allele fields (skip index 0 if it was GT)
            std::vector<std::string> df4_fields(df4.samp_string.size(), "");
            if(start_df4 != -1){
                size_t j = static_cast<size_t>(start_df4);
                while(j < df4.var_id.size() && df4.var_id[j] == current_var_id){
                    if(df4.samp_id[j] == static_cast<unsigned short>(i)){
                        for(size_t col = 0; col < df4.samp_string.size(); col++){
                            if(!df4_fields[col].empty()) df4_fields[col] += ",";
                            df4_fields[col] += df4.samp_string[col].i_string[j];
                        }
                    }
                    j++;
                }
            }

            for(size_t col = df4_start; col < df4_fields.size(); col++){
                sample_data += ":";
                sample_data += (df4_fields[col].empty() ? "." : df4_fields[col]);
            }

            ss << sample_data;
        }
    }


    return ss.str();
}