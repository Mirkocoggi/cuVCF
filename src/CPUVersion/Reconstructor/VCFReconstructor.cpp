#include "VCFReconstructor.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <stdexcept>
#include <omp.h>
#include <cmath>
#include <algorithm>
#include <numeric>


VCFReconstructor::VCFReconstructor(const std::string& output_vcf_path,
                                   const std::string& header_text)
    : output_vcf_path(output_vcf_path),
      header_text(header_text) {
}

/*Base ID of a FORMAT column. A declared ID like PS1 stays as it is. A split column (AD0, AD1, ...)
  loses one digit at a time until it matches an ID declared with a fixed Number greater than 1,
  so AF20 gives AF2 and not AF. Without a match all trailing digits are removed*/
static std::string formatBase(const std::string& name, const std::map<std::string, std::string>& numbers){
    if (numbers.count(name)) return name;
    std::string r = name;
    while (!r.empty() && std::isdigit(static_cast<unsigned char>(r.back()))) {
        r.pop_back();
        auto it = numbers.find(r);
        if (it != numbers.end() && std::atoi(it->second.c_str()) > 1) return r;
    }
    return r;
}

/*
 * Recovers the original string values for CHROM and FILTER.
 * In the DataFrames, these are stored as compact codes to save RAM.
 * We invert the maps (code -> string) so we can quickly print the actual text to the VCF.
 */
void VCFReconstructor::buildInverseMaps(const var_columns_df& df1) {
    for (const auto& pair : df1.chrom_map) {
        inv_chrom_map[pair.second] = pair.first;
    }
    for (const auto& pair : df1.filter_map) {
        inv_filter_map[pair.second] = pair.first;
    }
}

/*
 * Builds a 0-indexed vector of sample names. 
 * This ensures the columns in the VCF are printed in the exact order of their internal samp_id.
 */
void VCFReconstructor::buildSampleNames(const sample_columns_df& df3) {
    ordered_samp_names.resize(df3.sampNames.size());
    for (const auto& pair : df3.sampNames) {
        ordered_samp_names[pair.second] = pair.first;
    }
}

std::map<std::string, std::string> VCFReconstructor::parseFormatNumbers(const std::string& header_text) {
    std::map<std::string, std::string> result;
    std::istringstream stream(header_text);
    std::string line;
    
    while (std::getline(stream, line)) {
        if (line.substr(0, 8) != "##FORMAT") continue;
        
        // Extract ID
        size_t id_pos = line.find("ID=");
        size_t num_pos = line.find("Number=");
        if (id_pos == std::string::npos || num_pos == std::string::npos) continue;
        
        size_t id_start = id_pos + 3;
        size_t id_end = line.find_first_of(",", id_start);
        std::string id = line.substr(id_start, id_end - id_start);
        
        size_t num_start = num_pos + 7;
        size_t num_end = line.find_first_of(",", num_start);
        std::string number = line.substr(num_start, num_end - num_start);
        
        result[id] = number;
    }
    return result;
}

std::vector<VCFReconstructor::ChunkCursor> VCFReconstructor::computeChunkCursor(
    const var_columns_df& df1,
    const alt_columns_df& df2,
    const alt_format_df& df4,
    const std::vector<unsigned int>& rows,
    int num_threads
) {
    size_t total_vars = rows.size();
    size_t chunk_size = (total_vars + num_threads - 1) / num_threads;

    std::vector<ChunkCursor> chunks(num_threads);

    for(int c = 0; c < num_threads; c++){
        size_t var_start = std::min(c * chunk_size, total_vars);
        size_t var_end = std::min(var_start + chunk_size, total_vars);

        chunks[c].var_start = var_start;
        chunks[c].var_end = var_end;
        chunks[c].df2_start = 0;
        chunks[c].df4_start = 0;

        if(var_start == var_end) continue;       //empty chunk: nothing to write

        // DF2 and DF4 are sorted by var_id: start at the rows of the chunk's first variant.
        // formatVariant() then skips forward over the rows of the variants not in the selection.
        unsigned int first_var_id = df1.var_number[rows[var_start]];
        chunks[c].df2_start = std::lower_bound(df2.var_id.begin(), df2.var_id.end(), first_var_id) - df2.var_id.begin();
        chunks[c].df4_start = std::lower_bound(df4.var_id.begin(), df4.var_id.end(), first_var_id) - df4.var_id.begin();
    }
    return chunks;
}

void VCFReconstructor::run(const var_columns_df& df1,
                            const alt_columns_df& df2,
                            const sample_columns_df& df3,
                            const alt_format_df& df4,
                            int num_threads)
{
    // Every variant, in file order
    std::vector<unsigned int> rows(df1.var_number.size());
    std::iota(rows.begin(), rows.end(), 0u);
    run(df1, df2, df3, df4, num_threads, rows);
}

void VCFReconstructor::run(const var_columns_df& df1,
                            const alt_columns_df& df2,
                            const sample_columns_df& df3,
                            const alt_format_df& df4,
                            int num_threads,
                            const std::vector<unsigned int>& rows)
{

    if(num_threads < 1) num_threads = 1;    

    std::ofstream vcf_file(output_vcf_path);
    if (!vcf_file.is_open()) {
        throw std::runtime_error("Error opening output VCF file: " + output_vcf_path);
    }

    // 1. Setup lookup tables
    buildInverseMaps(df1);
    buildSampleNames(df3);
    format_numbers = parseFormatNumbers(header_text);

    // 2. Write the static header and standard columns
    vcf_file << header_text;
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    // Append sample names to the header if they exist
    if (df3.numSample > 0) {
        vcf_file << "\tFORMAT";
        for (const auto& name : ordered_samp_names) {
            vcf_file << "\t" << name;
        }
    }
    vcf_file << "\n";

    // 3. Determine where the Genotype (GT) data lives
    // It can be in DF3 or DF4
    gt_in_df3 = false;
    gt_in_df4 = false;
    if (df3.numSample > 0) {
        gt_in_df3 = !df3.sample_GT.empty();
        if (!gt_in_df3) {
            gt_in_df4 = !df4.sample_GT.GT.empty();
        }
    }
    has_gt = gt_in_df3 || gt_in_df4;

    // 4. Pre-build the FORMAT string (e.g., "GT:AD:DP:GQ")
    // We only need to build this once because cuVCF enforces a uniform FORMAT across all variants  
    format_str = "";
    if (df3.numSample > 0) {
        if (has_gt) format_str = "GT";

        // Group DF3 integer fields
        // E.g., if we see AD0, AD1, we only append "AD" to the format string once
        std::string prev_base = "";
        for (size_t i = 0; i < df3.samp_int.size(); i++) {
            std::string base = formatBase(df3.samp_int[i].name, format_numbers);
            if (base != prev_base) {
                if (!format_str.empty()) format_str += ":";
                format_str += base;
                prev_base = base;
            }
        }

        // Group DF3 float fields
        prev_base = "";
        for (size_t i = 0; i < df3.samp_float.size(); i++) {
            std::string base = formatBase(df3.samp_float[i].name, format_numbers);
            if (base != prev_base) {
                if (!format_str.empty()) format_str += ":";
                format_str += base;
                prev_base = base;
            }
        }

        // DF3 string fields are usually singular (no grouping needed)
        for (size_t i = 0; i < df3.samp_string.size(); i++) {
            if (!format_str.empty()) format_str += ":";
            format_str += df3.samp_string[i].name;
        }

        // Group DF4 integer fields (per-allele format fields)
        prev_base = "";
        for (size_t i = 0; i < df4.samp_int.size(); i++) {
            std::string base = formatBase(df4.samp_int[i].name, format_numbers);
            if (base != prev_base) {
                if (!format_str.empty()) format_str += ":";
                format_str += base;
                prev_base = base;
            }
        }

        // Group DF4 float fields
        prev_base = "";
        for (size_t i = 0; i < df4.samp_float.size(); i++) {
            std::string base = formatBase(df4.samp_float[i].name, format_numbers);
            if (base != prev_base) {
                if (!format_str.empty()) format_str += ":";
                format_str += base;
                prev_base = base;
            }
        }

        // DF4 string fields
        for (size_t i = 0; i < df4.samp_string.size(); i++) {
            if (!format_str.empty()) format_str += ":";
            format_str += df4.samp_string[i].name;
        }
    }

    // 5. Pre-compute chunk cursors 
    std::vector<ChunkCursor> chunks = computeChunkCursor(df1, df2, df4, rows, num_threads);

    // 6. Parallel processing - each thread writes to its own temp file.
    // An exception must not leave the OpenMP region (std::terminate): keep the first one, rethrow below.
    std::string error;
    #pragma omp parallel for num_threads(num_threads) schedule(static)
    for (int c = 0; c < num_threads; c++) {
        try {
            std::string tmp_path = output_vcf_path + ".tmp_" + std::to_string(c);
            std::ofstream tmp_file(tmp_path);
            if (!tmp_file.is_open())
                throw std::runtime_error("Error opening temporary file: " + tmp_path);

            // Each thread has its own local cursors — no shared state
            size_t local_df2_cursor = chunks[c].df2_start;
            size_t local_df4_cursor = chunks[c].df4_start;

            for(size_t p = chunks[c].var_start; p < chunks[c].var_end; p++){
                tmp_file << formatVariant(rows[p], df1, df2, df3, df4,
                                          local_df2_cursor,
                                          local_df4_cursor) << "\n";
            }
            tmp_file.close();
            if (!tmp_file)
                throw std::runtime_error("Error writing temporary file: " + tmp_path);
        } catch (const std::exception& e) {
            #pragma omp critical(reconstructor_error)
            if (error.empty()) error = e.what();
        }
    }

    // 7. Concatenate temp files in order. An empty chunk leaves an empty file: skip it, inserting
    // the rdbuf() of an empty file sets failbit on vcf_file.
    for(int c = 0; c < num_threads && error.empty(); c++){
        std::string tmp_path = output_vcf_path + ".tmp_" + std::to_string(c);
        std::ifstream tmp_file(tmp_path);
        if (!tmp_file.is_open()) error = "Error reading temporary file: " + tmp_path;
        else if (tmp_file.peek() != std::ifstream::traits_type::eof()) vcf_file << tmp_file.rdbuf();
    }
    for(int c = 0; c < num_threads; c++)
        std::remove((output_vcf_path + ".tmp_" + std::to_string(c)).c_str());

    vcf_file.close();
    if (error.empty() && !vcf_file)
        error = "Error writing output VCF file (disk full?): " + output_vcf_path;
    if (!error.empty())
        throw std::runtime_error(error + " - the output VCF is incomplete");
} 

std::string VCFReconstructor::formatVariant(int index,
                                             const var_columns_df& df1,
                                             const alt_columns_df& df2,
                                             const sample_columns_df& df3,
                                             const alt_format_df& df4,
                                             size_t& df2_cursor,
                                             size_t& df4_cursor) {
    std::stringstream ss;
    unsigned int current_var_id = df1.var_number[index];

    // --- 1. CHROM, POS, ID, REF ---
    // Look up the char code in the inverse map. Default to "." if not found
    auto chrom_it = inv_chrom_map.find(df1.chrom[index]);
    ss << (chrom_it != inv_chrom_map.end() ? chrom_it->second : ".") << "\t";
    ss << df1.pos[index] << "\t";
    ss << df1.id[index] << "\t";
    ss << df1.ref[index] << "\t";

    
    // --- 2. ALT ---
    // Fetch all alternate alleles for this variant from DF2.
    std::string alt_string = "";
    bool first_alt = true;
    std::vector<size_t> alt_indices;

    // Advance the cursor until we reach the current variant ID.
    // This ignores old IDs without doing a full array scan.
    while (df2_cursor < df2.var_id.size() && df2.var_id[df2_cursor] < current_var_id) {
        df2_cursor++;
    }

    // Collect all alleles matching this variant ID
    size_t j = df2_cursor;
    while (j < df2.var_id.size() && df2.var_id[j] == current_var_id) {
        if (!first_alt) alt_string += ",";
        alt_string += df2.alt[j];
        first_alt = false;

        // Save the index so we can quickly grab per-allele INFO fields later
        alt_indices.push_back(j);
        j++;
    }
    if (alt_string.empty()) alt_string = ".";
    ss << alt_string << "\t";


    // --- 3. QUAL ---
    float qual_val = (float)df1.qual[index];
    if (qual_val == -1.0f) {
        ss << ".\t";
    } else {
        ss << qual_val << "\t";
    }

    // --- 4. FILTER ---
    auto filter_it = inv_filter_map.find(df1.filter[index]);
    if (filter_it != inv_filter_map.end()) {
        ss << filter_it->second << "\t";
    } else {
        ss << ".\t";
    }

    // --- 5. INFO ---
    bool first_info = true;

    // 5a-5d: Iterate through standard INFO fields from DF1.
    // - Flags are appended as a simple KEY if true (e.g., "SOMATIC").
    // - Integers, Floats, and Strings are appended as KEY=VALUE if present.

    // 5a. Flags
    for (size_t col = 0; col < df1.in_flag.size(); col++) {
        if (df1.in_flag[col].i_flag[index]) {
            if (!first_info) ss << ";";
            ss << df1.in_flag[col].name;
            first_info = false;
        }
    }

    // 5b. Integers
    for (size_t col = 0; col < df1.in_int.size(); col++) {
        int val = df1.in_int[col].i_int[index];
        if (val != MISSING_INT && !(df1.in_int[col].name == "TSA" && val == TSA_UNKNOWN)) {
            if (!first_info) ss << ";";
            const std::string& name = df1.in_int[col].name;
            if(name == "TSA" && val >= 0 && val < static_cast<int>(std::size(tsa_names)))
                ss << name << "=" << tsa_names[val];
            else
                ss << name << "=" << val;
            first_info = false;
        }
    }

    // 5c. Floats
    for (size_t col = 0; col < df1.in_float.size(); col++) {
        float val = (float)df1.in_float[col].i_float[index];
        if (!std::isnan(val)) {
            if (!first_info) ss << ";";
            ss << df1.in_float[col].name << "=" << val;
            first_info = false;
        }
    }

    // 5d. Strings
    for (size_t col = 0; col < df1.in_string.size(); col++) {
        const std::string& val = df1.in_string[col].i_string[index];
        if (!val.empty() && val != "." && val != "\0") {
            if (!first_info) ss << ";";
            ss << df1.in_string[col].name << "=" << val;
            first_info = false;
        }
    }

    // 5e. Per-allele floats from DF2
    // Fetch float data that is specific to each alternate allele (using alt_indices).
    // Multiple values are joined by commas (e.g., "AF=0.5,0.2"). Missing values become ".".
    for (size_t col = 0; col < df2.alt_float.size(); col++) {
        std::string val_str = "";
        bool has_data = false;
        for (size_t idx = 0; idx < alt_indices.size(); idx++) {
            float v = (float)df2.alt_float[col].i_float[alt_indices[idx]];
            if (idx > 0) val_str += ",";
            if (!std::isnan(v)) {
                val_str += std::to_string(v);
                has_data = true;
            } else {
                val_str += ".";
            }
        }
        if (has_data) {
            if (!first_info) ss << ";";
            ss << df2.alt_float[col].name << "=" << val_str;
            first_info = false;
        }
    }

    // 5f. Per-allele ints from DF2
    // Fetch integer data specific to each alternate allele. Join with commas.
    for (size_t col = 0; col < df2.alt_int.size(); col++) {
        std::string val_str = "";
        bool has_data = false;
        for (size_t idx = 0; idx < alt_indices.size(); idx++) {
            int v = df2.alt_int[col].i_int[alt_indices[idx]];
            if (idx > 0) val_str += ",";
            if (v != MISSING_INT) {
                val_str += std::to_string(v);
                has_data = true;
            } else {
                val_str += ".";
            }
        }
        if (has_data) {
            if (!first_info) ss << ";";
            ss << df2.alt_int[col].name << "=" << val_str;
            first_info = false;
        }
    }

    // 5g. Per-allele strings from DF2
    // Fetch string data specific to each alternate allele. Join with commas.
    for (size_t col = 0; col < df2.alt_string.size(); col++) {
        std::string val_str = "";
        bool has_data = false;
        for (size_t idx = 0; idx < alt_indices.size(); idx++) {
            const std::string& v = df2.alt_string[col].i_string[alt_indices[idx]];
            if (idx > 0) val_str += ",";
            if (!v.empty() && v != "\0") {
                val_str += v;
                has_data = true;
            } else {
                val_str += ".";
            }
        }
        if (has_data) {
            if (!first_info) ss << ";";
            ss << df2.alt_string[col].name << "=" << val_str;
            first_info = false;
        }
    }

    // If no INFO fields were present for this variant, write a dot as required by VCF specs.
    if (first_info) ss << ".";

    // --- 6-7. FORMAT + SAMPLES ---
    // Process sample columns only if samples exist in the dataset.
    if (df3.numSample > 0) {
        ss << "\t" << format_str << "\t";

        // Advance DF4 cursor
        // Synchronize df4_cursor with the current variant ID
        while (df4_cursor < df4.var_id.size() && df4.var_id[df4_cursor] < current_var_id) {
            df4_cursor++;
        }

        // Save the starting index of DF4 data for this variant (if it exists).
        int start_df4 = (df4_cursor < df4.var_id.size() && df4.var_id[df4_cursor] == current_var_id)
                        ? static_cast<int>(df4_cursor) : -1;

        // Iterate horizontally over every sample
        for (int i = 0; i < df3.numSample; i++) {
            if (i > 0) ss << "\t";

            // Calculate the absolute 1D index for DF3 (which is flattened: variant * numSamples + sample_id)
            size_t df3_idx = (static_cast<size_t>(index) * df3.numSample) + i;
            std::string sample_data = "";

            // --- GT ---
            // Retrieve Genotype string. Fallback to "./." if the genotype is not stored.
            if (has_gt) {
                bool found = false;
                char gt_code = 0;
                if(gt_in_df3) {
                    gt_code = df3.sample_GT[0].GT[df3_idx];
                    found = true;
                } else if (gt_in_df4&& start_df4 != -1) {
                    size_t k = static_cast<size_t>(start_df4);

                    //Search for the specific sample_id within the current variant's DF4 chunk
                    while(k < df4.var_id.size() && df4.var_id[k] == current_var_id) {
                        if(df4.samp_id[k] == static_cast<unsigned short>(i)){
                            gt_code = df4.sample_GT.GT[k];
                            found = true;
                            break; 
                        }
                        k++;
                    }
                }

                sample_data = found ? df3.getGTStringFromChar(gt_code) : "./.";
            }

            // --- DF3 int fields ---
            // Re-group array fields (e.g., AD0, AD1) separated by commas.
            {
                size_t col = 0;
                while (col < df3.samp_int.size()) {
                    std::string base = formatBase(df3.samp_int[col].name, format_numbers);
                    if (!sample_data.empty()) sample_data += ":";

                    // Each column is one position of the array: a missing value keeps its place
                    std::string group_val = "";
                    bool any = false;
                    while (col < df3.samp_int.size() &&
                        formatBase(df3.samp_int[col].name, format_numbers) == base) {
                        if (!group_val.empty()) group_val += ",";
                        int v = df3.samp_int[col].i_int[df3_idx];
                        if (v != MISSING_INT) {
                            group_val += std::to_string(v);
                            any = true;
                        } else {
                            group_val += ".";
                        }
                        col++;
                    }
                    sample_data += any ? group_val : ".";
                }
            }

            // --- DF3 float fields ---
            // Re-group float array fields.
            {
                size_t col = 0;
                while (col < df3.samp_float.size()) {
                    std::string base = formatBase(df3.samp_float[col].name, format_numbers);
                    if (!sample_data.empty()) sample_data += ":";

                    // Each column is one position of the array: a missing value keeps its place
                    std::string group_val = "";
                    bool any = false;
                    while(col < df3.samp_float.size() && formatBase(df3.samp_float[col].name, format_numbers) == base){
                        if(!group_val.empty()) group_val += ",";

                        float v = (float)df3.samp_float[col].i_float[df3_idx];

                        if(!std::isnan(v)){
                            group_val += std::to_string(v);
                            any = true;
                        }else{
                            group_val += ".";
                        }
                        col++;
                    }
                    sample_data+=any ? group_val : ".";
                }
            }

            // --- DF3 string fields ---
            // Append single string fields
            for (size_t col = 0; col < df3.samp_string.size(); col++) {
                if (!sample_data.empty()) sample_data += ":";
                const std::string& s_val = df3.samp_string[col].i_string[df3_idx];
                sample_data += s_val.empty() ? "." : s_val;
            }

            // --- DF4 per-allele fields ---
            // Extract sample-specific and allele-specific data
            if (start_df4 != -1) {
                // DF4 ints (grouped)
                {
                    size_t col = 0;
                    while (col < df4.samp_int.size()) {
                        std::string base = formatBase(df4.samp_int[col].name, format_numbers);
                        std::string field_val = "";

                        // Identify the bounds of the current group (e.g., AD0, AD1)
                        size_t group_start = col;
                        while (col < df4.samp_int.size() && formatBase(df4.samp_int[col].name, format_numbers) == base) {
                            col++;
                        }

                        // Iterate DF4 rows specifically for this variant and sample
                        size_t k = static_cast<size_t>(start_df4);
                        while (k < df4.var_id.size() && df4.var_id[k] == current_var_id) {
                            if (df4.samp_id[k] == static_cast<unsigned short>(i)) {
                                for (size_t c = group_start; c < col; c++) {
                                    if (!field_val.empty()) field_val += ",";
                                    int v = df4.samp_int[c].i_int[k];
                                    field_val += (v != MISSING_INT) ? std::to_string(v) : ".";
                                }
                            }
                            k++;
                        }
                        if (!sample_data.empty()) sample_data += ":";
                        sample_data += field_val.empty() ? "." : field_val;
                    }
                }

                // DF4 floats (grouped)
                {
                    size_t col = 0;
                    while (col < df4.samp_float.size()) {
                        std::string base = formatBase(df4.samp_float[col].name, format_numbers);
                        std::string field_val = "";

                        size_t group_start = col;
                        while (col < df4.samp_float.size() && formatBase(df4.samp_float[col].name, format_numbers) == base) {
                            col++;
                        }
                        size_t k = static_cast<size_t>(start_df4);
                        while (k < df4.var_id.size() && df4.var_id[k] == current_var_id) {
                            if (df4.samp_id[k] == static_cast<unsigned short>(i)) {
                                for (size_t c = group_start; c < col; c++) {
                                    if (!field_val.empty()) field_val += ",";
                                    float v = (float)df4.samp_float[c].i_float[k];
                                    field_val += (!std::isnan(v)) ? std::to_string(v) : ".";
                                }
                            }
                            k++;
                        }
                        if (!sample_data.empty()) sample_data += ":";
                        sample_data += field_val.empty() ? "." : field_val;
                    }
                }

                // DF4 strings (no grouping)
                for (size_t col = 0; col < df4.samp_string.size(); col++) {
                    std::string field_val = "";
                    size_t k = static_cast<size_t>(start_df4);
                    while (k < df4.var_id.size() && df4.var_id[k] == current_var_id) {
                        if (df4.samp_id[k] == static_cast<unsigned short>(i)) {
                            if (!field_val.empty()) field_val += ",";
                            const std::string& s_val = df4.samp_string[col].i_string[k];
                            field_val += s_val.empty() ? "." : s_val;
                        }
                        k++;
                    }
                    if (!sample_data.empty()) sample_data += ":";
                    sample_data += field_val.empty() ? "." : field_val;
                }
            }


            // Trailing missing fields can be dropped (VCF spec): "./.:.:." is written as "./."
            while(sample_data.size() >= 2  && sample_data.compare(sample_data.size() - 2,2,":.") == 0)
                sample_data.resize(sample_data.size() - 2);

            // Append the fully constructed sample format string to the line
            ss << sample_data;
        }
    }

    return ss.str();
}