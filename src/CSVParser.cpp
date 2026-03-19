#include "CSVParser.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/algorithm/string.hpp>

// Parse DF1: Core variant info and INFO fields
void CSVParser::parseDF1(var_columns_df& df1) {
    std::ifstream file(df1_path);
    if (!file) throw std::runtime_error("File not found: " + df1_path);
    
    std::string line;
    if (!std::getline(file, line)) return;
    std::vector<std::string> headers = splitLine(line, ',');

    // Init dynamic INFO containers (starts at col 7)
    for(size_t i = 7; i < headers.size(); i++) {
        info_string tmp_info;
        tmp_info.name = headers[i];
        df1.in_string.push_back(tmp_info);
    }    

    // Read variant records
    while (std::getline(file, line)) {
        std::vector<std::string> row = splitLine(line, ',');
        if (row.size() != headers.size()) throw std::runtime_error("Line size mismatch in DF1");
        
        // Populate core fixed fields
        df1.var_number.push_back(std::stoul(row[0]));              
        df1.chrom.push_back(static_cast<char>(std::stoi(row[1])));  
        df1.pos.push_back(std::stoul(row[2]));                     
        
        // Handle ID, REF, FILTER, QUAL
        std::string id_val = (row[3].empty() || row[3] == " ") ? "." : row[3];
        df1.id.push_back(id_val);
        
        df1.ref.push_back(row[4]);                                 

        char filter_val = (row[5].empty() || row[5] == " ") ? -1 : row[5][0];
        df1.filter.push_back(filter_val);

        float qual_val = (row[6].empty() || row[6] == " ") ? -1.0f : std::stof(row[6]);
        df1.qual.push_back(qual_val);

        // Distribute dynamic INFO fields
        for(size_t i = 7; i < row.size(); i++) {   
            df1.in_string[i-7].i_string.push_back(row[i]);
        }
    }
}

// Parse DF2: Alternative alleles
void CSVParser::parseDF2(alt_columns_df& df2) {
    std::ifstream file(df2_path);
    if (!file) throw std::runtime_error("File not found: " + df2_path);

    std::string line;
    if (!std::getline(file, line)) return;
    std::vector<std::string> headers = splitLine(line, ',');

    // Init dynamic ALT string fields (starts at col 3)
    for(size_t i = 3; i < headers.size(); i++) {
        info_string tmp_info_alt;
        tmp_info_alt.name = headers[i];
        df2.alt_string.push_back(tmp_info_alt);
    }    

    // Read allele records
    while (std::getline(file, line)) {
        std::vector<std::string> row = splitLine(line, ',');
        if (row.size() != headers.size()) throw std::runtime_error("Line size mismatch in DF2");

        // Populate fixed ALT fields
        df2.var_id.push_back(std::stoul(row[0]));
        df2.alt_id.push_back(static_cast<char>(std::stoi(row[1])));
        df2.alt.push_back(row[2]);

        // Distribute dynamic ALT string fields
        for(size_t i = 3; i < row.size(); i++) {   
            df2.alt_string[i-3].i_string.push_back(row[i]);
        }
    }
}

// Parse DF3: Sample core info (per-variant, per-sample)
void CSVParser::parseDF3(sample_columns_df& df3) {
    std::ifstream file(df3_path);
    if (!file) throw std::runtime_error("File not found: " + df3_path);

    std::string line;
    if (!std::getline(file, line)) return;
    std::vector<std::string> headers = splitLine(line, ',');

    // Ensure sample_GT vector has at least one element for indexing
    if (df3.sample_GT.empty()) {
        samp_GT initial_gt;
        df3.sample_GT.push_back(initial_gt);
    }

    // Init dynamic format containers (starts at col 3)
    for(size_t i = 3; i < headers.size(); i++) {
        samp_String tmp_format;
        tmp_format.name = headers[i];
        df3.samp_string.push_back(tmp_format);
    }    

    df3.numSample = 0; 

    // Read sample records
    while (std::getline(file, line)) {
        std::vector<std::string> row = splitLine(line, ',');
        if (row.size() != headers.size()) throw std::runtime_error("Line size mismatch in DF3");

        df3.var_id.push_back(std::stoul(row[0]));
        
        unsigned short current_samp_id = static_cast<unsigned short>(std::stoi(row[1]));
        df3.samp_id.push_back(current_samp_id);

        // Update total sample count
        if (current_samp_id >= df3.numSample) {
            df3.numSample = current_samp_id + 1;
        }

        // Extract raw byte for Genotype (GT0) at col 2
        char gt_char = row[2].empty() ? (char)0 : row[2][0];
        df3.sample_GT[0].GT.push_back(gt_char); 

        // Distribute remaining dynamic fields
        for(size_t i = 3; i < row.size(); i++) {
            df3.samp_string[i - 3].i_string.push_back(row[i]);
        }
    }
}

// Parse DF4: Sample info per allele
void CSVParser::parseDF4(alt_format_df& df4) {
    std::ifstream file(df4_path);
    if (!file) throw std::runtime_error("File not found: " + df4_path);

    std::string line;
    if (!std::getline(file, line)) return;
    std::vector<std::string> headers = splitLine(line, ',');

    // Init dynamic format containers (starts at col 3, no GT present)
    for(size_t i = 3; i < headers.size(); i++) {
        samp_String tmp_format;
        tmp_format.name = headers[i];
        df4.samp_string.push_back(tmp_format);
    }    

    // Read sample-allele records
    while (std::getline(file, line)) {
        std::vector<std::string> row = splitLine(line, ',');
        if (row.size() != headers.size()) throw std::runtime_error("Line size mismatch in DF4");

        // Populate 3 fixed key columns
        df4.var_id.push_back(std::stoul(row[0]));
        df4.samp_id.push_back(static_cast<unsigned short>(std::stoi(row[1])));
        df4.alt_id.push_back(static_cast<char>(std::stoi(row[2])));

        // Distribute dynamic fields
        for(size_t i = 3; i < row.size(); i++) {
            df4.samp_string[i - 3].i_string.push_back(row[i]);
        }
    }
}

// Split CSV line, removing trailing Windows carriage returns
std::vector<std::string> CSVParser::splitLine(const std::string& line, char delimiter) {
    std::string clean_line = line;
    
    if (!clean_line.empty() && clean_line.back() == '\r') {
        clean_line.pop_back(); 
    }

    std::vector<std::string> tokens;
    boost::split(tokens, clean_line, boost::is_any_of(std::string(1, delimiter)));
    return tokens;
}

// Constructor
CSVParser::CSVParser(
    const std::string& df1_path,
    const std::string& df2_path,
    const std::string& df3_path,
    const std::string& df4_path
) : df1_path(df1_path),
    df2_path(df2_path),
    df3_path(df3_path),
    df4_path(df4_path) {}

// Orchestrator: Init maps and load all DFs
void CSVParser::loadAll(var_columns_df& df1,
                        alt_columns_df& df2,
                        sample_columns_df& df3,
                        alt_format_df& df4){
    
    // Init hardcoded genotype dictionaries
    df3.initMapGT();
    df4.initMapGT();

    parseDF1(df1);
    parseDF2(df2);
    parseDF3(df3);
    parseDF4(df4);
}