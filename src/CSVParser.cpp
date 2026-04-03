#include "CSVParser.h"
#include "Utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>

// Parse the VCF header file: accumulate raw text and build the field type map
void CSVParser::parseHeader(){
    std::ifstream file(header_path);
    if(!file) throw std::runtime_error("Header not found: " + header_path);

    std::string line;
    while (std::getline(file, line)){
        // Accumulate every line for verbatim output
        header_text += line + "\n";

        // Extract ID and Type from ##INFO and ##FORMAT meta-lines
        if(line.find("##INFO=")==0 || line.find("##FORMAT=")==0){
            size_t id_pos = line.find("ID=");
            size_t id_end = line.find(",", id_pos);
            std::string id = line.substr(id_pos + 3, id_end - id_pos - 3);

            size_t type_pos = line.find("Type=");
            size_t type_end = line.find(",", type_pos);
            std::string type = line.substr(type_pos + 5, type_end - type_pos - 5);

            field_types[id] = type;
        }
    }
}

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
        tmp_info.name = stripTypePrefix(headers[i]);
        df1.in_string.push_back(tmp_info);
    }    

    // Read variant records
    while (std::getline(file, line)) {
        std::vector<std::string> row = splitLine(line, ',');
        if (row.size() != headers.size()) 
            throw std::runtime_error("Line size mismatch in DF1 at line " + std::to_string(df1.var_number.size() + 2) 
            + " (expected " + std::to_string(headers.size()) + ", got " + std::to_string(row.size()) + ")");
        
        // Populate core fixed fields
        df1.var_number.push_back(std::stoul(row[0]));              
        df1.chrom.push_back(row[1]);
        df1.pos.push_back(std::stoul(row[2]));                     
        
        // Handle ID, REF, QUAL, FILTER
        std::string id_val = (row[3].empty() || row[3] == " ") ? "." : row[3];
        df1.id.push_back(id_val);
        
        df1.ref.push_back(row[4]);                                 

        float qual_val = -1.0f;
        if (!row[5].empty() && row[5] != " " && row[5] != ".") {
            try { qual_val = std::stof(row[5]); }
            catch (const std::exception&) { qual_val = -1.0f; } // valore non numerico (es. "NA")
        }
        df1.qual.push_back(qual_val);

        df1.filter.push_back(row[6]);

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
        tmp_info_alt.name = stripTypePrefix(headers[i]);
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

    // Init dynamic format containers (starts at col 3)
    size_t data_start = (headers.size() > 2 && headers[2] == "sample_name") ? 3 : 2;
    for(size_t i = data_start; i < headers.size(); i++) {
        samp_String tmp_format;
        tmp_format.name = stripTypePrefix(headers[i]);
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

        // Collect unique sample names in order
        if(headers[2] == "sample_name"){
            if(current_samp_id >= static_cast<unsigned short>(samp_names.size())){
                samp_names.push_back(row[2]);
            }
        }

        // Update total sample count
        if (current_samp_id >= df3.numSample) {
            df3.numSample = current_samp_id + 1;
        }

        // Distribute remaining dynamic fields
        for(size_t i = data_start; i < row.size(); i++) {
            df3.samp_string[i - data_start].i_string.push_back(row[i]);
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

    size_t data_start = 0;
    for(size_t i = 0; i < headers.size(); i++){
        if(headers[i] == "var_id" || headers[i] == "samp_id" || 
           headers[i] == "sample_name" || headers[i] == "alt_id"){
            data_start = i + 1;
        } else {
            break;
        }
    }

    // Find key column positions by name
    size_t samp_id_col = 0, alt_id_col = 0;
    for(size_t i = 0; i < headers.size(); i++){
        if(headers[i] == "samp_id")    samp_id_col = i;
        if(headers[i] == "alt_id")     alt_id_col  = i;
    }

    // Init dynamic format containers (starts after key columns)
    for(size_t i = data_start; i < headers.size(); i++) {
        samp_String tmp_format;
        tmp_format.name = stripTypePrefix(headers[i]);
        df4.samp_string.push_back(tmp_format);
    }    

    // Read sample-allele records
    while (std::getline(file, line)) {
        std::vector<std::string> row = splitLine(line, ',');
        if (row.size() != headers.size()) throw std::runtime_error("Line size mismatch in DF4");

        // Populate 3 fixed key columns
        df4.var_id.push_back(std::stoul(row[0]));
        df4.samp_id.push_back(static_cast<unsigned short>(std::stoi(row[samp_id_col])));
        df4.alt_id.push_back(static_cast<char>(std::stoi(row[alt_id_col])));

        // Distribute dynamic fields
        for(size_t i = data_start; i < row.size(); i++) {
            df4.samp_string[i - data_start].i_string.push_back(row[i]);
        }
    }
}

// Constructor
CSVParser::CSVParser(
    const std::string& df1_path,
    const std::string& df2_path,
    const std::string& df3_path,
    const std::string& df4_path,
    const std::string& header_path
) : df1_path(df1_path),
    df2_path(df2_path),
    df3_path(df3_path),
    df4_path(df4_path),
    header_path(header_path) {}

// Orchestrator: Init maps and load all DFs
void CSVParser::loadAll(var_columns_df& df1,
                        alt_columns_df& df2,
                        sample_columns_df& df3,
                        alt_format_df& df4){
    

    parseHeader();
    parseDF1(df1);
    parseDF2(df2);
    parseDF3(df3);
    mergeFields(df3);
    moveGTFirst(df3.samp_string);
    parseDF4(df4);
    moveGTFirst(df4.samp_string);
}