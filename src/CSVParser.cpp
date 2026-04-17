#include "CSVParser.h"
#include "Utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>

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


// Parse VCF header file
void CSVParser::parseHeader() {
    std::ifstream file(header_path);
    if (!file) throw std::runtime_error("Header not found: " + header_path);

    std::string line;
    while (std::getline(file, line)) {
        header_text += line + "\n";

        if (line.find("##INFO=") == 0 || line.find("##FORMAT=") == 0) {
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

// Parse DF1: variant core + typed INFO fields
void CSVParser::parseDF1(var_columns_df& df1) {
    std::ifstream file(df1_path);
    if (!file) throw std::runtime_error("File not found: " + df1_path);

    std::string line;
    if (!std::getline(file, line)) return;
    std::vector<std::string> headers = splitLine(line, ',');

    // Build column mappings for dynamic fields (col 7+)
    std::vector<ColumnMapping> col_map;
    for (size_t i = 7; i < headers.size(); i++) {
        std::string name = stripTypePrefix(headers[i]);         // "int_AD0" → "AD0" (nome nel container)
        std::string baseName = stripTrailingDigits(name);        // "AD0" → "AD" (per lookup in field_types)

        std::string type;
        if (field_types.count(baseName)) {
            type = field_types[baseName];
        } else if (field_types.count(name)) {
            type = field_types[name];
        } else {
            type = "String";
        }

        ColumnMapping cm;
        cm.type = type;

        if (type == "Flag") {
            info_flag tmp; 
            tmp.name = name;
            df1.in_flag.push_back(tmp);
            cm.index = df1.in_flag.size() - 1;
        } else if (type == "Integer") {
            info_int tmp; 
            tmp.name = name;
            df1.in_int.push_back(tmp);
            cm.index = df1.in_int.size() - 1;
        } else if (type == "Float") {
            info_float tmp; 
            tmp.name = name;
            df1.in_float.push_back(tmp);
            cm.index = df1.in_float.size() - 1;
        } else {
            info_string tmp; 
            tmp.name = name;
            df1.in_string.push_back(tmp);
            cm.index = df1.in_string.size() - 1;
        }
        col_map.push_back(cm);
    }

    // Read data rows
    while (std::getline(file, line)) {
        std::vector<std::string> row = splitLine(line, ',');
        if (row.size() != headers.size())
            throw std::runtime_error("Line size mismatch in DF1");

        df1.var_number.push_back(std::stoul(row[0]));

        // CHROM: store as char via chrom_map
        std::string chrom_str = row[1];
        if (df1.chrom_map.find(chrom_str) == df1.chrom_map.end()) {
            df1.chrom_map[chrom_str] = static_cast<char>(df1.chrom_map.size());
        }
        df1.chrom.push_back(df1.chrom_map[chrom_str]);

        df1.pos.push_back(std::stoul(row[2]));

        std::string id_val = (row[3].empty() || row[3] == " ") ? "." : row[3];
        df1.id.push_back(id_val);

        df1.ref.push_back(row[4]);

        // QUAL: store as half
        half qual_val = (half)(-1.0f);
        if (!row[5].empty() && row[5] != " " && row[5] != ".") {
            try { qual_val = (half)std::stof(row[5]); }
            catch (const std::exception&) { qual_val = (half)(-1.0f); }
        }
        df1.qual.push_back(qual_val);

        // FILTER: store as char via filter_map
        std::string filter_str = row[6];
        if (filter_str.empty()) filter_str = ".";
        if (df1.filter_map.find(filter_str) == df1.filter_map.end()) {
            df1.filter_map[filter_str] = static_cast<char>(df1.filter_map.size());
        }
        df1.filter.push_back(df1.filter_map[filter_str]);

        // Route dynamic fields to typed vectors
        for (size_t i = 7; i < row.size(); i++) {
            const auto& cm = col_map[i - 7];
            const std::string& val = row[i];

            if (cm.type == "Flag") {
                bool flag_val = (!val.empty() && val != "0" && val != ".");
                df1.in_flag[cm.index].i_flag.push_back(flag_val);
            } else if (cm.type == "Integer") {
                int int_val = -1;
                if (!val.empty() && val != "." && val != " ") {
                    try { int_val = std::stoi(val); }
                    catch (const std::exception&) { int_val = -1; }
                }
                df1.in_int[cm.index].i_int.push_back(int_val);
            } else if (cm.type == "Float") {
                half float_val = (half)-1.0f;
                if (!val.empty() && val != "." && val != " ") {
                    try { float_val = (half)std::stof(val); }
                    catch (const std::exception&) { float_val = (half)-1.0f; }
                }
                df1.in_float[cm.index].i_float.push_back(float_val);
            } else {
                df1.in_string[cm.index].i_string.push_back(val);
            }
        }
    }
}

// Parse DF2: alternative alleles + typed per-allele INFO
void CSVParser::parseDF2(alt_columns_df& df2) {
    std::ifstream file(df2_path);
    if (!file) throw std::runtime_error("File not found: " + df2_path);

    std::string line;
    if (!std::getline(file, line)) return;
    std::vector<std::string> headers = splitLine(line, ',');

    // Build column mappings (col 3+)
    std::vector<ColumnMapping> col_map;
    for (size_t i = 3; i < headers.size(); i++) {
        std::string name = stripTypePrefix(headers[i]);         // "int_AD0" → "AD0" (nome nel container)
        std::string baseName = stripTrailingDigits(name);        // "AD0" → "AD" (per lookup in field_types)

        std::string type;
        if (field_types.count(baseName)) {
            type = field_types[baseName];
        } else if (field_types.count(name)) {
            type = field_types[name];
        } else {
            type = "String";
        }

        ColumnMapping cm;
        cm.type = type;

        if (type == "Float") {
            info_float tmp; 
            tmp.name = name;
            df2.alt_float.push_back(tmp);
            cm.index = df2.alt_float.size() - 1;
        } else if (type == "Integer") {
            info_int tmp; 
            tmp.name = name;
            df2.alt_int.push_back(tmp);
            cm.index = df2.alt_int.size() - 1;
        } else if (type == "Flag") {
            info_flag tmp; 
            tmp.name = name;
            df2.alt_flag.push_back(tmp);
            cm.index = df2.alt_flag.size() - 1;
        } else {
            info_string tmp; 
            tmp.name = name;
            df2.alt_string.push_back(tmp);
            cm.index = df2.alt_string.size() - 1;
        }
        col_map.push_back(cm);
    }

    // Read data rows
    while (std::getline(file, line)) {
        std::vector<std::string> row = splitLine(line, ',');
        if (row.size() != headers.size())
            throw std::runtime_error("Line size mismatch in DF2");

        df2.var_id.push_back(std::stoul(row[0]));
        df2.alt_id.push_back(static_cast<char>(std::stoi(row[1])));
        df2.alt.push_back(row[2]);

        for (size_t i = 3; i < row.size(); i++) {
            const auto& cm = col_map[i - 3];
            const std::string& val = row[i];

            if (cm.type == "Float") {
                half v = (half)-1.0f;
                if (!val.empty() && val != ".") {
                    try { v = (half)std::stof(val); }
                    catch (const std::exception&) {}
                }
                df2.alt_float[cm.index].i_float.push_back(v);
            } else if (cm.type == "Integer") {
                int v = -1;
                if (!val.empty() && val != ".") {
                    try { v = std::stoi(val); }
                    catch (const std::exception&) {}
                }
                df2.alt_int[cm.index].i_int.push_back(v);
            } else if (cm.type == "Flag") {
                bool v = (!val.empty() && val != "0" && val != ".");
                df2.alt_flag[cm.index].i_flag.push_back(v);
            } else {
                df2.alt_string[cm.index].i_string.push_back(val);
            }
        }
    }
}

// Parse DF3: sample core + GT encoding
void CSVParser::parseDF3(sample_columns_df& df3) {
    std::ifstream file(df3_path);
    if (!file) throw std::runtime_error("File not found: " + df3_path);

    std::string line;
    if (!std::getline(file, line)) return;
    std::vector<std::string> headers = splitLine(line, ',');

    // Determine data start (skip var_id, samp_id, sample_name)
    size_t data_start = (headers.size() > 2 && headers[2] == "sample_name") ? 3 : 2;

    // Init GT map
    df3.initMapGT();
    df3.numSample = 0;

    // Build column mappings
    std::vector<ColumnMapping> col_map;
    for (size_t i = data_start; i < headers.size(); i++) {
        std::string name = stripTypePrefix(headers[i]);         // "int_AD0" → "AD0" (nome nel container)
        std::string baseName = stripTrailingDigits(name);        // "AD0" → "AD" (per lookup in field_types)

        std::string type;
        if (isGTField(baseName)) {
            type = "GT";
        } else {
            if (field_types.count(baseName)) {
                type = field_types[baseName];
            } else if (field_types.count(name)) {
                type = field_types[name];
            } else {
                type = "String";
            }
        }

        ColumnMapping cm;
        cm.type = type;

        if (type == "GT") {
            samp_GT tmp; 
            tmp.numb = 1;
            df3.sample_GT.push_back(tmp);
            cm.index = df3.sample_GT.size() - 1;
        } else if (type == "Float") {
            samp_Float tmp; 
            tmp.name = name; 
            tmp.numb = 1;
            df3.samp_float.push_back(tmp);
            cm.index = df3.samp_float.size() - 1;
        } else if (type == "Integer") {
            samp_Int tmp; 
            tmp.name = name; 
            tmp.numb = 1;
            df3.samp_int.push_back(tmp);
            cm.index = df3.samp_int.size() - 1;
        } else if (type == "Flag") {
            samp_Flag tmp; 
            tmp.name = name; 
            tmp.numb = 1;
            df3.samp_flag.push_back(tmp);
            cm.index = df3.samp_flag.size() - 1;
        } else {
            samp_String tmp; 
            tmp.name = name; 
            tmp.numb = 1;
            df3.samp_string.push_back(tmp);
            cm.index = df3.samp_string.size() - 1;
        }
        col_map.push_back(cm);
    }

    // Read data rows
    while (std::getline(file, line)) {
        std::vector<std::string> row = splitLine(line, ',');
        if (row.size() != headers.size())
            throw std::runtime_error("Line size mismatch in DF3");

        df3.var_id.push_back(std::stoul(row[0]));

        unsigned short current_samp_id = static_cast<unsigned short>(std::stoi(row[1]));
        df3.samp_id.push_back(current_samp_id);

        // Collect sample names
        if (headers[2] == "sample_name") {
            std::string sname = row[2];
            if (df3.sampNames.find(sname) == df3.sampNames.end()) {
                df3.sampNames[sname] = current_samp_id;
            }
        }

        if (current_samp_id >= df3.numSample) {
            df3.numSample = current_samp_id + 1;
        }

        // Route data to typed vectors
        for (size_t i = data_start; i < row.size(); i++) {
            const auto& cm = col_map[i - data_start];
            const std::string& val = row[i];

            if (cm.type == "GT") {
                // Encode GT string to char via GTMap
                char gt_code = -1;
                if (df3.GTMap.find(val) != df3.GTMap.end()) {
                    gt_code = df3.GTMap[val];
                }
                df3.sample_GT[cm.index].GT.push_back(gt_code);
            } else if (cm.type == "Float") {
                half v = (half)-1.0f;
                if (!val.empty() && val != ".") {
                    try { v = (half)std::stof(val); }
                    catch (const std::exception&) {}
                }
                df3.samp_float[cm.index].i_float.push_back(v);
            } else if (cm.type == "Integer") {
                int v = -1;
                if (!val.empty() && val != ".") {
                    try { v = std::stoi(val); }
                    catch (const std::exception&) {}
                }
                df3.samp_int[cm.index].i_int.push_back(v);
            } else if (cm.type == "Flag") {
                bool v = (!val.empty() && val != "0" && val != ".");
                df3.samp_flag[cm.index].i_flag.push_back(v);
            } else {
                df3.samp_string[cm.index].i_string.push_back(val);
            }
        }
    }
}

// Parse DF4: sample × allele + GT encoding
void CSVParser::parseDF4(alt_format_df& df4) {
    std::ifstream file(df4_path);
    if (!file) throw std::runtime_error("File not found: " + df4_path);

    std::string line;
    if (!std::getline(file, line)) return;
    std::vector<std::string> headers = splitLine(line, ',');

    // Find data_start and key column positions
    size_t data_start = 0;
    size_t samp_id_col = 0, alt_id_col = 0;
    for (size_t i = 0; i < headers.size(); i++) {
        if (headers[i] == "var_id" || headers[i] == "samp_id" ||
            headers[i] == "sample_name" || headers[i] == "alt_id") {
            data_start = i + 1;
        } else {
            break;
        }
    }
    for (size_t i = 0; i < headers.size(); i++) {
        if (headers[i] == "samp_id") samp_id_col = i;
        if (headers[i] == "alt_id")  alt_id_col = i;
    }

    // Init GT map
    df4.initMapGT();

    // Build column mappings
    std::vector<ColumnMapping> col_map;
    for (size_t i = data_start; i < headers.size(); i++) {
        std::string name = stripTypePrefix(headers[i]);
        std::string baseName = stripTrailingDigits(name);
        std::string type;
        if (isGTField(baseName)) {
            type = "GT";
        } else {
            if (field_types.count(baseName)) {
                type = field_types[baseName];
            } else if (field_types.count(name)) {
                type = field_types[name];
            } else {
                type = "String";
            }
        }

        ColumnMapping cm;
        cm.type = type;

        if (type == "GT") {
            // DF4 has a single samp_GT, not a vector
            cm.index = 0;
        } else if (type == "Float") {
            samp_Float tmp; 
            tmp.name = name; 
            tmp.numb = 1;
            df4.samp_float.push_back(tmp);
            cm.index = df4.samp_float.size() - 1;
        } else if (type == "Integer") {
            samp_Int tmp; 
            tmp.name = name;
            tmp.numb = 1;
            df4.samp_int.push_back(tmp);
            cm.index = df4.samp_int.size() - 1;
        } else if (type == "Flag") {
            samp_Flag tmp; 
            tmp.name = name; 
            tmp.numb = 1;
            df4.samp_flag.push_back(tmp);
            cm.index = df4.samp_flag.size() - 1;
        } else {
            samp_String tmp; 
            tmp.name = name; 
            tmp.numb = 1;
            df4.samp_string.push_back(tmp);
            cm.index = df4.samp_string.size() - 1;
        }
        col_map.push_back(cm);
    }

    // Read data rows
    while (std::getline(file, line)) {
        std::vector<std::string> row = splitLine(line, ',');
        if (row.size() != headers.size())
            throw std::runtime_error("Line size mismatch in DF4");

        df4.var_id.push_back(std::stoul(row[0]));
        df4.samp_id.push_back(static_cast<unsigned short>(std::stoi(row[samp_id_col])));
        df4.alt_id.push_back(static_cast<char>(std::stoi(row[alt_id_col])));

        for (size_t i = data_start; i < row.size(); i++) {
            const auto& cm = col_map[i - data_start];
            const std::string& val = row[i];

            if (cm.type == "GT") {
                char gt_code = -1;
                if (df4.GTMap.find(val) != df4.GTMap.end()) {
                    gt_code = df4.GTMap[val];
                }
                df4.sample_GT.GT.push_back(gt_code);
            } else if (cm.type == "Float") {
                half v = (half)-1.0f;
                if (!val.empty() && val != ".") {
                    try { v = (half)std::stof(val); }
                    catch (const std::exception&) {}
                }
                df4.samp_float[cm.index].i_float.push_back(v);
            } else if (cm.type == "Integer") {
                int v = -1;
                if (!val.empty() && val != ".") {
                    try { v = std::stoi(val); }
                    catch (const std::exception&) {}
                }
                df4.samp_int[cm.index].i_int.push_back(v);
            } else if (cm.type == "Flag") {
                bool v = (!val.empty() && val != "0" && val != ".");
                df4.samp_flag[cm.index].i_flag.push_back(v);
            } else {
                df4.samp_string[cm.index].i_string.push_back(val);
            }
        }
    }
}

// Orchestrator
void CSVParser::loadAll(var_columns_df& df1,
                              alt_columns_df& df2,
                              sample_columns_df& df3,
                              alt_format_df& df4) {
    parseHeader();
    parseDF1(df1);
    parseDF2(df2);
    parseDF3(df3);
    parseDF4(df4);
}