#include "Utils.h"
#include <string>
#include <vector>
#include <cctype>
#include <stdexcept>

// Split CSV line, removing trailing Windows carriage returns
std::vector<std::string> splitLine(const std::string& line, char delimiter) {
    std::string clean_line = line;
    
    if (!clean_line.empty() && clean_line.back() == '\r') {
        clean_line.pop_back(); 
    }

    std::vector<std::string> tokens;
    std::string current;
    bool in_quotes = false;

    for (size_t i = 0; i < clean_line.size(); i++) {
        char c = clean_line[i];
        if (c == '"') {
            in_quotes = !in_quotes;
        } else if (c == delimiter && !in_quotes) {
            tokens.push_back(current);
            current.clear();
        } else {
            current += c;
        }
    }
    tokens.push_back(current);

    return tokens;
}

// Strips trailing digits from a field name (e.g. "AD0" -> "AD", "PL12" -> "PL")
std::string stripTrailingDigits(const std::string& field_name){
    std::string result = field_name;
    while (!result.empty() && std::isdigit(result.back())) {
        result.pop_back();
    }
    return result;
}

// Strips type prefix from CSV field name (e.g. "flag_SOMATIC" -> "SOMATIC")
std::string stripTypePrefix(const std::string& name) {
    size_t underscore = name.find('_');
    if(underscore != std::string::npos){
        std::string prefix = name.substr(0, underscore);
        if(prefix == "flag" || prefix == "float" || prefix == "int" ||
           prefix == "string" || prefix == "char"){
            return name.substr(underscore + 1);
        }
    }
    return name;
}

// Merges split FORMAT fields (e.g. AD0, AD1 → AD) by concatenating values with comma
void mergeFields(sample_columns_df& df3){

    if(df3.samp_string.empty()) return;
    
    // First, clean all names
    for(size_t i = 0; i < df3.samp_string.size(); i++){
        df3.samp_string[i].name = stripTrailingDigits(df3.samp_string[i].name);
    }

    // Then merge consecutive fields with same name, backwards
    for(size_t i = df3.samp_string.size() - 1; i > 0; i--){
        if(df3.samp_string[i].name == df3.samp_string[i-1].name){
            // Concatenate each value with comma
            for(size_t j = 0; j < df3.samp_string[i-1].i_string.size(); j++){
                if(df3.samp_string[i].i_string[j] == "."){
                    continue;
                }
                df3.samp_string[i-1].i_string[j] += "," + df3.samp_string[i].i_string[j];
            }
            // Remove the merged field
            df3.samp_string.erase(df3.samp_string.begin() + i);
        }
    }
}

void moveGTFirst(std::vector<samp_String>& fields){
    for(size_t i = 0; i < fields.size(); i++){
        if(isGTField(fields[i].name)){
            if(i > 0) std::swap(fields[0], fields[i]);
            return;
        }
    }
}

bool isGTField(const std::string& name){
    return name == "GT" || name == "GT0";
}