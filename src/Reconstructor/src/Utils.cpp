#include "Utils.h"
#include <string>
#include <vector>
#include <cctype>
#include <stdexcept>

// Split CSV line, removing trailing Windows carriage returns
std::vector<std::string> splitLine(const std::string& line, char delimiter) {
    std::string clean_line = line;
    
    // Remove trailing \r from Windows line endings
    if (!clean_line.empty() && clean_line.back() == '\r') {
        clean_line.pop_back(); 
    }

    std::vector<std::string> tokens;
    std::string current;
    bool in_quotes = false;

    for (size_t i = 0; i < clean_line.size(); i++) {
        char c = clean_line[i];
        if (c == '"') {
            // Toggle quote mode — delimiters inside quotes are ignored
            in_quotes = !in_quotes;
        } else if (c == delimiter && !in_quotes) {
            // End of field — push and reset
            tokens.push_back(current);
            current.clear();
        } else {
            current += c;
        }
    }
    // Push the last field (after the final delimiter)
    tokens.push_back(current);

    return tokens;
}

// Strips trailing digits from a field name (e.g. "AD0" -> "AD")
std::string stripTrailingDigits(const std::string& field_name){
    std::string result = field_name;
    // Remove digits from the end one by one
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
        // Only strip known cuVCF type prefixes
        if(prefix == "flag" || prefix == "float" || prefix == "int" ||
           prefix == "string" || prefix == "char"){
            return name.substr(underscore + 1);
        }
    }
    // No known prefix found — return unchanged
    return name;
}

bool isGTField(const std::string& name){
    // cuVCF may name GT as "GT" or "GT0" depending on the dataset
    return name == "GT" || name == "GT0";
}


