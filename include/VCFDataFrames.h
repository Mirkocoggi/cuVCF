#ifndef VCF_DATAFRAMES_H
#define VCF_DATAFRAMES_H

#include <vector>
#include <string>
#include <map>

// ===================================================================
// BASE STRUCTURES 
// ===================================================================

struct info_string {
    std::vector<std::string> i_string;
    std::string name;
};

struct samp_String {
    std::vector<std::string> i_string;
    std::string name;
    int numb;
};

struct samp_GT {
    std::vector<char> GT;
    int numb=0;    
};

// ===================================================================
// DATAFRAME CLASSES
// ===================================================================

// -------------------------------------------------------------------
// DF1: Variant core
// -------------------------------------------------------------------
class var_columns_df {
public:
    std::vector<unsigned int> var_number;
    std::vector<char> chrom;                         
    std::vector<unsigned int> pos;
    std::vector<std::string> id;
    std::vector<std::string> ref;
    std::vector<float> qual;                         

    std::vector<char> filter;
    std::vector<info_string> in_string;

    /*// Dictionaries for decoding (ID to String)
    std::vector<std::string> chrom_names;            // index = char code
    std::vector<std::string> filter_names;           // index = char code
    */
};

// -------------------------------------------------------------------
// DF2: Alternative alleles detail
// -------------------------------------------------------------------
class alt_columns_df {
public:
    std::vector<unsigned int> var_id;    
    std::vector<char> alt_id;             
    std::vector<std::string> alt;         
    std::vector<info_string> alt_string;
};

// -------------------------------------------------------------------
// DF3: Sample core 
// -------------------------------------------------------------------
class sample_columns_df {
public:
    int numSample = 0;                                         
    std::vector<unsigned int> var_id;                      
    std::vector<unsigned short> samp_id;                   
    std::vector<samp_GT> sample_GT;                         
    std::vector<samp_String> samp_string;

    /*// Dictionary for sample decoding (ID to String)
    std::vector<std::string> samp_names;                   
    */

    std::map<std::string, char> GTMap;                      

    // Initialize the genotype mapping table
    void initMapGT() {
        int value = 0;
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "|" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "/" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
    }

    // Retrieve original genotype string from char ID
    std::string getGTStringFromChar(char gtChar) const {
        for (const auto& pair : GTMap) {
            if (pair.second == gtChar) return pair.first;
        }
        return ".";
    }
};

// -------------------------------------------------------------------
// DF4: Sample x Allele 
// -------------------------------------------------------------------
class alt_format_df {
public:
    std::vector<unsigned int> var_id;        
    std::vector<unsigned short> samp_id;     
    std::vector<char> alt_id;                 
    samp_GT sample_GT;                         
    std::vector<samp_String> samp_string;

    std::map<std::string, char> GTMap;                      
    int numSample = 0;

    // Initialize the genotype mapping table
    void initMapGT() {
        int value = 0;
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "|" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "/" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
    }

    // Retrieve original genotype string from char ID
    std::string getGTStringFromChar(char gtChar) const {
        for (const auto& pair : GTMap) {
            if (pair.second == gtChar) return pair.first;
        }
        return ".";
    }
};

#endif