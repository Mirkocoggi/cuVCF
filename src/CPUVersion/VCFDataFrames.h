
#ifndef VCF_DATAFRAMES_H
#define VCF_DATAFRAMES_H

#include <Imath/half.h>
#include <vector>
#include <string>
#include <map>


/**
 * @brief Constants defining the types of VCF fields
 * 
 * These constants are used to identify different types of fields in the VCF file,
 * including basic types (FLAG, INT, FLOAT, STRING) and their variants for ALT
 * and FORMAT sections.
 */
const int FLAG = 0;
const int INT = 1;
const int FLOAT = 2;
const int STRING = 3;
const int INT_ALT = 4;
const int FLOAT_ALT = 5;
const int STRING_ALT = 6;
const int STRING_FORMAT = 8;
const int INT_FORMAT = 9;
const int FLOAT_FORMAT = 10;
const int STRING_FORMAT_ALT = 11;
const int INT_FORMAT_ALT = 12;
const int FLOAT_FORMAT_ALT = 13;

using namespace std;

struct info_flag
{
    vector<bool> i_flag;
    string name;
};

struct info_string
{
    vector<string> i_string;
    string name;
};

struct info_float
{
    vector<half> i_float;
    string name;
};

struct info_int
{
    vector<int> i_int;
    string name;
};

struct samp_Flag
{
    vector<bool> i_flag;
    string name;
    int numb;
};

struct samp_String
{
    vector<string> i_string;
    string name;
    int numb;
};

struct samp_Float
{
    vector<half> i_float;
    string name;
    int numb;
};

struct samp_Int
{
    vector<int> i_int;
    string name;
    int numb;
};

struct samp_GT 
{
    vector<char> GT;
    int numb;    
};

struct header_element
{
    vector<string> ID;
    vector<string> Number;
    vector<string> Type;
    int total_values=0;
    int alt_values=0;
    int no_alt_values=0;
    int ints_alt=0;
    int floats_alt=0;
    int strings_alt=0;
    int flags_alt=0;
    int ints=0;
    int floats=0;
    int strings=0;
    int flags=0;
    bool hasGT = false;
    char numGT = 0;
};

/**
 * @class var_columns_df
 * @brief Column-oriented data frame for VCF variant data
 * 
 * Stores VCF data in a column-oriented format for efficient processing and memory usage.
 * Each field from the VCF file (CHROM, POS, ID, etc.) is stored in its own vector,
 * allowing for better cache utilization and SIMD processing.
 */

class var_columns_df
{
public:
    /** @brief Vector storing variant numbers for each entry */
    vector<unsigned int> var_number;
    
    /** 
     * @brief Maps chromosome names to compact char codes
     * @details Provides efficient storage by mapping chromosome names (e.g., "chr1") 
     *          to single character codes
     */
    std::map<std::string, char> chrom_map;
    
    /** @brief Vector storing chromosome codes sequentially */
    vector<char> chrom;
    
    /** @brief Vector storing genomic positions */
    vector<unsigned int> pos;
    
    /** @brief Vector storing variant identifiers */
    vector<string> id;
    
    /** @brief Vector storing reference alleles */
    vector<string> ref;
    
    /** @brief Vector storing quality scores in half precision */
    vector<half> qual;

    /** @brief Maps filter names to compact char representations */
    std::map<std::string, char> filter_map;
    
    /** @brief Vector storing filter status codes */
    vector<char> filter;
    
    /** @brief Vector storing float-type INFO fields */
    vector<info_float> in_float;
    
    /** @brief Vector storing flag-type INFO fields */
    vector<info_flag> in_flag;
    
    /** @brief Vector storing string-type INFO fields */
    vector<info_string> in_string;
    
    /** @brief Vector storing integer-type INFO fields */
    vector<info_int> in_int;
    
    /** @brief Maps INFO field names to their indices */
    map<string,int> info_map1;
 
};

class alt_columns_df
{
    public:
    vector<unsigned int> var_id;
    vector<char> alt_id;
    vector<string> alt;
    vector<info_float> alt_float;
    vector<info_flag> alt_flag; //non gestite per ora
    vector<info_string> alt_string;
    vector<info_int> alt_int;
    int numAlt;
};

class sample_columns_df //aka df3
{
    public:
    vector<unsigned int> var_id;
    vector<unsigned short> samp_id;
    vector<samp_Float> samp_float;
    vector<samp_Flag> samp_flag;
    vector<samp_String> samp_string;
    vector<samp_Int> samp_int;
    std::map<std::string, unsigned short> sampNames;
    map<string, char> GTMap;
    vector<samp_GT> sample_GT;
    int numSample; //numero di sample per riga

    void initMapGT(){
        int value = 0;
        // First half of the map from 0|0 to 10|10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "|" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
        // Second half of the map from 0/0 to 10/10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "/" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
    }

    std::string getGTStringFromChar(char gtChar) const {
        for (const auto& pair : GTMap) {
            if (pair.second == gtChar) {
                return pair.first;
            }
        }
        return "Not found";
    }
};

class alt_format_df //aka df4 in progress
{
    public:
    vector<unsigned int> var_id;
    vector<unsigned short> samp_id;
    vector<char> alt_id;
    vector<samp_Float> samp_float;
    vector<samp_Flag> samp_flag;
    vector<samp_String> samp_string;
    vector<samp_Int> samp_int;
    std::map<std::string, unsigned short> sampNames;
    samp_GT sample_GT;
    map<string, char> GTMap;
    int numSample; 


    void initMapGT(){
        int value = 0;
        // First half of the map from 0|0 to 10|10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "|" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
        // Second half of the map from 0/0 to 10/10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "/" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
    }

    std::string getGTStringFromChar(char gtChar) const {
        for (const auto& pair : GTMap) {
            if (pair.second == gtChar) {
                return pair.first;
            }
        }
        return "Not found";
    }
};

#endif