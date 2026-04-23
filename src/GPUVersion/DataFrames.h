/**
 * @file DataFrames.h
 * @brief Contains data frame classes for handling VCF column data.
 *
 * The file defines several classes:
 * - var_columns_df (DF1): Manages variant columns such as chromosome, position, etc.
 * - alt_columns_df (DF2): Handles alternative allele columns.
 * - sample_columns_df (DF3): Manages sample-related columns.
 * - alt_format_df (DF4): Handles formatted alternative allele data for samples.
 *
 * Each class provides methods for initialization, cloning, and printing of data.
 */

#ifndef DATA_FRAMES_H
#define DATA_FRAMES_H

#include "DataStructures.h"
#include <chrono>
#include <boost/algorithm/string.hpp>
#include <cuda_runtime.h>     
#include <cuda_fp16.h>  
#include <map>
#include <iostream>
#include <vector>

using namespace std;

/**
 * @class var_columns_df
 * @brief Data frame for variant columns.
 *
 * This class manages variant-specific data such as variant numbers, chromosome,
 * positions, IDs, reference alleles, quality scores, filter information, and additional info fields.
 */
class var_columns_df //DF1
{
public:
    /// Vector of variant numbers.
    vector<unsigned int> var_number;
    /// Map from chromosome string to a unique unsigned char code.
    std::map<std::string, unsigned char> chrom_map;
    /// Vector of chromosome codes.
    vector<char> chrom;
    /// Vector of variant positions.
    vector<unsigned int> pos;
    /// Vector of variant IDs.
    vector<string> id;
    /// Vector of reference alleles.
    vector<string> ref;
    /// Vector of quality scores in half precision.
    vector<__half> qual;
    /// Map from filter string to a unique char code.
    std::map<std::string, char> filter_map;
    /// Vector of filter codes.
    vector<char> filter;
    /// Vector of float info fields.
    vector<info_float> in_float;
    /// Vector of flag info fields.
    vector<info_flag> in_flag;
    /// Vector of string info fields.
    vector<info_string> in_string;
    /// Vector of integer info fields.
    vector<info_int> in_int;
    /// Map for additional info, mapping field names to integer codes.
    map<string, int> info_map1;
};

/**
 * @class alt_columns_df
 * @brief Data frame for alternative allele columns.
 *
 * This class stores variant IDs, alternative allele IDs, allele strings, and
 * associated information (floats, flags, strings, and integers) for alternative alleles.
 */
class alt_columns_df //DF2
{
    public:
    /// Vector of variant IDs corresponding to each alternative allele entry.
    vector<unsigned int> var_id;
    /// Vector of alternative allele IDs.
    vector<unsigned char> alt_id;
    /// Vector of alternative allele strings.
    vector<string> alt;
    /// Vector of float information related to alternative alleles.
    vector<info_float> alt_float;
    /// Vector of flag information (not yet handled) for alternative alleles.
    vector<info_flag> alt_flag;
    /// Vector of string information related to alternative alleles.
    vector<info_string> alt_string;
    /// Vector of integer information related to alternative alleles.
    vector<info_int> alt_int;
    /// Number of alternative alleles.
    int numAlt;
};

/**
 * @class sample_columns_df
 * @brief Data frame for sample-related columns.
 *
 * This class manages sample information including variable IDs, sample IDs,
 * float, flag, string, and integer data for samples, along with sample names and genotype mapping.
 */

class sample_columns_df //aka df3
{
    public:
    /// Vector of variant IDs for sample data.
    vector<unsigned int> var_id;
    /// Vector of sample IDs.
    vector<unsigned short> samp_id;
    /// Vector of sample float data.
    vector<samp_Float> samp_float;
    /// Vector of sample flag data.
    vector<samp_Flag> samp_flag;
    /// Vector of sample string data.
    vector<samp_String> samp_string;
    /// Vector of sample integer data.
    vector<samp_Int> samp_int;
    /// Map from sample name to sample ID.
    std::map<std::string, unsigned short> sampNames;
    /// Map from genotype string to a char code.
    map<string, char> GTMap;
    /// Vector of sample genotype data.
    vector<samp_GT> sample_GT;
    /// Number of samples per row.
    int numSample;

    /**
     * @brief Initializes the genotype map for samples.
     *
     * The function builds a map for genotype strings:
     * - First half: keys from "0|0" to "10|10"
     * - Second half: keys from "0/0" to "10/10"
     * It also maps missing genotype values.
     */
    void initMapGT(){
        GTMap[".|."] = static_cast<char>(254);
        GTMap["./."] = static_cast<char>(255);
        int value = 0;
        // First __half of the map from 0|0 to 10|10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "|" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
        // Second __half of the map from 0/0 to 10/10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "/" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }       
    }

    /**
     * @brief Retrieves the genotype string corresponding to a genotype character.
     *
     * Performs a reverse lookup in the GTMap.
     * @param gtChar Genotype character code.
     * @return The genotype string if found, otherwise "Not found".
     */
    std::string getGTStringFromChar(char gtChar) const {
        for (const auto& pair : GTMap) {
            if (pair.second == gtChar) {
                return pair.first;
            }
        }
        return "Not found";
    }
};

/**
 * @class alt_format_df
 * @brief Data frame for formatted alternative allele data for samples.
 * 
 * @details This class stores formatted alternative allele information for samples,
 * including variant IDs, sample IDs, alternative allele IDs, and associated
 * sample data (float, flag, string, int). It provides methods for initialization,
 * data management and printing.
 *
 * @note The class name "df4" is deprecated and should not be used
 */
class alt_format_df //aka df4 in progress
{
    public:
    /// Vector of variant IDs.
    vector<unsigned int> var_id;
    /// Vector of sample IDs.
    vector<unsigned short> samp_id;
    /// Vector of alternative allele IDs.
    vector<char> alt_id;
    /// Vector of sample float data.
    vector<samp_Float> samp_float;
    /// Vector of sample flag data.
    vector<samp_Flag> samp_flag;
    /// Vector of sample string data.
    vector<samp_String> samp_string;
    /// Vector of sample integer data.
    vector<samp_Int> samp_int;
    /// Map from sample name to sample ID.
    std::map<std::string, unsigned short> sampNames;
    /// Sample genotype data.
    samp_GT sample_GT;
    /// Map from genotype string to a char code.
    map<string, char> GTMap;
    /// Number of samples.
    int numSample;  

    /**
     * @brief Initializes the genotype map for formatted alternative alleles.
     *
     * Builds the GTMap for genotype lookup.
     */
    void initMapGT(){
        int value = 0;
        // First __half of the map from 0|0 to 10|10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "|" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
        // Second __half of the map from 0/0 to 10/10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "/" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
    }

    /**
     * @brief Retrieves the genotype string corresponding to a genotype character.
     *
     * Performs a reverse lookup in the GTMap.
     * @param gtChar Genotype character code.
     * @return The corresponding genotype string if found; otherwise, "Not found".
     */
    std::string getGTStringFromChar(char gtChar) const {
        for (const auto& pair : GTMap) {
            if (pair.second == gtChar) {
                return pair.first;
            }
        }
        return "Not found";
    }     
    }   
};



#endif