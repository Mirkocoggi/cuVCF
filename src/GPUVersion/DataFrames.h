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

    /**
     * @brief Prints the variant columns data frame.
     *
     * Outputs the first num_lines entries including variant number, chromosome,
     * position, ID, reference, and filter information.
     *
     * @param num_lines Number of lines to print.
     */
    void print(long num_lines) {
        std::cout << "VarID\tChrom\tPos\tID\tRef\tQUAL\tFilter\tFlag\t\tInt\t\tFloat\t\tString" << std::endl;
        long iter = (num_lines > static_cast<long>(id.size())) ? id.size() : num_lines;
        for (long i = 0; i < iter; i++) {
            // Stampa VarID (presumibilmente var_number)
            if (i < static_cast<long>(var_number.size()))
                std::cout << var_number[i] << "\t";
            else
                std::cout << "nan\t";
            
            // Stampa Chrom: costruiamo la chiave da un singolo carattere
            if (i < static_cast<long>(chrom.size())) {
                unsigned char code = chrom[i];
                bool found = false;
                for (const auto& pair : chrom_map) {
                    if (pair.second == code) {
                        std::cout << pair.first << "\t";
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    std::cout << "nan\t";
                }
            } else {
                std::cout << "nan\t";
            }
            
            
            // Stampa Pos
            if (i < static_cast<long>(pos.size()))
                std::cout << pos[i] << "\t";
            else
                std::cout << "nan\t";
            
            // Stampa ID
            if (i < static_cast<long>(id.size()))
                std::cout << id[i] << "\t";
            else
                std::cout << "nan\t";
            
            // Stampa Ref
            if (i < static_cast<long>(ref.size()))
                std::cout << ref[i] << "\t";
            else
                std::cout << "nan\t";

            // Stampa Qual
            if (i < static_cast<long>(qual.size()))
                if(__half2float(qual[i]) != -1.0f){
                    std::cout << __half2float(qual[i]) << "\t";
                }else{
                    std::cout << "." << "\t";
                }
            else
                std::cout << "nan\t";
            
            // Stampa Filter: costruiamo la chiave dal carattere in filter
            if (i < static_cast<long>(filter.size())) {
                char code = filter[i];
                bool found = false;
                for (const auto& pair : filter_map) {
                    if (pair.second == code) {
                        std::cout << pair.first << "\t";
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    std::cout << "nan\t";
                }
            } else {
                std::cout << "nan\t";
            }
            
            
            // Stampa i campi dei vettori info: Flag, Int, Float e String
            // Si assume che ogni vettore info_xxx abbia la stessa dimensione (almeno per indice i)
            // e che il campo 'name' contenga il nome della colonna.
            
            // Flag
            for (size_t j = 0; j < in_flag.size(); j++) {
                if (i < static_cast<long>(in_flag[j].i_flag.size())){
                    if(in_flag[j].i_flag[i]==1) {
                        std::cout << in_flag[j].name << ";";
                    }
                }else{
                    std::cout << in_flag[j].name << "nan ";
                }
            }
            std::cout << "\t";
            
            // Int
            for (size_t j = 0; j < in_int.size(); j++) {
                if (i < static_cast<long>(in_int[j].i_int.size()))
                    std::cout << in_int[j].name << ": " << in_int[j].i_int[i] << " ";
                else
                    std::cout << in_int[j].name << ": nan ";
            }
            std::cout << "\t";
            
            // Float
            for (size_t j = 0; j < in_float.size(); j++) {
                if (i < static_cast<long>(in_float[j].i_float.size()))
                    std::cout << in_float[j].name << ": " << static_cast<float>(in_float[j].i_float[i]) << " ";
                else
                    std::cout << in_float[j].name << ": nan ";
            }
            std::cout << "\t";
            
            // String
            for (size_t j = 0; j < in_string.size(); j++) {
                if (i < static_cast<long>(in_string[j].i_string.size()))
                    std::cout << in_string[j].name << ": " << in_string[j].i_string[i] << " ";
                else
                    std::cout << in_string[j].name << ": nan ";
            }
            std::cout << std::endl;
        }
    }
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

    /**
     * @brief Initializes the alternative columns data frame.
     *
     * The function resizes the internal vectors based on 2 * batch_size and populates
     * the float, int, and string information using a reference data frame and header info.
     *
     * @param ref Reference alt_columns_df used to copy names.
     * @param INFO Header element containing alternative info counts.
     * @param batch_size The batch size used for processing.
     */
    void init(alt_columns_df ref, header_element INFO, long batch_size){
        int numAlt = 2*batch_size;
        var_id.resize(numAlt, 0);
        alt.resize(numAlt, "\0");
        alt_id.resize(numAlt, (char)0);
        int tmp = INFO.floats_alt;
        if(tmp>0){
            info_float tmpInfoFloat;
            for(int i = 0; i<tmp; i++){
                tmpInfoFloat.name = ref.alt_float[i].name;
                tmpInfoFloat.i_float.resize(numAlt, 0.0f);
                alt_float.push_back(tmpInfoFloat);
            }
            alt_float.resize(tmp);
        }
    
        tmp = INFO.ints_alt;

        if(tmp>0){
            info_int tmpInfoInt;
            for(int i = 0; i<tmp; i++){
                tmpInfoInt.name = ref.alt_int[i].name;
                tmpInfoInt.i_int.resize(numAlt, 0);
                alt_int.push_back(tmpInfoInt);

            }
            alt_int.resize(tmp);
        }
        tmp = INFO.strings_alt;
        if(tmp>0){
            info_string tmpInfoString;
            for(int i = 0; i<tmp; i++){
                tmpInfoString.name = ref.alt_string[i].name;
                tmpInfoString.i_string.resize(numAlt, "\0");
                alt_string.push_back(tmpInfoString);

            }
            alt_string.resize(tmp);
        }       
    }

    /**
     * @brief Clones the alternative columns data frame from a reference.
     *
     * The function clones the data frame by resizing the vectors based on the
     * number of alternative values specified in the header.
     *
     * @param ref Reference alt_columns_df to clone from.
     * @param INFO Header element containing alternative info counts.
     * @param batch_size The batch size used for processing.
     */
    void clone(alt_columns_df ref, header_element INFO, long batch_size){
        int numAlt = INFO.alt_values;//sigsegv
        var_id.resize(numAlt, 0);
        alt.resize(numAlt, "\0");
        alt_id.resize(numAlt, (char)0);
        int tmp = INFO.floats_alt;
        if(tmp>0){
            info_float tmpInfoFloat;
            for(int i = 0; i<tmp; i++){
                tmpInfoFloat.name = ref.alt_float[i].name;
                tmpInfoFloat.i_float.resize(ref.alt_float[i].i_float.size(), 0.0f);
                alt_float.push_back(tmpInfoFloat);
            }
            alt_float.resize(tmp);
        }
    
        tmp = INFO.ints_alt;

        if(tmp>0){
            info_int tmpInfoInt;
            for(int i = 0; i<tmp; i++){
                tmpInfoInt.name = ref.alt_int[i].name;
                tmpInfoInt.i_int.resize(ref.alt_int[i].i_int.size(), 0);
                alt_int.push_back(tmpInfoInt);

            }
            alt_int.resize(tmp);
        }
        tmp = INFO.strings_alt;
        if(tmp>0){
            info_string tmpInfoString;
            for(int i = 0; i<tmp; i++){
                tmpInfoString.name = ref.alt_string[i].name;
                tmpInfoString.i_string.resize(ref.alt_string[i].i_string.size(), "\0");
                alt_string.push_back(tmpInfoString);

            }
            alt_string.resize(tmp);
        }       
    }
    
    /**
     * @brief Prints the alternative columns data frame.
     *
     * Outputs the first n entries including variant ID, alternative ID, alternative allele,
     * and associated float, int, and string information.
     *
     * @param n Number of entries to print.
     */
    void print(int n){
        cout << "VarID\tAltID\tAlt\tFloat\t\tInt\t\tStr" << endl;
        int iter = (n>var_id.size()) ? var_id.size() : n;
        for(int i=0; i<iter; i++){

            cout << var_id[i] << "\t" << (int)alt_id[i] << "\t" << alt[i] << "\t";
            
            for(int j=0; j < alt_float.size(); j++){
                cout << alt_float[j].name << "=" << static_cast<float>(alt_float[j].i_float[i]) << ";";
            }
        
            cout << "\t";

            for(int j=0; j < alt_int.size(); j++){
                cout << alt_int[j].name << "=" << alt_int[j].i_int[i] << ";";
            }

            cout << "\t";

            for(int j=0; j < alt_string.size(); j++){
                cout << alt_string[j].name << "=" << alt_string[j].i_string[i] << ";";
            }
            
            cout << endl;
        }
    }
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

    /**
     * @brief Prints the sample columns data frame.
     *
     * Outputs the first n entries including variant ID, sample ID (resolved to name),
     * and associated sample float, integer, string, and genotype data.
     *
     * @param n Number of entries to print.
     */
   void print(int n){
        cout << "VarID\tSampID\tFloat\t\tInt\t\tStr\t\tGT" << endl;
        
        int iter = (n>var_id.size()) ? var_id.size() : n;

        for(int i=0; i<iter; i++){
            cout << var_id[i] << "\t";
            // Reverse lookup in sampNames map to get the string associated with samp_id
            for (const auto& pair : sampNames) {
                if (pair.second == samp_id[i]) {
                    cout << pair.first << "\t";
                    break;
                }
            }

            for(int j=0; j < samp_float.size(); j++){
                cout << samp_float[j].name << "=" << static_cast<float>(samp_float[j].i_float[i]) << ";";
            }
            cout << "\t";
            
            for(int j=0; j < samp_int.size(); j++){
                cout << samp_int[j].name << "=" << samp_int[j].i_int[i] << ";";
            }
            cout << "\t";
            
            for(int j=0; j < samp_string.size(); j++){
                cout << samp_string[j].name << "=" << samp_string[j].i_string[i] << ";";
            }
            cout << "\t";

            for(int j=0; j < sample_GT[0].numb; j++){
                cout << "GT"<< j << "=" << getGTStringFromChar(sample_GT[j].GT[i]) << ";";
            }

            cout << endl;
        }
    }
};

/**
 * @class alt_format_df
 * @brief Data frame for formatted alternative allele data for samples.
 *
 * This class stores formatted alternative allele information for samples,
 * including variant IDs, sample IDs, alternative allele IDs, and associated
 * sample float, flag, string, and integer data.
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
     * @brief Initializes the formatted alternative allele data frame.
     *
     * Resizes and initializes internal vectors based on batch size and number of samples.
     *
     * @param ref Reference alt_format_df to copy names and sample info.
     * @param FORMAT Header element containing format information.
     * @param batch_size Batch size used for processing.
     */
    void init(alt_format_df ref, header_element FORMAT, long batch_size){
        numSample = ref.numSample;
        sampNames = ref.sampNames;
        int numAlt = batch_size*numSample*2;
        var_id.resize(numAlt, 0);
        samp_id.resize(numAlt, 0);
        alt_id.resize(numAlt, (char)0);
        
        int tmp = FORMAT.floats_alt;
        if(tmp>0){
            samp_Float tmpFloat;
            for(int i = 0; i<tmp; i++){
                tmpFloat.name = ref.samp_float[i].name;
                tmpFloat.numb = ref.samp_float[i].numb;
                tmpFloat.i_float.resize(numAlt, 0);
                samp_float.push_back(tmpFloat);
            }
            samp_float.resize(tmp);
        }
        tmp = FORMAT.ints_alt;
        if(tmp>0){
            samp_Int tmpInt;
            for(int i = 0; i<tmp; i++){
                tmpInt.name = ref.samp_int[i].name;
                tmpInt.numb = ref.samp_int[i].numb;
                tmpInt.i_int.resize(numAlt, 0);
                samp_int.push_back(tmpInt);
            }
            samp_int.resize(tmp);
        }
        tmp = FORMAT.strings_alt;
        if(tmp>0){
            samp_String tmpString;
            for(int i = 0; i<tmp; i++){
                tmpString.name = ref.samp_string[i].name;
                tmpString.numb = ref.samp_string[i].numb;
                tmpString.i_string.resize(numAlt, "\0");
                samp_string.push_back(tmpString);
            }
            samp_string.resize(tmp);
        }       
    }   

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

    /**
     * @brief Clones the formatted alternative allele data frame from a reference.
     *
     * This method is not currently used and needs to be updated.
     *
     * @param ref Reference alt_format_df to clone from.
     * @param FORMAT Header element containing format information.
     */
    void clone(alt_format_df ref, header_element FORMAT){
        int numAlt = FORMAT.alt_values;
        var_id.resize(numAlt, 0);
        samp_id.resize(numAlt, 0);
        alt_id.resize(numAlt, (char)0);
        numSample = ref.numSample;
        sampNames = ref.sampNames;

        int tmp = FORMAT.floats_alt;
        if(tmp>0){
            samp_Float tmpFloat;
            for(int i = 0; i<tmp; i++){
                tmpFloat.name = ref.samp_float[i].name;
                tmpFloat.numb = ref.samp_float[i].numb;
                tmpFloat.i_float.resize(ref.samp_float[i].i_float.size(), 0);
                samp_float.push_back(tmpFloat);
            }
            samp_float.resize(tmp);
        }
    
        tmp = FORMAT.ints_alt;

        if(tmp>0){
            samp_Int tmpInt;
            for(int i = 0; i<tmp; i++){
                tmpInt.name = ref.samp_int[i].name;
                tmpInt.numb = ref.samp_int[i].numb;
                tmpInt.i_int.resize(ref.samp_int[i].i_int.size(), 0);
                samp_int.push_back(tmpInt);
            }
            samp_int.resize(tmp);
        }
        tmp = FORMAT.strings_alt;
        if(tmp>0){
            samp_String tmpString;
            for(int i = 0; i<tmp; i++){
                tmpString.name = ref.samp_string[i].name;
                tmpString.numb = ref.samp_string[i].numb;
                tmpString.i_string.resize(ref.samp_string[i].i_string.size(), "\0");
                samp_string.push_back(tmpString);
            }
            samp_string.resize(tmp);
        }       
    }   
 
    /**
     * @brief Prints the formatted alternative allele data frame.
     *
     * Outputs the first n entries including variant ID, sample ID (resolved to name),
     * alternative allele ID, and associated sample float, int, string, and genotype data.
     *
     * @param n Number of entries to print.
     */
    void print(int n){
        cout << "VarID\tSampID\talt_id\tFloat\t\tInt\t\tStr\t\tGT" << endl;

        int iter = (n>samp_id.size()) ? samp_id.size() : n;

        for(int i=0; i<iter; i++){

            cout << var_id[i] << "\t";

            // Reverse lookup in sampNames map to get the string associated with samp_id
            for (const auto& pair : sampNames) {
                if (pair.second == samp_id[i]) {
                    cout << pair.first << "\t";
                    break;
                }
            }

            cout << (int)alt_id[i] << "\t";

            for(int j=0; j < samp_float.size(); j++){
                cout << samp_float[j].name << "=" << static_cast<float>(samp_float[j].i_float[i]) << ";";
            }
            cout << "\t";

            for(int j=0; j < samp_int.size(); j++){
                cout << samp_int[j].name << "=" << samp_int[j].i_int[i] << ";";
            }

            cout << "\t";
            for(int j=0; j < samp_string.size(); j++){
                cout << samp_string[j].name << "=" << samp_string[j].i_string[i] << ";";
            }

            cout << "\t";
            cout << getGTStringFromChar(sample_GT.GT[i]) << ";";
            cout << endl;
        }
    }
};



#endif