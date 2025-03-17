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

class alt_columns_df //DF2
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

        GTMap[".|."] = static_cast<char>(254);
        GTMap["./."] = static_cast<char>(255);
        
    }

    std::string getGTStringFromChar(char gtChar) const {
        for (const auto& pair : GTMap) {
            if (pair.second == gtChar) {
                return pair.first;
            }
        }
        return "Not found";
    }

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

            //TODO - il segfault è qui dentro
            for(int j=0; j < sample_GT[0].numb; j++){
                cout << "GT"<< j << "=" << getGTStringFromChar(sample_GT[j].GT[i]) << ";";
            }

            cout << endl;
        }
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

    std::string getGTStringFromChar(char gtChar) const {
        for (const auto& pair : GTMap) {
            if (pair.second == gtChar) {
                return pair.first;
            }
        }
        return "Not found";
    }

    //Not used, need to be updated
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
 
    void print(int n){ //TODO
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

class var_columns_df //DF1
{
public:
    vector<unsigned int> var_number;
    std::map<std::string, unsigned char> chrom_map;
    vector<char> chrom;
    vector<unsigned int>pos;
    vector<string> id;
    vector<string> ref;
    vector<__half> qual;
    std::map<std::string, char> filter_map;
    vector<char> filter;
    vector<info_float> in_float;
    vector<info_flag> in_flag;
    vector<info_string> in_string;
    vector<info_int> in_int;
    map<string,int> info_map1;

    void print(long num_lines) {
        std::cout << "VarID\tChrom\tPos\tID\tRef\tFilter\tFlag\t\tInt\t\tFloat\t\tString" << std::endl;
        long iter = (num_lines > static_cast<long>(id.size())) ? id.size() : num_lines;
        
        for (long i = 0; i < iter; i++) {
            // Stampa VarID (presumibilmente var_number)
            if (i < static_cast<long>(var_number.size()))
                std::cout << var_number[i] << "\t";
            else
                std::cout << "nan\t";
            
            // Stampa Chrom: costruiamo la chiave da un singolo carattere
            if (i < static_cast<long>(chrom.size())) {
                std::string chromKey(1, chrom[i]);
                auto itChrom = chrom_map.find(chromKey);
                if (itChrom != chrom_map.end()) {
                    std::cout << itChrom->first << "\t";
                } else {
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
            
            // Stampa Filter: costruiamo la chiave dal carattere in filter
            if (i < static_cast<long>(filter.size())) {
                std::string filterKey(1, filter[i]);
                auto itFilter = filter_map.find(filterKey);
                if (itFilter != filter_map.end()) {
                    std::cout << itFilter->first << "\t";
                } else {
                    std::cout << "nan\t";
                }
            } else {
                std::cout << "nan\t";
            }
            
            /* Stampa i campi dei vettori info: Flag, Int, Float e String
            // Si assume che ogni vettore info_xxx abbia la stessa dimensione (almeno per indice i)
            // e che il campo 'name' contenga il nome della colonna.
            
            // Flag
            for (size_t j = 0; j < in_flag.size(); j++) {
                if (i < static_cast<long>(in_flag[j].i_flag.size()))
                    std::cout << in_flag[j].name << ": " << in_flag[j].i_flag[i] << " ";
                else
                    std::cout << in_flag[j].name << ": nan ";
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
            */
            std::cout << std::endl;
        }
    }
};


#endif