#ifndef VCF_STRUCTS_H
#define VCF_STRUCTS_H

#include <chrono>
#include <boost/algorithm/string.hpp>
#include <cuda_runtime.h>     
#include <cuda_fp16.h>  
#include <iostream>
#include <map>

using namespace std;

struct KernelParams {
    char *line;
    unsigned int *var_number;
    unsigned int *pos;
    __half *qual;
    __half *in_float;
    bool *in_flag;
    int *in_int;
    char* float_name;
    char* flag_name;
    char* int_name;
    unsigned int *new_lines_index;
    unsigned int *samp_var_id;
    unsigned short *samp_id;
    __half *samp_float;
    bool *samp_flag;
    int *samp_int;
    char* samp_float_name;
    char* samp_flag_name;
    char* samp_int_name;
    int* samp_float_numb;
    int* samp_flag_numb;
    int* samp_int_numb;
    char *sample_GT;
    int numSample;
    unsigned int numLines;
    int numGT;
};

struct info_flag
{
    vector<uint8_t> i_flag;
    string name;
};

struct info_flag_d
{
    bool* i_flag;
    char* name;
};

struct info_string
{
    vector<string> i_string;
    string name;
};

struct info_float
{
    vector<__half> i_float;
    string name;
};

struct info_float_d
{
    __half *i_float;
    char *name;
};

struct info_int
{
    vector<int> i_int;
    string name;
};

struct info_int_d
{
    int *i_int;
    char *name;
};

struct samp_Flag
{
    vector<uint8_t> i_flag;
    string name;
    int numb;
};

struct samp_Flag_d
{   
    int *tot;
    bool *i_flag;
    char *name;
    int *numb;
};

struct samp_String
{
    vector<string> i_string;
    string name;
    int numb;
};

struct samp_Float
{
    vector<__half> i_float;
    string name;
    int numb;
};

struct samp_Float_d
{   
    int *tot;
    __half *i_float;
    char *name;
    int *numb;
};

struct samp_Int
{
    vector<int> i_int;
    string name;
    int numb;
};

struct samp_Int_d
{
    int *tot;
    int *i_int;
    char *name;
    int *numb;
};

struct samp_GT 
{
    vector<char> GT;
    int numb;    
};

struct samp_GT_d
{
    char *GT;
    int *numb;    
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

            for(int j=0; j < sample_GT.size(); j++){
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
#endif