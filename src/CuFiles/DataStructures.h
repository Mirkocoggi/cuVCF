#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#include <chrono>
#include <boost/algorithm/string.hpp>
#include <cuda_runtime.h>     
#include <cuda_fp16.h>  
#include <iostream>
#include <map>
#include <vector>

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

#endif