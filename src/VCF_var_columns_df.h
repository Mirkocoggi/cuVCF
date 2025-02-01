#ifndef VCF_VARCOLUMNS_H
#define VCF_VARCOLUMNS_H
#include "VCFparser_mt_col_struct.h"
#include "VCF_parsed.h"
#include "VCF_var.h"
#include <chrono>
#include <boost/algorithm/string.hpp>
#include <Imath/half.h>

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

class var_columns_df
{
public:
    vector<unsigned int> var_number;
    std::map<std::string, char> chrom_map;
    vector<char> chrom; // TODO si riempie in modo sequenziale, ha gli stessi valori di var_number
    vector<unsigned int>pos;
    vector<string> id;
    vector<string> ref;
    vector<half> qual;
    std::map<std::string, char> filter_map;
    vector<char> filter;
    vector<info_float> in_float;
    vector<info_flag> in_flag;
    vector<info_string> in_string;
    vector<info_int> in_int;
    map<string,int> info_map1;
    
    void get_vcf_line_in_var_columns(char *line, long start, long end, long i, alt_columns_df* tmp_alt, int *tmp_num_alt)
    { 
        bool find1 = false;
        long iter=0;
        int local_alt = 1;
        string tmp="\0";
        vector<string> tmp_split;
        vector<string> tmp_format_split;
        
        //Chromosome
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(chrom_map.find(tmp) == chrom_map.end()){
                    chrom_map.insert(std::make_pair(tmp, (char)chrom_map.size()));
                }
                chrom[i] = chrom_map[tmp];
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Position
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                pos[i] = stoul(tmp);
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //ID
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                id[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Reference
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                ref[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Alternative
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                boost::split(tmp_split, tmp, boost::is_any_of(","));
                local_alt = tmp_split.size();
                for(int y = 0; y<local_alt; y++){
                    (*tmp_alt).alt[(*tmp_num_alt)+y] = tmp_split[y];
                    (*tmp_alt).alt_id[(*tmp_num_alt)+y] = (char)y;
                    (*tmp_alt).var_id[(*tmp_num_alt)+y] = var_number[i];
                }
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Quality
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(strcmp(&tmp[0], ".")==0){
                    qual[i] = 0.0f;
                }else{
                    try{
                        qual[i] = (half)stof(tmp);
                    }catch (const std::exception& e){
                        qual[i] = 0;
                    }
                }
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        
        //Filter
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(filter_map.find(tmp) == filter_map.end()){
                    filter_map.insert(std::make_pair(tmp, (char)filter_map.size()));
                }
                filter[i] = filter_map[tmp];
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Info
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
                find1 = true;
                iter++;
                vector<string> tmp_el;
                boost::split(tmp_el, tmp, boost::is_any_of(";")); //Info arguments separation
                vector<string> tmp_elems;
                for(int ii=0; ii<tmp_el.size(); ii++){
                    boost::split(tmp_elems, tmp_el[ii], boost::is_any_of("=")); //info_id separation from contents
                    bool find_info_type = false;
                    bool find_info_elem = false;
                    if(tmp_elems.size()==2){
                        while(!find_info_type){
                            if(info_map1[tmp_elems[0]]==INT){
                                //Int
                                bool isAlt = false;
                                int el=0;
                                while(!find_info_elem){
                                    if(in_int[el].name == tmp_elems[0]){
                                        in_int[el].i_int[i] = stoi(tmp_elems[1]);
                                        find_info_elem = true;
                                    }
                                    el++; 
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==FLOAT){
                                //Float
                                int el=0;
                                while(!find_info_elem){
                                    if(in_float[el].name == tmp_elems[0]){
                                        try{
                                            in_float[el].i_float[i] = (half)stof(tmp_elems[1]);
                                        }catch (const std::exception& e){
                                            in_float[el].i_float[i] = 0;
                                        }
                                        find_info_elem = true;
                                    } 
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==STRING){
                                //String
                                int el=0;
                                while(!find_info_elem){
                                    if(in_string[el].name == tmp_elems[0]){
                                        in_string[el].i_string[i] = tmp_elems[1];
                                        find_info_elem = true;
                                    } 
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==INT_ALT){
                                //Int Alt
                                int el=0;
                                while(!find_info_elem){
                                    if((*tmp_alt).alt_int[el].name == tmp_elems[0]){
                                        boost::split(tmp_split, tmp_elems[1], boost::is_any_of(","));
                                        for(int y = 0; y<local_alt; y++){
                                            (*tmp_alt).alt_int[el].i_int[(*tmp_num_alt)+y] = stoi(tmp_split[y]);
                                        }
                                        find_info_elem = true;
                                    }
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==FLOAT_ALT){
                                //Float Alt
                                int el=0;
                                while(!find_info_elem){
                                    if((*tmp_alt).alt_float[el].name == tmp_elems[0]){
                                        boost::split(tmp_split, tmp_elems[1], boost::is_any_of(","));
                                        
                                        for(int y = 0; y<local_alt; y++){
                                            try{
                                                (*tmp_alt).alt_float[el].i_float[(*tmp_num_alt)+y] = (half)stof(tmp_split[y]);
                                            }catch (const std::exception& e){
                                                (*tmp_alt).alt_float[el].i_float[(*tmp_num_alt)+y] = 0;
                                            }
                                        }
                                        find_info_elem = true;
                                    }
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==STRING_ALT){
                                //String Alt
                                int el=0;
                                while(!find_info_elem){
                                    if((*tmp_alt).alt_string[el].name == tmp_elems[0]){
                                        boost::split(tmp_split, tmp_elems[1], boost::is_any_of(","));
                                        for(int y = 0; y<local_alt; y++){
                                            (*tmp_alt).alt_string[el].i_string[(*tmp_num_alt)+y] = tmp_split[y];
                                        }
                                        find_info_elem = true;
                                    }
                                    el++;
                                }
                                find_info_type = true;                            
                            }else{
                                find_info_type = true;
                            }
                        }
                    }else{
                        if((info_map1[tmp_elems[0]]==FLAG) && strcmp(&tmp_elems[0][0],"")){
                            //Flag
                            int el=0;
                            while(!find_info_elem){
                                if(in_flag[el].name == tmp_elems[0]){
                                    in_flag[el].i_flag[i] = 1;
                                    find_info_elem = true;
                                } 
                                el++;
                            }
                                
                        }
                        tmp_elems.clear();
                    }
                }
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        (*tmp_num_alt) = (*tmp_num_alt)+local_alt;
    }

    void get_vcf_line_in_var_columns_format(char *line, long start, long end, long i, alt_columns_df* tmp_alt, int *tmp_num_alt, sample_columns_df* sample, header_element* FORMAT, int *tmp_num_alt_format, alt_format_df* tmp_alt_format)
    {
        bool find1 = false;
        long iter=0;
        int local_alt = 1;
        string tmp="\0";
        vector<string> tmp_split;
        vector<string> tmp_format_split;
        vector<string> tmp_subSplit;

        //Chromosome
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(chrom_map.find(tmp) == chrom_map.end()){
                    chrom_map.insert(std::make_pair(tmp, (char)chrom_map.size()));
                }
                chrom[i] = chrom_map[tmp];
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Position
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                pos[i] = stoul(tmp);
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //ID
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                id[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Reference
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                ref[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Alternatives
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                boost::split(tmp_split, tmp, boost::is_any_of(","));
                local_alt = tmp_split.size();
                for(int y = 0; y<local_alt; y++){
                    (*tmp_alt).alt[(*tmp_num_alt)+y] = tmp_split[y];
                    (*tmp_alt).alt_id[(*tmp_num_alt)+y] = (char)y;
                    (*tmp_alt).var_id[(*tmp_num_alt)+y] = var_number[i];
                }
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Quality
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(strcmp(&tmp[0], ".")==0){
                    qual[i] = (half)0.0f;
                }else{
                    try{
                        qual[i] = (half)stof(tmp);
                    }catch (const std::exception& e){
                        qual[i] = 0;
                    }
                }
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Filter
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(filter_map.find(tmp) == filter_map.end()){
                    filter_map.insert(std::make_pair(tmp, (char)filter_map.size()));
                }
                filter[i] = filter_map[tmp];
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Info
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
                find1 = true;
                iter++;
                vector<string> tmp_el;
                boost::split(tmp_el, tmp, boost::is_any_of(";")); //Info arguments separation
                vector<string> tmp_elems;
                for(int ii=0; ii<tmp_el.size(); ii++){
                    boost::split(tmp_elems, tmp_el[ii], boost::is_any_of("=")); //info_id separation from contents
                    bool find_info_type = false;
                    bool find_info_elem = false;
                    if(tmp_elems.size()==2){
                        while(!find_info_type){
                            if(info_map1[tmp_elems[0]]==INT){
                                //Int
                                bool isAlt = false;
                                int el=0;
                                while(!find_info_elem){
                                    if(in_int[el].name == tmp_elems[0]){
                                        in_int[el].i_int[i] = stoi(tmp_elems[1]);
                                        find_info_elem = true;
                                    }
                                    el++; 
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==FLOAT){
                                //Float                  
                                int el=0;
                                while(!find_info_elem){
                                    if(in_float[el].name == tmp_elems[0]){
                                        try{
                                            in_float[el].i_float[i] = (half)stof(tmp_elems[1]);
                                        }catch (const std::exception& e){
                                            in_float[el].i_float[i] = 0;
                                        }
                                        find_info_elem = true;
                                    } 
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==STRING){
                                //String                                
                                int el=0;
                                while(!find_info_elem){
                                    if(in_string[el].name == tmp_elems[0]){
                                        in_string[el].i_string[i] = tmp_elems[1];
                                        find_info_elem = true;
                                    } 
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==INT_ALT){
                                //Int Alternatives
                                int el=0;
                                while(!find_info_elem){
                                    if((*tmp_alt).alt_int[el].name == tmp_elems[0]){
                                        boost::split(tmp_split, tmp_elems[1], boost::is_any_of(","));
                                        for(int y = 0; y<local_alt; y++){
                                            (*tmp_alt).alt_int[el].i_int[(*tmp_num_alt)+y] = stoi(tmp_split[y]);
                                        }
                                        find_info_elem = true;
                                    }
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==FLOAT_ALT){
                                //Float Alt
                                int el=0;
                                while(!find_info_elem){
                                    if((*tmp_alt).alt_float[el].name == tmp_elems[0]){
                                        boost::split(tmp_split, tmp_elems[1], boost::is_any_of(","));
                                        for(int y = 0; y<local_alt; y++){
                                            try{
                                                (*tmp_alt).alt_float[el].i_float[(*tmp_num_alt)+y] = (half)stof(tmp_split[y]);
                                            }catch (const std::exception& e){
                                                (*tmp_alt).alt_float[el].i_float[(*tmp_num_alt)+y] = 0;
                                            }
                                        }
                                        find_info_elem = true;
                                    }
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==STRING_ALT){
                                //String Alt
                                int el=0;
                                while(!find_info_elem){
                                    if((*tmp_alt).alt_string[el].name == tmp_elems[0]){
                                        boost::split(tmp_split, tmp_elems[1], boost::is_any_of(","));
                                        for(int y = 0; y<local_alt; y++){
                                            (*tmp_alt).alt_string[el].i_string[(*tmp_num_alt)+y] = tmp_split[y];
                                        }
                                        find_info_elem = true;
                                    }
                                    el++;
                                }
                                find_info_type = true;   
                            }else{
                                find_info_type = true;
                            }
                        }
                    }else{
                        if((info_map1[tmp_elems[0]]==FLAG) && strcmp(&tmp_elems[0][0],"")){
                            //Flag
                            int el=0;
                            while(!find_info_elem){
                                if(in_flag[el].name == tmp_elems[0]){
                                    in_flag[el].i_flag[i] = 1;
                                    find_info_elem = true;
                                } 
                                el++;
                            }
                        }
                        tmp_elems.clear();
                    }
                }
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        (*tmp_num_alt) = (*tmp_num_alt)+local_alt;

        //Format decomposition
        tmp="\0";
        find1=false;

        //Format's template
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                boost::split(tmp_format_split, tmp, boost::is_any_of(":"));
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        
        //(*sample).print();

        //TODO - Da gestire il GT qui dentro
        int samp;
        for(samp = 0; samp < (*sample).numSample; samp++){
            tmp="\0";
            find1=false;
            while(!find1){
                if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
                    find1 = true;
                    iter++;
                    boost::split(tmp_split, tmp, boost::is_any_of(":"));
                    vector<string> tmp_sub;
                    for(int j = 0; j < tmp_split.size(); j++){
                        bool find_type = false;
                        bool find_elem = false;
                        while(!find_type){
                            if(!strcmp(tmp_format_split[j].c_str(), "GT")){
                                // TODO - da considerare come rifare la struttura visto che se numb > 1 servono + array
                                (*sample).var_id[i*(*sample).numSample + samp] = var_number[i];
                                
                                (*sample).samp_id[i*(*sample).numSample + samp] =  static_cast<unsigned short>(samp);
                                if((*sample).sample_GT.size() > 1){ //sample_columns_df* sample, alt_format_df* tmp_alt_format
                                    tmp_sub;
                                    boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                    for(int k=0; k < (*sample).sample_GT[0].numb; k++){ 
                                        (*sample).sample_GT[k].GT[i*(*sample).numSample + samp] = (*sample).GTMap[tmp_sub[k]];
                                    }
                                }else if((*sample).sample_GT.size() == 1){
                                    (*sample).sample_GT[0].GT[i*(*sample).numSample + samp] = (*sample).GTMap[tmp_split[j]];
                                }else{
                                    boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                    local_alt = tmp_sub.size();
                                    for(int y = 0; y<local_alt; y++){
                                        //Fill a tuple for each alternatives
                                        (*tmp_alt_format).var_id[(*tmp_num_alt_format) + y] = var_number[i];
                                        (*tmp_alt_format).samp_id[(*tmp_num_alt_format) + y] = samp;
                                        (*tmp_alt_format).alt_id[(*tmp_num_alt_format) + y] = (char)y;
                                        (*tmp_alt_format).sample_GT.GT[(*tmp_num_alt_format) + y] = (*tmp_alt_format).GTMap[tmp_sub[y]];
                                    }
                                    (*tmp_num_alt_format) = (*tmp_num_alt_format) + local_alt;
                                }
                                find_type = true;
                            }else if(info_map1[tmp_format_split[j]] == STRING_FORMAT || info_map1[tmp_format_split[j] + std::to_string(1)] == STRING_FORMAT){
                                //String - deterministic
                                (*sample).var_id[i*(*sample).numSample + samp] = var_number[i];
                                (*sample).samp_id[i*(*sample).numSample + samp] =  static_cast<unsigned short>(samp);
                                int el = 0;
                                
                                while(!find_elem){
                                    if(!(*sample).samp_string[el].name.compare(0, tmp_format_split[j].length(), tmp_format_split[j], 0, tmp_format_split[j].length())){
                                        if((*sample).samp_string[el].numb==1){ //String with numb = 1
                                            //Update the corresponing cell
                                            (*sample).samp_string[el].i_string[i*(*sample).numSample + samp] = tmp_split[j];
                                        }else{
                                            //String with numb > 1 (separated by commas)
                                            //Iterate over the alternatives (lists with the same name + ascending number, e.g., el1, el2, ...)
                                            //referred with 'el + number'
                                            vector<string> tmp_sub;
                                            boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                            for(int k = 0; k<(*sample).samp_string[el].numb; k++){
                                                (*sample).samp_string[el+k].i_string[i*(*sample).numSample + samp] = tmp_sub[k];
                                            }
                                        }
                                        find_elem = true;
                                    }
                                    el++;
                                }
                                find_type = true;
                            }else if(info_map1[tmp_format_split[j]] == INT_FORMAT || info_map1[tmp_format_split[j] + std::to_string(1)] == INT_FORMAT){
                                //Integer - deterministic
                                (*sample).var_id[i*(*sample).numSample + samp] = var_number[i];
                                (*sample).samp_id[i*(*sample).numSample + samp] = samp;
                                int el = 0;
                                while(!find_elem){
                                    if(!(*sample).samp_int[el].name.compare(0, tmp_format_split[j].length(), tmp_format_split[j], 0, tmp_format_split[j].length())){
                                        if((*sample).samp_int[el].numb==1){
                                            //Integer with numb = 1
                                            (*sample).samp_int[el].i_int[i*(*sample).numSample + samp] = std::stoi(tmp_split[j]);
                                        }else{
                                            //Integer with numb > 1
                                            vector<string> tmp_sub;
                                            boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                            for(int i = 0; i<(*sample).samp_int[el].numb; i++){
                                                (*sample).samp_int[el+i].i_int[i*(*sample).numSample + samp] = std::stoi(tmp_sub[i]);
                                            }
                                        }
                                        find_elem = true;
                                    }
                                    el++;
                                }
                                find_type = true;
                            }else if(info_map1[tmp_format_split[j]] == FLOAT_FORMAT || info_map1[tmp_format_split[j] + std::to_string(1)] == FLOAT_FORMAT){
                                //Float - deterministic
                                (*sample).var_id[i*(*sample).numSample + samp] = var_number[i];
                                (*sample).samp_id[i*(*sample).numSample + samp] = samp;
                                int el = 0;
                                while(!find_elem){
                                    if(!(*sample).samp_float[el].name.compare(0, tmp_format_split[j].length(), tmp_format_split[j], 0, tmp_format_split[j].length())){
                                        if((*sample).samp_float[el].numb==1){
                                            //Float with numb = 1
                                            try{ //Check if format is 10E-40 to crop to 0
                                                (*sample).samp_float[el].i_float[i*(*sample).numSample + samp] = (half)std::stof(tmp_split[j]);
                                            }catch (const std::exception& e){
                                                (*sample).samp_float[el].i_float[i*(*sample).numSample + samp] = 0;
                                            }
                                        }else{
                                            //Float with numb > 1
                                            vector<string> tmp_sub;
                                            boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                            for(int i = 0; i<(*sample).samp_float[el].numb; i++){
                                                try{ //Check if format is 10E-40 to crop to 0
                                                    (*sample).samp_float[el+i].i_float[i*(*sample).numSample + samp] = (half)std::stof(tmp_sub[i]);
                                                }catch (const std::exception& e){
                                                    (*sample).samp_float[el+i].i_float[i*(*sample).numSample + samp] = 0;
                                                }
                                            }
                                        }
                                        find_elem = true;
                                    }
                                    el++;
                                }
                                find_type = true;
                            }else if(info_map1[tmp_format_split[j]] == STRING_FORMAT_ALT){
                                //String alternatives
                                // TODO - Da gestire i GT con alternatives
                                boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                local_alt = tmp_sub.size();
                                int el = 0;
                                while(!find_elem){
                                    //Search the corresponding element
                                    if(!(*tmp_alt_format).samp_string[el].name.compare(tmp_format_split[j])){
                                        for(int y = 0; y<local_alt; y++){
                                            //Fill a tuple for each alternatives
                                            (*tmp_alt_format).var_id[(*tmp_num_alt_format) + y] = var_number[i];
                                            (*tmp_alt_format).samp_id[(*tmp_num_alt_format) + y] = samp;
                                            (*tmp_alt_format).alt_id[(*tmp_num_alt_format) + y] = (char)y;
                                            (*tmp_alt_format).samp_string[el].i_string[(*tmp_num_alt_format) + y] = tmp_sub[y];
                                        }
                                        find_elem = true;
                                        (*tmp_num_alt_format) = (*tmp_num_alt_format) + local_alt;
                                    }
                                    el++;
                                }
                                find_type = true;
                            }else if(info_map1[tmp_format_split[j]] == INT_FORMAT_ALT){
                                //Integer alternatives
                                boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                local_alt = tmp_sub.size();
                                int el = 0;
                                while(!find_elem){
                                    //Search the corresponding element
                                    if(!(*tmp_alt_format).samp_int[el].name.compare(tmp_format_split[j])){
                                        //Fill a tuple for each alternatives
                                        for(int y = 0; y<local_alt; y++){
                                            (*tmp_alt_format).var_id[(*tmp_num_alt_format) + y] = var_number[i];
                                            (*tmp_alt_format).samp_id[(*tmp_num_alt_format) + y] = samp;
                                            (*tmp_alt_format).alt_id[(*tmp_num_alt_format) + y] = (char)y;
                                            (*tmp_alt_format).samp_int[el].i_int[(*tmp_num_alt_format) + y] = std::stoi(tmp_sub[y]);
                                        }
                                        find_elem = true;
                                        (*tmp_num_alt_format) = (*tmp_num_alt_format) + local_alt;
                                    }
                                    el++;
                                }
                                find_type = true;
                            }else if(info_map1[tmp_format_split[j]] == FLOAT_FORMAT_ALT){
                                //Float alternatives
                                boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                local_alt = tmp_sub.size();
                                int el = 0;
                                while(!find_elem){ 
                                    //Search the corresponding element
                                    if(!(*tmp_alt_format).samp_float[el].name.compare(tmp_format_split[j])){
                                        //Fill a tuple for each alternatives
                                        for(int y = 0; y<local_alt; y++){
                                            (*tmp_alt_format).var_id[(*tmp_num_alt_format) + y] = var_number[i];
                                            (*tmp_alt_format).samp_id[(*tmp_num_alt_format) + y] = samp;
                                            (*tmp_alt_format).alt_id[(*tmp_num_alt_format) + y] = (char)y;
                                            try{
                                                (*tmp_alt_format).samp_float[el].i_float[(*tmp_num_alt_format) + y] = (half)std::stof(tmp_sub[y]);
                                            }catch (const std::exception& e){
                                                (*tmp_alt_format).samp_float[el].i_float[(*tmp_num_alt_format) + y] = 0;
                                            }
                                        }
                                        find_elem = true;
                                        (*tmp_num_alt_format) = (*tmp_num_alt_format) + local_alt;
                                    }
                                    el++;
                                }
                                find_type = true;
                            }
                        }
                    }
                }else{
                    tmp += line[start+iter];
                    iter++;
                }
            }
        }
        
    }

    void print(long num_lines){
        int iter = (num_lines>var_number.size()) ? var_number.size() : num_lines;
        for(long i=0; i<num_lines; i++){
            cout << "Var" << var_number[i] << ":\t";
            cout << chrom_map.find(std::string(1, chrom[i]))->first << "\t";
            cout << to_string(pos[i]) << "\t";
            cout << id[i] << "\t";
            cout << ref[i] << "\t";
            cout << filter_map.find(std::string(1, filter[i]))->first << "\t";

            cout << "AAAAAAAA 1.1" << endl;

            for(int j=0; j<in_flag.size(); j++){
                cout<<in_flag[j].name<<": "<<in_flag[j].i_flag[i]<<", ";
            }
            for(int j=0; j<in_int.size(); j++){
                cout<<in_int[j].name<<": "<<in_int[j].i_int[i]<<", ";
            }
            for(int j=0; j<in_float.size(); j++){
                cout<<in_float[j].name<<": "<<in_float[j].i_float[i]<<", ";
            }
            for(int j=0; j<in_string.size(); j++){
                cout<<in_string[j].name<<": "<<in_string[j].i_string[i]<<", ";
            }
            cout<<endl;

        }
    }    
};

#endif