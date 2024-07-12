#ifndef VCF_STRUCTS_H
#define VCF_STRUCTS_H
#include <chrono>
#include <boost/algorithm/string.hpp>

using namespace std;

class info_flag
{
    public:
    vector<bool> i_flag;
    string name;
};

class info_string
{
    public:
    vector<string> i_string;
    string name;
};

class info_float
{
    public:
    vector<float> i_float;
    string name;
};

class info_int
{
    public:
    vector<int> i_int;
    string name;
};

class samp_Flag
{
    public:
    vector<bool> i_flag;
    string name;
    int numb;
};

class samp_String
{
    public:
    vector<string> i_string;
    string name;
    int numb;
};

class samp_Float
{
    public:
    vector<float> i_float;
    string name;
    int numb;
};

class samp_Int
{
    public:
    vector<int> i_int;
    string name;
    int numb;
};

class header_element
{
public:
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

};

class alt_columns_df
{
    public:
    vector<string> var_id;
    vector<int> alt_id;
    vector<string> alt;
    vector<info_float> alt_float;
    vector<info_flag> alt_flag; //non gestite per ora
    vector<info_string> alt_string;
    vector<info_int> alt_int;
    int numAlt;

    void clone(alt_columns_df ref, header_element INFO){
        int numAlt = INFO.alt_values;//sigsegv
        var_id.resize(numAlt, "\0");
        alt.resize(numAlt, "\0");
        alt_id.resize(numAlt, 0);
        int tmp = INFO.floats_alt;
        if(tmp>0){
            info_float tmpInfoFloat;
            for(int i = 0; i<tmp; i++){
                tmpInfoFloat.name = ref.alt_float[i].name;
                tmpInfoFloat.i_float.resize(ref.alt_float[i].i_float.size(), 0);
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
    
    void print(){
        cout << "VarID\tAltID\tAlt\tFloat\t\tInt\t\tStr" << endl;
        for(int i=0; i<numAlt; i++){
            cout << var_id[i] << "\t" << alt_id[i] << "\t" << alt[i] << "\t";
            for(int j=0; j < alt_float.size(); j++){
                cout << alt_float[j].name << "=" << alt_float[j].i_float[i] << ";";
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
    vector<string> var_id;
    vector<string> samp_id;
    vector<samp_Float> samp_float;
    vector<samp_Flag> samp_flag;
    vector<samp_String> samp_string;
    vector<samp_Int> samp_int;
    vector<string> sampNames;
    int numSample; //numero di sample per riga

   void print(){
        cout << "VarID\tSampID\tFloat\t\tInt\t\tStr" << endl;
        int iter = samp_id.size();
        for(int i=0; i<iter; i++){
            cout << var_id[i] << "\t" << samp_id[i] << "\t";
            for(int j=0; j < samp_float.size(); j++){
                cout << samp_float[j].name << "=" << samp_float[j].i_float[i] << ";";
            }
            cout << "\t";
            for(int j=0; j < samp_int.size(); j++){
                cout << samp_int[j].name << "=" << samp_int[j].i_int[i] << ";";
            }
            cout << "\t";
            for(int j=0; j < samp_string.size(); j++){
                cout << samp_string[j].name << "=" << samp_string[j].i_string[i] << ";";
            }
            cout << endl;
        }
    }
};

class alt_format_df //aka df4 in progress
{
    public:
    vector<string> var_id;
    vector<string> samp_id;
    vector<int> alt_id;
    vector<samp_Float> samp_float;
    vector<samp_Flag> samp_flag;
    vector<samp_String> samp_string;
    vector<samp_Int> samp_int;
    vector<string> sampNames;
    int numSample;

    void clone(alt_format_df ref, header_element FORMAT){
        int numAlt = FORMAT.alt_values;
        var_id.resize(numAlt, "\0");
        samp_id.resize(numAlt, "\0");
        alt_id.resize(numAlt, 0);
        numSample = ref.numSample;
        sampNames.resize(numSample, "\0");

        for(int i=0; i<numSample; i++){
            sampNames[i] = ref.sampNames[i];
        }

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
 
    void print(){ //in progress
        cout << "VarID\tSampID\talt_id\tFloat\t\tInt\t\tStr" << endl;
        int iter = samp_id.size();
        for(int i=0; i<iter; i++){
            cout << var_id[i] << "\t" << samp_id[i] << "\t" << alt_id[i] << "\t";
            for(int j=0; j < samp_float.size(); j++){
                cout << samp_float[j].name << "=" << samp_float[j].i_float[i] << ";";
            }
            cout << "\t";
            for(int j=0; j < samp_int.size(); j++){
                cout << samp_int[j].name << "=" << samp_int[j].i_int[i] << ";";
            }
            cout << "\t";
            for(int j=0; j < samp_string.size(); j++){
                cout << samp_string[j].name << "=" << samp_string[j].i_string[i] << ";";
            }
            cout << endl;
        }
    }
};

class var_columns_df
{
public:
    vector<long> var_number;
    vector<string> chrom;
    vector<long>pos;
    vector<string> id;
    vector<string> ref;
    vector<string> alt;
    vector<float> qual;
    vector<string> filter;
    vector<string> info;
    vector<info_float> in_float;
    vector<info_flag> in_flag;
    vector<info_string> in_string;
    vector<info_int> in_int;
    map<string,int> info_map1;
    void get_vcf_line_in_var_columns(char *line, long start, long end, long i, alt_columns_df* tmp_alt, int *tmp_num_alt)
    {
        //cout<<i<<"\t";
        bool find1 = false;
        long iter=0;
        int local_alt = 1;
        string tmp="\0";
        vector<string> tmp_split;
        vector<string> tmp_format_split;
        
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                chrom[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
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
                    (*tmp_alt).alt_id[(*tmp_num_alt)+y] = y;
                    (*tmp_alt).var_id[(*tmp_num_alt)+y] = id[i]; 
                }
            }else{
                tmp += line[start+iter]; // per ora le salvo tutte come un unico array of char, andrebbe cambiato per salvarle separatamente
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(strcmp(&tmp[0], ".")==0){
                    qual[i] = 0.0f;
                }else{
                    qual[i] = stof(tmp); // da cambiare, in futuro
                }
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                filter[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
                find1 = true;
                iter++;
                info[i] = tmp;
                vector<string> tmp_el;
                boost::split(tmp_el, tmp, boost::is_any_of(";")); //separo gli argomenti di info
                vector<string> tmp_elems;
                for(int ii=0; ii<tmp_el.size(); ii++){
                    boost::split(tmp_elems, tmp_el[ii], boost::is_any_of("=")); //separo id dell'info dal contenuto
                    bool find_info_type = false;
                    bool find_info_elem = false;
                    if(tmp_elems.size()==2){
                        while(!find_info_type){
                            if(info_map1[tmp_elems[0]]==1){
                                //Int - divisione
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
                            }else if(info_map1[tmp_elems[0]]==2){
                                //Float
                                
                                int el=0;
                                while(!find_info_elem){
                                    if(in_float[el].name == tmp_elems[0]){
                                        in_float[el].i_float[i] = stof(tmp_elems[1]);
                                        find_info_elem = true;
                                    } 
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==3){
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
                            }else if(info_map1[tmp_elems[0]]==4){
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
                            }else if(info_map1[tmp_elems[0]]==5){
                                //Float Alt
                                int el=0;
                                while(!find_info_elem){
                                    if((*tmp_alt).alt_float[el].name == tmp_elems[0]){
                                        boost::split(tmp_split, tmp_elems[1], boost::is_any_of(","));
                                        for(int y = 0; y<local_alt; y++){
                                            (*tmp_alt).alt_float[el].i_float[(*tmp_num_alt)+y] = stof(tmp_split[y]);
                                        }
                                        find_info_elem = true;
                                    }
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==6){
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
                                //cout<<"\n\nhere\n\n";
                                find_info_type = true;
                            }
                        }
                    }else{
                        if((info_map1[tmp_elems[0]]==0) && strcmp(&tmp_elems[0][0],"")){
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
                tmp += line[start+iter]; // anche qui andrebbero separate le info e creati dizionari con keys e values
                iter++;
                
            }
        }
        (*tmp_num_alt) = (*tmp_num_alt)+local_alt;
    }

    void get_vcf_line_in_var_columns_format(char *line, long start, long end, long i, alt_columns_df* tmp_alt, int *tmp_num_alt, sample_columns_df* sample, header_element* FORMAT, int *tmp_num_alt_format, alt_format_df* tmp_alt_format)
    {
        //cout<<i<<"\t";
        bool find1 = false;
        long iter=0;
        int local_alt = 1;
        string tmp="\0";
        vector<string> tmp_split;
        vector<string> tmp_format_split;
        vector<string> tmp_subSplit;
        
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                chrom[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
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
                    (*tmp_alt).alt_id[(*tmp_num_alt)+y] = y;
                    (*tmp_alt).var_id[(*tmp_num_alt)+y] = id[i];
                }
            }else{
                tmp += line[start+iter]; // per ora le salvo tutte come un unico array of char, andrebbe cambiato per salvarle separatamente
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(strcmp(&tmp[0], ".")==0){
                    qual[i] = 0.0f;
                }else{
                    qual[i] = stof(tmp); // da cambiare, in futuro
                }
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                filter[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
                find1 = true;
                iter++;
                info[i] = tmp;
                vector<string> tmp_el;
                boost::split(tmp_el, tmp, boost::is_any_of(";")); //separo gli argomenti di info
                vector<string> tmp_elems;
                for(int ii=0; ii<tmp_el.size(); ii++){
                    boost::split(tmp_elems, tmp_el[ii], boost::is_any_of("=")); //separo id dell'info dal contenuto
                    bool find_info_type = false;
                    bool find_info_elem = false;
                    if(tmp_elems.size()==2){
                        while(!find_info_type){
                            if(info_map1[tmp_elems[0]]==1){
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
                            }else if(info_map1[tmp_elems[0]]==2){
                                //Float                                
                                int el=0;
                                while(!find_info_elem){
                                    if(in_float[el].name == tmp_elems[0]){
                                        in_float[el].i_float[i] = stof(tmp_elems[1]);
                                        find_info_elem = true;
                                    } 
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==3){
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
                            }else if(info_map1[tmp_elems[0]]==4){
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
                            }else if(info_map1[tmp_elems[0]]==5){
                                //Float Alt
                                int el=0;
                                while(!find_info_elem){
                                    if((*tmp_alt).alt_float[el].name == tmp_elems[0]){
                                        boost::split(tmp_split, tmp_elems[1], boost::is_any_of(","));
                                        for(int y = 0; y<local_alt; y++){
                                            (*tmp_alt).alt_float[el].i_float[(*tmp_num_alt)+y] = stof(tmp_split[y]);
                                        }
                                        find_info_elem = true;
                                    }
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==6){
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
                                //cout<<"\n\nhere\n\n";
                                find_info_type = true;
                            }
                        }
                    }else{
                        if((info_map1[tmp_elems[0]]==0) && strcmp(&tmp_elems[0][0],"")){
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
                tmp += line[start+iter]; // anche qui andrebbero separate le info e creati dizionari con keys e values
                iter++;
                
            }
        }
        (*tmp_num_alt) = (*tmp_num_alt)+local_alt;
        tmp="\0";
        find1=false; 
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
                            if(info_map1[tmp_format_split[j]] == 8 || info_map1[tmp_format_split[j] + std::to_string(1)] == 8){
                                //String
                                (*sample).var_id[i*(*sample).numSample + samp] = id[i];
                                (*sample).samp_id[i*(*sample).numSample + samp] = (*sample).sampNames[samp];
                                int el = 0;
                                while(!find_elem){
                                    if(!(*sample).samp_string[el].name.compare(0, tmp_format_split[j].length(), tmp_format_split[j], 0, tmp_format_split[j].length())){
                                        if((*sample).samp_string[el].numb==1){
                                            //aggiorna solo la cella corrispondente
                                            (*sample).samp_string[el].i_string[i*(*sample).numSample + samp] = tmp_split[j];
                                        }else{
                                            //itera anche sulle liste che hanno il nome con in aggiunta il numero che saranno ad el + numero
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
                            }else if(info_map1[tmp_format_split[j]] == 9 || info_map1[tmp_format_split[j] + std::to_string(1)] == 9){
                                //Integer
                                (*sample).var_id[i*(*sample).numSample + samp] = id[i];
                                (*sample).samp_id[i*(*sample).numSample + samp] = (*sample).sampNames[samp];
                                int el = 0;
                                while(!find_elem){
                                    if(!(*sample).samp_int[el].name.compare(0, tmp_format_split[j].length(), tmp_format_split[j], 0, tmp_format_split[j].length())){
                                        if((*sample).samp_int[el].numb==1){
                                            //aggiorna solo la cella corrispondente
                                            (*sample).samp_int[el].i_int[i*(*sample).numSample + samp] = std::stoi(tmp_split[j]);
                                        }else{
                                            //itera anche sulle liste che hanno il nome con in aggiunta il numero che saranno ad el + numero
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
                            }else if(info_map1[tmp_format_split[j]] == 10 || info_map1[tmp_format_split[j] + std::to_string(1)] == 10){
                                //Float
                                (*sample).var_id[i*(*sample).numSample + samp] = id[i];
                                (*sample).samp_id[i*(*sample).numSample + samp] = (*sample).sampNames[samp];
                                int el = 0;
                                while(!find_elem){
                                    if(!(*sample).samp_float[el].name.compare(0, tmp_format_split[j].length(), tmp_format_split[j], 0, tmp_format_split[j].length())){
                                        if((*sample).samp_float[el].numb==1){
                                            //aggiorna solo la cella corrispondente
                                            (*sample).samp_float[el].i_float[i*(*sample).numSample + samp] = std::stof(tmp_split[j]);
                                        }else{
                                            //itera anche sulle liste che hanno il nome con in aggiunta il numero che saranno ad el + numero
                                            vector<string> tmp_sub;
                                            boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                            for(int i = 0; i<(*sample).samp_float[el].numb; i++){
                                                (*sample).samp_float[el+i].i_float[i*(*sample).numSample + samp] = std::stof(tmp_sub[i]);
                                            }
                                        }
                                        find_elem = true;
                                    }
                                    el++;
                                }
                                find_type = true;
                            }else if(info_map1[tmp_format_split[j]] == 11){
                                //String alt
                                boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                local_alt = tmp_sub.size();
                                int el = 0;
                                while(!find_elem){ //in progress da cambiare
                                    if(!(*tmp_alt_format).samp_string[el].name.compare(tmp_format_split[j])){
                                        //da fare su misura per alternativesp
                                        for(int y = 0; y<local_alt; y++){ //in progress da inizializzare local_alt
                                            (*tmp_alt_format).var_id[(*tmp_num_alt_format) + y] = id[i];
                                            (*tmp_alt_format).samp_id[(*tmp_num_alt_format) + y] = (*tmp_alt_format).sampNames[samp];
                                            (*tmp_alt_format).alt_id[(*tmp_num_alt_format) + y] = y;
                                            (*tmp_alt_format).samp_string[el].i_string[(*tmp_num_alt_format) + y] = tmp_sub[y];
                                        }
                                        find_elem = true;
                                        (*tmp_num_alt_format) = (*tmp_num_alt_format) + local_alt;
                                    }
                                    el++;
                                }
                                find_type = true;
                            }else if(info_map1[tmp_format_split[j]] == 12){
                                //int alt
                                boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                local_alt = tmp_sub.size();
                                int el = 0;
                                while(!find_elem){ //in progress da cambiare
                                    if(!(*tmp_alt_format).samp_int[el].name.compare(tmp_format_split[j])){
                                        //da fare su misura per alternatives
                                        for(int y = 0; y<local_alt; y++){ //in progress da inizializzare local_alt
                                            (*tmp_alt_format).var_id[(*tmp_num_alt_format) + y] = id[i];
                                            (*tmp_alt_format).samp_id[(*tmp_num_alt_format) + y] = (*tmp_alt_format).sampNames[samp];
                                            (*tmp_alt_format).alt_id[(*tmp_num_alt_format) + y] = y;
                                            (*tmp_alt_format).samp_int[el].i_int[(*tmp_num_alt_format) + y] = std::stoi(tmp_sub[y]);
                                        }
                                        find_elem = true;
                                        (*tmp_num_alt_format) = (*tmp_num_alt_format) + local_alt;
                                    }
                                    el++;
                                }
                                find_type = true;
                            }else if(info_map1[tmp_format_split[j]] == 13){
                                //float alt
                                boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                local_alt = tmp_sub.size();
                                int el = 0;
                                while(!find_elem){ //in progress da cambiare
                                    if(!(*tmp_alt_format).samp_float[el].name.compare(tmp_format_split[j])){
                                        //da fare su misura per alternatives
                                        for(int y = 0; y<local_alt; y++){ //in progress da inizializzare local_alt
                                            (*tmp_alt_format).var_id[(*tmp_num_alt_format) + y] = id[i];
                                            (*tmp_alt_format).samp_id[(*tmp_num_alt_format) + y] = (*tmp_alt_format).sampNames[samp];
                                            (*tmp_alt_format).alt_id[(*tmp_num_alt_format) + y] = y;
                                            (*tmp_alt_format).samp_float[el].i_float[(*tmp_num_alt_format) + y] = std::stof(tmp_sub[y]);
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

    void print_var_columns(long num_lines){
        for(long i=0; i<num_lines; i++){
            cout << "Var" << var_number[i] << ":\t";
            cout << chrom[i] << "\t";
            cout << to_string(pos[i]) << "\t";
            cout << id[i] << "\t";
            cout << ref[i] << "\t";
            cout << alt[i] << "\t";
            cout << to_string(qual[i]) << "\t";
            cout << filter[i] << "\t";
            cout << info[i] << "\t";
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

class var
{
public:
    long var_number;
	string chrom="\0"; // su quale cromosoma della reference
    long pos; // posizione nel chromosoma
    string id="\0"; // id della variation, spesso .
    string ref="\0"; // nucleotide della reference
    string alt="\0"; // possibili alternative (possono essere più di una)
    float qual; // quality score
    string filter="\0"; // filtro penso usato durante il sequenziamento
    string info="\0"; // info varie deducibili dall'header, potrebbe essere utile averle in collegamento
    string format="\0"; // formato dei samples, info variabiliti, mi dice come sono ordinate
    string samples="\0"; // una colonna per ogni sample in cui ognuno descrive i valori indicati in format
    void get_vcf_line(char *line, long start, long end)
    {
        bool find1 = false;
        long iter=0;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                chrom += line[start+iter];
                iter++;
            }
        }
        string tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                pos = stoul(tmp); // da cambiare, in futuro
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                id += line[start+iter];
                iter++;
            }
        }
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                ref += line[start+iter];
                iter++;
            }
        }
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                alt += line[start+iter]; // per ora le salvo tutte come un unico array of char, andrebbe cambiato per salvarle separatamente
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(strcmp(&tmp[0], ".")==0){
                    qual = 0.0;
                }else{
                    qual = stof(tmp); // da cambiare, in futuro
                }
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                filter += line[start+iter];
                iter++;
            }
        }
        find1=false;
        while(!find1 && (start+iter)<(end-1)){
            if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
                find1 = true;
                iter++;
            }else{
                info += line[start+iter]; // anche qui andrebbero separate le info e creati dizionari con keys e values
                iter++;
            }
        }
        /*find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
                find1 = true;
                iter++;
            }else{
                format += line[start+iter]; // vanno separati e soprattutti connessi ai valori nei samples (tipo se nel format ho ..:AF e nel sample ..:2, devo poter fare la query AF=2 in tutte le var)
                iter++;
            }
        }
        
        while((start+iter)<(end-1)){ //segmentation fault -> perchè non (start+iter)<(end-1)?
            samples += line[start+iter]; // sicuramente da dividere i vari sample e eventualmente creare un dizionario per ogni sample cosi da associare con il format
            iter++;
        }*/
        
    }
    void print_var()
    {
        cout << "Var" << var_number << ":\t";
        cout << chrom << "\t";
        cout << to_string(pos) << "\t";
        cout << id << "\t";
        cout << ref << "\t";
        cout << alt << "\t";
        cout << to_string(qual) << "\t";
        cout << filter << "\t";
        cout << info << endl;
        //cout << format << "\t";
        //cout << samples;
        //cout << endl;
    }
};

class vcf_parsed
{
public:
    int id;
    string filename;
    var *var_df;
    string header;
    header_element INFO;
    map<string,int> info_map; //Flag=0, Int=1, Float=2, String=3;
    header_element FORMAT; //aggiungere Filter e gli altri, magari con switch case
    char *filestring;
    int header_size=0;
    long filesize;
    long variants_size;
    long num_lines=0;
    long *new_lines_index;
    bool samplesON = false;
    var_columns_df var_columns;
    alt_columns_df alt_columns;
    sample_columns_df samp_columns;
    alt_format_df alt_sample;

    void get_filename(string path_to_filename){
        vector<string> line_el;
        boost::split(line_el, path_to_filename, boost::is_any_of("/"));
        filename = line_el[line_el.size()-1];
        //cout << "\nOnly file name: " << filename << endl;
        line_el.clear();
        boost::split(line_el, filename, boost::is_any_of("."));
        if(line_el[line_el.size()-1] == "gz")
        {
            //cout << "\nFile must be uncompressed!\n" << endl;
            exit(0);
        }else{
            //cout << "\nFile already uncompressed!\n" << endl;
        }
    }
    void get_file_size(string filename){
        filesize = filesystem::file_size(filename);
    }
    void get_header(ifstream *file){
        string line;
        //removing the header and storing it in vcf.header
        while (getline(*file, line) && line[0]=='#' && line[1]=='#'){
            header.append(line + '\n');
            header_size += line.length() + 1;
        }
        header_size += line.length() + 1;
        //cout << "\nheader char: " << to_string(header_size) << endl;
        variants_size = filesize - header_size; //ora che ho tolto l'header ho un file piu piccolo quindi una nuova size
        //cout<<"filesize: "<<filesize<<" variants_size: "<<variants_size<<endl;
    }
    void print_header(){
        //cout << "VCF header:\n" << header << endl;
    }
    void get_and_parse_header(ifstream *file){
        string line;
        vector<string> line_el; //all the characteristics together
        vector<string> line_el1; //each characteristic
        vector<string> line_el2; //keys and values
        //removing the header and storing it in vcf.header
        
        while (getline(*file, line) && line[0]=='#' && line[1]=='#'){
            header.append(line + '\n');
            header_size += line.length() + 1;
            bool Info = (line[2]=='I');
            bool Format = (line[2]=='F' && line[3]=='O');
            
            if(Info || Format){
                boost::split(line_el, line, boost::is_any_of("><"));
                boost::split(line_el1, line_el[1], boost::is_any_of(","));
                for(int i=0; i<3; i++){
                    boost::split(line_el2, line_el1[i], boost::is_any_of("="));
                    if(Info){
                        if(i==0) INFO.ID.push_back(line_el2[1]);
                        if(i==1) INFO.Number.push_back(line_el2[1]);
                        if(i==1 && line_el2[1] == "A") INFO.alt_values++;
                        if(i==2) INFO.Type.push_back(line_el2[1]);
                    }
                    if(Format){
                        if(i==0) FORMAT.ID.push_back(line_el2[1]);
                        if(i==1) FORMAT.Number.push_back(line_el2[1]);
                        if(i==1 && line_el2[1] == "A") FORMAT.alt_values++;
                        if(i==2) FORMAT.Type.push_back(line_el2[1]);
                    }
                }
            }
        }

        vector<string> tmp_split;
        boost::split(tmp_split, line, boost::is_any_of("\t"));
        if(tmp_split.size() > 9){
            samplesON = true;
            samp_columns.numSample = tmp_split.size() - 9;
            alt_sample.numSample = samp_columns.numSample;
            samp_columns.sampNames.resize(samp_columns.numSample, "\0");
            alt_sample.sampNames.resize(samp_columns.numSample, "\0");
            for(int i = 0; i < samp_columns.numSample; i++){
                samp_columns.sampNames[i] = tmp_split[9+i];
                alt_sample.sampNames[i] = tmp_split[9+i];
            }
        }else{
            samp_columns.numSample = 0;
        }
        
 
        // for(int i=0; i<INFO.ID.size(); i++){
        //    cout<<"INFO.ID["<<i<<"]: "<<INFO.ID[i]<<" | INFO.Number["<<i<<"]: "<<INFO.Number[i]<<" | INFO.Type["<<i<<"]: "<<INFO.Type[i]<<endl;
        // }
        INFO.total_values = INFO.ID.size();
        INFO.no_alt_values = INFO.total_values - INFO.alt_values;
        // cout<<"INFO.total_values: "<< INFO.total_values << " INFO.alt_values: "<<INFO.alt_values<<" INFO.no_alt_values: "<<INFO.no_alt_values<<endl;
        header_size += line.length() + 1;
        //cout << "\nheader char: " << to_string(header_size) << endl;
        variants_size = filesize - header_size; //ora che ho tolto l'header ho un file piu piccolo quindi una nuova size
        //cout<<"filesize: "<<filesize<<" variants_size: "<<variants_size<<endl;
    }   
    void allocate_filestring(){
        filestring = (char*)malloc(variants_size);
    }
    void find_new_lines_index(string w_filename, int num_threads){ //popola anche il filestring 
        new_lines_index = (long*)malloc(variants_size); //per ora ho esagerato con la dimensione (è come se permettessi tutti \n. Si puo ridurre (euristicamente), pero ipotizzarlo è meglio perche senno devo passare il file due volte solo per vedere dove iniziano le linee)
        new_lines_index[0] = 0; //il primo elemento lo metto a zero per indicare l'inizio della prima linea
        num_lines++;
        long tmp_num_lines[num_threads];
        
        
        auto before = chrono::system_clock::now();
        long batch_infile = (variants_size - 1 + num_threads)/num_threads; //numero di char che verrà processato da ogni thread        
#pragma omp parallel
        {
            int thr_ID = omp_get_thread_num();
            ifstream infile(w_filename); //apro lo stesso file con un ifstream in ogni thread, da considerare se ha senso creare prima un array di ifstream
            infile.seekg((header_size + thr_ID*batch_infile), ios::cur); //ogni thread parte a copiare inFile in filestring a una distanza di batch_size
            long start, end;
            start = thr_ID*batch_infile; // inizio del batch dello specifico thread
            end = start + batch_infile; // fine del batch
            
            tmp_num_lines[thr_ID] = 0;
            if(thr_ID==0){
                tmp_num_lines[0] = 1;
            } 

            for(long i=start; i<end && i<variants_size; i++){
                filestring[i] = infile.get();
                if(filestring[i]=='\n'){
                    tmp_num_lines[thr_ID] = tmp_num_lines[thr_ID] + 1;
                }
            }
        }
        while(filestring[variants_size-1]=='\n'){
            variants_size--;
        }
        filestring[variants_size] = '\n';
        variants_size++;
        auto after = chrono::system_clock::now();
        auto filestring_time = std::chrono::duration<double>(after - before).count();
        
        before = chrono::system_clock::now();
#pragma omp parallel
        {
            int thr_ID = omp_get_thread_num();
            long start, end;
            start = thr_ID*batch_infile; // inizio del batch dello specifico thread
            end = start + batch_infile; // fine del batch
            long startNLI = 1;
            if(thr_ID!=0){
                for(int i=0; i<thr_ID; i++){
                    startNLI = startNLI + tmp_num_lines[i];
                }
                startNLI--;
            }
            long lineCount = 0;
            for(long i=start; i<end && i<variants_size; i++){
                if(filestring[i]=='\n'&& filestring[i+1]!='\n'){
                    new_lines_index[startNLI+lineCount] = i+1; // PROBLEMA: funziona solo se l'ultimo char è uno /n, altrimenti si rompe
                    lineCount++;
                }
            }
        }
        num_lines = tmp_num_lines[0];
        for(int i=1; i<num_threads; i++){
            num_lines= num_lines + tmp_num_lines[i];
        }
        after = chrono::system_clock::now();
        auto f_new_lines = std::chrono::duration<double>(after - before).count(); 
        //cout << "\nFilestring time: " << filestring_time << " s " << "New lines time: " << f_new_lines << " s\n\n" << endl;
    }
    void create_sample_vectors(int num_threads){
        long batch_size = (num_lines-2+num_threads)/num_threads;
        samp_Flag samp_flag_tmp;
        samp_Float samp_float_tmp;
        samp_Int samp_int_tmp;
        samp_String samp_string_tmp;

        samp_Float samp_alt_float_tmp;
        samp_Int samp_alt_int_tmp;
        samp_String samp_alt_string_tmp;

        int numIter = FORMAT.ID.size();
        // cout << "numiter: " << numIter<<endl;
        for(int i = 0; i < numIter; i++){
            if(strcmp(&FORMAT.Number[i][0], "A") != 0){
                if(strcmp(&FORMAT.Number[i][0], "1")==0){
                    if(!strcmp(&FORMAT.Type[i][0], "String")){ 
                        samp_string_tmp.name = FORMAT.ID[i];
                        samp_columns.samp_string.push_back(samp_string_tmp);
                        samp_columns.samp_string.back().i_string.resize((num_lines-1)*samp_columns.numSample, "\0");
                        samp_columns.samp_string.back().numb = std::stoi(FORMAT.Number[i]);
                        info_map[FORMAT.ID[i]] = 8;
                        var_columns.info_map1[FORMAT.ID[i]] = 8;
                        FORMAT.strings++;
                    }else if(!strcmp(&FORMAT.Type[i][0], "Integer")){
                        samp_int_tmp.name = FORMAT.ID[i];
                        samp_columns.samp_int.push_back(samp_int_tmp);
                        samp_columns.samp_int.back().i_int.resize((num_lines-1)*samp_columns.numSample, 0);
                        samp_columns.samp_int.back().numb = std::stoi(FORMAT.Number[i]);
                        info_map[FORMAT.ID[i]] = 9;
                        var_columns.info_map1[FORMAT.ID[i]] = 9;
                        FORMAT.ints++;
                    }else if(!strcmp(&FORMAT.Type[i][0], "Float")){
                        samp_float_tmp.name = FORMAT.ID[i];
                        samp_columns.samp_float.push_back(samp_float_tmp);
                        samp_columns.samp_float.back().i_float.resize((num_lines-1)*samp_columns.numSample, 0);
                        samp_columns.samp_float.back().numb = std::stoi(FORMAT.Number[i]);
                        info_map[FORMAT.ID[i]] = 10;
                        var_columns.info_map1[FORMAT.ID[i]] = 10;
                        FORMAT.floats++;
                    }
                }else if(strcmp(&FORMAT.Number[i][0], "0")==0){
                    samp_flag_tmp.name = FORMAT.ID[i];
                    samp_columns.samp_flag.push_back(samp_flag_tmp);
                    samp_columns.samp_flag.back().i_flag.resize((num_lines-1)*samp_columns.numSample, 0);
                    samp_columns.samp_flag.back().numb = std::stoi(FORMAT.Number[i]);
                    info_map[FORMAT.ID[i]] = 11;
                    var_columns.info_map1[FORMAT.ID[i]] = 11;
                    FORMAT.flags++;
                }else{
                    if(!strcmp(&FORMAT.Type[i][0], "String")){
                        for(int j = 0; j < std::stoi(FORMAT.Number[i]); j++){
                            samp_string_tmp.name = FORMAT.ID[i] + std::to_string(j);
                            samp_columns.samp_string.push_back(samp_string_tmp);
                            samp_columns.samp_string.back().i_string.resize((num_lines-1)*samp_columns.numSample, "\0");
                            samp_columns.samp_string.back().numb = std::stoi(FORMAT.Number[i]);
                            info_map[FORMAT.ID[i]+std::to_string(j)] = 8;
                            var_columns.info_map1[FORMAT.ID[i]+std::to_string(j)] = 8;
                            FORMAT.strings++;
                        }
                    }else if(!strcmp(&FORMAT.Type[i][0], "Integer")){
                        for(int j = 0; j < std::stoi(FORMAT.Number[i]); j++){
                            samp_int_tmp.name = FORMAT.ID[i] + std::to_string(j);
                            samp_columns.samp_int.push_back(samp_int_tmp);
                            samp_columns.samp_int.back().i_int.resize((num_lines-1)*samp_columns.numSample, 0);
                            samp_columns.samp_int.back().numb = std::stoi(FORMAT.Number[i]);
                            info_map[FORMAT.ID[i]+std::to_string(j)] = 9;
                            var_columns.info_map1[FORMAT.ID[i]+std::to_string(j)] = 9;
                            FORMAT.ints++;
                        }
                    }else if(!strcmp(&FORMAT.Type[i][0], "Float")){
                        for(int j = 0; j < std::stoi(FORMAT.Number[i]); j++){
                            samp_float_tmp.name = FORMAT.ID[i] + std::to_string(j);
                            samp_columns.samp_float.push_back(samp_float_tmp);
                            samp_columns.samp_float.back().i_float.resize((num_lines-1)*samp_columns.numSample, 0);
                            samp_columns.samp_float.back().numb = std::stoi(FORMAT.Number[i]);
                            info_map[FORMAT.ID[i]+std::to_string(j)] = 10;
                            var_columns.info_map1[FORMAT.ID[i]+std::to_string(j)] = 10;
                            FORMAT.floats++;
                        }
                    }
                }
            }else{
                //In progress -> parte con alternatives
                if(!strcmp(&FORMAT.Type[i][0], "String")){ 
                    samp_alt_string_tmp.name = FORMAT.ID[i];
                    alt_sample.samp_string.push_back(samp_alt_string_tmp);
                    alt_sample.samp_string.back().i_string.resize(batch_size*samp_columns.numSample*2, "\0");
                    alt_sample.samp_string.back().numb = -1;
                    info_map[FORMAT.ID[i]] = 11;
                    var_columns.info_map1[FORMAT.ID[i]] = 11;
                    FORMAT.strings_alt++;
                }else if(!strcmp(&FORMAT.Type[i][0], "Integer")){
                    samp_alt_int_tmp.name = FORMAT.ID[i];
                    alt_sample.samp_int.push_back(samp_alt_int_tmp);
                    alt_sample.samp_int.back().i_int.resize(batch_size*samp_columns.numSample*2, 0);
                    alt_sample.samp_int.back().numb = -1;
                    info_map[FORMAT.ID[i]] = 12;
                    var_columns.info_map1[FORMAT.ID[i]] = 12;
                    FORMAT.ints_alt++;
                }else if(!strcmp(&FORMAT.Type[i][0], "Float")){
                    samp_alt_float_tmp.name = FORMAT.ID[i];
                    alt_sample.samp_float.push_back(samp_alt_float_tmp);
                    alt_sample.samp_float.back().i_float.resize(batch_size*samp_columns.numSample*2, 0);
                    alt_sample.samp_float.back().numb = -1;
                    info_map[FORMAT.ID[i]] = 13;
                    var_columns.info_map1[FORMAT.ID[i]] = 13;
                    FORMAT.floats_alt++;
                }
            }
        }
        
        // cout << "numiter: " << numIter<<endl;
        
        samp_columns.samp_flag.resize(FORMAT.flags);
        samp_columns.samp_int.resize(FORMAT.ints);
        samp_columns.samp_float.resize(FORMAT.floats);
        samp_columns.samp_string.resize(FORMAT.strings);
        samp_columns.var_id.resize((num_lines-1)*samp_columns.numSample, "\0");
        samp_columns.samp_id.resize((num_lines-1)*samp_columns.numSample, "\0");

        alt_sample.samp_flag.resize(FORMAT.flags_alt);
        alt_sample.samp_int.resize(FORMAT.ints_alt);
        alt_sample.samp_float.resize(FORMAT.floats_alt);
        alt_sample.samp_string.resize(FORMAT.strings_alt);
        alt_sample.var_id.resize((num_lines-1)* alt_sample.numSample, "\0");
        alt_sample.samp_id.resize((num_lines-1)*alt_sample.numSample, "\0");
    }
    void create_info_vectors(int num_threads){
        long batch_size = (num_lines-2+num_threads)/num_threads;
        info_flag info_flag_tmp;
        info_float info_float_tmp;
        info_int info_int_tmp;
        info_string info_string_tmp;
        
        info_float alt_float_tmp;
        info_int alt_int_tmp;
        info_string alt_string_tmp;
        for(int i=0; i<INFO.total_values; i++){
            if(strcmp(&INFO.Number[i][0], "A") == 0){ 
                if(strcmp(&INFO.Type[i][0], "Integer")==0){
                    INFO.ints_alt++;
                    alt_int_tmp.name = INFO.ID[i];
                    alt_int_tmp.i_int.resize(2*batch_size, 0);
                    alt_columns.alt_int.push_back(alt_int_tmp);
                    info_map[INFO.ID[i]] = 4;
                    var_columns.info_map1[INFO.ID[i]] = 4;
                }
                if(strcmp(&INFO.Type[i][0], "Float")==0){
                    INFO.floats_alt++;
                    alt_float_tmp.name = INFO.ID[i];
                    alt_float_tmp.i_float.resize(2*batch_size, 0);
                    alt_columns.alt_float.push_back(alt_float_tmp);
                    info_map[INFO.ID[i]] = 5;
                    var_columns.info_map1[INFO.ID[i]] = 5;
                }
                if(strcmp(&INFO.Type[i][0], "String")==0){
                    INFO.strings_alt++;
                    alt_string_tmp.name = INFO.ID[i];
                    alt_string_tmp.i_string.resize(2*batch_size, "\0");
                    alt_columns.alt_string.push_back(alt_string_tmp);
                    info_map[INFO.ID[i]] = 6;
                    var_columns.info_map1[INFO.ID[i]] = 6;
                }
                if(strcmp(&INFO.Type[i][0], "Flag")==0){ //Per ora non gestito
                    INFO.flags_alt++;
                    info_map[INFO.ID[i]] = 7;
                }
            }else if((strcmp(&INFO.Number[i][0], "1") == 0)||(strcmp(&INFO.Number[i][0], "0") == 0)){
                if(strcmp(&INFO.Type[i][0], "Integer")==0){
                    INFO.ints++;
                    info_int_tmp.name = INFO.ID[i];
                    info_int_tmp.i_int.resize(num_lines-1, 0);
                    var_columns.in_int.push_back(info_int_tmp);
                    info_map[INFO.ID[i]] = 1;
                    var_columns.info_map1[INFO.ID[i]] = 1;
                }
                if(strcmp(&INFO.Type[i][0], "Float")==0){
                    INFO.floats++;
                    info_float_tmp.name = INFO.ID[i];
                    info_float_tmp.i_float.resize(num_lines-1, 0);
                    var_columns.in_float.push_back(info_float_tmp);
                    info_map[INFO.ID[i]] = 2;
                    var_columns.info_map1[INFO.ID[i]] = 2;
                }
                if(strcmp(&INFO.Type[i][0], "String")==0){
                    INFO.strings++;
                    info_string_tmp.name = INFO.ID[i];
                    info_string_tmp.i_string.resize(num_lines-1, "\0");
                    var_columns.in_string.push_back(info_string_tmp);
                    info_map[INFO.ID[i]] = 3;
                    var_columns.info_map1[INFO.ID[i]] = 3;
                }
                if(strcmp(&INFO.Type[i][0], "Flag")==0){
                    INFO.flags++;
                    info_flag_tmp.name = INFO.ID[i];
                    info_flag_tmp.i_flag.resize(num_lines-1, 0);
                    var_columns.in_flag.push_back(info_flag_tmp);
                    info_map[INFO.ID[i]] = 0;
                    var_columns.info_map1[INFO.ID[i]] = 0;
                }else{
                    //in progress, se num > 1
                }
            }
        }
        
        var_columns.in_flag.resize(INFO.flags);
        var_columns.in_int.resize(INFO.ints);
        var_columns.in_float.resize(INFO.floats);
        var_columns.in_string.resize(INFO.strings);
        alt_columns.alt_int.resize(INFO.ints_alt);
        alt_columns.alt_float.resize(INFO.floats_alt);
        alt_columns.alt_string.resize(INFO.strings_alt);
    }
    void print_info_map(){
        for(const auto& element : info_map){
            cout<<element.first<<": "<<element.second<<endl;
        }
    }
    void print_info(){
        cout<<"Flags size: "<<var_columns.in_flag.size()<<endl;
        for(int i=0; i<var_columns.in_flag.size(); i++){
            cout<<var_columns.in_flag[i].name<<": ";
            for(int j=0; j<10; j++){
                cout<<var_columns.in_flag[i].i_flag[j]<<" ";
            }
            cout<<" size: "<<var_columns.in_flag[i].i_flag.size();
            cout<<endl;
        }
        cout<<endl;
        cout<<"Floats size: "<<var_columns.in_float.size()<<endl;
        for(int i=0; i<var_columns.in_float.size(); i++){
            cout<<var_columns.in_float[i].name<<": ";
            for(int j=0; j<10; j++){
                cout<<var_columns.in_float[i].i_float[j]<<" ";
            }
            cout<<" size: "<<var_columns.in_float[i].i_float.size();
            cout<<endl;
        }
        cout<<endl;
        cout<<"Strings size: "<<var_columns.in_string.size()<<endl;
        for(int i=0; i<var_columns.in_string.size(); i++){
            cout<<var_columns.in_string[i].name<<": ";
            for(int j=0; j<10; j++){
                cout<<var_columns.in_string[i].i_string[j]<<" ";
            }
            cout<<" size: "<<var_columns.in_string[i].i_string.size();
            cout<<endl;
        }
        cout<<endl;
        cout<<"Ints size: "<<var_columns.in_int.size()<<endl;
        for(int i=0; i<var_columns.in_int.size(); i++){
            cout<<var_columns.in_int[i].name<<": ";
            for(int j=0; j<10; j++){
                cout<<var_columns.in_int[i].i_int[j]<<" ";
            }
            cout<<" size: "<<var_columns.in_int[i].i_int.size();
            cout<<endl;
        }
        //cout<<endl;
        // cout<<"Flags: "<<endl;
        // for(int i=0; i<var_columns.in_flag.size(); i++){
        //     cout<<var_columns.in_flag[i].name<<": ";
        //     for(int j=0; j<var_columns.in_flag[i].i_flag.size(); j++){
        //         cout<<var_columns.in_flag[i].i_flag[j]<<" ";
        //     }
        //     cout<<endl;
        // }
        // cout<<endl;
        // cout<<"Floats: "<<endl;
        // for(int i=0; i<var_columns.in_float.size(); i++){
        //     cout<<var_columns.in_float[i].name<<": ";
        //     for(int j=0; j<var_columns.in_float[i].i_float.size(); j++){
        //         cout<<var_columns.in_float[i].i_float[j]<<" ";
        //     }
        //     cout<<endl;
        // }
        // cout<<endl;
        // cout<<"Strings: "<<endl;
        // for(int i=0; i<var_columns.in_string.size(); i++){
        //     cout<<var_columns.in_string[i].name<<": ";
        //     for(int j=0; j<var_columns.in_string[i].i_string.size(); j++){
        //         cout<<var_columns.in_string[i].i_string[j]<<" ";
        //     }
        //     cout<<endl;
        // }
        // cout<<endl;
        // cout<<"Ints: "<<endl;
        // for(int i=0; i<var_columns.in_int.size(); i++){
        //     cout<<var_columns.in_int[i].name<<": ";
        //     for(int j=0; j<var_columns.in_int[i].i_int.size(); j++){
        //         cout<<var_columns.in_int[i].i_int[j]<<" ";
        //     }
        //     cout<<endl;
        // }
        // cout<<endl;
    }
    void reserve_var_columns(){
        var_columns.var_number.resize(num_lines-1);
        var_columns.chrom.resize(num_lines-1);
        var_columns.id.resize(num_lines-1);
        var_columns.pos.resize(num_lines-1);
        var_columns.ref.resize(num_lines-1);
        var_columns.alt.resize(num_lines-1);
        var_columns.qual.resize(num_lines-1);
        var_columns.filter.resize(num_lines-1);
        var_columns.info.resize(num_lines-1);
        //cout<<"Finish resize!"<<endl;
    }
    void populate_var_columns(int num_threads){
        long batch_size = (num_lines-2+num_threads)/num_threads;
        alt_columns_df tmp_alt[num_threads];
        int tmp_num_alt[num_threads];
        //in progress
        alt_format_df tmp_alt_format[num_threads];
        int tmp_num_alt_format[num_threads];
        // cout<<"sqa"<<endl;
#pragma omp parallel
        {
            long start, end;
            int th_ID = omp_get_thread_num();
            //struttura temporanea del thread con alternatives.
            tmp_alt[th_ID].clone(alt_columns, INFO);
            // cout<<"sqae"<<endl;
            // tmp_alt_format[th_ID].clone(alt_sample, FORMAT); // si blocca qui
            // cout<<"sqar"<<endl;
            tmp_num_alt[th_ID] = 0;
            tmp_alt[th_ID].var_id.resize(batch_size*2, "\0");
            tmp_alt[th_ID].alt_id.resize(batch_size*2, 0);
            tmp_alt[th_ID].alt.resize(batch_size*2, "\0");

            // tmp_num_alt_format[th_ID] = 0;
            // tmp_alt_format[th_ID].var_id.resize(batch_size*2*samp_columns.numSample, "\0");
            // tmp_alt_format[th_ID].alt_id.resize(batch_size*2*samp_columns.numSample, 0);
            // tmp_alt_format[th_ID].samp_id.resize(batch_size*2*samp_columns.numSample, "\0");
            //cout << "\nThread id: "<<th_ID<<endl;
            start = th_ID*batch_size; // inizio del batch dello specifico thread
            end = start + batch_size; // fine del batch
            //cout<<"\nNum Lines: "<<num_lines-2<<endl;
            if(samplesON){ 
                // cout<<"sa"<<endl;
                tmp_alt_format[th_ID].clone(alt_sample, FORMAT);
                tmp_num_alt_format[th_ID] = 0;
                tmp_alt_format[th_ID].var_id.resize(batch_size*2*samp_columns.numSample, "\0");
                tmp_alt_format[th_ID].alt_id.resize(batch_size*2*samp_columns.numSample, 0);
                tmp_alt_format[th_ID].samp_id.resize(batch_size*2*samp_columns.numSample, "\0");
                for(long i=start; i<end && i<num_lines-1; i++){ 
                    var_columns.var_number[i] = i;
                    var_columns.get_vcf_line_in_var_columns_format(filestring, new_lines_index[i], new_lines_index[i+1], i, &(tmp_alt[th_ID]), &(tmp_num_alt[th_ID]), &samp_columns, &FORMAT, &(tmp_num_alt_format[th_ID]), &(tmp_alt_format[th_ID]));
                }

                tmp_alt[th_ID].var_id.resize(tmp_num_alt[th_ID]);
                tmp_alt[th_ID].alt_id.resize(tmp_num_alt[th_ID]);
                tmp_alt[th_ID].alt.resize(tmp_num_alt[th_ID]);
                for(int i=0; i<INFO.ints_alt; i++){
                    tmp_alt[th_ID].alt_int[i].i_int.resize(tmp_num_alt[th_ID]);
                }
                for(int i=0; i<INFO.floats_alt; i++){
                    tmp_alt[th_ID].alt_float[i].i_float.resize(tmp_num_alt[th_ID]);
                }
                for(int i=0; i<INFO.strings_alt; i++){
                    tmp_alt[th_ID].alt_string[i].i_string.resize(tmp_num_alt[th_ID]);
                }
                tmp_alt[th_ID].numAlt = tmp_num_alt[th_ID];

                tmp_alt_format[th_ID].var_id.resize(tmp_num_alt_format[th_ID]);
                tmp_alt_format[th_ID].alt_id.resize(tmp_num_alt_format[th_ID]);
                tmp_alt_format[th_ID].samp_id.resize(tmp_num_alt_format[th_ID]);
                for(int i=0; i<FORMAT.ints_alt; i++){
                   tmp_alt_format[th_ID].samp_int[i].i_int.resize(tmp_num_alt_format[th_ID]);
                }
                for(int i=0; i<FORMAT.floats_alt; i++){
                    tmp_alt_format[th_ID].samp_float[i].i_float.resize(tmp_num_alt_format[th_ID]);
                }
                for(int i=0; i<FORMAT.strings_alt; i++){
                    tmp_alt_format[th_ID].samp_string[i].i_string.resize(tmp_num_alt_format[th_ID]);
                }
                tmp_alt_format[th_ID].numSample = tmp_num_alt_format[th_ID]; 

            }else{
                // cout<<"no"<<endl;
                for(long i=start; i<end && i<num_lines-1; i++){
                    var_columns.var_number[i] = i;
                    var_columns.get_vcf_line_in_var_columns(filestring, new_lines_index[i], new_lines_index[i+1], i, &(tmp_alt[th_ID]), &(tmp_num_alt[th_ID]));
                }
                tmp_alt[th_ID].var_id.resize(tmp_num_alt[th_ID]);
                tmp_alt[th_ID].alt_id.resize(tmp_num_alt[th_ID]);
                tmp_alt[th_ID].alt.resize(tmp_num_alt[th_ID]);
                for(int i=0; i<INFO.ints_alt; i++){
                tmp_alt[th_ID].alt_int[i].i_int.resize(tmp_num_alt[th_ID]);
                }
                for(int i=0; i<INFO.floats_alt; i++){
                tmp_alt[th_ID].alt_float[i].i_float.resize(tmp_num_alt[th_ID]);
                }
                for(int i=0; i<INFO.strings_alt; i++){
                tmp_alt[th_ID].alt_string[i].i_string.resize(tmp_num_alt[th_ID]);
                }
                tmp_alt[th_ID].numAlt = tmp_num_alt[th_ID];
            }            
                       
        }
        int totAlt = 0;
        int totSampAlt = 0;
        // in progress Non troppo efficiente ma valido per ora
        for(int i=0; i<INFO.ints_alt; i++){
            alt_columns.alt_int[i].i_int.resize(0);
        }
        for(int i=0; i<INFO.floats_alt; i++){
            alt_columns.alt_float[i].i_float.resize(0);
        }
        for(int i=0; i<INFO.strings_alt; i++){
            alt_columns.alt_string[i].i_string.resize(0);
        }
        for(int i=0; i<num_threads; i++){
            alt_columns.var_id.insert(
                alt_columns.var_id.end(),
                tmp_alt[i].var_id.begin(),
                tmp_alt[i].var_id.end()
            );
            alt_columns.alt_id.insert(
                alt_columns.alt_id.end(),
                tmp_alt[i].alt_id.begin(),
                tmp_alt[i].alt_id.end()
            );
            alt_columns.alt.insert(
                alt_columns.alt.end(),
                tmp_alt[i].alt.begin(),
                tmp_alt[i].alt.end()
            );
            for(int j=0; j<INFO.ints_alt; j++){
                alt_columns.alt_int[j].i_int.insert(
                    alt_columns.alt_int[j].i_int.end(), 
                    tmp_alt[i].alt_int[j].i_int.begin(), 
                    tmp_alt[i].alt_int[j].i_int.end());
            }
            for(int j=0; j<INFO.floats_alt; j++){
                alt_columns.alt_float[j].i_float.insert(
                    alt_columns.alt_float[j].i_float.end(), 
                    tmp_alt[i].alt_float[j].i_float.begin(), 
                    tmp_alt[i].alt_float[j].i_float.end());
            }
            for(int j=0; j<INFO.strings_alt; j++){
                alt_columns.alt_string[j].i_string.insert(
                    alt_columns.alt_string[j].i_string.end(), 
                    tmp_alt[i].alt_string[j].i_string.begin(), 
                    tmp_alt[i].alt_string[j].i_string.end());
            }
            totAlt+=tmp_num_alt[i];
        }
        if(samplesON){
        // in progress Non troppo efficiente ma valido per ora
        for(int i=0; i<FORMAT.ints_alt; i++){
            alt_sample.samp_int[i].i_int.resize(0);
        }
        for(int i=0; i<FORMAT.floats_alt; i++){
            alt_sample.samp_float[i].i_float.resize(0);
        }
        for(int i=0; i<FORMAT.strings_alt; i++){
            alt_sample.samp_string[i].i_string.resize(0);
        }
        alt_sample.var_id.resize(0);
        alt_sample.alt_id.resize(0);
        alt_sample.samp_id.resize(0);    
        for(int i=0; i<num_threads; i++){
            alt_sample.var_id.insert(
                alt_sample.var_id.end(),
                tmp_alt_format[i].var_id.begin(),
                tmp_alt_format[i].var_id.end()
            );
            alt_sample.alt_id.insert(
                alt_sample.alt_id.end(),
                tmp_alt_format[i].alt_id.begin(),
                tmp_alt_format[i].alt_id.end()
            );
            alt_sample.samp_id.insert(
                alt_sample.samp_id.end(),
                tmp_alt_format[i].samp_id.begin(),
                tmp_alt_format[i].samp_id.end()
            );
            for(int j=0; j<FORMAT.ints_alt; j++){
                alt_sample.samp_int[j].i_int.insert(
                    alt_sample.samp_int[j].i_int.end(), 
                    tmp_alt_format[i].samp_int[j].i_int.begin(), 
                    tmp_alt_format[i].samp_int[j].i_int.end());
            }
            for(int j=0; j<FORMAT.floats_alt; j++){
                alt_sample.samp_float[j].i_float.insert(
                    alt_sample.samp_float[j].i_float.end(), 
                    tmp_alt_format[i].samp_float[j].i_float.begin(), 
                    tmp_alt_format[i].samp_float[j].i_float.end());
            }
            for(int j=0; j<FORMAT.strings_alt; j++){
                alt_sample.samp_string[j].i_string.insert(
                    alt_sample.samp_string[j].i_string.end(), 
                    tmp_alt_format[i].samp_string[j].i_string.begin(), 
                    tmp_alt_format[i].samp_string[j].i_string.end());
            }
            totSampAlt+=tmp_num_alt_format[i];
        }

        alt_columns.numAlt = totAlt;
        alt_columns.var_id.resize(totAlt);
        alt_columns.alt_id.resize(totAlt);
        alt_columns.alt.resize(totAlt);
        for(int j=0; j<INFO.ints_alt; j++){
            alt_columns.alt_int[j].i_int.resize(totAlt);
        }
        for(int j=0; j<INFO.floats_alt; j++){
            alt_columns.alt_float[j].i_float.resize(totAlt);
        }
        for(int j=0; j<INFO.strings_alt; j++){
            alt_columns.alt_string[j].i_string.resize(totAlt);
        }

        alt_sample.numSample = totSampAlt;
        alt_sample.var_id.resize(totSampAlt);
        alt_sample.samp_id.resize(totSampAlt);
        alt_sample.alt_id.resize(totSampAlt);
        for(int j=0; j<FORMAT.ints_alt; j++){
            alt_sample.samp_int[j].i_int.resize(totSampAlt);
        }
        for(int j=0; j<FORMAT.floats_alt; j++){
            alt_sample.samp_float[j].i_float.resize(totSampAlt);
        }
        for(int j=0; j<FORMAT.strings_alt; j++){
            alt_sample.samp_string[j].i_string.resize(totSampAlt);
        }
        }
        
    }

    void populate_var_struct(int num_threads){
        
        auto before = chrono::system_clock::now();
        
        var_df = (var*)calloc((num_lines-1), sizeof(var)); // allocating var_df    
        //cout << "\nBegin tmp: \n" <<"newlines: "<<num_lines<<" num threads: "<<num_threads<<endl;
        long batch_size = (num_lines-2+num_threads)/num_threads; //numero di lines che verrà processato da ogni thread
        
        //cout << "\nBatch size: "<<batch_size<<endl;
        auto after = chrono::system_clock::now();
        auto pre_pragma = std::chrono::duration<double>(after - before).count();
        //cout << "Pre_pragma: " << pre_pragma << " s" << endl;

#pragma omp parallel
        {
            long start, end;
            int th_ID = omp_get_thread_num();
            //cout << "\nThread id: "<<th_ID<<endl;

            start = th_ID*batch_size; // inizio del batch dello specifico thread
            end = start + batch_size; // fine del batch

            //cout << "\nstart: " << start << " end: " << end << endl; //segfault

            for(long i=start; i<end && i<num_lines-1; i++){ //start e end mi dicono l'intervallo di linee che eseguirà ogni thread, quindi la i rappresenta l'iesima linea (var) che inizia a new_lines_index[i] e finisce a new_lines_index[i+1] (escluso)
                var_df[i].get_vcf_line(filestring, new_lines_index[i], new_lines_index[i+1]); //qui traduco da char* a var structure la specifica line
                var_df[i].var_number = i; // dato aggiuntivo come se fosse un altro ID interno 
            }

        }
    }

};

#endif