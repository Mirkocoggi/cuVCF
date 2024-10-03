#ifndef VCF_PARCED_H
#define VCF_PARCED_H
#include "VCFparser_mt_col_struct.h"
#include "VCF_var.h"
#include "VCF_var_columns_df.h"
#include <chrono>
#include <boost/algorithm/string.hpp>
#include <OpenEXR/half.h>

class vcf_parsed
{
public:
    int id;
    string filename;
    var *var_df;
    string header;
    header_element INFO;
    map<string,int> info_map; //Flag=0, Int=1, Float=2, String=3;
    header_element FORMAT;
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

    void get_filename(string path_to_filename){ //TODO è ancora necessario il controllo sul .gz?
        vector<string> line_el;
        boost::split(line_el, path_to_filename, boost::is_any_of("/"));
        filename = line_el[line_el.size()-1];
        //cout << "\nOnly file name: " << filename << endl;
        line_el.clear();
        boost::split(line_el, filename, boost::is_any_of("."));
        if(line_el[line_el.size()-1] == "gz")
        {
            cout << "\nFile must be uncompressed!\n" << endl;
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
        variants_size = filesize - header_size; // New size without the header
        //cout<<"filesize: "<<filesize<<" variants_size: "<<variants_size<<endl;
    }
    
    void print_header(){
        cout << "VCF header:\n" << header << endl;
    }
    
    void get_and_parse_header(ifstream *file){
        string line;
        vector<string> line_el;     //all the characteristics together
        vector<string> line_el1;    //each characteristic
        vector<string> line_el2;    //keys and values
        // removing the header and storing it in vcf.header
        
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
                        if(i==0){
                            // TODO
                            //if(line_el2[1] == "GT"){
                            //    FORMAT.hasGT = true;
                            //    boost::split(line_el2, line_el1[1], boost::is_any_of("="));
                            //    FORMAT.numGT = line_el2[1];
                            //    i+=3;
                            //}else{
                                FORMAT.ID.push_back(line_el2[1]);
                            //}
                        } 
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

            for(int i = 0; i < samp_columns.numSample; i++){
                samp_columns.sampNames.insert(std::make_pair(tmp_split[9+i], i));
                alt_sample.sampNames.insert(std::make_pair(tmp_split[9+i], i));
            }

        }else{
            samp_columns.numSample = 0;
        }
        
        INFO.total_values = INFO.ID.size();
        INFO.no_alt_values = INFO.total_values - INFO.alt_values;

        header_size += line.length() + 1;

        variants_size = filesize - header_size; // New size without the header
    }   
    
    void allocate_filestring(){
        filestring = (char*)malloc(variants_size);
    }
  
    // This function identifies the indices of new line characters in the file and populates the `new_lines_index` array.
    // It also fills the 'filestring' variable (if applicable).
    void find_new_lines_index(string w_filename, int num_threads){
        // Allocate memory for the `new_lines_index` array. The size is exaggerated (assuming every character is a new line).
        new_lines_index = (long*)malloc(variants_size);
        new_lines_index[0] = 0; // The first element is set to 0, indicating the start of the first line.
        num_lines++; // Increment the line counter to account for the first line.
        long tmp_num_lines[num_threads]; // Temporary array to store the number of lines found by each thread.
        
        
        auto before = chrono::system_clock::now();
        long batch_infile = (variants_size - 1 + num_threads)/num_threads; // Number of characters each thread will process        
#pragma omp parallel
        {
            int thr_ID = omp_get_thread_num();
            ifstream infile(w_filename); // Open the same file with an ifstream in each thread
            infile.seekg((header_size + thr_ID*batch_infile), ios::cur); // Move each thread's file pointer to its starting position
            long start, end;
            start = thr_ID*batch_infile; // Starting position of the current thread’s batch
            end = start + batch_infile; // Ending position of the batch for the current thread
            
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
            start = thr_ID*batch_infile; // Starting position of the current thread’s batch
            end = start + batch_infile; // Ending position of the batch for the current thread
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
                    new_lines_index[startNLI+lineCount] = i+1;
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

        for(int i = 0; i < numIter; i++){
            if(strcmp(&FORMAT.Number[i][0], "A") != 0){
                // Without Alternatives
                if(strcmp(&FORMAT.Number[i][0], "1")==0){ 
                    //Number = 1
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
                    //Number = 0; so it's a flag
                    samp_flag_tmp.name = FORMAT.ID[i];
                    samp_columns.samp_flag.push_back(samp_flag_tmp);
                    samp_columns.samp_flag.back().i_flag.resize((num_lines-1)*samp_columns.numSample, 0);
                    samp_columns.samp_flag.back().numb = std::stoi(FORMAT.Number[i]);
                    info_map[FORMAT.ID[i]] = 11;
                    var_columns.info_map1[FORMAT.ID[i]] = 11;
                    FORMAT.flags++;
                }else{ 
                    //Number > 1
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
                //Alternatives
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
                
        samp_columns.samp_flag.resize(FORMAT.flags);
        samp_columns.samp_int.resize(FORMAT.ints);
        samp_columns.samp_float.resize(FORMAT.floats);
        samp_columns.samp_string.resize(FORMAT.strings);
        samp_columns.var_id.resize((num_lines-1)*samp_columns.numSample, "\0");
        samp_columns.samp_id.resize((num_lines-1)*samp_columns.numSample, static_cast<unsigned short>(0));

        alt_sample.samp_flag.resize(FORMAT.flags_alt);
        alt_sample.samp_int.resize(FORMAT.ints_alt);
        alt_sample.samp_float.resize(FORMAT.floats_alt);
        alt_sample.samp_string.resize(FORMAT.strings_alt);
        alt_sample.var_id.resize((num_lines-1)* alt_sample.numSample, "\0");
        alt_sample.samp_id.resize((num_lines-1)*alt_sample.numSample, static_cast<unsigned short>(0));

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
                // Alternatives
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
                // Without Alternatives and number = 1 or a flag
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
                    //in progress, se num > 1 TODO
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
    }
    
    void populate_var_columns(int num_threads){
        long batch_size = (num_lines-2+num_threads)/num_threads;
        
        alt_columns_df tmp_alt[num_threads];
        int tmp_num_alt[num_threads];
        
        alt_format_df tmp_alt_format[num_threads];
        int tmp_num_alt_format[num_threads];
        
#pragma omp parallel
        {
            long start, end;
            int th_ID = omp_get_thread_num();
            // Temporary structure of the thread with alternatives.
            tmp_alt[th_ID].clone(alt_columns, INFO);
            
            tmp_num_alt[th_ID] = 0;
            tmp_alt[th_ID].var_id.resize(batch_size*2, "\0");
            tmp_alt[th_ID].alt_id.resize(batch_size*2, 0);
            tmp_alt[th_ID].alt.resize(batch_size*2, "\0");
            
            tmp_num_alt_format[th_ID] = 0;
            tmp_alt_format[th_ID].var_id.resize(batch_size*2*samp_columns.numSample, "\0");
            tmp_alt_format[th_ID].alt_id.resize(batch_size*2*samp_columns.numSample, 0);
            tmp_alt_format[th_ID].samp_id.resize(batch_size*2*samp_columns.numSample, static_cast<unsigned short>(0));
            
            start = th_ID*batch_size; // Starting point of the thread's batch
            end = start + batch_size; // Ending point of the thread's batch
            
            if(samplesON){
                // There are samples in the dataset
                tmp_alt_format[th_ID].clone(alt_sample, FORMAT);
                tmp_num_alt_format[th_ID] = 0;
                tmp_alt_format[th_ID].var_id.resize(batch_size*2*samp_columns.numSample, "\0");
                tmp_alt_format[th_ID].alt_id.resize(batch_size*2*samp_columns.numSample, 0);
                tmp_alt_format[th_ID].samp_id.resize(batch_size*2*samp_columns.numSample, static_cast<unsigned short>(0));
                
                // For each line in the batch
                for(long i=start; i<end && i<num_lines-1; i++){ 
                    var_columns.var_number[i] = i;
                    var_columns.get_vcf_line_in_var_columns_format(filestring, new_lines_index[i], new_lines_index[i+1], i, &(tmp_alt[th_ID]), &(tmp_num_alt[th_ID]), &samp_columns, &FORMAT, &(tmp_num_alt_format[th_ID]), &(tmp_alt_format[th_ID]));
                }
                
                tmp_alt[th_ID].var_id.resize(tmp_num_alt[th_ID]);
                tmp_alt[th_ID].alt_id.resize(tmp_num_alt[th_ID]);
                tmp_alt[th_ID].alt.resize(tmp_num_alt[th_ID]);

                // For each integer variable
                for(int i=0; i<INFO.ints_alt; i++){
                    tmp_alt[th_ID].alt_int[i].i_int.resize(tmp_num_alt[th_ID]);
                }
                // For each float variable
                for(int i=0; i<INFO.floats_alt; i++){
                    tmp_alt[th_ID].alt_float[i].i_float.resize(tmp_num_alt[th_ID]);
                }
                // For each string variable
                for(int i=0; i<INFO.strings_alt; i++){
                    tmp_alt[th_ID].alt_string[i].i_string.resize(tmp_num_alt[th_ID]);
                }

                tmp_alt[th_ID].numAlt = tmp_num_alt[th_ID];
                tmp_alt_format[th_ID].var_id.resize(tmp_num_alt_format[th_ID]);
                tmp_alt_format[th_ID].alt_id.resize(tmp_num_alt_format[th_ID]);
                tmp_alt_format[th_ID].samp_id.resize(tmp_num_alt_format[th_ID]);

                // For each integer variable
                for(int i=0; i<FORMAT.ints_alt; i++){
                   tmp_alt_format[th_ID].samp_int[i].i_int.resize(tmp_num_alt_format[th_ID]);
                }
                // For each float variable
                for(int i=0; i<FORMAT.floats_alt; i++){
                    tmp_alt_format[th_ID].samp_float[i].i_float.resize(tmp_num_alt_format[th_ID]);
                }
                // For each string variable
                for(int i=0; i<FORMAT.strings_alt; i++){
                    tmp_alt_format[th_ID].samp_string[i].i_string.resize(tmp_num_alt_format[th_ID]);
                }

                tmp_alt_format[th_ID].numSample = tmp_num_alt_format[th_ID]; 

            }else{
                // There aren't samples in the dataset
                for(long i=start; i<end && i<num_lines-1; i++){
                    var_columns.var_number[i] = i;
                    var_columns.get_vcf_line_in_var_columns(filestring, new_lines_index[i], new_lines_index[i+1], i, &(tmp_alt[th_ID]), &(tmp_num_alt[th_ID]));
                }

                tmp_alt[th_ID].var_id.resize(tmp_num_alt[th_ID]);
                tmp_alt[th_ID].alt_id.resize(tmp_num_alt[th_ID]);
                tmp_alt[th_ID].alt.resize(tmp_num_alt[th_ID]);
                
                // For each integer variable
                for(int i=0; i<INFO.ints_alt; i++){
                tmp_alt[th_ID].alt_int[i].i_int.resize(tmp_num_alt[th_ID]);
                }
                // For each float variable
                for(int i=0; i<INFO.floats_alt; i++){
                tmp_alt[th_ID].alt_float[i].i_float.resize(tmp_num_alt[th_ID]);
                }
                // For each string variable
                for(int i=0; i<INFO.strings_alt; i++){
                tmp_alt[th_ID].alt_string[i].i_string.resize(tmp_num_alt[th_ID]);
                }

                tmp_alt[th_ID].numAlt = tmp_num_alt[th_ID];
            }         
        }

    //Here finish the parallel part and we merge the threads results

        int totAlt = 0;
        int totSampAlt = 0;
        // TODO Can be more efficient
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
            // TODO Can be more efficient
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

        long batch_size = (num_lines-2+num_threads)/num_threads; // number of lines passed to each thread
        
        auto after = chrono::system_clock::now();
        auto pre_pragma = std::chrono::duration<double>(after - before).count();

#pragma omp parallel
        {
            long start, end;
            int th_ID = omp_get_thread_num();

            start = th_ID*batch_size;
            end = start + batch_size;
            
            // start and end tells me the interval executed by each thread,
            //so the 'i' represents the i_th line (var) that starts at new_lines_index[i] and ends at new_lines_index[i+1] (excluded)
            for(long i=start; i<end && i<num_lines-1; i++){
                var_df[i].get_vcf_line(filestring, new_lines_index[i], new_lines_index[i+1]); //translation from char* to var structure of the line
                var_df[i].var_number = i; // dato aggiuntivo come se fosse un altro ID interno TODO, è necessario??
                                          // Il var_df alla fine di tutto ha una copia di tutti i dati dentro ma quali servono??
            }

        }
    }

};


#endif