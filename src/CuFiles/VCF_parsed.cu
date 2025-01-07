#ifndef VCF_PARCED_CU_H
#define VCF_PARCED_CU_H
#include <cuda_runtime.h>     
#include <cuda_fp16.h>  
#include <thrust/device_ptr.h> 
#include <thrust/sort.h>
#include <boost/algorithm/string.hpp>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <sys/wait.h>
#include <unistd.h>
#include <map>
#include <omp.h> 
#include "VCFparser_mt_col_struct_CU.h"
#include "VCF_CUDA_implementation.cu"
#include "VCF_var_columns_df_CU.h"

using namespace std;

class vcf_parsed
{
public:
    int id;
    string filename;
    string path_to_filename;
    string header;
    header_element INFO;
    map<string,int> info_map; //Flag=0, Int=1, Float=2, String=3;
    header_element FORMAT;
    char *filestring;
    char *d_filestring; 
    int header_size=0;
    long filesize;
    long variants_size;
    long num_lines=0;
    unsigned int *new_lines_index;
    unsigned int *d_new_lines_index;
    bool samplesON = false;
    bool hasDetSamples = false;
    var_columns_df var_columns;
    alt_columns_df alt_columns;
    sample_columns_df samp_columns;
    alt_format_df alt_sample;

    unsigned int *d_VC_var_number;
    unsigned int *d_VC_pos;
    __half *d_VC_qual;
    info_float_d *d_VC_in_float;
    info_flag_d *d_VC_in_flag;
    info_int_d *d_VC_in_int;

    unsigned int *d_SC_var_id;
    unsigned short *d_SC_samp_id;
    samp_Float_d *d_SC_samp_float;
    samp_Flag_d *d_SC_samp_flag;
    samp_Int_d *d_SC_samp_int;
    samp_GT_d *d_SC_sample_GT;

    void unzip_gz_file(char* vcf_filename) {
        // Check the extension ".gz"
        if (strcmp(vcf_filename + strlen(vcf_filename) - 3, ".gz") == 0) {
            pid_t pid = fork();
            if (pid == 0) {
                // Child proces
                execlp("gzip", "gzip", "-df", vcf_filename, nullptr);
                // If execlp fails
                cout<< "ERROR: Failed to execute gzip command" << std::endl;
                exit(EXIT_FAILURE);
            } else if (pid > 0) {
                // Parent proces waits for the child proces
                int status;
                waitpid(pid, &status, 0);
                if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
                    // Remove ".gz" from the filename
                    char* mutable_vcf_filename = const_cast<char*>(vcf_filename);
                    mutable_vcf_filename[strlen(vcf_filename) - 3] = '\0';
                } else {
                    cout<< "ERROR: cannot unzip file" << std::endl;
                }
            } else {
                // Error in the fork() call
                cout<< "ERROR: Failed to fork process" << std::endl;
            }
        }
    }

    void run(char* vcf_filename, int num_threadss){
        string filename = vcf_filename; 
        string line;
        vcf_parsed vcf;

        // TODO - definire parametri da usare per lanciare su GPU
            //lanciare una query per capire il numero di risorse della scheda e da li settare i parametro
        // Variables to hold device information
        size_t globalMemory = 0;        // Total global memory
        size_t sharedMemory = 0;       // Shared memory per block
        size_t constantMemory = 0;     // Constant memory
        size_t textureAlignment = 0;   // Texture alignment
        int maxThreadsPerBlock = 0;    // Maximum threads per block
        int threadsDim[3] = {0};       // Maximum threads per block dimensions
        int gridDim[3] = {0};          // Maximum grid dimensions

        // Query device properties
        cudaDeviceProp prop;
        cudaError_t err = cudaGetDeviceProperties(&prop, 0); // Query the first (and only) device

        if (err == cudaSuccess) {
            globalMemory = prop.totalGlobalMem;
            sharedMemory = prop.sharedMemPerBlock;
            constantMemory = prop.totalConstMem;
            textureAlignment = prop.textureAlignment;
            maxThreadsPerBlock = prop.maxThreadsPerBlock;

            // Threads and grid dimensions
            threadsDim[0] = prop.maxThreadsDim[0];
            threadsDim[1] = prop.maxThreadsDim[1];
            threadsDim[2] = prop.maxThreadsDim[2];

            gridDim[0] = prop.maxGridSize[0];
            gridDim[1] = prop.maxGridSize[1];
            gridDim[2] = prop.maxGridSize[2];

            // Print the saved values
            std::cout << "Device Information:" << std::endl;
            std::cout << "Global memory: " << globalMemory / (1024.0 * 1024.0) << " MB" << std::endl;
            std::cout << "Shared memory per block: " << sharedMemory / 1024.0 << " KB" << std::endl;
            std::cout << "Constant memory: " << constantMemory / 1024.0 << " KB" << std::endl;
            std::cout << "Texture alignment: " << textureAlignment << " bytes" << std::endl;
            std::cout << "Maximum threads per block: " << maxThreadsPerBlock << std::endl;
            std::cout << "Threads per block dimensions: "
                    << threadsDim[0] << " x " << threadsDim[1] << " x " << threadsDim[2] << std::endl;
            std::cout << "Grid dimensions: "
                    << gridDim[0] << " x " << gridDim[1] << " x " << gridDim[2] << std::endl;
        } else {
            std::cerr << "Failed to query device properties: " << cudaGetErrorString(err) << std::endl;
        }

        omp_set_num_threads(num_threadss);

        // Open input file, gzip -df compressed_file1.gz
        if(!strcmp((vcf_filename + strlen(vcf_filename) - 3), ".gz")){
            unzip_gz_file(vcf_filename);
        }
        
        ifstream inFile(filename);
        if(!inFile){
            cout << "ERROR: cannot open file " << filename << endl;
        }
        // Saving filename
        get_filename(filename);
        
        // Getting filesize (number of char in the file)
        auto before = chrono::system_clock::now();
        get_file_size(filename);
        auto after = chrono::system_clock::now();
        auto get_file_size = std::chrono::duration<double>(after - before).count();
        // Getting the header (Saving the header into a string and storing the header size )
        before = chrono::system_clock::now();
        get_and_parse_header(&inFile); //serve per separare l'header dal resto del file
        //vcf.print_header();
        after = chrono::system_clock::now();
        auto get_header = std::chrono::duration<double>(after - before).count();
        inFile.close();
        // Allocating the filestring (the variations as a big char*, the dimension is: filesize - header_size)
        allocate_filestring();
        // Populate filestring and getting the number of lines (num_lines), saving the starting char index of each lines
        before = chrono::system_clock::now();
        find_new_lines_index(filename, num_threadss); 
        after = chrono::system_clock::now();
        auto find_new_lines = std::chrono::duration<double>(after - before).count();

        auto populate_var_struct = std::chrono::duration<double>(after - before).count();
        before = chrono::system_clock::now();
        create_info_vectors(num_threadss);
        reserve_var_columns();
        create_sample_vectors(num_threadss);
        after = chrono::system_clock::now();
        auto reserve_var_columns = std::chrono::duration<double>(after - before).count();
        
         /*
        TODO:
            AAA - TODO 
            0 - capire se posso passare a python una stringa con il padding -> Va riparsato!!
            1 - allocare memoria var_column e samp_columns su GPU -> non più di 80% global, (string a 10 char + heap + atomic o mutex)
            2 - scrivere i kernel dedicati e lanciarli
            3 - scaricare i risultati da GPU sulle strutture dell'host:
                cudaMemcpy(cpuStruct.i_int.data(), gpuStruct.i_int, size * sizeof(int), cudaMemcpyDeviceToHost);
        */

       //Allocate and initialize device memory
        device_allocation();

        before = chrono::system_clock::now();
        populate_var_columns(num_threadss);
        after = chrono::system_clock::now();
        auto populate_var_columns = std::chrono::duration<double>(after - before).count();
        
        cout << "Get file size: " << get_file_size << " s" << endl;
        cout << "get_header: " << get_header << " s" << endl;
        cout << "find_new_lines: " << find_new_lines << " s" << endl;
        cout << "populate_var_struct: " << populate_var_struct << " s" << endl;
        cout << "reserve: " << reserve_var_columns << " s" << endl;
        cout << "populate_var_columns: " << populate_var_columns << " s" << endl;
        free(filestring);
        free(new_lines_index);
    }
    
    void copyMapToConstantMemory(const std::map<std::string, char>& map) {
        char h_keys[NUM_KEYS_GT][MAX_KEY_LENGTH_GT] = {0};
        char h_values[NUM_KEYS_GT] = {0};

        size_t index = 0;
        for (const auto& [key, value] : map) {
            if (index >= NUM_KEYS_GT) break;

            // Copia la chiave, con padding se necessario
            std::strncpy(h_keys[index], key.c_str(), MAX_KEY_LENGTH_GT - 1);

            // Copia il valore corrispondente
            h_values[index] = value;

            ++index;
        }
        cudaMemcpyToSymbol(d_keys_gt, h_keys, sizeof(h_keys));
        cudaMemcpyToSymbol(d_values_gt, h_values, sizeof(h_values));
    }

    void initialize_map1(const std::map<std::string, int> &my_map){
        
        if(my_map.size() > NUM_KEYS_MAP1) {
            std::cerr << "Too many keys." << std::endl;
            return;
        }

        // Host buffers
        char h_keys[NUM_KEYS_MAP1][MAX_KEY_LENGTH_MAP1] = {0};
        int h_values[NUM_KEYS_MAP1] = {0};

        // Copy keys and values into host buffers
        int i = 0;
        for (const auto& [key, value] : my_map) {
            strncpy(h_keys[i], key.c_str(), MAX_KEY_LENGTH_MAP1 - 1);
            h_values[i] = value;
            ++i;
        }

        // Copy to device memory
        cudaMemcpyToSymbol(d_keys_map1, h_keys, sizeof(h_keys));
        cudaMemcpyToSymbol(d_values_map1, h_values, sizeof(h_values));
    }

    void device_allocation(){

       /*
        Per ora:
        unsigned int *d_VC_var_number;
        unsigned int *d_VC_pos;
        __half *d_VC_qual;
        info_float_d *d_VC_in_float;
        info_flag_d *d_VC_in_flag;
        info_int_d *d_VC_in_int;

        TODO - serve la info_map1 su gpu - costant memory
       */

        cudaMalloc(&d_VC_var_number, (num_lines-1)*sizeof(unsigned int));
        cudaMalloc(&d_VC_pos, (num_lines-1)*sizeof(unsigned int));
        cudaMalloc(&d_VC_qual, (num_lines-1)*sizeof(__half));

        int tmp = var_columns.in_float.size();
        cudaMalloc(&d_VC_in_float, tmp*sizeof(info_float_d));
        for(int i=0; i < tmp; i++){
            cudaMalloc(&(d_VC_in_float[i].i_float), (num_lines-1)*sizeof(__half));
            cudaMalloc(&(d_VC_in_float[i].name), var_columns.in_float[i].name.size() + 1);
            cudaMemcpy(d_VC_in_float[i].name, var_columns.in_float[i].name.c_str(), var_columns.in_float[i].name.size() + 1, cudaMemcpyHostToDevice);
        }

        tmp = var_columns.in_flag.size();
        cudaMalloc(&d_VC_in_flag, tmp*sizeof(info_flag_d));
        for(int i=0; i < tmp; i++){
            cudaMalloc(&(d_VC_in_flag[i].i_flag), (num_lines-1)*sizeof(bool));
            cudaMalloc(&(d_VC_in_flag[i].name), var_columns.in_flag[i].name.size() + 1);
            cudaMemcpy(d_VC_in_flag[i].name, var_columns.in_flag[i].name.c_str(), var_columns.in_flag[i].name.size() + 1, cudaMemcpyHostToDevice);
        }

        tmp = var_columns.in_int.size();
        cudaMalloc(&d_VC_in_int, tmp*sizeof(info_int_d));
        for(int i=0; i < tmp; i++){
            cudaMalloc(&(d_VC_in_int[i].i_int), (num_lines-1)*sizeof(int));
            cudaMalloc(&(d_VC_in_int[i].name), var_columns.in_int[i].name.size() + 1);
            cudaMemcpy(d_VC_in_int[i].name, var_columns.in_int[i].name.c_str(), var_columns.in_int[i].name.size() + 1, cudaMemcpyHostToDevice);
        }

        /*
        unsigned int *d_SC_var_id;
        unsigned short *d_SC_samp_id;
        samp_Float *d_SC_samp_float;
        samp_Flag *d_SC_samp_flag;
        samp_Int *d_SC_samp_int;
        samp_GT *d_SC_sample_GT;
        map<string, char> GTMap;
        */

       if(hasDetSamples){ //TODO rialloca dimensione giusta
            copyMapToConstantMemory(samp_columns.GTMap);
            initialize_map1(var_columns.info_map1);
            cudaMalloc(&d_SC_var_id, (num_lines-1)*samp_columns.numSample*sizeof(unsigned int));
            cudaMemset(d_SC_var_id, 0, (num_lines-1)*sizeof(unsigned int));
            cudaMalloc(&d_SC_samp_id, (num_lines-1)*sizeof(unsigned short));
            cudaMemset(d_SC_samp_id, 0, (num_lines-1)*sizeof(unsigned short));

            tmp = samp_columns.samp_float.size();
            cudaMalloc(&d_SC_samp_float, tmp*sizeof(samp_Float_d));
            for(int i=0; i < tmp; i++){
                cudaMalloc(&(d_SC_samp_float[i].i_float), (num_lines-1)*sizeof(__half));
                cudaMalloc(&(d_SC_samp_float[i].name), samp_columns.samp_float[i].name.size() + 1);
                cudaMemcpy(d_SC_samp_float[i].name, samp_columns.samp_float[i].name.c_str(), samp_columns.samp_float[i].name.size() + 1, cudaMemcpyHostToDevice);
                cudaMemcpy(&(d_SC_samp_float[i].numb), &(samp_columns.samp_float[i].numb), sizeof(int), cudaMemcpyHostToDevice);
            }

            tmp = samp_columns.samp_flag.size();
            cudaMalloc(&d_SC_samp_flag, tmp*sizeof(samp_Flag_d));
            for(int i=0; i < tmp; i++){
                cudaMalloc(&(d_SC_samp_flag[i].i_flag), (num_lines-1)*sizeof(bool));
                cudaMalloc(&(d_SC_samp_flag[i].name), samp_columns.samp_flag[i].name.size() + 1);
                cudaMemcpy(d_SC_samp_flag[i].name, samp_columns.samp_flag[i].name.c_str(), samp_columns.samp_flag[i].name.size() + 1, cudaMemcpyHostToDevice);
                cudaMemcpy(&(d_SC_samp_flag[i].numb), &(samp_columns.samp_flag[i].numb), sizeof(int), cudaMemcpyHostToDevice);
            }

            tmp = samp_columns.samp_int.size();
            cudaMalloc(&d_SC_samp_int, tmp*sizeof(samp_Int_d));
            for(int i=0; i < tmp; i++){
                cudaMalloc(&(d_SC_samp_int[i].i_int), (num_lines-1)*sizeof(int));
                cudaMalloc(&(d_SC_samp_int[i].name), samp_columns.samp_int[i].name.size() + 1);
                cudaMemcpy(d_SC_samp_int[i].name, samp_columns.samp_int[i].name.c_str(), samp_columns.samp_int[i].name.size() + 1, cudaMemcpyHostToDevice);
                cudaMemcpy(&(d_SC_samp_int[i].numb), &(samp_columns.samp_int[i].numb), sizeof(int), cudaMemcpyHostToDevice);
            }

            tmp = samp_columns.sample_GT.size();
            cudaMalloc(&d_SC_sample_GT, tmp*sizeof(samp_GT_d));
            for(int i=0; i < tmp; i++){
                cudaMalloc(&(d_SC_sample_GT[i].GT), (num_lines-1)*sizeof(char));
                cudaMemcpy(&(d_SC_sample_GT[i].numb), &(samp_columns.sample_GT[i].numb), sizeof(int), cudaMemcpyHostToDevice);
            }
        }
    }

    void device_free(){
        cudaFree(&d_VC_var_number);
        cudaFree(&d_VC_pos);
        cudaFree(&d_VC_qual);

        int tmp = var_columns.in_float.size();
        
        for(int i=0; i < tmp; i++){
            cudaFree(&d_VC_in_float[i].i_float);
        }
        cudaFree(&d_VC_in_float);

        tmp = var_columns.in_flag.size();
        for(int i=0; i < tmp; i++){
            cudaFree(&d_VC_in_flag[i].i_flag);
        }
        cudaFree(&d_VC_in_flag);

        tmp = var_columns.in_int.size();
        for(int i=0; i < tmp; i++){
            cudaFree(&d_VC_in_int[i].i_int);
        }
        cudaFree(&d_VC_in_int);

        cudaFree(&d_filestring);
        cudaFree(&d_new_lines_index);

       if(hasDetSamples){
            cudaFree(&d_SC_var_id);
            cudaFree(&d_SC_samp_id);

            tmp = samp_columns.samp_float.size();
            for(int i=0; i < tmp; i++){
                cudaFree(&d_SC_samp_float[i].i_float);
            }
            cudaFree(&d_SC_samp_float);

            tmp = samp_columns.samp_flag.size();
            for(int i=0; i < tmp; i++){
                cudaFree(&d_SC_samp_flag[i].i_flag);
            }
            cudaFree(&d_SC_samp_flag);

            tmp = samp_columns.samp_int.size();
            for(int i=0; i < tmp; i++){
                cudaFree(&d_SC_samp_int[i].i_int);
            }
            cudaFree(&d_SC_samp_int);

            tmp = samp_columns.sample_GT.size();
            for(int i=0; i < tmp; i++){
                cudaFree(&d_SC_sample_GT[i].GT);
            }
            cudaFree(&d_SC_sample_GT);
        }
    }

    void find_new_lines_index(string w_filename, int num_threads){
        // Allocate memory for the `new_lines_index` array. The size is exaggerated (assuming every character is a new line).
        // The first element is set to 0, indicating the start of the first line.
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
        num_lines = tmp_num_lines[0];
        for(int i=1; i<num_threads; i++){
            num_lines= num_lines + tmp_num_lines[i];
        }
        new_lines_index = (unsigned int*)malloc(sizeof(unsigned int)*(num_lines+1));
        new_lines_index[0] = 0;

        cudaMalloc(&d_filestring, variants_size*sizeof(char));
        cudaMalloc(&d_new_lines_index, sizeof(unsigned int)*(num_lines+1));
        cudaMemcpy(d_filestring, filestring, sizeof(char)*variants_size, cudaMemcpyHostToDevice);
        unsigned int* d_count;
        cudaMalloc(&d_count, sizeof(unsigned int));
        cudaMemset(d_count, 0, sizeof(long int));
        dim3 threads = 1024;
        dim3 blocks = ceil(variants_size/1024);
        cu_find_new_lines_index<<<threads, blocks>>>(
            d_filestring,
            variants_size,
            d_new_lines_index,
            num_lines+1,
            d_count
        );

        //ordering with Thrust library
        thrust::device_ptr<unsigned int> d_new_lines_index_ptr = thrust::device_pointer_cast(d_new_lines_index);
        thrust::sort(d_new_lines_index_ptr, d_new_lines_index_ptr + (num_lines+1));
        cudaMemcpy(new_lines_index, d_new_lines_index, sizeof(unsigned int)*(num_lines+1), cudaMemcpyDeviceToHost);

        //da capire cosa lasciare su GPU e cosa liberare
        cudaFree(&d_count);
        cudaFree(&d_new_lines_index_ptr);
        cudaFree(&d_filestring);
        
        for(int i = 0; i <= num_lines; i++){
            cout << new_lines_index[i] << " | ";
        }
        cout << endl;
    }

    void get_filename(string path_filename){
        vector<string> line_el;
        path_to_filename = path_filename;
        boost::split(line_el, path_filename, boost::is_any_of("/"));
        filename = line_el[line_el.size()-1];
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
                            if(line_el2[1] == "GT"){
                                FORMAT.hasGT = true;
                                boost::split(line_el2, line_el1[1], boost::is_any_of("="));
                                FORMAT.numGT = line_el2[1][0];
                                i+=3;
                            }else{
                                FORMAT.ID.push_back(line_el2[1]);
                            }
                        } 
                        if(i==1) FORMAT.Number.push_back(line_el2[1]);
                        if(i==1 && line_el2[1] == "A"){
                            FORMAT.alt_values++;
                        }else{
                            hasDetSamples = true;
                        }
                        if(i==2) FORMAT.Type.push_back(line_el2[1]);
                    }
                }
            }
        }

        vector<string> tmp_split;
        boost::split(tmp_split, line, boost::is_any_of("\t "));
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
        
        cout << "Num Samples = " << samp_columns.numSample << endl;

        INFO.total_values = INFO.ID.size();
        INFO.no_alt_values = INFO.total_values - INFO.alt_values;

        header_size += line.length() + 1;

        variants_size = filesize - header_size; // New size without the header
    }   
    
    void allocate_filestring(){
        filestring = (char*)malloc(variants_size);
    }

    void create_sample_vectors(int num_threads){
        samp_Flag samp_flag_tmp;
        samp_Float samp_float_tmp;
        samp_Int samp_int_tmp;
        samp_String samp_string_tmp;

        samp_Float samp_alt_float_tmp;
        samp_Int samp_alt_int_tmp;
        samp_String samp_alt_string_tmp;

        int numIter = FORMAT.ID.size();

        if(FORMAT.hasGT && FORMAT.numGT == 'A'){
            alt_sample.initMapGT();
            alt_sample.sample_GT.numb = -1;
        }else if(FORMAT.hasGT){
            int iter = FORMAT.numGT - '0';//std::stoi(FORMAT.numGT); //number of vector needed TODO - da sistemare, se no non possiamo avere numb=10...
            samp_columns.initMapGT();
            for(int i=0; i<iter; i++){ //create a vector of vectors
                samp_GT tmp;
                tmp.numb = iter;
                tmp.GT.resize((num_lines-1)*samp_columns.numSample, (char)0);
                samp_columns.sample_GT.push_back(tmp);
            }
            samp_columns.sample_GT.resize(FORMAT.numGT);   
        }

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
                    //alt_sample.samp_string.back().i_string.resize(batch_size*samp_columns.numSample*2, "\0");
                    alt_sample.samp_string.back().numb = -1;
                    info_map[FORMAT.ID[i]] = 11;
                    var_columns.info_map1[FORMAT.ID[i]] = 11;
                    FORMAT.strings_alt++;
                }else if(!strcmp(&FORMAT.Type[i][0], "Integer")){
                    samp_alt_int_tmp.name = FORMAT.ID[i];
                    alt_sample.samp_int.push_back(samp_alt_int_tmp);
                    //alt_sample.samp_int.back().i_int.resize(batch_size*samp_columns.numSample*2, 0);
                    alt_sample.samp_int.back().numb = -1;
                    info_map[FORMAT.ID[i]] = 12;
                    var_columns.info_map1[FORMAT.ID[i]] = 12;
                    FORMAT.ints_alt++;
                }else if(!strcmp(&FORMAT.Type[i][0], "Float")){
                    samp_alt_float_tmp.name = FORMAT.ID[i];
                    alt_sample.samp_float.push_back(samp_alt_float_tmp);
                    //alt_sample.samp_float.back().i_float.resize(batch_size*samp_columns.numSample*2, 0);
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
        if(hasDetSamples){
            samp_columns.var_id.resize((num_lines-1)*samp_columns.numSample, 0);
            samp_columns.samp_id.resize((num_lines-1)*samp_columns.numSample, static_cast<unsigned short>(0));
        }    
        alt_sample.samp_flag.resize(FORMAT.flags_alt);
        alt_sample.samp_int.resize(FORMAT.ints_alt);
        alt_sample.samp_float.resize(FORMAT.floats_alt);
        alt_sample.samp_string.resize(FORMAT.strings_alt);
        alt_sample.var_id.resize((num_lines-1)* alt_sample.numSample, 0);
        alt_sample.samp_id.resize((num_lines-1)*alt_sample.numSample, static_cast<unsigned short>(0));

    }
    
    void create_info_vectors(int num_threads){
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
                    //alt_int_tmp.i_int.resize(2*batch_size, 0);
                    alt_columns.alt_int.push_back(alt_int_tmp);
                    info_map[INFO.ID[i]] = 4;
                    var_columns.info_map1[INFO.ID[i]] = 4;
                }
                if(strcmp(&INFO.Type[i][0], "Float")==0){
                    INFO.floats_alt++;
                    alt_float_tmp.name = INFO.ID[i];
                    //alt_float_tmp.i_float.resize(2*batch_size, 0);
                    alt_columns.alt_float.push_back(alt_float_tmp);
                    info_map[INFO.ID[i]] = 5;
                    var_columns.info_map1[INFO.ID[i]] = 5;
                }
                if(strcmp(&INFO.Type[i][0], "String")==0){
                    INFO.strings_alt++;
                    alt_string_tmp.name = INFO.ID[i];
                    //alt_string_tmp.i_string.resize(2*batch_size, "\0");
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
                } else if(strcmp(&INFO.Type[i][0], "Float")==0){
                    INFO.floats++;
                    info_float_tmp.name = INFO.ID[i];
                    info_float_tmp.i_float.resize(num_lines-1, 0);
                    var_columns.in_float.push_back(info_float_tmp);
                    info_map[INFO.ID[i]] = 2;
                    var_columns.info_map1[INFO.ID[i]] = 2;
                } else if(strcmp(&INFO.Type[i][0], "String")==0){
                    INFO.strings++;
                    info_string_tmp.name = INFO.ID[i];
                    info_string_tmp.i_string.resize(num_lines-1, "\0");
                    var_columns.in_string.push_back(info_string_tmp);
                    info_map[INFO.ID[i]] = 3;
                    var_columns.info_map1[INFO.ID[i]] = 3;
                } else if(strcmp(&INFO.Type[i][0], "Flag")==0){
                    INFO.flags++;
                    info_flag_tmp.name = INFO.ID[i];
                    info_flag_tmp.i_flag.resize(num_lines-1, 0);
                    var_columns.in_flag.push_back(info_flag_tmp);
                    info_map[INFO.ID[i]] = 0;
                    var_columns.info_map1[INFO.ID[i]] = 0;
                }
            }else if(/*fai il punto*/ false){
                
            }else{
                //in progress, se num > 1 TODO
                // Può avere solo come valori: 0, 1, R, A, G, .
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
                cout << static_cast<float>(var_columns.in_float[i].i_float[j]) << " ";
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
        /*cout<<endl;
        cout<<"Flags: "<<endl;
        for(int i=0; i<var_columns.in_flag.size(); i++){
            cout<<var_columns.in_flag[i].name<<": ";
            for(int j=0; j<var_columns.in_flag[i].i_flag.size(); j++){
                cout<<var_columns.in_flag[i].i_flag[j]<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
        cout<<"Floats: "<<endl;
        for(int i=0; i<var_columns.in_float.size(); i++){
            cout<<var_columns.in_float[i].name<<": ";
            for(int j=0; j<var_columns.in_float[i].i_float.size(); j++){
                cout<<var_columns.in_float[i].i_float[j]<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
        cout<<"Strings: "<<endl;
        for(int i=0; i<var_columns.in_string.size(); i++){
            cout<<var_columns.in_string[i].name<<": ";
            for(int j=0; j<var_columns.in_string[i].i_string.size(); j++){
                cout<<var_columns.in_string[i].i_string[j]<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
        cout<<"Ints: "<<endl;
        for(int i=0; i<var_columns.in_int.size(); i++){
            cout<<var_columns.in_int[i].name<<": ";
            for(int j=0; j<var_columns.in_int[i].i_int.size(); j++){
                cout<<var_columns.in_int[i].i_int[j]<<" ";
            }
            cout<<endl;
        }
        cout<<endl;*/
    }
    
    void reserve_var_columns(){
        var_columns.var_number.resize(num_lines-1);
        var_columns.chrom.resize(num_lines-1);
        var_columns.id.resize(num_lines-1);
        var_columns.pos.resize(num_lines-1);
        var_columns.ref.resize(num_lines-1); 
        var_columns.qual.resize(num_lines-1);
        var_columns.filter.resize(num_lines-1);
    }
    
    void populate_var_columns(int num_threads){
        int threadsPerBlock = 1024;
        int blocksPerGrid = std::ceil(num_lines/threadsPerBlock);

       if(samplesON){
            get_vcf_line_format_kernel<<<blocksPerGrid, threadsPerBlock>>>(
                d_filestring,
                d_VC_var_number,
                d_VC_pos,
                d_VC_qual,
                d_VC_in_float,
                d_VC_in_flag,
                d_VC_in_int,
                d_new_lines_index,
                d_SC_var_id,
                d_SC_samp_id,
                d_SC_samp_float,
                d_SC_samp_flag,
                d_SC_samp_int,
                d_SC_sample_GT,
                samp_columns.numSample,
                num_lines,
                FORMAT.numGT
            );
       }else{
            get_vcf_line_kernel<<<blocksPerGrid, threadsPerBlock>>>(
                d_filestring,
                d_VC_var_number,
                d_VC_pos,
                d_VC_qual,
                d_VC_in_float,
                d_VC_in_flag,
                d_VC_in_int,
                d_new_lines_index,
                num_lines
            );
       }

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
            tmp_alt[th_ID].init(alt_columns, INFO, batch_size);
            
            tmp_num_alt[th_ID] = 0;
            
            tmp_num_alt_format[th_ID] = 0;
            tmp_alt_format[th_ID].var_id.resize(batch_size*2*samp_columns.numSample, 0);
            tmp_alt_format[th_ID].alt_id.resize(batch_size*2*samp_columns.numSample, 0);
            tmp_alt_format[th_ID].samp_id.resize(batch_size*2*samp_columns.numSample, static_cast<unsigned short>(0));
            
            start = th_ID * batch_size; // Starting point of the thread's batch
            end = start + batch_size; // Ending point of the thread's batch
            
            if(samplesON){
                // There are samples in the dataset
                tmp_alt_format[th_ID].init(alt_sample, FORMAT, batch_size);
                tmp_num_alt_format[th_ID] = 0;
                tmp_alt_format[th_ID].var_id.resize(batch_size*2*samp_columns.numSample, 0);
                tmp_alt_format[th_ID].alt_id.resize(batch_size*2*samp_columns.numSample, 0);
                tmp_alt_format[th_ID].samp_id.resize(batch_size*2*samp_columns.numSample, static_cast<unsigned short>(0));
                if(FORMAT.hasGT && FORMAT.numGT == 'A'){
                    tmp_alt_format[th_ID].sample_GT.GT.resize(batch_size*2*samp_columns.numSample, (char)0),
                    tmp_alt_format[th_ID].initMapGT();
                }
                
                // For each line in the batch
                for(long i=start; i<end && i<num_lines-1; i++){ 
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

        cudaDeviceSynchronize();
        cudaMemcpy(var_columns.var_number.data(), d_VC_var_number, (num_lines-1)*sizeof(unsigned int), cudaMemcpyDeviceToHost);
        cudaMemcpy(var_columns.pos.data(), d_VC_pos, (num_lines-1)*sizeof(unsigned int), cudaMemcpyDeviceToHost);
        cudaMemcpy(var_columns.qual.data(), d_VC_qual, (num_lines-1)*sizeof(__half), cudaMemcpyDeviceToHost);

        for(int i=0; i<var_columns.in_float.size(); i++){
            cudaMemcpy(var_columns.in_float[i].i_float.data(), &(d_VC_in_float[i].i_float[0]), (num_lines-1)*sizeof(__half), cudaMemcpyDeviceToHost);
        }

        for(int i=0; i<var_columns.in_flag.size(); i++){
            cudaMemcpy(var_columns.in_flag[i].i_flag.data(), &(d_VC_in_flag[i].i_flag[0]), (num_lines-1)*sizeof(bool), cudaMemcpyDeviceToHost);
        }

        for(int i=0; i<var_columns.in_int.size(); i++){
            cudaMemcpy(var_columns.in_int[i].i_int.data(), &(d_VC_in_int[i].i_int[0]), (num_lines-1)*sizeof(int), cudaMemcpyDeviceToHost);
        }

        if(samplesON){ //TODO samp det sono numlines*numsamp e non solo num lines controlla allocazione
            cudaMemcpy(samp_columns.var_id.data(), d_SC_var_id, (num_lines-1)*samp_columns.numSample*sizeof(unsigned int), cudaMemcpyDeviceToHost);
            cudaMemcpy(samp_columns.samp_id.data(), d_SC_samp_id, (num_lines-1)*samp_columns.numSample*sizeof(unsigned short), cudaMemcpyDeviceToHost);

            for(int i=0; i<samp_columns.samp_float.size(); i++){
                cudaMemcpy(samp_columns.samp_float[i].i_float.data(), &(d_SC_samp_float[i].i_float[0]), (num_lines-1)*samp_columns.numSample*sizeof(__half), cudaMemcpyDeviceToHost);
            }
            
            for(int i=0; i<samp_columns.samp_flag.size(); i++){
                cudaMemcpy(samp_columns.samp_flag[i].i_flag.data(), d_VC_in_flag[i].i_flag, (num_lines - 1) * samp_columns.numSample * sizeof(bool), cudaMemcpyDeviceToHost);
            }

            for(int i=0; i<samp_columns.samp_int.size(); i++){
                cudaMemcpy(samp_columns.samp_int[i].i_int.data(), &(d_SC_samp_int[i].i_int[0]), (num_lines-1)*samp_columns.numSample*sizeof(int), cudaMemcpyDeviceToHost);
            }

            for(int i=0; i<samp_columns.sample_GT.size(); i++){
                cudaMemcpy(samp_columns.sample_GT[i].GT.data(), &(d_SC_sample_GT[i].GT[0]), (num_lines-1)*samp_columns.numSample*sizeof(char), cudaMemcpyDeviceToHost);
            }
            
        }

        device_free();

        int totAlt = 0;
        int totSampAlt = 0;
        for(int i=0; i<num_threads; i++){
            alt_columns.var_id.insert(
                alt_columns.var_id.end(),
                std::make_move_iterator(tmp_alt[i].var_id.begin()),
                std::make_move_iterator(tmp_alt[i].var_id.end())
            );
            alt_columns.alt_id.insert(
                alt_columns.alt_id.end(),
                std::make_move_iterator(tmp_alt[i].alt_id.begin()),
                std::make_move_iterator(tmp_alt[i].alt_id.end())
            );
            alt_columns.alt.insert(
                alt_columns.alt.end(),
                std::make_move_iterator(tmp_alt[i].alt.begin()),
                std::make_move_iterator(tmp_alt[i].alt.end())
            );
            for(int j=0; j<INFO.ints_alt; j++){
                alt_columns.alt_int[j].i_int.insert(
                    alt_columns.alt_int[j].i_int.end(), 
                    std::make_move_iterator(tmp_alt[i].alt_int[j].i_int.begin()), 
                    std::make_move_iterator(tmp_alt[i].alt_int[j].i_int.end()));
            }
            for(int j=0; j<INFO.floats_alt; j++){
                alt_columns.alt_float[j].i_float.insert(
                    alt_columns.alt_float[j].i_float.end(), 
                    std::make_move_iterator(tmp_alt[i].alt_float[j].i_float.begin()), 
                    std::make_move_iterator(tmp_alt[i].alt_float[j].i_float.end()));
            }
            for(int j=0; j<INFO.strings_alt; j++){
                alt_columns.alt_string[j].i_string.insert(
                    alt_columns.alt_string[j].i_string.end(), 
                    std::make_move_iterator(tmp_alt[i].alt_string[j].i_string.begin()), 
                    std::make_move_iterator(tmp_alt[i].alt_string[j].i_string.end()));
            }
            totAlt+=tmp_num_alt[i];
        }
        
        if(samplesON){
            for(int i=0; i<num_threads; i++){
                alt_sample.var_id.insert(
                    alt_sample.var_id.end(),
                    std::make_move_iterator(tmp_alt_format[i].var_id.begin()),
                    std::make_move_iterator(tmp_alt_format[i].var_id.end())
                );
                alt_sample.alt_id.insert(
                    alt_sample.alt_id.end(),
                    std::make_move_iterator(tmp_alt_format[i].alt_id.begin()),
                    std::make_move_iterator(tmp_alt_format[i].alt_id.end())
                );
                alt_sample.samp_id.insert(
                    alt_sample.samp_id.end(),
                    std::make_move_iterator(tmp_alt_format[i].samp_id.begin()),
                    std::make_move_iterator(tmp_alt_format[i].samp_id.end())
                );
                for(int j=0; j<FORMAT.ints_alt; j++){
                    alt_sample.samp_int[j].i_int.insert(
                        alt_sample.samp_int[j].i_int.end(), 
                        std::make_move_iterator(tmp_alt_format[i].samp_int[j].i_int.begin()), 
                        std::make_move_iterator(tmp_alt_format[i].samp_int[j].i_int.end()));
                }
                for(int j=0; j<FORMAT.floats_alt; j++){
                    alt_sample.samp_float[j].i_float.insert(
                        alt_sample.samp_float[j].i_float.end(), 
                        std::make_move_iterator(tmp_alt_format[i].samp_float[j].i_float.begin()), 
                        std::make_move_iterator(tmp_alt_format[i].samp_float[j].i_float.end()));
                }
                for(int j=0; j<FORMAT.strings_alt; j++){
                    alt_sample.samp_string[j].i_string.insert(
                        alt_sample.samp_string[j].i_string.end(), 
                        std::make_move_iterator(tmp_alt_format[i].samp_string[j].i_string.begin()), 
                        std::make_move_iterator(tmp_alt_format[i].samp_string[j].i_string.end()));
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

};

#endif