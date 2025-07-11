/**
 * @file Parser.cu
 * @brief Implements the VCF file parsing workflow.
 *
 * This file defines the vcf_parsed class, which encapsulates the entire VCF parsing process:
 * - Reading and parsing the VCF header.
 * - Allocating host and device memory.
 * - Launching CUDA kernels for parsing VCF lines.
 * - Merging GPU and CPU results into structured data.
 */

#ifndef PARSER_CU
#define PARSER_CU

#include "DataStructures.h"
#include "Kernels.cu"
#include "Utils.h"
#include "DataFrames.h"
#include "CUDAUtils.cuh"
#include "Parser.h"

#include <cuda_runtime.h>
#include <cuda_fp16.h>  
#include <thrust/device_ptr.h> 
#include <thrust/sort.h>

#include <boost/algorithm/string.hpp>
#include <chrono>
#include <fstream>
#include <iostream>  
#include <vector>    
#include <cstring>   
#include <cstdlib>   
#include <map>
#include <omp.h> 
#include <thread>
#include <functional>
#include <future>


using namespace std;

#define CUDA_CHECK_ERROR(call)                             \
    do {                                                   \
        cudaError_t err = call;                            \
        if (err != cudaSuccess) {                          \
            std::cerr << "CUDA error in " << #call         \
                      << " at " << __FILE__ << ":" << __LINE__ \
                      << " - " << cudaGetErrorString(err)    \
                      << std::endl;                          \
            exit(EXIT_FAILURE);                            \
        }                                                  \
    } while (0)



/**
    * @brief Runs the VCF parsing process.
    *
    * Performs the following steps:
    *  - Initializes CUDA device and queries device properties.
    *  - Opens the VCF file (uncompressing if needed) and extracts the header.
    *  - Allocates memory for the file content and identifies the start of each variant line.
    *  - Creates and reserves vectors for variant and sample data.
    *  - Allocates and initializes device memory and lookup maps.
    *  - Launches CUDA kernels to parse the VCF lines.
    *  - Copies the parsed results back to host memory and frees device memory.
    *
    * @param vcf_filename Path to the VCF file.
    * @param num_threadss Number of threads to use for OpenMP parallel processing.
    */
void vcf_parsed::run(char* vcf_filename, int num_threadss){
    string filename = vcf_filename; 
    string line;
    vcf_parsed vcf;

    // Variables to hold device information
    size_t globalMemory = 0;       // Total global memory
    size_t sharedMemory = 0;       // Shared memory per block
    size_t constantMemory = 0;     // Constant memory
    size_t textureAlignment = 0;   // Texture alignment
    int maxThreadsPerBlock = 0;    // Maximum threads per block
    int cudaCores;                 // Total number of CUDA cores
    int threadsDim[3] = {0};       // Maximum threads per block dimensions
    int gridDim[3] = {0};          // Maximum grid dimensions

    // Query device properties
    int deviceCount = 0;
    cudaError_t error = cudaGetDeviceCount(&deviceCount);
    if (error != cudaSuccess || deviceCount == 0) {
        std::cerr << "No CUDA-capable devices found. Exiting..." << std::endl;
        exit(1);
    }

    int deviceID = 0; // Usa il dispositivo 0 (il primo)
    cudaSetDevice(deviceID);

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

        // Determine number of CUDA cores per SM based on compute capability
        int coresPerSM = 0;
        if (prop.major == 1) {
            // Tesla architecture
            coresPerSM = 8;
        } else if (prop.major == 2) {
            // Fermi architecture
            coresPerSM = (prop.minor == 0 || prop.minor == 1) ? 32 : 48;
        } else if (prop.major == 3) {
            // Kepler architecture
            coresPerSM = 192;
        } else if (prop.major == 5) {
            // Maxwell architecture
            coresPerSM = 128;
        } else if (prop.major == 6) {
            // Pascal architecture
            coresPerSM = (prop.minor == 1 || prop.minor == 2) ? 128 : 64;
        } else if (prop.major == 7) {
            // Volta or Turing architecture
            coresPerSM = (prop.minor == 0) ? 64 : 64;  // Adjust if needed
        } else if (prop.major == 8) {
            // Ampere architecture
            coresPerSM = (prop.minor == 0) ? 64 : (prop.minor == 6 ? 128 : 64);
        } else {
            // Fallback assumption
            coresPerSM = 128;
        }

        cudaCores = coresPerSM * prop.multiProcessorCount;

        //Print the saved values
        //std::cout << "Device Information:" << std::endl;
        //std::cout << "Global memory: " << globalMemory / (1024.0 * 1024.0) << " MB" << std::endl;
        //std::cout << "Shared memory per block: " << sharedMemory / 1024.0 << " KB" << std::endl;
        //std::cout << "Constant memory: " << constantMemory / 1024.0 << " KB" << std::endl;
        //std::cout << "Texture alignment: " << textureAlignment << " bytes" << std::endl;
        //std::cout << "Maximum threads per block: " << maxThreadsPerBlock << std::endl;
        //std::cout << "Threads per block dimensions: " << threadsDim[0] << " x " << threadsDim[1] << " x " << threadsDim[2] << std::endl;
        //std::cout << "Grid dimensions: " << gridDim[0] << " x " << gridDim[1] << " x " << gridDim[2] << std::endl;
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
    filename = get_filename(filename, path_to_filename);
    
    // Getting filesize (number of char in the file)
    filesize = get_file_size(path_to_filename);
    // Getting the header (Saving the header into a string and storing the header size )
    get_and_parse_header(&inFile); //serve per separare l'header dal resto del file
    inFile.close();
    // Allocating the filestring (the variations as a big char*, the dimension is: filesize - header_size)
    allocate_filestring();
    // Populate filestring and getting the number of lines (num_lines), saving the starting char index of each lines
    find_new_lines_index(path_to_filename, num_threadss);
    create_info_vectors(num_threadss);
    reserve_var_columns();
    create_sample_vectors(num_threadss);
    //Allocate and initialize device memory
    device_allocation();
    populate_var_columns(num_threadss, cudaCores);
    device_free();

}
    
/**
    * @brief Copies a genotype map to device constant memory.
    *
    * Copies the key-value pairs from the provided host map into the device constant
    * memory arrays (d_keys_gt and d_values_gt).
    *
    * @param map Host map with genotype keys and corresponding char values.
    */
void vcf_parsed::copyMapToConstantMemory(const std::map<std::string, char>& map) {
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

/**
    * @brief Initializes the INFO field lookup map (Map1) in device constant memory.
    *
    * Copies keys and integer values from the host map to device constant memory.
    * Prints an error if the number of keys exceeds the allowed limit.
    *
    * @param my_map Host map containing keys and their corresponding integer values.
    */
void vcf_parsed::initialize_map1(const std::map<std::string, int> &my_map){
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
    CUDA_CHECK_ERROR(cudaMemcpyToSymbol(d_keys_map1, h_keys, sizeof(h_keys)));
    CUDA_CHECK_ERROR(cudaMemcpyToSymbol(d_values_map1, h_values, sizeof(h_values)));
}

/**
    * @brief Allocates device memory for VCF parsing data.
    *
    * Allocates memory on the GPU for variant numbers, positions, quality scores,
    * and INFO fields. If sample data is present, it also allocates memory for sample fields.
    */
void vcf_parsed::device_allocation(){
    cudaMalloc(&d_VC_var_number, (num_lines) * sizeof(unsigned int));
    cudaMalloc(&d_VC_pos, (num_lines) * sizeof(unsigned int));
    cudaMalloc(&d_VC_qual, (num_lines) * sizeof(__half));

    int tmp = var_columns.in_float.size();

    const int max_name_size = 16;
    // allocazione di tutti i campi float in successione
    cudaMalloc(&(d_VC_in_float->i_float), tmp * (num_lines) * sizeof(__half)); //Va in segfault qui
    cudaMalloc(&(d_VC_in_float->name), tmp * sizeof(char) * max_name_size); //allocazione in successione di tutti i nomi (assumo massimo 16 caratteri)

    for (int i = 0; i < tmp; i++) {
        cudaMemcpy(d_VC_in_float->name+i*max_name_size, var_columns.in_float[i].name.c_str(), var_columns.in_float[i].name.size() + 1, cudaMemcpyHostToDevice);
    }

    tmp = var_columns.in_flag.size();
    cudaMalloc(&(d_VC_in_flag->i_flag), tmp * (num_lines) * sizeof(bool));
    cudaMalloc(&(d_VC_in_flag->name), tmp * sizeof(char) * max_name_size);

    for (int i = 0; i < tmp; i++) {
        cudaMemcpy(d_VC_in_flag->name+i*max_name_size, var_columns.in_flag[i].name.c_str(), var_columns.in_flag[i].name.size() + 1, cudaMemcpyHostToDevice);
    }

    // Ripeti per in_int
    tmp = var_columns.in_int.size();
    cudaMalloc(&(d_VC_in_int->i_int), tmp * (num_lines) * sizeof(int));
    cudaMalloc(&(d_VC_in_int->name), tmp * sizeof(char) * max_name_size);
    for (int i = 0; i < tmp; i++) {
        cudaMemcpy(d_VC_in_int->name+i*max_name_size, var_columns.in_int[i].name.c_str(), var_columns.in_int[i].name.size() + 1, cudaMemcpyHostToDevice);
    }

    initialize_map1(var_columns.info_map1);
    
    if (hasDetSamples) {
        // Copy map to constant memory
        copyMapToConstantMemory(samp_columns.GTMap);
        // Allocate and initialize d_SC_var_id
        cudaMalloc(&d_SC_var_id, (num_lines) * samp_columns.numSample * sizeof(unsigned int));
        cudaMemset(d_SC_var_id, 0, (num_lines) * samp_columns.numSample * sizeof(unsigned int));
        
        // Allocate and initialize d_SC_samp_id
        cudaMalloc(&d_SC_samp_id, (num_lines) * samp_columns.numSample * sizeof(unsigned short));
        cudaMemset(d_SC_samp_id, 0, (num_lines) * sizeof(unsigned short));

        // Allocate and initialize samp_float
        tmp = samp_columns.samp_float.size();
        cudaMalloc(&(d_SC_samp_float->i_float), tmp * (num_lines * samp_columns.numSample) * sizeof(__half));
        cudaMalloc(&(d_SC_samp_float->name), tmp * sizeof(char) * max_name_size);
        cudaMalloc(&(d_SC_samp_float->numb), tmp * sizeof(int));

        for (int i = 0; i < tmp; i++) {
            cudaMemcpy(d_SC_samp_float->name+i*max_name_size, samp_columns.samp_float[i].name.c_str(), samp_columns.samp_float[i].name.size() + 1, cudaMemcpyHostToDevice);
            cudaMemcpy(d_SC_samp_float->numb+i, &(samp_columns.samp_float[i].numb), sizeof(int), cudaMemcpyHostToDevice);
        }

        // Allocate and initialize samp_flag
        tmp = samp_columns.samp_flag.size();
        cudaMalloc(&(d_SC_samp_flag->i_flag), tmp * (num_lines * samp_columns.numSample) * sizeof(bool));
        cudaMalloc(&(d_SC_samp_flag->name), tmp * sizeof(char) * max_name_size);
        cudaMalloc(&(d_SC_samp_flag->numb), tmp * sizeof(int));

        for (int i = 0; i < tmp; i++) {
            cudaMemcpy(d_SC_samp_flag->name+i*max_name_size, samp_columns.samp_flag[i].name.c_str(), samp_columns.samp_flag[i].name.size() + 1, cudaMemcpyHostToDevice);
            cudaMemcpy(d_SC_samp_flag->numb+i, &(samp_columns.samp_flag[i].numb), sizeof(int), cudaMemcpyHostToDevice);
        }

        // Allocate and initialize samp_int
        tmp = samp_columns.samp_int.size();
        cudaMalloc(&(d_SC_samp_int->i_int), tmp * (num_lines * samp_columns.numSample) * sizeof(int));
        cudaMalloc(&(d_SC_samp_int->name), tmp * sizeof(char) * max_name_size);
        cudaMalloc(&(d_SC_samp_int->numb), tmp * sizeof(int));

        for (int i = 0; i < tmp; i++) {
            cudaMemcpy(d_SC_samp_int->name+i*max_name_size, samp_columns.samp_int[i].name.c_str(), samp_columns.samp_int[i].name.size() + 1, cudaMemcpyHostToDevice);
            cudaMemcpy(d_SC_samp_int->numb+i, &(samp_columns.samp_int[i].numb), sizeof(int), cudaMemcpyHostToDevice);
        }

        // Allocate and initialize samp_GT
        tmp = samp_columns.sample_GT.size();
        cudaMalloc(&(d_SC_sample_GT->GT), tmp * (num_lines * samp_columns.numSample) * sizeof(char));
        cudaMalloc(&(d_SC_sample_GT->numb), sizeof(int));
        cudaMemcpy(d_SC_sample_GT->numb, &(samp_columns.sample_GT[0].numb), sizeof(int), cudaMemcpyHostToDevice);

        //TODO - Da rimuovere questa sinchronize utile per debugging:
        //cudaError_t err = cudaDeviceSynchronize();
        //if (err != cudaSuccess) {
        //    fprintf(stderr, "CUDA error in device allocation: %s\n", cudaGetErrorString(err));
        //    exit(EXIT_FAILURE);
        //}
    }

}

/**
    * @brief Frees all allocated device memory.
    *
    * Releases device memory allocated during the parsing process.
    */
void vcf_parsed::device_free() {
    cudaFree(d_VC_var_number);
    cudaFree(d_VC_pos);
    cudaFree(d_VC_qual);
    cudaFree(d_VC_in_float->i_float);
    cudaFree(d_VC_in_float->name);
    cudaFree(d_VC_in_flag->i_flag);
    cudaFree(d_VC_in_flag->name);
    cudaFree(d_VC_in_int->i_int);
    cudaFree(d_VC_in_int->name);
    cudaFree(d_filestring);
    cudaFree(d_new_lines_index);
    cudaFree(d_count);

    if (hasDetSamples) {
        cudaFree(d_SC_var_id);
        cudaFree(d_SC_samp_id);
        cudaFree(d_SC_samp_float->i_float);
        cudaFree(d_SC_samp_float->name);
        cudaFree(d_SC_samp_float->numb);
        cudaFree(d_SC_samp_flag->i_flag);
        cudaFree(d_SC_samp_flag->name);
        cudaFree(d_SC_samp_flag->numb);
        cudaFree(d_SC_samp_int->i_int);
        cudaFree(d_SC_samp_int->numb);
        cudaFree(d_SC_sample_GT->GT);
        cudaFree(d_SC_sample_GT->numb);
        cudaFree(d_params);
    }
}

/**
    * @brief Finds newline indices in the VCF file.
    *
    * Reads the VCF file in parallel using OpenMP to determine the starting index of each line.
    * The indices are stored in an array and copied to device memory for use by CUDA kernels.
    *
    * @param w_filename The path to the VCF file.
    * @param num_threads Number of threads to use for parallel processing.
    */
void vcf_parsed::find_new_lines_index(string w_filename, int num_threads){
    // Allocate memory for the `new_lines_index` array. The size is exaggerated (assuming every character is a new line).
    // The first element is set to 0, indicating the start of the first line.
    num_lines++; // Increment the line counter to account for the first line.
    long tmp_num_lines[num_threads]; // Temporary array to store the number of lines found by each thread.
    variants_size--;
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
    before = chrono::system_clock::now();
    num_lines = tmp_num_lines[0];
    for(int i=1; i<num_threads; i++){
        num_lines= num_lines + tmp_num_lines[i];
    }

    new_lines_index = (unsigned int*)malloc(sizeof(unsigned int)*(num_lines+1));
    new_lines_index[0] = 0;
    cudaMalloc(&d_filestring, (variants_size + 8)* sizeof(char));
    cudaMalloc(&d_new_lines_index, (num_lines + 1) * sizeof(unsigned int));
    cudaMemcpy(d_filestring, filestring, sizeof(char)*variants_size, cudaMemcpyHostToDevice);
    cudaMalloc(&d_count, sizeof(unsigned int));
    cudaMemset(d_count, 0, sizeof(unsigned int));
    
    dim3 threads = 1024;
    dim3 blocks((variants_size + threads.x - 1) / threads.x);
    cu_find_new_lines_index<<<blocks, threads>>>(
        d_filestring,
        variants_size,
        d_new_lines_index,
        num_lines+1,
        d_count
    );

    cudaError_t kernelErr = cudaGetLastError();
    if (kernelErr != cudaSuccess) {
        std::cerr << "Kernel launch error: " << cudaGetErrorString(kernelErr) << std::endl;
        return;
    }

    cudaDeviceSynchronize();

    //ordering with Thrust library
    thrust::device_ptr<unsigned int> d_new_lines_index_ptr(d_new_lines_index);
    if (!d_new_lines_index_ptr) {
        std::cerr << "Invalid device pointer for d_new_lines_index." << std::endl;
        return;
    }
    try {
        thrust::sort(d_new_lines_index_ptr, d_new_lines_index_ptr + (num_lines + 1));
    } catch (thrust::system_error &e) {
        std::cerr << "Thrust error: " << e.what() << std::endl;
    }


    cudaDeviceSynchronize();
    
    cudaMemcpy(new_lines_index, d_new_lines_index, sizeof(unsigned int)*(num_lines+1), cudaMemcpyDeviceToHost);
}
    
/**
    * @brief Reads the VCF header from the input file.
    *
    * Extracts header lines (starting with "##") from the VCF file,
    * storing them in the header string and updating the header size.
    *
    * @param file Pointer to the input file stream.
    */
void vcf_parsed::get_header(ifstream *file){
    string line;
    //removing the header and storing it in vcf.header
    while (getline(*file, line) && line[0]=='#' && line[1]=='#'){
        header.append(line + '\n');
        header_size += line.length() + 1;
    }
    header_size += line.length() + 1;
    variants_size = filesize - header_size; // New size without the header
}
    
/**
    * @brief Prints the VCF header to standard output.
    */
void vcf_parsed::print_header(){
    cout << "VCF header:\n" << header << endl;
}
    
/**
    * @brief Reads and parses the VCF header.
    *
    * Separates header and variant data, extracts INFO and FORMAT information,
    * and determines the number of samples present.
    *
    * @param file Pointer to the input file stream.
    */
void vcf_parsed::get_and_parse_header(ifstream *file){
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
    
    INFO.total_values = INFO.ID.size();
    INFO.no_alt_values = INFO.total_values - INFO.alt_values;

    header_size += line.length() + 1;

    variants_size = filesize - header_size; // New size without the header
}   
    
/**
    * @brief Allocates a character array to store the variant portion of the VCF file.
    *
    * The allocated size is based on the file size minus the header size.
    */
void vcf_parsed::allocate_filestring(){
    filestring = (char*)malloc(variants_size + 8);
    memset(filestring, '\0', variants_size + 8);
}

/**
    * @brief Creates and initializes vectors for sample data.
    *
    * Based on the FORMAT header, this method initializes vectors for sample genotype,
    * float, integer, and string data.
    *
    * @param num_threads Number of threads to use for parallel processing.
    */
void vcf_parsed::create_sample_vectors(int num_threads){
    samp_Flag samp_flag_tmp;
    samp_Float samp_float_tmp;
    samp_Int samp_int_tmp;
    samp_String samp_string_tmp;

    samp_Float samp_alt_float_tmp;
    samp_Int samp_alt_int_tmp;
    samp_String samp_alt_string_tmp;

    if(FORMAT.hasGT && FORMAT.numGT == 'A'){
        alt_sample.initMapGT();
        alt_sample.sample_GT.numb = -1;
    }else if(FORMAT.hasGT){
        int iter = FORMAT.numGT - '0';
        samp_columns.initMapGT();
        for(int i=0; i<iter; i++){ 
            samp_GT tmp;
            tmp.numb = iter;
            tmp.GT.resize((num_lines)*samp_columns.numSample, (char)0);
            samp_columns.sample_GT.push_back(tmp);
        }
        samp_columns.sample_GT.resize(FORMAT.numGT-'0');   
    }

    int numIter = FORMAT.ID.size();
    if(numIter == 0 ) return;

    for(int i = 0; i < numIter; i++){
        if(strcmp(&FORMAT.Number[i][0], "A") != 0){
            // Without Alternatives
            if(strcmp(&FORMAT.Number[i][0], "1")==0){ 
                //Number = 1
                if(!strcmp(&FORMAT.Type[i][0], "String")){
                    samp_string_tmp.name = FORMAT.ID[i];
                    samp_columns.samp_string.push_back(samp_string_tmp);
                    samp_columns.samp_string.back().i_string.resize((num_lines)*samp_columns.numSample, "\0");
                    samp_columns.samp_string.back().numb = std::stoi(FORMAT.Number[i]);
                    info_map[FORMAT.ID[i]] = 8;
                    var_columns.info_map1[FORMAT.ID[i]] = 8;
                    FORMAT.strings++;                        
                }else if(!strcmp(&FORMAT.Type[i][0], "Integer")){
                    samp_int_tmp.name = FORMAT.ID[i];
                    samp_columns.samp_int.push_back(samp_int_tmp);
                    samp_columns.samp_int.back().i_int.resize((num_lines)*samp_columns.numSample, 0);
                    samp_columns.samp_int.back().numb = std::stoi(FORMAT.Number[i]);
                    info_map[FORMAT.ID[i]] = 9;
                    var_columns.info_map1[FORMAT.ID[i]] = 9;
                    FORMAT.ints++;
                }else if(!strcmp(&FORMAT.Type[i][0], "Float")){
                    samp_float_tmp.name = FORMAT.ID[i];
                    samp_columns.samp_float.push_back(samp_float_tmp);
                    samp_columns.samp_float.back().i_float.resize((num_lines)*samp_columns.numSample, 0);
                    samp_columns.samp_float.back().numb = std::stoi(FORMAT.Number[i]);
                    info_map[FORMAT.ID[i]] = 10;
                    var_columns.info_map1[FORMAT.ID[i]] = 10;
                    FORMAT.floats++;
                }
            }else if(strcmp(&FORMAT.Number[i][0], "0")==0){ 
                //Number = 0; so it's a flag
                samp_flag_tmp.name = FORMAT.ID[i];
                samp_columns.samp_flag.push_back(samp_flag_tmp);
                samp_columns.samp_flag.back().i_flag.resize((num_lines)*samp_columns.numSample, 0);
                samp_columns.samp_flag.back().numb = std::stoi(FORMAT.Number[i]);
                info_map[FORMAT.ID[i]] = 11;
                var_columns.info_map1[FORMAT.ID[i]] = 11;
                FORMAT.flags++;
            }else if(strcmp(&FORMAT.Number[i][0], ".") == 0){
                int userNumber = -1;
                // Ciclo per richiedere un input corretto (intero non negativo)
                while(true) {
                    std::cout << "Select the number for the field " 
                            << FORMAT.ID[i] << " (non negative integer): ";
                    if (std::cin >> userNumber && userNumber >= 0) {
                        break; // Input corretto, esce dal ciclo
                    } else {
                        std::cout << "Valore non valido. Riprova." << std::endl;
                        std::cin.clear(); // Ripristina lo stato del flusso
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Scarta l'input errato
                    }
                }
            
                // Caso in cui l'utente inserisca 1: trattiamo il campo come un singolo elemento
                if(userNumber == 1){
                    if(!strcmp(&FORMAT.Type[i][0], "String")){
                        samp_string_tmp.name = FORMAT.ID[i];
                        samp_columns.samp_string.push_back(samp_string_tmp);
                        samp_columns.samp_string.back().i_string.resize((num_lines)*samp_columns.numSample, "\0");
                        samp_columns.samp_string.back().numb = userNumber;
                        info_map[FORMAT.ID[i]] = 8;
                        var_columns.info_map1[FORMAT.ID[i]] = 8;
                        FORMAT.strings++;                        
                    } else if(!strcmp(&FORMAT.Type[i][0], "Integer")){
                        samp_int_tmp.name = FORMAT.ID[i];
                        samp_columns.samp_int.push_back(samp_int_tmp);
                        samp_columns.samp_int.back().i_int.resize((num_lines)*samp_columns.numSample, 0);
                        samp_columns.samp_int.back().numb = userNumber;
                        info_map[FORMAT.ID[i]] = 9;
                        var_columns.info_map1[FORMAT.ID[i]] = 9;
                        FORMAT.ints++;
                    } else if(!strcmp(&FORMAT.Type[i][0], "Float")){
                        samp_float_tmp.name = FORMAT.ID[i];
                        samp_columns.samp_float.push_back(samp_float_tmp);
                        samp_columns.samp_float.back().i_float.resize((num_lines)*samp_columns.numSample, 0);
                        samp_columns.samp_float.back().numb = userNumber;
                        info_map[FORMAT.ID[i]] = 10;
                        var_columns.info_map1[FORMAT.ID[i]] = 10;
                        FORMAT.floats++;
                    }
                }else if(userNumber == 0){
                    samp_flag_tmp.name = FORMAT.ID[i];
                    samp_columns.samp_flag.push_back(samp_flag_tmp);
                    samp_columns.samp_flag.back().i_flag.resize((num_lines)*samp_columns.numSample, 0);
                    samp_columns.samp_flag.back().numb = userNumber;
                    info_map[FORMAT.ID[i]] = 11;
                    var_columns.info_map1[FORMAT.ID[i]] = 11;
                    FORMAT.flags++;
                }else { // userNumber>1
                    if(!strcmp(&FORMAT.Type[i][0], "String")){
                        for(int j = 0; j < userNumber; j++){
                            samp_string_tmp.name = FORMAT.ID[i] + std::to_string(j);
                            samp_columns.samp_string.push_back(samp_string_tmp);
                            samp_columns.samp_string.back().i_string.resize((num_lines)*samp_columns.numSample, "\0");
                            samp_columns.samp_string.back().numb = userNumber;
                            info_map[FORMAT.ID[i] + std::to_string(j)] = 8;
                            var_columns.info_map1[FORMAT.ID[i] + std::to_string(j)] = 8;
                            FORMAT.strings++;
                        }
                    } else if(!strcmp(&FORMAT.Type[i][0], "Integer")){
                        for(int j = 0; j < userNumber; j++){
                            samp_int_tmp.name = FORMAT.ID[i] + std::to_string(j);
                            samp_columns.samp_int.push_back(samp_int_tmp);
                            samp_columns.samp_int.back().i_int.resize((num_lines)*samp_columns.numSample, 0);
                            samp_columns.samp_int.back().numb = userNumber;
                            info_map[FORMAT.ID[i] + std::to_string(j)] = 9;
                            var_columns.info_map1[FORMAT.ID[i] + std::to_string(j)] = 9;
                            FORMAT.ints++;
                        }
                    } else if(!strcmp(&FORMAT.Type[i][0], "Float")){
                        for(int j = 0; j < userNumber; j++){
                            samp_float_tmp.name = FORMAT.ID[i] + std::to_string(j);
                            samp_columns.samp_float.push_back(samp_float_tmp);
                            samp_columns.samp_float.back().i_float.resize((num_lines)*samp_columns.numSample, 0);
                            samp_columns.samp_float.back().numb = userNumber;
                            info_map[FORMAT.ID[i] + std::to_string(j)] = 10;
                            var_columns.info_map1[FORMAT.ID[i] + std::to_string(j)] = 10;
                            FORMAT.floats++;
                        }
                    }
                }                    
            }else{ 
                //Number > 1
                if(!strcmp(&FORMAT.Type[i][0], "String")){
                    for(int j = 0; j < std::stoi(FORMAT.Number[i]); j++){
                        samp_string_tmp.name = FORMAT.ID[i] + std::to_string(j);
                        samp_columns.samp_string.push_back(samp_string_tmp);
                        samp_columns.samp_string.back().i_string.resize((num_lines)*samp_columns.numSample, "\0");
                        samp_columns.samp_string.back().numb = std::stoi(FORMAT.Number[i]);
                        info_map[FORMAT.ID[i]+std::to_string(j)] = 8;
                        var_columns.info_map1[FORMAT.ID[i]+std::to_string(j)] = 8;
                        FORMAT.strings++;
                    }
                }else if(!strcmp(&FORMAT.Type[i][0], "Integer")){
                    for(int j = 0; j < std::stoi(FORMAT.Number[i]); j++){
                        samp_int_tmp.name = FORMAT.ID[i] + std::to_string(j);
                        samp_columns.samp_int.push_back(samp_int_tmp);
                        samp_columns.samp_int.back().i_int.resize((num_lines)*samp_columns.numSample, 0);
                        samp_columns.samp_int.back().numb = std::stoi(FORMAT.Number[i]);
                        info_map[FORMAT.ID[i]+std::to_string(j)] = 9;
                        var_columns.info_map1[FORMAT.ID[i]+std::to_string(j)] = 9;
                        FORMAT.ints++;
                    }
                }else if(!strcmp(&FORMAT.Type[i][0], "Float")){
                    for(int j = 0; j < std::stoi(FORMAT.Number[i]); j++){
                        samp_float_tmp.name = FORMAT.ID[i] + std::to_string(j);
                        samp_columns.samp_float.push_back(samp_float_tmp);
                        samp_columns.samp_float.back().i_float.resize((num_lines)*samp_columns.numSample, 0);
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
        samp_columns.var_id.resize((num_lines)*samp_columns.numSample, 0);
        samp_columns.samp_id.resize((num_lines)*samp_columns.numSample, static_cast<unsigned short>(0));
    }    
    alt_sample.samp_flag.resize(FORMAT.flags_alt);
    alt_sample.samp_int.resize(FORMAT.ints_alt);
    alt_sample.samp_float.resize(FORMAT.floats_alt);
    alt_sample.samp_string.resize(FORMAT.strings_alt);
    alt_sample.var_id.resize((num_lines)* alt_sample.numSample, 0);
    alt_sample.samp_id.resize((num_lines)*alt_sample.numSample, static_cast<unsigned short>(0));

}
    
/**
    * @brief Creates and initializes vectors for INFO field data.
    *
    * Initializes vectors to store INFO field values (flag, integer, float, and string)
    * based on header information.
    *
    * @param num_threads Number of threads to use for parallel processing.
    */
void vcf_parsed::create_info_vectors(int num_threads){
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
                //TODO - implementazione temporanea start - TSA - per provare la metto negli int deterministici
                if(strcmp(INFO.ID[i].c_str(), "TSA")==0){
                    INFO.ints++;
                    info_int_tmp.name = INFO.ID[i];
                    info_int_tmp.i_int.resize(num_lines-1, 0);
                    var_columns.in_int.push_back(info_int_tmp);
                    info_map[INFO.ID[i]] = 1;
                    var_columns.info_map1[INFO.ID[i]] = 1;
                }else{ //TODO - implementazione temporanea end
                    INFO.strings++;
                    info_string_tmp.name = INFO.ID[i];
                    info_string_tmp.i_string.resize(num_lines-1, "\0");
                    var_columns.in_string.push_back(info_string_tmp);
                    info_map[INFO.ID[i]] = 3;
                    var_columns.info_map1[INFO.ID[i]] = 3;
                }
            } else if(strcmp(&INFO.Type[i][0], "Flag")==0){
                INFO.flags++;
                info_flag_tmp.name = INFO.ID[i];
                info_flag_tmp.i_flag.resize(num_lines-1, 0);
                var_columns.in_flag.push_back(info_flag_tmp);
                info_map[INFO.ID[i]] = 0;
                var_columns.info_map1[INFO.ID[i]] = 0;
            }
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
    
/**
    * @brief Prints the INFO field mapping.
    *
    * Outputs the mapping from INFO field names to their corresponding type codes.
    */
void vcf_parsed::print_info_map(){
    for(const auto& element : info_map){
        cout<<element.first<<": "<<element.second<<endl;
    }
}
    
/**
    * @brief Prints a summary of INFO field data.
    *
    * Displays a brief summary of the sizes and first few entries for each INFO field type.
    */
void vcf_parsed::print_info(){
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
}
    
/**
    * @brief Reserves space in the variant columns structure.
    *
    * Resizes the vectors in the var_columns_df structure based on the number of variants.
    */
void vcf_parsed::reserve_var_columns(){
    var_columns.var_number.resize(num_lines-1);
    var_columns.chrom.resize(num_lines-1);
    var_columns.id.resize(num_lines-1);
    var_columns.pos.resize(num_lines-1);
    var_columns.ref.resize(num_lines-1); 
    var_columns.qual.resize(num_lines-1);
    var_columns.filter.resize(num_lines-1);
}

/**
    * @brief Allocates device memory for KernelParams and copies the host structure.
    *
    * @param d_params Double pointer to the device KernelParams.
    * @param h_params Pointer to the host KernelParams structure.
    */
void vcf_parsed::allocParamPointers(KernelParams **d_params, KernelParams *h_params) {
    // Allocate memory for KernelParams on GPU
    cudaMalloc((void**)d_params, sizeof(KernelParams));

    // Copy the structure from host to device
    cudaMemcpy(*d_params, h_params, sizeof(KernelParams), cudaMemcpyHostToDevice);
}

/**
    * @brief Launches the CUDA kernel to parse VCF lines and merges the results.
    *
    * Sets up CUDA streams and events, launches the appropriate kernel (based on whether sample data is present),
    * and asynchronously copies the parsed data from device to host.
    */
void vcf_parsed::populate_runner(int numb_cores){
    int threadsPerBlock = 32;
    int blocksPerGrid = (numb_cores/threadsPerBlock) + 1; // TODO gestire dinamicamente in funzione del numero di core
    cudaEvent_t kernel_done;
    cudaEventCreate(&kernel_done);
    auto start = chrono::system_clock::now();
    char* my_mem;
    int batchSize = threadsPerBlock*blocksPerGrid;
    cudaMalloc(&my_mem, batchSize*MAX_TOKEN_LEN*MAX_TOKENS*3);

    cudaStream_t stream1, stream2;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);

    // Pass existing device pointers to h_params
    h_params.line = d_filestring;
    h_params.var_number = d_VC_var_number;
    h_params.pos = d_VC_pos;
    h_params.qual = d_VC_qual;
    h_params.in_float = d_VC_in_float->i_float;
    h_params.in_flag = d_VC_in_flag->i_flag;
    h_params.in_int = d_VC_in_int->i_int;
    h_params.float_name = d_VC_in_float->name;
    h_params.flag_name = d_VC_in_flag->name;
    h_params.int_name = d_VC_in_int->name;
    h_params.new_lines_index = d_new_lines_index;
    h_params.numLines = num_lines;

    if(samplesON){            
        
        h_params.samp_var_id = d_SC_var_id;
        h_params.samp_id = d_SC_samp_id;
        h_params.samp_float = d_SC_samp_float->i_float;
        h_params.samp_flag = d_SC_samp_flag->i_flag;
        h_params.samp_int = d_SC_samp_int->i_int;
        h_params.samp_float_name = d_SC_samp_float->name;
        h_params.samp_flag_name = d_SC_samp_flag->name;
        h_params.samp_int_name = d_SC_samp_int->name;
        h_params.samp_float_numb = d_SC_samp_float->numb;
        h_params.samp_flag_numb = d_SC_samp_flag->numb;
        h_params.samp_int_numb = d_SC_samp_int->numb;
        h_params.sample_GT = d_SC_sample_GT->GT;
        h_params.numSample = samp_columns.numSample;
        h_params.numGT = (int)(FORMAT.numGT - '0');

        // Allocate d_params and copy h_params to GPU
        allocParamPointers(&d_params, &h_params);
        //TODO - Da rimuovere questa sinchronize utile per debugging:
        //cudaError_t err = cudaDeviceSynchronize();
        //if (err != cudaSuccess) {
        //    fprintf(stderr, "CUDA error in param alloc: %s\n", cudaGetErrorString(err));
        //    exit(EXIT_FAILURE);
        //}

        // Launch kernel
        kernel<<<blocksPerGrid, threadsPerBlock, 0, stream1>>>(d_params, my_mem, batchSize, true);

        cudaEventRecord(kernel_done, stream1);
        // Check for errors
        //err = cudaGetLastError();
        //if (err != cudaSuccess) {
        //    fprintf(stderr, "CUDA error in kernel launch: %s\n", cudaGetErrorString(err));
        //    exit(EXIT_FAILURE);
        //}

        //TODO - Da rimuovere questa sinchronize utile per debugging:
        //err = cudaDeviceSynchronize();
        //if (err != cudaSuccess) {
        //    fprintf(stderr, "CUDA error in kernel execution: %s\n", cudaGetErrorString(err));
        //    exit(EXIT_FAILURE);
        //}
    }else{
        // Allocate d_params and copy h_params to GPU
        allocParamPointers(&d_params, &h_params);
        kernel<<<blocksPerGrid, threadsPerBlock, 0, stream1>>>(d_params, my_mem, batchSize, false);
        cudaEventRecord(kernel_done, stream1);
        
        // Controllo errori di lancio del kernel
        cudaError_t err = cudaGetLastError();
        if(err != cudaSuccess) {
            fprintf(stderr, "Errore nel lancio di get_vcf_line_kernel: %s\n", cudaGetErrorString(err));
            exit(EXIT_FAILURE);
        }
        
    }
    
    //TODO - memcpyAsinc may be broken 

    cudaStreamWaitEvent(stream2, kernel_done, 0);
    cudaMemcpyAsync(var_columns.var_number.data(), d_VC_var_number, (num_lines) * sizeof(unsigned int), cudaMemcpyDeviceToHost, stream2);
    cudaMemcpyAsync(var_columns.pos.data(), d_VC_pos, (num_lines) * sizeof(unsigned int), cudaMemcpyDeviceToHost, stream2);
    cudaMemcpyAsync(var_columns.qual.data(), d_VC_qual, (num_lines) * sizeof(__half), cudaMemcpyDeviceToHost, stream2);

    for(int i=0; i<var_columns.in_float.size(); i++){
        cudaMemcpyAsync(var_columns.in_float[i].i_float.data(), d_VC_in_float->i_float + i * (num_lines), (num_lines)*sizeof(__half), cudaMemcpyDeviceToHost, stream2);
    }

    for(int i=0; i<var_columns.in_flag.size(); i++){
        cudaMemcpyAsync(var_columns.in_flag[i].i_flag.data(), d_VC_in_flag->i_flag + i * (num_lines), (num_lines)*sizeof(bool), cudaMemcpyDeviceToHost, stream2);
    }

    for(int i=0; i<var_columns.in_int.size(); i++){
        cudaMemcpyAsync(var_columns.in_int[i].i_int.data(), d_VC_in_int->i_int + i * (num_lines), (num_lines) * sizeof(int), cudaMemcpyDeviceToHost, stream2);
    }

    if(samplesON){
        cudaMemcpyAsync(samp_columns.var_id.data(), d_SC_var_id, (num_lines) * samp_columns.numSample * sizeof(unsigned int), cudaMemcpyDeviceToHost, stream2);
        cudaMemcpyAsync(samp_columns.samp_id.data(), d_SC_samp_id, (num_lines) * samp_columns.numSample * sizeof(unsigned short), cudaMemcpyDeviceToHost, stream2);

        for (int i = 0; i < samp_columns.samp_float.size(); i++) {
            cudaMemcpyAsync(samp_columns.samp_float[i].i_float.data(), d_SC_samp_float->i_float + i * ((num_lines) * samp_columns.numSample), 
                (num_lines) * samp_columns.numSample * sizeof(__half), cudaMemcpyDeviceToHost, stream2);
        }

        for (int i = 0; i < samp_columns.samp_flag.size(); i++) {
            cudaMemcpyAsync(samp_columns.samp_flag[i].i_flag.data(), d_SC_samp_flag->i_flag + i * ((num_lines) * samp_columns.numSample), 
                (num_lines) * samp_columns.numSample * sizeof(bool), cudaMemcpyDeviceToHost, stream2);
        }

        for (int i = 0; i < samp_columns.samp_int.size(); i++) {
            cudaMemcpyAsync(samp_columns.samp_int[i].i_int.data(), d_SC_samp_int->i_int + (i * num_lines * samp_columns.numSample), 
                (num_lines) * samp_columns.numSample * sizeof(int), cudaMemcpyDeviceToHost, stream2);                
        }   

        for(int i=0; i<samp_columns.sample_GT.size(); i++){
            cudaMemcpyAsync(samp_columns.sample_GT[i].GT.data(), d_SC_sample_GT->GT + i * ((num_lines) * samp_columns.numSample), 
                (num_lines)*samp_columns.numSample*sizeof(char), cudaMemcpyDeviceToHost, stream2);
        } 
    }
    cudaStreamSynchronize(stream2);

    // Cleanup
    cudaEventDestroy(kernel_done);
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaFree(my_mem);
}

/**
    * @brief Populates variant columns by processing VCF lines in parallel.
    *
    * Spawns a worker thread to run the CUDA kernel for parsing and uses OpenMP to merge alternative allele
    * data from multiple threads into the final data structures (alt_columns_df and alt_format_df).
    *
    * @param num_threads Number of threads to use for parallel merging.
    */
void vcf_parsed::populate_var_columns(int num_threads, int numb_cores){

    std::thread worker_thread(&vcf_parsed::populate_runner, this, numb_cores);

    long batch_size = (num_lines-2+num_threads)/num_threads;
    
    std::vector<alt_columns_df> tmp_alt(num_threads);
    std::vector<int> tmp_num_alt(num_threads);

    std::vector<alt_format_df> tmp_alt_format(num_threads);
    std::vector<int> tmp_num_alt_format(num_threads);

    
    int totAlt = 0;
    int totSampAlt = 0;

    #pragma omp parallel //TODO lentissimo
    {
        long start, end;
        int th_ID = omp_get_thread_num();
        // Temporary structure of the thread with alternatives.
        tmp_alt[th_ID].init(alt_columns, INFO, batch_size);
        
        tmp_num_alt[th_ID] = 0;

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
            for(long i=start; i<end && i<num_lines-1; i++){ //TODO - da sistemare
                get_vcf_line_in_var_columns_format(filestring, new_lines_index[i], new_lines_index[i+1], i, &(tmp_alt[th_ID]), &(tmp_num_alt[th_ID]), &samp_columns, &FORMAT, &(tmp_num_alt_format[th_ID]), &(tmp_alt_format[th_ID]));
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
            for(long i=start; i<end && i<num_lines-1; i++){ //TODO - Da sistemare
                get_vcf_line_in_var_columns(filestring, new_lines_index[i], new_lines_index[i+1], i, &(tmp_alt[th_ID]), &(tmp_num_alt[th_ID]));
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

    // Merge results in parallel
    std::thread t1(merge_member_vector<alt_columns_df, unsigned int>, std::ref(tmp_alt),
            std::ref(alt_columns.var_id), num_threads, &alt_columns_df::var_id);

    std::thread t2(merge_member_vector<alt_columns_df, unsigned char>, std::ref(tmp_alt),
                std::ref(alt_columns.alt_id), num_threads, &alt_columns_df::alt_id);
    
    std::thread t3(merge_member_vector<alt_columns_df, string>, std::ref(tmp_alt),
                std::ref(alt_columns.alt), num_threads, &alt_columns_df::alt);

    std::thread t4(merge_nested_member_vector<alt_columns_df, info_int, int>, 
        std::ref(tmp_alt), std::ref(alt_columns.alt_int), num_threads, INFO.ints_alt, &alt_columns_df::alt_int, &info_int::i_int);

    std::thread t5(merge_nested_member_vector<alt_columns_df, info_float, __half>, 
        std::ref(tmp_alt), std::ref(alt_columns.alt_float), num_threads, INFO.floats_alt, &alt_columns_df::alt_float, &info_float::i_float);

    std::thread t6(merge_nested_member_vector<alt_columns_df, info_string, string>, 
        std::ref(tmp_alt), std::ref(alt_columns.alt_string), num_threads, INFO.strings_alt, &alt_columns_df::alt_string, &info_string::i_string);

    std::thread t_sum([&]() {
        int somma = 0;
        for (int i = 0; i < num_threads; i++) {
            somma += tmp_num_alt[i];
        }
        totAlt = somma;
    });

    if (samplesON) {
        std::thread t7(merge_member_vector<alt_format_df, unsigned int>, std::ref(tmp_alt_format),
                std::ref(alt_sample.var_id), num_threads, &alt_format_df::var_id);

        std::thread t8(merge_member_vector<alt_format_df, char>, std::ref(tmp_alt_format),
            std::ref(alt_sample.alt_id), num_threads, &alt_format_df::alt_id);

        std::thread t9(merge_member_vector<alt_format_df, unsigned short>, std::ref(tmp_alt_format),
                std::ref(alt_sample.samp_id), num_threads, &alt_format_df::samp_id);

        std::thread t10(merge_nested_member_vector<alt_format_df, samp_String, string>, std::ref(tmp_alt_format),
                std::ref(alt_sample.samp_string), num_threads, FORMAT.strings_alt, &alt_format_df::samp_string, 
                &samp_String::i_string);
        
        std::thread t11(merge_nested_member_vector<alt_format_df, samp_Int, int>, std::ref(tmp_alt_format),
                std::ref(alt_sample.samp_int), num_threads, FORMAT.ints_alt, &alt_format_df::samp_int, 
                &samp_Int::i_int);

        std::thread t12(merge_nested_member_vector<alt_format_df, samp_Float, __half>, std::ref(tmp_alt_format),
                std::ref(alt_sample.samp_float), num_threads, FORMAT.floats_alt, &alt_format_df::samp_float, 
                &samp_Float::i_float);

        std::thread t_sum_samp([&]() {
            int somma = 0;
            for (int i = 0; i < num_threads; i++) {
                somma += tmp_num_alt_format[i];
            }
            totSampAlt = somma;
        });

        t7.join();
        t8.join();
        t9.join();
        t10.join();
        t11.join();
        t12.join();
        t_sum_samp.join();
    }

    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();
    t_sum.join();

    //Here finish the parallel part and we merge the threads results

    alt_columns.numAlt = totAlt;
    alt_sample.numSample = totSampAlt;
    // Eseguiamo il resize in parallelo per alt_columns
    {
        // Task per ridimensionare le vector "piatte"
        auto fut1 = std::async(std::launch::async, [&]() {
            alt_columns.var_id.resize(totAlt);
            alt_columns.alt_id.resize(totAlt);
            alt_columns.alt.resize(totAlt);
        });
        
        // Task per ridimensionare le vector interne di alt_int
        auto fut2 = std::async(std::launch::async, [&]() {
            for (int j = 0; j < INFO.ints_alt; j++) {
                alt_columns.alt_int[j].i_int.resize(totAlt);
            }
        });
        
        // Task per ridimensionare le vector interne di alt_float
        auto fut3 = std::async(std::launch::async, [&]() {
            for (int j = 0; j < INFO.floats_alt; j++) {
                alt_columns.alt_float[j].i_float.resize(totAlt);
            }
        });
        
        // Task per ridimensionare le vector interne di alt_string
        auto fut4 = std::async(std::launch::async, [&]() {
            for (int j = 0; j < INFO.strings_alt; j++) {
                alt_columns.alt_string[j].i_string.resize(totAlt);
            }
        });
        
        // Aspettiamo che tutti i task completino
        fut1.get();
        fut2.get();
        fut3.get();
        fut4.get();
    }

    // Se samplesON è attivo, facciamo la stessa cosa per alt_sample
    if (samplesON) {
        auto fut1 = std::async(std::launch::async, [&]() {
            alt_sample.var_id.resize(totSampAlt);
            alt_sample.samp_id.resize(totSampAlt);
            alt_sample.alt_id.resize(totSampAlt);
        });
        
        auto fut2 = std::async(std::launch::async, [&]() {
            for (int j = 0; j < FORMAT.ints_alt; j++) {
                alt_sample.samp_int[j].i_int.resize(totSampAlt);
            }
        });
        
        auto fut3 = std::async(std::launch::async, [&]() {
            for (int j = 0; j < FORMAT.floats_alt; j++) {
                alt_sample.samp_float[j].i_float.resize(totSampAlt);
            }
        });
        
        auto fut4 = std::async(std::launch::async, [&]() {
            for (int j = 0; j < FORMAT.strings_alt; j++) {
                alt_sample.samp_string[j].i_string.resize(totSampAlt);
            }
        });
        
        fut1.get();
        fut2.get();
        fut3.get();
        fut4.get();
    }
    worker_thread.join();    
}

/**
    * @brief Parses a VCF line and populates variant columns data.
    *
    * This function processes a single VCF line (from index @p start to @p end) by reading it character by character.
    * It extracts key variant fields such as chromosome, variant ID, reference allele, alternative alleles,
    * filter information, and the INFO field. The alternative allele field is split using commas and stored in
    * the provided alt_columns_df structure (@p tmp_alt). The count of alternative alleles processed is tracked by
    * the integer pointed to by @p tmp_num_alt.
    *
    * @param line Pointer to the VCF line as a C-string.
    * @param start The starting index of the line within the file.
    * @param end The ending index of the line within the file.
    * @param i The index (row number) corresponding to the current variant.
    * @param tmp_alt Pointer to an alt_columns_df structure for storing alternative allele data.
    * @param tmp_num_alt Pointer to an integer tracking the current number of alternative alleles processed.
    */
void vcf_parsed::get_vcf_line_in_var_columns(char *line, long start, long end, long i, alt_columns_df* tmp_alt, int *tmp_num_alt)
{ 
    bool find1 = false;
    long iter=0;
    int local_alt = 1;
    string tmp="\0";
    vector<string> tmp_split;
    vector<string> tmp_format_split;

    if(line[start+iter]=='\n'){
        iter++;
    } 

    //Chromosome
    while(!find1){
        if(line[start+iter]=='\t'||line[start+iter]==' '){
            find1 = true;
            iter++;
            if(var_columns.chrom_map.find(tmp) == var_columns.chrom_map.end()){
                var_columns.chrom_map.insert(std::make_pair(tmp, (unsigned char)var_columns.chrom_map.size()));
            }
            var_columns.chrom[i] = var_columns.chrom_map[tmp];
        }else{
            tmp += line[start+iter];
            iter++;
        }
    }

    //Position on device
    // Salta la sottostringa delimitata da '\t' o ' '
    while (line[start + iter] != '\t' && line[start + iter] != ' ') {
        iter++;
    }
    iter++; // Salta anche il delimitatore

    //ID
    tmp="\0";
    find1=false;
    while(!find1){
        if(line[start+iter]=='\t'||line[start+iter]==' '){
            find1 = true;
            iter++;
            var_columns.id[i] = tmp;
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
            var_columns.ref[i] = tmp;
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
                (*tmp_alt).alt[(*tmp_num_alt)] = tmp_split[y];
                (*tmp_alt).alt_id[(*tmp_num_alt)] = (char)y;
                (*tmp_alt).var_id[(*tmp_num_alt)] = i;
            }
        }else{
            tmp += line[start+iter];
            iter++;
        }
    }

    //Quality - on device
    while (line[start + iter] != '\t' && line[start + iter] != ' ') {
        iter++;
    }
    iter++;
    
    //Filter
    tmp="\0";
    find1=false;
    while(!find1){
        if(line[start+iter]=='\t'||line[start+iter]==' '){
            find1 = true;
            iter++;
            if(var_columns.filter_map.find(tmp) == var_columns.filter_map.end()){
                var_columns.filter_map.insert(std::make_pair(tmp, (char)var_columns.filter_map.size()));
            }
            var_columns.filter[i] = var_columns.filter_map[tmp];
        }else{
            tmp += line[start+iter];
            iter++;
        }
    }

    //Info
    tmp="\0";
    find1=false;
    bool find_info_type = false;
    bool find_info_elem = false;
    int el=0;
    while(!find1){
        if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
            find1 = true;
            iter++;
            vector<string> tmp_el;
            boost::split(tmp_el, tmp, boost::is_any_of(";")); //Info arguments separation
            vector<string> tmp_elems;
            for(int ii=0; ii<tmp_el.size(); ii++){
                boost::split(tmp_elems, tmp_el[ii], boost::is_any_of("=")); //info_id separation from contents
                find_info_type = false;
                find_info_elem = false;
                if(tmp_elems.size()==2){
                    while(!find_info_type){
                        if(var_columns.info_map1[tmp_elems[0]]==STRING){
                            //String
                            el=0;
                            while(!find_info_elem){
                                if(var_columns.in_string[el].name == tmp_elems[0]){
                                    var_columns.in_string[el].i_string[i] = tmp_elems[1];
                                    find_info_elem = true;
                                } 
                                el++;
                            }
                            find_info_type = true;
                        }else if(var_columns.info_map1[tmp_elems[0]]==INT_ALT){
                            //Int Alt
                            el=0;
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
                        }else if(var_columns.info_map1[tmp_elems[0]]==FLOAT_ALT){
                            //Float Alt
                            el=0;
                            while(!find_info_elem){
                                if((*tmp_alt).alt_float[el].name == tmp_elems[0]){
                                    boost::split(tmp_split, tmp_elems[1], boost::is_any_of(","));
                                    
                                    for(int y = 0; y<local_alt; y++){
                                        try{
                                            (*tmp_alt).alt_float[el].i_float[(*tmp_num_alt)+y] = (__half)stof(tmp_split[y]);
                                        }catch (const std::exception& e){
                                            (*tmp_alt).alt_float[el].i_float[(*tmp_num_alt)+y] = 0;
                                        }
                                    }
                                    find_info_elem = true;
                                }
                                el++;
                            }
                            find_info_type = true;
                        }else if(var_columns.info_map1[tmp_elems[0]]==STRING_ALT){
                            //String Alt
                            el=0;
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
                } // la flag sarà su device
            }
        }else{
            tmp += line[start+iter];
            iter++;
        }
    }
    (*tmp_num_alt) = (*tmp_num_alt)+local_alt;
}

/**
    * @brief Parses a VCF line with formatted variant and sample data.
    *
    * This function processes a single VCF line to extract both variant and sample-related information,
    * following a predefined FORMAT template. It parses the chromosome, variant ID, reference allele,
    * alternative alleles, and filter field, and then further splits the FORMAT field to extract per-sample
    * data (e.g., genotype, float, integer, and string values). The extracted sample data is stored in the
    * provided alt_format_df structure (@p tmp_alt_format), while variant data is updated in the global structures.
    *
    * @param line Pointer to the VCF line as a C-string.
    * @param start The starting index of the line within the file.
    * @param end The ending index of the line within the file.
    * @param i The index (row number) corresponding to the current variant.
    * @param tmp_alt Pointer to an alt_columns_df structure for storing alternative allele data.
    * @param tmp_num_alt Pointer to an integer tracking the number of alternative alleles processed.
    * @param sample Pointer to a sample_columns_df structure for storing sample-specific data.
    * @param FORMAT Pointer to a header_element structure describing the FORMAT fields.
    * @param tmp_num_alt_format Pointer to an integer tracking the number of formatted alternative entries processed.
    * @param tmp_alt_format Pointer to an alt_format_df structure for storing formatted sample data.
    */
void vcf_parsed::get_vcf_line_in_var_columns_format(char *line, long start, long end, long i, alt_columns_df* tmp_alt, int *tmp_num_alt, sample_columns_df* sample, header_element* FORMAT, int *tmp_num_alt_format, alt_format_df* tmp_alt_format)
{
    bool find1 = false;
    long iter=0;
    int local_alt = 1;
    string tmp="\0";
    vector<string> tmp_split;
    vector<string> tmp_format_split;
    vector<string> tmp_subSplit;

    if(line[start+iter]=='\n'){
        iter++;
    } 
    //Chromosome
    while(!find1){
        if(line[start+iter]=='\t'||line[start+iter]==' '){
            find1 = true;
            iter++;
            if(var_columns.chrom_map.find(tmp) == var_columns.chrom_map.end()){
                var_columns.chrom_map.insert(std::make_pair(tmp, (unsigned char)var_columns.chrom_map.size()));
            }
            var_columns.chrom[i] = var_columns.chrom_map[tmp];
        }else{
            tmp += line[start+iter];
            iter++;
        }
    }

    //Position on device
    // Salta la sottostringa delimitata da '\t' o ' '
    while (line[start + iter] != '\t' && line[start + iter] != ' ') {
        iter++;
    }
    iter++; // Salta anche il delimitatore

    //ID
    tmp="\0";
    find1=false;
    while(!find1){
        if(line[start+iter]=='\t'||line[start+iter]==' '){
            find1 = true;
            iter++;
            var_columns.id[i] = tmp;
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
            var_columns.ref[i] = tmp;
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
                (*tmp_alt).alt[(*tmp_num_alt) + y] = tmp_split[y];
                (*tmp_alt).alt_id[(*tmp_num_alt) + y] = (char)y;
                (*tmp_alt).var_id[(*tmp_num_alt) + y] = i;
            }
        }else{
            tmp += line[start+iter];
            iter++;
        }
    }
    
    //Quality - on device
    while (line[start + iter] != '\t' && line[start + iter] != ' ') {
        iter++;
    }
    iter++;

    //Filter
    tmp="\0";
    find1=false;
    while(!find1){
        if(line[start+iter]=='\t'||line[start+iter]==' '){
            find1 = true;
            iter++;
            if(var_columns.filter_map.find(tmp) == var_columns.filter_map.end()){
                var_columns.filter_map.insert(std::make_pair(tmp, (char)var_columns.filter_map.size()));
            }
            var_columns.filter[i] = var_columns.filter_map[tmp];
        }else{
            tmp += line[start+iter];
            iter++;
        }
    }

    //Info
    tmp="\0";
    find1=false;
    bool find_info_type = false;
    bool find_info_elem = false;
    int el=0;
    while(!find1){
        if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
            find1 = true;
            iter++;
            vector<string> tmp_el;
            boost::split(tmp_el, tmp, boost::is_any_of(";")); //Info arguments separation
            vector<string> tmp_elems;
            for(int ii=0; ii<tmp_el.size(); ii++){
                boost::split(tmp_elems, tmp_el[ii], boost::is_any_of("=")); //info_id separation from contents
                find_info_type = false;
                find_info_elem = false;
                if(tmp_elems.size()==2){
                    while(!find_info_type){
                        if(var_columns.info_map1[tmp_elems[0]]==STRING){
                            //String
                            el=0;
                            while(!find_info_elem){
                                if(var_columns.in_string[el].name == tmp_elems[0]){
                                    var_columns.in_string[el].i_string[i] = tmp_elems[1];
                                    find_info_elem = true;
                                } 
                                el++;
                            }
                            find_info_type = true;
                        }else if(var_columns.info_map1[tmp_elems[0]]==INT_ALT){
                            //Int Alt
                            el=0;
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
                        }else if(var_columns.info_map1[tmp_elems[0]]==FLOAT_ALT){
                            //Float Alt
                            el=0;
                            while(!find_info_elem){
                                if((*tmp_alt).alt_float[el].name == tmp_elems[0]){
                                    boost::split(tmp_split, tmp_elems[1], boost::is_any_of(","));
                                    
                                    for(int y = 0; y<local_alt; y++){
                                        try{
                                            (*tmp_alt).alt_float[el].i_float[(*tmp_num_alt)+y] = (__half)stof(tmp_split[y]);
                                        }catch (const std::exception& e){
                                            (*tmp_alt).alt_float[el].i_float[(*tmp_num_alt)+y] = 0;
                                        }
                                    }
                                    find_info_elem = true;
                                }
                                el++;
                            }
                            find_info_type = true;
                        }else if(var_columns.info_map1[tmp_elems[0]]==STRING_ALT){
                            //String Alt
                            el=0;
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
                } // la flag sarà su device
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

    int samp;
    bool find_type = false;
    bool find_elem = false;
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
                    find_type = false;
                    find_elem = false;
                    while(!find_type){
                        if(!strcmp(tmp_format_split[j].c_str(), "GT")){
                            if(!((*sample).sample_GT.size() >= 1)){
                                boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                                local_alt = tmp_sub.size();
                                for(int y = 0; y<local_alt; y++){
                                    //Fill a tuple for each alternatives
                                    (*tmp_alt_format).var_id[(*tmp_num_alt_format) + y] = var_columns.var_number[i];
                                    (*tmp_alt_format).samp_id[(*tmp_num_alt_format) + y] = samp;
                                    (*tmp_alt_format).alt_id[(*tmp_num_alt_format) + y] = (char)y;
                                    (*tmp_alt_format).sample_GT.GT[(*tmp_num_alt_format) + y] = (*tmp_alt_format).GTMap[tmp_sub[y]];
                                }
                                (*tmp_num_alt_format) = (*tmp_num_alt_format) + local_alt;
                            }
                            find_type = true;
                        }else if(var_columns.info_map1[tmp_format_split[j]] == STRING_FORMAT || var_columns.info_map1[tmp_format_split[j] + std::to_string(1)] == STRING_FORMAT){
                            //String - deterministic
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
                        }else if(var_columns.info_map1[tmp_format_split[j]] == INT_FORMAT || var_columns.info_map1[tmp_format_split[j] + std::to_string(1)] == INT_FORMAT){
                            //Integer - deterministic - on device
                            find_elem = true;
                            find_type = true;
                        }else if(var_columns.info_map1[tmp_format_split[j]] == FLOAT_FORMAT || var_columns.info_map1[tmp_format_split[j] + std::to_string(1)] == FLOAT_FORMAT){
                            //Float - deterministic - on device
                            find_elem = true;
                            find_type = true;
                        }else if(var_columns.info_map1[tmp_format_split[j]] == STRING_FORMAT_ALT){
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
                                        (*tmp_alt_format).var_id[(*tmp_num_alt_format) + y] = var_columns.var_number[i];
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
                        }else if(var_columns.info_map1[tmp_format_split[j]] == INT_FORMAT_ALT){
                            //Integer alternatives
                            boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                            local_alt = tmp_sub.size();
                            int el = 0;
                            while(!find_elem){
                                //Search the corresponding element
                                if(!(*tmp_alt_format).samp_int[el].name.compare(tmp_format_split[j])){
                                    //Fill a tuple for each alternatives
                                    for(int y = 0; y<local_alt; y++){
                                        (*tmp_alt_format).var_id[(*tmp_num_alt_format) + y] = var_columns.var_number[i];
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
                        }else if(var_columns.info_map1[tmp_format_split[j]] == FLOAT_FORMAT_ALT){
                            //Float alternatives
                            boost::split(tmp_sub, tmp_split[j], boost::is_any_of(","));
                            local_alt = tmp_sub.size();
                            int el = 0;
                            while(!find_elem){ 
                                //Search the corresponding element
                                if(!(*tmp_alt_format).samp_float[el].name.compare(tmp_format_split[j])){
                                    //Fill a tuple for each alternatives
                                    for(int y = 0; y<local_alt; y++){
                                        (*tmp_alt_format).var_id[(*tmp_num_alt_format) + y] = var_columns.var_number[i];
                                        (*tmp_alt_format).samp_id[(*tmp_num_alt_format) + y] = samp;
                                        (*tmp_alt_format).alt_id[(*tmp_num_alt_format) + y] = (char)y;
                                        try{
                                            (*tmp_alt_format).samp_float[el].i_float[(*tmp_num_alt_format) + y] = (__half)std::stof(tmp_sub[y]);
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


#endif