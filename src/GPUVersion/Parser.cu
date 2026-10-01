/**
 * @file Parser.cu
 * @brief CUDA-accelerated VCF file parser implementation
 * @date 2025-07-16
 *
 * @details Implements a parallel VCF parser using CUDA:
 *  - Header parsing and metadata extraction
 *  - Memory management (host and device)
 *  - CUDA kernel launches and synchronization
 *  - Multi-threaded data merging
 *  - Sample and variant data processing
 *
 * @note Requires CUDA toolkit 11.0 or later
 * @warning Memory allocation sizes must account for maximum VCF file size
 */

#ifndef PARSER_CU
#define PARSER_CU

#include "DataStructures.h"
#include "Kernels.h"
#include "Utils.h"
#include "DataFrames.h"
#include "Parser.h"

#include <cuda_runtime.h>
#include <cuda_fp16.h>  
#include <fcntl.h>
#include <unistd.h>
#include <atomic>
#include <cerrno>
#include <charconv>
#include <unordered_map>
#include <unordered_set>

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
#include <string_view>
#include <stdexcept>
#include <exception>
#include <algorithm>
#include <cctype>


using namespace std;

// Throws instead of exit(): exit() would kill the Python interpreter running GPUParser.
// ponytail: device buffers allocated before the error are not freed; process exit reclaims them.
#define CUDA_CHECK_ERROR(call)                             \
    do {                                                   \
        cudaError_t err = call;                            \
        if (err != cudaSuccess) {                          \
            throw std::runtime_error(std::string("CUDA error in ") + #call \
                      + " at " + __FILE__ + ":" + std::to_string(__LINE__) \
                      + " - " + cudaGetErrorString(err));  \
        }                                                  \
    } while (0)



// Header attribute helpers.
/**
 * @brief True if a header Number= value is a fixed count ("0", "1", "2", ...).
 */
static inline bool is_fixed_number(const std::string& number) {
    return !number.empty() && std::all_of(number.begin(), number.end(), [](unsigned char c){ return std::isdigit(c); });
}

/**
 * @brief Returns the text between '<' and the last '>' of a ##INFO/##FORMAT header line.
 */
static inline std::string_view header_attr_body(std::string_view line) {
    const auto l = line.find('<');
    const auto r = line.rfind('>');
    if (l == std::string_view::npos || r == std::string_view::npos || r <= l) return {};
    return line.substr(l + 1, r - l - 1);
}

/**
 * @brief Returns the value of attribute keyEq (e.g. "ID=") up to the next comma, whatever its position.
 */
static inline std::string header_attr(std::string_view body, std::string_view keyEq) {
    const auto pos0 = body.find(keyEq);
    if (pos0 == std::string_view::npos) return {};
    const auto pos = pos0 + keyEq.size();
    auto end = body.find(',', pos);
    if (end == std::string_view::npos) end = body.size();
    return std::string(body.substr(pos, end - pos));
}

// VCF field walking: fields end at a tab, a space or the end of the line
static inline bool is_field_end(char c){ return c == '\t' || c == ' ' || c == '\n'; }

// End of the field starting at p (a tab, a space or the end of the line)
static inline const char* field_end(const char* p, const char* e){
    while(p < e && !is_field_end(*p)) ++p;
    return p;
}

// Next field after a field ending at p
static inline const char* next_field(const char* p, const char* e){
    return (p < e && (*p == '\t' || *p == ' ')) ? p + 1 : p;
}

// FORMAT columns are named <ID> (Number=1) or <ID>0, <ID>1, ... (Number>1): match the first column
// of key exactly. A prefix match picked GQX for GQ when GQX was declared first.
static inline bool format_name_matches(const std::string& name, const std::string& key){
    return name == key || (name.size() == key.size() + 1 && name.back() == '0' && name.compare(0, key.size(), key) == 0);
}

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
    int cudaCores;                 // Total number of CUDA cores

    // Query device properties
    int deviceCount = 0;
    cudaError_t error = cudaGetDeviceCount(&deviceCount);
    if (error != cudaSuccess || deviceCount == 0) {
        throw std::runtime_error("no CUDA-capable device found");
    }

    int deviceID = 0; 
    CUDA_CHECK_ERROR(cudaSetDevice(deviceID));

    cudaDeviceProp prop;
    cudaError_t err = cudaGetDeviceProperties(&prop, 0); // Query the first (and only) device

    if (err == cudaSuccess) {
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

    } else {
        std::cerr << "Failed to query device properties: " << cudaGetErrorString(err) << std::endl;
    }

    omp_set_num_threads(num_threadss);

    // Open input file, gzip -df compressed_file1.gz
    unzip_gz_file(vcf_filename); // no-op unless the name ends in .gz
    filename = vcf_filename; // after unzip_gz_file, which strips the .gz
    
    ifstream inFile(filename);
    if(!inFile){
        throw std::runtime_error("cannot open file " + filename);
    }
    // Saving filename
    filename = get_filename(filename, path_to_filename);
    
    // Getting filesize (number of char in the file)
    filesize = get_file_size(path_to_filename);
    // Getting the header (Saving the header into a string and storing the header size )
    get_and_parse_header(&inFile);
    inFile.close();
    // Allocating the filestring (the variations as a big char*, the dimension is: filesize - header_size)
    allocate_filestring();
    // Populate filestring and getting the number of lines (num_lines), saving the starting char index of each lines
    find_new_lines_index(path_to_filename, num_threadss);
    create_info_vectors(num_threadss);
    reserve_var_columns();
    create_sample_vectors(num_threadss);
    if(num_lines == 0) return; // header-only file: the columns exist and are empty
    // Allocate and initialize device memory
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

        std::strncpy(h_keys[index], key.c_str(), MAX_KEY_LENGTH_GT - 1);

        h_values[index] = value;

        ++index;
    }
    CUDA_CHECK_ERROR(upload_gt_table(h_keys, h_values));
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
    CUDA_CHECK_ERROR(upload_map1_table(h_keys, h_values));
}

/**
    * @brief Copies a field name into a fixed-size MAX_NAME_SIZE device slot (truncated, zero-padded).
    */
static void copy_fixed_name(char* dst, const string& name){
    char buf[MAX_NAME_SIZE] = {0};
    strncpy(buf, name.c_str(), MAX_NAME_SIZE - 1);
    CUDA_CHECK_ERROR(cudaMemcpy(dst, buf, MAX_NAME_SIZE, cudaMemcpyHostToDevice));
}

/**
    * @brief Allocates device memory for VCF parsing data.
    *
    * Allocates memory on the GPU for variant numbers, positions, quality scores,
    * and INFO fields. If sample data is present, it also allocates memory for sample fields.
    */
void vcf_parsed::device_allocation(){
    // The kernel writes a column cell only when the record has that field: every column is zeroed here,
    // or the cells left unwritten would be copied back with whatever the reused device memory held.
    CUDA_CHECK_ERROR(cudaMalloc(&d_VC_var_number, (num_lines) * sizeof(unsigned int)));
    CUDA_CHECK_ERROR(cudaMalloc(&d_VC_pos, (num_lines) * sizeof(unsigned int)));
    CUDA_CHECK_ERROR(cudaMalloc(&d_VC_qual, (num_lines) * sizeof(__half)));

    int tmp = var_columns.in_float.size();

    CUDA_CHECK_ERROR(cudaMalloc(&(d_VC_in_float->i_float), tmp * (num_lines) * sizeof(__half)));
    CUDA_CHECK_ERROR(cudaMemset(d_VC_in_float->i_float, 0, tmp * (num_lines) * sizeof(__half)));
    CUDA_CHECK_ERROR(cudaMalloc(&(d_VC_in_float->name), tmp * sizeof(char) * MAX_NAME_SIZE)); 

    for (int i = 0; i < tmp; i++) {
        copy_fixed_name(d_VC_in_float->name + i*MAX_NAME_SIZE, var_columns.in_float[i].name);
    }

    tmp = var_columns.in_flag.size();
    CUDA_CHECK_ERROR(cudaMalloc(&(d_VC_in_flag->i_flag), tmp * (num_lines) * sizeof(uint8_t)));
    CUDA_CHECK_ERROR(cudaMemset(d_VC_in_flag->i_flag, 0, tmp * (num_lines) * sizeof(uint8_t)));
    CUDA_CHECK_ERROR(cudaMalloc(&(d_VC_in_flag->name), tmp * sizeof(char) * MAX_NAME_SIZE));

    for (int i = 0; i < tmp; i++) {
        copy_fixed_name(d_VC_in_flag->name + i*MAX_NAME_SIZE, var_columns.in_flag[i].name);
    }

    tmp = var_columns.in_int.size();
    CUDA_CHECK_ERROR(cudaMalloc(&(d_VC_in_int->i_int), tmp * (num_lines) * sizeof(int)));
    CUDA_CHECK_ERROR(cudaMemset(d_VC_in_int->i_int, 0, tmp * (num_lines) * sizeof(int)));
    CUDA_CHECK_ERROR(cudaMalloc(&(d_VC_in_int->name), tmp * sizeof(char) * MAX_NAME_SIZE));
    for (int i = 0; i < tmp; i++) {
        copy_fixed_name(d_VC_in_int->name + i*MAX_NAME_SIZE, var_columns.in_int[i].name);
    }

    initialize_map1(var_columns.info_map1);
    
    if (hasDetSamples) {
        // Copy map to constant memory
        copyMapToConstantMemory(samp_columns.GTMap);
        // Allocate and initialize d_SC_var_id
        CUDA_CHECK_ERROR(cudaMalloc(&d_SC_var_id, (num_lines) * samp_columns.numSample * sizeof(unsigned int)));
        CUDA_CHECK_ERROR(cudaMemset(d_SC_var_id, 0, (num_lines) * samp_columns.numSample * sizeof(unsigned int)));
        
        // Allocate and initialize d_SC_samp_id
        CUDA_CHECK_ERROR(cudaMalloc(&d_SC_samp_id, (num_lines) * samp_columns.numSample * sizeof(unsigned short)));
        CUDA_CHECK_ERROR(cudaMemset(d_SC_samp_id, 0, (num_lines) * samp_columns.numSample * sizeof(unsigned short)));

        // Allocate and initialize samp_float
        tmp = samp_columns.samp_float.size();
        CUDA_CHECK_ERROR(cudaMalloc(&(d_SC_samp_float->i_float), tmp * (num_lines * samp_columns.numSample) * sizeof(__half)));
        CUDA_CHECK_ERROR(cudaMemset(d_SC_samp_float->i_float, 0, tmp * (num_lines * samp_columns.numSample) * sizeof(__half)));
        CUDA_CHECK_ERROR(cudaMalloc(&(d_SC_samp_float->name), tmp * sizeof(char) * MAX_NAME_SIZE));
        CUDA_CHECK_ERROR(cudaMalloc(&(d_SC_samp_float->numb), tmp * sizeof(int)));

        for (int i = 0; i < tmp; i++) {
            copy_fixed_name(d_SC_samp_float->name + i*MAX_NAME_SIZE, samp_columns.samp_float[i].name);
            CUDA_CHECK_ERROR(cudaMemcpy(d_SC_samp_float->numb+i, &(samp_columns.samp_float[i].numb), sizeof(int), cudaMemcpyHostToDevice));
        }

        // Allocate and initialize samp_flag
        tmp = samp_columns.samp_flag.size();
        CUDA_CHECK_ERROR(cudaMalloc(&(d_SC_samp_flag->i_flag), tmp * (num_lines * samp_columns.numSample) * sizeof(uint8_t)));
        CUDA_CHECK_ERROR(cudaMemset(d_SC_samp_flag->i_flag, 0, tmp * (num_lines * samp_columns.numSample) * sizeof(uint8_t)));
        CUDA_CHECK_ERROR(cudaMalloc(&(d_SC_samp_flag->name), tmp * sizeof(char) * MAX_NAME_SIZE));
        CUDA_CHECK_ERROR(cudaMalloc(&(d_SC_samp_flag->numb), tmp * sizeof(int)));

        for (int i = 0; i < tmp; i++) {
            copy_fixed_name(d_SC_samp_flag->name + i*MAX_NAME_SIZE, samp_columns.samp_flag[i].name);
            CUDA_CHECK_ERROR(cudaMemcpy(d_SC_samp_flag->numb+i, &(samp_columns.samp_flag[i].numb), sizeof(int), cudaMemcpyHostToDevice));
        }

        // Allocate and initialize samp_int
        tmp = samp_columns.samp_int.size();
        CUDA_CHECK_ERROR(cudaMalloc(&(d_SC_samp_int->i_int), tmp * (num_lines * samp_columns.numSample) * sizeof(int)));
        CUDA_CHECK_ERROR(cudaMemset(d_SC_samp_int->i_int, 0, tmp * (num_lines * samp_columns.numSample) * sizeof(int)));
        CUDA_CHECK_ERROR(cudaMalloc(&(d_SC_samp_int->name), tmp * sizeof(char) * MAX_NAME_SIZE));
        CUDA_CHECK_ERROR(cudaMalloc(&(d_SC_samp_int->numb), tmp * sizeof(int)));

        for (int i = 0; i < tmp; i++) {
            copy_fixed_name(d_SC_samp_int->name + i*MAX_NAME_SIZE, samp_columns.samp_int[i].name);
            CUDA_CHECK_ERROR(cudaMemcpy(d_SC_samp_int->numb+i, &(samp_columns.samp_int[i].numb), sizeof(int), cudaMemcpyHostToDevice));
        }

        // Allocate and initialize samp_GT
        tmp = samp_columns.sample_GT.size();
        CUDA_CHECK_ERROR(cudaMalloc(&(d_SC_sample_GT->GT), tmp * (num_lines * samp_columns.numSample) * sizeof(char)));
        CUDA_CHECK_ERROR(cudaMemset(d_SC_sample_GT->GT, 0, tmp * (num_lines * samp_columns.numSample) * sizeof(char)));
        CUDA_CHECK_ERROR(cudaMalloc(&(d_SC_sample_GT->numb), sizeof(int)));
        // sample_GT is empty when GT is not declared (or is Number=A, parsed on the host)
        if(!samp_columns.sample_GT.empty())
            CUDA_CHECK_ERROR(cudaMemcpy(d_SC_sample_GT->numb, &(samp_columns.sample_GT[0].numb), sizeof(int), cudaMemcpyHostToDevice));

    }

}

/**
    * @brief Frees all allocated device memory.
    *
    * Releases device memory allocated during the parsing process.
    */
void vcf_parsed::device_free() {
    CUDA_CHECK_ERROR(cudaFree(d_VC_var_number));
    CUDA_CHECK_ERROR(cudaFree(d_VC_pos));
    CUDA_CHECK_ERROR(cudaFree(d_VC_qual));
    CUDA_CHECK_ERROR(cudaFree(d_VC_in_float->i_float));
    CUDA_CHECK_ERROR(cudaFree(d_VC_in_float->name));
    CUDA_CHECK_ERROR(cudaFree(d_VC_in_flag->i_flag));
    CUDA_CHECK_ERROR(cudaFree(d_VC_in_flag->name));
    CUDA_CHECK_ERROR(cudaFree(d_VC_in_int->i_int));
    CUDA_CHECK_ERROR(cudaFree(d_VC_in_int->name));
    CUDA_CHECK_ERROR(cudaFree(d_filestring));
    CUDA_CHECK_ERROR(cudaFree(d_new_lines_index));

    if (hasDetSamples) {
        CUDA_CHECK_ERROR(cudaFree(d_SC_var_id));
        CUDA_CHECK_ERROR(cudaFree(d_SC_samp_id));
        CUDA_CHECK_ERROR(cudaFree(d_SC_samp_float->i_float));
        CUDA_CHECK_ERROR(cudaFree(d_SC_samp_float->name));
        CUDA_CHECK_ERROR(cudaFree(d_SC_samp_float->numb));
        CUDA_CHECK_ERROR(cudaFree(d_SC_samp_flag->i_flag));
        CUDA_CHECK_ERROR(cudaFree(d_SC_samp_flag->name));
        CUDA_CHECK_ERROR(cudaFree(d_SC_samp_flag->numb));
        CUDA_CHECK_ERROR(cudaFree(d_SC_samp_int->i_int));
        CUDA_CHECK_ERROR(cudaFree(d_SC_samp_int->name));
        CUDA_CHECK_ERROR(cudaFree(d_SC_samp_int->numb));
        CUDA_CHECK_ERROR(cudaFree(d_SC_sample_GT->GT));
        CUDA_CHECK_ERROR(cudaFree(d_SC_sample_GT->numb));
    }
    CUDA_CHECK_ERROR(cudaFree(d_params));
}

/**
    * @brief Reads the variant body and indexes its records.
    *
    * Reads the body into filestring with parallel pread() calls, then builds new_lines_index on the host
    * (one entry per record end, blank lines skipped) and copies both to the device for the kernel.
    *
    * @param w_filename The path to the VCF file.
    * @param num_threads Number of threads to use for parallel processing.
    */
void vcf_parsed::find_new_lines_index(string w_filename, int num_threads){
    // Parallel pread of the variant body straight into filestring, then each chunk counts the record ends
    // and writes their positions at its prefix offset. Chunks are in file order and scanned forward,
    // so new_lines_index comes out sorted: no device kernel, sort or per-chunk buffers are needed.
    const long filestring_size = variants_size + 8; // as allocated by allocate_filestring
    const long batch_infile = (variants_size - 1 + num_threads)/num_threads; // Number of characters each chunk holds

    const int fd = open(w_filename.c_str(), O_RDONLY);
    if(fd < 0) throw std::runtime_error("cannot open file " + w_filename);
    posix_fadvise(fd, header_size, variants_size, POSIX_FADV_SEQUENTIAL); // read-ahead hint, best effort
    std::atomic<int> read_errno{0};

    #pragma omp parallel for schedule(static)
    for(int c = 0; c < num_threads; c++){
        const long start = std::min(c*batch_infile, variants_size);
        const size_t want = std::min(start + batch_infile, variants_size) - start;
        size_t got = 0;
        while(got < want){
            const ssize_t n = pread(fd, filestring + start + got, want - got, (off_t)header_size + start + got);
            if(n > 0){ got += n; continue; }
            if(n < 0 && errno == EINTR) continue;
            int none = 0;
            read_errno.compare_exchange_strong(none, n == 0 ? EIO : errno); // n == 0: the file got shorter
            break;
        }
    }
    close(fd);
    if(read_errno) throw std::runtime_error("cannot read " + w_filename + ": " + strerror(read_errno));

    // Trailing newlines are dropped and a single '\n' terminator closes the last record
    while(variants_size > 0 && filestring[variants_size-1]=='\n') variants_size--;
    if(variants_size == 0){ num_lines = 0; return; } // header-only file: no records, nothing for the device
    filestring[variants_size] = '\n';
    variants_size++;
    memset(filestring + variants_size, '\0', filestring_size - variants_size); // NUL tail after the terminator

    // A '\n' ends a record unless a blank line follows it (the CPU parser skips blank lines the same way);
    // the terminator, followed by the NUL tail, always does
    const long batch = (variants_size + num_threads - 1)/num_threads;
    auto chunk_start = [&](int c){ return std::min(c*batch, variants_size); };
    // Calls f(position) for every record end of chunk c
    auto for_each_record_end = [&](int c, auto&& f){
        const char* p = filestring + chunk_start(c);
        const char* const e = filestring + chunk_start(c + 1);
        while((p = static_cast<const char*>(memchr(p, '\n', e - p)))){
            if(p[1] != '\n') f(static_cast<unsigned long long>(p - filestring));
            ++p;
        }
    };
    std::vector<size_t> chunk_count(num_threads + 1, 0);
    #pragma omp parallel for schedule(static)
    for(int c = 0; c < num_threads; c++) for_each_record_end(c, [&](unsigned long long){ chunk_count[c + 1]++; });
    for(int c = 0; c < num_threads; c++) chunk_count[c + 1] += chunk_count[c]; // prefix offsets

    // new_lines_index = [0, every record end..., terminator]: record i spans new_lines_index[i]..
    // new_lines_index[i+1] (a '\n' then the record, or the record itself for i == 0)
    num_lines = chunk_count[num_threads];
    new_lines_index = (unsigned long long*)malloc(sizeof(unsigned long long)*(num_lines + 1));
    if(!new_lines_index) throw std::runtime_error("cannot allocate the line index (" + std::to_string(num_lines + 1) + " entries)");
    new_lines_index[0] = 0;
    #pragma omp parallel for schedule(static)
    for(int c = 0; c < num_threads; c++){
        unsigned long long* out = new_lines_index + 1 + chunk_count[c];
        for_each_record_end(c, [&](unsigned long long pos){ *out++ = pos; });
    }
    if(filestring[0] == '\n'){ // blank lines before the first record: its start is the first record end found
        memmove(new_lines_index, new_lines_index + 1, sizeof(unsigned long long)*num_lines);
        num_lines--;
    }
    CUDA_CHECK_ERROR(cudaMalloc(&d_filestring, (variants_size + 8)* sizeof(char)));
    CUDA_CHECK_ERROR(cudaMalloc(&d_new_lines_index, (num_lines + 1) * sizeof(unsigned long long)));
    CUDA_CHECK_ERROR(cudaMemcpy(d_filestring, filestring, sizeof(char)*variants_size, cudaMemcpyHostToDevice));
    CUDA_CHECK_ERROR(cudaMemcpy(d_new_lines_index, new_lines_index, sizeof(unsigned long long)*(num_lines+1), cudaMemcpyHostToDevice));
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
    variants_size = std::max(0L, filesize - header_size); // New size without the header (0 when the header has no final '\n')
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
    // removing the header and storing it in vcf.header
    
    while (getline(*file, line) && line[0]=='#' && line[1]=='#'){
        header.append(line + '\n');
        header_size += line.length() + 1;
        bool Info = (line.rfind("##INFO=<", 0) == 0);
        bool Format = (line.rfind("##FORMAT=<", 0) == 0);
        
        if(Info || Format){
            // Attributes are looked up by name: the VCF spec does not fix their order.
            const std::string_view body = header_attr_body(line);
            const string id = header_attr(body, "ID=");
            const string number = header_attr(body, "Number=");
            const string type = header_attr(body, "Type=");
            if(Info){
                INFO.ID.push_back(id);
                INFO.Number.push_back(number);
                if(number == "A") INFO.alt_values++;
                INFO.Type.push_back(type);
            }else if(id == "GT"){
                FORMAT.hasGT = true;
                // GT is Number=1 by spec: anything but A or a fixed count (e.g. '.') counts as 1
                FORMAT.numGT = (number == "A" || is_fixed_number(number)) ? number[0] : '1';
                hasDetSamples = true;
            }else{
                FORMAT.ID.push_back(id);
                FORMAT.Number.push_back(number);
                if(number == "A"){
                    FORMAT.alt_values++;
                }else{
                    hasDetSamples = true;
                }
                FORMAT.Type.push_back(type);
            }
        }
    }

    vector<string> tmp_split;
    split_on(tmp_split, line, '\t'); // tab only: sample names may contain spaces
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

    variants_size = std::max(0L, filesize - header_size); // New size without the header (0 when the header has no final '\n')
}   
    
/**
    * @brief Allocates a character array to store the variant portion of the VCF file.
    *
    * The allocated size is based on the file size minus the header size.
    */
void vcf_parsed::allocate_filestring(){
    // No memset: find_new_lines_index overwrites the whole body and zeroes the tail
    filestring = (char*)malloc(variants_size + 8);
    if(!filestring) throw std::runtime_error("cannot allocate " + std::to_string(variants_size + 8) + " bytes for the VCF body");
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
    if(numIter == 0 && !FORMAT.hasGT) return;

    for(int i = 0; i < numIter; i++){
        // Number=R, G and . are not supported yet: skip the field instead of throwing in
        // std::stoi (R, G) or prompting on stdin (.). The line parser then ignores it.
        if(strcmp(&FORMAT.Number[i][0], "A") != 0 && !is_fixed_number(FORMAT.Number[i])) continue;
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
                info_map[FORMAT.ID[i]] = FLAG_FORMAT;
                var_columns.info_map1[FORMAT.ID[i]] = FLAG_FORMAT;
                FORMAT.flags++;
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
    if(hasDetSamples){ // only Number=A FORMAT fields: everything goes to DF4, DF3 stays empty
        samp_columns.var_id.resize((num_lines)*samp_columns.numSample, 0);
        samp_columns.samp_id.resize((num_lines)*samp_columns.numSample, static_cast<unsigned short>(0));
    }    
    alt_sample.samp_flag.resize(FORMAT.flags_alt);
    alt_sample.samp_int.resize(FORMAT.ints_alt);
    alt_sample.samp_float.resize(FORMAT.floats_alt);
    alt_sample.samp_string.resize(FORMAT.strings_alt);

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
            if(strcmp(&INFO.Type[i][0], "Flag")==0){ // not handled yet
                INFO.flags_alt++;
                info_map[INFO.ID[i]] = 7;
            }
        }else if((strcmp(&INFO.Number[i][0], "1") == 0)||(strcmp(&INFO.Number[i][0], "0") == 0)){
            // Without Alternatives and number = 1 or a flag
            if(strcmp(&INFO.Type[i][0], "Integer")==0){
                INFO.ints++;
                info_int_tmp.name = INFO.ID[i];
                info_int_tmp.i_int.resize(num_lines, 0);
                var_columns.in_int.push_back(info_int_tmp);
                info_map[INFO.ID[i]] = 1;
                var_columns.info_map1[INFO.ID[i]] = 1;
            } else if(strcmp(&INFO.Type[i][0], "Float")==0){
                INFO.floats++;
                info_float_tmp.name = INFO.ID[i];
                info_float_tmp.i_float.resize(num_lines, 0);
                var_columns.in_float.push_back(info_float_tmp);
                info_map[INFO.ID[i]] = 2;
                var_columns.info_map1[INFO.ID[i]] = 2;
            } else if(strcmp(&INFO.Type[i][0], "String")==0){
                if(strcmp(INFO.ID[i].c_str(), "TSA")==0){
                    INFO.ints++;
                    info_int_tmp.name = INFO.ID[i];
                    info_int_tmp.i_int.resize(num_lines, 0);
                    var_columns.in_int.push_back(info_int_tmp);
                    info_map[INFO.ID[i]] = 1;
                    var_columns.info_map1[INFO.ID[i]] = 1;
                }else{ 
                    INFO.strings++;
                    info_string_tmp.name = INFO.ID[i];
                    info_string_tmp.i_string.resize(num_lines, "\0");
                    var_columns.in_string.push_back(info_string_tmp);
                    info_map[INFO.ID[i]] = 3;
                    var_columns.info_map1[INFO.ID[i]] = 3;
                }
            } else if(strcmp(&INFO.Type[i][0], "Flag")==0){
                INFO.flags++;
                info_flag_tmp.name = INFO.ID[i];
                info_flag_tmp.i_flag.resize(num_lines, 0);
                var_columns.in_flag.push_back(info_flag_tmp);
                info_map[INFO.ID[i]] = 0;
                var_columns.info_map1[INFO.ID[i]] = 0;
            }
        }else {
            /**
             * @todo Implement support for INFO fields with Number != 1
             * 
             * @details Valid Number values to handle:
             *  - 0: Flag type (already handled)
             *  - 1: Single value (already handled)
             *  - R: One value per allele (including reference)
             *  - A: One value per alternate allele
             *  - G: One value per possible genotype
             *  - .: Unknown number of values
             */
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
    var_columns.var_number.resize(num_lines);
    var_columns.chrom.resize(num_lines);
    var_columns.id.resize(num_lines);
    var_columns.pos.resize(num_lines);
    var_columns.ref.resize(num_lines); 
    var_columns.qual.resize(num_lines);
    var_columns.filter.resize(num_lines);
}

/**
    * @brief Allocates device memory for KernelParams and copies the host structure.
    *
    * @param d_params Double pointer to the device KernelParams.
    * @param h_params Pointer to the host KernelParams structure.
    */
void vcf_parsed::allocParamPointers(KernelParams **d_params, KernelParams *h_params) {
    // Allocate memory for KernelParams on GPU
    CUDA_CHECK_ERROR(cudaMalloc((void**)d_params, sizeof(KernelParams)));

    // Copy the structure from host to device
    CUDA_CHECK_ERROR(cudaMemcpy(*d_params, h_params, sizeof(KernelParams), cudaMemcpyHostToDevice));
}

/**
    * @brief Launches the CUDA kernel to parse VCF lines and merges the results.
    *
    * Sets up CUDA streams and events, launches the appropriate kernel (based on whether sample data is present),
    * and asynchronously copies the parsed data from device to host.
    */
void vcf_parsed::populate_runner(int numb_cores){
    int threadsPerBlock = 32;
    int blocksPerGrid = (numb_cores/threadsPerBlock) + 1; 
    cudaEvent_t kernel_done;
    CUDA_CHECK_ERROR(cudaEventCreate(&kernel_done));
    auto start = chrono::system_clock::now();
    char* my_mem;
    int batchSize = threadsPerBlock*blocksPerGrid;
    CUDA_CHECK_ERROR(cudaMalloc(&my_mem, batchSize*MAX_TOKEN_LEN*MAX_TOKENS*3));

    cudaStream_t stream1, stream2;
    CUDA_CHECK_ERROR(cudaStreamCreate(&stream1));
    CUDA_CHECK_ERROR(cudaStreamCreate(&stream2));

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
    h_params.numInfoFloat = var_columns.in_float.size();
    h_params.numInfoFlag = var_columns.in_flag.size();
    h_params.numInfoInt = var_columns.in_int.size();
    h_params.new_lines_index = d_new_lines_index;
    h_params.numLines = num_lines;

    // The sample buffers are allocated only when hasDetSamples (GT or a non Number=A FORMAT field)
    if(hasDetSamples){            
        
        h_params.samp_var_id = d_SC_var_id;
        h_params.samp_id = d_SC_samp_id;
        h_params.samp_float = d_SC_samp_float->i_float;
        h_params.samp_flag = d_SC_samp_flag->i_flag;
        h_params.samp_int = d_SC_samp_int->i_int;
        h_params.samp_float_name = d_SC_samp_float->name;
        h_params.samp_flag_name = d_SC_samp_flag->name;
        h_params.samp_int_name = d_SC_samp_int->name;
        h_params.numSampFloat = samp_columns.samp_float.size();
        h_params.numSampInt = samp_columns.samp_int.size();
        h_params.samp_float_numb = d_SC_samp_float->numb;
        h_params.samp_flag_numb = d_SC_samp_flag->numb;
        h_params.samp_int_numb = d_SC_samp_int->numb;
        h_params.sample_GT = d_SC_sample_GT->GT;
        h_params.numSample = samp_columns.numSample;
        h_params.numGT = (int)samp_columns.sample_GT.size(); // 0 without a Number=1 GT column
        h_params.hasGT = FORMAT.hasGT;

        // Allocate d_params and copy h_params to GPU
        allocParamPointers(&d_params, &h_params);

        // Launch kernel
        CUDA_CHECK_ERROR(launch_parse_kernel(blocksPerGrid, threadsPerBlock, stream1, d_params, my_mem, batchSize, true));

        CUDA_CHECK_ERROR(cudaEventRecord(kernel_done, stream1));

    }else{
        // Allocate d_params and copy h_params to GPU
        allocParamPointers(&d_params, &h_params);
        CUDA_CHECK_ERROR(launch_parse_kernel(blocksPerGrid, threadsPerBlock, stream1, d_params, my_mem, batchSize, false));
        CUDA_CHECK_ERROR(cudaEventRecord(kernel_done, stream1));
        
    }
    
    CUDA_CHECK_ERROR(cudaStreamWaitEvent(stream2, kernel_done, 0));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(var_columns.var_number.data(), d_VC_var_number, (num_lines) * sizeof(unsigned int), cudaMemcpyDeviceToHost, stream2));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(var_columns.pos.data(), d_VC_pos, (num_lines) * sizeof(unsigned int), cudaMemcpyDeviceToHost, stream2));
    CUDA_CHECK_ERROR(cudaMemcpyAsync(var_columns.qual.data(), d_VC_qual, (num_lines) * sizeof(__half), cudaMemcpyDeviceToHost, stream2));

    for(int i=0; i<var_columns.in_float.size(); i++){
        CUDA_CHECK_ERROR(cudaMemcpyAsync(var_columns.in_float[i].i_float.data(), d_VC_in_float->i_float + i * (num_lines), (num_lines)*sizeof(__half), cudaMemcpyDeviceToHost, stream2));
    }

    for(int i=0; i<var_columns.in_flag.size(); i++){
        CUDA_CHECK_ERROR(cudaMemcpyAsync(var_columns.in_flag[i].i_flag.data(), d_VC_in_flag->i_flag + i * (num_lines), (num_lines)*sizeof(uint8_t), cudaMemcpyDeviceToHost, stream2));
    }

    for(int i=0; i<var_columns.in_int.size(); i++){
        CUDA_CHECK_ERROR(cudaMemcpyAsync(var_columns.in_int[i].i_int.data(), d_VC_in_int->i_int + i * (num_lines), (num_lines) * sizeof(int), cudaMemcpyDeviceToHost, stream2));
    }

    if(hasDetSamples){
        CUDA_CHECK_ERROR(cudaMemcpyAsync(samp_columns.var_id.data(), d_SC_var_id, (num_lines) * samp_columns.numSample * sizeof(unsigned int), cudaMemcpyDeviceToHost, stream2));
        CUDA_CHECK_ERROR(cudaMemcpyAsync(samp_columns.samp_id.data(), d_SC_samp_id, (num_lines) * samp_columns.numSample * sizeof(unsigned short), cudaMemcpyDeviceToHost, stream2));

        for (int i = 0; i < samp_columns.samp_float.size(); i++) {
            CUDA_CHECK_ERROR(cudaMemcpyAsync(samp_columns.samp_float[i].i_float.data(), d_SC_samp_float->i_float + i * ((num_lines) * samp_columns.numSample), 
                (num_lines) * samp_columns.numSample * sizeof(__half), cudaMemcpyDeviceToHost, stream2));
        }

        for (int i = 0; i < samp_columns.samp_flag.size(); i++) {
            CUDA_CHECK_ERROR(cudaMemcpyAsync(samp_columns.samp_flag[i].i_flag.data(), d_SC_samp_flag->i_flag + i * ((num_lines) * samp_columns.numSample), 
                (num_lines) * samp_columns.numSample * sizeof(uint8_t), cudaMemcpyDeviceToHost, stream2));
        }

        for (int i = 0; i < samp_columns.samp_int.size(); i++) {
            CUDA_CHECK_ERROR(cudaMemcpyAsync(samp_columns.samp_int[i].i_int.data(), d_SC_samp_int->i_int + (i * num_lines * samp_columns.numSample), 
                (num_lines) * samp_columns.numSample * sizeof(int), cudaMemcpyDeviceToHost, stream2));                
        }   

        for(int i=0; i<samp_columns.sample_GT.size(); i++){
            CUDA_CHECK_ERROR(cudaMemcpyAsync(samp_columns.sample_GT[i].GT.data(), d_SC_sample_GT->GT + i * ((num_lines) * samp_columns.numSample), 
                (num_lines)*samp_columns.numSample*sizeof(char), cudaMemcpyDeviceToHost, stream2));
        } 
    }
    CUDA_CHECK_ERROR(cudaStreamSynchronize(stream2));

    // Cleanup
    CUDA_CHECK_ERROR(cudaEventDestroy(kernel_done));
    CUDA_CHECK_ERROR(cudaStreamDestroy(stream1));
    CUDA_CHECK_ERROR(cudaStreamDestroy(stream2));
    CUDA_CHECK_ERROR(cudaFree(my_mem));
}

/**
    * @brief Fills var_columns.chrom_map / filter_map single-threaded, before the parallel parse.
    *
    * The parsing threads used to insert into these std::map concurrently (a data race that
    * corrupted the trees). Codes are assigned in order of first appearance in the file.
    */
void vcf_parsed::prebuild_chrom_filter_maps(){
    var_columns.chrom_map.clear();
    var_columns.filter_map.clear();

    // Each chunk lists its distinct CHROM / FILTER names in order of first appearance; merging the
    // chunks in file order then gives the same codes as a single sequential pass.
    const int n_chunks = std::max(1, omp_get_max_threads());
    const long lines_per_chunk = (num_lines + n_chunks - 1)/n_chunks;
    std::vector<std::vector<std::string_view>> chunk_chroms(n_chunks), chunk_filters(n_chunks);
    std::exception_ptr chunk_error; // an exception may not leave an OpenMP region: kept, rethrown after it

    #pragma omp parallel for schedule(static)
    for(int c = 0; c < n_chunks; c++) try {
        std::unordered_set<std::string_view> seen_chrom, seen_filter;
        const long first = c*lines_per_chunk, last = std::min(num_lines, first + lines_per_chunk);
        for(long i = first; i < last; i++){
            const char* p = filestring + new_lines_index[i];
            const char* e = filestring + new_lines_index[i + 1];
            if(p < e && *p == '\n') p++;

            const char* q = field_end(p, e);
            const std::string_view chrom(p, q - p);
            if(seen_chrom.insert(chrom).second) chunk_chroms[c].push_back(chrom);
            p = next_field(q, e);

            for(int field = 2; field <= 6 && p < e; field++) p = next_field(field_end(p, e), e); // skip POS, ID, REF, ALT, QUAL

            q = field_end(p, e);
            const std::string_view filter(p, q - p);
            if(seen_filter.insert(filter).second) chunk_filters[c].push_back(filter);
        }
    } catch(...) {
        #pragma omp critical(prebuild_chrom_filter_error)
        if(!chunk_error) chunk_error = std::current_exception();
    }
    if(chunk_error) std::rethrow_exception(chunk_error);
    for(int c = 0; c < n_chunks; c++){
        for(auto name : chunk_chroms[c]) var_columns.chrom_map.emplace(std::string(name), static_cast<unsigned char>(var_columns.chrom_map.size()));
        for(auto name : chunk_filters[c]) var_columns.filter_map.emplace(std::string(name), static_cast<char>(var_columns.filter_map.size()));
    }
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
    prebuild_chrom_filter_maps();
    build_host_lookup();

    // The CUDA worker runs next to the host parse: an exception thrown there would call
    // std::terminate, so it is caught and rethrown on this thread after join().
    std::exception_ptr worker_error;
    std::thread worker_thread([this, numb_cores, &worker_error]{
        try{
            populate_runner(numb_cores);
        }catch(...){
            cudaDeviceSynchronize(); // let queued async copies finish before the host vectors can go away
            worker_error = std::current_exception();
        }
    });

    long batch_size = (num_lines-1+num_threads)/num_threads;
    
    std::vector<alt_columns_df> tmp_alt(num_threads);
    std::vector<int> tmp_num_alt(num_threads);

    std::vector<alt_format_df> tmp_alt_format(num_threads);
    std::vector<int> tmp_num_alt_format(num_threads);

    
    int totAlt = 0;
    int totSampAlt = 0;

    // One chunk per index (not per OpenMP thread): every chunk is parsed even with fewer threads
    #pragma omp parallel for schedule(static)
    for(int th_ID = 0; th_ID < num_threads; th_ID++)
    {
        long start, end;
        format_plan_cache plans; // FORMAT templates classified in this chunk
        // Temporary structure of the thread with alternatives.
        tmp_alt[th_ID].init(alt_columns, INFO, batch_size);
        
        tmp_num_alt[th_ID] = 0;

        start = th_ID * batch_size; // Starting point of the thread's batch
        end = start + batch_size; // Ending point of the thread's batch
        
        if(samplesON){
            // There are samples in the dataset
            tmp_alt_format[th_ID].init(alt_sample, FORMAT, batch_size);
            tmp_num_alt_format[th_ID] = 0;
            if(FORMAT.hasGT && FORMAT.numGT == 'A'){
                tmp_alt_format[th_ID].sample_GT.GT.resize(batch_size*2*samp_columns.numSample, (char)0),
                tmp_alt_format[th_ID].initMapGT();
            }

            // For each line in the batch
            for(long i=start; i<end && i<num_lines; i++){ 
                get_vcf_line_in_var_columns_format(filestring, new_lines_index[i], new_lines_index[i+1], i, &(tmp_alt[th_ID]), &(tmp_num_alt[th_ID]), &samp_columns, &FORMAT, &(tmp_num_alt_format[th_ID]), &(tmp_alt_format[th_ID]), &plans);
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
            for(long i=start; i<end && i<num_lines; i++){ 
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
    // Resize the alt_columns vectors in parallel
    {
        // Task resizing the flat vectors
        auto fut1 = std::async(std::launch::async, [&]() {
            alt_columns.var_id.resize(totAlt);
            alt_columns.alt_id.resize(totAlt);
            alt_columns.alt.resize(totAlt);
        });
        
        // Task resizing the inner vectors of alt_int
        auto fut2 = std::async(std::launch::async, [&]() {
            for (int j = 0; j < INFO.ints_alt; j++) {
                alt_columns.alt_int[j].i_int.resize(totAlt);
            }
        });
        
        // Task resizing the inner vectors of alt_float
        auto fut3 = std::async(std::launch::async, [&]() {
            for (int j = 0; j < INFO.floats_alt; j++) {
                alt_columns.alt_float[j].i_float.resize(totAlt);
            }
        });
        
        // Task resizing the inner vectors of alt_string
        auto fut4 = std::async(std::launch::async, [&]() {
            for (int j = 0; j < INFO.strings_alt; j++) {
                alt_columns.alt_string[j].i_string.resize(totAlt);
            }
        });
        
        // Wait for all the tasks to finish
        fut1.get();
        fut2.get();
        fut3.get();
        fut4.get();
    }

    // With samplesON, do the same for alt_sample
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
    if(worker_error) std::rethrow_exception(worker_error);
}

// Per-thread alternative buffers are pre-sized for ~2 ALTs per line; grow them when a chunk needs more
// (same approach as the CPU backend).
static void ensure_alt_capacity(alt_columns_df* tmp_alt, int needed){
    if(needed <= 0 || static_cast<int>(tmp_alt->alt.size()) >= needed) return;
    tmp_alt->var_id.resize(needed, 0);
    tmp_alt->alt.resize(needed, "\0");
    tmp_alt->alt_id.resize(needed, (char)0);
    for(auto& c : tmp_alt->alt_int) c.i_int.resize(needed, 0);
    for(auto& c : tmp_alt->alt_float) c.i_float.resize(needed, 0.0f);
    for(auto& c : tmp_alt->alt_string) c.i_string.resize(needed, "\0");
}

static void ensure_alt_format_capacity(alt_format_df* tmp_alt_format, int needed){
    if(needed <= 0 || static_cast<int>(tmp_alt_format->var_id.size()) >= needed) return;
    tmp_alt_format->var_id.resize(needed, 0);
    tmp_alt_format->alt_id.resize(needed, (char)0);
    tmp_alt_format->samp_id.resize(needed, static_cast<unsigned short>(0));
    for(auto& c : tmp_alt_format->samp_int) c.i_int.resize(needed, 0);
    for(auto& c : tmp_alt_format->samp_float) c.i_float.resize(needed, 0.0f);
    for(auto& c : tmp_alt_format->samp_string) c.i_string.resize(needed, "\0");
    if(!tmp_alt_format->sample_GT.GT.empty()) tmp_alt_format->sample_GT.GT.resize(needed, (char)0);
}

// ---- Host line parsing -------------------------------------------------------------------------
// The host side parses what the kernel does not: CHROM, ID, REF, ALT, FILTER, INFO String Number=1,
// INFO Number=A, FORMAT String, FORMAT Number=A and GT Number=A. Fields are walked with pointers
// (no std::string or split per field) and keys are resolved once, in build_host_lookup.

// Calls f(token_begin, token_end, index) for every sep-separated token of [b, e) (one token if none)
template <class F>
static inline int for_each_token(const char* b, const char* e, char sep, F&& f){
    int n = 0;
    for(const char* t = b;; ++n){
        const char* q = static_cast<const char*>(memchr(t, sep, e - t));
        if(!q) q = e;
        f(t, q, n);
        if(q == e) return n + 1;
        t = q + 1;
    }
}

static inline int count_tokens(const char* b, const char* e, char sep){
    int n = 1;
    for(const char* q = b; (q = static_cast<const char*>(memchr(q, sep, e - q))); ++q) ++n;
    return n;
}

// Same result as std::stoi, without throwing: leading spaces and an optional sign, 0 if nothing parses or out of range
static inline int parse_int_token(const char* b, const char* e){
    while(b < e && isspace(static_cast<unsigned char>(*b))) ++b;
    if(e - b > 1 && *b == '+' && isdigit(static_cast<unsigned char>(b[1]))) ++b;
    int v = 0;
    return std::from_chars(b, e, v).ec == std::errc() ? v : 0;
}

// Same result as std::stof, without throwing (stof is strtof: 0 if nothing parses or on ERANGE)
static inline float parse_float_token(const char* b, const char* e){
    char buf[64];
    std::string big;
    const size_t n = e - b;
    const char* s = buf;
    if(n < sizeof(buf)){ memcpy(buf, b, n); buf[n] = '\0'; }
    else{ big.assign(b, n); s = big.c_str(); }
    char* end;
    errno = 0;
    const float v = strtof(s, &end);
    return (end == s || errno == ERANGE) ? 0.0f : v;
}

/**
 * @brief Builds the read-only lookup tables used by the parallel host parse.
 *
 * Maps CHROM / FILTER names to their codes, and every INFO key to its type code and the index of its
 * host-side column (in_string, alt_int, alt_float or alt_string; -1 when the host has nothing to do).
 */
void vcf_parsed::build_host_lookup(){
    host_chrom.clear();
    host_filter.clear();
    host_info.clear();
    for(const auto& kv : var_columns.chrom_map) host_chrom.emplace(kv.first, kv.second);
    for(const auto& kv : var_columns.filter_map) host_filter.emplace(kv.first, kv.second);
    auto index_of = [](const auto& columns, const std::string& name){
        for(size_t el = 0; el < columns.size(); el++) if(columns[el].name == name) return static_cast<int>(el);
        return -1;
    };
    for(const auto& kv : var_columns.info_map1){
        int el = -1;
        switch(kv.second){
            case STRING:     el = index_of(var_columns.in_string, kv.first); break;
            case INT_ALT:    el = index_of(alt_columns.alt_int, kv.first); break;
            case FLOAT_ALT:  el = index_of(alt_columns.alt_float, kv.first); break;
            case STRING_ALT: el = index_of(alt_columns.alt_string, kv.first); break;
        }
        host_info.emplace(kv.first, host_key{kv.second, el});
    }
}

/**
 * @brief Classifies the fields of a FORMAT template once: what the host does for each position.
 */
format_plan vcf_parsed::make_format_plan(const char* b, const char* e){
    format_plan plan;
    for_each_token(b, e, ':', [&](const char* tb, const char* te, int){
        const std::string key(tb, te - tb);
        format_step step;
        if(key == "GT"){
            // GT with Number=A goes to DF4; GT Number=1 is parsed by the kernel; an undeclared GT is skipped
            if(samp_columns.sample_GT.empty() && FORMAT.hasGT) step.kind = format_step::GT_ALT;
        }else{
            const int code = var_columns.info_code(key);
            const int code1 = var_columns.info_code(key + "1");
            auto index_of = [&](const auto& columns, bool exact){
                for(size_t el = 0; el < columns.size(); el++)
                    if(exact ? columns[el].name == key : format_name_matches(columns[el].name, key)) return static_cast<int>(el);
                return -1;
            };
            if(code == STRING_FORMAT || code1 == STRING_FORMAT){
                step.el = index_of(samp_columns.samp_string, false);
                if(step.el >= 0){ step.kind = format_step::STR; step.numb = samp_columns.samp_string[step.el].numb; }
            }else if(code == INT_FORMAT || code1 == INT_FORMAT || code == FLOAT_FORMAT || code1 == FLOAT_FORMAT){
                // parsed by the kernel
            }else if(code == STRING_FORMAT_ALT){
                step.el = index_of(alt_sample.samp_string, true);
                if(step.el >= 0) step.kind = format_step::STR_ALT;
            }else if(code == INT_FORMAT_ALT){
                step.el = index_of(alt_sample.samp_int, true);
                if(step.el >= 0) step.kind = format_step::INT_ALT;
            }else if(code == FLOAT_FORMAT_ALT){
                step.el = index_of(alt_sample.samp_float, true);
                if(step.el >= 0) step.kind = format_step::FLT_ALT;
            }
        }
        plan.steps.push_back(step);
        if(step.kind != format_step::NONE) plan.used = plan.steps.size();
    });
    return plan;
}

/**
 * @brief Parses the fixed fields and INFO of a record (everything up to FORMAT).
 * @return Pointer to the field after INFO.
 */
const char* vcf_parsed::parse_record_head(const char* p, const char* e, long i, alt_columns_df* tmp_alt, int* tmp_num_alt){
    if(p < e && *p == '\n') ++p;

    // CHROM (chrom_map is filled before the parallel parse, read only here)
    const char* q = field_end(p, e);
    auto chrom_it = host_chrom.find(std::string_view(p, q - p));
    var_columns.chrom[i] = (chrom_it != host_chrom.end()) ? chrom_it->second : static_cast<unsigned char>(0);
    p = next_field(q, e);
    // POS: on device
    p = next_field(field_end(p, e), e);
    // ID
    q = field_end(p, e);
    var_columns.id[i].assign(p, q - p);
    p = next_field(q, e);
    // REF
    q = field_end(p, e);
    var_columns.ref[i].assign(p, q - p);
    p = next_field(q, e);
    // ALT: one DF2 row per allele
    q = field_end(p, e);
    const int base = *tmp_num_alt;
    const int local_alt = count_tokens(p, q, ',');
    ensure_alt_capacity(tmp_alt, base + local_alt);
    for_each_token(p, q, ',', [&](const char* tb, const char* te, int y){
        tmp_alt->alt[base + y].assign(tb, te - tb);
        tmp_alt->alt_id[base + y] = (char)y;
        tmp_alt->var_id[base + y] = i;
    });
    p = next_field(q, e);
    // QUAL: on device
    p = next_field(field_end(p, e), e);
    // FILTER
    q = field_end(p, e);
    auto filter_it = host_filter.find(std::string_view(p, q - p));
    var_columns.filter[i] = (filter_it != host_filter.end()) ? filter_it->second : static_cast<char>(0);
    p = next_field(q, e);

    // INFO: key=value entries separated by ';' (an entry with no or more than one '=' is ignored)
    q = field_end(p, e);
    for_each_token(p, q, ';', [&](const char* tb, const char* te, int){
        const char* eq = static_cast<const char*>(memchr(tb, '=', te - tb));
        if(!eq || memchr(eq + 1, '=', te - eq - 1)) return;
        auto it = host_info.find(std::string_view(tb, eq - tb));
        if(it == host_info.end() || it->second.el < 0) return;
        const int el = it->second.el;
        const char* vb = eq + 1;
        switch(it->second.code){
            case STRING:
                var_columns.in_string[el].i_string[i].assign(vb, te - vb);
                break;
            case INT_ALT:
            case FLOAT_ALT:
            case STRING_ALT: {
                // One value per ALT; a missing value ('.' is a single token) gives 0 / ""
                int y = 0;
                for_each_token(vb, te, ',', [&](const char* vtb, const char* vte, int k){
                    if(k >= local_alt) return;
                    if(it->second.code == INT_ALT) tmp_alt->alt_int[el].i_int[base + k] = parse_int_token(vtb, vte);
                    else if(it->second.code == FLOAT_ALT) tmp_alt->alt_float[el].i_float[base + k] = (__half)parse_float_token(vtb, vte);
                    else tmp_alt->alt_string[el].i_string[base + k].assign(vtb, vte - vtb);
                    y = k + 1;
                });
                for(; y < local_alt; y++){
                    if(it->second.code == INT_ALT) tmp_alt->alt_int[el].i_int[base + y] = 0;
                    else if(it->second.code == FLOAT_ALT) tmp_alt->alt_float[el].i_float[base + y] = (__half)0.0f;
                    else tmp_alt->alt_string[el].i_string[base + y] = "";
                }
                break;
            }
        }
    });
    *tmp_num_alt = base + local_alt;
    return next_field(q, e);
}

/**
    * @brief Parses a VCF line and populates variant columns data.
    *
    * This function processes a single VCF line (from index @p start to @p end). It extracts the host-side
    * fields (chromosome, variant ID, reference allele, alternative alleles, filter and the INFO fields the
    * kernel does not parse). The alternative alleles are stored in the provided alt_columns_df structure
    * (@p tmp_alt), and @p tmp_num_alt counts the alternative alleles processed.
    *
    * @param line Pointer to the VCF body.
    * @param start The starting index of the line within the file.
    * @param end The ending index of the line within the file.
    * @param i The index (row number) corresponding to the current variant.
    * @param tmp_alt Pointer to an alt_columns_df structure for storing alternative allele data.
    * @param tmp_num_alt Pointer to an integer tracking the current number of alternative alleles processed.
    */
void vcf_parsed::get_vcf_line_in_var_columns(char *line, long start, long end, long i, alt_columns_df* tmp_alt, int *tmp_num_alt)
{
    parse_record_head(line + start, line + end, i, tmp_alt, tmp_num_alt);
}

/**
    * @brief Parses a VCF line with sample data.
    *
    * Parses the record head like get_vcf_line_in_var_columns, then the FORMAT template and the samples.
    * The template is classified once (format_plan, cached per chunk in @p plans): when no FORMAT field
    * needs the host (e.g. only GT Number=1 and numeric fields, all parsed by the kernel) the samples are
    * not scanned at all, otherwise each sample is split only up to the last field the host needs.
    *
    * @param line Pointer to the VCF body.
    * @param start The starting index of the line within the file.
    * @param end The ending index of the line within the file.
    * @param i The index (row number) corresponding to the current variant.
    * @param tmp_alt Pointer to an alt_columns_df structure for storing alternative allele data.
    * @param tmp_num_alt Pointer to an integer tracking the number of alternative alleles processed.
    * @param sample Pointer to a sample_columns_df structure for storing sample-specific data.
    * @param FORMAT Pointer to a header_element structure describing the FORMAT fields.
    * @param tmp_num_alt_format Pointer to an integer tracking the number of formatted alternative entries processed.
    * @param tmp_alt_format Pointer to an alt_format_df structure for storing formatted sample data.
    * @param plans Per-chunk cache of the classified FORMAT templates.
    */
void vcf_parsed::get_vcf_line_in_var_columns_format(char *line, long start, long end, long i, alt_columns_df* tmp_alt, int *tmp_num_alt, sample_columns_df* sample, header_element* FORMAT, int *tmp_num_alt_format, alt_format_df* tmp_alt_format, format_plan_cache* plans)
{
    const char* e = line + end;
    const char* p = parse_record_head(line + start, e, i, tmp_alt, tmp_num_alt);

    // FORMAT template
    const char* q = field_end(p, e);
    const std::string_view template_key(p, q - p); // points into filestring, which outlives the chunk's cache
    auto plan_it = plans->find(template_key);
    if(plan_it == plans->end()){
        plan_it = plans->emplace(template_key, make_format_plan(p, q)).first;
    }
    const format_plan& plan = plan_it->second;
    if(plan.used == 0) return; // nothing for the host in this record's samples
    p = q < e ? q + 1 : e;

    // Every separator opens one more sample, so an empty last sample (a trailing tab) is still parsed;
    // sample columns missing at the end of the record read as empty samples, as on the CPU and the kernel
    const unsigned int n_samp = sample->numSample;
    for(unsigned int samp = 0; samp < n_samp; samp++){
        q = field_end(p, e);
        const char* tb = p;
        for(size_t j = 0; j < plan.used && tb <= q; j++){
            const char* te = static_cast<const char*>(memchr(tb, ':', q - tb));
            if(!te) te = q;
            const format_step& step = plan.steps[j];
            const size_t cell = i*sample->numSample + samp;
            switch(step.kind){
                case format_step::NONE: break;
                case format_step::STR:
                    if(step.numb == 1){
                        sample->samp_string[step.el].i_string[cell].assign(tb, te - tb);
                    }else{
                        // Fixed Number>1: columns <ID>0..<ID>k-1; a '.' value has a single token
                        int got = 0;
                        for_each_token(tb, te, ',', [&](const char* vb, const char* ve, int k){
                            if(k < step.numb){ sample->samp_string[step.el + k].i_string[cell].assign(vb, ve - vb); got = k + 1; }
                        });
                        for(int k = got; k < step.numb; k++) sample->samp_string[step.el + k].i_string[cell] = "";
                    }
                    break;
                default: {
                    // Number=A (and GT Number=A): one DF4 row per value
                    const int base = *tmp_num_alt_format;
                    const int local_alt = count_tokens(tb, te, ',');
                    ensure_alt_format_capacity(tmp_alt_format, base + local_alt);
                    for_each_token(tb, te, ',', [&](const char* vb, const char* ve, int y){
                        tmp_alt_format->var_id[base + y] = static_cast<unsigned int>(i); // var_number is still being copied back from the device
                        tmp_alt_format->samp_id[base + y] = samp;
                        tmp_alt_format->alt_id[base + y] = (char)y;
                        switch(step.kind){
                            case format_step::GT_ALT: {
                                auto gt = tmp_alt_format->GTMap.find(std::string(vb, ve - vb));
                                tmp_alt_format->sample_GT.GT[base + y] = (gt != tmp_alt_format->GTMap.end()) ? gt->second : (char)0;
                                break;
                            }
                            case format_step::STR_ALT: tmp_alt_format->samp_string[step.el].i_string[base + y].assign(vb, ve - vb); break;
                            case format_step::INT_ALT: tmp_alt_format->samp_int[step.el].i_int[base + y] = parse_int_token(vb, ve); break;
                            case format_step::FLT_ALT: tmp_alt_format->samp_float[step.el].i_float[base + y] = (__half)parse_float_token(vb, ve); break;
                            default: break;
                        }
                    });
                    *tmp_num_alt_format = base + local_alt;
                }
            }
            if(te == q) break;
            tb = te + 1;
        }
        p = q < e ? q + 1 : e;
    }
}


#endif