#include "VCFReconstructorGPU.h"
#include "CudaUtils.cuh"
#include "Utils.h"
#include "Maps.h"
#include <sstream>
#include <omp.h>
#include <chrono>
#include <cub/cub.cuh>

//Constructor
VCFReconstructorGPU::VCFReconstructorGPU(const std::string& output_vcf_path,
                                         const std::string& header_text)
    : output_vcf_path(output_vcf_path),
      header_text(header_text),
      d_output(nullptr),
      d_line_lens(nullptr)
{}

/**
 * @brief Safely frees device memory and nullifies the pointer.
 * 
 * This templated utility ensures that `cudaFree` is only called on valid 
 * (non-null) pointers. By taking the pointer by reference (`T*&`), it 
 * automatically sets the original pointer to `nullptr` after freeing, 
 * preventing dangling pointer bugs or double-free errors during the 
 * cleanup phase of the GPU pipeline.
 * 
 * @tparam T The data type of the device pointer.
 * @param d_ptr Reference to the device pointer to be freed.
 */
template <typename T>
inline void safe_cuda_free(T*& d_ptr) {
    if (d_ptr != nullptr) {
        gpuErrchk( cudaFree(d_ptr));
        d_ptr = nullptr;
    }
}

/**
 * This method iterates through all device pointers (d_*) associated with the 
 * sparse output buffer, stream compaction arrays, and the flattened shadow 
 * DataFrames (DF1, DF2, DF3). It relies on the `safe_cuda_free` utility to 
 * avoid double-free errors.
 */
void VCFReconstructorGPU::freeDevice(){
    // --- Sparse Output Buffers ---
    safe_cuda_free(d_output);
    safe_cuda_free(d_line_lens);
    safe_cuda_free(d_output_offsets); 

    // --- DF1 (Core Variants & INFO) ---
    safe_cuda_free(d_df1.var_number);
    safe_cuda_free(d_df1.chrom);
    safe_cuda_free(d_df1.pos);
    safe_cuda_free(d_df1.id_data);
    safe_cuda_free(d_df1.id_offsets);
    safe_cuda_free(d_df1.ref_data);
    safe_cuda_free(d_df1.ref_offsets);
    safe_cuda_free(d_df1.qual);
    safe_cuda_free(d_df1.filter);
    safe_cuda_free(d_df1.in_int);
    safe_cuda_free(d_df1.in_float);
    safe_cuda_free(d_df1.in_flag);
    safe_cuda_free(d_df1.in_string_data);
    safe_cuda_free(d_df1.in_string_offsets);
    
    // --- DF2 (Alternative Alleles & INFO-ALT) ---
    safe_cuda_free(d_df2.var_id);
    safe_cuda_free(d_df2.alt_data);
    safe_cuda_free(d_df2.alt_offsets);
    safe_cuda_free(d_df2.alt_start);
    safe_cuda_free(d_df2.alt_count);
    safe_cuda_free(d_df2.alt_int);
    safe_cuda_free(d_df2.alt_float);
    safe_cuda_free(d_df2.alt_string_data);
    safe_cuda_free(d_df2.alt_string_offsets);

    // --- DF3/DF4 (Packed Sample Data) ---
    safe_cuda_free(d_df3.sample_data);
    safe_cuda_free(d_df3.sample_offsets);

    // NOTE: d_compacted is already freed inside writeChunk() 
    // to keep its lifecycle local to the I/O operation.
}

//Destructor
VCFReconstructorGPU::~VCFReconstructorGPU() {
    freeDevice(); // Clean up VRAM

    // Clean up host-side dynamically allocated compaction buffer
    if (h_compacted) safe_cuda_free(h_compacted);
}

// Device Allocation & Data Transfer Setup
void VCFReconstructorGPU::allocateDevice(const var_columns_df& df1,
                                          const alt_columns_df& df2,
                                          const sample_columns_df& df3,
                                          int chunk_size,
                                          int chunk_start,
                                          int chunk_end,
                                          int df2_start) {

    // --- Output buffers ---
    // Allocate space for the sparse output buffer where each row takes up MAX_LINE_LEN chars.
    gpuErrchk( cudaMalloc((void**)&d_output,    chunk_size * MAX_LINE_LEN * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_line_lens, chunk_size * sizeof(unsigned int)));
    gpuErrchk(cudaMalloc((void**)&d_output_offsets, (chunk_size + 1) * sizeof(unsigned int)));

    // Note: d_compacted is dynamically allocated in writeChunk() because its size
    // relies on the sum of line_lens, which is known only after the prefix scan.

    // --- DF1 scalar fields ---
    gpuErrchk( cudaMalloc((void**)&d_df1.var_number, chunk_size * sizeof(unsigned int)));
    gpuErrchk( cudaMalloc((void**)&d_df1.chrom,      chunk_size * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_df1.pos,        chunk_size * sizeof(unsigned int)));
    gpuErrchk( cudaMalloc((void**)&d_df1.qual,       chunk_size * sizeof(__half)));
    gpuErrchk( cudaMalloc((void**)&d_df1.filter,     chunk_size * sizeof(char)));

    // --- DF1 ID (Variable-length strings) ---
    // Pre-calculate total memory required to hold all null-terminated ID strings in this chunk
    size_t id_total_chars = 0;
    for (int i = chunk_start; i < chunk_end; i++)
        id_total_chars += df1.id[i].size() + 1;
    gpuErrchk( cudaMalloc((void**)&d_df1.id_data,    id_total_chars * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_df1.id_offsets, chunk_size * sizeof(unsigned int)));

    // --- DF1 REF (Variable-length strings) ---
    size_t ref_total_chars = 0;
    for (int i = chunk_start; i < chunk_end; i++)
        ref_total_chars += df1.ref[i].size() + 1;
    gpuErrchk( cudaMalloc((void**)&d_df1.ref_data,    ref_total_chars * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_df1.ref_offsets, chunk_size * sizeof(unsigned int)));

    // --- DF1 INFO fields ---
    d_df1.num_int_fields   = df1.in_int.size();
    d_df1.num_float_fields = df1.in_float.size();
    d_df1.num_flag_fields  = df1.in_flag.size();
    gpuErrchk( cudaMalloc((void**)&d_df1.in_int,      d_df1.num_int_fields   * chunk_size * sizeof(int)));
    gpuErrchk( cudaMalloc((void**)&d_df1.in_float,    d_df1.num_float_fields * chunk_size * sizeof(__half)));
    gpuErrchk( cudaMalloc((void**)&d_df1.in_flag,     d_df1.num_flag_fields  * chunk_size * sizeof(bool)));

    // --- DF1 INFO Strings ---
    d_df1.num_string_fields = df1.in_string.size();
    size_t info_str_total_chars = 0;
    for (int f = 0; f < d_df1.num_string_fields; f++) {
        for (int i = chunk_start; i < chunk_end; i++) {
            info_str_total_chars += df1.in_string[f].i_string[i].size() + 1; // +1 for the '\0' terminator
        }
    }
    if (info_str_total_chars > 0) {
        gpuErrchk(cudaMalloc((void**)&d_df1.in_string_data, info_str_total_chars * sizeof(char)));
    }
    gpuErrchk(cudaMalloc((void**)&d_df1.in_string_offsets, d_df1.num_string_fields * chunk_size * sizeof(unsigned int)));


    // --- DF2 (Alternative alleles count & allocation) ---
    int df2_count = 0;
    size_t alt_total_chars = 0;
    for (size_t j = df2_start; j < df2.var_id.size() && (int)df2.var_id[j] < chunk_end; j++) {
        alt_total_chars += df2.alt[j].size() + 1;
        df2_count++;
    }

    d_df2.num_entries = df2_count;
    d_df2.num_alt_int_fields   = df2.alt_int.size();
    d_df2.num_alt_float_fields = df2.alt_float.size();

    gpuErrchk( cudaMalloc((void**)&d_df2.var_id,       df2_count * sizeof(unsigned int)));
    gpuErrchk( cudaMalloc((void**)&d_df2.alt_data,         alt_total_chars * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_df2.alt_offsets,      df2_count * sizeof(unsigned int)));

    gpuErrchk(cudaMalloc((void**)&d_df2.alt_start, chunk_size * sizeof(unsigned int)));
    gpuErrchk(cudaMalloc((void**)&d_df2.alt_count, chunk_size * sizeof(unsigned int)));

    gpuErrchk( cudaMalloc((void**)&d_df2.alt_int,          d_df2.num_alt_int_fields   * df2_count * sizeof(int)));
    gpuErrchk( cudaMalloc((void**)&d_df2.alt_float,        d_df2.num_alt_float_fields * df2_count * sizeof(__half)));

    // --- DF2 INFO-ALT Strings ---
    d_df2.num_alt_string_fields = df2.alt_string.size();
    size_t alt_info_str_total_chars = 0;
    for (int f = 0; f < d_df2.num_alt_string_fields; f++) {
        for (int i = df2_start; i < df2_start + df2_count; i++) {
            alt_info_str_total_chars += df2.alt_string[f].i_string[i].size() + 1;
        }
    }
    if (alt_info_str_total_chars > 0) {
        gpuErrchk(cudaMalloc((void**)&d_df2.alt_string_data, alt_info_str_total_chars * sizeof(char)));
    }
    gpuErrchk(cudaMalloc((void**)&d_df2.alt_string_offsets, d_df2.num_alt_string_fields * df2_count * sizeof(unsigned int)));

    // --- DF3 (Sample Strings upper bound) ---
    // Calculate a safe upper bound size for sample strings to prevent OOM errors.
    size_t max_sample_chars = (size_t)chunk_size * df3.numSample * MAX_SAMPLE_STRING_LEN;

    gpuErrchk( cudaMalloc(&d_df3.sample_data,    max_sample_chars * sizeof(char)));
    gpuErrchk( cudaMalloc(&d_df3.sample_offsets, chunk_size * df3.numSample * sizeof(unsigned int)));
    d_df3.num_samples = df3.numSample;

}

void VCFReconstructorGPU::uploadToDevice(const var_columns_df& df1,
                                          const alt_columns_df& df2,
                                          const HostBuffers& buffers) {
    int chunk_size  = buffers.chunk_size;
    int chunk_start = buffers.chunk_start;
    int df2_start   = buffers.df2_start;
    int df2_count   = buffers.df2_count;

    // --- DF1 Direct Copies (Fixed-size arrays) ---
    gpuErrchk(cudaMemcpy(d_df1.var_number, df1.var_number.data() + chunk_start,
                         chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df1.chrom, df1.chrom.data() + chunk_start,
                         chunk_size * sizeof(char), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df1.pos, df1.pos.data() + chunk_start,
                         chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df1.qual, df1.qual.data() + chunk_start,
                         chunk_size * sizeof(__half), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df1.filter, df1.filter.data() + chunk_start,
                         chunk_size * sizeof(char), cudaMemcpyHostToDevice));

    // --- DF1 Staging Copies (Variable-length flattened data) ---
    gpuErrchk(cudaMemcpy(d_df1.id_data, buffers.id_buffer.data(),
                         buffers.id_buffer.size() * sizeof(char), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df1.id_offsets, buffers.id_offsets.data(),
                         chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df1.ref_data, buffers.ref_buffer.data(),
                         buffers.ref_buffer.size() * sizeof(char), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df1.ref_offsets, buffers.ref_offsets.data(),
                         chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));

    // --- DF1 INFO Copies ---
    gpuErrchk(cudaMemcpy(d_df1.in_int, buffers.in_int_buffer.data(),
                         d_df1.num_int_fields * chunk_size * sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df1.in_float, buffers.in_float_buffer.data(),
                         d_df1.num_float_fields * chunk_size * sizeof(__half), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df1.in_flag, buffers.in_flag_buffer.data(),
                         d_df1.num_flag_fields * chunk_size * sizeof(uint8_t), cudaMemcpyHostToDevice));

    // --- DF1 String INFO Copies ---
    if (buffers.in_string_buffer.size() > 0) {
        gpuErrchk(cudaMemcpy(d_df1.in_string_data, buffers.in_string_buffer.data(),
                            buffers.in_string_buffer.size() * sizeof(char), cudaMemcpyHostToDevice));
    }
    gpuErrchk(cudaMemcpy(d_df1.in_string_offsets, buffers.in_string_offsets.data(),
                        d_df1.num_string_fields * chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));

    // --- DF2 Direct and Staging Copies ---
    gpuErrchk(cudaMemcpy(d_df2.var_id, df2.var_id.data() + df2_start,
                         df2_count * sizeof(unsigned int), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_df2.alt_data, buffers.alt_data_buffer.data(),
                         buffers.alt_data_buffer.size() * sizeof(char), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df2.alt_offsets, buffers.alt_data_offsets.data(),
                         df2_count * sizeof(unsigned int), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_df2.alt_start, buffers.alt_start_buf.data(),
                         chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df2.alt_count, buffers.alt_count_buf.data(),
                         chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_df2.alt_int, buffers.alt_int_buffer.data(),
                         d_df2.num_alt_int_fields * df2_count * sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df2.alt_float, buffers.alt_float_buffer.data(),
                         d_df2.num_alt_float_fields * df2_count * sizeof(__half), cudaMemcpyHostToDevice));

    // --- DF2 String INFO-ALT Copies ---
    if (buffers.alt_string_buffer.size() > 0) {
        gpuErrchk(cudaMemcpy(d_df2.alt_string_data, buffers.alt_string_buffer.data(),
                            buffers.alt_string_buffer.size() * sizeof(char), cudaMemcpyHostToDevice));
    }
    gpuErrchk(cudaMemcpy(d_df2.alt_string_offsets, buffers.alt_string_offsets.data(),
                        d_df2.num_alt_string_fields * df2_count * sizeof(unsigned int), cudaMemcpyHostToDevice));

    // --- DF3 Sample Data Copies ---
    gpuErrchk(cudaMemcpy(d_df3.sample_data, buffers.sample_buffer.data(),
                         buffers.sample_buffer.size() * sizeof(char), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df3.sample_offsets, buffers.sample_offsets.data(),
                         buffers.sample_offsets.size() * sizeof(unsigned int), cudaMemcpyHostToDevice));
}

/**
 * @brief Stream compaction kernel to consolidate sparse VCF string buffers.
 * 
 * `reconstructKernel` outputs strings into a sparse 2D buffer (`d_output`) where 
 * each row occupies a fixed, maximum length (`max_line_len`) to prevent data races. 
 * This kernel reads the actual length of each generated string (`line_lens`) and 
 * its designated starting position (`offsets`, calculated via prefix sum) to pack 
 * everything into a contiguous, dense 1D array (`dst`).
 * 
 * @param src         Pointer to the sparse input buffer (d_output).
 * @param line_lens   Array holding the actual string length of each row.
 * @param offsets     Prefix-summed array indicating where each compacted row should start.
 * @param dst         Pointer to the dense destination buffer (d_compacted).
 * @param chunk_size  Total number of variants (rows) processed in this block.
 * @param max_line_len The fixed stride used in the sparse `src` buffer.
 */
__global__ void compactKernel(
    const char* __restrict__ src,       
    const unsigned int* __restrict__ line_lens,
    const unsigned int* __restrict__ offsets, 
    char* __restrict__ dst,               
    int chunk_size,
    int max_line_len)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;

    // Grid-stride loop handles chunks larger than the total number of launched threads
    for (int i = tid; i < chunk_size; i += stride) {
        unsigned int len = line_lens[i];
        unsigned int dst_off = offsets[i];
        const char* s = src + i * max_line_len;
        char* d = dst + dst_off;
        for (unsigned int k = 0; k < len; k++) {
            d[k] = s[k];
        }
    }
}

/**
 * @brief Main execution kernel that constructs VCF formatted strings per variant.
 * 
 * This kernel performs the core text reconstruction using a grid-stride loop. 
 * While each iteration processes a single variant row, a single CUDA thread 
 * may process multiple variants if the chunk size exceeds the total number of 
 * launched threads. It reads the raw data from the shadowed DataFrame structures, 
 * translates integer codes back to strings using the mapped dictionaries, handles 
 * multi-allelic expansions, and concatenates the final VCF string into a sparse buffer.
 * 
 * @param maps       Dictionaries for translating codes to strings (Chrom, Filter).
 * @param df1        Device memory view of the core variants DataFrame.
 * @param df2        Device memory view of the alternative alleles DataFrame.
 * @param df3        Device memory view containing pre-packed format/sample data.
 * @param output     The sparse destination buffer for the generated text.
 * @param line_lens  Output array where the thread records the final length of each generated string.
 * @param chunk_size The total number of variants to process in this execution block.
 */
__global__ void reconstructKernel(
    DeviceMaps maps,
    DeviceVarColumns df1,
    DeviceAltColumns df2,
    DeviceSampleColumns df3,
    char* output, unsigned int* line_lens,
    int chunk_size
){
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = tid; i < chunk_size; i += stride){
        // Pointer to the start of this thread's designated row in the sparse matrix
        char* line = output + i * MAX_LINE_LEN;
        int pos = 0;

        // --- 1. CHROM ---
        unsigned char chrom_code = df1.chrom[i];
        pos += device_strcpy(line + pos, maps.chrom_strings + maps.chrom_offsets[chrom_code]);
        line[pos++] = '\t';

        // --- 2. POS ---
        pos += device_itoa(df1.pos[i], line + pos);
        line[pos++] = '\t';

        // --- 3. ID ---
        pos += device_strcpy(line + pos, df1.id_data + df1.id_offsets[i]);
        line[pos++] = '\t';

        // --- 4. REF ---
        pos += device_strcpy(line + pos, df1.ref_data + df1.ref_offsets[i]);
        line[pos++] = '\t';

        // --- 5. ALT ---
        // Handle multiple alternative alleles separated by commas
        int alt_begin = df2.alt_start[i];
        int alt_end   = alt_begin + df2.alt_count[i];
        bool first_alt = true;
        for (int j = alt_begin; j < alt_end; j++){
            if (!first_alt) line[pos++] = ',';
            first_alt = false;
            pos += device_strcpy(line + pos, df2.alt_data + df2.alt_offsets[j]);
        }
        line[pos++] = '\t';

        // --- 6. QUAL ---
        float q = df1.qual[i];
        pos += device_ftoa(q, line + pos);
        line[pos++] = '\t';

        // --- 7. FILTER ---
        unsigned char filter_code = df1.filter[i];
        pos += device_strcpy(line + pos, maps.filter_strings + maps.filter_offsets[filter_code]);
        line[pos++] = '\t';

        // --- 8. INFO ---
        bool first_info = true;
        
        // Append variant-level Flag fields (e.g., SOMATIC)
        for (int f = 0; f < df1.num_flag_fields; f++){
            if (df1.in_flag[f * chunk_size + i]){
                if(!first_info){
                    line[pos++] = ';';
                }
                first_info = false;
                pos += device_strcpy(line + pos, df1.flag_names + f * MAX_NAME_LEN);
            }
        }

        // Append variant-level Integer fields
        for (int f = 0; f < df1.num_int_fields; f++){
            if (df1.in_int[f * chunk_size + i] != -1){
                if(!first_info){
                    line[pos++] = ';';
                }
                first_info = false;
                pos += device_strcpy(line + pos, df1.int_names + f * MAX_NAME_LEN);
                line[pos++] = '=';
                pos += device_itoa(df1.in_int[f * chunk_size + i], line + pos);
            }
        }

        // Append variant-level Float fields
        for (int f = 0; f < df1.num_float_fields; f++){
            if (__half2float(df1.in_float[f * chunk_size + i]) != -1.0f){
                if(!first_info){
                    line[pos++] = ';';
                }
                first_info = false;
                pos += device_strcpy(line + pos, df1.float_names + f * MAX_NAME_LEN);
                line[pos++] = '=';
                pos += device_ftoa(__half2float(df1.in_float[f * chunk_size + i]), line + pos);
            }
        }

        // --- Append variant-level String fields ---
        for (int f = 0; f < df1.num_string_fields; f++){
            int flat_idx = f * chunk_size + i;
            unsigned int offset = df1.in_string_offsets[flat_idx];
            const char* str_val = df1.in_string_data + offset;
            
            // Ensure the string is neither empty nor the missing value sentinel "."
            if (str_val[0] != '\0' && !(str_val[0] == '.' && str_val[1] == '\0')){
                if(!first_info) line[pos++] = ';';
                first_info = false;
                
                const char* field_name = df1.string_names + f * MAX_NAME_LEN;
                pos += device_strcpy(line + pos, field_name);
                line[pos++] = '=';
                
                bool mapped = false;

                // Decode specific mapped fields (TSA, PolyPhen, CSQ)
                if (device_strcmp(field_name, "TSA") && device_is_numeric(str_val)) {
                    int val = device_atoi(str_val);
                    if (val >= 0 && val < 256) { // Bounds check matching 'it != map.end()'
                        const char* mapped_str = maps.tsa_strings + maps.tsa_offsets[val];
                        if (mapped_str[0] != '\0') {
                            pos += device_strcpy(line + pos, mapped_str);
                            mapped = true;
                        }
                    }
                } else if (device_strcmp(field_name, "PolyPhen") && device_is_numeric(str_val)) {
                    int val = device_atoi(str_val);
                    if (val >= 0 && val < 256) {
                        const char* mapped_str = maps.polyphen_strings + maps.polyphen_offsets[val];
                        if (mapped_str[0] != '\0') {
                            pos += device_strcpy(line + pos, mapped_str);
                            mapped = true;
                        }
                    }
                } else if (device_strcmp(field_name, "CSQ") && device_is_numeric(str_val)) {
                    int val = device_atoi(str_val);
                    if (val >= 0 && val < 256) {
                        const char* mapped_str = maps.csq_strings + maps.csq_offsets[val];
                        if (mapped_str[0] != '\0') {
                            pos += device_strcpy(line + pos, mapped_str);
                            mapped = true;
                        }
                    }
                }
                
                // Fallback (simulating the 'catch' block or unmapped string)
                if (!mapped) {
                    pos += device_strcpy(line + pos, str_val);
                }
            }
        }

        // Append allele-specific Integer fields from DF2 (comma-separated if multi-allelic)
        for (int f = 0; f < df2.num_alt_int_fields; f++){
            // Validate if the group contains any actual data
            bool any_valid = false;
            for (int j = alt_begin; j < alt_end; j++){
                if (df2.alt_int[f * df2.num_entries + j] != -1){
                    any_valid = true;
                    break;
                }
            }
            if (!any_valid) continue;

            if (!first_info) line[pos++] = ';';
            first_info = false;

            pos += device_strcpy(line + pos, df2.alt_int_names + f * MAX_NAME_LEN);
            line[pos++] = '=';

            bool first_val = true;
            for (int j = alt_begin; j < alt_end; j++){
                if (!first_val) line[pos++] = ',';
                first_val = false;
                int v = df2.alt_int[f * df2.num_entries + j];
                if (v == -1){
                    line[pos++] = '.';
                } else {
                    pos += device_itoa(v, line + pos);
                }
            }
        }

        // Append allele-specific Float fields from DF2
        for (int f = 0; f < df2.num_alt_float_fields; f++){
            bool any_valid = false;
            for (int j = alt_begin; j < alt_end; j++){
                if (__half2float(df2.alt_float[f * df2.num_entries + j]) != -1.0f){
                    any_valid = true;
                    break;
                }
            }
            if (!any_valid) continue;

            if (!first_info) line[pos++] = ';';
            first_info = false;

            pos += device_strcpy(line + pos, df2.alt_float_names + f * MAX_NAME_LEN);
            line[pos++] = '=';

            bool first_val = true;
            for (int j = alt_begin; j < alt_end; j++){
                if (!first_val) line[pos++] = ',';
                first_val = false;
                float v = __half2float(df2.alt_float[f * df2.num_entries + j]);
                if (v == -1.0f){
                    line[pos++] = '.';
                } else {
                    pos += device_ftoa(v, line + pos);
                }
            }
        }

        // --- Append allele-specific String fields from DF2 ---
        for (int f = 0; f < df2.num_alt_string_fields; f++){
            bool any_valid = false;
            
            // First pass: check if there is at least one valid data point in the multi-allelic group
            for (int j = alt_begin; j < alt_end; j++){
                int flat_idx = f * df2.num_entries + j;
                const char* str_val = df2.alt_string_data + df2.alt_string_offsets[flat_idx];
                if (str_val[0] != '\0' && !(str_val[0] == '.' && str_val[1] == '\0')){
                    any_valid = true;
                    break;
                }
            }
            if (!any_valid) continue;

            if (!first_info) line[pos++] = ';';
            first_info = false;

            pos += device_strcpy(line + pos, df2.alt_string_names + f * MAX_NAME_LEN);
            line[pos++] = '=';

            bool first_val = true;
            for (int j = alt_begin; j < alt_end; j++){
                if (!first_val) line[pos++] = ',';
                first_val = false;
                
                int flat_idx = f * df2.num_entries + j;
                const char* str_val = df2.alt_string_data + df2.alt_string_offsets[flat_idx];
                
                // Write the missing indicator if the string is empty or already a sentinel
                if (str_val[0] == '\0' || (str_val[0] == '.' && str_val[1] == '\0')){
                    line[pos++] = '.';
                } else {
                    pos += device_strcpy(line + pos, str_val);
                }
            }
        }

        // If no INFO fields were present, write the missing indicator
        if(first_info){
            line[pos++] = '.';
        }

        // --- 9. FORMAT & SAMPLES ---
        if (df3.num_samples > 0) {
            line[pos++] = '\t';
            pos += device_strcpy(line + pos, df3.format_str);

            for (int s = 0; s < df3.num_samples; s++) {
                line[pos++] = '\t';
                int idx = i * df3.num_samples + s;
                // Copy the pre-packed string for this specific variant-sample combination
                pos += device_strcpy(line + pos, df3.sample_data + df3.sample_offsets[idx]);
            }
        }

        // Finalize the row and record its total length
        line[pos++] = '\n';
        line_lens[i] = pos;
    }
}

void VCFReconstructorGPU::writeChunk(int num_variants, std::ofstream& out){
    // 1. Prefix scan of line_lens on d_line_lens -> d_output_offsets
    //    We use CUB ExclusiveSum: offset[0]=0, offset[i] = sum(line_lens[0..i-1])
    void* d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;

    // First call: only to get temp_storage_bytes required by CUB
    cub::DeviceScan::ExclusiveSum(
        d_temp_storage, temp_storage_bytes,
        d_line_lens, d_output_offsets, num_variants);

    gpuErrchk(cudaMalloc(&d_temp_storage, temp_storage_bytes));

    // Second call: actual prefix scan
    cub::DeviceScan::ExclusiveSum(
        d_temp_storage, temp_storage_bytes,
        d_line_lens, d_output_offsets, num_variants);

    // 2. Calculate total_bytes = last_offset + last_len
    //    Retrieve the last two elements from the device to compute the exact buffer size
    unsigned int last_offset, last_len;
    gpuErrchk(cudaMemcpy(&last_offset, d_output_offsets + (num_variants - 1),
                         sizeof(unsigned int), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(&last_len, d_line_lens + (num_variants - 1),
                         sizeof(unsigned int), cudaMemcpyDeviceToHost));
    size_t total_bytes = last_offset + last_len;

    // 3. Allocate d_compacted to the exact required size
    gpuErrchk(cudaMalloc((void**)&d_compacted, total_bytes * sizeof(char)));

    // 3. Allocate d_compacted to the exact required size
    int block_size = 32;
    int num_blocks = (num_variants + block_size - 1) / block_size;
    compactKernel<<<num_blocks, block_size>>>(
        d_output, d_line_lens, d_output_offsets,
        d_compacted, num_variants, MAX_LINE_LEN);
    gpuErrchk(cudaPeekAtLastError());

    // 5. Grow the host buffer if necessary
    if (total_bytes > h_compacted_capacity) {
        if (h_compacted) free(h_compacted);
        h_compacted_capacity = total_bytes + total_bytes / 4; // 25% margin
        h_compacted = (char*)malloc(h_compacted_capacity);
    }

    // 6. Single D2H copy of the compacted buffer
    gpuErrchk(cudaMemcpy(h_compacted, d_compacted,
                         total_bytes * sizeof(char), cudaMemcpyDeviceToHost));

    // 7. Single write to file
    out.write(h_compacted, total_bytes);

    // 8. Cleanup temporary buffers
    gpuErrchk(cudaFree(d_temp_storage));
    gpuErrchk(cudaFree(d_compacted));
    d_compacted = nullptr;
}

/**
 * @brief Parses and groups format fields to determine their cardinality.
 * 
 * Groups fields based on their base name (stripping trailing digits) and 
 * cross-references them with the parsed FORMAT numbers from the VCF header.
 * 
 * @tparam FieldVec The type of the vector containing the fields.
 * @param fields The vector of fields to group.
 * @param format_numbers Map associating FORMAT IDs to their 'Number' value.
 * @return std::vector<VCFReconstructorGPU::GroupInfo> A vector describing the grouped fields and their rules.
 */
template <typename FieldVec>
static std::vector<VCFReconstructorGPU::GroupInfo> buildGroups(
    const FieldVec& fields,
    const std::map<std::string, std::string>& format_numbers)
{
    std::vector<VCFReconstructorGPU::GroupInfo> out;
    size_t col = 0;
    while (col < fields.size()) {
        std::string base = stripTrailingDigits(fields[col].name);
        size_t end = col;
        while (end < fields.size() &&
               stripTrailingDigits(fields[end].name) == base) {
            end++;
        }
        VCFReconstructorGPU::GroupInfo g;
        g.col_start = col;
        g.col_end   = end;
        auto it = format_numbers.find(base);
        std::string num = (it != format_numbers.end()) ? it->second : "R";
        if      (num == "A") g.number_kind = 1;
        else if (num == "R") g.number_kind = 0;
        else                 g.number_kind = 2;
        out.push_back(g);
        col = end;
    }
    return out;
}

void VCFReconstructorGPU::run(const var_columns_df& df1,
         const alt_columns_df& df2,
         const sample_columns_df& df3,
         const alt_format_df& df4){

    // Initialize mappings and formatting rules based on the header
    buildInverseMaps(df1);
    buildSampleNames(df3);
    format_numbers = parseFormatNumbers(header_text);

    df3_int_groups   = buildGroups(df3.samp_int,   format_numbers);
    df3_float_groups = buildGroups(df3.samp_float, format_numbers);
    df4_int_groups   = buildGroups(df4.samp_int,   format_numbers);
    df4_float_groups = buildGroups(df4.samp_float, format_numbers);

    // Open file and write initial VCF headers
    std::ofstream vcf_file(output_vcf_path);
    if (!vcf_file.is_open()) {
        throw std::runtime_error("Error opening output VCF file: " + output_vcf_path);
    }

    vcf_file << header_text;
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    if (df3.numSample > 0) {
        vcf_file << "\tFORMAT";
        for (const auto& name : ordered_samp_names) {
            vcf_file << "\t" << name;
        }
    }
    vcf_file << "\n";

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Construct the FORMAT string
    gt_in_df3 = !df3.sample_GT.empty();
    gt_in_df4 = !gt_in_df3 && !df4.sample_GT.GT.empty();
    has_gt    = gt_in_df3 || gt_in_df4;

    format_str = "";
    if (df3.numSample > 0) {
        if (has_gt) format_str = "GT";

        // Group DF3 integer fields
        // E.g., if we see AD0, AD1, we only append "AD" to the format string once
        std::string prev_base = "";
        for (size_t i = 0; i < df3.samp_int.size(); i++) {
            std::string base = stripTrailingDigits(df3.samp_int[i].name);
            if (base != prev_base) {
                if (!format_str.empty()) format_str += ":";
                format_str += base;
                prev_base = base;
            }
        }

        // Group DF3 float fields
        prev_base = "";
        for (size_t i = 0; i < df3.samp_float.size(); i++) {
            std::string base = stripTrailingDigits(df3.samp_float[i].name);
            if (base != prev_base) {
                if (!format_str.empty()) format_str += ":";
                format_str += base;
                prev_base = base;
            }
        }

        // DF3 string fields are usually singular (no grouping needed)
        for (size_t i = 0; i < df3.samp_string.size(); i++) {
            if (!format_str.empty()) format_str += ":";
            format_str += df3.samp_string[i].name;
        }

        // Group DF4 integer fields (per-allele format fields)
        prev_base = "";
        for (size_t i = 0; i < df4.samp_int.size(); i++) {
            std::string base = stripTrailingDigits(df4.samp_int[i].name);
            if (base != prev_base) {
                if (!format_str.empty()) format_str += ":";
                format_str += base;
                prev_base = base;
            }
        }

        // Group DF4 float fields
        prev_base = "";
        for (size_t i = 0; i < df4.samp_float.size(); i++) {
            std::string base = stripTrailingDigits(df4.samp_float[i].name);
            if (base != prev_base) {
                if (!format_str.empty()) format_str += ":";
                format_str += base;
                prev_base = base;
            }
        }

        // DF4 string fields
        for (size_t i = 0; i < df4.samp_string.size(); i++) {
            if (!format_str.empty()) format_str += ":";
            format_str += df4.samp_string[i].name;
        }
    }

    gpuErrchk( cudaMalloc(&d_df3.format_str, (format_str.size() + 1) * sizeof(char)));
    gpuErrchk( cudaMemcpy(d_df3.format_str, format_str.c_str(), format_str.size() + 1, cudaMemcpyHostToDevice));
    d_df3.format_str_len = format_str.size();
    d_df3.num_samples    = df3.numSample;

    // ====== Dataset Constants (Allocated and copied once) ======
    {
        // --- Helper per allocazione sicura ---
        auto safe_alloc_and_copy = [](void** device_ptr, const std::vector<char>& buffer) {
            if (!buffer.empty()) {
                gpuErrchk(cudaMalloc(device_ptr, buffer.size()));
                gpuErrchk(cudaMemcpy(*device_ptr, buffer.data(), buffer.size(), cudaMemcpyHostToDevice));
            } else {
                *device_ptr = nullptr; // Fondamentale per evitare crash in free
            }
        };

        // Preparazione buffer nomi (DF1)
        std::vector<char> int_names_buf(df1.in_int.size() * MAX_NAME_LEN, '\0');
        for (size_t f = 0; f < df1.in_int.size(); f++) strncpy(&int_names_buf[f * MAX_NAME_LEN], df1.in_int[f].name.c_str(), MAX_NAME_LEN - 1);
        
        std::vector<char> float_names_buf(df1.in_float.size() * MAX_NAME_LEN, '\0');
        for (size_t f = 0; f < df1.in_float.size(); f++) strncpy(&float_names_buf[f * MAX_NAME_LEN], df1.in_float[f].name.c_str(), MAX_NAME_LEN - 1);
        
        std::vector<char> flag_names_buf(df1.in_flag.size() * MAX_NAME_LEN, '\0');
        for (size_t f = 0; f < df1.in_flag.size(); f++) strncpy(&flag_names_buf[f * MAX_NAME_LEN], df1.in_flag[f].name.c_str(), MAX_NAME_LEN - 1);

        std::vector<char> string_names_buf(df1.in_string.size() * MAX_NAME_LEN, '\0');
        for (size_t f = 0; f < df1.in_string.size(); f++) strncpy(&string_names_buf[f * MAX_NAME_LEN], df1.in_string[f].name.c_str(), MAX_NAME_LEN - 1);

        // Preparazione buffer nomi (DF2)
        std::vector<char> alt_int_names_buf(df2.alt_int.size() * MAX_NAME_LEN, '\0');
        for (size_t f = 0; f < df2.alt_int.size(); f++) strncpy(&alt_int_names_buf[f * MAX_NAME_LEN], df2.alt_int[f].name.c_str(), MAX_NAME_LEN - 1);
        
        std::vector<char> alt_float_names_buf(df2.alt_float.size() * MAX_NAME_LEN, '\0');
        for (size_t f = 0; f < df2.alt_float.size(); f++) strncpy(&alt_float_names_buf[f * MAX_NAME_LEN], df2.alt_float[f].name.c_str(), MAX_NAME_LEN - 1);

        std::vector<char> alt_string_names_buf(df2.alt_string.size() * MAX_NAME_LEN, '\0');
        for (size_t f = 0; f < df2.alt_string.size(); f++) strncpy(&alt_string_names_buf[f * MAX_NAME_LEN], df2.alt_string[f].name.c_str(), MAX_NAME_LEN - 1);

        // Esecuzione allocazioni protette
        safe_alloc_and_copy((void**)&d_df1.int_names,   int_names_buf);
        safe_alloc_and_copy((void**)&d_df1.float_names, float_names_buf);
        safe_alloc_and_copy((void**)&d_df1.flag_names,  flag_names_buf);
        safe_alloc_and_copy((void**)&d_df1.string_names, string_names_buf);

        safe_alloc_and_copy((void**)&d_df2.alt_int_names,   alt_int_names_buf);
        safe_alloc_and_copy((void**)&d_df2.alt_float_names, alt_float_names_buf);
        safe_alloc_and_copy((void**)&d_df2.alt_string_names, alt_string_names_buf);

        // Update counters
        d_df1.num_int_fields    = df1.in_int.size();
        d_df1.num_float_fields  = df1.in_float.size();
        d_df1.num_flag_fields   = df1.in_flag.size();
        d_df1.num_string_fields = df1.in_string.size();
        
        d_df2.num_alt_int_fields    = df2.alt_int.size();
        d_df2.num_alt_float_fields  = df2.alt_float.size();
        d_df2.num_alt_string_fields = df2.alt_string.size();
    }

    // Inverse maps buffers
    {
        auto safe_map_alloc = [](void** str_ptr, void** off_ptr, const std::vector<char>& s_buf, const std::vector<unsigned int>& o_buf) {
            gpuErrchk(cudaMalloc(off_ptr, 256 * sizeof(unsigned int)));
            gpuErrchk(cudaMemcpy(*off_ptr, o_buf.data(), 256 * sizeof(unsigned int), cudaMemcpyHostToDevice));

            if (!s_buf.empty()) {
                gpuErrchk(cudaMalloc(str_ptr, s_buf.size()));
                gpuErrchk(cudaMemcpy(*str_ptr, s_buf.data(), s_buf.size(), cudaMemcpyHostToDevice));
            } else {
                *str_ptr = nullptr;
            }
        };

        std::vector<char> chrom_strings_buf;
        std::vector<unsigned int> chrom_offsets_buf(256, 0);
        unsigned int offset = 0;
        for (const auto& pair : inv_chrom_map) {
            int idx = static_cast<int>(pair.first);
            if (idx >= 0 && idx < 256) {
                chrom_offsets_buf[idx] = offset;
                for (char c : pair.second) chrom_strings_buf.push_back(c);
                chrom_strings_buf.push_back('\0');
                offset += pair.second.size() + 1;
            }
        }

        std::vector<char> filter_strings_buf;
        std::vector<unsigned int> filter_offsets_buf(256, 0);
        offset = 0;
        for (const auto& pair : inv_filter_map) {
            int idx = static_cast<int>(pair.first);
            if (idx >= 0 && idx < 256) {
                filter_offsets_buf[idx] = offset;
                for (char c : pair.second) filter_strings_buf.push_back(c);
                filter_strings_buf.push_back('\0');
                offset += pair.second.size() + 1;
            }
        }

        // --- PolyPhen Map ---
        std::vector<char> polyphen_strings_buf;
        std::vector<unsigned int> polyphen_offsets_buf(256, 0);
        unsigned int offset_poly = 0;
        for (const auto& pair : inv_polyphen_map) {
            int idx = static_cast<int>(pair.first);
            if (idx >= 0 && idx < 256) {
                polyphen_offsets_buf[idx] = offset_poly;
                for (char c : pair.second) polyphen_strings_buf.push_back(c);
                polyphen_strings_buf.push_back('\0');
                offset_poly += pair.second.size() + 1;
            }
        }

        // --- CSQ Map ---
        std::vector<char> csq_strings_buf;
        std::vector<unsigned int> csq_offsets_buf(256, 0);
        unsigned int offset_csq = 0;
        for (const auto& pair : inv_csq_map) {
            int idx = static_cast<int>(pair.first);
            if (idx >= 0 && idx < 256) {
                csq_offsets_buf[idx] = offset_csq;
                for (char c : pair.second) csq_strings_buf.push_back(c);
                csq_strings_buf.push_back('\0');
                offset_csq += pair.second.size() + 1;
            }
        }

        // --- TSA Map ---
        std::vector<char> tsa_strings_buf;
        std::vector<unsigned int> tsa_offsets_buf(256, 0);
        unsigned int offset_tsa = 0;
        for (const auto& pair : inv_tsa_map) {
            int idx = static_cast<int>(pair.first);
            if (idx >= 0 && idx < 256) {
                tsa_offsets_buf[idx] = offset_tsa;
                for (char c : pair.second) tsa_strings_buf.push_back(c);
                tsa_strings_buf.push_back('\0');
                offset_tsa += pair.second.size() + 1;
            }
        }

        safe_map_alloc((void**)&d_maps.chrom_strings,  (void**)&d_maps.chrom_offsets,  chrom_strings_buf,  chrom_offsets_buf);
        safe_map_alloc((void**)&d_maps.filter_strings, (void**)&d_maps.filter_offsets, filter_strings_buf, filter_offsets_buf);
        safe_map_alloc((void**)&d_maps.polyphen_strings, (void**)&d_maps.polyphen_offsets, polyphen_strings_buf, polyphen_offsets_buf);
        safe_map_alloc((void**)&d_maps.csq_strings,    (void**)&d_maps.csq_offsets,    csq_strings_buf,    csq_offsets_buf);
        safe_map_alloc((void**)&d_maps.tsa_strings,    (void**)&d_maps.tsa_offsets,    tsa_strings_buf,    tsa_offsets_buf);
    }

    // Main execution loop: Process data in chunks
    int df2_start = 0;
    for (int chunk_start = 0; chunk_start < (int)df1.var_number.size(); chunk_start += CHUNK_SIZE){
        int chunk_end = std::min(chunk_start + CHUNK_SIZE, (int)df1.var_number.size());
        int chunk_size = chunk_end - chunk_start;

        while(df2_start < (int)df2.var_id.size() && 
                df2.var_id[df2_start] < df1.var_number[chunk_start]){
            df2_start++;
        }

        float ms_alloc, ms_u2d, ms_kernel, ms_write, ms_free;

        // Allocation phase
        cudaEventRecord(start);
        allocateDevice(df1, df2, df3, chunk_size, chunk_start, chunk_end, df2_start);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&ms_alloc, start, stop);


        host_buffers.chunk_size  = chunk_size;
        host_buffers.chunk_start = chunk_start;
        host_buffers.chunk_end   = chunk_end;
        host_buffers.df2_start   = df2_start;

        // Host preparation phase
        auto t0 = std::chrono::high_resolution_clock::now();
        prepareHostBuffers(df1, df2, df3, df4, host_buffers);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms_prep = std::chrono::duration<double, std::milli>(t1 - t0).count();
           
        // Upload phase
        cudaEventRecord(start);
        uploadToDevice(df1, df2, host_buffers);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&ms_u2d, start, stop);

        int block_size = 32;
        int num_blocks = (chunk_size + block_size - 1) / block_size;

        // Kernel execution phase
        cudaEventRecord(start);
        reconstructKernel<<<num_blocks, block_size>>>(
            d_maps,
            d_df1,
            d_df2,
            d_df3,
            d_output,
            d_line_lens,
            chunk_size
        );
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&ms_kernel, start, stop);

        gpuErrchk(cudaPeekAtLastError());
        gpuErrchk(cudaDeviceSynchronize());

        // Stream compaction and write phase
        cudaEventRecord(start);
        writeChunk(chunk_size, vcf_file);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&ms_write, start, stop);

        // Cleanup phase
        cudaEventRecord(start);
        freeDevice();
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&ms_free, start, stop);

        printf("alloc: %.2f | u2d: %.2f | prep: %.2f | kernel: %.2f | write: %.2f | free: %.2f\n",
               ms_alloc, ms_u2d, ms_prep, ms_kernel, ms_write, ms_free);
    }

    // Final cleanup of global constants
    safe_cuda_free(d_df1.int_names);
    safe_cuda_free(d_df1.float_names);
    safe_cuda_free(d_df1.flag_names);
    safe_cuda_free(d_df2.alt_int_names);
    safe_cuda_free(d_df2.alt_float_names);
    safe_cuda_free(d_maps.chrom_strings);
    safe_cuda_free(d_maps.chrom_offsets);
    safe_cuda_free(d_maps.filter_strings);
    safe_cuda_free(d_maps.filter_offsets);
    safe_cuda_free(d_df3.format_str);
    safe_cuda_free(d_maps.tsa_strings);
    safe_cuda_free(d_maps.tsa_offsets);
    safe_cuda_free(d_maps.polyphen_strings);
    safe_cuda_free(d_maps.polyphen_offsets);
    safe_cuda_free(d_maps.csq_strings);
    safe_cuda_free(d_maps.csq_offsets);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    vcf_file.close();
}

void VCFReconstructorGPU::buildSampleNames(const sample_columns_df& df3) {
    ordered_samp_names.resize(df3.sampNames.size());
    for (const auto& pair : df3.sampNames) {
        ordered_samp_names[pair.second] = pair.first;
    }
}

std::map<std::string, std::string> VCFReconstructorGPU::parseFormatNumbers(const std::string& header_text) {
    std::map<std::string, std::string> result;
    std::istringstream stream(header_text);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.substr(0, 8) != "##FORMAT") continue;
        size_t id_pos  = line.find("ID=");
        size_t num_pos = line.find("Number=");
        if (id_pos == std::string::npos || num_pos == std::string::npos) continue;
        size_t id_start = id_pos + 3;
        size_t id_end   = line.find_first_of(",", id_start);
        std::string id  = line.substr(id_start, id_end - id_start);
        size_t num_start = num_pos + 7;
        size_t num_end   = line.find_first_of(",", num_start);
        std::string number = line.substr(num_start, num_end - num_start);
        result[id] = number;
    }
    return result;
}

/**
 * @param df1         Core variants DataFrame.
 * @param df2         Alternative alleles DataFrame.
 * @param df3         Sample core DataFrame.
 * @param df4         Intersected sample-allele DataFrame.
 * @param chunk_size  Number of variants in the current block.
 * @param chunk_start Starting index in DF1.
 * @param chunk_end   Ending index in DF1.
 * @param df2_start   Starting index in DF2.
 * @param buffer      Output vector that will hold the flattened character data.
 * @param offsets     Output vector that will hold the prefix-summed offsets for the character buffer.
 */
void VCFReconstructorGPU::buildSampleStrings(const var_columns_df& df1,
                                             const alt_columns_df& df2,
                                             const sample_columns_df& df3,
                                             const alt_format_df& df4,
                                             int chunk_size,
                                             int chunk_start,
                                             int chunk_end,
                                             int df2_start,
                                             std::vector<char>& buffer,
                                             std::vector<unsigned int>& offsets) {
    int num_samples = df3.numSample;
    int total_cells = chunk_size * num_samples;
    offsets.resize(total_cells);

    // Sequential pre-calculation per variant
    std::vector<size_t> per_var_num_alts(chunk_size, 0);
    std::vector<long>   per_var_start_df4(chunk_size, -1);
    {
        size_t df2_cursor = (size_t)df2_start;
        size_t df4_cursor = 0;
        while (df4_cursor < df4.var_id.size() &&
               df4.var_id[df4_cursor] < df1.var_number[chunk_start]) df4_cursor++;
               
        for (int i = chunk_start; i < chunk_end; i++) {
            unsigned int current_var_id = df1.var_number[i];
            
            while (df2_cursor < df2.var_id.size() &&
                   df2.var_id[df2_cursor] < current_var_id) df2_cursor++;
                   
            size_t num_alts = 0;
            {
                size_t j = df2_cursor;
                while (j < df2.var_id.size() && df2.var_id[j] == current_var_id) { num_alts++; j++; }
            }
            per_var_num_alts[i - chunk_start] = num_alts;
            
            while (df4_cursor < df4.var_id.size() &&
                   df4.var_id[df4_cursor] < current_var_id) df4_cursor++;
                   
            per_var_start_df4[i - chunk_start] =
                (df4_cursor < df4.var_id.size() &&
                 df4.var_id[df4_cursor] == current_var_id) ? (long)df4_cursor : -1;
        }
    }

    // Parallel Pass 1: Each thread writes to its private scratch buffer
    std::vector<unsigned int> cell_lens(total_cells, 0);
    int num_threads = omp_get_max_threads();
    std::vector<std::vector<char>> scratch(num_threads);
    std::vector<std::vector<unsigned int>> cell_offsets_in_scratch(num_threads);
    std::vector<std::vector<int>> global_indices(num_threads);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int my_var_count = (chunk_size + num_threads - 1) / num_threads;
        scratch[tid].resize((size_t)my_var_count * num_samples * MAX_SAMPLE_STRING_LEN);
        cell_offsets_in_scratch[tid].reserve(my_var_count * num_samples);
        global_indices[tid].reserve(my_var_count * num_samples);

        size_t local_pos = 0;

        #pragma omp for schedule(static)
        for (int i_rel = 0; i_rel < chunk_size; i_rel++) {
            int i = chunk_start + i_rel;
            unsigned int current_var_id = df1.var_number[i];
            size_t num_alts  = per_var_num_alts[i_rel];
            long   start_df4 = per_var_start_df4[i_rel];

            for (int s = 0; s < num_samples; s++) {
                size_t df3_idx = (size_t)i * num_samples + s;

                char* w  = scratch[tid].data() + local_pos;
                char* w0 = w;
                bool first_field = true;
                auto write_sep = [&]() { if (!first_field) *w++ = ':'; first_field = false; };

                // FLAG PER TRONCARE I CAMPI SE GT E' MANCANTE
                bool skip_rest = false;

                // Handle Genotype (GT)
                if (has_gt) {
                    char gt_code = -1;
                    if (gt_in_df3) {
                        gt_code = df3.sample_GT[0].GT[df3_idx];
                    } else if (gt_in_df4 && start_df4 != -1) {
                        long k = start_df4;
                        while (k < (long)df4.var_id.size() && df4.var_id[k] == current_var_id) {
                            if (df4.samp_id[k] == (unsigned short)s) { gt_code = df4.sample_GT.GT[k]; break; }
                            k++;
                        }
                    }
                    write_sep();
                    if (gt_code == -1) { 
                        *w++ = '.'; *w++ = '|'; *w++ = '.'; 
                        skip_rest = true; 
                    }
                    else {
                        const std::string& gt = df3.getGTStringFromChar(gt_code);
                        memcpy(w, gt.data(), gt.size()); w += gt.size();
                        // Protezione extra se il mapping originario conteneva la stringa missing
                        if (gt == ".|." || gt == "./.") {
                            skip_rest = true;
                        }
                    }
                }

                // PROCEDI CON IL RESTO SOLO SE IL GENOTIPO E' VALIDO
                if (!skip_rest) {
                    // Process DF3 Integer fields
                    for (const auto& g : df3_int_groups) {
                        write_sep();
                        size_t group_size = g.col_end - g.col_start;
                        size_t limit = (group_size > 1)
                            ? std::min(group_size, (g.number_kind == 1) ? num_alts : num_alts + 1) : 1;
                        bool any = false, first_v = true;
                        for (size_t cc = g.col_start; cc < g.col_start + limit; cc++) {
                            int v = df3.samp_int[cc].i_int[df3_idx];
                            if (v == -1) continue;
                            if (!first_v) *w++ = ','; first_v = false;
                            w += snprintf(w, 16, "%d", v);
                            any = true;
                        }
                        if (!any) *w++ = '.';
                    }

                    // Process DF3 Float fields
                    for (const auto& g : df3_float_groups) {
                        write_sep();
                        size_t group_size = g.col_end - g.col_start;
                        size_t limit = (group_size > 1)
                            ? std::min(group_size, (g.number_kind == 1) ? num_alts : num_alts + 1) : 1;
                        bool any = false, first_v = true;
                        for (size_t cc = g.col_start; cc < g.col_start + limit; cc++) {
                            float v = (float)df3.samp_float[cc].i_float[df3_idx];
                            if (v == -1.0f) continue;
                            if (!first_v) *w++ = ','; first_v = false;
                            w += snprintf(w, 32, "%g", v);
                            any = true;
                        }
                        if (!any) *w++ = '.';
                    }

                    // Process DF3 String fields
                    for (size_t col = 0; col < df3.samp_string.size(); col++) {
                        write_sep();
                        const std::string& s_val = df3.samp_string[col].i_string[df3_idx];
                        if (s_val.empty()) *w++ = '.';
                        else { memcpy(w, s_val.data(), s_val.size()); w += s_val.size(); }
                    }

                    // Process DF4 Allele-specific fields
                    if (start_df4 != -1) {
                        for (const auto& g : df4_int_groups) {
                            write_sep();
                            bool any = false;
                            long k = start_df4;
                            while (k < (long)df4.var_id.size() && df4.var_id[k] == current_var_id) {
                                if (df4.samp_id[k] == (unsigned short)s) {
                                    for (size_t cc = g.col_start; cc < g.col_end; cc++) {
                                        if (any) *w++ = ',';
                                        int v = df4.samp_int[cc].i_int[k];
                                        if (v == -1) *w++ = '.';
                                        else w += snprintf(w, 16, "%d", v);
                                        any = true;
                                    }
                                }
                                k++;
                            }
                            if (!any) *w++ = '.';
                        }

                        for (const auto& g : df4_float_groups) {
                            write_sep();
                            bool any = false;
                            long k = start_df4;
                            while (k < (long)df4.var_id.size() && df4.var_id[k] == current_var_id) {
                                if (df4.samp_id[k] == (unsigned short)s) {
                                    for (size_t cc = g.col_start; cc < g.col_end; cc++) {
                                        if (any) *w++ = ',';
                                        float v = (float)df4.samp_float[cc].i_float[k];
                                        if (v == -1.0f) *w++ = '.';
                                        else w += snprintf(w, 32, "%g", v);
                                        any = true;
                                    }
                                }
                                k++;
                            }
                            if (!any) *w++ = '.';
                        }

                        for (size_t col = 0; col < df4.samp_string.size(); col++) {
                            write_sep();
                            bool any = false;
                            long k = start_df4;
                            while (k < (long)df4.var_id.size() && df4.var_id[k] == current_var_id) {
                                if (df4.samp_id[k] == (unsigned short)s) {
                                    const std::string& s_val = df4.samp_string[col].i_string[k];
                                    if (any) *w++ = ',';
                                    if (s_val.empty()) *w++ = '.';
                                    else { memcpy(w, s_val.data(), s_val.size()); w += s_val.size(); }
                                    any = true;
                                }
                                k++;
                            }
                            if (!any) *w++ = '.';
                        }
                    }
                } // FINE BLOCCO if (!skip_rest)

                *w++ = '\0';
                size_t cell_len = w - w0;

                int global_idx = i_rel * num_samples + s;
                cell_lens[global_idx] = (unsigned int)cell_len;
                cell_offsets_in_scratch[tid].push_back((unsigned int)local_pos);
                global_indices[tid].push_back(global_idx);
                local_pos += cell_len;
            }
        }
    }

    // Prefix sum to calculate global offsets
    size_t total_size = 0;
    for (int k = 0; k < total_cells; k++) {
        offsets[k] = (unsigned int)total_size;
        total_size += cell_lens[k];
    }

    // Parallel Pass 2: Copy from thread-local scratch buffers to the final contiguous buffer
    buffer.resize(total_size);
    #pragma omp parallel for schedule(static)
    for (int tid = 0; tid < num_threads; tid++) {
        const auto& src_buf  = scratch[tid];
        const auto& src_offs = cell_offsets_in_scratch[tid];
        const auto& gidx     = global_indices[tid];
        for (size_t k = 0; k < gidx.size(); k++) {
            int g = gidx[k];
            unsigned int len = cell_lens[g];
            memcpy(buffer.data() + offsets[g],
                   src_buf.data() + src_offs[k],
                   len);
        }
    }
}

void VCFReconstructorGPU::prepareHostBuffers(const var_columns_df& df1,
                                              const alt_columns_df& df2,
                                              const sample_columns_df& df3,
                                              const alt_format_df& df4,
                                              HostBuffers& buffers) {
    int chunk_size  = buffers.chunk_size;
    int chunk_start = buffers.chunk_start;
    int chunk_end   = buffers.chunk_end;
    int df2_start   = buffers.df2_start;

    // --- Prepare DF1 ID buffer ---
    buffers.id_buffer.clear();
    buffers.id_offsets.resize(chunk_size);
    {
        unsigned int offset = 0;
        for (int i = chunk_start; i < chunk_end; i++) {
            buffers.id_offsets[i - chunk_start] = offset;
            for (char c : df1.id[i]) buffers.id_buffer.push_back(c);
            buffers.id_buffer.push_back('\0');
            offset += df1.id[i].size() + 1;
        }
    }

    // --- Prepare DF1 REF buffer ---
    buffers.ref_buffer.clear();
    buffers.ref_offsets.resize(chunk_size);
    {
        unsigned int offset = 0;
        for (int i = chunk_start; i < chunk_end; i++) {
            buffers.ref_offsets[i - chunk_start] = offset;
            for (char c : df1.ref[i]) buffers.ref_buffer.push_back(c);
            buffers.ref_buffer.push_back('\0');
            offset += df1.ref[i].size() + 1;
        }
    }

    // --- Prepare DF1 INFO integers ---
    // Flatten multi-dimensional array into a 1D contiguous block
    buffers.in_int_buffer.assign(d_df1.num_int_fields * chunk_size, 0);
    for (int f = 0; f < d_df1.num_int_fields; f++) {
        for (int i = chunk_start; i < chunk_end; i++) {
            buffers.in_int_buffer[f * chunk_size + (i - chunk_start)] = df1.in_int[f].i_int[i];
        }
    }

    // --- Prepare DF1 INFO floats ---
    buffers.in_float_buffer.assign(d_df1.num_float_fields * chunk_size, __half(0));
    for (int f = 0; f < d_df1.num_float_fields; f++) {
        for (int i = chunk_start; i < chunk_end; i++) {
            buffers.in_float_buffer[f * chunk_size + (i - chunk_start)] = df1.in_float[f].i_float[i];
        }
    }

    // --- Prepare DF1 INFO flags ---
    buffers.in_flag_buffer.assign(d_df1.num_flag_fields * chunk_size, 0);
    for (int f = 0; f < d_df1.num_flag_fields; f++) {
        for (int i = chunk_start; i < chunk_end; i++) {
            buffers.in_flag_buffer[f * chunk_size + (i - chunk_start)] = df1.in_flag[f].i_flag[i];
        }
    }

    // --- Prepare DF1 INFO Strings ---
    buffers.in_string_buffer.clear();
    buffers.in_string_offsets.assign(d_df1.num_string_fields * chunk_size, 0);
    unsigned int current_offset = 0;

    for (int f = 0; f < d_df1.num_string_fields; f++) {
        for (int i = chunk_start; i < chunk_end; i++) {
            int flat_idx = f * chunk_size + (i - chunk_start);
            buffers.in_string_offsets[flat_idx] = current_offset;
            
            const std::string& val = df1.in_string[f].i_string[i];
            for (char c : val) buffers.in_string_buffer.push_back(c);
            buffers.in_string_buffer.push_back('\0'); // C-style string terminator
            
            current_offset += val.size() + 1;
        }
    }

    // --- Prepare DF2 counts ---
    int df2_count = 0;
    for (size_t j = df2_start; j < df2.var_id.size() && (int)df2.var_id[j] < chunk_end; j++) df2_count++;
    buffers.df2_count = df2_count;

    // --- Prepare DF2 ALT strings ---
    buffers.alt_data_buffer.clear();
    buffers.alt_data_offsets.resize(df2_count);
    {
        unsigned int offset = 0;
        for (int i = df2_start; i < df2_start + df2_count; i++) {
            buffers.alt_data_offsets[i - df2_start] = offset;
            for (char c : df2.alt[i]) buffers.alt_data_buffer.push_back(c);
            buffers.alt_data_buffer.push_back('\0');
            offset += df2.alt[i].size() + 1;
        }
    }

    // --- Setup alt_start and alt_count mapping arrays ---
    buffers.alt_start_buf.assign(chunk_size, 0);
    buffers.alt_count_buf.assign(chunk_size, 0);
    {
        int j = df2_start;
        for (int i = 0; i < chunk_size; i++) {
            unsigned int var_num = df1.var_number[chunk_start + i];
            buffers.alt_start_buf[i] = j - df2_start;
            while (j < df2_start + df2_count && df2.var_id[j] == var_num) j++;
            buffers.alt_count_buf[i] = (j - df2_start) - buffers.alt_start_buf[i];
        }
    }

    // --- Prepare DF2 ALT INFO fields ---
    buffers.alt_int_buffer.assign(d_df2.num_alt_int_fields * df2_count, 0);
    for (int f = 0; f < d_df2.num_alt_int_fields; f++) {
        for (int i = df2_start; i < df2_start + df2_count; i++) {
            buffers.alt_int_buffer[f * df2_count + (i - df2_start)] = df2.alt_int[f].i_int[i];
        }
    }

    buffers.alt_float_buffer.assign(d_df2.num_alt_float_fields * df2_count, __half(0));
    for (int f = 0; f < d_df2.num_alt_float_fields; f++) {
        for (int i = df2_start; i < df2_start + df2_count; i++) {
            buffers.alt_float_buffer[f * df2_count + (i - df2_start)] = df2.alt_float[f].i_float[i];
        }
    }

    // --- Prepare DF2 INFO-ALT Strings ---
    buffers.alt_string_buffer.clear();
    buffers.alt_string_offsets.assign(d_df2.num_alt_string_fields * df2_count, 0);
    current_offset = 0;

    for (int f = 0; f < d_df2.num_alt_string_fields; f++) {
        for (int i = df2_start; i < df2_start + df2_count; i++) {
            int flat_idx = f * df2_count + (i - df2_start);
            buffers.alt_string_offsets[flat_idx] = current_offset;
            
            const std::string& val = df2.alt_string[f].i_string[i];
            for (char c : val) buffers.alt_string_buffer.push_back(c);
            buffers.alt_string_buffer.push_back('\0');
            
            current_offset += val.size() + 1;
        }
    }

    // --- Package complex sample strings (delegates to buildSampleStrings) ---
    buildSampleStrings(df1, df2, df3, df4,
                       chunk_size, chunk_start, chunk_end, df2_start,
                       buffers.sample_buffer, buffers.sample_offsets);
}

/*
* The parser inherently builds forward dictionaries (e.g., String -> char ID) to 
* compress the data in RAM. This method flips those dictionaries (char ID -> String) 
* so that the GPU kernel can retrieve the original textual representation.
*/
void VCFReconstructorGPU::buildInverseMaps(const var_columns_df& df1) {
    for (const auto& pair : df1.chrom_map) {
        inv_chrom_map[pair.second] = pair.first;
    }
    for (const auto& pair : df1.filter_map) {
        inv_filter_map[pair.second] = pair.first;
    }
    for (const auto& pair : polyphenCharMap) {
        inv_polyphen_map[pair.second] = pair.first;
    }
    for (const auto& pair : csqCharMap) {
        inv_csq_map[pair.second] = pair.first;
    }
    for (const auto& pair : tsaCharMap) {
        inv_tsa_map[pair.second] = pair.first;
    }
}