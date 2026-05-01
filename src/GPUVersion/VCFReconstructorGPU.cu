#include "VCFReconstructorGPU.h"
#include "CudaUtils.cuh"
#include "Utils.h"
#include "Maps.h"
#include <sstream>

VCFReconstructorGPU::VCFReconstructorGPU(const std::string& output_vcf_path,
                                         const std::string& header_text)
    : output_vcf_path(output_vcf_path),
      header_text(header_text),
      //DF1
      d_output(nullptr),
      d_line_lens(nullptr)
{}

template <typename T>
inline void safe_cuda_free(T*& d_ptr) {
    if (d_ptr != nullptr) {
        gpuErrchk( cudaFree(d_ptr));
        d_ptr = nullptr;
    }
}

void VCFReconstructorGPU::freeDevice(){
    //DF1
      safe_cuda_free(d_output);
      safe_cuda_free(d_line_lens);
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
      safe_cuda_free(d_df1.int_names);
      safe_cuda_free(d_df1.float_names);
      safe_cuda_free(d_df1.flag_names);
    //DF2
      safe_cuda_free(d_df2.var_id);
      safe_cuda_free(d_df2.alt_data);
      safe_cuda_free(d_df2.alt_offsets);
      safe_cuda_free(d_df2.alt_start);
      safe_cuda_free(d_df2.alt_count);
      safe_cuda_free(d_df2.alt_int);
      safe_cuda_free(d_df2.alt_float);
      safe_cuda_free(d_df2.alt_int_names);
      safe_cuda_free(d_df2.alt_float_names);
    //Inverse maps device
      safe_cuda_free(d_maps.chrom_strings);
      safe_cuda_free(d_maps.chrom_offsets);
      safe_cuda_free(d_maps.filter_strings);
      safe_cuda_free(d_maps.filter_offsets);
    //DF3
        safe_cuda_free(d_df3.sample_data);
        safe_cuda_free(d_df3.sample_offsets);
}

VCFReconstructorGPU::~VCFReconstructorGPU() {
    freeDevice();
}

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

void VCFReconstructorGPU::allocateDevice(const var_columns_df& df1,
                                          const alt_columns_df& df2,
                                          const sample_columns_df& df3,
                                          int chunk_size,
                                          int chunk_start,
                                          int chunk_end,
                                          int df2_start) {
    // Output buffer
    gpuErrchk( cudaMalloc((void**)&d_output,    chunk_size * MAX_LINE_LEN * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_line_lens, chunk_size * sizeof(unsigned int)));

    // DF1 scalar
    gpuErrchk( cudaMalloc((void**)&d_df1.var_number, chunk_size * sizeof(unsigned int)));
    gpuErrchk( cudaMalloc((void**)&d_df1.chrom,      chunk_size * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_df1.pos,        chunk_size * sizeof(unsigned int)));
    gpuErrchk( cudaMalloc((void**)&d_df1.qual,       chunk_size * sizeof(__half)));
    gpuErrchk( cudaMalloc((void**)&d_df1.filter,     chunk_size * sizeof(char)));

    // DF1 id (stringhe variabili)
    size_t id_total_chars = 0;
    for (int i = chunk_start; i < chunk_end; i++)
        id_total_chars += df1.id[i].size() + 1;
    gpuErrchk( cudaMalloc((void**)&d_df1.id_data,    id_total_chars * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_df1.id_offsets, chunk_size * sizeof(unsigned int)));

    // DF1 ref (stringhe variabili)
    size_t ref_total_chars = 0;
    for (int i = chunk_start; i < chunk_end; i++)
        ref_total_chars += df1.ref[i].size() + 1;
    gpuErrchk( cudaMalloc((void**)&d_df1.ref_data,    ref_total_chars * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_df1.ref_offsets, chunk_size * sizeof(unsigned int)));

    // DF1 INFO fields
    d_df1.num_int_fields   = df1.in_int.size();
    d_df1.num_float_fields = df1.in_float.size();
    d_df1.num_flag_fields  = df1.in_flag.size();
    gpuErrchk( cudaMalloc((void**)&d_df1.in_int,      d_df1.num_int_fields   * chunk_size * sizeof(int)));
    gpuErrchk( cudaMalloc((void**)&d_df1.in_float,    d_df1.num_float_fields * chunk_size * sizeof(__half)));
    gpuErrchk( cudaMalloc((void**)&d_df1.in_flag,     d_df1.num_flag_fields  * chunk_size * sizeof(bool)));
    gpuErrchk( cudaMalloc((void**)&d_df1.int_names,   d_df1.num_int_fields   * MAX_NAME_LEN * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_df1.float_names, d_df1.num_float_fields * MAX_NAME_LEN * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_df1.flag_names,  d_df1.num_flag_fields  * MAX_NAME_LEN * sizeof(char)));

    // DF2 - conta entries nel chunk
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
    gpuErrchk( cudaMalloc((void**)&d_df2.alt_int_names,    d_df2.num_alt_int_fields   * MAX_NAME_LEN * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_df2.alt_float_names,  d_df2.num_alt_float_fields * MAX_NAME_LEN * sizeof(char)));

    // Mappe inverse
    size_t chrom_total_chars = 0;
    for (const auto& pair : inv_chrom_map)
        chrom_total_chars += pair.second.size() + 1;
    gpuErrchk( cudaMalloc((void**)&d_maps.chrom_strings, chrom_total_chars * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_maps.chrom_offsets, 256 * sizeof(unsigned int)));

    size_t filter_total_chars = 0;
    for (const auto& pair : inv_filter_map)
        filter_total_chars += pair.second.size() + 1;
    gpuErrchk( cudaMalloc((void**)&d_maps.filter_strings, filter_total_chars * sizeof(char)));
    gpuErrchk( cudaMalloc((void**)&d_maps.filter_offsets, 256 * sizeof(unsigned int)));

    // Stima upper-bound: chunk_size * num_samples * MAX_SAMPLE_STRING_LEN
    // MAX_SAMPLE_STRING_LEN dipende dal tuo dataset, sull'IRBT direi 256 byte abbondanti
    size_t max_sample_chars = (size_t)chunk_size * df3.numSample * MAX_SAMPLE_STRING_LEN;

    gpuErrchk( cudaMalloc(&d_df3.sample_data,    max_sample_chars * sizeof(char)));
    gpuErrchk( cudaMalloc(&d_df3.sample_offsets, chunk_size * df3.numSample * sizeof(unsigned int)));
    d_df3.num_samples = df3.numSample;

}

void VCFReconstructorGPU::hostToDevice(const var_columns_df& df1,
                                          const alt_columns_df& df2,
                                          const sample_columns_df& df3,
                                          const alt_format_df& df4,
                                          int chunk_size,
                                          int chunk_start,
                                          int chunk_end,
                                          int df2_start) {
    // DF1 scalar
    gpuErrchk( cudaMemcpy(d_df1.var_number, df1.var_number.data() + chunk_start, chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));
    gpuErrchk( cudaMemcpy(d_df1.chrom,      df1.chrom.data()      + chunk_start, chunk_size * sizeof(char),     cudaMemcpyHostToDevice));
    gpuErrchk( cudaMemcpy(d_df1.pos,        df1.pos.data()        + chunk_start, chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));
    gpuErrchk( cudaMemcpy(d_df1.qual,       df1.qual.data()       + chunk_start, chunk_size * sizeof(__half),   cudaMemcpyHostToDevice));
    gpuErrchk( cudaMemcpy(d_df1.filter,     df1.filter.data()     + chunk_start, chunk_size * sizeof(char),     cudaMemcpyHostToDevice));

    // Build buffer host for id
    std::vector<char> id_buffer;
    std::vector<unsigned int> id_offsets(chunk_size);
    unsigned int offset = 0;

    for (int i = chunk_start; i < chunk_end; i++) {
        id_offsets[i - chunk_start] = offset;
        // copy the characters of the string into the buffer
        for (char c : df1.id[i]) {
            id_buffer.push_back(c);
        }
        id_buffer.push_back('\0');  // terminator
        offset += df1.id[i].size() + 1;
    }

    // Now copy to GPU
    gpuErrchk( cudaMemcpy(d_df1.id_data,    id_buffer.data(),  id_buffer.size() * sizeof(char),         cudaMemcpyHostToDevice));
    gpuErrchk( cudaMemcpy(d_df1.id_offsets, id_offsets.data(), chunk_size       * sizeof(unsigned int),  cudaMemcpyHostToDevice));

    // Build buffer host for ref
    std::vector<char> ref_buffer;
    std::vector<unsigned int> ref_offsets(chunk_size);
    offset = 0;
    // copy the characters of the string into the buffer
    for (int i = chunk_start; i < chunk_end; i++){
        ref_offsets[i - chunk_start] = offset;

        for (char c: df1.ref[i]){
            ref_buffer.push_back(c);
        }
        ref_buffer.push_back('\0'); // terminator
        offset += df1.ref[i].size() + 1;
    }

    // Now copy to GPU
    gpuErrchk( cudaMemcpy(d_df1.ref_data, ref_buffer.data(), ref_buffer.size() * sizeof(char), cudaMemcpyHostToDevice));
    gpuErrchk( cudaMemcpy(d_df1.ref_offsets, ref_offsets.data(), chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));


    // DF1 INFO fields
    std::vector<int> in_int_buffer(d_df1.num_int_fields * chunk_size);
    for (int f = 0; f < d_df1.num_int_fields; f++) {
        for (int i = chunk_start; i < chunk_end; i++) {
            in_int_buffer[f * chunk_size + (i - chunk_start)] = df1.in_int[f].i_int[i];
        }
    }

    gpuErrchk( cudaMemcpy(d_df1.in_int, in_int_buffer.data(), d_df1.num_int_fields * chunk_size * sizeof(int), cudaMemcpyHostToDevice));

    std::vector<__half> in_float_buffer(d_df1.num_float_fields * chunk_size);
    for (int f = 0; f < d_df1.num_float_fields; f++) {
        for (int i = chunk_start; i < chunk_end; i++) {
            in_float_buffer[f * chunk_size + (i - chunk_start)] = df1.in_float[f].i_float[i];
        }
    }

    gpuErrchk( cudaMemcpy(d_df1.in_float, in_float_buffer.data(), d_df1.num_float_fields * chunk_size * sizeof(__half), cudaMemcpyHostToDevice));

    std::vector<uint8_t> in_flag_buffer(d_df1.num_flag_fields * chunk_size);
    for (int f = 0; f < d_df1.num_flag_fields; f++) {
        for (int i = chunk_start; i < chunk_end; i++) {
            in_flag_buffer[f * chunk_size + (i - chunk_start)] = df1.in_flag[f].i_flag[i];
        }
    }

    gpuErrchk( cudaMemcpy(d_df1.in_flag, in_flag_buffer.data(), d_df1.num_flag_fields * chunk_size * sizeof(uint8_t), cudaMemcpyHostToDevice));

    std::vector<char> int_names_buffer(d_df1.num_int_fields * MAX_NAME_LEN, '\0');
    for (int f = 0; f < d_df1.num_int_fields; f++) {
        const std::string& name = df1.in_int[f].name;
        strncpy(&int_names_buffer[f * MAX_NAME_LEN], name.c_str(), MAX_NAME_LEN - 1);
    }
    gpuErrchk( cudaMemcpy(d_df1.int_names, int_names_buffer.data(), d_df1.num_int_fields * MAX_NAME_LEN * sizeof(char), cudaMemcpyHostToDevice));

    std::vector<char> float_names_buffer(d_df1.num_float_fields * MAX_NAME_LEN, '\0');
    for (int f = 0; f < d_df1.num_float_fields; f++) {
        const std::string& name = df1.in_float[f].name;
        strncpy(&float_names_buffer[f * MAX_NAME_LEN], name.c_str(), MAX_NAME_LEN - 1);
    }
    gpuErrchk( cudaMemcpy(d_df1.float_names, float_names_buffer.data(), d_df1.num_float_fields * MAX_NAME_LEN * sizeof(char), cudaMemcpyHostToDevice));

    std::vector<char> flag_names_buffer(d_df1.num_flag_fields * MAX_NAME_LEN, '\0');
    for (int f = 0; f < d_df1.num_flag_fields; f++) {
        const std::string& name = df1.in_flag[f].name;
        strncpy(&flag_names_buffer[f * MAX_NAME_LEN], name.c_str(), MAX_NAME_LEN - 1);
    }
    gpuErrchk( cudaMemcpy(d_df1.flag_names, flag_names_buffer.data(), d_df1.num_flag_fields * MAX_NAME_LEN * sizeof(char), cudaMemcpyHostToDevice));

    //DF2
    int df2_count = 0;
    for (size_t j = df2_start; j < df2.var_id.size() && (int)df2.var_id[j] < chunk_end; j++){
        df2_count++;
    }
    gpuErrchk( cudaMemcpy(d_df2.var_id, df2.var_id.data() + df2_start, d_df2.num_entries * sizeof(unsigned int), cudaMemcpyHostToDevice));

    // Build buffer host for id
    std::vector<char> alt_data_buffer;
    std::vector<unsigned int> alt_data_offsets(d_df2.num_entries);
    offset = 0;

    for (int i = df2_start; i < df2_start + df2_count; i++) {
        alt_data_offsets[i - df2_start] = offset;
        // copy the characters of the string into the buffer
        for (char c : df2.alt[i]) {
            alt_data_buffer.push_back(c);
        }
        alt_data_buffer.push_back('\0');  // terminator
        offset += df2.alt[i].size() + 1;
    }

    gpuErrchk( cudaMemcpy(d_df2.alt_data, alt_data_buffer.data(), alt_data_buffer.size() * sizeof(char), cudaMemcpyHostToDevice));
    gpuErrchk( cudaMemcpy(d_df2.alt_offsets, alt_data_offsets.data(), d_df2.num_entries * sizeof(unsigned int), cudaMemcpyHostToDevice));

    std::vector<unsigned int> alt_start_buf(chunk_size, 0);
    std::vector<unsigned int> alt_count_buf(chunk_size, 0);

    int j = df2_start;
    for (int i = 0; i < chunk_size; i++) {
        unsigned int var_num = df1.var_number[chunk_start + i];
        alt_start_buf[i] = j - df2_start;
        while (j < df2_start + df2_count && df2.var_id[j] == var_num) {
            j++;
        }
        alt_count_buf[i] = (j - df2_start) - alt_start_buf[i];
    }

    gpuErrchk(cudaMemcpy(d_df2.alt_start, alt_start_buf.data(), chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df2.alt_count, alt_count_buf.data(), chunk_size * sizeof(unsigned int), cudaMemcpyHostToDevice));


    std::vector<int> alt_int_buffer(d_df2.num_alt_int_fields * df2_count);
    for (int f = 0; f < d_df2.num_alt_int_fields; f++) {
        for (int i = df2_start; i < df2_start + d_df2.num_entries; i++) {
            alt_int_buffer[f * d_df2.num_entries + (i - df2_start)] = df2.alt_int[f].i_int[i];
        }
    }

    gpuErrchk( cudaMemcpy(d_df2.alt_int, alt_int_buffer.data(), d_df2.num_alt_int_fields * df2_count * sizeof(int), cudaMemcpyHostToDevice));

    std::vector<__half> alt_float_buffer(d_df2.num_alt_float_fields * df2_count);
    for (int f = 0; f < d_df2.num_alt_float_fields; f++) {
        for (int i = df2_start; i < df2_start + d_df2.num_entries; i++) {
            alt_float_buffer[f * d_df2.num_entries + (i - df2_start)] = df2.alt_float[f].i_float[i];
        }
    }

    gpuErrchk( cudaMemcpy(d_df2.alt_float, alt_float_buffer.data(), d_df2.num_alt_float_fields * df2_count * sizeof(__half), cudaMemcpyHostToDevice));

    std::vector<char> alt_int_names_buffer(d_df2.num_alt_int_fields * MAX_NAME_LEN, '\0');
    for (int f = 0; f < d_df2.num_alt_int_fields; f++) {
        const std::string& name = df2.alt_int[f].name;
        strncpy(&alt_int_names_buffer[f * MAX_NAME_LEN], name.c_str(), MAX_NAME_LEN - 1);
    }
    gpuErrchk( cudaMemcpy(d_df2.alt_int_names, alt_int_names_buffer.data(), d_df2.num_alt_int_fields * MAX_NAME_LEN * sizeof(char), cudaMemcpyHostToDevice));

    std::vector<char> alt_float_names_buffer(d_df2.num_alt_float_fields * MAX_NAME_LEN, '\0');
    for (int f = 0; f < d_df2.num_alt_float_fields; f++) {
        const std::string& name = df2.alt_float[f].name;
        strncpy(&alt_float_names_buffer[f * MAX_NAME_LEN], name.c_str(), MAX_NAME_LEN - 1);
    }
    gpuErrchk( cudaMemcpy(d_df2.alt_float_names, alt_float_names_buffer.data(), d_df2.num_alt_float_fields * MAX_NAME_LEN * sizeof(char), cudaMemcpyHostToDevice));

    // Inverse maps
    std::vector<char> chrom_strings_buffer;
    std::vector<unsigned int> chrom_offsets_buffer(256,0);
    offset = 0;
    for (const auto& pair : inv_chrom_map) {
        int idx = static_cast<int>(pair.first); 
        if (idx >= 0 && idx < 256) {
            chrom_offsets_buffer[idx] = offset; 
            
            for (char c : pair.second) {
                chrom_strings_buffer.push_back(c);
            }
            chrom_strings_buffer.push_back('\0');
            offset += pair.second.size() + 1;
        }
    }
    gpuErrchk( cudaMemcpy(d_maps.chrom_strings, chrom_strings_buffer.data(), chrom_strings_buffer.size() * sizeof(char), cudaMemcpyHostToDevice));
    gpuErrchk( cudaMemcpy(d_maps.chrom_offsets, chrom_offsets_buffer.data(), 256 * sizeof(unsigned int), cudaMemcpyHostToDevice));

    std::vector<char> filter_strings_buffer;
    std::vector<unsigned int> filter_offsets_buffer(256,0);
    offset = 0;
    for (const auto& pair : inv_filter_map) {
        int idx = static_cast<int>(pair.first);
        if (idx >= 0 && idx < 256) {
            filter_offsets_buffer[idx] = offset; 
            
            for (char c : pair.second) {
                filter_strings_buffer.push_back(c);
            }
            filter_strings_buffer.push_back('\0');
            offset += pair.second.size() + 1;
        }
    }
    gpuErrchk( cudaMemcpy(d_maps.filter_strings, filter_strings_buffer.data(), filter_strings_buffer.size() * sizeof(char), cudaMemcpyHostToDevice));
    gpuErrchk( cudaMemcpy(d_maps.filter_offsets, filter_offsets_buffer.data(), 256 * sizeof(unsigned int), cudaMemcpyHostToDevice));

    // Costruzione dei sample pre-formattati
    std::vector<char> sample_buffer;
    std::vector<unsigned int> sample_offsets;
    buildSampleStrings(df1, df2, df3, df4,
                    chunk_size, chunk_start, chunk_end, df2_start,
                    sample_buffer, sample_offsets);

    size_t max_chars = (size_t)chunk_size * df3.numSample * MAX_SAMPLE_STRING_LEN;
    if (sample_buffer.size() > max_chars) {
        throw std::runtime_error("sample_buffer overflow: aumentare MAX_SAMPLE_STRING_LEN");
    }

    gpuErrchk(cudaMemcpy(d_df3.sample_data, sample_buffer.data(), sample_buffer.size() * sizeof(char), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_df3.sample_offsets, sample_offsets.data(), sample_offsets.size() * sizeof(unsigned int), cudaMemcpyHostToDevice));
}

__device__ int device_strcpy(char* dst, const char* src) {
    int i = 0;
    while (src[i] != '\0') {
        dst[i] = src[i];
        i++;
    }
    return i;
}
   
__device__ int device_itoa(int n, char* dst){
    char tmp[20];
    int dst_pos = 0;
    int tmp_len = 0;

    if(n==0){
        dst[0] = '0';
        return 1;
    }

    if(n<0){
        dst[0] = '-';
        n = -n;
        dst_pos = 1;
    }

    while(n>0){
        tmp[tmp_len] = '0' + n % 10;
        n /= 10;
        tmp_len++;
    }

    for (int j = 0; j < tmp_len; j++){
        dst[dst_pos + j] = tmp[tmp_len - j - 1];
    }
    
    return dst_pos + tmp_len;
}

__device__ int device_ftoa(float f, char* dst){

    if(f == -1.0f){
        dst[0] = '.';
        return 1;
    }

    int int_part = (int)f;
    int int_len = device_itoa(int_part, dst);
    dst[int_len] = '.';
    float frac_part = f - int_part;

    if(frac_part < 0){
        frac_part = -frac_part;
    }

    frac_part *= 1000000;
    int frac_len = device_itoa((int)frac_part, dst + int_len + 1);
    return int_len + frac_len + 1;
}   

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
        char* line = output + i * MAX_LINE_LEN;
        int pos = 0;

        //CHROM
        unsigned char chrom_code = df1.chrom[i];
        pos += device_strcpy(line + pos, maps.chrom_strings + maps.chrom_offsets[chrom_code]);
        line[pos++] = '\t';

        //POS
        pos += device_itoa(df1.pos[i], line + pos);
        line[pos++] = '\t';

        // ID
        pos += device_strcpy(line + pos, df1.id_data + df1.id_offsets[i]);
        line[pos++] = '\t';

        //REF
        pos += device_strcpy(line + pos, df1.ref_data + df1.ref_offsets[i]);
        line[pos++] = '\t';

        //ALT
        int alt_begin = df2.alt_start[i];
        int alt_end   = alt_begin + df2.alt_count[i];
        bool first_alt = true;
        for (int j = alt_begin; j < alt_end; j++){
            if (!first_alt) line[pos++] = ',';
            first_alt = false;
            pos += device_strcpy(line + pos, df2.alt_data + df2.alt_offsets[j]);
        }
        line[pos++] = '\t';

        //QUAL
        float q = df1.qual[i];
        pos += device_ftoa(q, line + pos);
        line[pos++] = '\t';

        //FILTER
        unsigned char filter_code = df1.filter[i];
        pos += device_strcpy(line + pos, maps.filter_strings + maps.filter_offsets[filter_code]);
        line[pos++] = '\t';

        //INFO
        bool first_info = true;
        
        //Flag
        for (int f = 0; f < df1.num_flag_fields; f++){
            if (df1.in_flag[f * chunk_size + i]){
                if(!first_info){
                    line[pos++] = ';';
                }
                first_info = false;
                pos += device_strcpy(line + pos, df1.flag_names + f * MAX_NAME_LEN);
            }
        }

        //Int
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

        //Float
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

        //Alt Int (DF2)
        for (int f = 0; f < df2.num_alt_int_fields; f++){
            // Verifica se almeno un valore nel gruppo è diverso da -1
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

        //Alt Float (DF2)
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

        if(first_info){
            line[pos++] = '.';
        }

        //FORMAT + SAMPLES
        if (df3.num_samples > 0) {
            line[pos++] = '\t';
            pos += device_strcpy(line + pos, df3.format_str);

            for (int s = 0; s < df3.num_samples; s++) {
                line[pos++] = '\t';
                int idx = i * df3.num_samples + s;
                pos += device_strcpy(line + pos,
                                    df3.sample_data + df3.sample_offsets[idx]);
            }
        }

        line[pos++] = '\n';
        line_lens[i] = pos;
    }
}

void VCFReconstructorGPU::writeChunk(int num_variants, std::ofstream& out){
    char* h_output = new char[num_variants * MAX_LINE_LEN];
    unsigned int* h_line_lens = new unsigned int[num_variants];

    gpuErrchk( cudaMemcpy(h_output, d_output, num_variants * MAX_LINE_LEN * sizeof(char), cudaMemcpyDeviceToHost));
    gpuErrchk( cudaMemcpy(h_line_lens, d_line_lens, num_variants * sizeof(unsigned int), cudaMemcpyDeviceToHost));

    for (int i = 0; i < num_variants; i++){
        out.write(h_output + i * MAX_LINE_LEN, h_line_lens[i]);
    }

    delete[] h_output;
    delete[] h_line_lens;
}

void VCFReconstructorGPU::run(const var_columns_df& df1,
         const alt_columns_df& df2,
         const sample_columns_df& df3,
         const alt_format_df& df4){

    buildInverseMaps(df1);
    buildSampleNames(df3);
    format_numbers = parseFormatNumbers(header_text);

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

    //Costruzione Format string
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


    int df2_start = 0;
    for (int chunk_start = 0; chunk_start < (int)df1.var_number.size(); chunk_start += CHUNK_SIZE){
        int chunk_end = std::min(chunk_start + CHUNK_SIZE, (int)df1.var_number.size());
        int chunk_size = chunk_end - chunk_start;

        while(df2_start < (int)df2.var_id.size() && 
                df2.var_id[df2_start] < df1.var_number[chunk_start]){
            df2_start++;
        }

        float ms_alloc, ms_h2d, ms_kernel, ms_write, ms_free;

        cudaEventRecord(start);
        allocateDevice(df1, df2, df3, chunk_size, chunk_start, chunk_end, df2_start);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&ms_alloc, start, stop);

        cudaEventRecord(start);
        hostToDevice(df1, df2, df3, df4, chunk_size, chunk_start, chunk_end, df2_start);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&ms_h2d, start, stop);

        int block_size = 32;
        int num_blocks = (chunk_size + block_size - 1) / block_size;

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

        cudaEventRecord(start);
        writeChunk(chunk_size, vcf_file);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&ms_write, start, stop);

        cudaEventRecord(start);
        freeDevice();
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&ms_free, start, stop);

        printf("alloc: %.2f | h2d: %.2f | kernel: %.2f | write: %.2f | free: %.2f\n",
               ms_alloc, ms_h2d, ms_kernel, ms_write, ms_free);
    }

    safe_cuda_free(d_df3.format_str);

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
    offsets.resize(chunk_size * num_samples);
    buffer.clear();
    buffer.reserve((size_t)chunk_size * num_samples * 64); // hint iniziale

    // Cursore DF4 sequenziale sul chunk
    size_t df4_cursor = 0;
    while (df4_cursor < df4.var_id.size() &&
           df4.var_id[df4_cursor] < df1.var_number[chunk_start]) {
        df4_cursor++;
    }
    // Nota: invece di una variabile membro, ti conviene passare il var_number dell'inizio
    // chunk come parametro, oppure ricalcolare dal df1 (vedi sotto).

    for (int i = chunk_start; i < chunk_end; i++) {
        unsigned int current_var_id = df1.var_number[i]; // <-- vedi sotto

        size_t df2_cursor = (size_t)df2_start;  // o passalo in modo pulito
        while (df2_cursor < df2.var_id.size() && df2.var_id[df2_cursor] < current_var_id) {
            df2_cursor++;
        }
        size_t num_alts = 0;
        {
            size_t j = df2_cursor;
            while (j < df2.var_id.size() && df2.var_id[j] == current_var_id) {
                num_alts++;
                j++;
            }
        }

        // Avanza df4_cursor fino al primo record con var_id == current_var_id
        while (df4_cursor < df4.var_id.size() &&
               df4.var_id[df4_cursor] < current_var_id) {
            df4_cursor++;
        }
        int start_df4 = (df4_cursor < df4.var_id.size() &&
                         df4.var_id[df4_cursor] == current_var_id)
                        ? (int)df4_cursor : -1;

        for (int s = 0; s < num_samples; s++) {
            size_t df3_idx = (size_t)i * num_samples + s;
            std::string sample_data;

            // ---------- GT ----------
            if (has_gt) {
                char gt_code = -1;
                if (gt_in_df3) {
                    gt_code = df3.sample_GT[0].GT[df3_idx];
                } else if (gt_in_df4 && start_df4 != -1) {
                    size_t k = (size_t)start_df4;
                    while (k < df4.var_id.size() && df4.var_id[k] == current_var_id) {
                        if (df4.samp_id[k] == (unsigned short)s) {
                            gt_code = df4.sample_GT.GT[k];
                            break;
                        }
                        k++;
                    }
                }
                if (gt_code == -1) {
                    sample_data = ".|.";
                } else {
                    sample_data = df3.getGTStringFromChar(gt_code);
                }
            }

            // ---------- DF3 int (raggruppati) ----------
            {
                size_t col = 0;
                while (col < df3.samp_int.size()) {
                    std::string base = stripTrailingDigits(df3.samp_int[col].name);
                    if (!sample_data.empty()) sample_data += ":";

                    size_t group_end = col;
                    while (group_end < df3.samp_int.size() &&
                           stripTrailingDigits(df3.samp_int[group_end].name) == base) {
                        group_end++;
                    }
                    size_t group_size = group_end - col;

                    size_t limit;
                    if (group_size > 1) {
                        std::string num = format_numbers.count(base) ? format_numbers[base] : "R";
                        size_t real = (num == "A") ? num_alts : num_alts + 1;
                        limit = std::min(group_size, real);
                    } else {
                        limit = 1;
                    }

                    std::string group_val;
                    size_t printed = 0;
                    while (col < df3.samp_int.size() &&
                           stripTrailingDigits(df3.samp_int[col].name) == base) {
                        if (printed < limit) {
                            int v = df3.samp_int[col].i_int[df3_idx];
                            if (v != -1) {
                                if (!group_val.empty()) group_val += ",";
                                group_val += std::to_string(v);
                            }
                            printed++;
                        }
                        col++;
                    }
                    sample_data += group_val.empty() ? "." : group_val;
                }
            }

            // ---------- DF3 float (raggruppati) ----------
            {
                size_t col = 0;
                while (col < df3.samp_float.size()) {
                    std::string base = stripTrailingDigits(df3.samp_float[col].name);
                    if (!sample_data.empty()) sample_data += ":";

                    size_t group_end = col;
                    while (group_end < df3.samp_float.size() &&
                           stripTrailingDigits(df3.samp_float[group_end].name) == base) {
                        group_end++;
                    }
                    size_t group_size = group_end - col;

                    size_t limit;
                    if (group_size > 1) {
                        std::string num = format_numbers.count(base) ? format_numbers[base] : "R";
                        size_t real = (num == "A") ? num_alts : num_alts + 1;
                        limit = std::min(group_size, real);
                    } else {
                        limit = 1;
                    }

                    std::string group_val;
                    size_t printed = 0;
                    while (col < df3.samp_float.size() &&
                           stripTrailingDigits(df3.samp_float[col].name) == base) {
                        if (printed < limit) {
                            float v = (float)df3.samp_float[col].i_float[df3_idx];
                            if (v != -1.0f) {
                                if (!group_val.empty()) group_val += ",";
                                group_val += std::to_string(v);
                            }
                            printed++;
                        }
                        col++;
                    }
                    sample_data += group_val.empty() ? "." : group_val;
                }
            }

            // ---------- DF3 string ----------
            for (size_t col = 0; col < df3.samp_string.size(); col++) {
                if (!sample_data.empty()) sample_data += ":";
                sample_data += df3.samp_string[col].i_string[df3_idx];
            }

            // ---------- DF4 (raggruppati int/float, poi string) ----------
            if (start_df4 != -1) {
                // int
                size_t col = 0;
                while (col < df4.samp_int.size()) {
                    std::string base = stripTrailingDigits(df4.samp_int[col].name);
                    std::string field_val;

                    size_t group_start = col;
                    while (col < df4.samp_int.size() &&
                           stripTrailingDigits(df4.samp_int[col].name) == base) {
                        col++;
                    }
                    size_t k = (size_t)start_df4;
                    while (k < df4.var_id.size() && df4.var_id[k] == current_var_id) {
                        if (df4.samp_id[k] == (unsigned short)s) {
                            for (size_t c = group_start; c < col; c++) {
                                if (!field_val.empty()) field_val += ",";
                                int v = df4.samp_int[c].i_int[k];
                                field_val += (v != -1) ? std::to_string(v) : ".";
                            }
                        }
                        k++;
                    }
                    sample_data += ":";
                    sample_data += field_val.empty() ? "." : field_val;
                }

                // float
                col = 0;
                while (col < df4.samp_float.size()) {
                    std::string base = stripTrailingDigits(df4.samp_float[col].name);
                    std::string field_val;

                    size_t group_start = col;
                    while (col < df4.samp_float.size() &&
                           stripTrailingDigits(df4.samp_float[col].name) == base) {
                        col++;
                    }
                    size_t k = (size_t)start_df4;
                    while (k < df4.var_id.size() && df4.var_id[k] == current_var_id) {
                        if (df4.samp_id[k] == (unsigned short)s) {
                            for (size_t c = group_start; c < col; c++) {
                                if (!field_val.empty()) field_val += ",";
                                float v = (float)df4.samp_float[c].i_float[k];
                                field_val += (v != -1.0f) ? std::to_string(v) : ".";
                            }
                        }
                        k++;
                    }
                    sample_data += ":";
                    sample_data += field_val.empty() ? "." : field_val;
                }

                // string
                for (size_t cc = 0; cc < df4.samp_string.size(); cc++) {
                    std::string field_val;
                    size_t k = (size_t)start_df4;
                    while (k < df4.var_id.size() && df4.var_id[k] == current_var_id) {
                        if (df4.samp_id[k] == (unsigned short)s) {
                            if (!field_val.empty()) field_val += ",";
                            field_val += df4.samp_string[cc].i_string[k];
                        }
                        k++;
                    }
                    sample_data += ":";
                    sample_data += field_val.empty() ? "." : field_val;
                }
            }

            // ---------- Pulizia finali vuoti (VCF Spec) ----------
            // Rimuove tutti i ":." finali finché ce ne sono
            while (sample_data.size() >= 2 && sample_data.substr(sample_data.size() - 2) == ":.") {
                sample_data.erase(sample_data.size() - 2);
            }
            // Rimuove un eventuale ":" finale rimasto appeso
            if (!sample_data.empty() && sample_data.back() == ':') {
                sample_data.pop_back();
            }

            // ---------- Append nel buffer ----------
            size_t out_idx = (size_t)(i - chunk_start) * num_samples + s;
            offsets[out_idx] = (unsigned int)buffer.size();
            buffer.insert(buffer.end(), sample_data.begin(), sample_data.end());
            buffer.push_back('\0');
        }
    }
}