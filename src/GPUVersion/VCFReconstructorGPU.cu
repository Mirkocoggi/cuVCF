#include "VCFReconstructorGPU.h"
#include "CudaUtils.cuh"
#include "Utils.h"
#include "Maps.h"

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
      safe_cuda_free(d_df2.alt_int);
      safe_cuda_free(d_df2.alt_float);
      safe_cuda_free(d_df2.alt_int_names);
      safe_cuda_free(d_df2.alt_float_names);
    //Inverse maps device
      safe_cuda_free(d_maps.chrom_strings);
      safe_cuda_free(d_maps.chrom_offsets);
      safe_cuda_free(d_maps.filter_strings);
      safe_cuda_free(d_maps.filter_offsets);
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
}

void VCFReconstructorGPU::hostToDevice(const var_columns_df& df1,
                                          const alt_columns_df& df2,
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
    char* output, unsigned int* line_lens,
    int chunk_size
){
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= chunk_size) return;

    char* line = output + tid * MAX_LINE_LEN;
    int pos = 0;

    //CHROM
    unsigned char chrom_code = df1.chrom[tid];
    pos += device_strcpy(line + pos, maps.chrom_strings + maps.chrom_offsets[chrom_code]);
    line[pos++] = '\t';

    //POS
    pos += device_itoa(df1.pos[tid], line + pos);
    line[pos++] = '\t';

    // ID
    pos += device_strcpy(line + pos, df1.id_data + df1.id_offsets[tid]);
    line[pos++] = '\t';

    //REF
    pos += device_strcpy(line + pos, df1.ref_data + df1.ref_offsets[tid]);
    line[pos++] = '\t';

    //ALT
    bool first_alt = true;
    for (int j = 0; j < df2.num_entries; j++){
        if (df2.var_id[j] == df1.var_number[tid]){
            if (!first_alt){
                line[pos++] = ',';
            }
            first_alt = false;
            pos += device_strcpy(line + pos, df2.alt_data + df2.alt_offsets[j]);
        }
    }
    line[pos++] = '\t';

    //QUAL
    float q = df1.qual[tid];
    pos += device_ftoa(q, line + pos);
    line[pos++] = '\t';

    //FILTER
    unsigned char filter_code = df1.filter[tid];
    pos += device_strcpy(line + pos, maps.filter_strings + maps.filter_offsets[filter_code]);
    line[pos++] = '\t';

    //INFO
    bool first_info = true;
    
    //Flag
    for (int f = 0; f < df1.num_flag_fields; f++){
        if (df1.in_flag[f * chunk_size + tid]){
            if(!first_info){
                line[pos++] = ';';
            }
            first_info = false;
            pos += device_strcpy(line + pos, df1.flag_names + f * MAX_NAME_LEN);
        }
    }

    //Int
    for (int f = 0; f < df1.num_int_fields; f++){
        if (df1.in_int[f * chunk_size + tid] != -1){
            if(!first_info){
                line[pos++] = ';';
            }
            first_info = false;
            pos += device_strcpy(line + pos, df1.int_names + f * MAX_NAME_LEN);
            line[pos++] = '=';
            pos += device_itoa(df1.in_int[f * chunk_size + tid], line + pos);
        }
    }

    //Float
    for (int f = 0; f < df1.num_float_fields; f++){
        if (__half2float(df1.in_float[f * chunk_size + tid]) != -1.0f){
            if(!first_info){
                line[pos++] = ';';
            }
            first_info = false;
            pos += device_strcpy(line + pos, df1.float_names + f * MAX_NAME_LEN);
            line[pos++] = '=';
            pos += device_ftoa(__half2float(df1.in_float[f * chunk_size + tid]), line + pos);
        }
    }

    if(first_info){
        line[pos++] = '.';
    }
    line[pos++] = '\n';
    line_lens[tid] = pos;
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

void VCFReconstructorGPU::run(const var_columns_df& df1, const alt_columns_df& df2){
    buildInverseMaps(df1);

    std::ofstream vcf_file(output_vcf_path);
    if (!vcf_file.is_open()) {
        throw std::runtime_error("Error opening output VCF file: " + output_vcf_path);
    }

    vcf_file << header_text;
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    int df2_start = 0;
    for (int chunk_start = 0; chunk_start < (int)df1.var_number.size(); chunk_start += CHUNK_SIZE){
        int chunk_end = std::min(chunk_start + CHUNK_SIZE, (int)df1.var_number.size());
        int chunk_size = chunk_end - chunk_start;

        while(df2_start < (int)df2.var_id.size() && 
                df2.var_id[df2_start] < df1.var_number[chunk_start]){
            df2_start++;
        }

        allocateDevice(df1, df2, chunk_size, chunk_start, chunk_end, df2_start);
        hostToDevice(df1, df2, chunk_size, chunk_start, chunk_end, df2_start);

        int block_size = 256;
        int num_blocks = (chunk_size + block_size -1) / block_size;
        
        int df2_count = 0;
        for(size_t j = df2_start; j < df2.var_id.size() && (int)df2.var_id[j] < chunk_end; j++){
            df2_count++;
        }

        reconstructKernel<<<num_blocks, block_size>>>(
            d_maps,
            d_df1,
            d_df2,
            d_output,
            d_line_lens,
            chunk_size
        );

        gpuErrchk( cudaPeekAtLastError() );

        gpuErrchk( cudaDeviceSynchronize() );

        writeChunk(chunk_size, vcf_file);
        freeDevice();
    }

    vcf_file.close();
}