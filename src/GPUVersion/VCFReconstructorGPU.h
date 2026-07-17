#ifndef VCF_RECONSTRUCTOR_GPU_H
#define VCF_RECONSTRUCTOR_GPU_H

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include "VCFDataFrames.h"
#include "GPUStructs.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>

#define MAX_LINE_LEN 756
#define MAX_NAME_LEN 16
#define MAX_SAMPLE_STRING_LEN 64
#define CHUNK_SIZE   10000
#define BLOCK_SIZE   256
#define NUM_BLOCKS   256

/**
 * @brief Staging structure residing in system RAM (Host).
 * 
 * Buffers and prepares data in chunks before sending them to the CUDA device.
 * Optimizes writes and avoids bottlenecks in PCIe transfers.
 */
struct HostBuffers {
    // DF1
    std::vector<char>          id_buffer;       ///< Temporary host buffer for variant IDs
    std::vector<unsigned int>  id_offsets;      ///< Temporary offsets for variant IDs
    std::vector<char>          ref_buffer;      ///< Temporary host buffer for REF
    std::vector<unsigned int>  ref_offsets;     ///< Temporary offsets for REF
    std::vector<int>           in_int_buffer;   ///< Temporary buffer for integer INFO
    std::vector<__half>        in_float_buffer; ///< Temporary buffer for float INFO
    std::vector<uint8_t>       in_flag_buffer;  ///< Temporary buffer for flag INFO
    std::vector<char>         in_string_buffer;   ///< Temporary flat host buffer for INFO strings.
    std::vector<unsigned int> in_string_offsets;  ///< Temporary offsets for INFO strings.

    // DF2
    int                        df2_count = 0;      ///< Number of alleles in the current chunk
    std::vector<char>          alt_data_buffer;    ///< Temporary host buffer for ALT
    std::vector<unsigned int>  alt_data_offsets;   ///< Offsets for ALT
    std::vector<unsigned int>  alt_start_buf;      ///< Start mapping variants->alleles
    std::vector<unsigned int>  alt_count_buf;      ///< Count mapping variants->alleles
    std::vector<int>           alt_int_buffer;     ///< Temporary buffer for INFO-ALT int
    std::vector<__half>        alt_float_buffer;   ///< Temporary buffer for INFO-ALT float
    std::vector<char>         alt_string_buffer;  ///< Temporary flat host buffer for INFO-ALT strings.
    std::vector<unsigned int> alt_string_offsets; ///< Temporary offsets for INFO-ALT strings.

    // DF3
    std::vector<char>          sample_buffer;      ///< Pre-calculated patient data buffer
    std::vector<unsigned int>  sample_offsets;     ///< Offsets for patient strings

    // Chunk Metadata
    int chunk_size  = 0;  ///< Size of the variant chunk
    int chunk_start = 0;  ///< Starting index in DF1
    int chunk_end   = 0;  ///< Ending index in DF1
    int df2_start   = 0;  ///< Synchronized starting index in DF2
};

struct TimingResult {
    double setup_ms  = 0;
    double alloc_ms  = 0;
    double u2d_ms    = 0;
    double prep_ms   = 0;
    double kernel_ms = 0;
    double write_ms  = 0;
    double free_ms   = 0;
    double drain_ms  = 0;
};

/**
 * @brief Manages the VCF reconstruction pipeline leveraging GPU acceleration.
 *
 * This class loads the DataFrames extracted from CSVs, splits them into chunks,
 * and instructs the CUDA device to build the VCF file strings in a highly parallel manner.
 */
class VCFReconstructorGPU {
public:

    /**
     * @brief Defines the grouping logic for FORMAT fields (e.g., AD, DP).
     */
    struct GroupInfo {
        size_t col_start;   ///< Starting column
        size_t col_end;     ///< Ending column
        char number_kind;   ///< Cardinality rule: 0='R' (Ref+Alt), 1='A' (Alt), 2=Fixed
    };

    /**
     * @brief Constructor for the GPU pipeline.
     * 
     * @param output_vcf_path Path to the final .vcf file on disk.
     * @param header_text     Full header text to write at the top of the file.
     */
    VCFReconstructorGPU(const std::string& output_vcf_path, const std::string& header_text);

    /**
    * @brief Class destructor. 
    * 
    * Ensures all allocated memory on both the CUDA device and the host is properly 
    * released when the VCFReconstructorGPU instance goes out of scope.
    */
    ~VCFReconstructorGPU();

    /**
     * @brief Main entry point for the reconstruction process.
     * 
     * Triggers map creation, allocation, chunk processing (GPU), 
     * and the final write to disk for each block.
     * 
     * @param df1 Base DataFrame (Core variants)
     * @param df2 DataFrame for multi-allelic variants (ALT)
     * @param df3 Base DataFrame for samples
     * @param df4 Intersected DataFrame for samples x alleles
     */
    TimingResult run(const var_columns_df& df1,
                     const alt_columns_df& df2,
                     const sample_columns_df& df3,
                     const alt_format_df& df4);

private:
    std::string output_vcf_path; ///< Path to the generated output VCF file.
    std::string header_text;     ///< The full VCF header string to be written at the top of the file.

    // Inverse maps (host side)
    std::map<char, std::string> inv_chrom_map;    ///< Maps chromosome char codes back to original strings (e.g., '1' -> "chr1").
    std::map<char, std::string> inv_filter_map;   ///< Maps filter char codes back to original strings (e.g., '0' -> "PASS").
    std::map<int,  std::string> inv_tsa_map;      ///< Maps TSA integer codes to original strings.
    std::map<char, std::string> inv_polyphen_map; ///< Maps PolyPhen char codes to original strings.
    std::map<char, std::string> inv_csq_map;      ///< Maps CSQ char codes to original strings.

    // Device output buffer
    char*         d_output;       ///< Sparse device buffer for VCF text. Size: [CHUNK_SIZE * MAX_LINE_LEN].
    unsigned int* d_line_lens;    ///< Device array storing the exact string length of each generated row.

    DeviceVarColumns d_df1;       ///< Shadow structure holding device pointers for DF1 (core variants).
    DeviceAltColumns d_df2;       ///< Shadow structure holding device pointers for DF2 (alternative alleles).
    DeviceSampleColumns d_df3;    ///< Shadow structure holding device pointers for DF3/DF4 (samples).
    DeviceMaps d_maps;            ///< Shadow structure holding device pointers for encoding dictionaries.

    // State and Format tracking
    bool has_gt;                  ///< Flag indicating if Genotype (GT) data exists in the dataset.
    bool gt_in_df3;               ///< Flag indicating if Genotype (GT) data is stored in DF3.
    bool gt_in_df4;               ///< Flag indicating if Genotype (GT) data is stored in DF4.
    std::string format_str;       ///< The dynamically built FORMAT string for the current VCF (e.g., "GT:AD:DP:GQ").
    std::vector<std::string> ordered_samp_names; ///< Sample names ordered by their internal integer ID.

    /**
     * @brief Parses the VCF header to extract the 'Number' attribute for each FORMAT ID.
     * 
     * @param header_text The raw VCF header string.
     * @return std::map<std::string, std::string> Map associating FORMAT IDs to their 'Number' value (e.g., "AD" -> "R").
     */
    std::map<std::string, std::string> parseFormatNumbers(const std::string& header_text);
    std::map<std::string, std::string> format_numbers; ///< Cached map of FORMAT ID cardinalities.

    // Host staging and grouping
    HostBuffers host_buffers;     ///< Struct managing staging buffers on the host before transferring to the device.

    std::vector<GroupInfo> df3_int_groups;   ///< Grouping logic for DF3 integer fields (handles multi-value keys).
    std::vector<GroupInfo> df3_float_groups; ///< Grouping logic for DF3 float fields.
    std::vector<GroupInfo> df4_int_groups;   ///< Grouping logic for DF4 integer fields.
    std::vector<GroupInfo> df4_float_groups; ///< Grouping logic for DF4 float fields.

    // Compaction and Disk I/O variables
    char* d_compacted;               ///< Device buffer for the compacted, contiguous VCF text.
    unsigned int* d_output_offsets;  ///< Device array storing prefix-summed offsets for stream compaction.

    //size_t h_compacted_capacity = 0; ///< Current capacity of the host-side compacted buffer.
    //char* h_compacted = nullptr;     ///< Host buffer receiving the compacted chunk from the device before disk write.   

    // Writer thread infrastructure
    /**
     * @brief Descriptor for a pending disk-write job dispatched to the writer thread.
     */
    struct WriteJob {
        int    buffer_idx;  ///< Index into h_compacted_pool identifying the source buffer.
        size_t total_bytes; ///< Number of bytes to write from the chosen buffer.
    };

    static constexpr int NUM_WRITE_BUFFERS = 2;
    char*  h_compacted_pool[NUM_WRITE_BUFFERS]          = {nullptr, nullptr}; ///< Double-buffered host memory pool for compacted VCF text.
    size_t h_compacted_pool_capacity[NUM_WRITE_BUFFERS] = {0, 0};             ///< Current allocated capacity of each pool slot in bytes.

    std::queue<WriteJob>    write_queue;        ///< Queue of pending write jobs consumed by the writer thread.
    std::mutex              write_mutex;        ///< Protects write_queue, write_buffer_busy, and writer_should_stop.
    std::condition_variable cv_job_available;  ///< Notified when a new job is pushed onto write_queue.
    std::condition_variable cv_buffer_free;    ///< Notified when a buffer slot becomes available for reuse.
    bool write_buffer_busy[NUM_WRITE_BUFFERS] = {false, false}; ///< Tracks which pool slots are currently held by the writer thread.
    bool writer_should_stop = false;           ///< Set to true to signal the writer thread to exit after draining the queue.

    std::thread   writer_thread; ///< Background thread responsible for writing compacted chunks to disk.
    std::ofstream* writer_out = nullptr; ///< Non-owning pointer to the output file stream managed by run().

    /**
     * @brief Background thread entry point; drains the write queue and flushes chunks to disk.
     */
    void writerLoop();

    /**
     * @brief Reconstructs the inverse maps from the DF1 dictionaries for host-side lookups.
     * 
     * @param df1 The parsed var_columns_df containing the forward dictionaries.
     */
    void buildInverseMaps(const var_columns_df& df1);

    /**
     * @brief Allocates all required memory buffers on the CUDA device for the current chunk.
     * 
     * @param df1 Core variants DataFrame.
     * @param df2 Alternative alleles DataFrame.
     * @param df3 Sample core DataFrame.
     * @param chunk_size Number of variants to process in the current execution block.
     * @param chunk_start Starting index in DF1.
     * @param chunk_end Ending index in DF1.
     * @param df2_start Starting index in DF2 corresponding to chunk_start.
     */
    void allocateDevice(const var_columns_df& df1,
                        const alt_columns_df& df2,
                        const sample_columns_df& df3,
                        int chunk_size,
                        int chunk_start,
                        int chunk_end,
                        int df2_start);

    /**
     * @brief Safely deallocates all CUDA device memory to prevent memory leaks.
     */
    void freeDevice();

    /**
     * @brief Populates the HostBuffers struct with flattened, contiguous data ready for optimal GPU transfer.
     * 
     * @param df1 Core variants DataFrame.
     * @param df2 Alternative alleles DataFrame.
     * @param df3 Sample core DataFrame.
     * @param df4 Intersected sample-allele DataFrame.
     * @param buffers The HostBuffers struct to be populated.
     */
    void prepareHostBuffers(const var_columns_df& df1,
                            const alt_columns_df& df2,
                            const sample_columns_df& df3,
                            const alt_format_df& df4,
                            HostBuffers& buffers);

    /**
     * @brief Transfers staging buffers (HostBuffers) from system RAM to GPU VRAM using cudaMemcpy.
     * 
     * @param df1 Core variants DataFrame.
     * @param df2 Alternative alleles DataFrame.
     * @param buffers The fully populated HostBuffers struct.
     */
    void uploadToDevice(const var_columns_df& df1,
                        const alt_columns_df& df2,
                        const HostBuffers& buffers);

    /**
     * @brief Performs stream compaction using Prefix Sum (CUB) and writes the processed, contiguous chunk to disk.
     * 
     * @param num_variants The number of variants processed in the current chunk.
     */
    void writeChunk(int num_variants);

    /**
     * @brief Populates the ordered_samp_names vector by mapping sample IDs to their string names.
     * 
     * @param df3 Sample core DataFrame containing the sample name dictionary.
     */
    void buildSampleNames(const sample_columns_df& df3);

    /**
     * @brief Pre-processes and flattens complex, multi-allelic sample data into a contiguous character buffer using OpenMP.
     * 
     * @param df1 Core variants DataFrame.
     * @param df2 Alternative alleles DataFrame.
     * @param df3 Sample core DataFrame.
     * @param df4 Intersected sample-allele DataFrame.
     * @param chunk_size Number of variants in the current block.
     * @param chunk_start Starting index in DF1.
     * @param chunk_end Ending index in DF1.
     * @param df2_start Starting index in DF2.
     * @param buffer Output vector that will hold the flattened character data.
     * @param offsets Output vector that will hold the prefix-summed offsets for the character buffer.
     */
    void buildSampleStrings(const var_columns_df& df1,
                        const alt_columns_df& df2,
                        const sample_columns_df& df3,
                        const alt_format_df& df4,
                        int chunk_size,
                        int chunk_start,
                        int chunk_end,
                        int df2_start,
                        std::vector<char>& buffer,
                        std::vector<unsigned int>& offsets);
    
};

#endif