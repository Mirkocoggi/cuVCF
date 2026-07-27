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
#include "PinnedAllocator.h"
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

// Vectors backed by pinned (page-locked) host memory: every one of these is
// used as a cudaMemcpy source in uploadToDevice(), so pinning them is what
// will let those transfers become truly asynchronous (cudaMemcpyAsync) later.
using PinnedCharVec = std::vector<char, PinnedAllocator<char>>;
using PinnedUintVec = std::vector<unsigned int, PinnedAllocator<unsigned int>>;
using PinnedIntVec  = std::vector<int, PinnedAllocator<int>>;
using PinnedHalfVec = std::vector<__half, PinnedAllocator<__half>>;
using PinnedU8Vec   = std::vector<uint8_t, PinnedAllocator<uint8_t>>;

/**
 * @brief Staging structure residing in system RAM (Host).
 *
 * Buffers and prepares data in chunks before sending them to the CUDA device.
 * Optimizes writes and avoids bottlenecks in PCIe transfers.
 */
struct HostBuffers {
    // DF1
    PinnedCharVec id_buffer;       ///< Temporary host buffer for variant IDs
    PinnedUintVec id_offsets;      ///< Temporary offsets for variant IDs
    PinnedCharVec ref_buffer;      ///< Temporary host buffer for REF
    PinnedUintVec ref_offsets;     ///< Temporary offsets for REF
    PinnedIntVec  in_int_buffer;   ///< Temporary buffer for integer INFO
    PinnedHalfVec in_float_buffer; ///< Temporary buffer for float INFO
    PinnedU8Vec   in_flag_buffer;  ///< Temporary buffer for flag INFO
    PinnedCharVec in_string_buffer;   ///< Temporary flat host buffer for INFO strings.
    PinnedUintVec in_string_offsets;  ///< Temporary offsets for INFO strings.

    // DF2
    int           df2_count = 0;      ///< Number of alleles in the current chunk
    PinnedCharVec alt_data_buffer;    ///< Temporary host buffer for ALT
    PinnedUintVec alt_data_offsets;   ///< Offsets for ALT
    PinnedUintVec alt_start_buf;      ///< Start mapping variants->alleles
    PinnedUintVec alt_count_buf;      ///< Count mapping variants->alleles
    PinnedIntVec  alt_int_buffer;     ///< Temporary buffer for INFO-ALT int
    PinnedHalfVec alt_float_buffer;   ///< Temporary buffer for INFO-ALT float
    PinnedCharVec alt_string_buffer;  ///< Temporary flat host buffer for INFO-ALT strings.
    PinnedUintVec alt_string_offsets; ///< Temporary offsets for INFO-ALT strings.

    // DF3
    PinnedCharVec sample_buffer;      ///< Pre-calculated patient data buffer
    PinnedUintVec sample_offsets;     ///< Offsets for patient strings

    // Chunk Metadata
    int chunk_size  = 0;  ///< Size of the variant chunk
    int chunk_start = 0;  ///< Starting index in DF1
    int chunk_end   = 0;  ///< Ending index in DF1
    int df2_start   = 0;  ///< Synchronized starting index in DF2
};

// Per-phase blocking timings (alloc/u2d/kernel/write/free) are gone: once
// chunk N+1's prep genuinely overlaps chunk N's GPU work, measuring each
// phase in isolation would require re-inserting the very blocking syncs the
// pipeline is designed to avoid. setup/prep/drain remain meaningful because
// they are either purely host-side (chrono, no GPU sync) or already
// one-time/sequential by construction. For phase-level GPU inspection, use
// the NVTX ranges still present in run() with Nsight Systems.
struct TimingResult {
    double setup_ms = 0; ///< One-time pre-loop initialization (chrono).
    double prep_ms  = 0; ///< Sum of per-chunk host prep time (chrono, host-only).
    double loop_ms  = 0; ///< Wall-clock time of the pipelined chunk loop (chrono).
    double drain_ms = 0; ///< Final chunk finalization + writer thread join (chrono).
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

    // Device buffers, in two slots (slot = chunk index % 2) so that chunk N+1
    // can be prepared and uploaded into one slot while chunk N is still being
    // read by the GPU from the other slot. d_maps is dataset-constant and
    // read-only after setup, so it is not duplicated: both slots' kernels can
    // safely read it concurrently.
    char*         d_output[2]       = {nullptr, nullptr}; ///< Sparse device buffer for VCF text. Size: [CHUNK_SIZE * MAX_LINE_LEN].
    unsigned int* d_line_lens[2]    = {nullptr, nullptr}; ///< Device array storing the exact string length of each generated row.

    DeviceVarColumns d_df1[2];    ///< Shadow structure holding device pointers for DF1 (core variants).
    DeviceAltColumns d_df2[2];    ///< Shadow structure holding device pointers for DF2 (alternative alleles).
    DeviceSampleColumns d_df3[2]; ///< Shadow structure holding device pointers for DF3/DF4 (samples).
    DeviceMaps d_maps;            ///< Shadow structure holding device pointers for encoding dictionaries. Not slotted (read-only dataset constant).

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

    // Host staging and grouping (slotted the same way as the device buffers).
    HostBuffers host_buffers[2];  ///< Struct managing staging buffers on the host before transferring to the device.

    std::vector<GroupInfo> df3_int_groups;   ///< Grouping logic for DF3 integer fields (handles multi-value keys).
    std::vector<GroupInfo> df3_float_groups; ///< Grouping logic for DF3 float fields.
    std::vector<GroupInfo> df4_int_groups;   ///< Grouping logic for DF4 integer fields.
    std::vector<GroupInfo> df4_float_groups; ///< Grouping logic for DF4 float fields.

    // Compaction and Disk I/O variables (slotted; see d_output/d_df1 above).
    char* d_compacted[2] = {nullptr, nullptr};    ///< Device buffer for the compacted, contiguous VCF text. Persistent, grown on demand.
    unsigned int* d_output_offsets[2] = {nullptr, nullptr}; ///< Device array storing prefix-summed offsets for stream compaction.

    // --- Persistent device buffer capacities (one per slot) ---
    // allocateDevice() no longer runs cudaMalloc/cudaFree every chunk. Fixed-size
    // buffers (bounded by CHUNK_SIZE) are allocated once in initDeviceBuffers().
    // Variable-length buffers below persist across chunks and only grow (never
    // shrink) when a chunk needs more space than currently allocated.
    size_t cap_id_data[2] = {0, 0}, cap_ref_data[2] = {0, 0}, cap_in_string_data[2] = {0, 0};
    size_t cap_df2_var_id[2] = {0, 0}, cap_alt_data[2] = {0, 0}, cap_alt_offsets[2] = {0, 0};
    size_t cap_alt_int[2] = {0, 0}, cap_alt_float[2] = {0, 0}, cap_alt_string_data[2] = {0, 0}, cap_alt_string_offsets[2] = {0, 0};
    size_t cap_compacted[2] = {0, 0};
    size_t cap_temp_storage[2] = {0, 0};
    unsigned char* d_temp_storage[2] = {nullptr, nullptr}; ///< Persistent scratch storage for cub::DeviceScan, grown on demand.

    // Explicit CUDA streams (created in run(), destroyed at teardown). Chunk i
    // is prepared, uploaded and launched on streams[i % 2] and device/host
    // buffer slot i % 2, letting chunk i+1's host prep run while chunk i-1's
    // GPU work (launched on the other stream, one iteration ago) finishes up.
    cudaStream_t streams[2] = {nullptr, nullptr};

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
     * @brief Allocates the device buffers whose size only depends on CHUNK_SIZE and
     * per-dataset field counts (chromosome/INFO/FORMAT field counts, sample count),
     * never on a specific chunk's content. Called once before the main chunk loop.
     *
     * @param df1 Core variants DataFrame (used for INFO field counts).
     * @param df2 Alternative alleles DataFrame (used for INFO-ALT field counts).
     * @param df3 Sample core DataFrame (used for the sample count).
     */
    void initDeviceBuffers(const var_columns_df& df1,
                           const alt_columns_df& df2,
                           const sample_columns_df& df3);

    /**
     * @brief Ensures the variable-length device buffers (whose required size depends
     * on the current chunk's content: ID/REF/ALT/INFO string totals, allele count)
     * are large enough for the current chunk, growing them only when needed.
     *
     * @param df1 Core variants DataFrame.
     * @param df2 Alternative alleles DataFrame.
     * @param chunk_start Starting index in DF1.
     * @param chunk_end Ending index in DF1.
     * @param df2_start Starting index in DF2 corresponding to chunk_start.
     * @param slot Which of the two buffer slots (0 or 1) to prepare.
     */
    void allocateDevice(const var_columns_df& df1,
                        const alt_columns_df& df2,
                        int chunk_start,
                        int chunk_end,
                        int df2_start,
                        int slot);

    /**
     * @brief Safely deallocates all CUDA device memory to prevent memory leaks.
     * Called once after the main chunk loop (and from the destructor as a safety net).
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
     * @param slot Which of the two buffer slots (0 or 1) to upload into; also selects streams[slot].
     */
    void uploadToDevice(const var_columns_df& df1,
                        const alt_columns_df& df2,
                        const HostBuffers& buffers,
                        int slot);

    /**
     * @brief Performs stream compaction using Prefix Sum (CUB) and writes the processed, contiguous chunk to disk.
     *
     * @param num_variants The number of variants processed in the current chunk.
     * @param slot Which of the two buffer slots (0 or 1) holds this chunk's reconstructed data.
     */
    void writeChunk(int num_variants, int slot);

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
                        PinnedCharVec& buffer,
                        PinnedUintVec& offsets);

    /**
     * @brief Pins (page-locks) the DF1/DF2 host arrays that are copied directly
     * to the device every chunk in uploadToDevice() (var_number, chrom, pos,
     * qual, filter, var_id). Called once at the start of run(); paired with
     * unregisterHostInputBuffers() at the end.
     *
     * @param df1 Core variants DataFrame.
     * @param df2 Alternative alleles DataFrame.
     */
    void registerHostInputBuffers(const var_columns_df& df1, const alt_columns_df& df2);

    /**
     * @brief Undoes registerHostInputBuffers(); must be called before df1/df2's
     * underlying memory could be freed or reallocated by the caller.
     */
    void unregisterHostInputBuffers(const var_columns_df& df1, const alt_columns_df& df2);

};

#endif