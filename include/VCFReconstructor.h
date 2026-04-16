/**
 * @file VCFReconstructor.h
 * @brief Defines the class responsible for reconstructing a VCF file from separated DataFrames.
 *
 * This file contains the VCFReconstructor class, which handles the merging and 
 * formatting of variant core data (DF1), alternate alleles (DF2), sample-level 
 * data (DF3), and per-allele sample data (DF4) back into a standard VCF text format.
 * Reconstruction is parallelized using OpenMP: variants are split into chunks,
 * each processed by a separate thread writing to a temporary file, which are then
 * concatenated in order into the final output.
 */

#ifndef VCF_RECONSTRUCTOR_H
#define VCF_RECONSTRUCTOR_H

#include <string>
#include <map>
#include <vector>
#include "VCFDataFrames.h"

/**
 * @brief Reconstructs a standard VCF (Variant Call Format) file from parsed column-oriented DataFrames.
 *
 * Parallelizes reconstruction via OpenMP by splitting the variant set into chunks.
 * Each thread processes its own chunk independently using pre-computed cursors into
 * DF2 and DF4, writing results to a temporary file. The main thread then concatenates
 * the temporary files in order to produce the final VCF output.
 */
class VCFReconstructor {
public:
    /**
     * @brief Constructs a new VCFReconstructor object.
     *
     * @param output_vcf_path Path to the destination VCF file.
     * @param header_text The raw VCF header text to be written at the beginning of the file.
     */
    VCFReconstructor(const std::string& output_vcf_path,
                     const std::string& header_text);

    /**
     * @brief Executes the VCF reconstruction process.
     *
     * Opens the output file, writes the header, builds necessary inverse mapping
     * dictionaries, pre-computes per-chunk cursors for DF2 and DF4, launches parallel
     * OpenMP workers each writing to a temporary file, and finally concatenates the
     * temporary files in order into the final output.
     *
     * @param df1 DataFrame containing variant core fields and basic INFO fields.
     * @param df2 DataFrame containing alternate alleles and per-allele INFO fields.
     * @param df3 DataFrame containing core sample formats and genotypes.
     * @param df4 DataFrame containing sample data specific to alternate alleles.
     */
    void run(const var_columns_df& df1,
             const alt_columns_df& df2,
             const sample_columns_df& df3,
             const alt_format_df& df4);

private:
    std::string output_vcf_path; ///< Path where the reconstructed VCF will be saved.
    std::string header_text;     ///< Raw VCF header string to prepend to the output file.
    std::string format_str;      ///< Cached FORMAT string built from sample data fields (e.g., "GT:AD:DP").

    bool has_gt;                 ///< True if Genotype (GT) data exists anywhere in the dataset.
    bool gt_in_df3;              ///< True if Genotype (GT) data is stored in the core sample DataFrame (DF3).
    bool gt_in_df4;              ///< True if Genotype (GT) data is stored in the per-allele DataFrame (DF4).

    std::map<char, std::string> inv_chrom_map;             ///< Decodes compact chromosome char keys back to original strings.
    std::map<char, std::string> inv_filter_map;            ///< Decodes compact filter char keys back to original strings.
    std::vector<std::string> ordered_samp_names;           ///< List of sample names, strictly ordered by their numerical IDs.
    std::map<std::string, std::string> format_numbers;     ///< Maps FORMAT field IDs to their Number attribute (e.g., "AD" -> "R").

    /**
     * @brief Holds the pre-computed start and end positions for a single thread's chunk.
     *
     * Cursors into DF2 and DF4 are pre-computed sequentially before parallel execution
     * so that each thread can start traversal at the correct offset without scanning
     * from the beginning of the DataFrame.
     */
    struct ChunkCursor {
        size_t var_start;  ///< Index of the first variant in DF1 assigned to this chunk.
        size_t var_end;    ///< Index one past the last variant in DF1 assigned to this chunk.
        size_t df2_start;  ///< Index in DF2 where this chunk's variants begin.
        size_t df4_start;  ///< Index in DF4 where this chunk's variants begin.
    };

    /**
     * @brief Pre-computes the DF2 and DF4 starting offsets for each parallel chunk.
     *
     * Performs a single sequential O(n) scan over DF2 and DF4 to determine where
     * each chunk begins. This avoids redundant scanning inside parallel workers and
     * ensures that each thread can traverse its portion of DF2 and DF4 independently
     * without shared state.
     *
     * @param df1 DataFrame containing variant core fields (used to determine chunk boundaries).
     * @param df2 DataFrame containing alternate alleles (scanned to find per-chunk offsets).
     * @param df4 DataFrame containing per-allele sample data (scanned to find per-chunk offsets).
     * @param num_threads Number of chunks to generate, one per thread.
     * @return A vector of ChunkCursor structs, one per thread, in ascending order.
     */
    std::vector<ChunkCursor> computeChunkCursor(
        const var_columns_df& df1,
        const alt_columns_df& df2,
        const alt_format_df& df4,
        int num_threads);

    /**
     * @brief Builds inverse mapping dictionaries for encoded fields (CHROM and FILTER).
     *
     * Extracts the forward maps from DF1 and reverses them for fast char-to-string 
     * lookup during formatting.
     *
     * @param df1 Reference to DF1 containing the original forward maps.
     */
    void buildInverseMaps(const var_columns_df& df1);

    /**
     * @brief Reconstructs the ordered list of sample names from DF3.
     *
     * Inverts the sampNames map (name -> id) into a vector indexed by id,
     * ensuring sample columns are written in the correct order in the VCF output.
     *
     * @param df3 Reference to DF3 containing the sample name to ID mapping.
     */
    void buildSampleNames(const sample_columns_df& df3);

    /**
     * @brief Parses the VCF header to extract FORMAT field Number declarations.
     *
     * Scans the header text for ##FORMAT lines and extracts the ID and Number
     * attributes. The resulting map is used during reconstruction to determine
     * how many values to print for split fields (e.g., Number=R prints num_alts+1
     * values, Number=A prints num_alts values).
     *
     * @param header_text The raw VCF header text.
     * @return A map where keys are FORMAT field IDs (e.g., "AD") and values are
     *         their Number attribute (e.g., "R", "A", "1", ".").
     */
    std::map<std::string, std::string> parseFormatNumbers(const std::string& header_text);

    /**
     * @brief Formats a single variant row into a standard VCF string.
     *
     * Merges data from all four DataFrames for a specific variant index. Assembles
     * CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, and per-sample FORMAT columns
     * into a single tab-separated VCF record. Uses thread-local cursors passed by
     * reference to traverse DF2 and DF4 without shared state across threads.
     *
     * @param index The positional index of the variant in DF1.
     * @param df1 DataFrame containing variant core fields.
     * @param df2 DataFrame containing alternate alleles.
     * @param df3 DataFrame containing core sample formats.
     * @param df4 DataFrame containing per-allele sample formats.
     * @param df2_cursor Thread-local cursor into DF2, advanced as ALT rows are consumed.
     * @param df4_cursor Thread-local cursor into DF4, advanced as per-allele rows are consumed.
     * @return A fully formatted string representing one tab-separated line of the VCF body.
     */
    std::string formatVariant(int index,
                              const var_columns_df& df1,
                              const alt_columns_df& df2,
                              const sample_columns_df& df3,
                              const alt_format_df& df4,
                              size_t& df2_cursor,
                              size_t& df4_cursor);
};

#endif // VCF_RECONSTRUCTOR_H