/**
 * @file VCFReconstructor.h
 * @brief Defines the class responsible for reconstructing a VCF file from separated DataFrames.
 *
 * This file contains the VCFReconstructor class, which handles the merging and 
 * formatting of variant core data (DF1), alternate alleles (DF2), sample-level 
 * data (DF3), and per-allele sample data (DF4) back into a standard VCF text format.
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
 * It uses internal tracking cursors (df2_cursor, df4_cursor) to efficiently 
 * traverse the DataFrames and align variant IDs without relying on deep nested loops.
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
     * Opens the output file, writes the header, builds necessary inverse mapping dictionaries, 
     * and iterates through all variants to format and write each VCF record to disk.
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
    
    size_t df2_cursor;           ///< Running index to track the current positional offset in DF2 without resetting.
    size_t df4_cursor;           ///< Running index to track the current positional offset in DF4 without resetting.

    std::map<char, std::string> inv_chrom_map;  ///< Decodes compact chromosome char keys back to original strings.
    std::map<char, std::string> inv_filter_map; ///< Decodes compact filter char keys back to original strings.
    std::vector<std::string> ordered_samp_names;///< List of sample names, strictly ordered by their numerical IDs.

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
     * @param df3 Reference to DF3 containing the sample name to ID mapping.
     */
    void buildSampleNames(const sample_columns_df& df3);

    /**
     * @brief Formats a single variant row into a standard VCF string format.
     *
     * Merges data from all four DataFrames for a specific variant index. It aligns 
     * multiple alternate alleles, merges heterogeneous INFO fields, and reconstructs 
     * the sample FORMAT columns.
     *
     * @param index The positional index of the variant in DF1.
     * @param df1 DataFrame containing variant core fields.
     * @param df2 DataFrame containing alternate alleles.
     * @param df3 DataFrame containing core sample formats.
     * @param df4 DataFrame containing per-allele sample formats.
     * @return A fully formatted string representing one tab-separated line of the VCF body.
     */
    std::string formatVariant(int index,
                              const var_columns_df& df1,
                              const alt_columns_df& df2,
                              const sample_columns_df& df3,
                              const alt_format_df& df4);
};

#endif // VCF_RECONSTRUCTOR_H