#ifndef CSV_PARSER_H
#define CSV_PARSER_H

#include <string>
#include <vector>
#include <map>
#include "VCFDataFrames.h"

/**
 * @brief Mapping from a CSV column to a typed vector in the dataframe.
 */
struct ColumnMapping {
    std::string type;   // "flag", "int", "float", "string", "gt"
    size_t index;       // index within the corresponding typed vector
};

/**
 * @brief Parses CSV files exported by cuVCF into cuVCF's typed DataFrames.
 *
 * Routes each CSV column to the correct typed vector (in_flag, in_int, in_float, etc.)
 * based on the type prefix in the column header (flag_, int_, float_, string_).
 */
class CSVParser {
public:
    std::string header_text; ///< Raw VCF header text loaded from header_path.
    std::map<std::string, std::string> field_types; ///< Maps field names to their type prefix (e.g., "AN" -> "int").

    /**
     * @brief Constructs a CSVParser and loads the VCF header from disk.
     *
     * @param df1_path    Path to the CSV file for core variant data (DF1).
     * @param df2_path    Path to the CSV file for alternate allele data (DF2).
     * @param df3_path    Path to the CSV file for core sample FORMAT data (DF3).
     * @param df4_path    Path to the CSV file for per-allele sample FORMAT data (DF4).
     * @param header_path Path to the file containing the raw VCF header text.
     */
    CSVParser(const std::string& df1_path,
                   const std::string& df2_path,
                   const std::string& df3_path,
                   const std::string& df4_path,
                   const std::string& header_path);

    /**
     * @brief Parses all four CSV DataFrames and populates the corresponding structures.
     *
     * @param df1 Output: core variant DataFrame.
     * @param df2 Output: alternate allele DataFrame.
     * @param df3 Output: core sample FORMAT DataFrame.
     * @param df4 Output: per-allele sample FORMAT DataFrame.
     */
    void loadAll(var_columns_df& df1,
                 alt_columns_df& df2,
                 sample_columns_df& df3,
                 alt_format_df& df4);

    /**
     * @brief Loads encoding dictionaries (CHROM, FILTER, etc.) from a maps file.
     *
     * @param maps_path Path to the serialized maps file produced by cuVCF.
     * @param df1       Output: receives the chrom_map and filter_map dictionaries.
     * @param df3       Output: receives the sample name map.
     * @param df4       Output: receives the sample name map for DF4.
     */
    void parseMaps(const std::string& maps_path,
                   var_columns_df& df1,
                   sample_columns_df& df3,
                   alt_format_df& df4);

private:
    std::string df1_path, df2_path, df3_path, df4_path, header_path;

    /** @brief Reads the header file and stores its content in header_text. */
    void parseHeader();
    /** @brief Parses the DF1 CSV and populates df1 with core variant fields. */
    void parseDF1(var_columns_df& df1);
    /** @brief Parses the DF2 CSV and populates df2 with alternate allele fields. */
    void parseDF2(alt_columns_df& df2);
    /** @brief Parses the DF3 CSV and populates df3 with core sample FORMAT fields. */
    void parseDF3(sample_columns_df& df3);
    /** @brief Parses the DF4 CSV and populates df4 with per-allele sample FORMAT fields. */
    void parseDF4(alt_format_df& df4);

    /**
     * @brief Extracts the type prefix from a CSV column name.
     * @return pair of (prefix, field_name), e.g. ("flag", "EVA_4")
     *         If no known prefix, returns ("string", original_name)
     */
    std::pair<std::string, std::string> extractTypeAndName(const std::string& header);
};

#endif