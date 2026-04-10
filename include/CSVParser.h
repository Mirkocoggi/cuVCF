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
    std::string header_text;
    std::map<std::string, std::string> field_types;
    header_element INFO;
    header_element FORMAT;

    CSVParser(const std::string& df1_path,
                   const std::string& df2_path,
                   const std::string& df3_path,
                   const std::string& df4_path,
                   const std::string& header_path);

    void loadAll(var_columns_df& df1,
                 alt_columns_df& df2,
                 sample_columns_df& df3,
                 alt_format_df& df4);

private:
    std::string df1_path, df2_path, df3_path, df4_path, header_path;

    void parseHeader();
    void parseDF1(var_columns_df& df1);
    void parseDF2(alt_columns_df& df2);
    void parseDF3(sample_columns_df& df3);
    void parseDF4(alt_format_df& df4);

    /**
     * @brief Extracts the type prefix from a CSV column name.
     * @return pair of (prefix, field_name), e.g. ("flag", "EVA_4")
     *         If no known prefix, returns ("string", original_name)
     */
    std::pair<std::string, std::string> extractTypeAndName(const std::string& header);
};

#endif