
#ifndef VCF_DATAFRAMES_H
#define VCF_DATAFRAMES_H

#include <Imath/half.h>
#include <vector>
#include <string>
#include <map>


/**
 * @brief Constants defining the types of VCF fields
 * 
 * These constants are used to identify different types of fields in the VCF file,
 * including basic types (FLAG, INT, FLOAT, STRING) and their variants for ALT
 * and FORMAT sections.
 */
const int FLAG = 0;
const int INT = 1;
const int FLOAT = 2;
const int STRING = 3;
const int INT_ALT = 4;
const int FLOAT_ALT = 5;
const int STRING_ALT = 6;
const int STRING_FORMAT = 8;
const int INT_FORMAT = 9;
const int FLOAT_FORMAT = 10;
const int STRING_FORMAT_ALT = 11;
const int INT_FORMAT_ALT = 12;
const int FLOAT_FORMAT_ALT = 13;

using namespace std;

/** @brief Stores a named boolean (flag) INFO column across all variants. */
struct info_flag
{
    vector<bool> i_flag; ///< Per-variant flag values.
    string name;         ///< Field name (e.g., "SOMATIC").
};

/** @brief Stores a named string INFO column across all variants. */
struct info_string
{
    vector<string> i_string; ///< Per-variant string values.
    string name;             ///< Field name.
};

/** @brief Stores a named half-precision float INFO column across all variants. */
struct info_float
{
    vector<half> i_float; ///< Per-variant float values (fp16).
    string name;          ///< Field name.
};

/** @brief Stores a named integer INFO column across all variants. */
struct info_int
{
    vector<int> i_int; ///< Per-variant integer values. -1 encodes a missing value.
    string name;       ///< Field name.
};

/** @brief Stores a named boolean FORMAT column for a specific allele group across all (variant, sample) pairs. */
struct samp_Flag
{
    vector<bool> i_flag; ///< Per-(variant, sample) flag values.
    string name;         ///< Field name.
    int numb;            ///< Allele index or column group number.
};

/** @brief Stores a named string FORMAT column for a specific allele group across all (variant, sample) pairs. */
struct samp_String
{
    vector<string> i_string; ///< Per-(variant, sample) string values.
    string name;             ///< Field name.
    int numb;                ///< Allele index or column group number.
};

/** @brief Stores a named half-precision float FORMAT column for a specific allele group across all (variant, sample) pairs. */
struct samp_Float
{
    vector<half> i_float; ///< Per-(variant, sample) float values (fp16).
    string name;          ///< Field name.
    int numb;             ///< Allele index or column group number.
};

/** @brief Stores a named integer FORMAT column for a specific allele group across all (variant, sample) pairs. */
struct samp_Int
{
    vector<int> i_int; ///< Per-(variant, sample) integer values. -1 encodes a missing value.
    string name;       ///< Field name.
    int numb;          ///< Allele index or column group number.
};

/** @brief Stores encoded genotype (GT) values across all (variant, sample) pairs. */
struct samp_GT
{
    vector<char> GT; ///< Encoded GT values; each char maps to a GT string via GTMap.
    int numb;        ///< Number of samples or allele group index.
};

/**
 * @brief Aggregated metadata parsed from the VCF ##INFO and ##FORMAT header lines.
 *
 * Holds the IDs, Numbers, and Types declared in the header together with
 * pre-computed counts of each field category, used to guide DataFrame
 * allocation and reconstruction.
 */
struct header_element
{
    vector<string> ID;
    vector<string> Number;
    vector<string> Type;
    int total_values=0;
    int alt_values=0;
    int no_alt_values=0;
    int ints_alt=0;
    int floats_alt=0;
    int strings_alt=0;
    int flags_alt=0;
    int ints=0;
    int floats=0;
    int strings=0;
    int flags=0;
    bool hasGT = false;
    char numGT = 0;
};

/**
 * @class var_columns_df
 * @brief Column-oriented data frame for VCF variant data
 * 
 * Stores VCF data in a column-oriented format for efficient processing and memory usage.
 * Each field from the VCF file (CHROM, POS, ID, etc.) is stored in its own vector,
 * allowing for better cache utilization and SIMD processing.
 */

class var_columns_df
{
public:
    /** @brief Vector storing variant numbers for each entry */
    vector<unsigned int> var_number;
    
    /** 
     * @brief Maps chromosome names to compact char codes
     * @details Provides efficient storage by mapping chromosome names (e.g., "chr1") 
     *          to single character codes
     */
    std::map<std::string, char> chrom_map;
    
    /** @brief Vector storing chromosome codes sequentially */
    vector<char> chrom;
    
    /** @brief Vector storing genomic positions */
    vector<unsigned int> pos;
    
    /** @brief Vector storing variant identifiers */
    vector<string> id;
    
    /** @brief Vector storing reference alleles */
    vector<string> ref;
    
    /** @brief Vector storing quality scores in half precision */
    vector<half> qual;

    /** @brief Maps filter names to compact char representations */
    std::map<std::string, char> filter_map;
    
    /** @brief Vector storing filter status codes */
    vector<char> filter;
    
    /** @brief Vector storing float-type INFO fields */
    vector<info_float> in_float;
    
    /** @brief Vector storing flag-type INFO fields */
    vector<info_flag> in_flag;
    
    /** @brief Vector storing string-type INFO fields */
    vector<info_string> in_string;
    
    /** @brief Vector storing integer-type INFO fields */
    vector<info_int> in_int;
    
    /** @brief Maps INFO field names to their indices */
    map<string,int> info_map1;
 
};

/**
 * @brief Column-oriented data frame for alternate allele data (DF2).
 *
 * Each row corresponds to one alternate allele entry.  Multi-allelic variants
 * produce multiple consecutive rows sharing the same var_id.
 */
class alt_columns_df
{
    public:
    vector<unsigned int> var_id;    ///< Variant ID this allele belongs to.
    vector<char>         alt_id;    ///< Allele index within its variant (0-based).
    vector<string>       alt;       ///< ALT allele string (e.g., "A", "<DEL>").
    vector<info_float>   alt_float; ///< Allele-specific float INFO fields.
    vector<info_flag>    alt_flag;  ///< Allele-specific flag INFO fields (not yet handled).
    vector<info_string>  alt_string;///< Allele-specific string INFO fields.
    vector<info_int>     alt_int;   ///< Allele-specific integer INFO fields.
    int numAlt;                     ///< Total number of allele rows in this DataFrame.
};

/**
 * @brief Column-oriented data frame for core per-sample FORMAT fields (DF3).
 *
 * Stores FORMAT data that is the same across all alternate alleles for a given
 * (variant, sample) pair. The flat layout is indexed as
 * [variant_index * numSample + sample_index].
 */
class sample_columns_df
{
    public:
    vector<unsigned int>  var_id;      ///< Variant ID for each row.
    vector<unsigned short> samp_id;   ///< Sample ID for each row.
    vector<samp_Float>    samp_float; ///< Named float FORMAT columns.
    vector<samp_Flag>     samp_flag;  ///< Named flag FORMAT columns.
    vector<samp_String>   samp_string;///< Named string FORMAT columns.
    vector<samp_Int>      samp_int;   ///< Named integer FORMAT columns.
    std::map<std::string, unsigned short> sampNames; ///< Maps sample name to its integer ID.
    map<string, char>     GTMap;       ///< Maps GT strings (e.g., "0|1") to compact char codes.
    vector<samp_GT>       sample_GT;  ///< Encoded genotype column; present when GT is in DF3.
    int numSample; ///< Number of samples per variant row.

    void initMapGT(){
        int value = 0;
        // First half of the map from 0|0 to 10|10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "|" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
        // Second half of the map from 0/0 to 10/10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "/" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
    }

    std::string getGTStringFromChar(char gtChar) const {
        for (const auto& pair : GTMap) {
            if (pair.second == gtChar) {
                return pair.first;
            }
        }
        return "Not found";
    }
};

/**
 * @brief Column-oriented data frame for per-allele sample FORMAT fields (DF4).
 *
 * Stores FORMAT data that varies per alternate allele for each (variant, sample)
 * pair. Each row is uniquely identified by the combination of var_id, samp_id,
 * and alt_id.
 */
class alt_format_df
{
    public:
    vector<unsigned int>  var_id;     ///< Variant ID for each row.
    vector<unsigned short> samp_id;  ///< Sample ID for each row.
    vector<char>           alt_id;   ///< Allele index for each row.
    vector<samp_Float>     samp_float; ///< Named float FORMAT columns (allele-specific).
    vector<samp_Flag>      samp_flag;  ///< Named flag FORMAT columns (allele-specific).
    vector<samp_String>    samp_string;///< Named string FORMAT columns (allele-specific).
    vector<samp_Int>       samp_int;   ///< Named integer FORMAT columns (allele-specific).
    std::map<std::string, unsigned short> sampNames; ///< Maps sample name to its integer ID.
    samp_GT               sample_GT; ///< Encoded genotype column; present when GT is in DF4.
    map<string, char>     GTMap;     ///< Maps GT strings to compact char codes.
    int numSample; ///< Number of samples per variant row.


    void initMapGT(){
        int value = 0;
        // First half of the map from 0|0 to 10|10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "|" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
        // Second half of the map from 0/0 to 10/10
        for (int i = 0; i < 11; ++i) {
            for (int j = 0; j < 11; ++j) {
                std::string key = std::to_string(i) + "/" + std::to_string(j);
                GTMap[key] = value;
                value++;
            }
        }
    }

    std::string getGTStringFromChar(char gtChar) const {
        for (const auto& pair : GTMap) {
            if (pair.second == gtChar) {
                return pair.first;
            }
        }
        return "Not found";
    }
};

#endif