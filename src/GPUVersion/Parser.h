#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <map>
#include "DataStructures.h"
#include "DataFrames.h"

/**
 * @class vcf_parsed
 * @brief Encapsulates the VCF file parsing workflow.
 *
 * This class manages reading the VCF file, extracting the header and variant/sample data,
 * allocating memory on both host and device, launching CUDA kernels, and merging the results.
 */

class vcf_parsed
{
public:
    /// Unique identifier for this VCF parsing instance.
    int id;
    /// Name of the VCF file.
    string filename;
    /// Full path to the VCF file.
    string path_to_filename;
    /// The header content of the VCF file.
    string header;
    /// Parsed header information for INFO fields.
    header_element INFO;
    /// Mapping of INFO field names to type codes (Flag=0, Int=1, Float=2, String=3).
    map<string,int> info_map;
    /// Parsed header information for FORMAT fields.
    header_element FORMAT;
    /// Host-side storage for variant data as a large character array.
    char *filestring;
    /// Device-side storage for variant data.
    char *d_filestring;
    /// Device-side counter used for various kernel operations.
    unsigned int *d_count;
    /// Size of the VCF header (in bytes).
    int header_size = 0;
    /// Total file size (in bytes).
    long filesize;
    /// Size of the variant portion of the file (filesize minus header size).
    long variants_size;
    /// Number of variant lines (excluding the header).
    long num_lines = 0;
    /// Host-side array storing the starting index of each variant line.
    unsigned int *new_lines_index;
    /// Device-side array storing the starting index of each variant line.
    unsigned int *d_new_lines_index;
    /// Flag indicating whether sample data is present.
    bool samplesON = false;
    /// Flag indicating whether detailed sample data is available.
    bool hasDetSamples = false;
    /// Data frame containing variant column information (e.g., chromosome, position).
    var_columns_df var_columns;
    /// Data frame containing alternative allele information.
    alt_columns_df alt_columns;
    /// Data frame containing sample column information.
    sample_columns_df samp_columns;
    /// Data frame for alternative allele sample format data.
    alt_format_df alt_sample;
    /// Host-side structure containing kernel parameters.
    KernelParams h_params;
    /// Pointer to the device-side kernel parameters structure.
    KernelParams *d_params;
    /// Device memory for storing variant numbers.
    unsigned int *d_VC_var_number;
    /// Device memory for storing variant positions.
    unsigned int *d_VC_pos;
    /// Device memory for storing quality scores (in half precision).
    __half *d_VC_qual;
    /// Device-side structure for float INFO field values.
    info_float_d *d_VC_in_float = (info_float_d*)malloc(sizeof(info_float_d));
    /// Device-side structure for flag INFO field values.
    info_flag_d *d_VC_in_flag = (info_flag_d*)malloc(sizeof(info_flag_d));
    /// Device-side structure for integer INFO field values.
    info_int_d *d_VC_in_int = (info_int_d*)malloc(sizeof(info_int_d));
    /// Device memory for storing sample variant IDs.
    unsigned int *d_SC_var_id;
    /// Device memory for storing sample IDs.
    unsigned short *d_SC_samp_id;
    /// Device-side structure for sample float values.
    samp_Float_d *d_SC_samp_float = (samp_Float_d*)malloc(sizeof(samp_Float_d));
    /// Device-side structure for sample flag values.
    samp_Flag_d *d_SC_samp_flag = (samp_Flag_d*)malloc(sizeof(samp_Flag_d));
    /// Device-side structure for sample integer values.
    samp_Int_d *d_SC_samp_int = (samp_Int_d*)malloc(sizeof(samp_Int_d));
    /// Device-side structure for sample genotype values.
    samp_GT_d *d_SC_sample_GT = (samp_GT_d*)malloc(sizeof(samp_GT_d));

    /**
    * @brief Runs the VCF parsing process.
    *
    * Performs the following steps:
    *  - Initializes CUDA device and queries device properties.
    *  - Opens the VCF file (uncompressing if needed) and extracts the header.
    *  - Allocates memory for the file content and identifies the start of each variant line.
    *  - Creates and reserves vectors for variant and sample data.
    *  - Allocates and initializes device memory and lookup maps.
    *  - Launches CUDA kernels to parse the VCF lines.
    *  - Copies the parsed results back to host memory and frees device memory.
    *
    * @param vcf_filename Path to the VCF file.
    * @param num_threadss Number of threads to use for OpenMP parallel processing.
    */
    void run(char* vcf_filename, int num_threadss);
    
    /**
    * @brief Copies a genotype map to device constant memory.
    *
    * Copies the key-value pairs from the provided host map into the device constant
    * memory arrays (d_keys_gt and d_values_gt).
    *
    * @param map Host map with genotype keys and corresponding char values.
    */
    void copyMapToConstantMemory(const std::map<std::string, char>& map);

    /**
    * @brief Initializes the INFO field lookup map (Map1) in device constant memory.
    *
    * Copies keys and integer values from the host map to device constant memory.
    * Prints an error if the number of keys exceeds the allowed limit.
    *
    * @param my_map Host map containing keys and their corresponding integer values.
    */
    void initialize_map1(const std::map<std::string, int> &my_map);

    /**
    * @brief Allocates device memory for VCF parsing data.
    *
    * Allocates memory on the GPU for variant numbers, positions, quality scores,
    * and INFO fields. If sample data is present, it also allocates memory for sample fields.
    */
    void device_allocation();

    /**
    * @brief Frees all allocated device memory.
    *
    * Releases device memory allocated during the parsing process.
    */
    void device_free();

    /**
    * @brief Finds newline indices in the VCF file.
    *
    * Reads the VCF file in parallel using OpenMP to determine the starting index of each line.
    * The indices are stored in an array and copied to device memory for use by CUDA kernels.
    *
    * @param w_filename The path to the VCF file.
    * @param num_threads Number of threads to use for parallel processing.
    */
    void find_new_lines_index(string w_filename, int num_threads);
    
    /**
    * @brief Reads the VCF header from the input file.
    *
    * Extracts header lines (starting with "##") from the VCF file,
    * storing them in the header string and updating the header size.
    *
    * @param file Pointer to the input file stream.
    */
    void get_header(ifstream *file);
    
    /**
    * @brief Prints the VCF header to standard output.
    */
    void print_header();
    
    /**
    * @brief Reads and parses the VCF header.
    *
    * Separates header and variant data, extracts INFO and FORMAT information,
    * and determines the number of samples present.
    *
    * @param file Pointer to the input file stream.
    */
    void get_and_parse_header(ifstream *file);
    
    /**
    * @brief Allocates a character array to store the variant portion of the VCF file.
    *
    * The allocated size is based on the file size minus the header size.
    */
    void allocate_filestring();

    /**
    * @brief Creates and initializes vectors for sample data.
    *
    * Based on the FORMAT header, this method initializes vectors for sample genotype,
    * float, integer, and string data.
    *
    * @param num_threads Number of threads to use for parallel processing.
    */
    void create_sample_vectors(int num_threads);
    
    /**
    * @brief Creates and initializes vectors for INFO field data.
    *
    * Initializes vectors to store INFO field values (flag, integer, float, and string)
    * based on header information.
    *
    * @param num_threads Number of threads to use for parallel processing.
    */
    void create_info_vectors(int num_threads);
    
    /**
    * @brief Prints the INFO field mapping.
    *
    * Outputs the mapping from INFO field names to their corresponding type codes.
    */
    void print_info_map();
    
    /**
    * @brief Prints a summary of INFO field data.
    *
    * Displays a brief summary of the sizes and first few entries for each INFO field type.
    */
    void print_info();
    
    /**
    * @brief Reserves space in the variant columns structure.
    *
    * Resizes the vectors in the var_columns_df structure based on the number of variants.
    */
    void reserve_var_columns();

    /**
    * @brief Allocates device memory for KernelParams and copies the host structure.
    *
    * @param d_params Double pointer to the device KernelParams.
    * @param h_params Pointer to the host KernelParams structure.
    */
    void allocParamPointers(KernelParams **d_params, KernelParams *h_params);

    /**
    * @brief Launches the CUDA kernel to parse VCF lines and merges the results.
    *
    * Sets up CUDA streams and events, launches the appropriate kernel (based on whether sample data is present),
    * and asynchronously copies the parsed data from device to host.
    */
    void populate_runner(int numb_cores);

    /**
    * @brief Populates variant columns by processing VCF lines in parallel.
    *
    * Spawns a worker thread to run the CUDA kernel for parsing and uses OpenMP to merge alternative allele
    * data from multiple threads into the final data structures (alt_columns_df and alt_format_df).
    *
    * @param num_threads Number of threads to use for parallel merging.
    */
    void populate_var_columns(int num_threads, int numb_cores);

    /**
    * @brief Parses a VCF line and populates variant columns data.
    *
    * This function processes a single VCF line (from index @p start to @p end) by reading it character by character.
    * It extracts key variant fields such as chromosome, variant ID, reference allele, alternative alleles,
    * filter information, and the INFO field. The alternative allele field is split using commas and stored in
    * the provided alt_columns_df structure (@p tmp_alt). The count of alternative alleles processed is tracked by
    * the integer pointed to by @p tmp_num_alt.
    *
    * @param line Pointer to the VCF line as a C-string.
    * @param start The starting index of the line within the file.
    * @param end The ending index of the line within the file.
    * @param i The index (row number) corresponding to the current variant.
    * @param tmp_alt Pointer to an alt_columns_df structure for storing alternative allele data.
    * @param tmp_num_alt Pointer to an integer tracking the current number of alternative alleles processed.
    */
    void get_vcf_line_in_var_columns(char *line, long start, long end, long i, alt_columns_df* tmp_alt, int *tmp_num_alt);

    /**
    * @brief Parses a VCF line with formatted variant and sample data.
    *
    * This function processes a single VCF line to extract both variant and sample-related information,
    * following a predefined FORMAT template. It parses the chromosome, variant ID, reference allele,
    * alternative alleles, and filter field, and then further splits the FORMAT field to extract per-sample
    * data (e.g., genotype, float, integer, and string values). The extracted sample data is stored in the
    * provided alt_format_df structure (@p tmp_alt_format), while variant data is updated in the global structures.
    *
    * @param line Pointer to the VCF line as a C-string.
    * @param start The starting index of the line within the file.
    * @param end The ending index of the line within the file.
    * @param i The index (row number) corresponding to the current variant.
    * @param tmp_alt Pointer to an alt_columns_df structure for storing alternative allele data.
    * @param tmp_num_alt Pointer to an integer tracking the number of alternative alleles processed.
    * @param sample Pointer to a sample_columns_df structure for storing sample-specific data.
    * @param FORMAT Pointer to a header_element structure describing the FORMAT fields.
    * @param tmp_num_alt_format Pointer to an integer tracking the number of formatted alternative entries processed.
    * @param tmp_alt_format Pointer to an alt_format_df structure for storing formatted sample data.
    */
    void get_vcf_line_in_var_columns_format(char *line, long start, long end, long i, alt_columns_df* tmp_alt, int *tmp_num_alt, sample_columns_df* sample, header_element* FORMAT, int *tmp_num_alt_format, alt_format_df* tmp_alt_format);
    
};
 

#endif // PARSER_H
