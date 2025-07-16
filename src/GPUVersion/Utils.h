/**
 * @class vcf_parsed
 * @brief Encapsulates the VCF file parsing workflow
 * @ingroup Parser
 *
 * @details This class manages:
 *  - VCF file reading and header extraction
 *  - Host and device memory allocation
 *  - CUDA kernel execution
 *  - Results merging and cleanup
 *
 * The parser supports:
 *  - Standard VCF fields (CHROM, POS, etc.)
 *  - INFO field parsing
 *  - FORMAT field parsing
 *  - Sample data processing
 *  - Compressed (.gz) input files
 *
 * @note All device memory is automatically managed
 * @warning Requires sufficient GPU memory for file size
 */

#ifndef UTILS_H
#define UTILS_H

#include <sys/wait.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>
#include <boost/algorithm/string.hpp>

/// Constant representing a flag type.
const int FLAG = 0;
/// Constant representing an integer type.
const int INT = 1;
/// Constant representing a float type.
const int FLOAT = 2;
/// Constant representing a string type.
const int STRING = 3;
/// Constant representing an alternative integer type.
const int INT_ALT = 4;
/// Constant representing an alternative float type.
const int FLOAT_ALT = 5;
/// Constant representing an alternative string type.
const int STRING_ALT = 6;
/// Constant representing a formatted string type.
const int STRING_FORMAT = 8;
/// Constant representing a formatted integer type.
const int INT_FORMAT = 9;
/// Constant representing a formatted float type.
const int FLOAT_FORMAT = 10;
/// Constant representing an alternative formatted string type.
const int STRING_FORMAT_ALT = 11;
/// Constant representing an alternative formatted integer type.
const int INT_FORMAT_ALT = 12;
/// Constant representing an alternative formatted float type.
const int FLOAT_FORMAT_ALT = 13;

// Mappa per PolyPhen
const std::unordered_map<std::string, char> polyphenCharMap = {
    {"", 0},
    {"benign", 1},
    {"possibly_damaging", 2},
    {"probably_damaging", 3}
};

// Mappa per CSQ
const std::unordered_map<std::string, char> csqCharMap = {
    {"synonymous_variant", 0},
    {"missense_variant", 1},
    {"nonsense_variant", 2},
    {"frameshift_variant", 3},
    {"inframe_insertion", 4},
    {"inframe_deletion", 5},
    {"5_prime_UTR_variant", 6},
    {"3_prime_UTR_variant", 7},
    {"splice_region_variant", 8},
    {"splice_acceptor_variant", 9},
    {"splice_donor_variant", 10},
    {"regulatory_region_variant", 11},
    {"upstream_gene_variant", 12},
    {"downstream_gene_variant", 13},
    {"coding_sequence_variant", 14},
    {"non_coding_transcript_variant", 15},
    {"mature_miRNA_variant", 16},
    {"ncRNA_exon_variant", 17},
    {"ncRNA_intronic_variant", 18},
    {"protein_altering_variant", 19},
    {"transcript_ablation", 20},
    {"transcript_amplification", 21},
    {"intron_variant", 22},
    {"intergenic_variant", 23},
    {"enhancer_variant", 24},
    {"regulatory_region_variant", 25},
    {"non_coding_exon_variant", 26},
    {"pseudogene_variant", 27}
};

/**
 * @brief Securely unzips a .gz file using gzip
 *
 * @param vcf_filename [in,out] Pointer to filename, .gz extension removed on success
 * @throw std::runtime_error If fork() or exec() fails
 * @note Uses fork() and execlp() for secure execution
 * @warning Modifies the input filename string on successful decompression
 */
void unzip_gz_file(char* vcf_filename) {
    // Check if the filename ends with ".gz"
    if (strcmp(vcf_filename + strlen(vcf_filename) - 3, ".gz") == 0) {
        pid_t pid = fork();
        if (pid == 0) {
            // Child process: execute gzip command
            execlp("gzip", "gzip", "-df", vcf_filename, nullptr);
            // If execlp fails, output error and exit.
            std::cout << "ERROR: Failed to execute gzip command" << std::endl;
            exit(EXIT_FAILURE);
        } else if (pid > 0) {
            // Parent process: wait for the child process to finish.
            int status;
            waitpid(pid, &status, 0);
            if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
                // On success, remove the ".gz" extension from the filename.
                char* mutable_vcf_filename = const_cast<char*>(vcf_filename);
                mutable_vcf_filename[strlen(vcf_filename) - 3] = '\0';
            } else {
                std::cout << "ERROR: cannot unzip file" << std::endl;
            }
        } else {
            // Fork failed: output error message.
            std::cout << "ERROR: Failed to fork process" << std::endl;
        }
    }
}

/**
 * @brief Extracts the filename from a full file path.
 *
 * This function splits the provided file path using "/" as a delimiter and returns
 * the last token, which is assumed to be the filename. The full path is also assigned
 * to the reference parameter.
 *
 * @param path_filename The full file path as a string.
 * @param path_to_filename Reference to a string that will be assigned the full file path.
 * @return The extracted filename.
 */
string get_filename(string path_filename, string &path_to_filename){
    vector<string> line_el;
    path_to_filename = path_filename;
    boost::split(line_el, path_filename, boost::is_any_of("/"));
    return line_el[line_el.size()-1];
}

/**
 * @brief Gets the size of a file
 *
 * @param filename Full path to the file
 * @return long File size in bytes
 * @throw std::filesystem::filesystem_error If file doesn't exist or is inaccessible
 * @note Uses std::filesystem::file_size
 */
long get_file_size(string filename){
    return std::filesystem::file_size(filename);
}

/**
 * @brief Merges member vectors from temporary data structures
 * @tparam T Type of temporary objects containing vectors to merge
 * @tparam U Type of elements in the vectors
 * 
 * @param tmp_alt [in] Vector of temporary objects
 * @param dest [out] Destination vector for merged elements
 * @param num_threads Number of temporary objects to process
 * @param member_ptr Pointer to member vector within T
 * 
 * @note Uses move semantics to avoid copying
 * @warning Source vectors are cleared after merging
 */
template <typename T, typename U>
void merge_member_vector(
    std::vector<T>& tmp_alt,
    std::vector<U>& dest,
    int num_threads,
    std::vector<U> T::* member_ptr
) {
    for (int i = 0; i < num_threads; i++) {
        dest.insert(
            dest.end(),
            std::make_move_iterator((tmp_alt[i].*member_ptr).begin()),
            std::make_move_iterator((tmp_alt[i].*member_ptr).end())
        );
        (tmp_alt[i].*member_ptr).clear();
    }
}

/**
 * @brief Merges nested member vectors from temporary data structures into corresponding destination vectors.
 *
 * This function handles cases where each temporary object contains an outer vector (e.g., representing
 * multiple alternative fields) and each element of that outer vector is itself a vector that needs to be merged.
 * For every temporary object and for each element in the outer vector (up to num_nested elements), the function
 * moves the contents of the nested vector into the corresponding nested vector in the destination object and clears
 * the source nested vector afterwards.
 *
 * @tparam T The type of the temporary objects.
 * @tparam S The type of the elements in the outer vector (e.g., a structure representing a field group).
 * @tparam V The type of the elements in the inner (nested) vectors.
 * @param tmp_alt A vector of temporary objects containing nested member vectors.
 * @param dest The destination vector (global) where the nested vectors will be merged.
 * @param num_threads The number of temporary objects (typically equal to tmp_alt.size()).
 * @param num_nested The number of elements in the outer vector (e.g., the number of alternative fields).
 * @param outer_member_ptr Pointer to the outer vector member within the temporary object T.
 * @param inner_member_ptr Pointer to the inner vector member within the outer element S.
 */
template <typename T, typename S, typename V>
void merge_nested_member_vector(
    std::vector<T>& tmp_alt,
    std::vector<S>& dest,
    int num_threads,
    int num_nested,
    std::vector<S> T::* outer_member_ptr,
    std::vector<V> S::* inner_member_ptr
) {
    for (int i = 0; i < num_threads; i++) {
        for (int j = 0; j < num_nested; j++) {
            (dest[j].*inner_member_ptr).insert(
                (dest[j].*inner_member_ptr).end(),
                std::make_move_iterator(((tmp_alt[i].*outer_member_ptr)[j].*inner_member_ptr).begin()),
                std::make_move_iterator(((tmp_alt[i].*outer_member_ptr)[j].*inner_member_ptr).end())
            );
            ((tmp_alt[i].*outer_member_ptr)[j].*inner_member_ptr).clear();
        }
    }
}
 

#endif
