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

#include <zlib.h>
#include <stdexcept>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

// Splits s at every sep into out: empty tokens are kept ("a,,b" gives 3) and an empty s gives one
// empty token
inline void split_on(std::vector<std::string>& out, std::string_view s, char sep){
    out.clear();
    size_t b = 0;
    for(size_t e; (e = s.find(sep, b)) != std::string_view::npos; b = e + 1) out.emplace_back(s.substr(b, e - b));
    out.emplace_back(s.substr(b));
}

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
const int FLAG_FORMAT = 17; // was 11, which collided with STRING_FORMAT_ALT


/**
 * @brief Decompresses a .gz file in place with zlib (like "gzip -df")
 *
 * @param vcf_filename [in,out] Pointer to filename, .gz extension removed on success
 * @warning Modifies the input filename string on successful decompression
 */
inline void unzip_gz_file(char* vcf_filename) {
    // Decompress in process with zlib (no external gzip binary, no fork),
    // with the same effect as "gzip -df": the .gz is replaced by the plain file.
    const size_t len = strlen(vcf_filename);
    if (len < 4 || strcmp(vcf_filename + len - 3, ".gz") != 0) return;
    const std::string out_name(vcf_filename, len - 3);

    gzFile in = gzopen(vcf_filename, "rb");
    FILE* out = in ? fopen(out_name.c_str(), "wb") : nullptr;
    bool ok = in && out;
    if (ok) {
        std::vector<char> buf(1 << 16);
        int n;
        while ((n = gzread(in, buf.data(), buf.size())) > 0) {
            if (fwrite(buf.data(), 1, n, out) != static_cast<size_t>(n)) { ok = false; break; }
        }
        // A truncated stream ends with gzread() == 0 like a clean EOF: only gzerror() reports it
        int zerr = Z_OK;
        gzerror(in, &zerr);
        if (n < 0 || zerr != Z_OK || gzdirect(in)) ok = false; // gzdirect: not gzip data
    }
    if (out && fclose(out) != 0) ok = false;
    if (in && gzclose(in) != Z_OK) ok = false;
    if (!ok) {
        if (out) remove(out_name.c_str()); // no partial output, and the .gz is kept
        throw std::runtime_error(std::string("cannot decompress ") + vcf_filename);
    }
    remove(vcf_filename);
    vcf_filename[len - 3] = '\0'; // continue with the decompressed file
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
inline std::string get_filename(std::string path_filename, std::string &path_to_filename){
    path_to_filename = path_filename;
    return path_filename.substr(path_filename.rfind('/') + 1); // the whole string when there is no '/'
}

/**
 * @brief Gets the size of a file
 *
 * @param filename Full path to the file
 * @return long File size in bytes
 * @throw std::filesystem::filesystem_error If file doesn't exist or is inaccessible
 * @note Uses std::filesystem::file_size
 */
inline long get_file_size(std::string filename){
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
    size_t total = dest.size();
    for (int i = 0; i < num_threads; i++) total += (tmp_alt[i].*member_ptr).size();
    dest.reserve(total);
    for (int i = 0; i < num_threads; i++) {
        dest.insert(
            dest.end(),
            std::make_move_iterator((tmp_alt[i].*member_ptr).begin()),
            std::make_move_iterator((tmp_alt[i].*member_ptr).end())
        );
        std::vector<U>().swap(tmp_alt[i].*member_ptr); // releases the memory; clear() keeps the capacity
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
    for (int j = 0; j < num_nested; j++) {
        size_t total = (dest[j].*inner_member_ptr).size();
        for (int i = 0; i < num_threads; i++) total += ((tmp_alt[i].*outer_member_ptr)[j].*inner_member_ptr).size();
        (dest[j].*inner_member_ptr).reserve(total);
    }
    for (int i = 0; i < num_threads; i++) {
        for (int j = 0; j < num_nested; j++) {
            (dest[j].*inner_member_ptr).insert(
                (dest[j].*inner_member_ptr).end(),
                std::make_move_iterator(((tmp_alt[i].*outer_member_ptr)[j].*inner_member_ptr).begin()),
                std::make_move_iterator(((tmp_alt[i].*outer_member_ptr)[j].*inner_member_ptr).end())
            );
            std::vector<V>().swap((tmp_alt[i].*outer_member_ptr)[j].*inner_member_ptr);
        }
    }
}
 

#endif
