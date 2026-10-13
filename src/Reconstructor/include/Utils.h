/**
 * @file Utils.h
 * @brief Utility functions for VCF field name manipulation and CSV parsing.
 *
 * Provides helper functions used by both the CSV parser and the VCF reconstructor
 * for handling field name conventions from cuVCF's CSV exports.
 */

#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

/**
 * @brief Splits a CSV line into fields using the given delimiter.
 *
 * Handles quoted fields (double quotes) so that delimiters inside
 * quotes are not treated as separators. Strips trailing carriage
 * returns for Windows line endings.
 *
 * @param line The raw CSV line to split.
 * @param delimiter The character used to separate fields (typically ',').
 * @return A vector of strings, one per field.
 */
std::vector<std::string> splitLine(const std::string& line, char delimiter);

/**
 * @brief Strips trailing digits from a field name.
 *
 * Used to recover the base field name from cuVCF's numbered columns.
 * Example: "AD0" → "AD", "DP" → "DP".
 *
 * @param field_name The field name potentially ending with digits.
 * @return The field name without trailing digits.
 */
std::string stripTrailingDigits(const std::string& field_name);

/**
 * @brief Strips the type prefix from a CSV column header.
 *
 * cuVCF prefixes column names with their data type.
 * Example: "flag_SOMATIC" → "SOMATIC", "int_AN" → "AN".
 * If no known prefix is found, returns the name unchanged.
 *
 * @param name The CSV column header.
 * @return The field name without the type prefix.
 */
std::string stripTypePrefix(const std::string& name);

/**
 * @brief Checks whether a field name refers to the genotype (GT) field.
 *
 * Matches both "GT" and "GT0" (the latter used by cuVCF when GT
 * has a numeric suffix).
 *
 * @param name The field name to check.
 * @return True if the field is a GT field.
 */
bool isGTField(const std::string& name);

#endif // UTILS_H