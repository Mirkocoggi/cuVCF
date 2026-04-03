#ifndef UTILS_H
#define UTILS_H

#include <string>
#include "VCFDataFrames.h"

   
    /*
     * Splits a CSV line into fields using the given delimiter.
     * Returns a vector of strings, one for each field.
     * This simple splitter does not handle quotes or escaped delimiters,
     * which is fine for the CSV exports we expect.
     */
    std::vector<std::string> splitLine(const std::string& line, char delimiter);


   /*
     * Strips trailing digits from a field name (e.g. "AD0" -> "AD", "PL12" -> "PL")
     */
    std::string stripTrailingDigits(const std::string& field_name);

    /*
     * Strips type prefix from CSV field name (e.g. "flag_SOMATIC" -> "SOMATIC")
     */
    std::string stripTypePrefix(const std::string& name);

    /*
     * Merges split FORMAT fields (e.g. AD0, AD1 → AD) by concatenating values with comma
     */
    void mergeFields(sample_columns_df& df3);  

    /*
     * Moves the GT field to the first position in the sample columns
     */
    void moveGTFirst(std::vector<samp_String>& fields);

    /*
     * Checks if a field is the GT field
     */
    bool isGTField(const std::string& name);

#endif