/**
 * @file VCF_var.h
 * @brief Class definition for VCF variant records
 * 
 * This header defines the var class which represents a single variant record
 * from a VCF file. It provides storage and parsing capabilities for all standard
 * VCF fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, and samples).
 * 
 * @note This is part of the CPU-only implementation and requires no CUDA dependencies
 */

#ifndef VCF_VAR_H
#define VCF_VAR_H
#include "VCFparser_mt_col_struct.h"
#include "VCF_parsed.h"
#include "VCF_var_columns_df.h"
#include <chrono>
#include <boost/algorithm/string.hpp>
#include <Imath/half.h>

/**
 * @class var
 * @brief Represents a single variant record from a VCF file
 * 
 * Stores and manages all fields from a VCF record line including both
 * required fields (CHROM, POS, etc.) and optional fields (FORMAT, samples).
 */
class var
{
public:

    /** @brief Unique identifier for the variant */
    long var_number;

    /** @brief Chromosome identifier where the variant is located */
    string chrom="\0";
    
    /** @brief 1-based position of the variant on the chromosome */
    long pos;
    
    /** @brief Semi-colon separated list of unique identifiers */
    string id="\0";
    
    /** @brief Reference sequence at the given position */
    string ref="\0";
    
    /** @brief Comma-separated list of alternative alleles */
    string alt="\0";
    
    /** @brief Phred-scaled quality score for the assertion */
    half qual;
    
    /** @brief Filter status, PASS or filter-specific codes */
    string filter="\0";
    
    /** @brief Additional annotation information */
    string info="\0";
    
    /** @brief Format string for the sample genotype fields */
    string format="\0";
    
    /** @brief Sample-specific genotype information */
    string samples="\0";

    /**
     * @brief Parses a VCF line into its constituent fields
     * 
     * @param line Pointer to the raw VCF line buffer
     * @param start Starting position in the line buffer
     * @param end Ending position in the line buffer
     * 
     * @details Processes a single line from a VCF file, parsing each field
     * according to the VCF specification. Fields are extracted in order:
     * CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, and samples.
     * 
     * @note Special handling is implemented for:
     *  - QUAL field: Converts "." to 0.0 and handles invalid floats
     *  - Position: Converts string to unsigned long
     */
    void get_vcf_line(char *line, long start, long end)
    {
        bool find1 = false;
        long iter=0;
        string tmp="\0";
        
        //Chromosome
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                chrom += line[start+iter];
                iter++;
            }
        }
        
        //Position
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                pos = stoul(tmp); 
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //ID
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                id += line[start+iter];
                iter++;
            }
        }

        //Reference
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                ref += line[start+iter];
                iter++;
            }
        }

        //Alternatives
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                alt += line[start+iter];
                iter++;
            }
        }

        //Quality
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(strcmp(&tmp[0], ".")==0){
                    qual = 0.0;
                }else{
                    try{
                        qual = (half)stof(tmp);
                    }catch (const std::exception& e){
                        qual = 0;
                    }
                }
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }

        //Filter
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                filter += line[start+iter];
                iter++;
            }
        }

        //Info
        find1=false;
        while(!find1 && (start+iter)<(end-1)){
            if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
                find1 = true;
                iter++;
            }else{
                info += line[start+iter];
                iter++;
            }
        }

    }
    
    /**
     * @brief Prints the variant record in VCF format
     * 
     * @details Outputs all fields in tab-separated format following VCF specification:
     * VAR_NUMBER CHROM POS ID REF ALT QUAL FILTER INFO
     * Each field is separated by tabs, and the line ends with a newline.
     */
    void print_var()
    {
        cout << "Var" << var_number << ":\t";
        cout << chrom << "\t";
        cout << to_string(pos) << "\t";
        cout << id << "\t";
        cout << ref << "\t";
        cout << alt << "\t";
        cout << to_string(qual) << "\t";
        cout << filter << "\t";
        cout << info << endl;
    }
};

#endif