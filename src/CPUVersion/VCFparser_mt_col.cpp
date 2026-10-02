/**
 * @file VCFparser_mt_col.cpp
 * @brief Multi-threaded VCF parser implementation with column-oriented storage
 * 
 * This implementation provides the main entry point for the VCF parser and
 * handles command-line argument processing and file decompression.
 * 
 * @note This is part of the CPU-only implementation and requires no CUDA dependencies
 */

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <map>
#include <tuple>
#include <zlib.h>
#include <queue>
#include <chrono>
#include <filesystem>
#include <Imath/half.h>
#include <omp.h>
#include <unistd.h>
#include "VCFparser_mt_col_struct.h"
#include "VCF_parsed.h"
#include "VCF_var_columns_df.h"

using namespace std;

/**
 * @brief Program entry point
 * 
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return int Exit status (0 for success)
 * 
 * @details Processes command line arguments:
 *   -v <filename>: VCF file to process
 *   -t <threads>:  Number of threads to use
 * 
 * Creates and runs a vcf_parsed instance to process the input file.
 */
int main(int argc, char *argv[]){

    int opt, num_threadss = 4;      // -t is optional
    char *vcf_filename = nullptr;   // -v is required

    while ((opt = getopt(argc, argv, "v:t:")) != -1)
    {
        switch (opt)
        {
        case 'v':
            vcf_filename = optarg;
            break;
        case 't':
            num_threadss = atoi(optarg);
            break;
        default: // '?': getopt already printed the error
            cerr << "Usage: " << argv[0] << " -v <file.vcf[.gz]> [-t <threads>]" << endl;
            return 1;
        }
    }
    if (vcf_filename == nullptr || num_threadss < 1) {
        cerr << "Usage: " << argv[0] << " -v <file.vcf[.gz]> [-t <threads>]" << endl;
        return 1;
    }

    if (num_threadss == 1)
    {
        cout << "Single thread execution, sequential process!!" << endl;
    }
    else
    {
        cout << "Multithreading execution, parallelization on " << num_threadss << " threads!!" << endl;
    }

    vcf_parsed vcf;
    auto s = std::chrono::steady_clock::now();
    try {
        vcf.run(vcf_filename, num_threadss);
    } catch (const std::exception& ex) {
        cerr << "ERROR: " << ex.what() << endl;
        return 1;
    }
    auto e = std::chrono::steady_clock::now();
    cerr << "vcf_parsed::run: " << std::chrono::duration<double, std::milli>(e - s).count() << " ms" << endl;

    // Preview of the four DataFrames
    cout << "------------------------------" << endl;
    vcf.var_columns.print(10);
    cout << "------------------------------" << endl;
    vcf.alt_columns.print(10);
    cout << "------------------------------" << endl;
    vcf.samp_columns.print(10);
    cout << "------------------------------" << endl;
    vcf.alt_sample.print(10);
    cout << "------------------------------" << endl;

    return 0;
}
