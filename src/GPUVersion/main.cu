/**
 * @file main.cu
 * @brief Entry point for the GPU-accelerated VCF parser
 * @date 2025-07-16
 *
 * @details Main application that:
 *  - Sets up CUDA device and environment
 *  - Processes command-line arguments
 *  - Initializes and runs the VCF parser
 *  - Manages resource cleanup
 *
 * Usage:
 *   ./VCFparser -v <vcf_filename> -t <num_threads>
 *
 * Arguments:
 *   -v : Path to input VCF file (required)
 *   -t : Number of CPU threads for parallel processing (required)
 *
 * @note CUDA device 0 is used by default
 */

#include "Parser.h"        

#include <cuda_runtime.h>   
#include <cuda_fp16.h>      

#include <chrono>           
#include <fstream>          
#include <filesystem>       
#include <unistd.h>         
#include <map>              
#include <omp.h>            
#include <iostream>         
#include <cstdlib>

using namespace std;


/**
 * @brief Program entry point
 * 
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return int Exit status (0 for success, non-zero for errors)
 *
 * @details Program workflow:
 *  1. Sets up CUDA device (in vcf_parsed::run)
 *  2. Processes command-line arguments
 *  3. Initializes VCF parser
 *  4. Runs parsing operation
 *  5. Cleans up resources
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
