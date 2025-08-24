/**
 * @file main.cu
 * @brief Entry point for the GPU-accelerated VCF parser
 * @author Your Name
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

#include "Parser.cu"        

#include <cuda_runtime.h>   
#include <cuda_fp16.h>      
#include <thrust/device_ptr.h> 
#include <thrust/sort.h>

#include <boost/algorithm/string.hpp> 
#include <chrono>           
#include <fstream>          
#include <filesystem>       
#include <sys/wait.h>       
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
 *  1. Sets up CUDA device
 *  2. Processes command-line arguments
 *  3. Initializes VCF parser
 *  4. Runs parsing operation
 *  5. Cleans up resources
 */
int main(int argc, char *argv[]){

    cudaSetDevice(0);  // Use device 0
   
    int opt, num_threadss;
    char *vcf_filename; 
   
    while ((opt = getopt(argc, argv, "v:t:")) != -1)
    {
        switch (opt)
        {
        case 'v':
            vcf_filename = optarg;
            break;
        case 't':
            num_threadss = atoi(optarg);
            if (num_threadss == 1)
            {
                cout << "Single thread execution, sequential process!!" << endl;
            }
            else
            {
                cout << "Multithreading execution, parallelization on " << num_threadss << " threads!!" << endl;
            }
            break;
        case '?':
            cout << "Unknown option: " << optopt << endl;
            break;
        }
    }
    
    vcf_parsed vcf;
    vcf.run(vcf_filename, num_threadss);

    return 0;
}