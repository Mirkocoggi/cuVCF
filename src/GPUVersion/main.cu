/**
 * @file main.cu
 * @brief Entry point for the VCF parser application.
 *
 * This file sets up the CUDA device, processes command-line arguments to obtain the
 * VCF file path (-v) and the number of threads (-t) for parallel processing, and
 * invokes the VCF parser workflow through the vcf_parsed class.
 *
 * Usage:
 *   ./VCFparser -v <vcf_filename> -t <num_threads>
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

    //vcf.print_header();
    //cout << "Start printing:" << endl;
    //vcf.var_columns.print(10);
    //vcf.alt_columns.print(10);
    //vcf.samp_columns.print(16); //TOTO: controlla come riprendi i dati da GPU
    //vcf.alt_sample.print(10);

    return 0;
}