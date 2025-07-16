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
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <filesystem>
#include <Imath/half.h>
#include <omp.h>
#include <sys/wait.h>
#include <unistd.h>
#include "VCFparser_mt_col_struct.h"
#include "VCF_parsed.h"
#include "VCF_var.h"
#include "VCF_var_columns_df.h"

using namespace std;

/**
 * @brief Safely decompresses a gzipped VCF file
 * 
 * @param vcf_filename Path to the gzipped VCF file
 * 
 * @details Uses fork() and execlp() to safely execute gzip for decompression.
 * The function handles process management and error conditions:
 *   1. Checks for .gz extension
 *   2. Forks a child process for decompression
 *   3. Parent process waits for completion
 *   4. Updates filename after successful decompression
 * 
 * @note Uses system calls that are POSIX-compliant
 */
void unzip_gz_file(char* vcf_filename) {
    // Check the extension ".gz"
    if (strcmp(vcf_filename + strlen(vcf_filename) - 3, ".gz") == 0) {
        pid_t pid = fork();
        if (pid == 0) {
            // Child process
            execlp("gzip", "gzip", "-df", vcf_filename, nullptr);
            // If execlp fails
            cout<< "ERROR: Failed to execute gzip command" << std::endl;
            exit(EXIT_FAILURE);
        } else if (pid > 0) {
            // Parent process waits for the child process
            int status;
            waitpid(pid, &status, 0);
            if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
                // Remove ".gz" from the filename
                char* mutable_vcf_filename = const_cast<char*>(vcf_filename);
                mutable_vcf_filename[strlen(vcf_filename) - 3] = '\0';
            } else {
                cout<< "ERROR: cannot unzip file" << std::endl;
            }
        } else {
            // Error in the fork() call
            cout<< "ERROR: Failed to fork process" << std::endl;
        }
    }
}

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