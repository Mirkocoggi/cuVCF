//#include "VCFparser_mt_col_struct_CU.h"
//#include "VCF_CUDA_implementation.cu"
#include "VCF_parsed.cu"
//#include "VCF_var_columns_df_CU.h"

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

//Securely unzip the file
void unzip_gz_file(char* vcf_filename) {
    // Check the extension ".gz"
    if (strcmp(vcf_filename + strlen(vcf_filename) - 3, ".gz") == 0) {
        pid_t pid = fork();
        if (pid == 0) {
            // Child proces
            execlp("gzip", "gzip", "-df", vcf_filename, nullptr);
            // If execlp fails
            cout<< "ERROR: Failed to execute gzip command" << std::endl;
            exit(EXIT_FAILURE);
        } else if (pid > 0) {
            // Parent proces waits for the child proces
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
    cout << "AAAAAA 1 " << endl;
    vcf.alt_columns.print(10);
    cout << "AAAAAA 2 " << endl;
    vcf.samp_columns.print(10); //TOTO: controlla come riprendi i dati da GPU
    cout << "AAAAAA 3 " << endl;
    vcf.alt_sample.print(10);
    cout << "AAAAAA 4 " << endl;

    return 0;
}