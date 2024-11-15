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
#include <OpenEXR/half.h>
#include <omp.h>
#include <sys/wait.h>
#include <unistd.h>
#include "VCFparser_mt_col_struct.h"
#include "VCF_parsed.h"
#include "VCF_var.h"
#include "VCF_var_columns_df.h"

using namespace std;

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
    //vcf.var_columns.print(10); //TODO
    vcf.alt_columns.print(10);
    vcf.samp_columns.print(10);
    vcf.alt_sample.print(10); //TODO
    /*
// Setting number of threads
    omp_set_num_threads(num_threadss);

// Open input file
    // gzip -df compressed_file1.gz
    if(!strcmp((vcf_filename + strlen(vcf_filename) - 3), ".gz")){
        unzip_gz_file(vcf_filename);
    }

    string filename = vcf_filename; 
    ifstream inFile(filename);
    if(!inFile){
        cout << "ERROR: cannot open file " << filename << endl;
        return -1;
    }
    //inFile.seekg(5033); //lunghezza header, stampa un po' di righe e basta, per vedere come sono fatte
    //std::string test;
    //getline(inFile, test);
    //cout << test << endl;
    //return -1;
    //---------------------
   
    int cont=0;
    string line;
    vcf_parsed vcf;
   
// Saving filename
    
    vcf.get_filename(filename);
    
// Getting filesize (number of char in the file)
    auto before = chrono::system_clock::now();
    vcf.get_file_size(filename);
    auto after = chrono::system_clock::now();
    auto get_file_size = std::chrono::duration<double>(after - before).count();

// Getting the header (Saving the header into a string and storing the header size )
    before = chrono::system_clock::now();
    vcf.get_and_parse_header(&inFile); //serve per separare l'header dal resto del file
    //vcf.print_header();
    after = chrono::system_clock::now();
    auto get_header = std::chrono::duration<double>(after - before).count();
    inFile.close();

// Allocating the filestring (the variations as a big char*, the dimension is: filesize - header_size)
    before = chrono::system_clock::now();
    vcf.allocate_filestring();

   
// Populate filestring and getting the number of lines (num_lines), saving the starting char index of each lines
    // before = chrono::system_clock::now();
    vcf.find_new_lines_index(filename, num_threadss);
    // after = chrono::system_clock::now();
    // auto find_new_lines = std::chrono::duration<double>(after - before).count();


// PRINT FROM FILESTRING
    
    // cout<<"\nnum_char: "<<vcf.variants_size<<endl;
    // for(long i=0; i<vcf.variants_size; i++){
    //     cout << vcf.filestring[i];
    // }

// PRINT NEW LINES INDEX

    // cout << "\nnew_line_inde: \n";
    // for(long i=0; i<vcf.num_lines; i++){
    //     cout << vcf.new_lines_index[i] << endl;
    // }

// Populating var structure:     

    // before = chrono::system_clock::now();
    // vcf.populate_var_struct(num_threadss);
    // after = chrono::system_clock::now();
    // auto populate_var_struct = std::chrono::duration<double>(after - before).count();

    // before = chrono::system_clock::now();
    vcf.create_info_vectors(num_threadss);
    //vcf.print_info_map();
    //vcf.print_info();
    vcf.reserve_var_columns();
    vcf.create_sample_vectors(num_threadss);
    after = chrono::system_clock::now();
    auto reserve_var_columns = std::chrono::duration<double>(after - before).count();
    before = chrono::system_clock::now();
    vcf.populate_var_columns(num_threadss);

    after = chrono::system_clock::now();
    auto populate_var_columns = std::chrono::duration<double>(after - before).count();
/* cout << "\nPrint from var_df: \n";
    for(int i=0; i<vcf.num_lines-1; i++){
       vcf.var_df[i].print_var();
    }

    cout << "\nPrint from var_columns: \n";
    vcf.var_columns.print_var_columns(10);

    cout << "\nPrint from var_df: \n";
    for(int i=0; i<5; i++){
        vcf.var_df[i].print_var();
    }
    cout << "\n\n\nLAST 100: \n\n\n";
    for(int i=vcf.num_lines - 101; i<vcf.num_lines - 1; i++){
        vcf.var_df[i].print_var();
    }
    cout << endl;

    vcf.samp_columns.print();
    vcf.alt_sample.print();
    vcf.alt_columns.print();
    
    cout << "Get file size: " << get_file_size << " s" << endl;
    cout << "get_header: " << get_header << " s" << endl;
    cout << "find_new_lines: " << find_new_lines << " s" << endl;
    cout << "populate_var_struct: " << populate_var_struct << " s" << endl;
    cout << "reserve: " << reserve_var_columns << " s" << endl;
*/
    //cout << "parsing time: " << populate_var_columns << " s" << endl;

    //free(vcf.var_df);
    //free(vcf.filestring);
    //free(vcf.new_lines_index);


    return 0;
}