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
#include <omp.h>
#include "VCFparser_mt_struct.h"

using namespace std;

int main(int argc, char *argv[]){
   
    int opt, num_threadss;
    char *vcf_filename; 
   
    while ((opt = getopt(argc, argv, "v:t:")) != -1)
    {
        switch (opt)
        {
        case 'v':
            vcf_filename = optarg;
            cout << "vcf_filename: " << vcf_filename << endl;
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

    omp_set_num_threads(num_threadss);

    string filename = vcf_filename;
    ifstream inFile(filename);
    if(!inFile){
        cout << "ERROR: cannot open file " << filename << endl;
        return -1;
    }
    int cont=0;
    string line;
    vcf_parsed vcf;
    vcf.get_filename(filename);
    vcf.get_file_size(filename);
    vcf.get_header(&inFile); //serve per separare l'header dal resto del file
    vcf.allocate_filestring();
    vcf.find_new_lines_index(&inFile);

    cout<<"\nnum_char: "<<vcf.variants_size<<endl;
    for(long i=0; i<vcf.variants_size; i++){
        cout << vcf.filestring[i];
    }
    cout << "\nnew_line_inde: \n";
    for(long i=0; i<vcf.num_lines; i++){
        cout << vcf.new_lines_index[i] << endl;
    }

    vcf.populate_var_struct(num_threadss);
    cout << "\nPrint from var_df: \n";
    for(int i=0; i<vcf.num_lines-1; i++){
        vcf.var_df[i].print_var();
    }
    inFile.close();
    
    
    return 0;
}