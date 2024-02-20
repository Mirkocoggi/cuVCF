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
#include "VCFparser_serial_struct.h"

using namespace std;

int main(int argc, char *argv[]){
   
    int opt;
    char *vcf_filename; 
   
    while ((opt = getopt(argc, argv, "v:")) != -1)
    {
        switch (opt)
        {
        case 'v':
            vcf_filename = optarg;
            cout << "vcf_filename: " << vcf_filename << endl;
            break;
        case '?':
            cout << "Unknown option: " << optopt << endl;
            break;
        }
    }

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
    vcf.get_header(&inFile);
    //vcf.print_header();
    vcf.create_var_struct(&inFile);
    
    
    
    return 0;
}