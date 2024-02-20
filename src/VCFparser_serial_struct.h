#ifndef VCF_STRUCTS_H
#define VCF_STRUCTS_H
#include <chrono>
#include <boost/algorithm/string.hpp>


using namespace std;

class var
{
public:
    int var_number;
	string chrom;
    int pos;
    string id;
    string ref;
    string alt;
    int qual;
    string filter;
    string info;
    string format;
    vector<string> samples;
    void get_vcf_line(string line)
    {
        vector<string> line_el;
        boost::split(line_el, line, boost::is_any_of("\t "));
        chrom = line_el[0];
        pos = stoi(line_el[1]);
        id = line_el[2];
        ref = line_el[3];
        alt = line_el[4];
        qual = stoi(line_el[5]);
        filter = line_el[6];
        info = line_el[7];
        format = line_el[8];
        string samples_string = line_el[9];
        line_el.clear();
        boost::split(line_el, samples_string, boost::is_any_of("\t "));
        for (int i=0; i<line_el.size(); i++){
            samples.push_back(line_el[i]);
        }
    }
    void print_var()
    {
        cout << "Var " << var_number << " : ";
        cout << chrom << "\t";
        cout << to_string(pos) << "\t";
        cout << id << "\t";
        cout << ref << "\t";
        cout << alt << "\t";
        cout << to_string(qual) << "\t";
        cout << filter << "\t";
        cout << info << "\t";
        cout << format << "\t";
        for (int i = 0; i < samples.size(); i++)
        {
            cout << samples[i] << "\t";
        }
        cout << endl;
    }
};

class vcf_parsed
{
public:
    int id;
    string filename;
    vector<var> var_df;
    string header;
    int num_var = 0;
    void get_filename(string path_to_filename)
    {
        vector<string> line_el;
        boost::split(line_el, path_to_filename, boost::is_any_of("/"));
        filename = line_el[line_el.size()-1];
        cout << "\nOnly file name: " << filename << endl;
        line_el.clear();
        boost::split(line_el, filename, boost::is_any_of("."));
        if(line_el[line_el.size()-1] == "gz")
        {
            cout << "\nFile must be uncompressed!\n" << endl;
            exit(0);
        }else{
            cout << "\nFile already uncompressed!\n" << endl;
        }
    }
    void get_header(ifstream *file)
    {
        string line;
        //removing the header and storing it in vcf.header
        while (getline(*file, line) && line[0]=='#' && line[1]=='#'){
            header.append(line + '\n');
        }
        cout << "\nLast line: "<< line<< endl;
        
    }
    void print_header()
    {
        cout << "VCF header:\n" << header << endl;
    }
    void create_var_struct(ifstream *file)
    {
        string line;
        var var_tmp;
        while (getline(*file, line))
        {
            var_tmp.get_vcf_line(line);
            num_var++;
            var_tmp.var_number = num_var;
            var_tmp.print_var();
            var_df.push_back(var_tmp);
            var_tmp.samples.clear();
        }
    }
};

#endif