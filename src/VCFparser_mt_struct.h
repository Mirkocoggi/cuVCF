#ifndef VCF_STRUCTS_H
#define VCF_STRUCTS_H
#include <chrono>
#include <boost/algorithm/string.hpp>


using namespace std;

class var
{
public:
    int var_number;
	string chrom="\0";
    int pos;
    string id="\0";
    string ref="\0";
    string alt="\0";
    int qual;
    string filter="\0";
    string info="\0";
    string format="\0";
    string samples;
    void get_vcf_line(char *line, int start, int end)
    {
        //da cambirare fare con char*
        bool find1 = false;
        int iter=0;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                chrom += line[start+iter];
                iter++;
            }
        }
        string tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                pos = stoi(tmp); // da cambiare, ci sarà un modo più semplice??!!
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                id += line[start+iter];
                iter++;
            }
        }
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                ref += line[start+iter];
                iter++;
            }
        }
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                alt += line[start+iter];
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                qual = stoi(tmp); // da cambiare, ci sarà un modo più semplice??!!
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                filter += line[start+iter];
                iter++;
            }
        }
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                info += line[start+iter];
                iter++;
            }
        }
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                format += line[start+iter];
                iter++;
            }
        }
        
        while((start+iter)!=(end-1)){
            samples += line[start+iter];
            iter++;
        }
        
        
    }
    void print_var()
    {
        cout << "Var" << var_number << ":\t";
        cout << chrom << "\t";
        cout << to_string(pos) << "\t";
        cout << id << "\t";
        cout << ref << "\t";
        cout << alt << "\t";
        cout << to_string(qual) << "\t";
        cout << filter << "\t";
        cout << info << "\t";
        cout << format << "\t";
        cout << samples;
        cout << endl;
    }
};

class vcf_parsed
{
public:
    int id;
    string filename;
    var *var_df;
    string header;
    char *filestring;
    int header_char=0;
    long filesize;
    long variants_size;
    long num_lines=0;
    long *new_lines_index;
    void get_filename(string path_to_filename){
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
    void get_file_size(string filename){
        filesize = filesystem::file_size(filename);
    }
    void get_header(ifstream *file){
        string line;
        //removing the header and storing it in vcf.header
        while (getline(*file, line) && line[0]=='#' && line[1]=='#'){
            header.append(line + '\n');
            header_char += line.length() + 1;
        }
        header_char += line.length() + 1;
        cout << "\nheader char: " << to_string(header_char) << endl;
        variants_size = filesize - header_char; //ora che ho tolto l'header ho un file piu piccolo quindi una nuova size
        cout<<"filesize: "<<filesize<<" variants_size: "<<variants_size<<endl;
    }
    void print_header(){
        cout << "VCF header:\n" << header << endl;
    }
    void allocate_filestring(){
        filestring = (char*)malloc(variants_size);
    }
    void find_new_lines_index(ifstream *inFile){ //popola anche il filestring
        long num_char=0;
        new_lines_index = (long*)malloc(variants_size); //per ora ho esagerato con la dimensione (è come se permettessi tutti \n. Si puo ridurre, pero ipotizzarlo è meglio perche senno devo passare il file due volte solo per vedere dove iniziano le linee)
        new_lines_index[0] = 0;
        num_lines++;
        while(num_char!=variants_size){
            filestring[num_char] = (*inFile).get(); 
            if(filestring[num_char]=='\n'){
                new_lines_index[num_lines] = num_char+1;
                num_lines++;
            }
            num_char++;
        }
    }
    void populate_var_struct(int num_threads){
        auto before = chrono::system_clock::now();
        var_df = (var*)calloc((num_lines-1), sizeof(var));
        cout << "\nBegin tmp: \n" <<"newlines: "<<num_lines<<" num threads: "<<num_threads<<endl;
        int batch_size = (num_lines-2+num_threads)/num_threads;
        cout << "\nBatch size: "<<batch_size<<endl;
        auto after = chrono::system_clock::now();
        auto pre_pragma = std::chrono::duration<double>(after - before).count();
        cout << "Pre_pragma: " << pre_pragma << " s" << endl;
#pragma omp parallel
        {
            int start, end;
            //string tmp ="\0";
            int th_ID = omp_get_thread_num();
            cout << "\nThread id: "<<th_ID<<endl;
            start = th_ID*batch_size;
            end = start + batch_size;
            cout << "\nstart: " << start << " end: " << end << endl;
            for(int i=start; i<end && i<num_lines-1; i++){
                //cout << "\nnew_lines_index[num_lines]: "<<new_lines_index[i] <<" new_lines_index[num_lines+1]: "<<new_lines_index[i+1]<<endl;
                // for(int j=new_lines_index[i]; j<new_lines_index[i+1]-1; j++){
                //     tmp += filestring[j]; //questo è un passaggio in piu che si puo togliere passando direttamente la sottostringa di filestring a get_vcf_line
                //     //cout<<filestring[j];
                // }
                // cout << "\nNew lines i: "<<new_lines_index[i]<<" New line i+1: "<<new_lines_index[i+1];
                var_df[i].get_vcf_line(filestring, new_lines_index[i], new_lines_index[i+1]);
                var_df[i].var_number = i;
                //var_df[i].print_var();
                //tmp.clear();
            
            }
        
        
        }
    }

};

#endif