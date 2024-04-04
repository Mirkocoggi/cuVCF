#ifndef VCF_STRUCTS_H
#define VCF_STRUCTS_H
#include <chrono>
#include <boost/algorithm/string.hpp>


using namespace std;

class info_flag
{
    public:
    vector<bool> i_flag;
    string name;
};

class info_string
{
    public:
    vector<string> i_string;
    string name;
};

class info_float
{
    public:
    vector<float> i_float;
    string name;
};

class info_int
{
    public:
    vector<int> i_int;
    string name;
};


class var_columns_df
{
public:
    vector<long> var_number;
    vector<string> chrom;
    vector<long>pos;
    vector<string> id;
    vector<string> ref;
    vector<string> alt;
    vector<float> qual;
    vector<string> filter;
    vector<string> info;
    vector<info_float> in_float;
    vector<info_flag> in_flag;
    vector<info_string> in_string;
    vector<info_int> in_int;
    map<string,int> info_map1;
    void get_vcf_line_in_var_columns(char *line, long start, long end, long i)
    {
        //cout<<i<<"\t";
        bool find1 = false;
        long iter=0;
        string tmp="\0";
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                chrom[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                pos[i] = stoul(tmp);
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                id[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                ref[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                alt[i] = tmp;
            }else{
                tmp += line[start+iter]; // per ora le salvo tutte come un unico array of char, andrebbe cambiato per salvarle separatamente
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(strcmp(&tmp[0], ".")==0){
                    qual[i] = 0.0f;
                }else{
                    qual[i] = stof(tmp); // da cambiare, in futuro
                }
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                filter[i] = tmp;
            }else{
                tmp += line[start+iter];
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
                find1 = true;
                iter++;
                info[i] = tmp;
                vector<string> tmp_el;
                boost::split(tmp_el, tmp, boost::is_any_of(";"));
                vector<string> tmp_elems;
                for(int ii=0; ii<tmp_el.size(); ii++){
                    boost::split(tmp_elems, tmp_el[ii], boost::is_any_of("="));
                    bool find_info_type = false;
                    bool find_info_elem = false;
                    if(tmp_elems.size()==2){
                        
                        while(!find_info_type){
                            
                            if(info_map1[tmp_elems[0]]==1){
                                //Int
                                
                                int el=0;
                                while(!find_info_elem){
                                    if(in_int[el].name == tmp_elems[0]){
                                        in_int[el].i_int[i] = stoi(tmp_elems[1]);
                                        find_info_elem = true;
                                    }
                                    el++; 
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==2){
                                //Float
                                
                                int el=0;
                                while(!find_info_elem){
                                    if(in_float[el].name == tmp_elems[0]){
                                        in_float[el].i_float[i] = stof(tmp_elems[1]);
                                        find_info_elem = true;
                                    } 
                                    el++;
                                }
                                find_info_type = true;
                            }else if(info_map1[tmp_elems[0]]==3){
                                //String
                                
                                int el=0;
                                while(!find_info_elem){
                                    if(in_string[el].name == tmp_elems[0]){
                                        in_string[el].i_string[i] = tmp_elems[1];
                                        find_info_elem = true;
                                    } 
                                    el++;
                                }
                                find_info_type = true;
                            }else{
                                cout<<"\n\nhere\n\n";
                                find_info_type = true;
                            }
                        }
                    }else{
                        if((info_map1[tmp_elems[0]]==0) && strcmp(&tmp_elems[0][0],"")){
                                //Flag
                                
                                int el=0;
                                while(!find_info_elem){
                                    if(in_flag[el].name == tmp_elems[0]){
                                        in_flag[el].i_flag[i] = 1;
                                        find_info_elem = true;
                                    } 
                                    el++;
                                }
                                
                    }
                    tmp_elems.clear();
                    }
                }
            }else{
                tmp += line[start+iter]; // anche qui andrebbero separate le info e creati dizionari con keys e values
                iter++;
                
            }
        }
    }

    
    void print_var_columns(long num_lines){
        for(long i=0; i<num_lines; i++){
            cout << "Var" << var_number[i] << ":\t";
            cout << chrom[i] << "\t";
            cout << to_string(pos[i]) << "\t";
            cout << id[i] << "\t";
            cout << ref[i] << "\t";
            cout << alt[i] << "\t";
            cout << to_string(qual[i]) << "\t";
            cout << filter[i] << "\t";
            cout << info[i] << "\t";
            for(int j=0; j<in_flag.size(); j++){
                cout<<in_flag[j].name<<": "<<in_flag[j].i_flag[i]<<", ";
            }
            for(int j=0; j<in_int.size(); j++){
                cout<<in_int[j].name<<": "<<in_int[j].i_int[i]<<", ";
            }
            for(int j=0; j<in_float.size(); j++){
                cout<<in_float[j].name<<": "<<in_float[j].i_float[i]<<", ";
            }
            for(int j=0; j<in_string.size(); j++){
                cout<<in_string[j].name<<": "<<in_string[j].i_string[i]<<", ";
            }
            cout<<endl;

        }
    }    
};

class var
{
public:
    long var_number;
	string chrom="\0"; // su quale cromosoma della reference
    long pos; // posizione nel chromosoma
    string id="\0"; // id della variation, spesso .
    string ref="\0"; // nucleotide della reference
    string alt="\0"; // possibili alternative (possono essere più di una)
    float qual; // quality score
    string filter="\0"; // filtro penso usato durante il sequenziamento
    string info="\0"; // info varie deducibili dall'header, potrebbe essere utile averle in collegamento
    string format="\0"; // formato dei samples, info variabiliti, mi dice come sono ordinate
    string samples="\0"; // una colonna per ogni sample in cui ognuno descrive i valori indicati in format
    void get_vcf_line(char *line, long start, long end)
    {
        bool find1 = false;
        long iter=0;
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
                pos = stoul(tmp); // da cambiare, in futuro
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
                alt += line[start+iter]; // per ora le salvo tutte come un unico array of char, andrebbe cambiato per salvarle separatamente
                iter++;
            }
        }
        tmp="\0";
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
                if(strcmp(&tmp[0], ".")==0){
                    qual = 0.0;
                }else{
                    qual = stof(tmp); // da cambiare, in futuro
                }
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
        while(!find1 && (start+iter)<(end-1)){
            if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
                find1 = true;
                iter++;
            }else{
                info += line[start+iter]; // anche qui andrebbero separate le info e creati dizionari con keys e values
                iter++;
            }
        }
        /*find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '||line[start+iter]=='\n'){
                find1 = true;
                iter++;
            }else{
                format += line[start+iter]; // vanno separati e soprattutti connessi ai valori nei samples (tipo se nel format ho ..:AF e nel sample ..:2, devo poter fare la query AF=2 in tutte le var)
                iter++;
            }
        }
        
        while((start+iter)<(end-1)){ //segmentation fault -> perchè non (start+iter)<(end-1)?
            samples += line[start+iter]; // sicuramente da dividere i vari sample e eventualmente creare un dizionario per ogni sample cosi da associare con il format
            iter++;
        }*/
        
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
        cout << info << endl;
        //cout << format << "\t";
        //cout << samples;
        //cout << endl;
    }
};

class header_element
{
public:
    vector<string> ID;
    vector<string> Number;
    vector<string> Type;
    int total_values=0;
    int alt_values=0;
    int no_alt_values=0;
    int ints_alt=0;
    int floats_alt=0;
    int strings_alt=0;
    int flags_alt=0;
    int ints=0;
    int floats=0;
    int strings=0;
    int flags=0;

};

class vcf_parsed
{
public:
    int id;
    string filename;
    var *var_df;
    string header;
    header_element INFO;
    map<string,int> info_map; //Flag=0, Int=1, Float=2, String=3;
    header_element FORMAT; //aggiungere Filter e gli altri, magari con switch case
    char *filestring;
    int header_size=0;
    long filesize;
    long variants_size;
    long num_lines=0;
    long *new_lines_index;
    var_columns_df var_columns;
    void get_filename(string path_to_filename){
        vector<string> line_el;
        boost::split(line_el, path_to_filename, boost::is_any_of("/"));
        filename = line_el[line_el.size()-1];
        //cout << "\nOnly file name: " << filename << endl;
        line_el.clear();
        boost::split(line_el, filename, boost::is_any_of("."));
        if(line_el[line_el.size()-1] == "gz")
        {
            //cout << "\nFile must be uncompressed!\n" << endl;
            exit(0);
        }else{
            //cout << "\nFile already uncompressed!\n" << endl;
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
            header_size += line.length() + 1;
        }
        header_size += line.length() + 1;
        //cout << "\nheader char: " << to_string(header_size) << endl;
        variants_size = filesize - header_size; //ora che ho tolto l'header ho un file piu piccolo quindi una nuova size
        //cout<<"filesize: "<<filesize<<" variants_size: "<<variants_size<<endl;
    }
    void print_header(){
        //cout << "VCF header:\n" << header << endl;
    }
    void get_and_parse_header(ifstream *file){
        string line;
        vector<string> line_el; //all the characteristics together
        vector<string> line_el1; //each characteristic
        vector<string> line_el2; //keys and values
        //removing the header and storing it in vcf.header
        
        while (getline(*file, line) && line[0]=='#' && line[1]=='#'){
            header.append(line + '\n');
            header_size += line.length() + 1;
            bool Info = (line[2]=='I');
            bool Format = (line[2]=='F' && line[3]=='O');
            
            if(Info || Format){
                boost::split(line_el, line, boost::is_any_of("><"));
                boost::split(line_el1, line_el[1], boost::is_any_of(","));
                for(int i=0; i<3; i++){
                    boost::split(line_el2, line_el1[i], boost::is_any_of("="));
                    if(Info){
                        if(i==0) INFO.ID.push_back(line_el2[1]);
                        if(i==1) INFO.Number.push_back(line_el2[1]);
                        if(i==1 && line_el2[1] == "A") INFO.alt_values++;
                        if(i==2) INFO.Type.push_back(line_el2[1]);
                    }
                    if(Format){
                        if(i==0) FORMAT.ID.push_back(line_el2[1]);
                        if(i==1) FORMAT.Number.push_back(line_el2[1]);
                        if(i==2) FORMAT.Type.push_back(line_el2[1]);
                    }
                }
            }
        }
        for(int i=0; i<INFO.ID.size(); i++){
            cout<<"INFO.ID["<<i<<"]: "<<INFO.ID[i]<<" | INFO.Number["<<i<<"]: "<<INFO.Number[i]<<" | INFO.Type["<<i<<"]: "<<INFO.Type[i]<<endl;
        }
        INFO.total_values = INFO.ID.size();
        INFO.no_alt_values = INFO.total_values - INFO.alt_values;
        cout<<"INFO.total_values: "<< INFO.total_values << " INFO.alt_values: "<<INFO.alt_values<<" INFO.no_alt_values: "<<INFO.no_alt_values<<endl;
        header_size += line.length() + 1;
        cout << "\nheader char: " << to_string(header_size) << endl;
        variants_size = filesize - header_size; //ora che ho tolto l'header ho un file piu piccolo quindi una nuova size
        cout<<"filesize: "<<filesize<<" variants_size: "<<variants_size<<endl;
    }
    
    void allocate_filestring(){
        filestring = (char*)malloc(variants_size);
    }
    void find_new_lines_index(string w_filename, int num_threads){ //popola anche il filestring
        long num_char=0;
        new_lines_index = (long*)malloc(variants_size); //per ora ho esagerato con la dimensione (è come se permettessi tutti \n. Si puo ridurre (euristicamente), pero ipotizzarlo è meglio perche senno devo passare il file due volte solo per vedere dove iniziano le linee)
        new_lines_index[0] = 0; //il primo elemento lo metto a zero per indicare l'inizio della prima linea
        num_lines++;
        auto before = chrono::system_clock::now();
        long batch_infile = (variants_size - 1 + num_threads)/num_threads; //numero di char che verrà processato da ogni thread

        /*-----------------------------------------------------
        std::ifstream file_streams[num_threads]; //std::vector<std::ifstream> file_streams(num_threads); // std::ifstream file_streams[num_threads] ??
        for (int i = 0; i < num_threads; i++) {
            // Set the offset for this thread
            int start = header_size + i*batch_infile;
            file_streams[i].open(w_filename); // Open the file stream for this thread
            file_streams[i].seekg(start, ios::cur); // Set the file stream's read position to the start offset
        }
        */        
#pragma omp parallel
        {
            int thr_ID = omp_get_thread_num();
            ifstream infile(w_filename); //apro lo stesso file con un ifstream in ogni thread, da considerare se ha senso creare prima un array di ifstream
            infile.seekg((header_size + thr_ID*batch_infile), ios::cur); //ogni thread parte a copiare inFile in filestring a una distanza di batch_size
            long start, end;
            start = thr_ID*batch_infile; // inizio del batch dello specifico thread
            end = start + batch_infile; // fine del batch
            //cout << "in find start: "<<start<<" end: "<<end<<endl;
            for(long i=start; i<end && i<variants_size; i++){
                filestring[i] = infile.get(); //file_streams[thr_ID].get(); //-----------------
            }
        }
        while(filestring[variants_size-1]=='\n'){
            variants_size--;
        }
        filestring[variants_size] = '\n';
        variants_size++;
        auto after = chrono::system_clock::now();
        auto filestring_time = std::chrono::duration<double>(after - before).count();
        before = chrono::system_clock
        ::now();
        num_char = 0;
        while(num_char!=variants_size){ //conto il numero di lines e ogni volta che trovo \n salvo il char successivo come inizio della riga successiva
            if(filestring[num_char]=='\n'&& filestring[num_char+1]!='\n'){
                new_lines_index[num_lines] = num_char+1; // PROBLEMA: funziona solo se l'ultimo char è uno /n, altrimenti si rompe
                num_lines++;
            }
            num_char++;
        } // si potrebbe fare multithreading ma vediamo se si ha beneficio
        after = chrono::system_clock::now();

        auto f_new_lines = std::chrono::duration<double>(after - before).count(); 
        //cout << /*"\nFilestring time: " <<*/ filestring_time /*<< " s " << "New lines time: " << f_new_lines << " s\n\n"*/ << endl; //da printare
    }
    void create_info_vectors(){
        info_flag info_flag_tmp;
        info_float info_float_tmp;
        info_int info_int_tmp;
        info_string info_string_tmp;
        for(int i=0; i<INFO.total_values; i++){
            if(strcmp(&INFO.Number[i][0], "A") == 0){
                if(strcmp(&INFO.Type[i][0], "Integer")==0){
                    INFO.ints_alt++;
                    info_map[INFO.ID[i]] = 1;
                }
                if(strcmp(&INFO.Type[i][0], "Float")==0){
                    INFO.floats_alt++;
                    info_map[INFO.ID[i]] = 2;
                }
                if(strcmp(&INFO.Type[i][0], "String")==0){
                    INFO.strings_alt++;
                    info_map[INFO.ID[i]] = 3;
                }
                if(strcmp(&INFO.Type[i][0], "Flag")==0){
                    INFO.flags_alt++;
                    info_map[INFO.ID[i]] = 0;
                }
            }else if((strcmp(&INFO.Number[i][0], "1") == 0)||(strcmp(&INFO.Number[i][0], "0") == 0)){
                if(strcmp(&INFO.Type[i][0], "Integer")==0){
                    INFO.ints++;
                    info_int_tmp.name = INFO.ID[i];
                    info_int_tmp.i_int.resize(num_lines-1, 0);
                    var_columns.in_int.push_back(info_int_tmp);
                    info_map[INFO.ID[i]] = 1;
                    var_columns.info_map1[INFO.ID[i]] = 1;
                }
                if(strcmp(&INFO.Type[i][0], "Float")==0){
                    INFO.floats++;
                    info_float_tmp.name = INFO.ID[i];
                    info_float_tmp.i_float.resize(num_lines-1, 0);
                    var_columns.in_float.push_back(info_float_tmp);
                    info_map[INFO.ID[i]] = 2;
                    var_columns.info_map1[INFO.ID[i]] = 2;
                }
                if(strcmp(&INFO.Type[i][0], "String")==0){
                    INFO.strings++;
                    info_string_tmp.name = INFO.ID[i];
                    info_string_tmp.i_string.resize(num_lines-1, "\0");
                    var_columns.in_string.push_back(info_string_tmp);
                    info_map[INFO.ID[i]] = 3;
                    var_columns.info_map1[INFO.ID[i]] = 3;
                }
                if(strcmp(&INFO.Type[i][0], "Flag")==0){
                    INFO.flags++;
                    info_flag_tmp.name = INFO.ID[i];
                    info_flag_tmp.i_flag.resize(num_lines-1, 0);
                    var_columns.in_flag.push_back(info_flag_tmp);
                    info_map[INFO.ID[i]] = 0;
                    var_columns.info_map1[INFO.ID[i]] = 0;
                }
            }
        }
        
        var_columns.in_flag.resize(INFO.flags);
        var_columns.in_int.resize(INFO.ints);
        var_columns.in_float.resize(INFO.floats);
        var_columns.in_string.resize(INFO.strings);

        

        
    }
    void print_info_map(){
        for(const auto& element : info_map){
            cout<<element.first<<": "<<element.second<<endl;
        }
    }
    void print_info(){
        cout<<"Flags size: "<<var_columns.in_flag.size()<<endl;
        for(int i=0; i<var_columns.in_flag.size(); i++){
            cout<<var_columns.in_flag[i].name<<": ";
            for(int j=0; j<10; j++){
                cout<<var_columns.in_flag[i].i_flag[j]<<" ";
            }
            cout<<" size: "<<var_columns.in_flag[i].i_flag.size();
            cout<<endl;
        }
        cout<<endl;
        cout<<"Floats size: "<<var_columns.in_float.size()<<endl;
        for(int i=0; i<var_columns.in_float.size(); i++){
            cout<<var_columns.in_float[i].name<<": ";
            for(int j=0; j<10; j++){
                cout<<var_columns.in_float[i].i_float[j]<<" ";
            }
            cout<<" size: "<<var_columns.in_float[i].i_float.size();
            cout<<endl;
        }
        cout<<endl;
        cout<<"Strings size: "<<var_columns.in_string.size()<<endl;
        for(int i=0; i<var_columns.in_string.size(); i++){
            cout<<var_columns.in_string[i].name<<": ";
            for(int j=0; j<10; j++){
                cout<<var_columns.in_string[i].i_string[j]<<" ";
            }
            cout<<" size: "<<var_columns.in_string[i].i_string.size();
            cout<<endl;
        }
        cout<<endl;
        cout<<"Ints size: "<<var_columns.in_int.size()<<endl;
        for(int i=0; i<var_columns.in_int.size(); i++){
            cout<<var_columns.in_int[i].name<<": ";
            for(int j=0; j<10; j++){
                cout<<var_columns.in_int[i].i_int[j]<<" ";
            }
            cout<<" size: "<<var_columns.in_int[i].i_int.size();
            cout<<endl;
        }
        cout<<endl;
        // cout<<"Flags: "<<endl;
        // for(int i=0; i<var_columns.in_flag.size(); i++){
        //     cout<<var_columns.in_flag[i].name<<": ";
        //     for(int j=0; j<var_columns.in_flag[i].i_flag.size(); j++){
        //         cout<<var_columns.in_flag[i].i_flag[j]<<" ";
        //     }
        //     cout<<endl;
        // }
        // cout<<endl;
        // cout<<"Floats: "<<endl;
        // for(int i=0; i<var_columns.in_float.size(); i++){
        //     cout<<var_columns.in_float[i].name<<": ";
        //     for(int j=0; j<var_columns.in_float[i].i_float.size(); j++){
        //         cout<<var_columns.in_float[i].i_float[j]<<" ";
        //     }
        //     cout<<endl;
        // }
        // cout<<endl;
        // cout<<"Strings: "<<endl;
        // for(int i=0; i<var_columns.in_string.size(); i++){
        //     cout<<var_columns.in_string[i].name<<": ";
        //     for(int j=0; j<var_columns.in_string[i].i_string.size(); j++){
        //         cout<<var_columns.in_string[i].i_string[j]<<" ";
        //     }
        //     cout<<endl;
        // }
        // cout<<endl;
        // cout<<"Ints: "<<endl;
        // for(int i=0; i<var_columns.in_int.size(); i++){
        //     cout<<var_columns.in_int[i].name<<": ";
        //     for(int j=0; j<var_columns.in_int[i].i_int.size(); j++){
        //         cout<<var_columns.in_int[i].i_int[j]<<" ";
        //     }
        //     cout<<endl;
        // }
        // cout<<endl;
    }
    void reserve_var_columns(){
        var_columns.var_number.resize(num_lines-1);
        var_columns.chrom.resize(num_lines-1);
        var_columns.id.resize(num_lines-1);
        var_columns.pos.resize(num_lines-1);
        var_columns.ref.resize(num_lines-1);
        var_columns.alt.resize(num_lines-1);
        var_columns.qual.resize(num_lines-1);
        var_columns.filter.resize(num_lines-1);
        var_columns.info.resize(num_lines-1);
        cout<<"Finish resize!"<<endl;
    }
    void populate_var_columns(int num_threads){
        long batch_size = (num_lines-2+num_threads)/num_threads;

#pragma omp parallel
        {
            long start, end;
            int th_ID = omp_get_thread_num();
            //cout << "\nThread id: "<<th_ID<<endl;

            start = th_ID*batch_size; // inizio del batch dello specifico thread
            end = start + batch_size; // fine del batch
            cout<<"\nNum Lines: "<<num_lines-2<<endl;
            for(long i=start; i<end && i<num_lines-1; i++){ //start e end mi dicono l'intervallo di linee che eseguirà ogni thread, quindi la i rappresenta l'iesima linea (var) che inizia a new_lines_index[i] e finisce a new_lines_index[i+1] (escluso)
                var_columns.var_number[i] = i;
                //cout<<"Var_number[i]: "<<var_columns.var_number[i]<<endl;
                var_columns.get_vcf_line_in_var_columns(filestring, new_lines_index[i], new_lines_index[i+1], i);
                //cout<<"enf"<<endl;
            }

        }
    }






    void populate_var_struct(int num_threads){
        auto before = chrono::system_clock::now();
        
        var_df = (var*)calloc((num_lines-1), sizeof(var)); // allocating var_df    
        //cout << "\nBegin tmp: \n" <<"newlines: "<<num_lines<<" num threads: "<<num_threads<<endl;
        long batch_size = (num_lines-2+num_threads)/num_threads; //numero di lines che verrà processato da ogni thread
        
        //cout << "\nBatch size: "<<batch_size<<endl;
        auto after = chrono::system_clock::now();
        auto pre_pragma = std::chrono::duration<double>(after - before).count();
        //cout << "Pre_pragma: " << pre_pragma << " s" << endl;

#pragma omp parallel
        {
            long start, end;
            int th_ID = omp_get_thread_num();
            //cout << "\nThread id: "<<th_ID<<endl;

            start = th_ID*batch_size; // inizio del batch dello specifico thread
            end = start + batch_size; // fine del batch

            //cout << "\nstart: " << start << " end: " << end << endl; //segfault

            for(long i=start; i<end && i<num_lines-1; i++){ //start e end mi dicono l'intervallo di linee che eseguirà ogni thread, quindi la i rappresenta l'iesima linea (var) che inizia a new_lines_index[i] e finisce a new_lines_index[i+1] (escluso)
                var_df[i].get_vcf_line(filestring, new_lines_index[i], new_lines_index[i+1]); //qui traduco da char* a var structure la specifica line
                var_df[i].var_number = i; // dato aggiuntivo come se fosse un altro ID interno 
            }

        }
    }

};

#endif