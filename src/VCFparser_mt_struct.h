#ifndef VCF_STRUCTS_H
#define VCF_STRUCTS_H
#include <chrono>
#include <boost/algorithm/string.hpp>


using namespace std;

class var
{
public:
    int var_number;
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
    void get_vcf_line(char *line, int start, int end)
    {
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
                qual = stof(tmp); // da cambiare, in futuro
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
                info += line[start+iter]; // anche qui andrebbero separate le info e creati dizionari con keys e values
                iter++;
            }
        }
        find1=false;
        while(!find1){
            if(line[start+iter]=='\t'||line[start+iter]==' '){
                find1 = true;
                iter++;
            }else{
                format += line[start+iter]; // vanno separati e soprattutti connessi ai valori nei samples (tipo se nel format ho ..:AF e nel sample ..:2, devo poter fare la query AF=2 in tutte le var)
                iter++;
            }
        }
        
        while((start+iter)!=(end-1)){
            samples += line[start+iter]; // sicuramente da dividere i vari sample e eventualmente creare un dizionario per ogni sample cosi da associare con il format
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
    int header_size=0;
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
            header_size += line.length() + 1;
        }
        header_size += line.length() + 1;
        cout << "\nheader char: " << to_string(header_size) << endl;
        variants_size = filesize - header_size; //ora che ho tolto l'header ho un file piu piccolo quindi una nuova size
        cout<<"filesize: "<<filesize<<" variants_size: "<<variants_size<<endl;
    }
    void print_header(){
        cout << "VCF header:\n" << header << endl;
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
            cout << "in find start: "<<start<<" end: "<<end<<endl;
            for(long i=start; i<end && i<variants_size; i++){
                filestring[i] = infile.get(); //file_streams[thr_ID].get(); //-----------------
            }
        }
        auto after = chrono::system_clock::now();
        auto filestring_time = std::chrono::duration<double>(after - before).count();
        before = chrono::system_clock
        ::now();
        num_char = 0;
        while(num_char!=variants_size){ //conto il numero di lines e ogni volta che trovo \n salvo il char successivo come inizio della riga successiva
            if(filestring[num_char]=='\n'){
                new_lines_index[num_lines] = num_char+1; // PROBLEMA: funziona solo se l'ultimo char è uno /n, altrimenti si rompe
                num_lines++;
            }
            num_char++;
        } // si potrebbe fare multithreading ma vediamo se si ha beneficio
        after = chrono::system_clock::now();

        auto f_new_lines = std::chrono::duration<double>(after - before).count();
        cout << "\nFilestring time: " << filestring_time << " s " << "New lines time: " << f_new_lines << " s\n\n";
    }
    void populate_var_struct(int num_threads){
        auto before = chrono::system_clock::now();
        
        //var_df = (var*)malloc((num_lines-1) * sizeof(var)); // Usando la malloc migliora in alcuni punti peggiora in altri, il bug era il \0
        var_df = (var*)calloc((num_lines-1), sizeof(var)); // allocating var_df
        cout << "\nBegin tmp: \n" <<"newlines: "<<num_lines<<" num threads: "<<num_threads<<endl;
        long batch_size = (num_lines-2+num_threads)/num_threads; //numero di lines che verrà processato da ogni thread
        
        cout << "\nBatch size: "<<batch_size<<endl;
        auto after = chrono::system_clock::now();
        auto pre_pragma = std::chrono::duration<double>(after - before).count();
        cout << "Pre_pragma: " << pre_pragma << " s" << endl;

#pragma omp parallel
        {
            long start, end;
            int th_ID = omp_get_thread_num();
            cout << "\nThread id: "<<th_ID<<endl;

            start = th_ID*batch_size; // inizio del batch dello specifico thread
            end = start + batch_size; // fine del batch

            cout << "\nstart: " << start << " end: " << end << endl;

            for(long i=start; i<end && i<num_lines-1; i++){ //start e end mi dicono l'intervallo di linee che eseguirà ogni thread, quindi la i rappresenta l'iesima linea (var) che inizia a new_lines_index[i] e finisce a new_lines_index[i+1] (escluso)

                var_df[i].get_vcf_line(filestring, new_lines_index[i], new_lines_index[i+1]); //qui traduco da char* a var structure la specifica line
                var_df[i].var_number = i; // dato aggiuntivo come se fosse un altro ID interno 
            
            }
        
        
        }
    }

};

#endif