#include <iostream>
#include <algorithm> // per std::min
#include "CSVParser.h"
#include "VCFDataFrames.h"
#include "VCFReconstructor.h"

int main() {
    // Inizializza il parser con i percorsi dei file
    CSVParser parser(
        "data/test1/df1.csv", 
        "data/test1/df2.csv", 
        "data/test1/df3.csv", 
        "data/test1/df4.csv"
    );

    // Crea le strutture dati vuote
    var_columns_df df1;
    alt_columns_df df2;
    sample_columns_df df3;
    alt_format_df df4;

    std::cout << "Caricamento di tutti i DataFrame in corso..." << std::endl;
    
    try {
        parser.loadAll(df1, df2, df3, df4);
        std::cout << "Caricamento completato con successo!\n" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Errore durante il caricamento: " << e.what() << std::endl;
        return 1;
    }

    /*
    // ==========================================
    // VERIFICA DATI DF1 (Variant Core)
    // ==========================================
    std::cout << "--- VERIFICA DATI DF1 (Core Varianti) ---" << std::endl;
    int righe_df1 = std::min(3, (int)df1.var_number.size()); 
    for (int i = 0; i < righe_df1; i++) {
        std::cout << "Riga " << i 
                  << ": VAR_NUM=" << df1.var_number[i]
                  << " CHROM_ID=" << (int)df1.chrom[i]
                  << " POS=" << df1.pos[i] 
                  << " ID=" << df1.id[i]
                  << " REF=" << df1.ref[i]
                  << " FILTER_ID=" << (int)df1.filter[i]
                  << " QUAL=" << df1.qual[i];
        
        if (!df1.in_string.empty()) {
            std::cout << "\n    -> INFO: ";
            for (const auto& campo : df1.in_string) {
                std::cout << campo.name << "=" << campo.i_string[i] << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << "Totale righe DF1: " << df1.var_number.size() << "\n" << std::endl;

    // ==========================================
    // VERIFICA DATI DF2 (Alleli Alternativi)
    // ==========================================
    std::cout << "--- VERIFICA DATI DF2 (Alleli Alternativi) ---" << std::endl;
    int righe_df2 = std::min(3, (int)df2.var_id.size()); 
    for (int i = 0; i < righe_df2; i++) {
        std::cout << "Riga " << i 
                  << ": VAR_ID=" << df2.var_id[i]
                  << " ALT_ID=" << (int)df2.alt_id[i]
                  << " ALT=" << df2.alt[i];
        
        if (!df2.alt_string.empty()) {
            std::cout << "\n    -> INFO ALT: ";
            for (const auto& campo : df2.alt_string) {
                std::cout << campo.name << "=" << campo.i_string[i] << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << "Totale righe DF2: " << df2.var_id.size() << "\n" << std::endl;

    // ==========================================
    // VERIFICA DATI DF3 (Campioni Core)
    // ==========================================
    std::cout << "--- VERIFICA DATI DF3 (Campioni Core) ---" << std::endl;
    int righe_df3 = std::min(3, (int)df3.var_id.size()); 
    for (int i = 0; i < righe_df3; i++) {
        std::cout << "Riga " << i 
                  << ": VAR_ID=" << df3.var_id[i]
                  << " SAMP_ID=" << df3.samp_id[i]
                  // Usiamo la funzione di decodifica per stampare il Genotipo leggibile (es. "0/1")
                  << " GT=" << df3.getGTStringFromChar(df3.sample_GT[0].GT[i]);
        
        if (!df3.samp_string.empty()) {
            std::cout << "\n    -> FORMAT: ";
            for (const auto& campo : df3.samp_string) {
                std::cout << campo.name << "=" << campo.i_string[i] << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << "Totale righe DF3: " << df3.var_id.size() << "\n" << std::endl;

    // ==========================================
    // VERIFICA DATI DF4 (Campioni x Alleli)
    // ==========================================
    std::cout << "--- VERIFICA DATI DF4 (Campioni x Alleli) ---" << std::endl;
    int righe_df4 = std::min(3, (int)df4.var_id.size()); 
    for (int i = 0; i < righe_df4; i++) {
        std::cout << "Riga " << i 
                  << ": VAR_ID=" << df4.var_id[i]
                  << " SAMP_ID=" << df4.samp_id[i]
                  << " ALT_ID=" << (int)df4.alt_id[i]
                  << " GT=" << df4.getGTStringFromChar(df4.sample_GT.GT[i]);
        
        if (!df4.samp_string.empty()) {
            std::cout << "\n    -> FORMAT ALT: ";
            for (const auto& campo : df4.samp_string) {
                std::cout << campo.name << "=" << campo.i_string[i] << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << "Totale righe DF4: " << df4.var_id.size() << std::endl;
    */
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "Avvio ricostruzione file VCF..." << std::endl;
    
    // 1. Creiamo l'oggetto Reconstructor e gli diciamo dove salvare il file
    VCFReconstructor reconstructor("build/output_ricostruito.vcf");
    
    // 2. Facciamo partire la magia passandogli i 4 DataFrame
    try {
        reconstructor.run(df1, df2, df3, df4);
        std::cout << "File VCF ricostruito con successo!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Errore fatale nel Reconstructor: " << e.what() << std::endl;
    }

    return 0;
}