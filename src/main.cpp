#include <iostream>
#include <algorithm>
#include <chrono>
#include "CSVParser.h"
#include "VCFDataFrames.h"
#include "VCFReconstructor.h"

int main() {
    // Initialize the parser with the input file paths
    CSVParser parser(
        "data/IRBT/df1.csv", 
        "data/IRBT/df2.csv",    
        "data/IRBT/df3.csv", 
        "data/IRBT/df4.csv",
        "data/IRBT/irbt_header.txt"
    );

    // Create empty dataframe structures
    var_columns_df df1;
    alt_columns_df df2;
    sample_columns_df df3;
    alt_format_df df4;

    std::cout << "Loading all DataFrames..." << std::endl;
    
    try {
        parser.loadAll(df1, df2, df3, df4);
        std::cout << "Loading completed successfully.\n" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error during loading: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "===========================================" << std::endl;
    std::cout << "Starting VCF reconstruction..." << std::endl;
    
    // Create the reconstructor and set the output path
    VCFReconstructor reconstructor(
        "build/output_ricostruito.vcf",
        parser.header_text,
        parser.field_types,
        parser.samp_names
    );
    
    // Run the reconstruction pipeline
    try {

        auto start = std::chrono::high_resolution_clock::now();

        reconstructor.run(df1, df2, df3, df4);

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);


        std::cout << "VCF file reconstructed successfully in " << duration.count() << " ms." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Fatal error in Reconstructor: " << e.what() << std::endl;
    }

    return 0;
}