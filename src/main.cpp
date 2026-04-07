#include <iostream>
#include <algorithm>
#include <chrono>
#include <sys/resource.h>

#include "CSVParser.h"
#include "VCFDataFrames.h"
#include "VCFReconstructor.h"

int main() {
    // Initialize the parser with the input file paths
    CSVParser parser(
        "data/IRBT2/df1.csv", 
        "data/IRBT2/df2.csv",    
        "data/IRBT2/df3.csv", 
        "data/IRBT2/df4.csv",
        "data/IRBT2/irbt_header.txt"
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

        struct rusage usage_before;
        getrusage(RUSAGE_SELF, &usage_before);

        auto start = std::chrono::high_resolution_clock::now();

        reconstructor.run(df1, df2, df3, df4);

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        struct rusage usage_after;
        getrusage(RUSAGE_SELF, &usage_after);

        std::cout << "VCF file reconstructed successfully in " << duration.count() << " ms." << std::endl;
        std::cout <<"Peak RSS: " << usage_after.ru_maxrss << " KB" << std::endl;
        std::cout << "Memory delta: " << (usage_after.ru_maxrss - usage_before.ru_maxrss) << " KB" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Fatal error in Reconstructor: " << e.what() << std::endl;
    }

    return 0;
}