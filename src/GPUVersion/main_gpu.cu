#include <iostream>
#include <chrono>
#include "../CPUVersion/CSVParser.h"
#include "VCFReconstructorGPU.h"

int main() {
    var_columns_df df1;
    alt_columns_df df2;
    sample_columns_df df3;
    alt_format_df df4;

    CSVParser parser(
        "data/IRBT2/df1.csv",
        "data/IRBT2/df2.csv",
        "data/IRBT2/df3.csv",
        "data/IRBT2/df4.csv",
        "data/IRBT/irbt_header.txt"
    );

    parser.parseMaps("data/IRBT/maps_used_IRBT.csv", df1, df3, df4);
    parser.loadAll(df1, df2, df3, df4);

    auto start_time = std::chrono::high_resolution_clock::now();

    VCFReconstructorGPU reconstructor("build/output_gpu.vcf", parser.header_text);

    try {
        reconstructor.run(df1, df2);
    } catch (const std::exception& e) {
        std::cerr << "Reconstruction Error: " << e.what() << std::endl;
        return 1;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "\n[Time]\n";
    std::cout << "  - Total Execution: " << duration.count() << " ms\n";

    return 0;
}