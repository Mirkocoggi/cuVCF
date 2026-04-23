#include <iostream>
#include <chrono>
#include <sys/resource.h>
#include "CSVParser.h" 
#include "VCFDataFrames.h"
#include "VCFReconstructor.h"

inline long getPeakRSS() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss; // Note: KB on Linux, Bytes on macOS
}

int main() {
    // ==========================================
    // PHASE 0: MOCK SETUP (DFs already in RAM)
    // ==========================================
    // You will remove this parsing section during integration.
    // It is only here to populate the DataFrames for the standalone test.
    CSVParser parser(
        "data/bos_taurus/df1.csv",
        "data/bos_taurus/df2.csv",
        "data/bos_taurus/df3.csv",
        "data/bos_taurus/df4.csv",
        "data/bos_taurus/bos_taurus_header.txt"
    );

    var_columns_df df1;
    alt_columns_df df2;
    sample_columns_df df3;
    alt_format_df df4;

    parser.loadAll(df1, df2, df3, df4);
    
    // ==========================================
    // PHASE 1: VCF RECONSTRUCTION PROFILING
    // ==========================================
    std::cout << "--- Starting VCF Reconstruction Profiling ---\n";
    
    // Track baseline memory and start time right before the reconstructor comes into play
    long mem_before = getPeakRSS();
    auto start_time = std::chrono::high_resolution_clock::now();
    
    VCFReconstructor reconstructor("build/output_ricostruito.vcf", parser.header_text);

    try {
        reconstructor.run(df1, df2, df3, df4);
    } catch (const std::exception& e) {
        std::cerr << "Reconstruction Error: " << e.what() << std::endl;
        return 1;
    }
    
    // Stop tracking right after the run finishes
    auto end_time = std::chrono::high_resolution_clock::now();
    long mem_after = getPeakRSS();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    // ==========================================
    // PHASE 2: REPORT
    // ==========================================
    std::cout << "\n[Time]\n";
    std::cout << "  - Total Execution (Init + Run) : " << duration.count() << " ms\n";
    
    std::cout << "\n[Peak Memory]\n";
    std::cout << "  - Baseline (Pre-Reconstructor) : " << mem_before << " KB\n";
    std::cout << "  - Peak (Post-Reconstructor)    : " << mem_after << " KB\n";
    std::cout << "  - Reconstructor RAM Delta      : +" << (mem_after - mem_before) << " KB\n";
    std::cout << "---------------------------------------------\n";

    return 0;
}