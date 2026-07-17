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
        "data/IRBT3/df1.csv",
        "data/IRBT3/df2.csv",
        "data/IRBT3/df3.csv",
        "data/IRBT3/df4.csv",
        "data/IRBT/irbt_header.txt"
    );

    parser.parseMaps("data/IRBT/maps_used_IRBT.csv", df1, df3, df4);
    parser.loadAll(df1, df2, df3, df4);

    const int N_RUNS = 30;
    long long total_wall_ms = 0;
    TimingResult total_timing;

    for (int i = 0; i < N_RUNS; ++i) {
        auto start_time = std::chrono::high_resolution_clock::now();

        VCFReconstructorGPU reconstructor("build/output_gpu.vcf", parser.header_text);

        TimingResult t;
        try {
            t = reconstructor.run(df1, df2, df3, df4);
        } catch (const std::exception& e) {
            std::cerr << "Reconstruction Error: " << e.what() << std::endl;
            return 1;
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        total_wall_ms    += std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        total_timing.setup_ms  += t.setup_ms;
        total_timing.alloc_ms  += t.alloc_ms;
        total_timing.u2d_ms    += t.u2d_ms;
        total_timing.prep_ms   += t.prep_ms;
        total_timing.kernel_ms += t.kernel_ms;
        total_timing.write_ms  += t.write_ms;
        total_timing.free_ms   += t.free_ms;
        total_timing.drain_ms  += t.drain_ms;
    }

    std::cout << "\n[Time] (avg over " << N_RUNS << " runs)\n";
    std::cout << "  - Total Execution : " << (total_wall_ms / N_RUNS) << " ms\n";
    std::cout << "  - setup           : " << (total_timing.setup_ms  / N_RUNS) << " ms\n";
    std::cout << "  - alloc           : " << (total_timing.alloc_ms  / N_RUNS) << " ms\n";
    std::cout << "  - u2d             : " << (total_timing.u2d_ms    / N_RUNS) << " ms\n";
    std::cout << "  - prep            : " << (total_timing.prep_ms   / N_RUNS) << " ms\n";
    std::cout << "  - kernel          : " << (total_timing.kernel_ms / N_RUNS) << " ms\n";
    std::cout << "  - write           : " << (total_timing.write_ms  / N_RUNS) << " ms\n";
    std::cout << "  - free            : " << (total_timing.free_ms   / N_RUNS) << " ms\n";
    std::cout << "  - drain           : " << (total_timing.drain_ms  / N_RUNS) << " ms\n";

    return 0;
}