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
        "data/IRBT/df1.csv",
        "data/IRBT/df2.csv",
        "data/IRBT/df3.csv",
        "data/IRBT/df4.csv",
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
        total_timing.setup_ms += t.setup_ms;
        total_timing.prep_ms  += t.prep_ms;
        total_timing.loop_ms  += t.loop_ms;
        total_timing.drain_ms += t.drain_ms;
    }

    // Per-phase alloc/u2d/kernel/write timings are gone: chunk N+1's prep now
    // genuinely overlaps chunk N's GPU work, so those phases no longer have
    // isolated windows to measure without re-introducing blocking syncs. Use
    // Nsight Systems (NVTX ranges already in run()) to inspect phase-level
    // overlap; prep_ms/loop_ms below give the aggregate picture.
    std::cout << "\n[Time] (avg over " << N_RUNS << " runs)\n";
    std::cout << "  - Total Execution : " << (total_wall_ms / N_RUNS) << " ms\n";
    std::cout << "  - setup           : " << (total_timing.setup_ms / N_RUNS) << " ms\n";
    std::cout << "  - prep (host)     : " << (total_timing.prep_ms  / N_RUNS) << " ms\n";
    std::cout << "  - loop (pipeline) : " << (total_timing.loop_ms  / N_RUNS) << " ms\n";
    std::cout << "  - drain           : " << (total_timing.drain_ms / N_RUNS) << " ms\n";

    return 0;
}