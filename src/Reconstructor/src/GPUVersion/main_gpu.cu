#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <numeric>
#include "../CPUVersion/CSVParser.h"
#include "VCFReconstructorGPU.h"

struct TimeStats {
    double mean, min, max;
};

TimeStats computeStats(const std::vector<double>& times) {
    TimeStats s{};
    double sum = std::accumulate(times.begin(), times.end(), 0.0);
    s.mean = sum / times.size();
    s.min = *std::min_element(times.begin(), times.end());
    s.max = *std::max_element(times.begin(), times.end());
    return s;
}

int main() {
    var_columns_df df1;
    alt_columns_df df2;
    sample_columns_df df3;
    alt_format_df df4;

    CSVParser parser(
        "data/IRBT_LONG/df1.csv",
        "data/IRBT_LONG/df2.csv",
        "data/IRBT_LONG/df3.csv",
        "data/IRBT_LONG/df4.csv",
        "data/IRBT/irbt_header.txt"
    );

    parser.parseMaps("data/IRBT/maps_used_IRBT.csv", df1, df3, df4);
    parser.loadAll(df1, df2, df3, df4);

    const int N_RUNS = 3;
    std::vector<double> run_times_ms;
    run_times_ms.reserve(N_RUNS);
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
        run_times_ms.push_back(std::chrono::duration<double, std::milli>(end_time - start_time).count());
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
    TimeStats stats = computeStats(run_times_ms);

    std::cout << "\n[Time] (over " << N_RUNS << " runs)\n";
    std::cout << "  - Mean : " << stats.mean << " ms\n";
    std::cout << "  - Min  : " << stats.min  << " ms\n";
    std::cout << "  - Max  : " << stats.max  << " ms\n";
    std::cout << "\n[Phases] (avg over " << N_RUNS << " runs)\n";
    std::cout << "  - setup           : " << (total_timing.setup_ms / N_RUNS) << " ms\n";
    std::cout << "  - prep (host)     : " << (total_timing.prep_ms  / N_RUNS) << " ms\n";
    std::cout << "  - loop (pipeline) : " << (total_timing.loop_ms  / N_RUNS) << " ms\n";
    std::cout << "  - drain           : " << (total_timing.drain_ms / N_RUNS) << " ms\n";

    return 0;
}