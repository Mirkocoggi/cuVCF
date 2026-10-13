#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <numeric>
#include <sys/resource.h>
#include <omp.h>
#include "CSVParser.h"
#include "VCFDataFrames.h"
#include "VCFReconstructor.h"

inline long getPeakRSS() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
}

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
    //Phase 0

    CSVParser parser(
        "data/IRBT_LONG/df1.csv",
        "data/IRBT_LONG/df2.csv",
        "data/IRBT_LONG/df3.csv",
        "data/IRBT_LONG/df4.csv",
        "data/IRBT/irbt_header.txt"
    );

    var_columns_df df1;
    alt_columns_df df2;
    sample_columns_df df3;
    alt_format_df df4;

    parser.parseMaps("data/IRBT/maps_used_IRBT.csv", df1, df3, df4);
    parser.loadAll(df1, df2, df3, df4);

    //Phase 1
    int active_threads = omp_get_max_threads();
    std::cout << "--- Starting VCF Reconstruction Profiling ---\n";
    std::cout << "OpenMP threads (from OMP_NUM_THREADS or default): " << active_threads << "\n";

    const int N_RUNS = 1;
    std::vector<double> run_times_ms;
    run_times_ms.reserve(N_RUNS);

    long mem_before = getPeakRSS();

    for (int i = 0; i < N_RUNS; ++i) {
        auto start_time = std::chrono::high_resolution_clock::now();

        VCFReconstructor reconstructor("build/output_ricostruito.vcf", parser.header_text);

        try {
            reconstructor.run(df1, df2, df3, df4);
        } catch (const std::exception& e) {
            std::cerr << "Reconstruction Error: " << e.what() << std::endl;
            return 1;
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
        run_times_ms.push_back(ms);
    }

    long mem_after = getPeakRSS();

    // PHASE 2: REPORT
    TimeStats stats = computeStats(run_times_ms);

    std::cout << "\n[Config]\n";
    std::cout << "  - Threads used : " << active_threads << "\n";
    std::cout << "  - N_RUNS       : " << N_RUNS << "\n";

    std::cout << "\n[Time] (over " << N_RUNS << " runs)\n";
    std::cout << "  - Mean : " << stats.mean << " ms\n";
    std::cout << "  - Min  : " << stats.min  << " ms\n";
    std::cout << "  - Max  : " << stats.max  << " ms\n";

    std::cout << "\n[Peak Memory]\n";
    std::cout << "  - Baseline (Pre-Reconstructor) : " << mem_before << " KB\n";
    std::cout << "  - Peak (Post-Reconstructor)    : " << mem_after << " KB\n";
    std::cout << "  - Reconstructor RAM Delta      : +" << (mem_after - mem_before) << " KB\n";

    std::cout << "---------------------------------------------\n";

    return 0;
}