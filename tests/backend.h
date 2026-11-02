#pragma once
#include <iostream>
#include <map>
#include <algorithm>
#include <iomanip>
#include <type_traits>
#ifdef CUVCF_TEST_GPU
#include "GPUVersion/Parser.h"
#else
#include "CPUVersion/VCF_parsed.h"
#endif
#include "harness.h"
#include <cmath>
#include <limits>
#include <zlib.h>

namespace test {
inline void parse(vcf_parsed& parsed, std::string path, int threads = 1) {
    parsed.run(path.data(), threads); // Both implementations may strip .gz in place.
}
inline float number(float x) { return x; }
#ifdef CUVCF_TEST_GPU
inline float number(__half x) { return __half2float(x); }
inline __half half_value(float x) { return __float2half(x); }
#else
inline float number(half x) { return static_cast<float>(x); }
inline half half_value(float x) { return half(x); }
#endif
// Semantic expectations must distinguish a missing value from zero or a valid -1.
inline int missing_int() { return std::numeric_limits<int>::min(); }
inline float missing_float() { return std::numeric_limits<float>::quiet_NaN(); }
inline std::string declaration(const std::string& scope, const std::string& id,
                               const std::string& n, const std::string& type) {
    return "##" + scope + "=<ID=" + id + ",Number=" + n + ",Type=" + type + ",Description=\"test\">\n";
}
inline std::string header(const std::string& meta = "", const std::string& samples = "") {
    return "##fileformat=VCFv4.3\n##contig=<ID=1>\n##contig=<ID=2>\n"
           "##FILTER=<ID=q10,Description=\"low quality\">\n" + meta +
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" +
           (samples.empty() ? "\n" : "\tFORMAT\t" + samples + "\n");
}
struct Frames {
    var_columns_df d1{};
    alt_columns_df d2{};
    sample_columns_df d3{};
    alt_format_df d4{};
};
inline Frames base(const std::vector<std::string>& samples = {}) {
    Frames f;
    f.d1.var_number = {0}; f.d1.chrom_map = {{"1", 0}}; f.d1.chrom = {0};
    f.d1.pos = {10}; f.d1.id = {"rs1"}; f.d1.ref = {"A"};
    f.d1.qual = {half_value(30)}; f.d1.filter_map = {{"PASS", 0}}; f.d1.filter = {0};
    f.d2.var_id = {0, 0}; f.d2.alt_id = {0, 1}; f.d2.alt = {"C", "G"}; f.d2.numAlt = 2;
    f.d3.numSample = f.d4.numSample = static_cast<int>(samples.size());
    for (size_t i = 0; i < samples.size(); ++i) {
        f.d3.sampNames[samples[i]] = f.d4.sampNames[samples[i]] = i;
        f.d3.var_id.push_back(0); f.d3.samp_id.push_back(i);
    }
    f.d3.initMapGT();
    return f;
}
inline std::string row(const std::string& info = ".", const std::string& sample = "") {
    return "1\t10\trs1\tA\tC,G\t30\tPASS\t" + info + (sample.empty() ? "\n" : "\t" + sample + "\n");
}
inline void alleles(Frames& f) {
    for (int s = 0; s < f.d3.numSample; ++s) {
        f.d4.var_id.insert(f.d4.var_id.end(), {0, 0});
        f.d4.samp_id.insert(f.d4.samp_id.end(), {static_cast<unsigned short>(s), static_cast<unsigned short>(s)});
        f.d4.alt_id.insert(f.d4.alt_id.end(), {0, 1});
    }
}
template<class V> void ints(V& cols, const std::string& name, std::initializer_list<int> values) {
    typename V::value_type c{}; c.name = name; c.i_int = values; cols.push_back(c);
}
template<class V> void floats(V& cols, const std::string& name, std::initializer_list<float> values) {
    typename V::value_type c{}; c.name = name;
    for (auto v : values) c.i_float.push_back(half_value(v));
    cols.push_back(c);
}
template<class V> void strings(V& cols, const std::string& name, std::initializer_list<std::string> values) {
    typename V::value_type c{}; c.name = name; c.i_string = values; cols.push_back(c);
}
template<class A, class B> void sequence(const A& actual, const B& expected, const std::string& key) {
    expect_equal(actual.size(), expected.size(), key + "/size");
    for (size_t i = 0; i < expected.size(); ++i) {
        // Numeric promotion ensures char codes print as numbers rather than control characters.
        if constexpr (std::is_arithmetic_v<typename A::value_type>)
            expect_equal(+actual.at(i), +expected.at(i), key + "/" + std::to_string(i));
        else expect_equal(actual.at(i), expected.at(i), key + "/" + std::to_string(i));
    }
}
template<class A, class B> void float_sequence(const A& actual, const B& expected, const std::string& key) {
    expect_equal(actual.size(), expected.size(), key + "/size");
    for (size_t i = 0; i < expected.size(); ++i) {
        auto a = number(actual.at(i)), e = number(expected.at(i));
        if (std::isnan(e)) require(std::isnan(a), key + "/" + std::to_string(i), "expected missing (NaN), got " + std::to_string(a));
        else expect_equal(a, e, key + "/" + std::to_string(i)); // Values are already quantized to half.
    }
}
template<class V, class Compare> void columns(const V& a, const V& e, const std::string& key, Compare compare) {
    expect_equal(a.size(), e.size(), key + "/columns");
    for (const auto& col : e) {
        auto it = std::find_if(a.begin(), a.end(), [&](const auto& v) { return v.name == col.name; });
        require(it != a.end(), key + "/" + col.name, "column missing");
        compare(*it, col, key + "/" + col.name);
    }
}
inline void check_frames(const vcf_parsed& p, const Frames& f) {
    const auto& a = p.var_columns; const auto& b = p.alt_columns;
    const auto& c = p.samp_columns; const auto& d = p.alt_sample;
#define SEQ(x,y) sequence(x, y, #x)
    SEQ(a.var_number,f.d1.var_number); SEQ(a.pos,f.d1.pos); SEQ(a.id,f.d1.id); SEQ(a.ref,f.d1.ref);
    SEQ(a.chrom,f.d1.chrom); SEQ(a.filter,f.d1.filter);
    require(a.chrom_map == f.d1.chrom_map, "chrom_map", "unexpected chromosome mapping");
    require(a.filter_map == f.d1.filter_map, "filter_map", "unexpected filter mapping");
    float_sequence(a.qual,f.d1.qual,"qual");
    SEQ(b.var_id,f.d2.var_id); SEQ(b.alt_id,f.d2.alt_id); SEQ(b.alt,f.d2.alt);
    SEQ(c.var_id,f.d3.var_id); SEQ(c.samp_id,f.d3.samp_id);
    SEQ(d.var_id,f.d4.var_id); SEQ(d.samp_id,f.d4.samp_id); SEQ(d.alt_id,f.d4.alt_id);
    require(c.sampNames == f.d3.sampNames && d.sampNames == f.d4.sampNames, "samples", "sample order/name mismatch");
    auto ci = [](const auto& x,const auto& y,const auto& k){ sequence(x.i_int,y.i_int,k); };
    auto cf = [](const auto& x,const auto& y,const auto& k){ float_sequence(x.i_float,y.i_float,k); };
    auto cs = [](const auto& x,const auto& y,const auto& k){ sequence(x.i_string,y.i_string,k); };
    auto cb = [](const auto& x,const auto& y,const auto& k){ sequence(x.i_flag,y.i_flag,k); };
#define COL(x,y,cmp) columns(x,y,#x,cmp)
    COL(a.in_int,f.d1.in_int,ci); COL(a.in_float,f.d1.in_float,cf); COL(a.in_string,f.d1.in_string,cs); COL(a.in_flag,f.d1.in_flag,cb);
    COL(b.alt_int,f.d2.alt_int,ci); COL(b.alt_float,f.d2.alt_float,cf); COL(b.alt_string,f.d2.alt_string,cs); COL(b.alt_flag,f.d2.alt_flag,cb);
    COL(c.samp_int,f.d3.samp_int,ci); COL(c.samp_float,f.d3.samp_float,cf); COL(c.samp_string,f.d3.samp_string,cs); COL(c.samp_flag,f.d3.samp_flag,cb);
    COL(d.samp_int,f.d4.samp_int,ci); COL(d.samp_float,f.d4.samp_float,cf); COL(d.samp_string,f.d4.samp_string,cs); COL(d.samp_flag,f.d4.samp_flag,cb);
    expect_equal(c.sample_GT.size(),f.d3.sample_GT.size(),"GT/columns");
    for(size_t i=0;i<c.sample_GT.size();++i) sequence(c.sample_GT[i].GT,f.d3.sample_GT[i].GT,"GT/"+std::to_string(i));
    SEQ(d.sample_GT.GT,f.d4.sample_GT.GT);
#undef SEQ
#undef COL
}
struct Fixture {
    std::string id, input;
    Frames expected;
    bool roundtrip = true;
};
}
