#pragma once
#include "backend.h"
#include "vcf_oracle.h"
#ifdef CUVCF_TEST_GPU
#include "GPUVersion/Reconstructor/VCFReconstructorGPU.h"
#else
#include "CPUVersion/Reconstructor/VCFReconstructor.h"
#endif
namespace test {
inline std::string metadata(const std::string& input) { return input.substr(0,input.find("#CHROM\t")); }
inline void reconstruct(const Frames& frames, const std::string& input,
                        const std::string& output, int threads=1) {
#ifdef CUVCF_TEST_GPU
    VCFReconstructorGPU writer(output,metadata(input));
#else
    VCFReconstructor writer(output,metadata(input));
#endif
    writer.run(frames.d1,frames.d2,frames.d3,frames.d4,threads);
}
inline std::string read_file(const std::string& path) {
    std::ifstream in(path); if(!in) throw std::runtime_error("missing output "+path);
    return std::string(std::istreambuf_iterator<char>(in),{});
}
inline Frames frames_of(const vcf_parsed& p) { return {p.var_columns,p.alt_columns,p.samp_columns,p.alt_sample}; }
inline void add_reconstruction(std::vector<Case>& cases,const std::vector<Fixture>& data) {
    for(const auto& f:data) if(f.roundtrip) {
        cases.push_back({"write/"+f.id,[f] {
            TempDir tmp; auto path=tmp.file("output.vcf").string();
            reconstruct(f.expected,f.input,path); compare_vcf(read_file(path),f.input);
        }});
        cases.push_back({"roundtrip/"+f.id,[f] {
            TempDir tmp; vcf_parsed p{}; parse(p,tmp.write("input.vcf",f.input));
            // Keep parser assertions in parse/*; the independent text oracle checks the composed flow.
            auto path=tmp.file("output.vcf").string(); reconstruct(frames_of(p),f.input,path);
            compare_vcf(read_file(path),f.input);
        }});
    }

    cases.push_back({"write/missing/empty_string_cell",[] {
        auto f=base({"Z"}); strings(f.d3.samp_string,"XS",{""}); f.d3.samp_string[0].numb=1;
        auto input=header(declaration("FORMAT","XS","1","String"),"Z")+row(".","XS\t.");
        TempDir tmp; auto path=tmp.file("out.vcf").string(); reconstruct(f,input,path); compare_vcf(read_file(path),input);
    }});
    const auto base_case=data.front();
    cases.push_back({"write/edited_consistent_values",[base_case] {
        auto f=base_case.expected; f.d1.pos[0]=42; f.d2.alt[0]="T";
        auto expected=base_case.input; auto pos=expected.find("1\t10\trs1\tA\tC,G");
        expected.replace(pos,std::string("1\t10\trs1\tA\tC,G").size(),"1\t42\trs1\tA\tT,G");
        TempDir tmp; auto path=tmp.file("out.vcf").string(); reconstruct(f,expected,path); compare_vcf(read_file(path),expected);
    }});
    cases.push_back({"write/missing_output_directory",[base_case] {
        TempDir tmp; bool thrown=false;
        try { reconstruct(base_case.expected,base_case.input,tmp.file("absent/out.vcf").string()); } catch(const std::runtime_error&) { thrown=true; }
        require(thrown,"output_open","expected an output-open exception");
    }});
    cases.push_back({"write/header_only",[] {
        Frames f{}; auto input=header(); TempDir tmp; auto path=tmp.file("empty.vcf").string();
        reconstruct(f,input,path); compare_vcf(read_file(path),input);
    }});
    cases.push_back({"write/extended/two_threads",[base_case] {
        TempDir tmp; auto path=tmp.file("out.vcf").string(); reconstruct(base_case.expected,base_case.input,path,2); compare_vcf(read_file(path),base_case.input);
    },true});
    cases.push_back({"roundtrip/extended/chunk_boundary",[] {
        // GPU writer chunks at 10,000 records. Keep this outside the pre-commit suite.
        std::string input=header();
        for(int i=1;i<=10001;++i) input+="1\t"+std::to_string(i)+"\t.\tA\tC\t30\tPASS\t.\n";
        TempDir tmp; vcf_parsed p{}; parse(p,tmp.write("in.vcf",input),2);
        auto path=tmp.file("out.vcf").string(); reconstruct(frames_of(p),input,path,2); compare_vcf(read_file(path),input);
    },true});
}
}
