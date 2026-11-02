#include "cases/fixtures.h"
#include "cases/reconstruction.h"
#include "support/expected_failures.h"

int main(int argc, char** argv) {
    using namespace test;
    std::vector<Case> cases;
#ifdef CUVCF_TEST_GPU
    // Listing is allowed without a device; executing GPU tests must never silently skip.
    bool list = false; for(int i=1;i<argc;++i) if(std::string(argv[i])=="--list") list=true;
    if(!list) {
        int count=0; auto error=cudaGetDeviceCount(&count);
        if(error!=cudaSuccess || count==0) { std::cerr<<"GPU suite requires an NVIDIA device: "<<cudaGetErrorString(error)<<'\n'; return 2; }
    }
#endif
    const auto data=fixtures();
    for(const auto& fixture:data) {
        cases.push_back({"parse/"+fixture.id,[fixture] {
            TempDir tmp; vcf_parsed parsed{};
            parse(parsed,tmp.write("input.vcf",fixture.input));
            expect_equal(parsed.header,fixture.input.substr(0,fixture.input.find("#CHROM\t")),"header/preservation");
            check_frames(parsed,fixture.expected);
        }});
    }
    auto baseline=data.front();
    cases.push_back({"parse/input/no_final_newline",[baseline] {
        TempDir tmp; vcf_parsed p{}; auto input=baseline.input; input.pop_back();
        parse(p,tmp.write("input.vcf",input)); check_frames(p,baseline.expected);
    }});
    cases.push_back({"parse/input/gzip",[baseline] {
        TempDir tmp; auto path=tmp.file("input.vcf.gz").string();
        gzFile file=gzopen(path.c_str(),"wb"); if(!file) throw std::runtime_error("gzip fixture open failed");
        auto wrote=gzwrite(file,baseline.input.data(),static_cast<unsigned int>(baseline.input.size()));
        auto closed=gzclose(file);
        if(wrote!=static_cast<int>(baseline.input.size()) || closed!=Z_OK) throw std::runtime_error("gzip fixture write failed");
        vcf_parsed p{}; parse(p,path); check_frames(p,baseline.expected);
    }});
    cases.push_back({"parse/input/header_only",[] {
        TempDir tmp; vcf_parsed p{}; parse(p,tmp.write("input.vcf",header()));
        require(p.var_columns.pos.empty() && p.alt_columns.var_id.empty() && p.samp_columns.var_id.empty() && p.alt_sample.var_id.empty(),"empty_rows","header-only input created rows");
    }});
    cases.push_back({"parse/input/missing_file",[] {
        TempDir tmp; vcf_parsed p{}; bool thrown=false;
        try { parse(p,tmp.file("absent.vcf").string()); } catch(const std::runtime_error&) { thrown=true; }
        require(thrown,"missing_file","expected an input-open exception");
    }});
    cases.push_back({"parse/extended/two_threads",[baseline] {
        TempDir tmp; vcf_parsed p{}; parse(p,tmp.write("input.vcf",baseline.input),2); check_frames(p,baseline.expected);
    },true});
    add_reconstruction(cases,data);
    try { apply_expected_failures(cases); }
    catch(const std::exception& e) { std::cerr << "Test setup failed: " << e.what() << '\n'; return 2; }
    return run(cases,argc,argv);
}
