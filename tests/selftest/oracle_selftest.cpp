#include "oracle/vcf_oracle.h"
using namespace test;
const std::string meta="##fileformat=VCFv4.3\n##INFO=<ID=I,Number=A,Type=Integer,Description=\"test\">\n##FORMAT=<ID=X,Number=2,Type=Integer,Description=\"test\">\n##FORMAT=<ID=F,Number=1,Type=Float,Description=\"test\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tZ\n";
const std::string row="1\t10\t.\tA\tC,G\t30\tPASS\tI=1,2\tX:F\t.,7:0.1\n";
void rejects(const std::string& bad,const std::string& key) {
    bool rejected=false;
    try { compare_vcf(meta+bad,meta+row); } catch(const Assertion& e) { rejected=true; expect_equal(e.key,key,"rejection/key"); }
    require(rejected,"rejection/missing","corrupt VCF accepted");
}
int main(int argc,char** argv) {
    auto document=[](const std::string& type,const std::string& n,const std::string& info) {
        return "##fileformat=VCFv4.3\n##INFO=<ID=I,Number="+n+",Type="+type+",Description=\"test\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t10\t.\tA\tC,G\t30\tPASS\t"+info+"\n";
    };
    std::vector<Case> cases{
        {"oracle/missing_INFO_equivalence",[document]{
            for(const auto& type:{"Integer","Float","String"}) for(const auto& n:{"1","A"}) {
                auto omitted=document(type,n,"."), explicit_missing=document(type,n,"I=.");
                compare_vcf(explicit_missing,omitted); compare_vcf(omitted,explicit_missing);
                if(std::string(n)=="A") { auto array=document(type,n,"I=.,."); compare_vcf(array,omitted); compare_vcf(omitted,array); }
            }
        }},
        {"oracle/missing_INFO_preserves_values",[document]{
            for(const auto& value:{"I=.,7","I=7,.","I=0,0"}) {
                bool rejected=false;
                try { compare_vcf(document("Integer","A",value),document("Integer","A",".")); }
                catch(const Assertion&) { rejected=true; }
                require(rejected,"missing_INFO/value","nonmissing INFO was discarded");
            }
        }},
        {"oracle/missing_INFO_preserves_validation",[document]{
            for(const auto& value:{"I=.,.,.","I=.;I=.","U=.","I="}) {
                bool rejected=false;
                try { compare_vcf(document("Integer","A",value),document("Integer","A",".")); }
                catch(const Assertion&) { rejected=true; }
                require(rejected,"missing_INFO/validation","invalid INFO accepted: "+std::string(value));
            }
            bool rejected=false;
            try { compare_vcf(document("Flag","0","I"),document("Flag","0",".")); }
            catch(const Assertion&) { rejected=true; }
            require(rejected,"missing_INFO/flag","flag presence was discarded");
        }},
        {"oracle/format_order_and_float",[]{compare_vcf(meta+"1\t10\t.\tA\tC,G\t30.000000\tPASS\tI=1,2\tF:X\t0.099976:.,7\n",meta+row);}},
        {"oracle/duplicate_info",[]{rejects("1\t10\t.\tA\tC,G\t30\tPASS\tI=1,2;I=3,4\tX:F\t.,7:0.1\n","info/duplicate");}},
        {"oracle/leading_colon",[]{rejects("1\t10\t.\tA\tC,G\t30\tPASS\tI=1,2\tX:F\t:.,7:0.1\n","format/width");}},
        {"oracle/empty_value",[]{rejects("1\t10\t.\tA\tC,G\t30\tPASS\tI=1,2\tX:F\t:0.1\n","format/empty");}},
        {"oracle/renamed_field",[]{rejects("1\t10\t.\tA\tC,G\t30\tPASS\tI=1,2\tY:F\t.,7:0.1\n","format/undeclared");}},
        {"oracle/truncated_array",[]{rejects("1\t10\t.\tA\tC,G\t30\tPASS\tI=1,2\tX:F\t7:0.1\n","record/0/sample/0/X/cardinality");}},
        {"oracle/missing_position",[]{rejects("1\t10\t.\tA\tC,G\t30\tPASS\tI=1,2\tX:F\t7,.:0.099976\n","record/0/sample/0/X/0");}},
        {"oracle/changed_float",[]{rejects("1\t10\t.\tA\tC,G\t30\tPASS\tI=1,2\tX:F\t.,7:0.11\n","record/0/sample/0/F/0");}},
        {"oracle/missing_trailing_field",[]{
            compare_vcf(meta+"1\t10\t.\tA\tC,G\t30\tPASS\tI=1,2\tX:F\t.,7\n",meta+"1\t10\t.\tA\tC,G\t30\tPASS\tI=1,2\tX:F\t.,7:.\n");
        }}
    };
    return run(cases,argc,argv);
}
