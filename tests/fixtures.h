#pragma once
#include "backend.h"
namespace test {
inline std::vector<Fixture> fixtures() {
    std::vector<Fixture> out;
    out.push_back({"core/multiallelic", header() + row(), base()});
    for (const auto& type : {std::string("Integer"), std::string("Float"), std::string("String")}) {
        for (const auto& count : {std::string("1"), std::string("A")}) {
            auto f = base(); const bool alt = count == "A";
            std::string value;
            if (type == "Integer") {
                value = alt ? "0,-7" : "-7";
                if (alt) ints(f.d2.alt_int,"X",{0,-7}); else ints(f.d1.in_int,"X",{-7});
            } else if (type == "Float") {
                value = alt ? "0.25,-0.5" : "0.25";
                if (alt) floats(f.d2.alt_float,"X",{0.25f,-0.5f}); else floats(f.d1.in_float,"X",{0.25f});
            } else {
                value = alt ? "alpha,beta" : "alpha";
                if (alt) strings(f.d2.alt_string,"X",{"alpha","beta"}); else strings(f.d1.in_string,"X",{"alpha"});
            }
            out.push_back({"info/"+type+"/"+count,header(declaration("INFO","X",count,type))+row("X="+value),f});
        }
        for (const auto& count : {std::string("1"),std::string("3"),std::string("A")}) {
            for (bool prefix : {false,true}) {
                auto f = base({"Z"}); const bool alt = count == "A", fixed = count == "3";
                if(alt) { alleles(f); if(!prefix) { f.d3.var_id.clear(); f.d3.samp_id.clear(); } }
                std::string value;
                if(type == "Integer") {
                    value = alt ? "7,9" : fixed ? "7,9,11" : "7";
                    if(alt) ints(f.d4.samp_int,"X",{7,9});
                    else if(fixed) { ints(f.d3.samp_int,"X0",{7}); ints(f.d3.samp_int,"X1",{9}); ints(f.d3.samp_int,"X2",{11}); }
                    else ints(f.d3.samp_int,"X",{7});
                } else if(type == "Float") {
                    value = alt ? "0.25,0.5" : fixed ? "0.25,0.5,0.75" : "0.25";
                    if(alt) floats(f.d4.samp_float,"X",{0.25f,0.5f});
                    else if(fixed) { floats(f.d3.samp_float,"X0",{0.25f}); floats(f.d3.samp_float,"X1",{0.5f}); floats(f.d3.samp_float,"X2",{0.75f}); }
                    else floats(f.d3.samp_float,"X",{0.25f});
                } else {
                    value = alt ? "a,b" : fixed ? "a,b,c" : "a";
                    if(alt) strings(f.d4.samp_string,"X",{"a","b"});
                    else if(fixed) { strings(f.d3.samp_string,"X0",{"a"}); strings(f.d3.samp_string,"X1",{"b"}); strings(f.d3.samp_string,"X2",{"c"}); }
                    else strings(f.d3.samp_string,"X",{"a"});
                }
                for(auto& c:f.d3.samp_int) c.numb=fixed?3:1;
                for(auto& c:f.d3.samp_float) c.numb=fixed?3:1;
                for(auto& c:f.d3.samp_string) c.numb=fixed?3:1;
                for(auto& c:f.d4.samp_int) c.numb=-1;
                for(auto& c:f.d4.samp_float) c.numb=-1;
                for(auto& c:f.d4.samp_string) c.numb=-1;
                if(fixed) { f.d2.var_id={0}; f.d2.alt_id={0}; f.d2.alt={"C"}; f.d2.numAlt=1; }
                auto meta=declaration("FORMAT","X",count,type);
                if(prefix) { f.d3.sample_GT.push_back({{static_cast<char>(122)},1}); meta=declaration("FORMAT","GT","1","String")+meta; }
                auto input=header(meta,"Z")+row(".",(prefix?"GT:X\t0/1:":"X\t")+value);
                if(fixed) input.replace(input.find("\tC,G\t"),5,"\tC\t");
                out.push_back({"format/"+type+"/"+count+(prefix?"/with_GT":"/alone"),input,f});
            }
        }
    }
    {
        auto f=base(); info_flag flag{}; flag.name="DB"; flag.i_flag={1}; f.d1.in_flag.push_back(flag);
        out.push_back({"info/flag_present",header(declaration("INFO","DB","0","Flag"))+row("DB"),f});
        f.d1.in_flag[0].i_flag={0};
        out.push_back({"info/flag_absent",header(declaration("INFO","DB","0","Flag"))+row(),f});
    }
    {
        auto f=base({"Z","A"});
        f.d3.sample_GT.push_back({{static_cast<char>(1),static_cast<char>(133)},1});
        ints(f.d3.samp_int,"DP",{7,0}); f.d3.samp_int[0].numb=1;
        alleles(f); ints(f.d4.samp_int,"XA",{1,2,3,4}); f.d4.samp_int[0].numb=-1;
        out.push_back({"samples/order_and_alleles",header(declaration("FORMAT","GT","1","String")+declaration("FORMAT","DP","1","Integer")+declaration("FORMAT","XA","A","Integer"),"Z\tA")+row(".","GT:DP:XA\t0|1:7:1,2\t1/1:0:3,4"),f});
    }
    for(const auto& gt: {std::string("./."),std::string(".|.")}) {
        auto f=base({"Z"}); f.d3.sample_GT.push_back({{static_cast<char>(gt=="./."?255:254)},1});
        ints(f.d3.samp_int,"DP",{12}); f.d3.samp_int[0].numb=1;
        out.push_back({"missing/GT_"+(gt=="./."?std::string("unphased"):std::string("phased")),header(declaration("FORMAT","GT","1","String")+declaration("FORMAT","DP","1","Integer"),"Z")+row(".","GT:DP\t"+gt+":12"),f});
    }
    for(const auto& type:{std::string("Integer"),std::string("Float")}) {
        auto f=base({"Z"});
        if(type=="Integer") { ints(f.d3.samp_int,"X0",{missing_int()}); ints(f.d3.samp_int,"X1",{7}); for(auto& c:f.d3.samp_int)c.numb=2; }
        else { floats(f.d3.samp_float,"X0",{missing_float()}); floats(f.d3.samp_float,"X1",{0.5f}); for(auto& c:f.d3.samp_float)c.numb=2; }
        out.push_back({"missing/array_"+type,header(declaration("FORMAT","X","2",type),"Z")+row(".","X\t.,"+(type=="Integer"?std::string("7"):std::string("0.5"))),f});
    }
    {
        auto f=base({"Z"}); ints(f.d3.samp_int,"DP",{12}); f.d3.samp_int[0].numb=1;
        strings(f.d3.samp_string,"XS",{"."}); f.d3.samp_string[0].numb=1;
        out.push_back({"missing/string",header(declaration("FORMAT","DP","1","Integer")+declaration("FORMAT","XS","1","String"),"Z")+row(".","DP:XS\t12:."),f});
    }
    for(const auto& id:{std::string("DP2"),std::string("X1")}) {
        auto f=base({"Z"}); ints(f.d3.samp_int,id,{7}); f.d3.samp_int[0].numb=1;
        std::string meta=declaration("FORMAT",id,"1","Integer"), sample=id+"\t7";
        if(id=="X1") { ints(f.d3.samp_int,"X2",{9}); f.d3.samp_int[1].numb=1; meta+=declaration("FORMAT","X2","1","Integer"); sample="X1:X2\t7:9"; }
        out.push_back({"identifiers/"+id,header(meta,"Z")+row(".",sample),f});
    }
    {
        auto f=base(); floats(f.d2.alt_float,"AF",{0.25f,0.5f}); floats(f.d2.alt_float,"AF2",{0.5f,0.75f});
        out.push_back({"identifiers/AF_and_AF2",header(declaration("INFO","AF","A","Float")+declaration("INFO","AF2","A","Float"))+row("AF=0.25,0.5;AF2=0.5,0.75"),f});
        auto g=base({"Z"}); ints(g.d3.samp_int,"AF20",{7}); ints(g.d3.samp_int,"AF21",{9}); for(auto& c:g.d3.samp_int)c.numb=2;
        out.push_back({"identifiers/array_AF2",header(declaration("FORMAT","AF2","2","Integer"),"Z")+row(".","AF2\t7,9"),g});
    }
    for(const auto& tsa: std::vector<std::pair<std::string,int>>{{"SNV",0},{"insertion",1},{"deletion",2},{"inversion",3}}) {
        auto f=base(); ints(f.d1.in_int,"TSA",{tsa.second});
        out.push_back({"special/TSA_"+tsa.first,header(declaration("INFO","TSA","1","String"))+row("TSA="+tsa.first),f});
    }
    for(const auto& n:{std::string("R"),std::string("G"),std::string(".")}) {
        auto f=base({"Z"}); ints(f.d1.in_int,"DP",{12}); ints(f.d3.samp_int,"Y",{9}); f.d3.samp_int[0].numb=1;
        std::string vals=n=="R"?"1,2,3":n=="G"?"1,2,3,4,5,6":"1,2";
        out.push_back({"unsupported/isolation_"+(n=="."?std::string("variable"):n),header(declaration("INFO","U",n,"Integer")+declaration("INFO","DP","1","Integer")+declaration("FORMAT","X",n,"Integer")+declaration("FORMAT","Y","1","Integer"),"Z")+row("U="+vals+";DP=12","X:Y\t"+vals+":9"),f,false});
    }
    {
        auto f=base(); ints(f.d1.in_int,"DP",{12});
        out.push_back({"unsupported/info_fixed",header(declaration("INFO","U","2","Integer")+declaration("INFO","DP","1","Integer"))+row("U=1,2;DP=12"),f,false});
    }
    {
        auto f=base();
        f.d1.var_number={0,1}; f.d1.chrom={0,1}; f.d1.chrom_map={{"1",0},{"2",1}};
        f.d1.pos={10,20}; f.d1.id={"rs1","."}; f.d1.ref={"A","T"}; f.d1.qual={half_value(30),half_value(-1)};
        f.d1.filter={0,1}; f.d1.filter_map={{"PASS",0},{"q10",1}};
        f.d2.var_id={0,0,1}; f.d2.alt_id={0,1,0}; f.d2.alt={"C","G","A"}; f.d2.numAlt=3;
        out.push_back({"core/multiple_records_missing_QUAL",header()+row()+"2\t20\t.\tT\tA\t.\tq10\t.\n",f});
    }

    for(const auto& type:{std::string("Integer"),std::string("Float")}) {
        for(bool alt:{false,true}) {
            auto f=base({"Z"}); if(alt) { alleles(f); f.d3.var_id.clear(); f.d3.samp_id.clear(); }
            const auto n=alt?"A":"1";
            if(type=="Integer") {
                if(alt) { ints(f.d2.alt_int,"I",{missing_int(),7}); ints(f.d4.samp_int,"F",{missing_int(),9}); f.d4.samp_int[0].numb=-1; }
                else { ints(f.d1.in_int,"I",{missing_int()}); ints(f.d3.samp_int,"F",{missing_int()}); f.d3.samp_int[0].numb=1; }
            } else {
                if(alt) { floats(f.d2.alt_float,"I",{missing_float(),0.5f}); floats(f.d4.samp_float,"F",{missing_float(),0.25f}); f.d4.samp_float[0].numb=-1; }
                else { floats(f.d1.in_float,"I",{missing_float()}); floats(f.d3.samp_float,"F",{missing_float()}); f.d3.samp_float[0].numb=1; }
            }
            // Missing scalar INFO may be omitted; use omission as the expected textual representation.
            auto info=alt?(type=="Integer"?"I=.,7":"I=.,0.5"):".";
            auto sample=alt?(type=="Integer"?"F\t.,9":"F\t.,0.25"):"F\t.";
            out.push_back({"missing/"+type+(alt?"_A":"_scalar"),header(declaration("INFO","I",n,type)+declaration("FORMAT","F",n,type),"Z")+row(info,sample),f});
        }
    }
    {
        auto f=base({"Z"}); floats(f.d3.samp_float,"F",{0.1f}); f.d3.samp_float[0].numb=1;
        out.push_back({"float/half_precision",header(declaration("FORMAT","F","1","Float"),"Z")+row(".","F\t0.1"),f});
    }
    {
        auto f=base({"Z"}); ints(f.d3.samp_int,"GQX",{4}); ints(f.d3.samp_int,"GQ",{9});
        for(auto& c:f.d3.samp_int)c.numb=1;
        out.push_back({"identifiers/prefix_overlap",header(declaration("FORMAT","GQX","1","Integer")+declaration("FORMAT","GQ","1","Integer"),"Z")+row(".","GQ:GQX\t9:4"),f});
    }
    {
        auto f=base({"Z","A"}); f.d1.var_number={0,1}; f.d1.chrom={0,0}; f.d1.pos={10,20};
        f.d1.id={"rs1","rs2"}; f.d1.ref={"A","T"}; f.d1.qual={half_value(30),half_value(40)}; f.d1.filter={0,0};
        f.d2.var_id={0,0,1}; f.d2.alt_id={0,1,0}; f.d2.alt={"C","G","A"}; f.d2.numAlt=3;
        f.d3.var_id={0,0,1,1}; f.d3.samp_id={0,1,0,1};
        f.d3.sample_GT.push_back({{static_cast<char>(122),static_cast<char>(1),static_cast<char>(255),static_cast<char>(255)},1});
        ints(f.d3.samp_int,"DP",{7,9,11,13}); f.d3.samp_int[0].numb=1;
        out.push_back({"samples/changing_FORMAT",header(declaration("FORMAT","GT","1","String")+declaration("FORMAT","DP","1","Integer"),"Z\tA")+row(".","GT:DP\t0/1:7\t0|1:9")+"1\t20\trs2\tT\tA\t40\tPASS\t.\tDP\t11\t13\n",f});
    }
    {
        auto f=base(); f.d2.alt={"<DEL>","*"};
        auto input=header()+row(); input.replace(input.find("\tC,G\t"),5,"\t<DEL>,*\t");
        out.push_back({"core/symbolic_and_spanning_ALT",input,f});
    }

    {
        auto f=base({"Z"}); ints(f.d1.in_int,"X",{7}); strings(f.d3.samp_string,"X",{"hello"}); f.d3.samp_string[0].numb=1;
        out.push_back({"identifiers/INFO_FORMAT_same_ID",header(declaration("INFO","X","1","Integer")+declaration("FORMAT","X","1","String"),"Z")+row("X=7","X\thello"),f});
    }
    {
        auto f=base(); ints(f.d1.in_int,"DP",{7});
        out.push_back({"header/reordered_attributes",header("##INFO=<Number=1,Type=Integer,ID=DP,Description=\"depth, test\">\n")+row("DP=7"),f});
    }
    return out;
}
}
