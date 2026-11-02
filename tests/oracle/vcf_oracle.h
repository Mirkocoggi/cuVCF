#pragma once
#include "support/harness.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

// A deliberately small checker for our supported synthetic VCF subset. It does not
// use cuVCF to parse reconstructed output, so parser/writer bugs cannot cancel out.
namespace test {
inline std::vector<std::string> tokens(const std::string& s, char delim) {
    std::vector<std::string> out; size_t start=0;
    for(size_t end; (end=s.find(delim,start))!=std::string::npos; start=end+1) out.push_back(s.substr(start,end-start));
    out.push_back(s.substr(start)); return out;
}
struct Field { std::string number, type; };
struct Vcf {
    std::set<std::string> metadata;
    std::vector<std::string> columns;
    std::map<std::string,Field> info, format;
    std::vector<std::vector<std::string>> rows;
};
inline std::string attribute(const std::string& text, const std::string& key) {
    auto start=text.find(key+"="); require(start!=std::string::npos,"header/attribute","missing "+key);
    start+=key.size()+1; auto end=text.find_first_of(",>",start); return text.substr(start,end-start);
}
inline Vcf read_vcf(const std::string& text) {
    Vcf v; std::istringstream input(text); std::string line; bool columns=false;
    while(std::getline(input,line)) {
        require(!line.empty(),"syntax/blank_line","empty VCF line");
        if(line.rfind("##",0)==0) {
            require(!columns,"header/order","metadata after column header");
            require(v.metadata.insert(line).second,"header/duplicate","duplicate metadata line");
            for(const auto& scope:{std::string("INFO"),std::string("FORMAT")}) if(line.rfind("##"+scope+"=<",0)==0) {
                auto& fields=scope=="INFO"?v.info:v.format;
                require(fields.emplace(attribute(line,"ID"),Field{attribute(line,"Number"),attribute(line,"Type")}).second,"header/duplicate_id","duplicate "+scope+" ID");
            }
        } else if(line.rfind("#CHROM\t",0)==0) {
            require(!columns,"header/columns","repeated column header"); columns=true; v.columns=tokens(line,'\t');
            require(v.columns.size()==8 || v.columns.size()>=10,"header/width","invalid header width");
        } else {
            require(columns,"header/missing","record before column header");
            auto row=tokens(line,'\t'); expect_equal(row.size(),v.columns.size(),"record/width");
            for(const auto& cell:row) require(!cell.empty(),"syntax/empty_cell","empty VCF cell");
            v.rows.push_back(std::move(row));
        }
    }
    require(columns,"header/missing","missing column header"); return v;
}
inline std::map<std::string,std::string> info_values(const std::string& cell, const Vcf& v) {
    std::map<std::string,std::string> values; if(cell==".")return values;
    for(const auto& token:tokens(cell,';')) {
        auto pos=token.find('='); auto name=token.substr(0,pos);
        require(v.info.count(name),"info/undeclared","undeclared INFO ID "+name);
        bool flag=v.info.at(name).type=="Flag";
        require(flag ? pos==std::string::npos : pos!=std::string::npos,"info/assignment","invalid assignment for "+name);
        require(values.emplace(name,flag?"<flag>":token.substr(pos+1)).second,"info/duplicate","duplicate INFO ID "+name);
    }
    return values;
}
inline std::map<std::string,std::string> sample_values(const std::string& format, const std::string& cell, const Vcf& v) {
    auto keys=tokens(format,':'), vals=tokens(cell,':');
    require(vals.size()<=keys.size(),"format/width","more sample values than FORMAT IDs");
    auto gt=std::find(keys.begin(),keys.end(),"GT");
    require(gt==keys.end() || gt==keys.begin(),"format/GT_order","GT must be first when present");
    std::map<std::string,std::string> result;
    for(size_t i=0;i<keys.size();++i) {
        require(v.format.count(keys[i]),"format/undeclared","undeclared FORMAT ID "+keys[i]);
        auto val=i<vals.size()?vals[i]:".";
        require(!val.empty(),"format/empty","empty FORMAT value for "+keys[i]);
        require(result.emplace(keys[i],val).second,"format/duplicate","duplicate FORMAT ID "+keys[i]);
    }
    return result;
}
inline double numeric(const std::string& s, const std::string& key) {
    size_t consumed=0; double n=0;
    try { n=std::stod(s,&consumed); } catch(const std::exception&) { throw Assertion(key,"not numeric: "+s); }
    require(consumed==s.size() && std::isfinite(n),key,"invalid finite number: "+s); return n;
}
inline double quantize_half(double x) {
    if(x==0) return x;
    int exponent=0; std::frexp(std::abs(x),&exponent);
    double step=std::ldexp(1.0,std::max(-24,exponent-11));
    return std::nearbyint(x/step)*step;
}
inline void value_equal(const std::string& actual, const std::string& expected,
                        const Field& field, size_t alt_count, const std::string& key) {
    if(field.type=="Flag") { expect_equal(actual,expected,key); return; }
    auto a=tokens(actual,','), e=tokens(expected,',');
    const size_t count=field.number=="A"?alt_count:static_cast<size_t>(std::stoul(field.number));
    // A single dot represents an entirely missing field, including an array.
    if(actual==".") a.assign(count,".");
    if(expected==".") e.assign(count,".");
    expect_equal(a.size(),count,key+"/cardinality"); expect_equal(e.size(),count,key+"/expected_cardinality");
    for(size_t i=0;i<count;++i) {
        const auto at=key+"/"+std::to_string(i);
        if(a[i]=="." || e[i]=="." || field.type=="String") expect_equal(a[i],e[i],at);
        else if(field.type=="Integer") {
            require(a[i].find_first_not_of("-+0123456789")==std::string::npos,at,"invalid integer token");
            expect_equal(numeric(a[i],at),numeric(e[i],at),at);
        } else {
            const auto av=numeric(a[i],at), ev=quantize_half(numeric(e[i],at));
            // Writers print six decimal places; allow text rounding after half conversion.
            require(quantize_half(av)==ev || std::abs(av-ev)<=0.00000051,at,"expected half-quantized "+std::to_string(ev)+", got "+a[i]);
        }
    }
}
inline void compare_vcf(const std::string& actual_text,const std::string& expected_text) {
    auto a=read_vcf(actual_text), e=read_vcf(expected_text);
    require(a.metadata==e.metadata,"header/metadata","metadata changed");
    require(a.columns==e.columns,"header/samples","column/sample order changed");
    expect_equal(a.rows.size(),e.rows.size(),"records/count");
    for(size_t i=0;i<e.rows.size();++i) {
        const auto& ar=a.rows[i]; const auto& er=e.rows[i]; const auto key="record/"+std::to_string(i);
        for(size_t col:{0u,1u,2u,3u,4u,6u}) expect_equal(ar[col],er[col],key+"/column/"+std::to_string(col));
        value_equal(ar[5],er[5],{"1","Float"},1,key+"/QUAL");
        auto ai=info_values(ar[7],a), ei=info_values(er[7],e);
        auto alts=tokens(er[4],',').size();
        // Validate first, then canonicalize wholly missing non-flag values. This
        // preserves duplicate/undeclared-key errors and partial-array positions.
        auto normalize=[&](auto& values,const auto& fields) {
            for(auto it=values.begin();it!=values.end();) {
                const auto& field=fields.at(it->first);
                value_equal(it->second,it->second,field,alts,key+"/INFO/"+it->first);
                auto parts=tokens(it->second,',');
                if(field.type!="Flag" && std::all_of(parts.begin(),parts.end(),[](const auto& part){return part==".";}))
                    it=values.erase(it);
                else ++it;
            }
        };
        normalize(ai,a.info); normalize(ei,e.info);
        expect_equal(ai.size(),ei.size(),key+"/info/count");
        for(const auto& kv:ei) {
            require(ai.count(kv.first),key+"/info/key","missing INFO "+kv.first);
            value_equal(ai.at(kv.first),kv.second,e.info.at(kv.first),alts,key+"/INFO/"+kv.first);
        }
        if(er.size()>8) for(size_t s=9;s<er.size();++s) {
            auto as=sample_values(ar[8],ar[s],a), es=sample_values(er[8],er[s],e);
            // The writer may emit a dataset-wide FORMAT union; added values must be missing.
            std::set<std::string> names; for(const auto& kv:as)names.insert(kv.first); for(const auto& kv:es)names.insert(kv.first);
            for(const auto& name:names) {
                auto av=as.count(name)?as.at(name):".", ev=es.count(name)?es.at(name):".";
                if(name=="GT") {
                    if(ev==".") require(av=="." || av=="./." || av==".|.",key+"/GT","invented genotype");
                    else expect_equal(av,ev,key+"/sample/"+std::to_string(s-9)+"/GT");
                } else value_equal(av,ev,e.format.at(name),alts,key+"/sample/"+std::to_string(s-9)+"/"+name);
            }
        }
    }
}
}
