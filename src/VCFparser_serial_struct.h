#ifndef VCF_STRUCTS_H
#define VCF_STRUCTS_H
#include <chrono>


using namespace std;

class var
{
public:
	string chrom;
    int pos;
    int id;
    string ref;
    string alt;
    int qual;
    string filter;
    string info;
    string format;
    vector<string> samples;
};

#endif