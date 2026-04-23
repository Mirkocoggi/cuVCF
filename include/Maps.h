#ifndef MAPS_H
#define MAPS_H

#include <map>
#include <string>

const std::map<std::string, char> polyphenCharMap = {
    {"", 0},
    {"benign", 1},
    {"possibly_damaging", 2},
    {"probably_damaging", 3}
};

const std::map<std::string, char> csqCharMap = {
    {"synonymous_variant", 0},
    {"missense_variant", 1},
    {"nonsense_variant", 2},
    {"frameshift_variant", 3},
    {"inframe_insertion", 4},
    {"inframe_deletion", 5},
    {"5_prime_UTR_variant", 6},
    {"3_prime_UTR_variant", 7},
    {"splice_region_variant", 8},
    {"splice_acceptor_variant", 9},
    {"splice_donor_variant", 10},
    {"regulatory_region_variant", 11},
    {"upstream_gene_variant", 12},
    {"downstream_gene_variant", 13},
    {"coding_sequence_variant", 14},
    {"non_coding_transcript_variant", 15},
    {"mature_miRNA_variant", 16},
    {"ncRNA_exon_variant", 17},
    {"ncRNA_intronic_variant", 18},
    {"protein_altering_variant", 19},
    {"transcript_ablation", 20},
    {"transcript_amplification", 21},
    {"intron_variant", 22},
    {"intergenic_variant", 23},
    {"enhancer_variant", 24},
    {"regulatory_region_variant", 25},
    {"non_coding_exon_variant", 26},
    {"pseudogene_variant", 27}
};

const std::map<std::string, int> tsaCharMap = {
    {"SNV",       0},
    {"INS",       1},
    {"insertion", 1},
    {"DEL",       2},
    {"deletion",  2},
    {"INV",       3},
    {"inversion", 3}
};

#endif