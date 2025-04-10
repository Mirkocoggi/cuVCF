#!/bin/bash

filter_expr="AF > 0.5"
echo "Esecuzione filtro: $filter_expr"
start=$(date +%s%N)
vcffilter -f "AF > 0.5" data/IRBT.vcf > "output.vcf"
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo per filtro ($output_file): ${runtime_ms} ms"
echo ""


filter_expr="AC >= 10"
echo "Esecuzione filtro: $filter_expr"
start=$(date +%s%N)
vcffilter -f "AF > 0.5" data/IRBT.vcf > "output.vcf"
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo per filtro ($output_file): ${runtime_ms} ms"
echo ""

# Filtro 1: Frequenza Allelica (AF[0]>=0.5) – filtro INFO
run_filter 'AF[0]>=0.5' output_AF.vcf -f

# Filtro 2: Allele Count (AC>=10) – filtro INFO
run_filter 'AC>=10' output_AC.vcf -f

# Filtro 3: Completezza del Genotipo:
#    (count(gt())-count(gt("./.")))/count(gt())>=0.75
# Usa il filtro Genotype (-g)
run_filter '(count(gt())-count(gt("./.")))/count(gt())>=0.75' output_complete.vcf -g

# Filtro 4: Profondità dell’Allele Alternativo (AD):
#    (sum(AD[1]>=8))/count(AD[1])>=0.5
# Usa il filtro Genotype (-g)
run_filter '(sum(AD[1]>=8))/count(AD[1])>=0.5' output_AD.vcf -g

# Filtro 5: Intervallo Genomico (Window Position):
#    CHROM=="1" && POS>=150 && POS<=3000
# Usa il filtro INFO (-f)
run_filter 'CHROM=="1" && POS>=150 && POS<=3000' output_window.vcf -f
