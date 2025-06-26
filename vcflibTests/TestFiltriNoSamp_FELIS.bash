#!/bin/bash

# Funzione per eseguire il filtro, cronometrare l'esecuzione e salvare l'output.
# Il terzo parametro, flag, indica se si utilizza -f (INFO filter) oppure -g (Genotype filter).

filter_expr="EVA_4 - felis_catus"
echo "Esecuzione filtro: $filter_expr" > ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "EVA_4" ../data/felis_catus.vcf   > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt


filter_expr="E_Multiple_observations - felis_catus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "E_Multiple_observations" ../data/felis_catus.vcf   > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt


filter_expr="TSA=SNV - felis_catus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "TSA = SNV" ../data/felis_catus.vcf   > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt

#CREAZIONE FILE
#echo "INIZIO CREAZIONE FILE"
#grep -v "^#" ../data/felis_catus.vcf | awk '{print $1"\t"($2-1)"\t"$2"\t"$2}' > myannotations.bed
#vcfannotate -b myannotations.bed -k mypos ../data/felis_catus.vcf > data/felis_catus_annotated.vcf
#sed '/^##INFO=<ID=mypos,/ s/Type=String/Type=Integer/' data/felis_catus_annotated.vcf > data/felis_catus_annotated_new.vcf
#sed -E 's/(mypos=[0-9]+)(:[0-9]+)+/\1/g' data/felis_catus_annotated_new.vcf > data/felis_catus_annotated_new_fixed.vcf 
#echo "FILE CREATO"

filter_expr="POS > 140680000 - felis_catus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "mypos > 140680000" data/felis_catus_annotated_new_fixed.vcf   > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt

filter_expr="POS > 140680000 & POS < 160680000 - felis_catus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "mypos > 140680000 & mypos < 160680000" data/felis_catus_annotated_new_fixed.vcf   > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt

filter_expr="! E_Multiple_observations and TSA=SNV - felis_catus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "! E_Multiple_observations & TSA = SNV" ../data/felis_catus.vcf > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt

filter_expr="POS > 140680000 & E_Multiple_observations - felis_catus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "mypos > 140680000 & E_Multiple_observations" data/felis_catus_annotated_new_fixed.vcf   > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt

filter_expr="POS > 140680000 & POS < 160680000 & E_Multiple_observations - felis_catus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "mypos > 140680000 & mypos < 160680000 & E_Multiple_observations" data/felis_catus_annotated_new_fixed.vcf   > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt

filter_expr="POS > 140680000 & TSA=SNV - felis_catus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "mypos > 140680000 & TSA = SNV" data/felis_catus_annotated_new_fixed.vcf   > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt

filter_expr="POS > 140680000 & POS < 160680000 & TSA=SNV - felis_catus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "mypos > 140680000 & mypos < 160680000 & TSA = SNV" data/felis_catus_annotated_new_fixed.vcf   > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt


filter_expr="EVA_4 & POS > 140680000 & TSA=SNV - felis_catus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "EVA_4 & mypos > 140680000 & TSA = SNV" data/felis_catus_annotated_new_fixed.vcf   > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt

filter_expr="EVA_4 & POS > 140680000 & POS < 160680000 & TSA=SNV - felis_catus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_results_felis.txt
start=$(date +%s%N)
vcffilter -f "EVA_4 & mypos > 140680000 & mypos < 160680000 & TSA = SNV" data/felis_catus_annotated_new_fixed.vcf   > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_results_felis.txt
echo "" >> ../result/Vcflib_results_felis.txt