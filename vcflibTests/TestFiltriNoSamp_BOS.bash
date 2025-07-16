#!/bin/bash

filter_expr="EVA_4 - bos_taurus"
echo "Esecuzione filtro: $filter_expr" > ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" > ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "EVA_4" ../data/bos_taurus.vcf  > output.vcf
vcffilter -f "EVA_4" ../data/bos_taurus.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt


filter_expr="E_Multiple_observations - bos_taurus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "E_Multiple_observations" ../data/bos_taurus.vcf  > output.vcf
vcffilter -f "E_Multiple_observations" ../data/bos_taurus.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt


filter_expr="TSA=SNV - bos_taurus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "TSA = SNV" ../data/bos_taurus.vcf  > output.vcf
vcffilter -f "TSA = SNV" ../data/bos_taurus.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt

#CREAZIONE FILE
grep -v "^#" ../data/bos_taurus.vcf | awk '{print $1"\t"($2-1)"\t"$2"\t"$2}' > myannotations.bed
vcfannotate -b myannotations.bed -k mypos ../data/bos_taurus.vcf > ../data/bos_taurus_annotated.vcf
sed '/^##INFO=<ID=mypos,/ s/Type=String/Type=Integer/' ../data/bos_taurus_annotated.vcf > ../data/bos_taurus_annotated_new.vcf
sed -E 's/(mypos=[0-9]+)(:[0-9]+)+/\1/g' ../data/bos_taurus_annotated_new.vcf > ../data/bos_taurus_annotated_new_fixed.vcf

filter_expr="POS > 200000 - bos_taurus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "mypos > 200000" ../data/bos_taurus_annotated_new_fixed.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt

filter_expr="POS > 200000 & POS < 300000 - bos_taurus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "mypos > 200000 & mypos < 300000" ../data/bos_taurus_annotated_new_fixed.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt

filter_expr="! E_Multiple_observations and TSA=SNV - bos_taurus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "! E_Multiple_observations & TSA = SNV" ../data/bos_taurus.vcf  > output.vcf
vcffilter -f "! E_Multiple_observations & TSA = SNV" ../data/bos_taurus.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt

filter_expr="POS > 200000 & E_Multiple_observations - bos_taurus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "mypos > 200000 & E_Multiple_observations" ../data/bos_taurus_annotated_new_fixed.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt

filter_expr="POS > 200000 & POS < 300000 & E_Multiple_observations - bos_taurus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "mypos > 200000 & mypos < 300000 & E_Multiple_observations" ../data/bos_taurus_annotated_new_fixed.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt

filter_expr="POS > 200000 & TSA=SNV - bos_taurus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "mypos > 200000 & TSA = SNV" ../data/bos_taurus_annotated_new_fixed.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt

filter_expr="POS > 200000 & POS < 300000 & TSA=SNV - bos_taurus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "mypos > 200000 & mypos < 300000 & TSA = SNV" ../data/bos_taurus_annotated_new_fixed.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt

filter_expr="EVA_4 & POS > 200000 & TSA=SNV - bos_taurus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "EVA_4 & mypos > 200000 & TSA = SNV" ../data/bos_taurus_annotated_new_fixed.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt

filter_expr="EVA_4 & POS > 200000 & POS < 300000 & TSA=SNV - bos_taurus"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_bos.txt
start=$(date +%s%N)
vcffilter -f "EVA_4 & mypos > 200000 & mypos < 300000 & TSA = SNV" ../data/bos_taurus_annotated_new_fixed.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt
echo "Tempo: ${runtime_ms} ms" >> ../result/Vcflib_result_bos.txt
echo "" >> ../result/Vcflib_result_bos.txt