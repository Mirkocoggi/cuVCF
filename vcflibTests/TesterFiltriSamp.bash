#!/bin/bash

#CREAZIONE FILE
grep -v "^#" ../data/IRBT.vcf | awk '{print $1"\t"($2-1)"\t"$2"\t"$2}' > myannotations.bed
vcfannotate -b myannotations.bed -k mypos ../data/IRBT.vcf > ../data/IRBT_annotated.vcf
sed '/^##INFO=<ID=mypos,/ s/Type=String/Type=Integer/' ../data/IRBT_annotated.vcf > ../data/IRBT_annotated_new.vcf
sed -E 's/(mypos=[0-9]+)(:[0-9]+)+/\1/g' ../data/IRBT_annotated_new.vcf > ../data/IRBT_annotated_new_fixed.vcf

filter_expr="AC > 8"
echo "Esecuzione filtro: $filter_expr" > ../result/Vcflib_result_IRBT.txt
echo "Esecuzione filtro: $filter_expr" > ../result/Vcflib_result_IRBT.txt
start=$(date +%s%N)
vcffilter -f "AC > 8" ../data/IRBT.vcf  > output.vcf
vcffilter -f "AC > 8" ../data/IRBT.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_sec=$(echo "scale=3; $runtime_ns/1000000000" | bc -l)
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt


filter_expr="AC > 8 & AF > 0.5"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
start=$(date +%s%N)
vcffilter -f "AC > 8 & AF > 0.5" ../data/IRBT.vcf  > output.vcf
vcffilter -f "AC > 8 & AF > 0.5" ../data/IRBT.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_sec=$(echo "scale=3; $runtime_ns/1000000000" | bc -l)
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt


filter_expr="GT = 1|1"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
start=$(date +%s%N)
vcffilter -g "GT = 1|1" ../data/IRBT.vcf  > output.vcf
vcffilter -g "GT = 1|1" ../data/IRBT.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_sec=$(echo "scale=3; $runtime_ns/1000000000" | bc -l)
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt

#NOTA - AD[0] > 3 - fatto con escamotage
filter_expr="AD > 3"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
start=$(date +%s%N)
vcffilter -g "AD > 3" ../data/IRBT.vcf  > output.vcf
vcffilter -g "AD > 3" ../data/IRBT.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_sec=$(echo "scale=3; $runtime_ns/1000000000" | bc -l)
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt

#NOTA - AD[0] > 3 & AD[1] < 10 - test failed

filter_expr="AC > 8 & GT = 1|1"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
start=$(date +%s%N)
vcffilter -f "AC > 8" ../data/IRBT.vcf > output.vcf
vcffilter -f "AC > 8" ../data/IRBT.vcf > output.vcf
vcffilter -g "GT = 1|1" output.vcf  > output1.vcf
end=$(date +%s%N)
rm output1.vcf
runtime_ns=$((end - start))
runtime_sec=$(echo "scale=3; $runtime_ns/1000000000" | bc -l)
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt

filter_expr="AC > 8 & AF > 0.5 & GT = 1|1"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
start=$(date +%s%N)
vcffilter -f "AC > 8 & AF > 0.5" ../data/IRBT.vcf  > output.vcf
vcffilter -f "AC > 8 & AF > 0.5" ../data/IRBT.vcf  > output.vcf
vcffilter -g "GT = 1|1" output.vcf  > output1.vcf
end=$(date +%s%N)
rm output1.vcf
runtime_ns=$((end - start))
runtime_sec=$(echo "scale=3; $runtime_ns/1000000000" | bc -l)
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt

#NOTA - AD[0] > 3 - fatto con escamotage
filter_expr="AC > 8 & AD > 3"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
start=$(date +%s%N)
vcffilter -f "AC > 8" ../data/IRBT.vcf > output.vcf
vcffilter -f "AC > 8" ../data/IRBT.vcf > output.vcf
vcffilter -g "AD > 3" output.vcf  > output1.vcf
end=$(date +%s%N)
rm output1.vcf
runtime_ns=$((end - start))
runtime_sec=$(echo "scale=3; $runtime_ns/1000000000" | bc -l)
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt

#NOTA - AC > 8 & AD[0] > 3 & AD[1] < 10 - test failed

#NOTA - AD[0] > 3 - fatto con escamotage
filter_expr="AC > 8 & AF > 0.5 & AD > 3"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
start=$(date +%s%N)
vcffilter -f "AC > 8 & AF > 0.5" ../data/IRBT.vcf > output.vcf
vcffilter -f "AC > 8 & AF > 0.5" ../data/IRBT.vcf > output.vcf
vcffilter -g "AD > 3" output.vcf  > output1.vcf
end=$(date +%s%N)
rm output1.vcf
runtime_ns=$((end - start))
runtime_sec=$(echo "scale=3; $runtime_ns/1000000000" | bc -l)
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt

#NOTA - AC > 8 & AF > 0.5 & AD[0] > 3 & AD[1] < 10 - test failed

#NOTA - AD[0] > 3 - fatto con escamotage - output non corretto ma runna FAILED?
filter_expr="AD > 3 & GT = 1|1"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
start=$(date +%s%N)
vcffilter -g "GT = 1|1 & AD > 3" ../data/IRBT.vcf  > output.vcf
vcffilter -g "GT = 1|1 & AD > 3" ../data/IRBT.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_sec=$(echo "scale=3; $runtime_ns/1000000000" | bc -l)
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt

#NOTA - AD[0] > 3 & AD[1] < 10 & GT = 1|1 - test failed

#NOTA - AD[0] > 3 - fatto con escamotage - output non corretto ma runna FAILED?
filter_expr="AC > 8 & GT = 1|1 & AD > 3"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
start=$(date +%s%N)
vcffilter -f "AC > 8" ../data/IRBT.vcf > output.vcf
vcffilter -f "AC > 8" ../data/IRBT.vcf > output.vcf
vcffilter -g "GT = 1|1 & AD > 3" output.vcf  > output1.vcf
end=$(date +%s%N)
rm output1.vcf
runtime_ns=$((end - start))
runtime_sec=$(echo "scale=3; $runtime_ns/1000000000" | bc -l)
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt

#NOTA - AC > 8 & AD[0] > 3 & AD[1] < 10 & GT = 1|1 - test failed

#NOTA - AD[0] > 3 - fatto con escamotage - output non corretto ma runna FAILED?
filter_expr="AC > 8 & AF > 0.5 & GT = 1|1 & AD > 3"
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
echo "Esecuzione filtro: $filter_expr" >> ../result/Vcflib_result_IRBT.txt
start=$(date +%s%N)
vcffilter -f "AC > 8 & AF > 0.5" ../data/IRBT.vcf > output.vcf
vcffilter -f "AC > 8 & AF > 0.5" ../data/IRBT.vcf > output.vcf
vcffilter -g "GT = 1|1 & AD > 3" output.vcf  > output1.vcf
end=$(date +%s%N)
rm output1.vcf
runtime_ns=$((end - start))
runtime_sec=$(echo "scale=3; $runtime_ns/1000000000" | bc -l)
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt
echo "Tempo: ${runtime_sec} s" >> ../result/Vcflib_result_IRBT.txt
echo "" >> ../result/Vcflib_result_IRBT.txt

#NOTA - AC > 8 & AF > 0.5 & GT = 1|1 & AD[0] > 3 & AD[1] < 10 - test failed