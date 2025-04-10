#!/bin/bash

# Funzione per eseguire il filtro, cronometrare l'esecuzione e salvare l'output.
# Il terzo parametro, flag, indica se si utilizza -f (INFO filter) oppure -g (Genotype filter).

filter_expr="EVA_4 - danio_rerio"
echo "Esecuzione filtro: $filter_expr" > result_danio.txt
start=$(date +%s%N)
vcffilter -f "EVA_4" data/danio_rerio.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt


filter_expr="E_Multiple_observations - danio_rerio"
echo "Esecuzione filtro: $filter_expr" >> result_danio.txt
start=$(date +%s%N)
vcffilter -f "E_Multiple_observations" data/danio_rerio.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt


filter_expr="TSA=SNV - danio_rerio"
echo "Esecuzione filtro: $filter_expr" >> result_danio.txt
start=$(date +%s%N)
vcffilter -f "TSA = SNV" data/danio_rerio.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt

#CREAZIONE FILE
#echo "INIZIO CREAZIONE FILE"
#grep -v "^#" data/danio_rerio.vcf | awk '{print $1"\t"($2-1)"\t"$2"\t"$2}' > myannotations.bed
#vcfannotate -b myannotations.bed -k mypos data/danio_rerio.vcf > data/danio_rerio_annotated.vcf
#sed '/^##INFO=<ID=mypos,/ s/Type=String/Type=Integer/' data/danio_rerio_annotated.vcf > data/danio_rerio_annotated_new.vcf
#sed -E 's/(mypos=[0-9]+)(:[0-9]+)+/\1/g' data/danio_rerio_annotated_new.vcf > data/danio_rerio_annotated_new_fixed.vcf 
#echo "FILE CREATO"

filter_expr="POS > 1780000 - danio_rerio"
echo "Esecuzione filtro: $filter_expr" >> result_danio.txt
start=$(date +%s%N)
vcffilter -f "mypos > 1780000" data/danio_rerio_annotated_new_fixed.vcf > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt

filter_expr="POS > 1780000 & POS < 1800000 - danio_rerio"
echo "Esecuzione filtro: $filter_expr" >> result_danio.txt
start=$(date +%s%N)
vcffilter -f "mypos > 1780000 & mypos < 1800000" data/danio_rerio_annotated_new_fixed.vcf > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt

filter_expr="! E_Multiple_observations and TSA=SNV - danio_rerio"
echo "Esecuzione filtro: $filter_expr" >> result_danio.txt
start=$(date +%s%N)
vcffilter -f "! E_Multiple_observations & TSA = SNV" data/danio_rerio.vcf  > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt

filter_expr="POS > 1780000 & E_Multiple_observations - danio_rerio"
echo "Esecuzione filtro: $filter_expr" >> result_danio.txt
start=$(date +%s%N)
vcffilter -f "mypos > 1780000 & E_Multiple_observations" data/danio_rerio_annotated_new_fixed.vcf > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt

filter_expr="POS > 1780000 & POS < 1800000 & E_Multiple_observations - danio_rerio"
echo "Esecuzione filtro: $filter_expr" >> result_danio.txt
start=$(date +%s%N)
vcffilter -f "mypos > 1780000 & mypos < 1800000 & E_Multiple_observations" data/danio_rerio_annotated_new_fixed.vcf > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt

filter_expr="POS > 1780000 & TSA=SNV - danio_rerio"
echo "Esecuzione filtro: $filter_expr" >> result_danio.txt
start=$(date +%s%N)
vcffilter -f "mypos > 1780000 & TSA = SNV" data/danio_rerio_annotated_new_fixed.vcf > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt

filter_expr="POS > 1780000 & POS < 1800000 & TSA=SNV - danio_rerio"
echo "Esecuzione filtro: $filter_expr" >> result_danio.txt
start=$(date +%s%N)
vcffilter -f "mypos > 1780000 & mypos < 1800000 & TSA = SNV" data/danio_rerio_annotated_new_fixed.vcf > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt


filter_expr="EVA_4 & POS > 1780000 & TSA=SNV - danio_rerio"
echo "Esecuzione filtro: $filter_expr" >> result_danio.txt
start=$(date +%s%N)
vcffilter -f "EVA_4 & mypos > 1780000 & TSA = SNV" data/danio_rerio_annotated_new_fixed.vcf > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt

filter_expr="EVA_4 & POS > 1780000 & POS < 1800000 & TSA=SNV - danio_rerio"
echo "Esecuzione filtro: $filter_expr" >> result_danio.txt
start=$(date +%s%N)
vcffilter -f "EVA_4 & mypos > 1780000 & mypos < 1800000 & TSA = SNV" data/danio_rerio_annotated_new_fixed.vcf > output.vcf
end=$(date +%s%N)
runtime_ns=$((end - start))
runtime_ms=$((runtime_ns / 1000000))
echo "Tempo: ${runtime_ms} ms" >> result_danio.txt
echo "" >> result_danio.txt