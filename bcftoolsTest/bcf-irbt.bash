#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f ../data/IRBT.vcf.gz ]]; then
    echo "⇢ Compressing and indexing IRBT.vcf ..."
    bgzip -c ../data/IRBT.vcf > ../data/IRBT.vcf.gz
    tabix -p vcf ../data/IRBT.vcf.gz
fi

VCF=../data/IRBT.vcf.gz
RES=../result/Bcftools_result_irbt.txt
mkdir -p ../result
: > "$RES"

run () {
    local label=$1 expr=$2
    echo "Esecuzione filtro: $label" | tee -a "$RES"

    local t0=$(date +%s%N)
    bcftools view -i "$expr" -Ou "$VCF" > /dev/null
    local ns=$(( $(date +%s%N) - t0 ))
    printf "Tempo: %.9f s\n\n" "$(awk "BEGIN{print $ns/1e9}")" \
        | tee -a "$RES"
}

run "AC > 8"                             'INFO/AC>8'
run "AC > 8 && AF > 0.5"                 'INFO/AC>8 && INFO/AF>0.5'

run "GT = 1|1"                           'FMT/GT="1|1"'
run "AD[0] > 3  (sample0)"               'FMT/AD[0:0]>3'
run "AD[0] > 3 && AD[1] < 10  (sample0)" 'FMT/AD[0:0]>3 && FMT/AD[0:1]<10'

run "GT=1|1 && AC>8"                     'FMT/GT="1|1" && INFO/AC>8'
run "GT=1|1 && AC>8 && AF>0.5"           'FMT/GT="1|1" && INFO/AC>8 && INFO/AF>0.5'

run "AD0>3(sample0) && AC>8"             'FMT/AD[0:0]>3 && INFO/AC>8'
run "AC>8 && AD0>3 && AD1<10 (sample0)"  'INFO/AC>8 && FMT/AD[0:0]>3 && FMT/AD[0:1]<10'

run "AC>8 && AF>0.5 && AD0>3 (sample0)"  'INFO/AC>8 && INFO/AF>0.5 && FMT/AD[0:0]>3'
run "AC>8 && AF>0.5 && AD0>3 && AD1<10 (sample0)" \
    'INFO/AC>8 && INFO/AF>0.5 && FMT/AD[0:0]>3 && FMT/AD[0:1]<10'

run "AD0>3 && GT=1|1 (sample0)"          'FMT/AD[0:0]>3 && FMT/GT="1|1"'
run "AD0>3 && AD1<10 && GT=1|1 (sample0)"\
    'FMT/AD[0:0]>3 && FMT/AD[0:1]<10 && FMT/GT="1|1"'

run "AC>8 && GT=1|1 && AD0>3 (sample0)"  'INFO/AC>8 && FMT/GT="1|1" && FMT/AD[0:0]>3'
run "AC>8 && AD0>3 && AD1<10 && GT=1|1 (sample0)" \
    'INFO/AC>8 && FMT/AD[0:0]>3 && FMT/AD[0:1]<10 && FMT/GT="1|1"'

run "AC>8 && AF>0.5 && GT=1|1 && AD0>3 (sample0)" \
    'INFO/AC>8 && INFO/AF>0.5 && FMT/GT="1|1" && FMT/AD[0:0]>3'
run "AC>8 && AF>0.5 && GT=1|1 && AD0>3 && AD1<10 (sample0)" \
    'INFO/AC>8 && INFO/AF>0.5 && FMT/GT="1|1" && FMT/AD[0:0]>3 && FMT/AD[0:1]<10'

echo "Risultati salvati in $RES"
