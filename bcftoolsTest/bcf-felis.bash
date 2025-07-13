#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 0.  Prepara il VCF compresso e indicizzato (solo la prima volta)
###############################################################################
if [[ ! -f ../data/felis_catus.vcf.gz ]]; then
    echo "⇢ Compressing and indexing felis_catus.vcf ..."
    bgzip -c ../data/felis_catus.vcf > ../data/felis_catus.vcf.gz
    tabix -p vcf ../data/felis_catus.vcf.gz
fi

VCF=../data/felis_catus.vcf.gz
RES=../result/Bcftools_result_felis.txt
mkdir -p ../result
: > "$RES"                     # svuota file risultato

###############################################################################
# Funzione di benchmark
###############################################################################
run () {
    local label=$1 expr=$2
    echo "Esecuzione filtro: $label" | tee -a "$RES"

    local t0=$(date +%s%N)              # timestamp iniziale in ns
    bcftools view -i "$expr" -Ou "$VCF" > /dev/null
    local dt_ns=$(( $(date +%s%N) - t0 ))

    # secondi con 9 cifre decimali
    local dt_s
    dt_s=$(awk 'BEGIN{printf "%.9f", ARGV[1]/1e9}' "$dt_ns")

    echo "Tempo: ${dt_s} s" | tee -a "$RES"
    echo >> "$RES"
}

###############################################################################
# 1.  Filtri INFO semplici
###############################################################################
run "EVA_4 - felis_catus"                'INFO/EVA_4=1'
run "E_Multiple_observations - felis_catus" \
                                          'INFO/E_Multiple_observations=1'
run "TSA=SNV - felis_catus"              'INFO/TSA=="SNV"'

###############################################################################
# 2.  Filtri basati su coordinate POS
###############################################################################
run "POS > 140680000"                    'POS>140680000'
run "140.68M < POS < 160.68M"            'POS>140680000 && POS<160680000'

###############################################################################
# 3.  Combinazioni INFO + POS
###############################################################################
run "!E_Multiple_observations && TSA=SNV" \
    'INFO/E_Multiple_observations=0 && INFO/TSA=="SNV"'

run "POS>140.68M && E_Multiple_observations" \
    'POS>140680000 && INFO/E_Multiple_observations=1'

run "140.68M<POS<160.68M && E_Multiple_observations" \
    'POS>140680000 && POS<160680000 && INFO/E_Multiple_observations=1'

run "POS>140.68M && TSA=SNV" \
    'POS>140680000 && INFO/TSA=="SNV"'

run "140.68M<POS<160.68M && TSA=SNV" \
    'POS>140680000 && POS<160680000 && INFO/TSA=="SNV"'

run "EVA_4 && POS>140.68M && TSA=SNV" \
    'INFO/EVA_4=1 && POS>140680000 && INFO/TSA=="SNV"'

run "EVA_4 && 140.68M<POS<160.68M && TSA=SNV" \
    'INFO/EVA_4=1 && POS>140680000 && POS<160680000 && INFO/TSA=="SNV"'

echo "Risultati salvati in $RES"
