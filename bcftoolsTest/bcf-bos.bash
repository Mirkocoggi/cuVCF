#!/usr/bin/env bash
set -euo pipefail
# prepara una sola volta
bgzip -c ../data/bos_taurus.vcf > ../data/bos_taurus.vcf.gz
tabix -p vcf ../data/bos_taurus.vcf.gz

# nel tuo script:
VCF=../data/bos_taurus.vcf.gz     # invece del .vcf

OUT=output.vcf.gz
RES=../result/Bcftools_result_bos.txt
: > "$RES"                  # svuota file risultato

run () {
    local label=$1 expr=$2
    echo "Esecuzione filtro: $label" | tee -a "$RES"

    local t0=$(date +%s%N)                     # nanosecondi iniziali
    bcftools view -i "$expr" -Ou "$VCF" > /dev/null
    local runtime_ns=$(( $(date +%s%N) - t0 )) # differenza in ns

    # converti in secondi con 9 cifre decimali (precisione nanosecondo)
    local runtime_s
    runtime_s=$(awk 'BEGIN{printf "%.9f", ARGV[1]/1e9}' "$runtime_ns")

    echo "Tempo: ${runtime_s} s" | tee -a "$RES"
    echo >> "$RES"                             # riga vuota separatrice
}


# … parti iniziali invariate …

run "EVA_4"                       'INFO/EVA_4=1'
run "E_Multiple_observations"     'INFO/E_Multiple_observations=1'
run "TSA=SNV"                     'INFO/TSA=="SNV"'

# POS-based
run "POS>200000"                  'POS>200000'
run "200k<POS<300k"               'POS>200000 && POS<300000'

# combinazioni
run "!E_Multiple_observations && TSA=SNV" \
    'INFO/E_Multiple_observations=0 && INFO/TSA=="SNV"'

run "POS>200k && E_Multiple_observations" \
    'POS>200000 && INFO/E_Multiple_observations=1'

run "200k<POS<300k && E_Multiple_observations" \
    'POS>200000 && POS<300000 && INFO/E_Multiple_observations=1'

run "POS>200k && TSA=SNV" \
    'POS>200000 && INFO/TSA=="SNV"'

run "200k<POS<300k && TSA=SNV" \
    'POS>200000 && POS<300000 && INFO/TSA=="SNV"'

run "EVA_4 && POS>200k && TSA=SNV" \
    'INFO/EVA_4=1 && POS>200000 && INFO/TSA=="SNV"'

run "EVA_4 && 200k<POS<300k && TSA=SNV" \
    'INFO/EVA_4=1 && POS>200000 && POS<300000 && INFO/TSA=="SNV"'

echo "Risultati salvati in $RES"
