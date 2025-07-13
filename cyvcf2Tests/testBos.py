#!/usr/bin/env python3
# cyvcf2 ≥ 0.30
import os, subprocess, time, pathlib
from cyvcf2 import VCF

RAW_VCF = "../data/bos_taurus.vcf"
GZ_VCF  = RAW_VCF + ".gz"
RES_TXT = "../result/cyvcf2_bos_times.txt"
pathlib.Path("../result").mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# 1) bgzip + tabix (una sola volta)
# ---------------------------------------------------------------------------
if not os.path.exists(GZ_VCF):
    print("⇢ bgzip bos_taurus.vcf …")
    with open(GZ_VCF, "wb") as gz:
        subprocess.check_call(["bgzip", "-c", RAW_VCF], stdout=gz)

if not os.path.exists(GZ_VCF + ".tbi"):
    print("⇢ tabix -p vcf bos_taurus.vcf.gz …")
    subprocess.check_call(["tabix", "-p", "vcf", GZ_VCF])

VCF_IN = GZ_VCF            # useremo sempre il file indicizzato

# ---------------------------------------------------------------------------
# 2) definizione dei filtri
# ---------------------------------------------------------------------------
def eva4(v):    return v.INFO.get("EVA_4") is not None
def emult(v):   return v.INFO.get("E_Multiple_observations") is not None
def tsa_snv(v): return (t := v.INFO.get("TSA")) is not None and t.strip() == "SNV"
def pos_gt(v, thr):              return v.POS > thr
def pos_between(v, lo, hi):      return lo < v.POS < hi

# lista di (etichetta, funzione-filtro)
TESTS = [
    ("EVA_4",                               eva4),
    ("E_Multiple_observations",             emult),
    ("TSA=SNV",                             tsa_snv),
    ("POS>200000",                          lambda v: pos_gt(v, 200_000)),
    ("200k<POS<300k",                       lambda v: pos_between(v, 200_000, 300_000)),
    ("!E_Multiple_observations && TSA=SNV", lambda v: (not emult(v)) and tsa_snv(v)),
    ("POS>200k && E_Multiple_observations", lambda v: pos_gt(v, 200_000) and emult(v)),
    ("200k<POS<300k && E_Multiple_observations",
                                             lambda v: pos_between(v, 200_000, 300_000) and emult(v)),
    ("POS>200k && TSA=SNV",                 lambda v: pos_gt(v, 200_000) and tsa_snv(v)),
    ("200k<POS<300k && TSA=SNV",
                                             lambda v: pos_between(v, 200_000, 300_000) and tsa_snv(v)),
    ("EVA_4 && POS>200k && TSA=SNV",
                                             lambda v: eva4(v) and pos_gt(v, 200_000) and tsa_snv(v)),
    ("EVA_4 && 200k<POS<300k && TSA=SNV",
                                             lambda v: eva4(v) and pos_between(v, 200_000, 300_000) and tsa_snv(v)),
]

# ---------------------------------------------------------------------------
# 3) benchmark
# ---------------------------------------------------------------------------
def benchmark(label, filt_fun):
    rdr = VCF(VCF_IN)               # nessun writer: iteriamo e basta
    t0 = time.perf_counter()
    for rec in rdr:
        _ = filt_fun(rec)           # valutazione filtro
    dt = time.perf_counter() - t0
    rdr.close()
    return dt

def main():
    with open(RES_TXT, "w") as fh:
        for lab, func in TESTS:
            elapsed = benchmark(lab, func)
            line = f"Esecuzione filtro: {lab}\nTempo: {elapsed:.9f} s\n\n"
            fh.write(line)
            print(line, end="")     # stampa anche a video

    print(f"\nTempistiche scritte in {RES_TXT}")

if __name__ == "__main__":
    main()
