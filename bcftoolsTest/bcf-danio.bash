#!/usr/bin/env python3
# cyvcf2 ≥ 0.30
import os, subprocess, time, pathlib
from cyvcf2 import VCF

RAW_VCF = "../data/danio_rerio.vcf"
GZ_VCF  = RAW_VCF + ".gz"
RES_TXT = "../result/cyvcf2_danio_times.txt"
pathlib.Path("../result").mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# 1) bgzip + tabix (solo la prima volta)
# ---------------------------------------------------------------------------
if not os.path.exists(GZ_VCF):
    print("⇢ bgzip danio_rerio.vcf …")
    with open(GZ_VCF, "wb") as gz:
        subprocess.check_call(["bgzip", "-c", RAW_VCF], stdout=gz)

if not os.path.exists(GZ_VCF + ".tbi"):
    print("⇢ tabix -p vcf danio_rerio.vcf.gz …")
    subprocess.check_call(["tabix", "-p", "vcf", GZ_VCF])

VCF_IN = GZ_VCF            # useremo sempre il file bgz + indice

# ---------------------------------------------------------------------------
# 2) predicati di filtro
# ---------------------------------------------------------------------------
def eva4(v):    return v.INFO.get("EVA_4") is not None
def emult(v):   return v.INFO.get("E_Multiple_observations") is not None
def tsa_snv(v): return (t := v.INFO.get("TSA")) is not None and t.strip() == "SNV"

def pos_gt(v, thr):          return v.POS > thr
def pos_between(v, lo, hi):  return lo < v.POS < hi

POS_LOW, POS_HIGH = 1_780_000, 1_800_000

TESTS = [
    ("EVA_4 - danio_rerio",                         eva4),
    ("E_Multiple_observations - danio_rerio",       emult),
    ("TSA=SNV - danio_rerio",                       tsa_snv),
    ("POS > 1 780 000",                             lambda v: pos_gt(v, POS_LOW)),
    ("1.78 M < POS < 1.80 M",                       lambda v: pos_between(v, POS_LOW, POS_HIGH)),
    ("!E_Multiple_observations && TSA=SNV",
                                                   lambda v: (not emult(v)) and tsa_snv(v)),
    ("POS>1.78 M && E_Multiple_observations",
                                                   lambda v: pos_gt(v, POS_LOW) and emult(v)),
    ("1.78 M<POS<1.80 M && E_Multiple_observations",
                                                   lambda v: pos_between(v, POS_LOW, POS_HIGH) and emult(v)),
    ("POS>1.78 M && TSA=SNV",
                                                   lambda v: pos_gt(v, POS_LOW) and tsa_snv(v)),
    ("1.78 M<POS<1.80 M && TSA=SNV",
                                                   lambda v: pos_between(v, POS_LOW, POS_HIGH) and tsa_snv(v)),
    ("EVA_4 && POS>1.78 M && TSA=SNV",
                                                   lambda v: eva4(v) and pos_gt(v, POS_LOW) and tsa_snv(v)),
    ("EVA_4 && 1.78 M<POS<1.80 M && TSA=SNV",
                                                   lambda v: eva4(v) and pos_between(v, POS_LOW, POS_HIGH) and tsa_snv(v)),
]

# ---------------------------------------------------------------------------
# 3) benchmark senza scrivere alcun VCF
# ---------------------------------------------------------------------------
def bench(label, pred):
    rdr = VCF(VCF_IN)          # lettura indicizzata: niente warning sui contig
    t0 = time.perf_counter()
    for rec in rdr:
        _ = pred(rec)          # valutiamo il predicato, ignoriamo l’output
    rdr.close()
    return time.perf_counter() - t0

def main():
    with open(RES_TXT, "w") as fh:
        for lab, fn in TESTS:
            t = bench(lab, fn)
            line = f"Esecuzione filtro: {lab}\nTempo: {t:.9f} s\n\n"
            fh.write(line)
            print(line, end="")

    print(f"\nTempistiche scritte in {RES_TXT}")

if __name__ == "__main__":
    main()
