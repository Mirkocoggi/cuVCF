#!/usr/bin/env python3

import os, subprocess, time, pathlib
from cyvcf2 import VCF

RAW_VCF = "../data/felis_catus.vcf"
GZ_VCF  = RAW_VCF + ".gz"
RES_TXT = "../result/cyvcf2_felis_times.txt"
pathlib.Path("../result").mkdir(exist_ok=True)

if not os.path.exists(GZ_VCF):
    print("⇢ bgzip felis_catus.vcf …")
    with open(GZ_VCF, "wb") as gz:
        subprocess.check_call(["bgzip", "-c", RAW_VCF], stdout=gz)
if not os.path.exists(GZ_VCF + ".tbi"):
    print("⇢ tabix -p vcf felis_catus.vcf.gz …")
    subprocess.check_call(["tabix", "-p", "vcf", GZ_VCF])

VCF_IN = GZ_VCF

def eva4(v):    return v.INFO.get("EVA_4") is not None
def emult(v):   return v.INFO.get("E_Multiple_observations") is not None
def tsa_snv(v): return (t := v.INFO.get("TSA")) is not None and t.strip() == "SNV"
def pos_gt(v,t):               return v.POS > t
def pos_between(v,lo,hi):      return lo < v.POS < hi

POS_LOW, POS_HIGH = 140_680_000, 160_680_000

TESTS = [
    ("EVA_4 - felis_catus",                        eva4),
    ("E_Multiple_observations - felis_catus",      emult),
    ("TSA=SNV - felis_catus",                      tsa_snv),
    ("POS > 140680000",                            lambda v: pos_gt(v, POS_LOW)),
    ("140.68M < POS < 160.68M",                    lambda v: pos_between(v, POS_LOW, POS_HIGH)),
    ("!E_Multiple_observations && TSA=SNV",        lambda v: (not emult(v)) and tsa_snv(v)),
    ("POS>140.68M && E_Multiple_observations",     lambda v: pos_gt(v, POS_LOW) and emult(v)),
    ("140.68M<POS<160.68M && E_Multiple_observations",
                                                  lambda v: pos_between(v, POS_LOW, POS_HIGH) and emult(v)),
    ("POS>140.68M && TSA=SNV",                    lambda v: pos_gt(v, POS_LOW) and tsa_snv(v)),
    ("140.68M<POS<160.68M && TSA=SNV",
                                                  lambda v: pos_between(v, POS_LOW, POS_HIGH) and tsa_snv(v)),
    ("EVA_4 && POS>140.68M && TSA=SNV",
                                                  lambda v: eva4(v) and pos_gt(v, POS_LOW) and tsa_snv(v)),
    ("EVA_4 && 140.68M<POS<160.68M && TSA=SNV",
                                                  lambda v: eva4(v) and pos_between(v, POS_LOW, POS_HIGH) and tsa_snv(v)),
]

def bench(label, filt):
    rdr = VCF(VCF_IN)
    t0 = time.perf_counter()
    for rec in rdr:
        _ = filt(rec)
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
