#!/usr/bin/env python3

import os, subprocess, time, pathlib, numpy as np
from cyvcf2 import VCF

RAW_VCF = "../data/IRBT.vcf"
GZ_VCF  = RAW_VCF + ".gz"
RES_TXT = "../result/cyvcf2_irbt_times.txt"
pathlib.Path("../result").mkdir(exist_ok=True)

if not os.path.exists(GZ_VCF):
    print("⇢ bgzip IRBT.vcf …")
    with open(GZ_VCF, "wb") as gz:
        subprocess.check_call(["bgzip", "-c", RAW_VCF], stdout=gz)

if not os.path.exists(GZ_VCF + ".tbi"):
    print("⇢ tabix -p vcf IRBT.vcf.gz …")
    subprocess.check_call(["tabix", "-p", "vcf", GZ_VCF])

VCF_IN = GZ_VCF

def _parse_info_numeric(val, cast):
    """
    Restituisce una lista di valori cast-ati:
      • val = None           → []
      • val = int/float      → [cast(val)]
      • val = '3,15,8'       → [3,15,8]
      • val = b'3,15,8'      → idem
    """
    if val is None:
        return []
    if isinstance(val, (int, float)):
        return [cast(val)]
    if isinstance(val, (bytes, str)):
        s = val.decode() if isinstance(val, bytes) else val
        return [cast(x) for x in s.split(",") if x != "."]

    try:
        return [cast(x) for x in val]
    except TypeError:
        return []

def ac_gt8(v):
    vals = _parse_info_numeric(v.INFO.get("AC"), int)
    return vals and max(vals) > 8

def af_gt05(v):
    vals = _parse_info_numeric(v.INFO.get("AF"), float)
    return vals and max(vals) > 0.5

def gt_1_1(v):
    g = v.genotypes[0]              
    return g[0] == 1 and g[1] == 1

def ad_values(v):
    """
    Ritorna (AD0, AD1) per il sample 0 oppure (None, None)
    AD è un ndarray shape (n_samples, n_subfields)
    """
    ad = v.format("AD")
    if ad is None or ad.size == 0:
        return None, None
    ad0 = int(ad[0][0])
    ad1 = int(ad[0][1]) if ad.shape[1] > 1 else None
    return ad0, ad1

def ad0_gt3(v):
    ad0, _ = ad_values(v)
    return ad0 is not None and ad0 > 3

def ad0_gt3_ad1_lt10(v):
    ad0, ad1 = ad_values(v)
    return (ad0 is not None and ad0 > 3) and (ad1 is not None and ad1 < 10)

TESTS = [

    ("AC > 8",                                    ac_gt8),
    ("AC > 8 && AF > 0.5",                        lambda v: ac_gt8(v) and af_gt05(v)),

    ("GT = 1|1",                                  gt_1_1),
    ("AD0 > 3 (sample0)",                         ad0_gt3),
    ("AD0 > 3 && AD1 < 10 (sample0)",             ad0_gt3_ad1_lt10),

    ("GT=1|1 && AC>8",                            lambda v: gt_1_1(v) and ac_gt8(v)),
    ("GT=1|1 && AC>8 && AF>0.5",
                                                 lambda v: gt_1_1(v) and ac_gt8(v) and af_gt05(v)),

    ("AD0>3 && AC>8 (sample0)",                   lambda v: ad0_gt3(v) and ac_gt8(v)),
    ("AC>8 && AD0>3 && AD1<10 (sample0)",
                                                 lambda v: ac_gt8(v) and ad0_gt3_ad1_lt10(v)),

    ("AC>8 && AF>0.5 && AD0>3 (sample0)",
                                                 lambda v: ac_gt8(v) and af_gt05(v) and ad0_gt3(v)),
    ("AC>8 && AF>0.5 && AD0>3 && AD1<10 (sample0)",
                                                 lambda v: ac_gt8(v) and af_gt05(v) and ad0_gt3_ad1_lt10(v)),

    ("AD0>3 && GT=1|1 (sample0)",                 lambda v: ad0_gt3(v) and gt_1_1(v)),
    ("AD0>3 && AD1<10 && GT=1|1 (sample0)",
                                                 lambda v: ad0_gt3_ad1_lt10(v) and gt_1_1(v)),

    ("AC>8 && GT=1|1 && AD0>3 (sample0)",
                                                 lambda v: ac_gt8(v) and gt_1_1(v) and ad0_gt3(v)),
    ("AC>8 && AD0>3 && AD1<10 && GT=1|1 (sample0)",
                                                 lambda v: ac_gt8(v) and gt_1_1(v) and ad0_gt3_ad1_lt10(v)),

    ("AC>8 && AF>0.5 && GT=1|1 && AD0>3 (sample0)",
                                                 lambda v: ac_gt8(v) and af_gt05(v) and gt_1_1(v) and ad0_gt3(v)),
    ("AC>8 && AF>0.5 && GT=1|1 && AD0>3 && AD1<10 (sample0)",
                                                 lambda v: ac_gt8(v) and af_gt05(v) and gt_1_1(v) and ad0_gt3_ad1_lt10(v)),
]

def bench(label, predicate):
    rdr = VCF(VCF_IN)
    t0 = time.perf_counter()
    for rec in rdr:
        _ = predicate(rec)
    rdr.close()
    return time.perf_counter() - t0

def main():
    with open(RES_TXT, "w") as fh:
        for lab, fn in TESTS:
            dt = bench(lab, fn)
            entry = f"Esecuzione filtro: {lab}\nTempo: {dt:.9f} s\n\n"
            fh.write(entry)
            print(entry, end="")

    print(f"\nTempistiche scritte in {RES_TXT}")

if __name__ == "__main__":
    main()
