"""Checks reconstruct_vcf(..., variants=...): writing only a selection of variants.

The expected output of a selection is the header and the selected records of the full
reconstruction (variants=None), byte for byte. The input files are the ones of
test_backends.py, built from roundtrip_cases.vcf (multi_chunk spans 2.5 GPU chunks).

Run from the repo root, after a CMake build into build/:
    PYTHONPATH=build python TestReconstructor/test_filter.py
"""
import os
import sys
import tempfile

import numpy as np
import pandas as pd

from test_backends import CHUNK, load, make_files
from test_roundtrip import compare, read_vcf

failures = []


def check(ok, what):
    print(f"  {'OK  ' if ok else 'FAIL'} {what}")
    if not ok:
        failures.append(what)


def split(path):
    with open(path) as f:
        lines = f.read().splitlines()
    return [l for l in lines if l.startswith("#")], [l for l in lines if not l.startswith("#")]


def info_values(cols, key):
    for item in cols[7].split(";"):
        k, _, v = item.partition("=")
        if k == key:
            return [float(x) for x in v.split(",") if x != "."]
    return []


def text_filter(path, keep):
    """Indices of the records of path for which keep(columns) holds, read from the VCF text."""
    _, _, rows = read_vcf(path)
    return [i for i, cols in enumerate(rows) if keep(cols)]


def selections(n):
    """Name -> var_number values to write, for a file of n variants."""
    sel = {
        "empty": [],
        "first variant": [0],
        "last variant": [n - 1],
        "all variants": list(range(n)),
        "one every 3": list(range(0, n, 3)),
        "random": sorted(np.random.default_rng(0).choice(n, size=min(1000, n), replace=False).tolist()),
    }
    if n > CHUNK + 10:
        sel["block across two GPU chunks"] = list(range(CHUNK - 10, CHUNK + 10))
    return sel


def check_file(module, case, path, tmp):
    res = module.vcf_parsed()
    res.run(path, 4)
    full = os.path.join(tmp, f"{module.__name__}_{case}_full.vcf")
    res.reconstruct_vcf(full, 4)
    header, records = split(full)
    n = len(records)

    def subset(name, variants):
        out = os.path.join(tmp, f"{module.__name__}_{case}_{name.replace(' ', '_')}.vcf")
        res.reconstruct_vcf(out, 4, variants=variants)
        return out

    for name, variants in selections(n).items():
        got = split(subset(name, variants))
        check(got == (header, [records[i] for i in variants]), f"{case}: {name} ({len(variants)} of {n})")

    # Same file as variants=None when every variant is selected
    with open(full) as a, open(subset("all again", np.arange(n))) as b:
        check(a.read() == b.read(), f"{case}: all variants, same file as variants=None")

    # Unsorted, with duplicates, and with other integer types: same as the sorted selection
    every3 = list(range(0, n, 3))
    expected = (header, [records[i] for i in every3])
    shuffled = np.random.default_rng(1).permutation(every3 + every3[:5])
    check(split(subset("unsorted", shuffled)) == expected, f"{case}: unsorted with duplicates")
    for dtype in (np.int32, np.int64, np.uint32):
        check(split(subset(f"dtype {dtype.__name__}", np.array(every3, dtype=dtype))) == expected,
              f"{case}: {dtype.__name__} array")

    # Values that are not a var_number
    for bad in ([n], [-1], [[0, 1]]):
        try:
            res.reconstruct_vcf(os.path.join(tmp, "bad.vcf"), 4, variants=bad)
            check(False, f"{case}: variants={bad} raises ValueError")
        except ValueError:
            check(True, f"{case}: variants={bad} raises ValueError")

    # Real filters on the tables, checked against the same filter on the VCF text
    df2 = pd.DataFrame(module.get_alt_columns_data(res.alt_columns))
    keep = np.unique(df2.var_id[df2.AC > 1])
    want = text_filter(path, lambda cols: any(v > 1 for v in info_values(cols, "AC")))
    got = split(subset("filter AC", keep))
    check(list(keep) == want and got == (header, [records[i] for i in want]),
          f"{case}: DF2 filter AC > 1 ({len(want)} variants)")

    df3 = pd.DataFrame(module.get_sample_columns_data(res.samp_columns))
    keep = np.unique(df3.var_id[df3.GT0 == module.GTMapGlobal["0/1"]])
    want = text_filter(path, lambda cols: "0/1" in [c.split(":")[0] for c in cols[9:]])
    got = split(subset("filter GT", keep))
    check(list(keep) == want and got == (header, [records[i] for i in want]),
          f"{case}: DF3 filter GT = 0/1 ({len(want)} variants)")


def main():
    backends = {n: m for n in ("CPUParser", "GPUParser") if (m := load(n))}
    with tempfile.TemporaryDirectory() as tmp:
        files = make_files(tmp)
        for name, module in backends.items():
            print(name)
            for case in ("cases", "multi_chunk"):
                check_file(module, case, files[case], tmp)

        if len(backends) == 2:
            # CPU and GPU write floats with different text: compare the selections by value
            print("CPU vs GPU")
            for case in ("cases", "multi_chunk"):
                cpu = os.path.join(tmp, f"CPUParser_{case}_random.vcf")
                gpu = os.path.join(tmp, f"GPUParser_{case}_random.vcf")
                errors = compare(gpu, cpu)
                check(not errors, f"{case}: same random selection, field by field")

    print("ALL OK" if not failures else f"{len(failures)} CHECKS FAILED")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
