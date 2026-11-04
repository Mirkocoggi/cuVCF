"""More round trip checks, on files built from TestReconstructor/roundtrip_cases.vcf:
  - the output does not depend on the number of threads;
  - CPU and GPU give the same tables and the same VCF;
  - a file with no variants, a file with no samples, a file larger than one GPU chunk;
  - the CLI (-o) writes the same file as reconstruct_vcf().

Run from the repo root, after a CMake build into build/:
    PYTHONPATH=build python TestReconstructor/test_backends.py
"""
import filecmp
import importlib
import os
import subprocess
import sys
import tempfile

import numpy as np

from test_roundtrip import CASES, compare

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHUNK = 10000  # CHUNK_SIZE of the GPU reconstructor
GETTERS = {
    "DF1": ("get_var_columns_data", "var_columns"),
    "DF2": ("get_alt_columns_data", "alt_columns"),
    "DF3": ("get_sample_columns_data", "samp_columns"),
    "DF4": ("get_alt_format_data", "alt_sample"),
}

failures = []


def check(ok, what):
    print(f"  {'OK  ' if ok else 'FAIL'} {what}")
    if not ok:
        failures.append(what)


def load(name):
    try:
        return importlib.import_module(name)
    except ImportError as e:
        print(f"{name}: skipped ({e})")
        return None


def parse(module, path, threads=4):
    res = module.vcf_parsed()
    res.run(path, threads)
    return res


def reconstruct(module, src, out, threads=4):
    parse(module, src, threads).reconstruct_vcf(out, threads)


def write_variants(path, records):
    header = [l for l in open(CASES) if l.startswith("#")]
    with open(path, "w") as f:
        f.writelines(header)
        f.writelines(records)


def make_files(tmp):
    """Build the derived inputs next to the output files."""
    records = [l for l in open(CASES) if not l.startswith("#")]
    files = {"cases": CASES}

    files["no_variants"] = os.path.join(tmp, "no_variants.vcf")
    write_variants(files["no_variants"], [])

    files["no_samples"] = os.path.join(tmp, "no_samples.vcf")
    with open(files["no_samples"], "w") as f:
        for l in open(CASES):
            f.write("\t".join(l.rstrip("\n").split("\t")[:8]) + "\n")

    # 2.5 GPU chunks: the hand-written records repeated, POS made unique
    files["multi_chunk"] = os.path.join(tmp, "multi_chunk.vcf")
    big = []
    for i in range(CHUNK * 5 // 2):
        cols = records[i % len(records)].split("\t")
        cols[1] = str(i + 1)
        big.append("\t".join(cols))
    write_variants(files["multi_chunk"], big)
    return files


def tables(module, res):
    return {df: getattr(module, getter)(getattr(res, attr)) for df, (getter, attr) in GETTERS.items()}


def same_tables(a, b):
    if a.keys() != b.keys():
        return False
    for df in a:
        if a[df].keys() != b[df].keys():
            return False
        for col in a[df]:
            x, y = a[df][col], b[df][col]
            if isinstance(x, np.ndarray) or isinstance(y, np.ndarray):
                if not np.array_equal(np.asarray(x), np.asarray(y)):
                    return False
            elif list(x) != list(y):
                return False
    return True


def main():
    backends = {n: m for n in ("CPUParser", "GPUParser") if (m := load(n))}
    with tempfile.TemporaryDirectory() as tmp:
        files = make_files(tmp)

        for name, module in backends.items():
            print(name)
            for case, path in files.items():
                out = os.path.join(tmp, f"{name}_{case}.vcf")
                reconstruct(module, path, out)
                errors = compare(out, path)
                for e in errors[:5]:
                    print(f"      {e}")
                check(not errors, f"{case}: field by field equal to the input")

            # Same bytes with 1, 3 and 16 threads: chunks of different sizes, empty ones included
            ref = os.path.join(tmp, f"{name}_t1.vcf")
            reconstruct(module, files["multi_chunk"], ref, 1)
            for threads in (3, 16):
                out = os.path.join(tmp, f"{name}_t{threads}.vcf")
                reconstruct(module, files["multi_chunk"], out, threads)
                check(filecmp.cmp(ref, out, shallow=False), f"multi_chunk: same output with 1 and {threads} threads")

            # The CLI writes what reconstruct_vcf() writes
            cli = os.path.join(REPO, "build", "VCFparser_cpu" if name == "CPUParser" else "VCFparser_gpu")
            if os.path.exists(cli):
                out = os.path.join(tmp, f"{name}_cli.vcf")
                subprocess.run([cli, "-v", CASES, "-t", "4", "-o", out],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                check(filecmp.cmp(os.path.join(tmp, f"{name}_cases.vcf"), out, shallow=False),
                      "CLI -o: same file as reconstruct_vcf()")

        if len(backends) == 2:
            cpu, gpu = backends["CPUParser"], backends["GPUParser"]
            print("CPU vs GPU")
            for case in ("cases", "multi_chunk"):
                same = same_tables(tables(cpu, parse(cpu, files[case])), tables(gpu, parse(gpu, files[case])))
                check(same, f"{case}: same DF1-DF4")
                # Same VCF except QUAL and float text: compare the GPU output with the CPU one, by value
                errors = compare(os.path.join(tmp, f"GPUParser_{case}.vcf"), os.path.join(tmp, f"CPUParser_{case}.vcf"))
                check(not errors, f"{case}: GPU VCF equal to the CPU one, field by field")

    print("ALL OK" if not failures else f"{len(failures)} CHECKS FAILED")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
