"""Round trip check: VCF -> cuVCF -> reconstructed VCF, compared field by field.

TestReconstructor/roundtrip_cases.vcf is written by hand, so it is the expected output too.
The comparison is on values, not on text. These differences are accepted:
  - order of the INFO keys and of the FORMAT keys;
  - float text and fp16 rounding (1e-3 relative tolerance);
  - a missing field written as "." or left out (the VCF spec allows both).

Run from the repo root, after a CMake build into build/:
    PYTHONPATH=build python TestReconstructor/test_roundtrip.py
"""
import importlib
import os
import re
import sys
import tempfile

CASES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "roundtrip_cases.vcf")


def read_vcf(path):
    """Return ({"INFO": {id: type}, "FORMAT": {id: type}}, header lines, records)."""
    types = {"INFO": {}, "FORMAT": {}}
    header, records = [], []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header.append(line)
                m = re.match(r"##(INFO|FORMAT)=<ID=([^,]+),.*Type=(\w+)", line)
                if m:
                    types[m.group(1)][m.group(2)] = m.group(3)
            else:
                records.append(line.split("\t"))
    return types, header, records


def number(text, kind):
    if text == ".":
        return None
    return float(text) if kind == "Float" else text


def values(text, kind):
    """Split a field: all-missing arrays count as a missing field."""
    vals = [number(v, kind) for v in text.split(",")]
    return None if all(v is None for v in vals) else vals


def parse_info(text, types):
    fields = {}
    if text == ".":
        return fields
    for item in text.split(";"):
        key, _, val = item.partition("=")
        parsed = True if not val else values(val, types.get(key, "String"))
        if parsed is not None:
            fields[key] = parsed
    return fields


def parse_cell(keys, cell, types):
    fields = {}
    for key, val in zip(keys, cell.split(":")):
        parsed = val if key == "GT" else values(val, types.get(key, "String"))
        if parsed is not None:
            fields[key] = parsed
    return fields


def same(a, b):
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, float) and isinstance(b, float):
        return abs(a - b) <= 1e-3 * max(abs(a), 1.0)
    return a == b


def compare_fields(where, expected, got, errors):
    for key in sorted(set(expected) | set(got)):
        if key not in got:
            errors.append(f"{where} {key}: expected {expected[key]}, missing in the output")
        elif key not in expected:
            errors.append(f"{where} {key}: not in the input, got {got[key]}")
        elif not same(expected[key], got[key]):
            errors.append(f"{where} {key}: expected {expected[key]}, got {got[key]}")


def compare(out_path, expected_path=CASES):
    types, header_in, rows_in = read_vcf(expected_path)
    _, header_out, rows_out = read_vcf(out_path)
    errors = []
    if header_in != header_out:
        errors.append("header differs")
    if len(rows_in) != len(rows_out):
        return errors + [f"expected {len(rows_in)} records, got {len(rows_out)}"]

    for exp, got in zip(rows_in, rows_out):
        where = f"{exp[0]}:{exp[1]}"
        for col, name in ((0, "CHROM"), (1, "POS"), (2, "ID"), (3, "REF"), (4, "ALT"), (6, "FILTER")):
            if exp[col] != got[col]:
                errors.append(f"{where} {name}: expected {exp[col]}, got {got[col]}")
        if not same(number(exp[5], "Float"), number(got[5], "Float")):
            errors.append(f"{where} QUAL: expected {exp[5]}, got {got[5]}")
        compare_fields(f"{where} INFO", parse_info(exp[7], types["INFO"]),
                       parse_info(got[7], types["INFO"]), errors)
        if len(exp) != len(got):
            errors.append(f"{where}: expected {len(exp)} columns, got {len(got)}")
            continue
        for s in range(9, len(exp)):
            compare_fields(f"{where} S{s - 8}",
                           parse_cell(exp[8].split(":"), exp[s], types["FORMAT"]),
                           parse_cell(got[8].split(":"), got[s], types["FORMAT"]), errors)
    return errors


def check_backend(name):
    try:
        module = importlib.import_module(name)
    except ImportError as e:
        print(f"{name}: skipped ({e})")
        return 0

    with tempfile.TemporaryDirectory() as tmp:
        out_path = os.path.join(tmp, "out.vcf")
        res = module.vcf_parsed()
        res.run(CASES, 4)
        res.reconstruct_vcf(out_path, 4)
        errors = compare(out_path)

        try:
            res.reconstruct_vcf(os.path.join(tmp, "missing_dir", "out.vcf"), 4)
            errors.append("unwritable path: no error raised")
        except RuntimeError:
            pass

        if os.path.exists("/dev/full"):  # Linux: every write fails with "no space left"
            try:
                res.reconstruct_vcf("/dev/full", 4)
                errors.append("disk full (/dev/full): no error raised")
            except RuntimeError:
                pass

    for e in errors:
        print(f"  {name} FAIL {e}")
    print(f"{name}: {'OK' if not errors else str(len(errors)) + ' errors'}")
    return len(errors)


def main():
    total = sum(check_backend(name) for name in ("CPUParser", "GPUParser"))
    sys.exit(1 if total else 0)


if __name__ == "__main__":
    main()
