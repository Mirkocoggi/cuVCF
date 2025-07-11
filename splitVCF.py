#!/usr/bin/env python3
"""
split_vcf.py  ────────────────────────────────────────────────────────────────
Crea due file VCF “metà & metà” a partire da un VCF originale,
mantenendo nei due file le **prime 254 righe di header**.

USO:
    python split_vcf.py data/chrx.vcf  out1.vcf  out2.vcf

Se i nomi di output non sono indicati verranno creati:
    <input>_part1.vcf  e  <input>_part2.vcf
───────────────────────────────────────────────────────────────────────────────
"""
import sys
from pathlib import Path
from math import ceil

def main(src, dst1=None, dst2=None):
    src_path = Path(src)
    if not src_path.is_file():
        sys.exit(f"File non trovato: {src}")

    # nomi di output di default
    if dst1 is None:
        dst1 = src_path.with_stem(src_path.stem + "_part1").with_suffix(".vcf")
    if dst2 is None:
        dst2 = src_path.with_stem(src_path.stem + "_part2").with_suffix(".vcf")

    # ── 1. leggi header (prime 254 righe) ────────────────────────────────
    header_lines = []
    variant_lines = []
    with src_path.open() as f:
        for _ in range(254):
            header_lines.append(f.readline())
        # le restanti righe sono varianti
        variant_lines = f.readlines()

    # ── 2. calcola il punto di split ────────────────────────────────────
    mid = len(variant_lines) // 2            # metà esatta
    first_half  = variant_lines[:mid]
    second_half = variant_lines[mid:]

    # ── 3. scrivi i due file ────────────────────────────────────────────
    for dst, body in ((dst1, first_half), (dst2, second_half)):
        with Path(dst).open("w") as out:
            out.writelines(header_lines)
            out.writelines(body)
        print(f"Creato: {dst}  ({len(body):,} record)")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Uso: python split_vcf.py <input.vcf> [output1.vcf output2.vcf]")
    main(*sys.argv[1:])
