# cuVCF

**cuVCF** is a **CPU/GPU-accelerated framework** that converts Variant Call Format (**VCF**) files into **columnar DataFrames**, directly usable with **pandas** and **cuDF**.  
Its **hybrid CPU+GPU pipeline** delivers substantial performance speedups over traditional tools.

## Requirements

### Build toolchain

* **g++** with **OpenMP** (`-fopenmp`)
* **CUDA Toolkit 12.x** (Makefile uses change to specify your architecture `-arch=sm_89`)
* **nvcc**
* **python3**, **python3-config**
* **pybind11** headers (`python3 -m pybind11 --includes`)
* **zlib** (`-lz`)
* **Imath** headers (Makefile includes `-I/usr/include/Imath`)
* **C++17**

### Python (runtime)

* **numpy**
* **pandas** (for CPU analytics)
* **cuDF** (optional, for GPU analytics - same sintax of pandas)
* Your module loads as `CPUParser` / `GPUParser` from `bin/`

### Tested environment

* **Ubuntu 22.04.5 LTS**
* **NVIDIA Ada (sm\_89)** class GPU (e.g., L40S).
  For other GPUs, change `-arch=sm_89` in `GPUFLAGS`/`DEBUG_GPUFLAGS`.

---

## Quick Start

### 1) Build

Clone and compile using `make`:

```bash
git clone https://github.com/<username>/cuVCF.git
cd cuVCF

# Build everything (CPU/GPU executables + Python bindings)
make all
````

Or build specific targets:

```bash
# CPU-only executable (columnar version)
make VARCOL

# GPU executable
make GPU

# Python modules (GPU + CPU)
make PYBIND   # -> bin/GPUParser.so
make CPUBIND  # -> bin/CPUParser.so

# Debug builds
make VARCOL_DEBUG
make DEBUG
make CPUBIND_DEBUG
```

Artifacts are produced in `bin/`:

* `bin/VCFparser` — CLI executable (CPU or GPU depending on the target used last)
* `bin/GPUParser.so` — pybind11 module (GPU backend)
* `bin/CPUParser.so` — pybind11 module (CPU backend)
* `bin/CPUParser_debug<pyext>` — debug pybind11 module with ASan/UBSan

> **Note:** The GPU compilation uses `-arch=sm_89` (Ada). Adjust `GPUFLAGS` if your GPU has a different compute architecture

### 2) Minimal CLI usage

The CLI executable is mainly intended for **debugging**.  
It parses a VCF file and converts it into the internal **columnar representation**, but it does **not** provide a complete analysis workflow.  
If you want to use the CLI instead of the Python interface, you’ll need to implement your own `main` function to process the parser’s output.

```bash
# After building with `make VARCOL` or `make GPU`:
./bin/VCFparser <path/to/input.vcf>
````

For full data access and analysis, the recommended entry point is the **Python bindings** (`CPUParser.so` / `GPUParser.so`).

---


### 3) Minimal Python usage

The project provides two Python extension modules built with pybind11:
- `GPUParser` → GPU backend (CUDA required)
- `CPUParser` → CPU backend

Add the `bin/` folder to `PYTHONPATH` (or to `sys.path`) and import one backend:

```python
import sys
sys.path.append("bin")

# Prefer GPU if available
try:
    import GPUParser as cuvcf
except ImportError:
    import CPUParser as cuvcf
````

---

#### Parse a VCF and build DataFrames

The typical workflow is:

1. Create a parser object
2. Run parsing on a VCF
3. Extract columnar tables (DF1–DF4) as Python dicts
4. Wrap them into pandas/cuDF DataFrames

```python
# 1) Create the parser
res = cuvcf.vcf_parsed()

# 2) Parse the VCF (second arg = CPU threads, adjust if needed)
res.run("data/example.vcf", 16)

# 3) Extract DF1–DF4 as dicts of numpy arrays
df1_dict = cuvcf.get_var_columns_data(res.var_columns)   # DF1: variant-level fields
df2_dict = cuvcf.get_alt_columns_data(res.alt_columns)   # DF2: allele-level annotations
df3_dict = cuvcf.get_sample_columns_data(res.samp_columns) # DF3: per-sample genotypes
df4_dict = cuvcf.get_alt_format_data(res.alt_sample)     # DF4: per-sample, per-allele metrics

# 4) Convert to DataFrames
import pandas as pd
DF1 = pd.DataFrame(df1_dict)
DF2 = pd.DataFrame(df2_dict)
DF3 = pd.DataFrame(df3_dict)
DF4 = pd.DataFrame(df4_dict)

# Or with GPU-accelerated cuDF
# import cudf
# DF1 = cudf.DataFrame(df1_dict)
# DF2 = cudf.DataFrame(df2_dict)
# DF3 = cudf.DataFrame(df3_dict)
# DF4 = cudf.DataFrame(df4_dict)
```

---

#### Example queries

```python
# Filter by flag (EVA_4 == 1)
eva = DF1[DF1["EVA_4"] == 1]

# Filter by category (TSA == SNV, encoded as 0)
snv = DF1[DF1["TSA"] == 0]

# Filter by position (POS > 200,000)
high_pos = DF1[DF1["pos"] > 200_000]
```

---

#### Export to CSV (optional)

```python
DF1.to_csv("df1.csv", index=False)
DF2.to_csv("df2.csv", index=False)
DF3.to_csv("df3.csv", index=False)
DF4.to_csv("df4.csv", index=False)
```

---

> **Notes**
>
> * DF1–DF4 are the normalized columnar representations of the VCF:
>
>   * **DF1** → variant-level fields
>   * **DF2** → allele-level annotations
>   * **DF3** → per-sample genotypes
>   * **DF4** → per-sample, per-allele metrics
> * Both `CPUParser` and `GPUParser` expose the same API.
> * Use pandas (CPU) or cuDF (GPU) depending on your environment.


## Repository Structure

```
cuVCF/
│── src/
│   ├── CPUVersion/ # CPU implementation and binding files
│   └── GPUVersion/ # GPU implementation and binding files
│   
│── bcftoolsTest/   # Test scripts for bcftool
│── cyvcf2Tests/    # Test scripts for cyvcf2
│── vcflibTest/     # Test scripts for vcflib
│── TestGPU/        # Test scripts for GPU implementation of cuVCF
│── TestCPU/        # Test scripts for CPU implementation of cuVCF
│── Makefile
│── README.md
│── Tester.bash             # Test script to test the parsing time of cuVCF from the CLI
│── TesterSanitizer.bash    # Test script to test the memory leaks of cuVCF from the CLI
│── README.md
```

---

## Notes & Tips

* Ensure `pybind11` and `python3-dev` (or equivalent) are installed so `python3-config` works.
* If your GPU is not **sm\_89**, update:

  * `GPUFLAGS` / `DEBUG_GPUFLAGS` in the Makefile (e.g., `-arch=sm_80` for A100).
* If importing from Python fails, check:

  * `sys.path` includes `bin/`
  * the extension name matches the module you import (`CPUParser` / `GPUParser`)
  * compatible Python version/ABI (see `PYEXT` in the Makefile)

* If you plan to run the **benchmarking/comparison scripts**, please check their dedicated documentation for the required Python libraries and external tools.  
* The datasets we used in our evaluation are publicly available at the following links:

  * [IRBT (Bos taurus population, EVA)](https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB6119/IRBT.population_sites.UMD3_1.20140322_EVA_ss_IDs.vcf.gz)  
  * [Felis catus (Ensembl release 112)](https://ftp.ensembl.org/pub/release-112/variation/vcf/felis_catus/felis_catus.vcf.gz)  
  * [Bos taurus (Ensembl release 113)](https://ftp.ensembl.org/pub/release-113/variation/vcf/bos_taurus/bos_taurus.vcf.gz)  
  * [Danio rerio (Ensembl release 113)](https://ftp.ensembl.org/pub/release-113/variation/vcf/danio_rerio/danio_rerio.vcf.gz)  
