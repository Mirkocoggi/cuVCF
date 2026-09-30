# cuVCF

**cuVCF** is a **CPU/GPU-accelerated framework** that converts Variant Call Format (**VCF**) files into **columnar DataFrames**, directly usable with **pandas** and **cuDF**.  
Its **hybrid CPU+GPU pipeline** delivers substantial performance speedups over traditional tools.

## Requirements

### Build toolchain

* **CMake ≥ 3.24** (installed automatically by `pip install .` if missing; on Ubuntu 22.04 apt ships 3.22, so plain CMake builds need a newer one, e.g. `pip install cmake`)
* **g++** with **OpenMP**, **C++17**
* **zlib** and **Imath** development packages (e.g. `zlib1g-dev libimath-dev`)
* **Python 3.8+** development headers (e.g. `python3-dev`); **pybind11** is fetched by pip, or for plain CMake builds it can come from the distro (`python3-pybind11` / `pybind11-dev`) or `pip install pybind11` into the same Python
* **CUDA Toolkit** ≥ 11.8 (tested with 12.x and 13.x) with **nvcc** — optional: without it only the CPU backend is built

### Python (runtime)

* **numpy**
* **pandas** (for CPU analytics)
* **cuDF** (optional, for GPU analytics - same syntax as pandas)

### Tested environment

* **Ubuntu 22.04.5 LTS**, **NVIDIA Ada (sm\_89)** class GPU (e.g., L40S)
* **Ubuntu 24.04.4 LTS**, **NVIDIA L4** (sm\_89), CUDA 13.1, CMake 3.28, GCC 13, Python 3.12

Other GPUs work too: by default the build targets the GPU found on the build machine (see *Notes & Tips*).

---

## Quick Start

### 1) Install the prerequisites

On Ubuntu:

```bash
sudo apt install build-essential cmake zlib1g-dev libimath-dev python3-dev python3-venv python3-pybind11
```

For the GPU backend, also install the [CUDA Toolkit](https://developer.nvidia.com/cuda-downloads) and make sure `nvcc` is on your `PATH`:

```bash
export PATH=/usr/local/cuda/bin:$PATH   # if nvcc is not found
nvcc --version
```

Without CUDA everything still builds, but only the CPU backend.

### 2) Build

#### Option A — Python modules with pip (recommended)

```bash
git clone https://github.com/Mirkocoggi/cuVCF.git
cd cuVCF
python3 -m venv .venv
source .venv/bin/activate
pip install .
```

This builds and installs `CPUParser` and, if CUDA is available, `GPUParser` into the virtual environment. The first build compiles the CUDA code and takes a few minutes.

* The virtual environment is required on recent Ubuntu/Debian releases, which refuse `pip install` into the system Python (`externally-managed-environment`).
* CPU-only build: `pip install . -C cmake.define.CUVCF_CUDA=OFF`
* Build for a specific GPU architecture: `pip install . -C cmake.define.CMAKE_CUDA_ARCHITECTURES=80`

Check the installation:

```bash
python -c "import CPUParser; print('CPU backend OK')"
python -c "import GPUParser; print('GPU backend OK')"   # CUDA builds only
```

#### Option B — CLIs and Python modules with CMake

```bash
git clone https://github.com/Mirkocoggi/cuVCF.git
cd cuVCF
cmake -S . -B build
cmake --build build -j
```

Artifacts in `build/`:

* `VCFparser_cpu` — CLI executable, CPU backend
* `VCFparser_gpu` — CLI executable, GPU backend (CUDA only)
* `CPUParser.cpython-*.so` — pybind11 module (CPU backend)
* `GPUParser.cpython-*.so` — pybind11 module (GPU backend, CUDA only)

Build options (`-D<option>=<value>` when configuring):

| Option | Default | Meaning |
|---|---|---|
| `CUVCF_CUDA` | `AUTO` | `AUTO`: build the GPU backend if nvcc is found; `ON`: require it; `OFF`: CPU only |
| `CMAKE_CUDA_ARCHITECTURES` | `native` (or `75;80;86;89;90` if no GPU is visible) | GPU architectures to compile for, e.g. `89` |
| `CMAKE_BUILD_TYPE` | `Release` | `Debug` adds `-g` (and `-G -lineinfo` for CUDA) and drops `-O3` (i.e. `-O0`) |
| `CUVCF_SANITIZE` | `OFF` | build the CPU targets with AddressSanitizer + UBSan |

Debug build example: `cmake -S . -B build-dbg -DCMAKE_BUILD_TYPE=Debug -DCUVCF_SANITIZE=ON && cmake --build build-dbg -j`

With `CUVCF_SANITIZE=ON` the `CPUParser` module is built with ASan, so the ASan runtime must be preloaded to import it:

```bash
LD_PRELOAD=$(gcc -print-file-name=libasan.so) ASAN_OPTIONS=detect_leaks=0 PYTHONPATH=build-dbg python3 -c "import CPUParser"
```

### 3) Run from the command line

The CLI executables are mainly intended for **debugging**.  
They parse a VCF file and convert it into the internal **columnar representation**, but they do **not** provide a complete analysis workflow.  
If you want to use the CLI instead of the Python interface, you’ll need to implement your own `main` function to process the parser’s output.

After an Option B build:

```bash
./build/VCFparser_gpu -v data/tiny.vcf -t 4
```

* `-v <file>` — input VCF, required (`.vcf`, or `.vcf.gz`: a gzipped input is decompressed **in place**, replacing the `.gz` file)
* `-t <n>` — number of CPU threads (default 4)
* `VCFparser_cpu` takes the same options.

The CLI prints the parse time on stderr and the first 10 rows of the four DataFrames on stdout.

For full data access and analysis, the recommended entry point is the **Python bindings** (`CPUParser` / `GPUParser`).

---

### 4) Use from Python

The project provides two Python extension modules built with pybind11:
- `GPUParser` → GPU backend (CUDA required)
- `CPUParser` → CPU backend

After Option A, activate the virtual environment and run your script as usual:

```bash
source .venv/bin/activate
python my_script.py
```

After Option B, point Python at the build directory instead:

```bash
PYTHONPATH=build python3 my_script.py
```

In your script, import one backend:

```python
# Prefer GPU if available
try:
    import GPUParser as cuvcf
except ImportError:
    import CPUParser as cuvcf
```

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
│── vcflibTests/    # Test scripts for vcflib
│── TestGPU/        # Test scripts for GPU implementation of cuVCF
│── TestCPU/        # Test scripts for CPU implementation of cuVCF
│── CMakeLists.txt  # Build definition (CLIs + Python modules)
│── pyproject.toml  # pip install . (scikit-build-core)
│── README.md
│── Tester.bash             # Test script to test the parsing time of cuVCF from the CLI
│── TesterSanitizer.bash    # Test script to test the memory leaks of cuVCF from the CLI
│── README.md
```

---

## Notes & Tips

* To build for a GPU other than the one on the build machine, pass `-DCMAKE_CUDA_ARCHITECTURES=<cc>` (e.g. `80` for A100), or `pip install . -C cmake.define.CMAKE_CUDA_ARCHITECTURES=80`. With the default `native` setting the build (including a pip-built wheel) only runs on the build machine's GPU generation.
* If importing from Python fails, check:

  * the package was installed with the same interpreter you run (`python -m pip install .`), or `build/` is on `PYTHONPATH` for plain CMake builds
  * `GPUParser` is missing: CUDA was not found at build time (the CMake log says `GPU backend skipped`)

* If you plan to run the **benchmarking/comparison scripts**, please check their dedicated documentation for the required Python libraries and external tools.  
* The datasets we used in our evaluation are publicly available at the following links:

  * [IRBT (Bos taurus population, EVA)](https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB6119/IRBT.population_sites.UMD3_1.20140322_EVA_ss_IDs.vcf.gz)  
  * [Felis catus (Ensembl release 112)](https://ftp.ensembl.org/pub/release-112/variation/vcf/felis_catus/felis_catus.vcf.gz)  
  * [Bos taurus (Ensembl release 113)](https://ftp.ensembl.org/pub/release-113/variation/vcf/bos_taurus/bos_taurus.vcf.gz)  
  * [Danio rerio (Ensembl release 113)](https://ftp.ensembl.org/pub/release-113/variation/vcf/danio_rerio/danio_rerio.vcf.gz)  
