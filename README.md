# VCF Reconstructor

Reconstructs a standard VCF file from the column-oriented DataFrames produced by **cuVCF**. Two backends are provided: a CPU backend parallelized with OpenMP, and a GPU backend built on CUDA.

---

## Overview

cuVCF compresses VCF data by splitting it into four typed, column-oriented DataFrames stored as CSV files. VCF Reconstructor reads those files back and reassembles them into a valid VCF.

| DataFrame | Contents |
|-----------|----------|
| **DF1** | Core variant fields: CHROM, POS, ID, REF, QUAL, FILTER, scalar INFO |
| **DF2** | Alternate alleles and per-allele INFO fields (one row per ALT) |
| **DF3** | Per-sample FORMAT fields constant across alleles (GT, DP, GQ, …) |
| **DF4** | Per-sample, per-allele FORMAT fields (AD, …) |

Both backends share the same parsing layer (`CSVParser`) and data structures (`VCFDataFrames.h`). The GPU backend further pre-formats all sample strings on the CPU using OpenMP before handing a single packed buffer to the CUDA kernel.

---

## Architecture

```
CSVParser ──► var_columns_df  (DF1) ──┐
             alt_columns_df  (DF2) ──┤
             sample_columns_df(DF3) ──┤──► VCFReconstructor     ──► output.vcf
             alt_format_df   (DF4) ──┘    (CPU: OpenMP chunks)

                                     or

                                          VCFReconstructorGPU   ──► output.vcf
                                          (GPU pipeline below)
```

### GPU Pipeline (per chunk)

```
allocateDevice  →  prepareHostBuffers  →  uploadToDevice (H2D)
     →  reconstructKernel  →  compactKernel (CUB prefix-sum)
     →  D2H copy  →  async disk write (background thread)
     →  freeDevice
```

- **`reconstructKernel`** — one CUDA thread per variant; grid-stride loop for large chunks. Writes each VCF line into a sparse 2D buffer (`chunk_size × MAX_LINE_LEN`).
- **`compactKernel`** — uses a CUB exclusive prefix-sum over line lengths to pack the sparse buffer into a contiguous byte stream.
- **Double-buffered writer thread** — disk I/O runs concurrently with the next chunk's GPU work, hiding write latency.

### CPU Pipeline

Variants are split into equal chunks, one per OpenMP thread. Each thread writes to a private temporary file; the main thread concatenates them in order into the final output.

---

## Dependencies

| Dependency | Purpose |
|------------|---------|
| C++17 | Required by both backends |
| OpenMP | CPU parallel reconstruction; also used inside the GPU backend for `buildSampleStrings` |
| CUDA Toolkit ≥ 11 | GPU backend (`nvcc`, CUB, NVTX) |
| [Imath](https://github.com/AcademySoftwareFoundation/Imath) | `Imath/half.h` — half-precision floats for QUAL and float INFO/FORMAT fields |

CUB and NVTX are bundled with the CUDA Toolkit and require no separate installation.

---

## Building

No build system file is included yet; compile directly with the commands below.

### CPU backend

```bash
g++ -std=c++17 -O2 -fopenmp \
    -I include \
    src/Utils.cpp \
    src/CPUVersion/CSVParser.cpp \
    src/CPUVersion/VCFReconstructor.cpp \
    src/CPUVersion/main_cpu.cpp \
    -o vcf_reconstructor_cpu \
    -lImath
```

### GPU backend

```bash
nvcc -std=c++17 -O2 -Xcompiler -fopenmp \
    -I include \
    src/Utils.cpp \
    src/CPUVersion/CSVParser.cpp \
    src/GPUVersion/VCFReconstructorGPU.cu \
    src/GPUVersion/main_gpu.cu \
    -o vcf_reconstructor_gpu \
    -lImath -lnvToolsExt
```

Output binaries are conventionally placed under `build/`.

---

## Input Files

cuVCF produces the following files for each dataset:

```
data/<dataset>/
├── df1.csv           # Core variant data
├── df2.csv           # Alternate alleles
├── df3.csv           # Core sample FORMAT data
├── df4.csv           # Per-allele sample FORMAT data
├── <name>_header.txt # Raw VCF header text
└── maps_used_<name>.csv  # Encoding dictionaries (CHROM, FILTER, …)
```

The maps file encodes compact char/int codes used internally to reduce memory footprint. It must be loaded before the DataFrames via `parseMaps`.

---

## Usage

Both `main_cpu.cpp` and `main_gpu.cu` demonstrate the full call sequence. In production code the `CSVParser` section is replaced by the DataFrames already held in memory by cuVCF.

```cpp
// 1. Parse encoding dictionaries and all four DataFrames
CSVParser parser("data/mydata/df1.csv",
                 "data/mydata/df2.csv",
                 "data/mydata/df3.csv",
                 "data/mydata/df4.csv",
                 "data/mydata/header.txt");

var_columns_df    df1;
alt_columns_df    df2;
sample_columns_df df3;
alt_format_df     df4;

parser.parseMaps("data/mydata/maps_used.csv", df1, df3, df4);
parser.loadAll(df1, df2, df3, df4);

// 2. Reconstruct — choose one backend
VCFReconstructor reconstructor("output.vcf", parser.header_text);   // CPU
// VCFReconstructorGPU reconstructor("output.vcf", parser.header_text); // GPU

reconstructor.run(df1, df2, df3, df4);
```

---

## Configuration (GPU backend)

The following compile-time constants in `VCFReconstructorGPU.h` control performance:

| Constant | Default | Description |
|----------|---------|-------------|
| `CHUNK_SIZE` | 100 000 | Variants processed per GPU launch |
| `BLOCK_SIZE` | 64 | CUDA threads per block |
| `NUM_BLOCKS` | 128 | CUDA blocks per grid |
| `MAX_LINE_LEN` | 756 | Maximum bytes per VCF line in the sparse buffer |
| `MAX_SAMPLE_STRING_LEN` | 64 | Maximum bytes per sample FORMAT string |
| `MAX_NAME_LEN` | 16 | Maximum characters for a field name on the device |

Increase `MAX_LINE_LEN` if your dataset has very long INFO fields or many samples.

---

## Data Encoding

cuVCF applies several space-saving encodings that the reconstructor reverses:

- **CHROM / FILTER** — mapped to single `char` codes; reversed via `inv_chrom_map` / `inv_filter_map`.
- **QUAL and float INFO/FORMAT fields** — stored as `half` (fp16); `–1.0f` is the missing-value sentinel, rendered as `.` in the output.
- **Integer INFO/FORMAT fields** — stored as `int`; `–1` is the missing-value sentinel.
- **Genotype (GT)** — encoded as a single `char` via `GTMap`; decoded with `getGTStringFromChar`.
- **TSA, PolyPhen, CSQ** — mapped to compact integer/char codes (see `include/Maps.h`); reversed on the GPU using the `DeviceMaps` dictionary buffers.

---

## Project Structure

```
VCF_Reconstructor/
├── include/
│   ├── VCFDataFrames.h      # All DataFrame struct/class definitions
│   ├── Utils.h              # String utility functions (splitLine, stripTypePrefix, …)
│   └── Maps.h               # Static encoding maps (TSA, PolyPhen, CSQ)
├── src/
│   ├── Utils.cpp
│   ├── CPUVersion/
│   │   ├── CSVParser.h/.cpp       # Parses cuVCF CSV exports into DataFrames
│   │   ├── VCFReconstructor.h/.cpp# OpenMP-based reconstruction
│   │   └── main_cpu.cpp
│   └── GPUVersion/
│       ├── GPUStructs.h           # Device-side shadow structures (DeviceVarColumns, …)
│       ├── CudaUtils.cuh          # Device utility functions (itoa, ftoa, strcpy, …)
│       ├── VCFReconstructorGPU.h/.cu  # CUDA reconstruction pipeline
│       └── main_gpu.cu
└── data/
    └── <dataset>/            # Input CSV files and headers (not tracked by git)
```
