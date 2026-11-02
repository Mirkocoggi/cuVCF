# Conversion regression tests

Parser stage: `test/parser-regressions`, based on master `6fe5a5f`. Full suite: `test/conversion-regressions`, rebased onto reconstructor `b2b04b5436ea72f819d3ffd4d1cafb5d6b2c5ee2`. Tests target the existing Linux C++17/OpenMP toolchain; no Apple Silicon build adaptations are provided. Native Linux and NVIDIA execution is pending. A test being present does **not** mean it has passed.


## Folder structure

```text
tests/
├── CMakeLists.txt          # Build targets and CTest registration
├── README.md               # Build/run instructions and failure ledger
├── TEST_CASES.md           # Complete expanded case index
├── cases/
│   ├── fixtures.h          # Small VCF inputs and independently specified DF expectations
│   └── reconstruction.h    # Writer and round-trip case registration
├── runners/
│   ├── conversion.cpp      # CPU entry point; parser case registration and suite assembly
│   └── conversion_gpu.cu   # Compiles the shared suite for the GPU backend
├── support/
│   ├── backend.h           # Backend adapter, DF builders and dataframe assertions
│   ├── harness.h           # Assertions, selection, reporting and temporary files
│   └── expected_failures.h # Confirmed failures, selected by backend and exact case
├── oracle/
│   └── vcf_oracle.h        # Independent semantic checker for reconstructed VCF
└── selftest/
    ├── harness_selftest.cpp
    └── oracle_selftest.cpp
```

Start with [the complete case list](TEST_CASES.md). The same synthetic fixture normally produces three cases: `parse/` compares parser output with expected frames; `write/` passes those independently built frames directly to the writer; `roundtrip/` checks the composed conversion. The independent checker validates writer output without using the production parser. Each backend runs its selected cases in one process.

Fixtures live in C++ so the literal VCF input and expected dataframe values stay together. To add coverage, add a fixture in `cases/fixtures.h`; it automatically receives all three directions unless marked parse-only. Add special writer cases in `cases/reconstruction.h`, and update `TEST_CASES.md`. The runtime `--list --extended` command is authoritative.

## Build once, run repeatedly

Install the normal project prerequisites (including pybind11). From the repository root:

```sh
cmake -S . -B build-tests -DCUVCF_CUDA=OFF -DCUVCF_BUILD_TESTS=ON
cmake --build build-tests --target cuvcf_test_cpu cuvcf_test_harness cuvcf_test_oracle -j2
ctest --test-dir build-tests -L 'cpu|harness' --output-on-failure --no-tests=error
```

Use `build-tests/tests/cuvcf_test_cpu --list` to list cases, `--filter substring` to select a case family, and `--extended` for extra threading and 10,001-record chunk-boundary checks. Empty selections and unknown arguments return nonzero. Rebuild incrementally after source changes; tests do not rebuild on every invocation. No pre-commit hook is installed automatically.

For an NVIDIA Linux machine, configure a separate build with `-DCUVCF_CUDA=ON`, build `cuvcf_test_gpu`, and run `ctest --test-dir build-gpu -L gpu --output-on-failure --no-tests=error`. A requested GPU run fails if the device is unavailable. Configure CUDA OFF for a CPU-only run; avoid bare CTest in a CUDA-enabled build when no device is available.

All ordinary cases run in one process per backend with one parser thread. There is no per-case shell command, Python interpreter, pandas, cuDF, dataset download or parallel test-worker pool. Temporary synthetic inputs are removed after each case, including failed assertions. The gzip test uses a disposable copy because parsing consumes the compressed file. CTest caps CPU/GPU suite time at 60/120 seconds. Runtime and memory measurements remain pending; the target is a few seconds for the default CPU suite after building.

## Supported subset and case map

| Feature | Case prefix | Contract |
|---|---|---|
| Core columns and alternate alleles | `core/` | Values, ordering, IDs and lookup maps; absent QUAL uses the documented -1 sentinel |
| INFO Integer/Float/String | `info/` | Number=1 and A |
| INFO flags | `info/flag_` | Number=0, presence and absence |
| FORMAT Integer/Float/String | `format/` | Number=1, fixed positive count, A; alone and after GT |
| Diploid GT and samples | `samples/`, `missing/GT_` | Phased/unphased, missing codes, sample order, allele relationships, changing per-record FORMAT |
| Missing values | `missing/` | Preserve positions and distinguish missing from zero and valid negative numbers |
| Original field IDs | `identifiers/` | Scalar digit suffixes, distinct X1/X2, AF/AF2, digit-ending array IDs |
| TSA special encoding | `special/` | SNV/insertion/deletion/inversion codes |
| Unsupported-field isolation | `unsupported/` | Ignored fields cannot shift or overwrite supported neighbours |
| Input handling | `input/` | Header-only, absent file, gzip, final newline |

INFO fixed arrays (>1), Number=R/G/., Character fields, non-diploid/partially missing GT, and nonstandard FORMAT Flag/GT Number=A are outside the supported lossless conversion contract. The implementation contains some nonstandard paths; their existence does not make them valid VCF features. Unsupported-field isolation tests are parse-only and do not assert lossless reconstruction.

Expectations are hand-authored native DF1–DF4 values, independent of the parser. Each test compares every populated data column and row association; missing integers use INT_MIN and missing floats NaN as the intended lossless contract. The old master encodes some missing values differently; tests intentionally expose that limitation. Finite floats are compared after half quantization. QUAL separately uses -1 for missing, consistent with the existing nonnegative QUAL contract.

## Failure ledger and triage

No native results from this revision are available yet. Do not pre-emptively mark suspected bugs as expected failures. Run on Linux, record the precise case and assertion, and then add an expected failure **only for that backend and direction**. Keep crash/setup failures active. The harness accepts an exact assertion key and treats an unexpected pass as failure; unrelated assertions and exceptions still fail.

| Candidate cases | Expected behavior / suspected issue | Status |
|---|---|---|
| `parse/missing/array_Integer`, `parse/missing/array_Float` | Preserve missing elements; master numeric fallback converts them to zero | Source-inspected; native reproduction pending |


Add confirmed failures in `support/expected_failures.h` and this ledger; the registry is currently empty. Each entry selects an exact backend, case ID and assertion key, and must include a reason. Stale case IDs and duplicate entries fail setup.

For each confirmed failure add: backend, base commit, case ID, exact assertion key, command, expected/actual result and issue/PR reference. Fix independently; remove the entry when its strict expected failure becomes an unexpected pass. Never weaken the oracle to make CI green.

## Harness checks

`cuvcf_test_harness` checks pass/fail classification, narrow expected failures, unexpected passes, setup-error propagation and temporary-file cleanup. It has no native parser dependency and can be compiled with a standard C++17 compiler. Suite exit codes are 0 for pass/expected-failure, 1 for test failure/unexpected-pass and 2 for invalid invocation or missing GPU.

## Writer and round-trip coverage

Every supported fixture has three separately reported cases: `parse/`, `write/` and `roundtrip/`. Unsupported-field isolation is parse-only. Writer cases construct native frames directly; round trips parse VCF then reconstruct it. Additional writer cases cover a consistent value edit, an empty string cell, header-only output and output-open errors. No filtering/selection algorithm or arbitrary inconsistent user dataframe is tested.

The standalone checker in `oracle/vcf_oracle.h` parses only the synthetic supported subset, not arbitrary VCF. It checks exact header metadata, sample/record order, declared keys, duplicate keys, cardinalities and missing positions. INFO and FORMAT ordering may differ. Trailing missing sample values may be omitted, and a dataset-wide FORMAT union is allowed only if added values are missing. Floating-point comparison allows half-precision quantization and six-decimal serialization rounding. It does not use the cuVCF parser, so matching parser/writer bugs cannot validate each other. Cases use ordinary finite half-range values; overflow/infinity serialization is outside this subset.

### Reconstructor failure candidates — native reproduction still required

These are review leads, **not confirmed test results or active expected failures**. New fixes for missing GT and empty strings are included in the rebased source; those regression tests remain ordinary assertions.

| Cases (CPU unless specified) | Intended result / reviewed risk | Previous review |
|---|---|---|
| `write/format/{Integer,Float}/3/*` | Preserve three values even with only one ALT; current grouping may clamp to ALT+REF | R02 |
| `write/missing/array_{Integer,Float}` | Keep leading missing array slots | R03 |
| `write/identifiers/{DP2,X1,AF_and_AF2,array_AF2}` | Preserve original IDs and distinct fields; also verify GPU grouping | R04 |
| `write/format/String/3/*` | Reassemble fixed-size string arrays under the original ID | R05 |
| `write/format/{Integer,Float,String}/A/alone` | Emit a sample without a leading colon | R06 |
| `parse/samples/changing_FORMAT`, corresponding writer/round-trip cases | Missing GT must not become 0\|0 on records with no GT | R18 |
| `parse/identifiers/INFO_FORMAT_same_ID` | INFO and FORMAT have separate ID namespaces; their types must not overwrite one another | New conversion regression |

The `roundtrip/` counterparts exercise the same risks across the composed flow. README patterns group cases; use `--list` for exact names. Known defects can make the first native run red until each is reproduced and recorded narrowly in the registry.

## Validation performed / pending

- Passed: three harness self-tests; nine independent VCF-checker self-tests.
- Passed: harness CLI selection, listing, invalid-argument and empty-selection exit-code checks.
- Passed: CPU suite C++17 syntax-only check against the repository headers. Existing production-header unused-expression warning remains unchanged.
- Pending: Linux CMake configure/build, actual CPU parsing/reconstruction, sanitizer execution, timing and memory measurement.
- Pending: CUDA compilation and all NVIDIA execution, including chunk boundaries.

No native conversion results are claimed from Apple Silicon. Build artifacts for the portable harness/checker checks are outside the repository. No platform-specific code, automatic hook or CI workflow was added.
