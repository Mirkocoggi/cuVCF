# Conversion regression tests

Parser stage based on master `6fe5a5f`. Tests target the existing Linux C++17/OpenMP toolchain; no Apple Silicon build adaptations are provided. Native Linux and NVIDIA execution is pending. A test being present does **not** mean it has passed.

## Build once, run repeatedly

Install the normal project prerequisites (including pybind11). From the repository root:

```sh
cmake -S . -B build-tests -DCUVCF_CUDA=OFF -DCUVCF_BUILD_TESTS=ON
cmake --build build-tests --target cuvcf_test_cpu cuvcf_test_harness -j2
ctest --test-dir build-tests -L 'cpu|harness' --output-on-failure --no-tests=error
```

Use `build-tests/tests/cuvcf_test_cpu --list` to list cases, `--filter substring` to select a case family, and `--extended` for the small extra threading checks. Empty selections and unknown arguments return nonzero. Rebuild incrementally after source changes; tests do not rebuild on every invocation. No pre-commit hook is installed automatically.

For an NVIDIA Linux machine, configure a separate build with `-DCUVCF_CUDA=ON`, build `cuvcf_test_gpu`, and run `ctest --test-dir build-gpu -L gpu --output-on-failure --no-tests=error`. A requested GPU run fails if the device is unavailable. Configure CUDA OFF for a CPU-only run; avoid bare CTest in a CUDA-enabled build when no device is available.

All ordinary cases run in one process per backend with one parser thread. There is no per-case shell command, Python interpreter, pandas, cuDF, dataset download or parallel test-worker pool. Temporary synthetic inputs are removed after each case, including failed assertions. The gzip test uses a disposable copy because parsing consumes the compressed file. CTest caps CPU/GPU suite time at 60/120 seconds. Runtime and memory measurements remain pending; the target is a few seconds for the default CPU suite after building.

## Supported subset and case map

| Feature | Case prefix | Contract |
|---|---|---|
| Core columns and alternate alleles | `core/` | Values, ordering, IDs and lookup maps; absent QUAL uses -1 |
| INFO Integer/Float/String | `info/` | Number=1 and A |
| INFO flags | `info/flag_` | Number=0, presence and absence |
| FORMAT Integer/Float/String | `format/` | Number=1, fixed positive count, A; alone and after GT |
| Diploid GT and samples | `samples/`, `missing/GT_` | Phased/unphased, missing codes, sample order, allele relationships |
| Missing values | `missing/` | Preserve positions and distinguish missing from zero and valid negative numbers |
| Original field IDs | `identifiers/` | Scalar digit suffixes, distinct X1/X2, AF/AF2, digit-ending array IDs |
| TSA special encoding | `special/` | SNV/insertion/deletion/inversion codes |
| Unsupported-field isolation | `unsupported/` | Ignored fields cannot shift or overwrite supported neighbours |
| Input handling | `input/` | Header-only, absent file, gzip, final newline |

INFO fixed arrays (>1), Number=R/G/., Character fields, non-diploid/partially missing GT, and nonstandard FORMAT Flag/GT Number=A are outside the supported lossless conversion contract. The implementation contains some nonstandard paths; their existence does not make them valid VCF features. Unsupported-field isolation tests are parse-only and do not assert lossless reconstruction.

Expectations are hand-authored native DF1–DF4 values, independent of the parser. Each test compares every populated data column and row association; missing integers use INT_MIN and missing floats NaN as the intended lossless contract. The old master encodes some missing values differently; tests intentionally expose that limitation. Finite floats are compared after half quantization. QUAL uses -1 for missing, consistent with its nonnegative valid range.

## Failure ledger and triage

No native results from this revision are available yet. Do not pre-emptively mark suspected bugs as expected failures. Run on Linux, record the precise case and assertion, and then add an expected failure **only for that backend and direction**. Keep crash/setup failures active. The harness accepts an exact assertion key and treats an unexpected pass as failure; unrelated assertions and exceptions still fail.

| Candidate cases | Expected behavior / suspected issue | Status |
|---|---|---|
| `parse/missing/array_Integer`, `parse/missing/array_Float` | Preserve missing elements; master numeric fallback converts them to zero | Source-inspected; native reproduction pending |


Record confirmed failures in `expected_failures.h` (currently empty) and this README.

For each confirmed failure add: backend, base commit, case ID, exact assertion key, command, expected/actual result and issue/PR reference. Fix independently; remove the entry when its strict expected failure becomes an unexpected pass. Never weaken the oracle to make CI green.

## Harness checks

`cuvcf_test_harness` checks pass/fail classification, narrow expected failures, unexpected passes, setup-error propagation and temporary-file cleanup. It has no native parser dependency and can be compiled with a standard C++17 compiler. Suite exit codes are 0 for pass/expected-failure, 1 for test failure/unexpected-pass and 2 for invalid invocation or missing GPU.

Harness self-checks and CLI exit-code checks passed. CPU test sources passed a syntax-only check. Native Linux execution and NVIDIA validation remain pending.
