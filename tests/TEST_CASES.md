# Complete test case index

This index lists the cases on `test/conversion-regressions`. Each conversion case is shared by the CPU and GPU suites; presence here does not mean native execution has passed.

## Conversion cases

Each row below lists the exact executable case IDs. `parse` checks VCF → DF, `write` checks independently constructed DF → VCF, and `roundtrip` checks VCF → DF → VCF. Unsupported-field isolation is parse-only.

| Parse case | Writer case | Round-trip case |
|---|---|---|
| `parse/core/multiallelic` | `write/core/multiallelic` | `roundtrip/core/multiallelic` |
| `parse/info/Integer/1` | `write/info/Integer/1` | `roundtrip/info/Integer/1` |
| `parse/info/Integer/A` | `write/info/Integer/A` | `roundtrip/info/Integer/A` |
| `parse/format/Integer/1/alone` | `write/format/Integer/1/alone` | `roundtrip/format/Integer/1/alone` |
| `parse/format/Integer/1/with_GT` | `write/format/Integer/1/with_GT` | `roundtrip/format/Integer/1/with_GT` |
| `parse/format/Integer/3/alone` | `write/format/Integer/3/alone` | `roundtrip/format/Integer/3/alone` |
| `parse/format/Integer/3/with_GT` | `write/format/Integer/3/with_GT` | `roundtrip/format/Integer/3/with_GT` |
| `parse/format/Integer/A/alone` | `write/format/Integer/A/alone` | `roundtrip/format/Integer/A/alone` |
| `parse/format/Integer/A/with_GT` | `write/format/Integer/A/with_GT` | `roundtrip/format/Integer/A/with_GT` |
| `parse/info/Float/1` | `write/info/Float/1` | `roundtrip/info/Float/1` |
| `parse/info/Float/A` | `write/info/Float/A` | `roundtrip/info/Float/A` |
| `parse/format/Float/1/alone` | `write/format/Float/1/alone` | `roundtrip/format/Float/1/alone` |
| `parse/format/Float/1/with_GT` | `write/format/Float/1/with_GT` | `roundtrip/format/Float/1/with_GT` |
| `parse/format/Float/3/alone` | `write/format/Float/3/alone` | `roundtrip/format/Float/3/alone` |
| `parse/format/Float/3/with_GT` | `write/format/Float/3/with_GT` | `roundtrip/format/Float/3/with_GT` |
| `parse/format/Float/A/alone` | `write/format/Float/A/alone` | `roundtrip/format/Float/A/alone` |
| `parse/format/Float/A/with_GT` | `write/format/Float/A/with_GT` | `roundtrip/format/Float/A/with_GT` |
| `parse/info/String/1` | `write/info/String/1` | `roundtrip/info/String/1` |
| `parse/info/String/A` | `write/info/String/A` | `roundtrip/info/String/A` |
| `parse/format/String/1/alone` | `write/format/String/1/alone` | `roundtrip/format/String/1/alone` |
| `parse/format/String/1/with_GT` | `write/format/String/1/with_GT` | `roundtrip/format/String/1/with_GT` |
| `parse/format/String/3/alone` | `write/format/String/3/alone` | `roundtrip/format/String/3/alone` |
| `parse/format/String/3/with_GT` | `write/format/String/3/with_GT` | `roundtrip/format/String/3/with_GT` |
| `parse/format/String/A/alone` | `write/format/String/A/alone` | `roundtrip/format/String/A/alone` |
| `parse/format/String/A/with_GT` | `write/format/String/A/with_GT` | `roundtrip/format/String/A/with_GT` |
| `parse/info/flag_present` | `write/info/flag_present` | `roundtrip/info/flag_present` |
| `parse/info/flag_absent` | `write/info/flag_absent` | `roundtrip/info/flag_absent` |
| `parse/samples/order_and_alleles` | `write/samples/order_and_alleles` | `roundtrip/samples/order_and_alleles` |
| `parse/missing/GT_unphased` | `write/missing/GT_unphased` | `roundtrip/missing/GT_unphased` |
| `parse/missing/GT_phased` | `write/missing/GT_phased` | `roundtrip/missing/GT_phased` |
| `parse/missing/array_Integer` | `write/missing/array_Integer` | `roundtrip/missing/array_Integer` |
| `parse/missing/array_Float` | `write/missing/array_Float` | `roundtrip/missing/array_Float` |
| `parse/missing/string` | `write/missing/string` | `roundtrip/missing/string` |
| `parse/identifiers/DP2` | `write/identifiers/DP2` | `roundtrip/identifiers/DP2` |
| `parse/identifiers/X1` | `write/identifiers/X1` | `roundtrip/identifiers/X1` |
| `parse/identifiers/AF_and_AF2` | `write/identifiers/AF_and_AF2` | `roundtrip/identifiers/AF_and_AF2` |
| `parse/identifiers/array_AF2` | `write/identifiers/array_AF2` | `roundtrip/identifiers/array_AF2` |
| `parse/special/TSA_SNV` | `write/special/TSA_SNV` | `roundtrip/special/TSA_SNV` |
| `parse/special/TSA_insertion` | `write/special/TSA_insertion` | `roundtrip/special/TSA_insertion` |
| `parse/special/TSA_deletion` | `write/special/TSA_deletion` | `roundtrip/special/TSA_deletion` |
| `parse/special/TSA_inversion` | `write/special/TSA_inversion` | `roundtrip/special/TSA_inversion` |
| `parse/unsupported/isolation_R` | — | — |
| `parse/unsupported/isolation_G` | — | — |
| `parse/unsupported/isolation_variable` | — | — |
| `parse/unsupported/info_fixed` | — | — |
| `parse/core/multiple_records_missing_QUAL` | `write/core/multiple_records_missing_QUAL` | `roundtrip/core/multiple_records_missing_QUAL` |
| `parse/missing/Integer_scalar` | `write/missing/Integer_scalar` | `roundtrip/missing/Integer_scalar` |
| `parse/missing/Integer_A` | `write/missing/Integer_A` | `roundtrip/missing/Integer_A` |
| `parse/missing/Float_scalar` | `write/missing/Float_scalar` | `roundtrip/missing/Float_scalar` |
| `parse/missing/Float_A` | `write/missing/Float_A` | `roundtrip/missing/Float_A` |
| `parse/float/half_precision` | `write/float/half_precision` | `roundtrip/float/half_precision` |
| `parse/identifiers/prefix_overlap` | `write/identifiers/prefix_overlap` | `roundtrip/identifiers/prefix_overlap` |
| `parse/samples/changing_FORMAT` | `write/samples/changing_FORMAT` | `roundtrip/samples/changing_FORMAT` |
| `parse/core/symbolic_and_spanning_ALT` | `write/core/symbolic_and_spanning_ALT` | `roundtrip/core/symbolic_and_spanning_ALT` |
| `parse/identifiers/INFO_FORMAT_same_ID` | `write/identifiers/INFO_FORMAT_same_ID` | `roundtrip/identifiers/INFO_FORMAT_same_ID` |
| `parse/header/reordered_attributes` | `write/header/reordered_attributes` | `roundtrip/header/reordered_attributes` |

## Additional default cases

- `parse/input/no_final_newline`
- `parse/input/gzip`
- `parse/input/header_only`
- `parse/input/missing_file`
- `write/missing/empty_string_cell`
- `write/edited_consistent_values`
- `write/missing_output_directory`
- `write/header_only`

## Extended cases

These require `--extended` and are excluded from the default fast suite.

- `parse/extended/two_threads`
- `write/extended/two_threads`
- `roundtrip/extended/chunk_boundary`

## Test-infrastructure self-checks

- `harness/outcomes`
- `harness/exit_status`
- `harness/cleanup`
- `oracle/format_order_and_float`
- `oracle/duplicate_info`
- `oracle/leading_colon`
- `oracle/empty_value`
- `oracle/renamed_field`
- `oracle/truncated_array`
- `oracle/missing_position`
- `oracle/changed_float`
- `oracle/missing_trailing_field`

## Authoritative runtime list

The C++ registrations are the source of truth. After building on Linux:

```sh
build-tests/tests/cuvcf_test_cpu --list --extended
build-tests/tests/cuvcf_test_harness --list
build-tests/tests/cuvcf_test_oracle --list
```

GPU uses the same case IDs: `build-gpu/tests/cuvcf_test_gpu --list --extended`. Listing does not require a GPU device, but the executable must have been built with CUDA.

Current inventory: **168 default conversion cases**, **3 extended conversion cases**, and **12 infrastructure self-checks**. Update this index when adding or renaming cases.
