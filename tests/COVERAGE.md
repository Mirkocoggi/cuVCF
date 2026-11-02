# Coverage audit

Audit of `test/conversion-regressions`, against the parser and reconstructor code inherited from `b2b04b5`. This is a source/fixture audit; native CPU and GPU results are still pending.

## Conclusion

The suite exercises the principal supported type/cardinality paths, but it does **not** cover all supported values or all combinations. Shared CPU/GPU case definitions do not establish backend parity until both execute successfully. The case count includes three directions per fixture and must not be read as a count of independent features.

| Dimension | Existing coverage | Audit additions / limits |
|---|---|---|
| INFO | Integer/Float/String × Number=1/A; Flag Number=0 | Combined all six payload classes and a flag in one record |
| FORMAT | Integer/Float/String × scalar/fixed/A, alone and following GT | Combined scalar, fixed and A fields of all three types |
| DF1–DF4 interaction | Mostly isolated types; integer sample/allele interaction | Three all-dataframe fixtures: no GT, present GT, missing GT |
| Missing arrays | Leading numeric missing values | Leading, middle and trailing slots for Integer/Float/String fixed arrays |
| Names | Digit suffixes, prefix overlap, INFO/FORMAT same ID | Already covered; not every name/type combination |
| Samples/records | Two sample order, multiple records, changing FORMAT | Mixed all-type fixture is one sample/record; fully mixed multi-record/multi-sample case remains absent |
| GT | Phased/unphased and fully missing, with scalar fields | Added missing GT with all three allele FORMAT types; encoding boundary allele 10 and the full diploid code space remain absent |
| Special strings | Four TSA codes; symbolic ALT and spanning deletion | Additional TSA codes/aliases and breakend ALT syntax remain absent |
| Input/output | Gzip, header-only, no final newline, open failures | CRLF, long strings and failed writes after opening are not exercised |
| Boundaries | Optional 10,001-record chunk case and two threads | High chromosome/sample counts, multi-digit fixed cardinality, large offsets and allocation failures remain absent |

## Deliberate limits

Number=R/G/., fixed-size INFO arrays, Character, non-diploid/partially missing GT and nonstandard FORMAT flags are outside the declared lossless subset. Unsupported-field isolation cases assert only that supported neighbours survive parsing. Overflow/nonfinite float serialization and arbitrary inconsistent dataframe mutations are outside this suite's current contract.

Supported does not mean defect-free: tests assert intended semantics and may expose existing bugs. No new expected-failure entries were added without native reproduction. See the README ledger for existing review candidates.

## Reading and extending the suite

The complete expanded IDs are in [TEST_CASES.md](TEST_CASES.md). The 12 fixtures added by this audit create 36 default conversion tests, bringing the inventory to 204 default cases plus three extended cases. The mixed fixtures expose interactions cheaply, but an early assertion failure can mask later defects in that same fixture; isolated cases remain necessary for diagnosis.

Prioritize the remaining gaps above based on intended datasets. Do not claim exhaustive combinations: scalar counts, field sets, sample counts, field ordering and missingness produce an unbounded space. A future pairwise matrix can make interaction coverage measurable without multiplying the fast suite indiscriminately.

Validation: updated CPU suite passes a C++17 syntax-only check. Native Linux execution, CUDA compilation and GPU execution remain pending. No production source changes were made by this audit.
