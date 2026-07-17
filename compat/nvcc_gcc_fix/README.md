# nvcc / gcc intrinsics workaround

`nvcc` 12.4's internal front-end (`cicc`/`cudafe`) does not recognize a handful of
very new x86 builtins (`__builtin_ia32_ldtilecfg`, `__builtin_ia32_cvtne2ps2bf16_*`, ...)
that recent gcc packages declare unconditionally in `<amxtileintrin.h>`,
`<avx512bf16intrin.h>` and `<avx512bf16vlintrin.h>` (pulled in transitively via
`<immintrin.h>`). Real gcc compiles these headers fine; only nvcc's parser chokes
on them.

This project never uses AMX-tile or AVX512-BF16 intrinsics, so the three headers
below are stub replacements (same include guard, empty body) that shadow the real
system headers for the nvcc build only. They are picked up because this directory
is added to the include path ahead of the system include dirs (see `-I` in the
Makefile's `NVCCFLAGS`).

If a future CUDA Toolkit release fixes this in nvcc's front-end, this directory
(and the `-I` flag referencing it) can be removed.
