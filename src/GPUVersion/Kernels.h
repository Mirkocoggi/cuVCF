/**
 * @file Kernels.h
 * @brief Host interface of the parsing kernel compiled in Kernels.cu
 *
 * @details The kernel, its device helpers (CUDAUtils.cuh) and the constant-memory tables it reads
 * all live in Kernels.cu: a translation unit that defined the tables itself would get its own copy,
 * which the kernel never sees. Other units fill the tables and launch the kernel through the
 * functions below.
 */

#ifndef KERNELS_H
#define KERNELS_H

#include <cuda_runtime.h>
#include "DataStructures.h"

/// Maximum number of keys in the GT map.
#define NUM_KEYS_GT 244

/// Maximum length for each key in the GT map.
#define MAX_KEY_LENGTH_GT 5

/// Maximum number of tokens when splitting strings.
#define MAX_TOKENS 16

/// Maximum length for each token after splitting.
#define MAX_TOKEN_LEN 32

/// Maximum length for temporary string buffers.
#define MAX_TMP_LEN 128

/// Maximum number of keys in Map1.
#define NUM_KEYS_MAP1 128
/// Maximum length for each key in Map1.
#define MAX_KEY_LENGTH_MAP1 32

/// Copies the GT codes into the kernel's constant memory (d_keys_gt, d_values_gt).
cudaError_t upload_gt_table(const char (&keys)[NUM_KEYS_GT][MAX_KEY_LENGTH_GT], const char (&values)[NUM_KEYS_GT]);

/// Copies the INFO/FORMAT name table into the kernel's constant memory (d_keys_map1, d_values_map1).
cudaError_t upload_map1_table(const char (&keys)[NUM_KEYS_MAP1][MAX_KEY_LENGTH_MAP1], const int (&values)[NUM_KEYS_MAP1]);

/// Launches the parsing kernel on stream; returns the launch error, if any.
cudaError_t launch_parse_kernel(int blocks, int threads, cudaStream_t stream, const KernelParams* params, char* my_mem, int batch_size, bool hasSamp);

#endif
