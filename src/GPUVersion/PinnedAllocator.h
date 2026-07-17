#ifndef PINNED_ALLOCATOR_H
#define PINNED_ALLOCATOR_H

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstddef>

/**
 * @brief Minimal std::allocator-compatible wrapper around cudaHostAlloc/cudaFreeHost.
 *
 * Lets a std::vector back its storage with page-locked (pinned) host memory
 * instead of regular pageable memory, while keeping all the usual vector
 * convenience methods (resize/assign/push_back/clear) unchanged. Pinned
 * memory is required for cudaMemcpyAsync to actually run asynchronously —
 * with pageable memory the driver silently falls back to a synchronous copy.
 *
 * Deliberately self-contained (no dependency on CudaUtils.cuh's gpuErrchk)
 * so it can be included from a plain header without risking a duplicate
 * (unguarded) inclusion of CudaUtils.cuh in the same translation unit.
 */
template <typename T>
struct PinnedAllocator {
    using value_type = T;

    PinnedAllocator() noexcept = default;
    template <typename U>
    PinnedAllocator(const PinnedAllocator<U>&) noexcept {}

    T* allocate(std::size_t n) {
        if (n == 0) return nullptr;
        void* p = nullptr;
        cudaError_t err = cudaHostAlloc(&p, n * sizeof(T), cudaHostAllocDefault);
        if (err != cudaSuccess) {
            fprintf(stderr, "PinnedAllocator: cudaHostAlloc failed: %s\n", cudaGetErrorString(err));
            exit(err);
        }
        return static_cast<T*>(p);
    }

    void deallocate(T* p, std::size_t) noexcept {
        if (p) cudaFreeHost(p);
    }
};

template <typename T, typename U>
bool operator==(const PinnedAllocator<T>&, const PinnedAllocator<U>&) noexcept { return true; }
template <typename T, typename U>
bool operator!=(const PinnedAllocator<T>&, const PinnedAllocator<U>&) noexcept { return false; }

#endif
