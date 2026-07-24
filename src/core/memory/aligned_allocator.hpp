#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <new>

namespace aleator::core {

/// Allocator that overaligns to `Alignment` bytes, for SIMD-friendly SoA
/// buffers (Highway vector loads/stores want at least the natural vector
/// width, commonly 64 bytes for AVX-512).
///
/// Uses sized, aligned `::operator new`/`::operator delete` (C++17) rather
/// than `posix_memalign` or `_aligned_malloc`, per the portability contract
/// in CLAUDE.md #6: no POSIX-only APIs, one code path for every target
/// compiler (GCC/Clang/MSVC) and platform.
template <typename T, std::size_t Alignment = 64>
class AlignedAllocator {
    static_assert(Alignment >= alignof(T), "Alignment must be at least alignof(T)");
    static_assert((Alignment & (Alignment - 1)) == 0, "Alignment must be a power of two");

public:
    using value_type = T;

    template <typename U>
    struct rebind {
        using other = AlignedAllocator<U, Alignment>;
    };

    AlignedAllocator() noexcept = default;

    template <typename U>
    explicit AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

    [[nodiscard]] T* allocate(std::size_t n) {
        if (n > (std::numeric_limits<std::size_t>::max)() / sizeof(T)) {
            throw std::bad_array_new_length();
        }
        const std::size_t bytes = n * sizeof(T);
        return static_cast<T*>(::operator new(bytes, std::align_val_t{Alignment}));
    }

    void deallocate(T* ptr, std::size_t) noexcept {
        ::operator delete(ptr, std::align_val_t{Alignment});
    }

    template <typename U>
    bool operator==(const AlignedAllocator<U, Alignment>&) const noexcept {
        return true;
    }

    template <typename U>
    bool operator!=(const AlignedAllocator<U, Alignment>&) const noexcept {
        return false;
    }
};

} // namespace aleator::core
