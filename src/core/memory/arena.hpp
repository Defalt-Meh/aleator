#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <stdexcept>

#include "core/memory/aligned_allocator.hpp"

namespace aleator::core {

/// Bump-pointer arena over a single aligned allocation. Intended for
/// per-step scratch buffers (neighbor-list scratch, MC trial-move buffers)
/// that are reset wholesale every iteration rather than freed piecemeal.
///
/// Not thread-safe: give each worker thread its own Arena.
class Arena {
public:
    explicit Arena(std::size_t capacityBytes, std::size_t alignment = 64)
        : capacity_(capacityBytes), alignment_(alignment) {
        if (alignment_ == 0 || (alignment_ & (alignment_ - 1)) != 0) {
            throw std::invalid_argument("Arena alignment must be a power of two");
        }
        if (capacity_ > 0) {
            buffer_ = static_cast<std::byte*>(
                ::operator new(capacity_, std::align_val_t{alignment_}));
        }
    }

    Arena(const Arena&) = delete;
    Arena& operator=(const Arena&) = delete;
    Arena(Arena&&) = delete;
    Arena& operator=(Arena&&) = delete;

    ~Arena() {
        if (buffer_ != nullptr) {
            ::operator delete(buffer_, std::align_val_t{alignment_});
        }
    }

    /// Allocates `bytes` from the arena, aligned to `alignment_`. Returns
    /// nullptr if the arena is exhausted; callers decide whether that is
    /// fatal (typical) or a signal to grow and retry.
    [[nodiscard]] std::byte* allocate(std::size_t bytes) noexcept {
        const std::size_t alignedOffset = (offset_ + alignment_ - 1) & ~(alignment_ - 1);
        if (alignedOffset + bytes > capacity_) {
            return nullptr;
        }
        offset_ = alignedOffset + bytes;
        return buffer_ + alignedOffset;
    }

    /// Resets the arena to empty without releasing the underlying
    /// allocation. Call once per simulation step.
    void reset() noexcept { offset_ = 0; }

    [[nodiscard]] std::size_t capacity() const noexcept { return capacity_; }
    [[nodiscard]] std::size_t used() const noexcept { return offset_; }

private:
    std::size_t capacity_;
    std::size_t alignment_;
    std::size_t offset_ = 0;
    std::byte* buffer_ = nullptr;
};

} // namespace aleator::core
