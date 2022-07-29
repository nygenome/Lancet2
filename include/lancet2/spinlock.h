#ifndef LANCET_SPINLOCK_H
#define LANCET_SPINLOCK_H

#include <atomic>

namespace lancet2::utils {
// https://rigtorp.se/spinlock/
class SpinLock {
 public:
  SpinLock() = default;

  void Lock() noexcept {
    for (;;) {
      // Optimistically assume the lock is free on first the try
      if (!alock_.exchange(true, std::memory_order_acquire)) {
        return;
      }

      // Wait for lock to be released without generating cache misses
      while (alock_.load(std::memory_order_relaxed)) {
        // Issue X86 PAUSE or ARM YIELD instruction to
        // reduce contention between hyper-threads
#if defined(__x86_64__)
        __builtin_ia32_pause();
#endif /* x64 */
#if defined(__arm__) || defined(__arm64__)
        __asm__ __volatile__("yield");
#endif
      }
    }
  }

  void Unlock() noexcept { alock_.store(false, std::memory_order_release); }

  auto TryLock() noexcept -> bool {
    // First do a relaxed load to check if lock is free in order to prevent
    // unnecessary cache misses if someone does while(!try_lock())
    return !alock_.load(std::memory_order_relaxed) && !alock_.exchange(true, std::memory_order_acquire);
  }

  inline void lock() noexcept { Lock(); }
  inline void unlock() noexcept { Unlock(); }

 private:
  std::atomic<bool> alock_;
};
}  // namespace lancet2::utils

#endif  // LANCET_SPINLOCK_H
