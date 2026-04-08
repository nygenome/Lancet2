#ifndef SRC_LANCET_CORE_VARIANT_STORE_H_
#define SRC_LANCET_CORE_VARIANT_STORE_H_

#include <array>
#include <iosfwd>
#include <memory>
#include <utility>
#include <vector>

#include "absl/base/thread_annotations.h"
#include "absl/container/flat_hash_map.h"
#include "absl/synchronization/mutex.h"
#include "lancet/base/types.h"
#include "lancet/caller/variant_call.h"
#include "lancet/core/window.h"

namespace lancet::core {

class VariantStore {
 public:
  using Key = caller::VariantID;
  using Value = std::unique_ptr<caller::VariantCall>;
  using Item = std::pair<const Key, Value>;

  VariantStore() = default;

  void AddVariants(std::vector<Value>&& variants);
  void FlushVariantsBeforeWindow(const Window& win, std::ostream& out);
  void FlushAllVariantsInStore(std::ostream& out);

 private:
  struct alignas(64) VariantBucket {
    mutable absl::Mutex mutex;
    absl::flat_hash_map<Key, Value> data ABSL_GUARDED_BY(mutex);
  };

  static constexpr usize NUM_SHARDS = 256;
  std::array<VariantBucket, NUM_SHARDS> mBuckets;
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_VARIANT_STORE_H_
