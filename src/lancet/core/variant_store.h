#ifndef SRC_LANCET_CORE_VARIANT_STORE_H_
#define SRC_LANCET_CORE_VARIANT_STORE_H_

#include <iosfwd>
#include <memory>
#include <utility>
#include <vector>

#include "absl/base/thread_annotations.h"
#include "absl/container/flat_hash_map.h"
#include "absl/synchronization/mutex.h"
#include "absl/types/span.h"
#include "lancet/caller/variant_call.h"
#include "lancet/core/window.h"

namespace lancet::core {

class VariantStore {
 public:
  using Key = caller::VariantID;
  using Value = std::unique_ptr<caller::VariantCall>;
  using Item = std::pair<const Key, Value>;

  VariantStore() = default;

  void AddVariants(std::vector<Value>&& variants) ABSL_LOCKS_EXCLUDED(mMutex);
  void FlushVariantsBeforeWindow(const Window& win, std::ostream& out) ABSL_LOCKS_EXCLUDED(mMutex);
  void FlushAllVariantsInStore(std::ostream& out) ABSL_LOCKS_EXCLUDED(mMutex);

 private:
  absl::Mutex mMutex;
  absl::flat_hash_map<Key, Value> mData ABSL_GUARDED_BY(mMutex);

  [[nodiscard]] ABSL_SHARED_LOCKS_REQUIRED(mMutex) auto KeysBeforeWindow(const Window& win) const -> std::vector<Key>;
  void ExtractKeysAndDumpToStream(absl::Span<const Key> keys, std::ostream& out) ABSL_EXCLUSIVE_LOCKS_REQUIRED(mMutex);
};

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_VARIANT_STORE_H_
