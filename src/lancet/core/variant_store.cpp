#include "lancet/core/variant_store.h"

#include "lancet/base/logging.h"
#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/core/window.h"

#include "absl/hash/hash.h"
#include "absl/synchronization/mutex.h"
#include "spdlog/fmt/bundled/ostream.h"

#include <algorithm>
#include <ostream>
#include <utility>
#include <vector>

namespace lancet::core {

void VariantStore::AddVariants(std::vector<Value>&& variants) {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (variants.empty())
    return;

  for (auto&& curr : variants) {
    auto const identifier = curr->Identifier();
    usize const shard_idx = absl::Hash<Key>()(identifier) & (NUM_SHARDS - 1);

    // Independent lock specifically for this shard
    // NOLINTBEGIN(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
    absl::MutexLock const lock(mBuckets[shard_idx].mUtex);
    auto& map = mBuckets[shard_idx].mData;
    // NOLINTEND(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)

    auto prev = map.find(identifier);
    if (prev == map.end()) {
      map.emplace(identifier, std::move(curr));
      continue;
    }

    // Duplicate found (same CHROM+POS+REF locus). Keep the variant with higher
    // total coverage — the better-covered window likely assembled a more complete
    // multi-allelic picture. See identity design note on VariantCall::operator<.
    if (prev->second->TotalCoverage() < curr->TotalCoverage()) {
      prev->second = std::move(curr);
    }
  }
}

void VariantStore::FlushVariantsBeforeWindow(Window const& win, std::ostream& out) {
  std::vector<Value> variants_to_write;

  // Extract variants from each bucket completely independently
  for (auto& bucket : mBuckets) {
    absl::MutexLock const lock(bucket.mUtex);

    std::vector<Key> keys_to_extract;
    for (auto const& item : bucket.mData) {
      auto const& var_ptr = item.second;
      auto const is_before_window = var_ptr->ChromIndex() != win.ChromIndex()
                                        ? var_ptr->ChromIndex() < win.ChromIndex()
                                        : var_ptr->StartPos1() < win.EndPos1();
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (is_before_window)
        keys_to_extract.emplace_back(item.first);
    }

    using caller::RawVariant::Type::REF;
    static auto const HAS_NO_SUPPORT = [](Value const& item) -> bool {
      return !item->HasAltSupport() ||
             std::ranges::all_of(item->Categories(), [](auto call) -> bool { return call == REF; });
    };

    for (auto const& key : keys_to_extract) {
      auto handle = bucket.mData.extract(key);
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (handle.empty() || HAS_NO_SUPPORT(handle.mapped()))
        continue;
      variants_to_write.emplace_back(std::move(handle.mapped()));
    }
  }

  // Complete the sort & output without holding any locks
  std::ranges::sort(variants_to_write,
                    [](Value const& lhs, Value const& rhs) -> bool { return *lhs < *rhs; });
  std::ranges::for_each(variants_to_write, [&out](Value const& item) -> void {
    fmt::print(out, "{}\n", item->AsVcfRecord());
  });

  if (!variants_to_write.empty()) {
    out.flush();
    LOG_DEBUG("Flushed {} variant(s) from VariantStore to output VCF file",
              variants_to_write.size())
  }
}

void VariantStore::FlushAllVariantsInStore(std::ostream& out) {
  std::vector<Value> variants_to_write;

  for (auto& bucket : mBuckets) {
    absl::MutexLock const lock(bucket.mUtex);

    using caller::RawVariant::Type::REF;
    static auto const HAS_NO_SUPPORT = [](Value const& item) -> bool {
      return !item->HasAltSupport() ||
             std::ranges::all_of(item->Categories(), [](auto call) -> bool { return call == REF; });
    };

    for (auto& item : bucket.mData) {
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (HAS_NO_SUPPORT(item.second))
        continue;
      variants_to_write.emplace_back(std::move(item.second));
    }
    bucket.mData.clear();
  }

  std::ranges::sort(variants_to_write,
                    [](Value const& lhs, Value const& rhs) -> bool { return *lhs < *rhs; });
  std::ranges::for_each(variants_to_write, [&out](Value const& item) -> void {
    fmt::print(out, "{}\n", item->AsVcfRecord());
  });

  if (!variants_to_write.empty()) {
    out.flush();
    LOG_DEBUG("Flushed {} variant(s) from VariantStore to output VCF file",
              variants_to_write.size())
  }
}

}  // namespace lancet::core
