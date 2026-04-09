#include "lancet/core/variant_store.h"

#include <algorithm>
#include <iterator>
#include <ostream>
#include <utility>
#include <vector>

#include "absl/hash/hash.h"
#include "absl/synchronization/mutex.h"
#include "absl/types/span.h"
#include "lancet/base/logging.h"
#include "lancet/caller/raw_variant.h"
#include "spdlog/fmt/bundled/ostream.h"
#include "spdlog/fmt/ostr.h"
#include "lancet/core/window.h"

namespace lancet::core {

void VariantStore::AddVariants(std::vector<Value> &&variants) {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (variants.empty()) return;

  for (auto &&curr : variants) {
    const auto identifier = curr->Identifier();
    const usize shard_idx = absl::Hash<Key>()(identifier) & (NUM_SHARDS - 1);
    
    // Independent lock specifically for this shard
    const absl::MutexLock lock(mBuckets[shard_idx].mutex);
    auto& map = mBuckets[shard_idx].data;

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

void VariantStore::FlushVariantsBeforeWindow(const Window &win, std::ostream &out) {
  std::vector<Value> variants_to_write;

  // Extract variants from each bucket completely independently
  for (auto& bucket : mBuckets) {
    const absl::MutexLock lock(bucket.mutex);
    
    std::vector<Key> keys_to_extract;
    for (const auto& item : bucket.data) {
      const auto &var_ptr = item.second;
      const auto is_before_window = var_ptr->ChromIndex() != win.ChromIndex() ? var_ptr->ChromIndex() < win.ChromIndex()
                                                                              : var_ptr->StartPos1() < win.EndPos1();
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (is_before_window) keys_to_extract.emplace_back(item.first);
    }

    using caller::RawVariant::Type::REF;
    static const auto has_no_support = [](const Value &item) { 
      return !item->HasAltSupport() || std::ranges::all_of(item->Categories(), [](auto c) { return c == REF; }); 
    };

    for (const auto& key : keys_to_extract) {
      auto handle = bucket.data.extract(key);
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (handle.empty() || has_no_support(handle.mapped())) continue;
      variants_to_write.emplace_back(std::move(handle.mapped()));
    }
  }

  // Complete the sort & output without holding any locks
  std::ranges::sort(variants_to_write, [](const Value &lhs, const Value &rhs) -> bool { return *lhs < *rhs; });
  std::ranges::for_each(variants_to_write, [&out](const Value &item) { fmt::print(out, "{}\n", item->AsVcfRecord()); });

  if (!variants_to_write.empty()) {
    out.flush();
    LOG_DEBUG("Flushed {} variant(s) from VariantStore to output VCF file", variants_to_write.size())
  }
}

void VariantStore::FlushAllVariantsInStore(std::ostream &out) {
  std::vector<Value> variants_to_write;

  for (auto& bucket : mBuckets) {
    const absl::MutexLock lock(bucket.mutex);
    
    using caller::RawVariant::Type::REF;
    static const auto has_no_support = [](const Value &item) { 
      return !item->HasAltSupport() || std::ranges::all_of(item->Categories(), [](auto c) { return c == REF; }); 
    };

    for (auto& item : bucket.data) {
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (has_no_support(item.second)) continue;
      variants_to_write.emplace_back(std::move(item.second));
    }
    bucket.data.clear();
  }

  std::ranges::sort(variants_to_write, [](const Value &lhs, const Value &rhs) -> bool { return *lhs < *rhs; });
  std::ranges::for_each(variants_to_write, [&out](const Value &item) { fmt::print(out, "{}\n", item->AsVcfRecord()); });

  if (!variants_to_write.empty()) {
    out.flush();
    LOG_DEBUG("Flushed {} variant(s) from VariantStore to output VCF file", variants_to_write.size())
  }
}

}  // namespace lancet::core
