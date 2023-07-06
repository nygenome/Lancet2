#include "lancet/core/variant_store.h"

#include <algorithm>
#include <utility>

#include "lancet/base/logging.h"
#include "spdlog/fmt/ostr.h"

namespace lancet::core {

void VariantStore::AddVariants(std::vector<Value> &&variants) {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (variants.empty()) return;

  const absl::MutexLock lock(&mMutex);
  for (auto &&curr : variants) {
    const auto identifier = curr->Identifier();
    auto prev = mData.find(identifier);
    if (prev == mData.end()) {
      mData.emplace(identifier, std::move(curr));
      continue;
    }

    if (prev->second->TotalCoverage() < curr->TotalCoverage() && prev->second->Quality() < curr->Quality()) {
      mData.erase(prev);
      mData.emplace(identifier, std::move(curr));
    }
  }
}

void VariantStore::FlushVariantsBeforeWindow(const Window &win, std::ostream &out) {
  const absl::MutexLock lock(&mMutex);
  const auto variant_keys_to_extract = KeysBeforeWindow(win);
  ExtractKeysAndDumpToStream(absl::MakeConstSpan(variant_keys_to_extract), out);
}

void VariantStore::FlushAllVariantsInStore(std::ostream &out) {
  const absl::MutexLock lock(&mMutex);
  std::vector<Key> variant_keys_to_extract;
  variant_keys_to_extract.reserve(mData.size());

  std::ranges::transform(mData, std::back_inserter(variant_keys_to_extract),
                         [](const Item &item) -> Key { return item.first; });

  ExtractKeysAndDumpToStream(absl::MakeConstSpan(variant_keys_to_extract), out);
}

auto VariantStore::KeysBeforeWindow(const Window &win) const -> std::vector<Key> {
  std::vector<Key> results;
  results.reserve(mData.size());

  std::ranges::for_each(mData, [&win, &results](const Item &item) {
    const auto &var_ptr = item.second;
    const auto is_before_window = var_ptr->ChromIndex() != win.ChromIndex() ? var_ptr->ChromIndex() < win.ChromIndex()
                                                                            : var_ptr->StartPos1() < win.EndPos1();

    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (is_before_window) results.emplace_back(item.first);
  });

  return results;
}

void VariantStore::ExtractKeysAndDumpToStream(absl::Span<const Key> keys, std::ostream &out) {
  std::vector<Value> variants;
  variants.reserve(keys.size());

  using caller::RawVariant::State::NONE;
  using caller::RawVariant::Type::REF;
  static const auto has_no_support = [](const Value &item) { return item->Category() == REF || item->State() == NONE; };
  std::ranges::for_each(keys, [&variants, this](const Key &key) {
    auto handle = this->mData.extract(key);
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (handle.empty() || has_no_support(handle.mapped())) return;
    variants.emplace_back(std::move(handle.mapped()));
  });

  std::ranges::sort(variants, [](const Value &lhs, const Value &rhs) -> bool { return *lhs < *rhs; });
  std::ranges::for_each(variants, [&out](const Value &item) { fmt::print(out, "{}\n", item->AsVcfRecord()); });

  if (!variants.empty()) {
    out.flush();
    LOG_DEBUG("Flushed {} variant(s) from VariantStore to output VCF file", variants.size())
  }
}

}  // namespace lancet::core
