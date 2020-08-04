#include "lancet/edge.h"

#include "absl/hash/hash.h"

namespace lancet {
auto Edge::ID() const -> std::uint64_t { return absl::Hash<Edge>()(*this); }
}  // namespace lancet
