#include "lancet2/edge.h"

#include "absl/hash/hash.h"

namespace lancet2 {
auto Edge::ID() const -> std::uint64_t { return absl::Hash<Edge>()(*this); }
}  // namespace lancet2
