#include "lancet2/edge.h"

#include "absl/hash/hash.h"

namespace lancet2 {
auto Edge::ID() const -> u64 { return absl::Hash<Edge>()(*this); }
}  // namespace lancet2
