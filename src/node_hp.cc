#include "lancet2/node_hp.h"

#include <algorithm>

#include "absl/types/span.h"
#include "lancet2/merge_node_info.h"

namespace lancet2 {
NodeHP::NodeHP(const NodeCov& node_cov) {
  tmrHPs.reserve(node_cov.Size());
  nmlHPs.reserve(node_cov.Size());

  const auto tmrBases = node_cov.BaseCovs(SampleLabel::TUMOR);
  const auto nmlBases = node_cov.BaseCovs(SampleLabel::NORMAL);

  for (const auto& bcov : tmrBases) tmrHPs.emplace_back(MakeDefaultHP(bcov));
  for (const auto& bcov : nmlBases) nmlHPs.emplace_back(MakeDefaultHP(bcov));
}

void NodeHP::MergeBuddy(const NodeHP& buddy, BuddyPosition dir, bool reverse_buddy, usize k) {
  MergeNodeInfo(&tmrHPs, absl::MakeConstSpan(buddy.tmrHPs), dir, reverse_buddy, k);
  MergeNodeInfo(&nmlHPs, absl::MakeConstSpan(buddy.nmlHPs), dir, reverse_buddy, k);
}

void NodeHP::Update(usize hp, SampleLabel label, const std::vector<bool>& bq_pass) {
  if (label == SampleLabel::TUMOR) {
    for (usize idx = 0; idx < tmrHPs.size(); ++idx) {
      tmrHPs[idx].at(hp).raw++;
      if (bq_pass[idx]) tmrHPs[idx][hp].bqPass++;
    }
  } else {
    for (usize idx = 0; idx < nmlHPs.size(); ++idx) {
      nmlHPs[idx].at(hp).raw++;
      if (bq_pass[idx]) nmlHPs[idx][hp].bqPass++;
    }
  }
}

void NodeHP::Update(usize hp, SampleLabel label, usize base_position) {
  if (label == SampleLabel::TUMOR) {
    tmrHPs[base_position][hp].bqPass += 1;
    std::for_each(tmrHPs.begin(), tmrHPs.end(), [&hp](auto& base) { base.at(hp).raw++; });
  } else {
    nmlHPs[base_position][hp].bqPass += 1;
    std::for_each(nmlHPs.begin(), nmlHPs.end(), [&hp](auto& base) { base.at(hp).raw++; });
  }
}

void NodeHP::Clear() {
  tmrHPs.clear();
  nmlHPs.clear();
}

void NodeHP::Reverse() {
  std::reverse(tmrHPs.begin(), tmrHPs.end());
  std::reverse(nmlHPs.begin(), nmlHPs.end());
}
}  // namespace lancet2
