#include "lancet/node_cov.h"

#include <algorithm>
#include <cmath>

#include "absl/types/span.h"
#include "lancet/merge_node_info.h"

namespace lancet {
void NodeCov::MergeBuddy(const NodeCov& buddy, BuddyPosition dir, bool reverse_buddy, std::size_t k) {
  const auto c1Ratio = static_cast<float>(tmrBases.size() - k + 1);
  const auto c2Ratio = static_cast<float>(buddy.Size() - k + 1);
  const auto combinedRatio = c1Ratio + c2Ratio;

  cntTumorFwd = std::ceil(((cntTumorFwd * c1Ratio) + (buddy.cntTumorFwd * c2Ratio)) / combinedRatio);     // NOLINT
  cntTumorRev = std::ceil(((cntTumorRev * c1Ratio) + (buddy.cntTumorRev * c2Ratio)) / combinedRatio);     // NOLINT
  cntNormalFwd = std::ceil(((cntNormalFwd * c1Ratio) + (buddy.cntNormalFwd * c2Ratio)) / combinedRatio);  // NOLINT
  cntNormalRev = std::ceil(((cntNormalRev * c1Ratio) + (buddy.cntNormalRev * c2Ratio)) / combinedRatio);  // NOLINT

  MergeNodeInfo(&tmrBases, absl::MakeConstSpan(buddy.tmrBases), dir, reverse_buddy, k);
  MergeNodeInfo(&nmlBases, absl::MakeConstSpan(buddy.nmlBases), dir, reverse_buddy, k);
}

auto NodeCov::StrandCov(SampleLabel label, Strand s) const -> std::uint16_t {
  if (label == SampleLabel::TUMOR) return s == Strand::FWD ? cntTumorFwd : cntTumorRev;
  return s == Strand::FWD ? cntNormalFwd : cntNormalRev;
}

auto NodeCov::TotalCov(SampleLabel label) const -> std::uint16_t {
  return StrandCov(label, Strand::FWD) + StrandCov(label, Strand::REV);
}

void NodeCov::Update(SampleLabel label, Strand s, const std::vector<std::size_t>& bq_pass) {
  if (label == SampleLabel::TUMOR) {
    if (s == Strand::FWD) {
      cntTumorFwd++;
      std::for_each(tmrBases.begin(), tmrBases.end(), [](auto& b) { b.fwdRaw++; });
      for (const auto& pos : bq_pass) tmrBases[pos].fwdBQPass += 1;
    } else {
      cntTumorRev++;
      std::for_each(tmrBases.begin(), tmrBases.end(), [](auto& b) { b.revRaw++; });
      for (const auto& pos : bq_pass) tmrBases[pos].revBQPass += 1;
    }
    return;
  }

  if (s == Strand::FWD) {
    cntNormalFwd++;
    std::for_each(nmlBases.begin(), nmlBases.end(), [](auto& b) { b.fwdRaw++; });
    for (const auto& pos : bq_pass) nmlBases[pos].fwdBQPass += 1;
  } else {
    cntNormalRev++;
    std::for_each(nmlBases.begin(), nmlBases.end(), [](auto& b) { b.revRaw++; });
    for (const auto& pos : bq_pass) nmlBases[pos].revBQPass += 1;
  }
}

void NodeCov::Update(std::uint16_t val, SampleLabel label, Strand s, const std::vector<std::size_t>& bq_pass) {
  if (label == SampleLabel::TUMOR) {
    if (s == Strand::FWD) {
      cntTumorFwd++;
      std::for_each(tmrBases.begin(), tmrBases.end(), [&val](auto& b) { b.fwdRaw = val; });
      for (auto& bcov : tmrBases) bcov.fwdRaw = val;
      for (const auto& pos : bq_pass) tmrBases[pos].fwdBQPass += 1;
    } else {
      cntTumorRev++;
      std::for_each(tmrBases.begin(), tmrBases.end(), [&val](auto& b) { b.revRaw = val; });
      for (const auto& pos : bq_pass) tmrBases[pos].revBQPass += 1;
    }
    return;
  }

  if (s == Strand::FWD) {
    cntNormalFwd++;
    std::for_each(nmlBases.begin(), nmlBases.end(), [&val](auto& b) { b.fwdRaw = val; });
    for (const auto& pos : bq_pass) nmlBases[pos].fwdBQPass += 1;
  } else {
    cntNormalRev++;
    std::for_each(nmlBases.begin(), nmlBases.end(), [&val](auto& b) { b.revRaw = val; });
    for (const auto& pos : bq_pass) nmlBases[pos].revBQPass += 1;
  }
}

void NodeCov::Update(SampleLabel label, Strand s, std::size_t pos) {
  if (label == SampleLabel::TUMOR) {
    if (s == Strand::FWD) {
      cntTumorFwd++;
      tmrBases[pos].fwdRaw += 1;
      tmrBases[pos].fwdBQPass += 1;
    } else {
      cntTumorRev++;
      tmrBases[pos].revRaw += 1;
      tmrBases[pos].revBQPass += 1;
    }
    return;
  }

  if (s == Strand::FWD) {
    cntNormalFwd++;
    nmlBases[pos].fwdRaw += 1;
    nmlBases[pos].fwdBQPass += 1;
  } else {
    cntNormalRev++;
    nmlBases[pos].revRaw += 1;
    nmlBases[pos].revBQPass += 1;
  }
}

void NodeCov::Clear() {
  cntTumorFwd = 0;
  cntTumorRev = 0;
  cntNormalFwd = 0;
  cntNormalRev = 0;
  tmrBases.clear();
  nmlBases.clear();
}

void NodeCov::Reverse() {
  std::reverse(tmrBases.begin(), tmrBases.end());
  std::reverse(nmlBases.begin(), nmlBases.end());
}
}  // namespace lancet
