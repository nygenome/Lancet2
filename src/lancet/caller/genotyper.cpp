#include "lancet/caller/genotyper.h"

#include <algorithm>
#include <cmath>
#include <string_view>

#include "absl/strings/numbers.h"
#include "lancet/base/assert.h"
#include "lancet/base/online_stats.h"
#include "mmpriv.h"
#include "spdlog/fmt/fmt.h"

namespace {

inline void FreeMinimap2Alignment(mm_reg1_t* regs, const int num_regs) {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (regs == nullptr) return;
  // NOLINTBEGIN(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc)
  for (int idx = 0; idx < num_regs; ++idx) {
    auto* curr_reg = &regs[idx];  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    std::free(curr_reg->p);
  }
  std::free(regs);
  // NOLINTEND(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc)
}

}  // namespace

namespace lancet::caller {

Genotyper::Genotyper() {
  // 0 -> no info, 1 -> error, 2 -> warning, 3 -> debug
  mm_verbose = 1;
  mm_set_opt(nullptr, mIndexingOpts.get(), mMappingOpts.get());
  mm_set_opt("sr", mIndexingOpts.get(), mMappingOpts.get());
  mMappingOpts->flag |= MM_F_CIGAR | MM_F_OUT_CS;
  mMappingOpts->best_n = 1;
}

auto Genotyper::Genotype(Haplotypes haplotypes, Reads reads, const VariantSet& vset) -> Result {
  ResetData(haplotypes);

  Result genotyped_variants;
  static constexpr usize DEFAULT_EXPECTED_SAMPLES_COUNT = 2;
  // Just in case caller of this method forgets to call SetNumSamples before
  genotyped_variants.reserve(std::max(mNumSamples, DEFAULT_EXPECTED_SAMPLES_COUNT));

  // so that we don't double count support
  // for REF and ALT alleles from same read
  AlnInfo::SupportsInfo read_supports;
  read_supports.reserve(vset.Count());

  static const auto by_descending_identity_and_score = [](const AlnInfo& lhs, const AlnInfo& rhs) {
    // Sort by gap compressed identity, then by alignment DP score and then by haplotype index
    // So if gap compressed identity and DP score are same, we pick ALT haplotype preferentially
    // NOLINTBEGIN(readability-braces-around-statements)
    if (lhs.mGcIden != rhs.mGcIden) return lhs.mGcIden > rhs.mGcIden;
    if (lhs.mDpScore != rhs.mDpScore) return lhs.mDpScore > rhs.mDpScore;
    return lhs.mHapIdx > rhs.mHapIdx;
    // NOLINTEND(readability-braces-around-statements)
  };

  for (const auto& read : reads) {
    read_supports.clear();
    auto alns_to_all_haps = AlignRead(read);
    std::ranges::sort(alns_to_all_haps, by_descending_identity_and_score);
    std::ranges::for_each(alns_to_all_haps, [&read_supports, &vset](const AlnInfo& item) {
      item.AddSupportingInfo(read_supports, vset);
    });

    // Add read supports to all the variants supported by it
    AddToTable(genotyped_variants, read, read_supports);
  }

  return genotyped_variants;
}

auto Genotyper::AlnInfo::IsEmpty() const noexcept -> bool {
  return mRefStart == -1 && mQryStart == -1 && mRefEnd == -1 && mQryEnd == -1 && mDpScore == -1 && mGcIden == 0.0 &&
         mCsTag.empty();
}

void Genotyper::AlnInfo::AddSupportingInfo(SupportsInfo& supports, const VariantSet& vset) const {
  const auto curr_allele = mHapIdx == REF_HAP_IDX ? Allele::REF : Allele::ALT;
  const auto hap_and_rd_identity_ranges = FindIdentityRanges();

  for (const auto& olap_var : vset.FindOverlappingVariants(mHapIdx, {mRefStart, mRefEnd})) {
    // If read has already been counted as support for this variant as support
    // skip the CS tag check to confirm that the read has exact match to the allele
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (supports.contains(olap_var)) continue;

    const auto al_start = olap_var->mHapStart0Idxs.at(mHapIdx);
    const auto al_len = curr_allele == Allele::REF ? olap_var->mRefAllele.length() : olap_var->mAltAllele.length();
    const auto allele_range = StartEndIndices({al_start, al_start + al_len - 1});
    const auto rd_start_idx_supporting_allele = FindQueryStart(hap_and_rd_identity_ranges, allele_range);
    if (rd_start_idx_supporting_allele) {
      auto qstart_strand = std::make_pair(rd_start_idx_supporting_allele.value(), curr_allele);
      supports.emplace(olap_var, std::move(qstart_strand));
    }
  }
}

auto Genotyper::AlnInfo::FindIdentityRanges() const -> RefQryIdentityRanges {
  std::vector<StartEndIndices> ref_iden_ranges;
  std::vector<StartEndIndices> qry_iden_ranges;

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mCsTag.empty()) return RefQryIdentityRanges({ref_iden_ranges, qry_iden_ranges});

  ref_iden_ranges.reserve(mCsTag.length());
  qry_iden_ranges.reserve(mCsTag.length());

  auto curr_ref_idx = static_cast<usize>(mRefStart);
  auto curr_qry_idx = static_cast<usize>(mQryStart);

  // Example CS Tag --> `:6-ata:10+gtc:4*at:3`
  std::string op_info;
  usize parsed_num = 0;
  std::vector<usize> idxs_with_ops;
  op_info.reserve(mCsTag.length());
  idxs_with_ops.reserve(mCsTag.length());

  // IDENTITY = ':', MISMATCH = '*', INSERTION = '+', DELETION = '-'
  const auto handle_cs_op = [&](const char token) {
    if (token == ':' && absl::SimpleAtoi(op_info, &parsed_num)) {
      ref_iden_ranges.emplace_back(StartEndIndices{curr_ref_idx, curr_ref_idx + parsed_num});
      qry_iden_ranges.emplace_back(StartEndIndices{curr_qry_idx, curr_qry_idx + parsed_num});
      curr_ref_idx += parsed_num;
      curr_qry_idx += parsed_num;
      return;
    }

    if (token == '*') {
      curr_ref_idx++;
      curr_qry_idx++;
      return;
    }

    if (token == '+') {
      curr_qry_idx += op_info.length();
      return;
    }

    if (token == '-') {
      curr_ref_idx += op_info.length();
      return;
    }
  };

  for (usize idx = 0; idx < mCsTag.length(); ++idx) {
    switch (mCsTag[idx]) {
      case ':':
      case '-':
      case '+':
      case '*':
        if (!idxs_with_ops.empty()) {
          handle_cs_op(mCsTag[idxs_with_ops.back()]);
          op_info.clear();
        }

        idxs_with_ops.emplace_back(idx);
        break;
      default:
        op_info.push_back(mCsTag[idx]);
        break;
    }
  }

  // Handle the final CS operation
  handle_cs_op(mCsTag[idxs_with_ops.back()]);
  op_info.clear();

  return RefQryIdentityRanges({ref_iden_ranges, qry_iden_ranges});
}

auto Genotyper::AlnInfo::FindQueryStart(const RefQryIdentityRanges& ref_qry_equal_ranges,
                                        const StartEndIndices& allele_span) const -> std::optional<usize> {
  // First, find if the variant allele range is contained within any of the ref/haplotype match ranges
  // Second, check if 100% of the qry is contained within the variant allele i.e variant longer than read
  // If either of these two scenarios happen, we return the start idx in read which supports variant

  const auto& [hap_identity_ranges, read_identity_ranges] = ref_qry_equal_ranges;
  const auto [var_allele_start, var_allele_end] = allele_span;
  LANCET_ASSERT(hap_identity_ranges.size() == read_identity_ranges.size())

  for (usize idx = 0; idx < hap_identity_ranges.size(); ++idx) {
    const auto [aln_hap_match_start, aln_hap_match_end] = hap_identity_ranges[idx];
    const auto [read_match_start, read_match_end] = read_identity_ranges[idx];
    // Check if alignment has exact match with variant allele within it
    if (aln_hap_match_start <= var_allele_start && aln_hap_match_end >= var_allele_end) {
      const auto allele_to_hap_match_start_diff = var_allele_start - aln_hap_match_start;
      return read_match_start + allele_to_hap_match_start_diff;
    }

    // Check if 100% of alignment length is contained within variant allele span
    // This will only happen if variant is longer than the read length itself
    const auto full_read_match = (read_match_end - read_match_start) == mQryLen && mGcIden == 1.0;
    if (full_read_match && var_allele_start <= aln_hap_match_start && var_allele_end >= aln_hap_match_end) {
      return read_identity_ranges[idx][0];
    }
  }

  return std::nullopt;
}

void Genotyper::ResetData(Haplotypes sequences) {
  mIndices.clear();
  mIndices.reserve(sequences.size() + 1);

  const auto* iopts = mIndexingOpts.get();
  std::ranges::for_each(sequences, [this, &iopts](const std::string& seq) {
    const char* raw_seq = seq.c_str();
    auto* idx_result = mm_idx_str(iopts->w, iopts->k, 0, iopts->bucket_bits, 1, &raw_seq, nullptr);
    this->mIndices.emplace_back(Minimap2Index(idx_result));
  });

  auto* mopts = mMappingOpts.get();
  std::ranges::for_each(mIndices, [&mopts](const Minimap2Index& mm2_idx) { mm_mapopt_update(mopts, mm2_idx.get()); });
}

auto Genotyper::AlignRead(const cbdg::Read& read) -> std::vector<AlnInfo> {
  std::vector<AlnInfo> results;
  results.reserve(mIndices.size());

  int nregs = 0;
  auto* tbuffer = mThreadBuffer.get();
  const auto* map_opts = mMappingOpts.get();
  const auto read_len = static_cast<int>(read.Length());

  for (usize idx = 0; idx < mIndices.size(); ++idx) {
    AlnInfo aln_info;

    const auto* hap_mm_idx = mIndices[idx].get();
    auto* regs = mm_map(hap_mm_idx, read_len, read.SeqPtr(), &nregs, tbuffer, map_opts, read.QnamePtr());
    if (regs == nullptr || nregs <= 0) {
      results.emplace_back(aln_info);
      FreeMinimap2Alignment(regs, nregs);
      continue;
    }

    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const mm_reg1_t* top_hit = &regs[0];
    aln_info.mRefStart = top_hit->rs;
    aln_info.mQryStart = top_hit->qs;
    aln_info.mRefEnd = top_hit->re;
    aln_info.mQryEnd = top_hit->qe;
    aln_info.mDpScore = top_hit->score;
    aln_info.mGcIden = mm_event_identity(top_hit);
    aln_info.mHapIdx = idx;
    aln_info.mQryLen = read.Length();

    int max_len = 0;
    char* cs_result_ptr = nullptr;
    const auto len_cs = mm_gen_cs(tbuffer->km, &cs_result_ptr, &max_len, hap_mm_idx, top_hit, read.SeqPtr(), 1);
    if (len_cs > 0 && cs_result_ptr != nullptr) {
      aln_info.mCsTag = std::string_view(cs_result_ptr, static_cast<usize>(len_cs));
      std::free(cs_result_ptr);  // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc)
    }

    results.emplace_back(std::move(aln_info));
    FreeMinimap2Alignment(regs, nregs);

    // If exact match with REF haplotype, skip aligning with ALTs
    const AlnInfo& added_aln_info = results.back();
    if (idx == REF_HAP_IDX && added_aln_info.mGcIden == 1.0) {
      break;
    }
  }

  return results;
}

void Genotyper::AddToTable(Result& result, const cbdg::Read& read, const SupportsInfo& read_supports) {
  const auto quals = read.QualView();
  const auto sample_name = read.SampleName();
  const auto read_strand = read.BitwiseFlag().IsFwdStrand() ? Strand::FWD : Strand::REV;

  for (const auto& [var_ptr, qry_start_and_allele] : read_supports) {
    auto& variant_evidence = result.try_emplace(var_ptr, PerSampleVariantEvidence()).first->second;
    auto& sample_variant = variant_evidence.try_emplace(sample_name, std::make_unique<VariantSupport>()).first->second;

    const auto [read_start_idx0, allele] = qry_start_and_allele;
    const auto allele_len = allele == Allele::REF ? var_ptr->mRefAllele.length() : var_ptr->mAltAllele.length();
    const auto allele_qualities = quals.subspan(read_start_idx0, allele_len);
    const auto mean_allele_qual = static_cast<u8>(std::round(OnlineStats::Mean(allele_qualities)));
    sample_variant->AddQual(mean_allele_qual, std::make_pair(allele, read_strand));
  }
}

}  // namespace lancet::caller
