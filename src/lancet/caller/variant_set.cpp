#include "lancet/caller/variant_set.h"

#include <algorithm>
#include <string>
#include <utility>

#include "lancet/base/assert.h"

static constexpr char ALIGN_GAP = '-';
static constexpr usize REF_HAP_IDX = 0;

namespace {

inline auto BuildAllele(const std::string_view seq, const std::array<usize, 2> &range) -> std::string {
  const auto [start, end] = range;
  const auto allele_with_gaps = seq.substr(start, end - start + 1);

  std::string final_allele;
  final_allele.reserve(allele_with_gaps.size());
  std::ranges::copy_if(allele_with_gaps, std::back_inserter(final_allele),
                       [](const char &base) -> bool { return base != ALIGN_GAP; });

  return final_allele;
}

using lancet::caller::RawVariant;
inline auto MakeVarType(const std::array<std::string_view, 2> ref_alt) -> RawVariant::Type {
  const auto [ref, alt] = ref_alt;
  const auto ref_len = static_cast<i64>(ref.length());
  const auto alt_len = static_cast<i64>(alt.length());
  const auto diff = alt_len - ref_len;
  // NOLINTBEGIN(readability-braces-around-statements)
  if (diff == 0 && ref.size() > 1 && alt.size() > 1) return RawVariant::Type::MNP;
  if (diff == 0 && ref.size() == 1 && alt.size() == 1) return RawVariant::Type::SNV;
  if (diff < 0 && ref.size() > 1) return RawVariant::Type::DEL;
  if (diff > 0 && alt.size() > 1) return RawVariant::Type::INS;
  // NOLINTEND(readability-braces-around-statements)
  return RawVariant::Type::REF;
}

inline auto GetAlleleLength(const std::array<std::string_view, 2> ref_alt, const RawVariant::Type vtype) -> i64 {
  const auto [ref, alt] = ref_alt;
  const auto ref_len = static_cast<i64>(ref.length());
  const auto alt_len = static_cast<i64>(alt.length());
  const auto diff = alt_len - ref_len;
  return vtype == RawVariant::Type::SNV ? 1 : diff == 0 ? alt_len : diff;
}

inline auto RemoveSuperfluousBases(std::string &ref, std::string &alt) -> usize {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (ref.length() == 1 || alt.length() == 1) return 0;

  // Only necessary because MSA of multiple paths with reference can cause nested indels. Example below
  // Here REF and HAP2 will generate superfluous CTC in both REF and ALT alleles after `FindVariationRanges`
  // 'CTC-----------------------------------------------------------------------------------------------------------'
  // 'CTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATTTATCTATTTATCTATCTATCTATCATCTACCTATCTCTCTATCTGTCATCTATCTCTCTGCCTT'
  // 'C----TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATTTATCTATTTATCTATCTATCTATCATCTACCTATCTCTCTATCTGTCATCTATCTCTCTGCCTT'

  usize idx = 1;
  while ((idx < ref.size()) && (idx < alt.size()) && ref[idx] == alt[idx]) {
    idx++;
  }

  // -1 because we want first base to match to left align REF and ALT
  const auto num_superfluous_bases = idx - 1;
  ref.erase(0, num_superfluous_bases);
  alt.erase(0, num_superfluous_bases);

  return num_superfluous_bases;
}

}  // namespace

namespace lancet::caller {

VariantSet::VariantSet(const MsaBuilder &bldr, const core::Window &win) {
  const auto msa = bldr.MultipleSequenceAlignment();
  const auto num_msa_seqs = msa.size();
  LANCET_ASSERT(num_msa_seqs > 1)
  LANCET_ASSERT(std::ranges::all_of(msa, [&msa](std::string_view seq) { return seq.size() == msa[0].size(); }))

  // Walk through each pairwise REF-ALT alignment in the MSA
  for (usize alt_hap_idx = 1; alt_hap_idx < num_msa_seqs; ++alt_hap_idx) {
    const auto ref_aln = msa[REF_HAP_IDX];
    const auto alt_aln = msa[alt_hap_idx];
    const auto alt_sequence = bldr.FetchHaplotypeSeqView(alt_hap_idx);

    for (const auto &mismatch : FindVariationRanges({ref_aln, alt_aln})) {
      auto ref_allele = std::move(BuildAllele(msa[REF_HAP_IDX], mismatch));
      auto alt_allele = std::move(BuildAllele(msa[alt_hap_idx], mismatch));
      const auto num_superfluous = RemoveSuperfluousBases(ref_allele, alt_allele);

      const auto nref_gaps = std::count(ref_aln.begin(), ref_aln.begin() + mismatch[0], ALIGN_GAP);
      const auto nalt_gaps = std::count(alt_aln.begin(), alt_aln.begin() + mismatch[0], ALIGN_GAP);
      const auto start_ref0 = mismatch[0] - nref_gaps + num_superfluous;
      const auto start_alt0 = mismatch[0] - nalt_gaps + num_superfluous;

      const auto var_type = MakeVarType({ref_allele, alt_allele});
      const auto allele_length = GetAlleleLength({ref_allele, alt_allele}, var_type);
      LANCET_ASSERT(var_type != RawVariant::Type::REF)
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (var_type == RawVariant::Type::REF) continue;

      RawVariant msa_variant;
      msa_variant.mChromIndex = win.ChromIndex();
      msa_variant.mGenomeStart1 = win.StartPos1() + start_ref0;
      msa_variant.mAlleleLength = allele_length;
      msa_variant.mType = var_type;
      msa_variant.mChromName = win.ChromName();
      msa_variant.mRefAllele = std::move(ref_allele);
      msa_variant.mAltAllele = std::move(alt_allele);

      auto itr = mResultVariants.find(msa_variant);
      if (itr == mResultVariants.end()) {
        msa_variant.mHapStart0Idxs.emplace(REF_HAP_IDX, start_ref0);
        msa_variant.mHapStart0Idxs.emplace(alt_hap_idx, start_alt0);
        msa_variant.mStrResult = FindStr(alt_sequence, start_alt0);
        mResultVariants.emplace(std::move(msa_variant));
      } else {
        itr->mHapStart0Idxs.emplace(alt_hap_idx, start_alt0);
        // NOLINTNEXTLINE(readability-braces-around-statements)
        if (!itr->mStrResult.mFoundStr) itr->mStrResult = FindStr(alt_sequence, start_alt0);
      }
    }
  }
}

auto VariantSet::FindOverlappingVariants(const usize hap_idx, const StartEndRange &interval) const -> VariantPtrs {
  const auto [aln_start, aln_end] = interval;

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (aln_start < 0 || aln_end < 0) return {};

  std::vector<const RawVariant *> results;
  results.reserve(mResultVariants.size());

  for (const auto &mvar : mResultVariants) {
    const auto hap_data_itr = mvar.mHapStart0Idxs.find(hap_idx);
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (hap_data_itr == mvar.mHapStart0Idxs.cend()) continue;

    const auto vstart = hap_data_itr->second;
    const auto vend = vstart + mvar.mAlleleLength;
    if (aln_start <= vstart && aln_end >= vend) {
      results.emplace_back(&mvar);
    }
  }

  return results;
}

auto VariantSet::FindVariationRanges(const PairwiseAlignment &alignment_view) -> VariationRanges {
  const auto [ref_aln, alt_aln] = alignment_view;
  LANCET_ASSERT(ref_aln.size() == alt_aln.size())

  // Skip all left and right flanking gaps from alignment
  const auto [gap_free_start, gap_free_end] = EndGapsFreeRange(alignment_view);
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (gap_free_start >= gap_free_end || gap_free_end > ref_aln.size()) return {};

  std::vector<StartEndIndices> mismatch_ranges;
  static constexpr usize ESTIMATED_NUM_RANGES_TO_ALLOCATE = 32;
  mismatch_ranges.reserve(ESTIMATED_NUM_RANGES_TO_ALLOCATE);

  const auto *itr_ref = ref_aln.begin() + gap_free_start;
  const auto *itr_alt = alt_aln.begin() + gap_free_start;
  const auto *const end_itr_ref = ref_aln.begin() + gap_free_end;
  const auto *const end_itr_alt = alt_aln.begin() + gap_free_end;

  while (itr_ref != end_itr_ref && itr_alt != end_itr_alt) {
    // std::mismatch returns a pair of iterators to the first mismatching elements of two ranges
    auto mis_pair = std::mismatch(itr_ref, end_itr_ref, itr_alt);

    // If the iterators are not at the end of gap free alignments, a mismatch was found
    if (mis_pair.first != end_itr_ref) {
      std::string_view::iterator mis_start = mis_pair.first;
      std::string_view::iterator mis_end = mis_start;

      // Find the end of this mismatch range by searching for the next matching pair
      while (*mis_end != *mis_pair.second && mis_end != end_itr_ref) {
        ++mis_end;
        ++mis_pair.second;
      }

      // Positions are given as distances from the start of the strings, starting from 0
      auto range_start = static_cast<usize>(std::distance(ref_aln.begin(), mis_start));
      const auto range_end = static_cast<usize>(std::distance(ref_aln.begin(), mis_end) - 1);

      while (range_start > 1 && ref_aln[range_start] == ALIGN_GAP || alt_aln[range_start] == ALIGN_GAP) {
        range_start--;
      }

      // Left align the variant if we are at an InDel/MNP, so that get normalized variant range
      const auto is_indel_or_mnp = ((range_end - range_start) > 1) || (range_start != range_end);
      while (range_start > 1 && is_indel_or_mnp && ref_aln[range_start] != alt_aln[range_start]) {
        range_start--;
      }

      LANCET_ASSERT(range_start < gap_free_end)
      LANCET_ASSERT(range_end < gap_free_end)
      // Check to make sure that range_start-- doesn't wrap below 0 to SIZE_MAX
      if (range_start < gap_free_end) {
        mismatch_ranges.emplace_back(StartEndIndices{range_start, range_end});
      }

      // Move the iterators past this mismatch range for the next loop
      itr_ref = mis_end;
      itr_alt = mis_pair.second;
    } else {
      break;
    }
  }

  return mismatch_ranges;
}

auto VariantSet::EndGapsFreeRange(const PairwiseAlignment &alignment_view) -> StartEndIndices {
  const auto [ref_aln, alt_aln] = alignment_view;
  LANCET_ASSERT(ref_aln.size() == alt_aln.size())

  usize start_idx = 0;
  while (ref_aln[start_idx] == ALIGN_GAP || alt_aln[start_idx] == ALIGN_GAP) {
    ++start_idx;
  }

  usize end_idx = ref_aln.length() - 1;
  while (ref_aln[end_idx] == ALIGN_GAP || alt_aln[end_idx] == ALIGN_GAP) {
    --end_idx;
  }

  return {start_idx, end_idx};
}

}  // namespace lancet::caller
