#include "lancet/caller/genotyper.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

extern "C" {
#include "minimap.h"
#include "mmpriv.h"
}

#include "absl/strings/numbers.h"
#include "lancet/base/assert.h"
#include "lancet/base/compute_stats.h"
#include "lancet/base/hash.h"
#include "lancet/base/types.h"
#include "lancet/caller/variant_set.h"
#include "lancet/caller/variant_support.h"
#include "lancet/hts/cigar_unit.h"

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

Genotyper::Genotyper(const Preset preset) {
  // 0 -> no info, 1 -> error, 2 -> warning, 3 -> debug
  // set default parameters first before seeting preset parameters
  mm_verbose = 1;
  mm_set_opt(nullptr, mIndexingOpts.get(), mMappingOpts.get());
  mm_set_opt(preset == Preset::ShortRead ? "sr" : "map-ont", mIndexingOpts.get(), mMappingOpts.get());
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
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (!read.PassesAlnFilters()) continue;

    read_supports.clear();
    auto alns_to_all_haps = AlignRead(read);
    std::ranges::sort(alns_to_all_haps, by_descending_identity_and_score);
    std::ranges::for_each(alns_to_all_haps, [&read_supports, &vset](const AlnInfo& item) {
      item.AddSupportingInfo(read_supports, vset);
    });

    AddToTable(genotyped_variants, read, read_supports);
  }

  return genotyped_variants;
}

auto Genotyper::AlnInfo::IsEmpty() const noexcept -> bool {
  return mRefStart == -1 && mQryStart == -1 && mRefEnd == -1 && mQryEnd == -1 && mDpScore == -1 && mGcIden == 0.0 &&
         mHapIdx == 0 && mQryLen == 0 && mCsTag.empty();
}

auto Genotyper::AlnInfo::IsFullQueryMatch() const noexcept -> bool {
  return (mQryEnd - mQryStart) == static_cast<i32>(mQryLen) && mGcIden == 1.0;
}

void Genotyper::AlnInfo::AddSupportingInfo(SupportsInfo& supports, const VariantSet& called_variants) const {
  const auto curr_allele = mHapIdx == REF_HAP_IDX ? Allele::REF : Allele::ALT;
  const auto identity_ranges = FindIdentityRanges();
  const auto non_indel_chunks = FindNonIndelChunks();

  for (const auto& variant : called_variants) {
    // If read has already been counted as support for this variant as support
    // skip the CS tag check to confirm that the read has exact match to the allele
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (!variant.mHapStart0Idxs.contains(mHapIdx) || supports.contains(&variant)) continue;

    const auto al_start = variant.mHapStart0Idxs.at(mHapIdx);
    const auto al_len = std::max(variant.mRefAllele.length(), variant.mAltAllele.length());
    const auto al_range = StartEndIndices({al_start, al_start + al_len - 1});
    const auto rd_start_idx_supporting_allele = FindQueryStartForAllele(identity_ranges, non_indel_chunks, al_range);
    if (rd_start_idx_supporting_allele) {
      auto qstart_strand = std::make_pair(rd_start_idx_supporting_allele.value(), curr_allele);
      supports.emplace(&variant, std::move(qstart_strand));
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
auto Genotyper::AlnInfo::FindNonIndelChunks() const -> NonIndelChunks {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mCsTag.empty()) return NonIndelChunks{};

  NonIndelChunks result_chunks;
  result_chunks.reserve(mCsTag.length());

  // Example CS Tag --> `:6-ata:10+gtc:4*at:3`
  // IDENTITY = ':', MISMATCH = '*', INSERTION = '+', DELETION = '-'
  usize parsed_number = 0;
  std::string op_len_data;
  std::vector<hts::CigarOp> cig_ops;
  std::vector<usize> cig_lens;
  cig_ops.reserve(mCsTag.length());
  cig_lens.reserve(mCsTag.length());

  const auto parse_operation_length = [&op_len_data, &parsed_number, &cig_lens]() {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (op_len_data.empty()) return;

    if (absl::SimpleAtoi(op_len_data, &parsed_number)) {
      cig_lens.push_back(parsed_number);
    } else {
      cig_lens.push_back(op_len_data.length());
    }
  };

  for (const char token : mCsTag) {
    switch (token) {
      case ':':
        parse_operation_length();
        cig_ops.emplace_back(hts::CigarOp::ALIGNMENT_MATCH);
        op_len_data.clear();
        break;

      case '-':
        parse_operation_length();
        cig_ops.emplace_back(hts::CigarOp::DELETION);
        op_len_data.clear();
        break;

      case '+':
        parse_operation_length();
        cig_ops.emplace_back(hts::CigarOp::INSERTION);
        op_len_data.clear();
        break;

      case '*':
        parse_operation_length();
        cig_ops.emplace_back(hts::CigarOp::SEQUENCE_MISMATCH);
        op_len_data.clear();
        break;

      default:
        op_len_data.push_back(token);
        break;
    }
  }

  // Parse the last operation length
  parse_operation_length();

  auto curr_ref_idx = static_cast<usize>(mRefStart);
  auto curr_qry_idx = static_cast<usize>(mQryStart);
  LANCET_ASSERT(cig_ops.size() == cig_lens.size())

  const auto update_ref_qry_idxs = [&curr_ref_idx, &curr_qry_idx](const hts::CigarOp cig_op, const usize len) {
    if (cig_op == hts::CigarOp::ALIGNMENT_MATCH || cig_op == hts::CigarOp::SEQUENCE_MISMATCH) {
      curr_ref_idx += len;
      curr_qry_idx += len;
      return std::array<StartEndIndices, 2>{StartEndIndices{curr_ref_idx - len, curr_ref_idx},
                                            StartEndIndices{curr_qry_idx - len, curr_qry_idx}};
    }

    if (cig_op == hts::CigarOp::DELETION) {
      curr_ref_idx += len;
      return std::array<StartEndIndices, 2>{StartEndIndices{curr_ref_idx - len, curr_ref_idx},
                                            StartEndIndices{curr_qry_idx, curr_qry_idx}};
    }

    if (cig_op == hts::CigarOp::INSERTION) {
      curr_qry_idx += len;
      return std::array<StartEndIndices, 2>{StartEndIndices{curr_ref_idx, curr_ref_idx},
                                            StartEndIndices{curr_qry_idx - len, curr_qry_idx}};
    }

    return std::array<StartEndIndices, 2>{StartEndIndices{curr_ref_idx, curr_ref_idx},
                                          StartEndIndices{curr_qry_idx, curr_qry_idx}};
  };

  for (usize idx = 0; idx < cig_ops.size(); ++idx) {
    // Skip adding indel chunks to the result chunks
    if (cig_ops[idx] == hts::CigarOp::INSERTION || cig_ops[idx] == hts::CigarOp::DELETION) {
      update_ref_qry_idxs(cig_ops[idx], cig_lens[idx]);
      continue;
    }

    const auto is_match = cig_ops[idx] == hts::CigarOp::ALIGNMENT_MATCH;
    const auto [ref_range, qry_range] = update_ref_qry_idxs(cig_ops[idx], cig_lens[idx]);

    if (result_chunks.empty()) {
      result_chunks.emplace_back(RefQryAlnChunk{ref_range, qry_range, is_match ? cig_lens[idx] : 0});
      continue;
    }

    // If previous operation is non indel, then add to previous existing chunk, instead of creating a new one
    if (cig_ops[idx - 1] == hts::CigarOp::ALIGNMENT_MATCH || cig_ops[idx - 1] == hts::CigarOp::SEQUENCE_MISMATCH) {
      result_chunks.back().mRefRange[1] = ref_range[1];
      result_chunks.back().mQryRange[1] = qry_range[1];
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (is_match) result_chunks.back().mNumExactMatches += cig_lens[idx];
      continue;
    }

    result_chunks.emplace_back(RefQryAlnChunk{ref_range, qry_range, is_match ? cig_lens[idx] : 0});
  }

  return result_chunks;
}

auto Genotyper::AlnInfo::FindQueryStartForAllele(const RefQryIdentityRanges& ref_qry_equal_ranges,
                                                 const NonIndelChunks& ref_qry_non_indel_ranges,
                                                 const StartEndIndices& allele_span) const -> std::optional<usize> {
  const auto& [hap_identity_ranges, read_identity_ranges] = ref_qry_equal_ranges;
  const auto [var_allele_start, var_allele_end] = allele_span;
  const auto one_third_read_length = static_cast<usize>(0.30 * f64(mQryLen));
  LANCET_ASSERT(hap_identity_ranges.size() == read_identity_ranges.size())

  for (usize idx = 0; idx < hap_identity_ranges.size(); ++idx) {
    const auto [aln_hap_match_start, aln_hap_match_end] = hap_identity_ranges[idx];
    const auto [read_match_start, read_match_end] = read_identity_ranges[idx];

    // For genotyping small variant alleles, where variant is within read matches
    // Check if variant allele is contained within the alignment match span
    if (aln_hap_match_start < var_allele_start && aln_hap_match_end > var_allele_end) {
      const auto allele_to_hap_match_start_diff = var_allele_start - aln_hap_match_start;
      return read_match_start + allele_to_hap_match_start_diff;
    }

    // For genotyping longer variant alleles, where read matches are within variant
    // Check if atleast 30% of the read matches the haplotype sequence exactly and
    // the entire match portion of the read lies within the variant allele span
    // This will usually only happen for long alleles or greater then read length.
    const auto partial_read_hap_match = (read_match_end - read_match_start) >= one_third_read_length;
    if (partial_read_hap_match && var_allele_start <= aln_hap_match_start && var_allele_end >= aln_hap_match_end) {
      return read_match_start;
    }
  }

  static constexpr usize LONG_ALLELE_THRESHOLD = 50;
  static constexpr f64 MIN_REQUIRED_MATCH_PERCENT = 0.9;

  // If long allele, allow fuzzy comparison of read & allele. i.e Match & Mismatch
  const auto var_length = var_allele_end - var_allele_start + 1;
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (var_length < LONG_ALLELE_THRESHOLD) return std::nullopt;

  for (const auto& non_indel_chunk : ref_qry_non_indel_ranges) {
    const auto chunk_len = static_cast<f64>(non_indel_chunk.mQryRange[1] - non_indel_chunk.mQryRange[0] + 1);
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (chunk_len < static_cast<f64>(one_third_read_length)) continue;

    const auto min_needed_matches = static_cast<usize>(std::ceil(chunk_len * MIN_REQUIRED_MATCH_PERCENT));
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (non_indel_chunk.mNumExactMatches < min_needed_matches) continue;

    const auto [aln_hap_non_indel_start, aln_hap_non_indel_end] = non_indel_chunk.mRefRange;
    const auto [read_non_indel_start, read_non_indel_end] = non_indel_chunk.mQryRange;

    // Check if variant allele is contained within the alignment non indel span
    if (aln_hap_non_indel_start < var_allele_start && aln_hap_non_indel_end > var_allele_end) {
      const auto allele_to_hap_non_indel_start_diff = var_allele_start - aln_hap_non_indel_start;
      return read_non_indel_start + allele_to_hap_non_indel_start_diff;
    }

    // Check if alignment non indel span is contained within variant allele
    if (var_allele_start <= aln_hap_non_indel_start && var_allele_end >= aln_hap_non_indel_end) {
      return read_non_indel_start;
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
    if (idx == REF_HAP_IDX && added_aln_info.IsFullQueryMatch()) {
      break;
    }
  }

  return results;
}

void Genotyper::AddToTable(Result& rslt, const cbdg::Read& read, const SupportsInfo& supports) {
  const auto quals = read.QualView();
  const auto sample_name = read.SampleName();
  const auto rname_hash = HashStr32(read.QnameView());
  const auto read_strand = read.BitwiseFlag().IsFwdStrand() ? Strand::FWD : Strand::REV;

  for (const auto& [var_ptr, qry_start_and_allele] : supports) {
    auto& variant_evidence = rslt.try_emplace(var_ptr, PerSampleVariantEvidence()).first->second;
    auto& sample_variant = variant_evidence.try_emplace(sample_name, std::make_unique<VariantSupport>()).first->second;

    const auto [read_start_idx0, allele] = qry_start_and_allele;
    const auto allele_len = allele == Allele::REF ? var_ptr->mRefAllele.length() : var_ptr->mAltAllele.length();
    const auto allele_qual = static_cast<u8>(Mean(quals.subspan(read_start_idx0, allele_len)));
    sample_variant->AddEvidence(rname_hash, allele, read_strand, allele_qual, read.MapQual(), read.PctAlnScoresDiff());
  }
}

}  // namespace lancet::caller
