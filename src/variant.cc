#include "lancet2/variant.h"

#include <vector>

#include "absl/hash/internal/city.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "lancet2/assert_macro.h"
#include "lancet2/fisher_exact.h"
#include "lancet2/utils.h"

namespace lancet2 {
Variant::Variant(const Transcript& transcript, usize kmer_size, VariantHpCov tmrCov, VariantHpCov nmlCov)
    : ChromName(transcript.ChromName()), Position(transcript.Position()), RefAllele(transcript.RefSeq()),
      AltAllele(transcript.AltSeq()), Kind(transcript.Code()), STRResult(transcript.STRResult()),
      Length(transcript.GetVariantLength()), KmerSize(kmer_size), TumorCov(tmrCov), NormalCov(nmlCov) {
  LANCET_ASSERT(Kind != TranscriptCode::REF_MATCH && transcript.IsFinalized());  // NOLINT
}

auto Variant::MakeVcfLine(const CliParams& params) const -> std::string {
  const auto somaticScore = PhredFisherScore(NormalCov.TotalRefCov(), TumorCov.TotalRefCov(), NormalCov.TotalAltCov(),
                                             TumorCov.TotalAltCov());

  const auto strandBiasScore = PhredFisherScore(TumorCov.RefAllele.FwdCov, TumorCov.RefAllele.RevCov,
                                                TumorCov.AltAllele.FwdCov, TumorCov.AltAllele.RevCov);

  const auto varState = ComputeState();
  LANCET_ASSERT(varState != VariantState::NONE);  // NOLINT

  auto info = absl::StrFormat("%s;FETS=%f;TYPE=%s;LEN=%d;KMERSIZE=%d;SB=%f;TMR_ALT_QUAL=%f", ToString(varState),
                              somaticScore, ToString(Kind), Length, KmerSize, strandBiasScore, TmrAltQual);

  if (!STRResult.empty()) info += absl::StrFormat(";MS=%s", STRResult);

  if (params.tenxMode) {
    const auto nmlHpScore = PhredFisherScore(NormalCov.RefHP(Haplotype::FIRST), NormalCov.RefHP(Haplotype::SECOND),
                                             NormalCov.AltHP(Haplotype::FIRST), NormalCov.AltHP(Haplotype::SECOND));

    const auto tmrHpScore = PhredFisherScore(TumorCov.RefHP(Haplotype::FIRST), TumorCov.RefHP(Haplotype::SECOND),
                                             TumorCov.AltHP(Haplotype::FIRST), TumorCov.AltHP(Haplotype::SECOND));

    const auto pairHpScore = PhredFisherScore(NormalCov.TotalHP(Haplotype::FIRST), NormalCov.TotalHP(Haplotype::SECOND),
                                              TumorCov.TotalHP(Haplotype::FIRST), TumorCov.TotalHP(Haplotype::SECOND));

    info += absl::StrFormat(";HPS=%f;HPSN=%f;HPST=%f", pairHpScore, nmlHpScore, tmrHpScore);
  }

  const auto lowSomaticScore =
      !STRResult.empty() ? somaticScore < params.minSTRFisher : somaticScore < params.minFisher;

  std::vector<std::string> filters{};
  if (lowSomaticScore && !STRResult.empty()) filters.emplace_back("LowFisherSTR");
  if (lowSomaticScore && STRResult.empty()) filters.emplace_back("LowFisherScore");
  if (NormalCov.TotalCov() < params.minNmlCov) filters.emplace_back("LowCovNormal");
  if (NormalCov.TotalCov() > params.maxNmlCov) filters.emplace_back("HighCovNormal");
  if (TumorCov.TotalCov() < params.minTmrCov) filters.emplace_back("LowCovTumor");
  if (TumorCov.TotalCov() > params.maxTmrCov) filters.emplace_back("HighCovTumor");
  if (TmrAltQual < static_cast<float>(params.minBaseQual)) filters.emplace_back("LowTmrAltQual");

  const auto tmrVaf = TumorCov.VAF();
  const auto nmlVaf = NormalCov.VAF();

  if (tmrVaf < params.minTmrVAF) filters.emplace_back("LowVafTumor");
  if (nmlVaf > params.maxNmlVAF) filters.emplace_back("HighVafNormal");
  if (TumorCov.TotalAltCov() < params.minTmrAltCnt) filters.emplace_back("LowAltCntTumor");
  if (NormalCov.TotalAltCov() > params.maxNmlAltCnt) filters.emplace_back("HighAltCntNormal");

  if (TumorCov.AltAllele.FwdCov < params.minStrandCnt || TumorCov.AltAllele.RevCov < params.minStrandCnt) {
    filters.emplace_back("StrandBias");
  }

  if (params.tenxMode && varState == VariantState::SOMATIC && TumorCov.AltHP(Haplotype::FIRST) > 0 &&
      TumorCov.AltHP(Haplotype::SECOND) > 0) {
    filters.emplace_back("MultiHP");
  }

  const std::string filter = filters.empty() ? "PASS" : absl::StrJoin(filters, ";");
  const auto* format = params.tenxMode ? "GT:AD:SR:SA:DP:HPR:HPA" : "GT:AD:SR:SA:DP";

  return absl::StrFormat("%s\t%d\t.\t%s\t%s\t%f\t%s\t%s\t%s\t%s\t%s\n", ChromName, Position, RefAllele, AltAllele,
                         somaticScore, filter, info, format, BuildSampleFormat(NormalCov, params.tenxMode),
                         BuildSampleFormat(TumorCov, params.tenxMode));
}

auto Variant::ID() const -> VariantID {
  const auto state = absl::StrFormat("%s|%d|%s|%s", ChromName, Position, RefAllele, AltAllele);
  return absl::hash_internal::CityHash64WithSeeds(state.c_str(), state.length(), utils::PRIME_0,  // NOLINT
                                                  utils::PRIME_1);
}

auto Variant::ComputeState() const -> VariantState {
  const auto nmlAlt = NormalCov.TotalAltCov();
  const auto tmrAlt = TumorCov.TotalAltCov();

  if (nmlAlt > 0 && tmrAlt > 0) return VariantState::SHARED;
  if (nmlAlt == 0 && tmrAlt > 0) return VariantState::SOMATIC;
  if (nmlAlt > 0 && tmrAlt == 0) return VariantState::NORMAL;
  return VariantState::NONE;
}

auto Variant::Genotype(int ref, int alt) -> std::string {
  if (ref > 0 && alt == 0) return "0/0";
  if (ref > 0 && alt > 0) return "0/1";
  if (ref == 0 && alt > 0) return "1/1";
  return "./.";
}

auto Variant::BuildSampleFormat(const VariantHpCov& v, bool is_tenx_mode) -> std::string {
  //  const auto FORMAT = params.tenxModeOn ? "GT:AD:SR:SA:DP:HPR:HPA" : "GT:AD:SR:SA:DP";
  auto result = absl::StrFormat("%s:%d,%d:%d,%d:%d,%d:%d", Genotype(v.TotalRefCov(), v.TotalAltCov()), v.TotalRefCov(),
                                v.TotalAltCov(), v.RefAllele.FwdCov, v.RefAllele.RevCov, v.AltAllele.FwdCov,
                                v.AltAllele.RevCov, v.TotalCov());

  if (is_tenx_mode) {
    result += absl::StrFormat(":%d,%d,%d:%d,%d,%d", v.RefHP(Haplotype::FIRST), v.RefHP(Haplotype::SECOND),
                              v.RefHP(Haplotype::UNASSIGNED), v.AltHP(Haplotype::FIRST), v.AltHP(Haplotype::SECOND),
                              v.AltHP(Haplotype::UNASSIGNED));
  }

  return result;
}
}  // namespace lancet2
