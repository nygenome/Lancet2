#include "lancet/cli/vcf_header_builder.h"

#include "lancet/base/types.h"
#include "lancet/base/version.h"
#include "lancet/core/sample_header_reader.h"
#include "lancet/hts/reference.h"

#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "absl/time/clock.h"
#include "absl/time/time.h"
#include "spdlog/fmt/bundled/base.h"
#include "spdlog/fmt/bundled/format.h"

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace {

using namespace std::string_view_literals;

// clang-format off

// ============================================================================
// VCF v4.5 header format template — named placeholders filled by fmt::format.
// Kept as a constexpr raw string to avoid line-by-line StrAppend noise.
// ============================================================================
constexpr auto FORMAT_STR_HEADER = R"raw(##fileformat=VCFv4.5
##fileDate={RUN_TIMESTAMP}
##source=Lancet_{FULL_VERSION_TAG}
##commandLine="{FULL_COMMAND_USED}"
##reference="{REFERENCE_PATH}"
{CONTIG_HDR_LINES}{CONDITIONAL_INFO_LINES}##INFO=<ID=TYPE,Number=A,Type=String,Description="Variant type (SNV, INS, DEL, MNP)">
##INFO=<ID=LENGTH,Number=A,Type=Integer,Description="Variant length in base pairs">
##INFO=<ID=MULTIALLELIC,Number=0,Type=Flag,Description="Indicates if the site has multiple ALT alleles">
{ANNOTATION_INFO_LINES}##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Forward strand allele depth">
##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Reverse strand allele depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=RMQ,Number=R,Type=Float,Description="RMS mapping quality per allele">
##FORMAT=<ID=NPBQ,Number=R,Type=Float,Description="Normalized posterior base quality per allele (raw PBQ / allele depth)">
##FORMAT=<ID=SB,Number=1,Type=Float,Description="Strand bias log odds ratio (Haldane-corrected, coverage-invariant)">
##FORMAT=<ID=SCA,Number=1,Type=Float,Description="Soft clip asymmetry (ALT minus REF)">
##FORMAT=<ID=FLD,Number=1,Type=Float,Description="Fragment length delta (signed mean ALT isize minus mean REF isize)">
##FORMAT=<ID=RPCD,Number=1,Type=Float,Description="Read position Cohen's D effect size (. if untestable)">
##FORMAT=<ID=BQCD,Number=1,Type=Float,Description="Base quality Cohen's D effect size (. if untestable)">
##FORMAT=<ID=MQCD,Number=1,Type=Float,Description="Mapping quality Cohen's D effect size (. if untestable)">
##FORMAT=<ID=ASMD,Number=1,Type=Float,Description="Allele-specific mismatch delta (mean ALT NM minus mean REF NM minus variant length)">
##FORMAT=<ID=SDFC,Number=1,Type=Float,Description="Site depth fold change (sample DP / per-sample window mean coverage)">
##FORMAT=<ID=PRAD,Number=1,Type=Float,Description="Polar radius: log10(1 + sqrt(AD_Ref^2 + AD_Alt^2))">
##FORMAT=<ID=PANG,Number=1,Type=Float,Description="Polar angle: allele identity ratio atan2(AD_Alt, AD_Ref) in radians">
##FORMAT=<ID=CMLOD,Number=A,Type=Float,Description="Continuous mixture log-odds score per ALT allele (base-quality-weighted LOD vs null)">
##FORMAT=<ID=FSSE,Number=1,Type=Float,Description="Fragment start Shannon entropy [0,1]. Measures spatial diversity of ALT read start positions. Catches PCR jackpot duplicates that survive MarkDuplicates via exonuclease fraying, alignment jitter, and representative read roulette. Missing (.) if fewer than 3 ALT-supporting reads.">
##FORMAT=<ID=AHDD,Number=1,Type=Float,Description="ALT-haplotype discordance delta. Mean ALT NM against own haplotype minus mean REF NM against REF. High values signal assembly hallucinations. Missing (.) if either group is empty.">
##FORMAT=<ID=HSE,Number=1,Type=Float,Description="Haplotype segregation entropy [0,1]. How concentrated ALT reads are on a single SPOA path. Near 0 = concentrated (true variant). Near 1 = scattered (noise). Missing (.) if fewer than 3 ALT reads or single haplotype.">
##FORMAT=<ID=PDCV,Number=1,Type=Float,Description="Path depth coefficient of variation. K-mer coverage uniformity along the ALT de Bruijn graph path. High values signal chimeric junctions with uneven support. Missing (.) if path has fewer than 2 nodes.">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods (Dirichlet-Multinomial model)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality (second-lowest PL from the Dirichlet-Multinomial model, capped at 99)">
)raw"sv;

constexpr auto CASE_CTRL_INFO_HDR_LINES = R"raw(##INFO=<ID=SHARED,Number=0,Type=Flag,Description="Variant ALT seen in both case & control sample(s)">
##INFO=<ID=CTRL,Number=0,Type=Flag,Description="Variant ALT seen only in control sample(s)">
##INFO=<ID=CASE,Number=0,Type=Flag,Description="Variant ALT seen only in case sample(s)">
)raw"sv;

constexpr auto ANNOTATION_INFO_HDR_LINES = R"raw(##INFO=<ID=GRAPH_CX,Number=3,Type=String,Description="Graph complexity metrics: GEI,TipToPathCovRatio,MaxSingleDirDegree">
##INFO=<ID=SEQ_CX,Number=11,Type=String,Description="Sequence complexity features: ContextHRun,ContextEntropy,ContextFlankLQ,ContextHaplotypeLQ,DeltaHRun,DeltaEntropy,DeltaFlankLQ,TrAffinity,TrPurity,TrPeriod,IsStutterIndel">
)raw"sv;

constexpr auto CONTIG_HDR_LINE_FORMAT = "##contig=<ID={},length={}>\n"sv;

// clang-format on

}  // namespace

namespace lancet::cli {

// ============================================================================
// BuildVcfHeader — assembles the complete VCF v4.5 header string
//
// Opens the reference FASTA to enumerate contigs.  Conditionally emits
// SHARED/CTRL/CASE INFO lines when case-control mode is active.
// Zero coupling to PipelineRunner instance state.
// ============================================================================
auto BuildVcfHeader(CliParams const& params) -> std::string {
  // --- Contig lines ---
  std::string contig_hdr_lines;
  static constexpr usize CONTIGS_BUFFER_SIZE = 524'288;
  contig_hdr_lines.reserve(CONTIGS_BUFFER_SIZE);
  hts::Reference const ref(params.mVariantBuilder.mRdCollParams.mRefPath);
  for (auto const& chrom : ref.ListChroms()) {
    absl::StrAppend(&contig_hdr_lines,
                    fmt::format(CONTIG_HDR_LINE_FORMAT, chrom.Name(), chrom.Length()));
  }

  // --- SHARED/CTRL/CASE INFO headers — only when case-control mode is active ---
  std::string conditional_info_lines;
  if (params.mIsCaseCtrlMode) {
    absl::StrAppend(&conditional_info_lines, CASE_CTRL_INFO_HDR_LINES);
  }

  // --- Complexity feature INFO headers — always emitted ---
  std::string annotation_info_lines;
  absl::StrAppend(&annotation_info_lines, ANNOTATION_INFO_HDR_LINES);

  auto full_hdr = fmt::format(
      // FORMAT_STR_HEADER is a constexpr `char const[N]`; fmt::format takes a string view and the
      // array decays to a pointer at the call site. The decay is required and unavoidable.
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      FORMAT_STR_HEADER,
      fmt::arg("RUN_TIMESTAMP",
               absl::FormatTime(absl::RFC3339_sec, absl::Now(), absl::LocalTimeZone())),
      fmt::arg("FULL_VERSION_TAG", lancet::base::LancetFullVersion()),
      fmt::arg("FULL_COMMAND_USED", params.mFullCmdLine),
      fmt::arg("REFERENCE_PATH", params.mVariantBuilder.mRdCollParams.mRefPath.string()),
      fmt::arg("CONTIG_HDR_LINES", contig_hdr_lines),
      fmt::arg("CONDITIONAL_INFO_LINES", conditional_info_lines),
      fmt::arg("ANNOTATION_INFO_LINES", annotation_info_lines));

  auto const rc_sample_list = core::BuildSampleNameList(params.mVariantBuilder.mRdCollParams);
  absl::StrAppend(&full_hdr,
                  fmt::format("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
                              absl::StrJoin(rc_sample_list, "\t")));

  return full_hdr;
}

}  // namespace lancet::cli
