#include "lancet2/variant_store.h"

#include <algorithm>
#include <filesystem>
#include <mutex>
#include <utility>

#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "absl/time/clock.h"
#include "absl/time/time.h"  // NOLINT
#include "generated/lancet2_version.h"
#include "lancet2/fasta_reader.h"

namespace lancet2 {
VariantStore::VariantStore(std::shared_ptr<const CliParams> p) : params(std::move(p)) {}

static inline auto GetContigHeaderLines(const std::filesystem::path& faPath) -> std::string {
  FastaReader refRdr(faPath);
  auto contigs = refRdr.GetContigs();
  const auto ctgIds = refRdr.GetContigIndexMap();
  std::sort(contigs.begin(), contigs.end(), [&ctgIds](const ContigInfo& c1, const ContigInfo& c2) -> bool {
    return ctgIds.at(c1.contigName) < ctgIds.at(c2.contigName);
  });

  std::vector<std::string> contigLines;
  contigLines.reserve(contigs.size());
  for (const auto& ctgInfo : contigs) {
    contigLines.emplace_back(absl::StrFormat("##contig=<ID=%s,length=%d>", ctgInfo.contigName, ctgInfo.contigLen));
  }
  return absl::StrJoin(contigLines, "\n");
}

auto VariantStore::GetHeader(const std::vector<std::string>& sample_names, const CliParams& p) -> std::string {
  // clang-format off
  const auto stdResult = absl::StrFormat(R"raw(##fileformat=VCFv4.3
##fileDate=%s
##source=lancet2_%s
##commandLine="%s"
##reference="%s"
%s
##FILTER=<ID=LowFisherSTR,Description="Fisher exact test score for tumor/normal STR allele counts less than %f">
##FILTER=<ID=LowFisherScore,Description="Fisher exact test score for tumor/normal allele counts less than %f">
##FILTER=<ID=LowCovNormal,Description="Allele coverage in normal less than %d">
##FILTER=<ID=HighCovNormal,Description="Allele coverage in normal greater than %d">
##FILTER=<ID=LowCovTumor,Description="Allele coverage in tumor less than %d">
##FILTER=<ID=HighCovTumor,Description="Allele coverage in tumor greater than %d}">
##FILTER=<ID=LowTmrRefQual,Description="Average base qualities for the reference allele in tumor less than %f">
##FILTER=<ID=LowTmrAltQual,Description="Average base qualities for the alternate allele in tumor less than %f">
##FILTER=<ID=LowNmlRefQual,Description="Average base qualities for the reference allele in normal less than %f">
##FILTER=<ID=LowNmlAltQual,Description="Average base qualities for the alternate allele in normal less than %f">
##FILTER=<ID=LowVafTumor,Description="Variant allele frequency in tumor less than %f">
##FILTER=<ID=HighVafNormal,Description="Variant allele frequency in normal greater than %f">
##FILTER=<ID=LowAltCntTumor,Description="Alternative allele count in tumor less than %d">
##FILTER=<ID=HighAltCntNormal,Description="Alternative allele count in normal greater than %d">
##FILTER=<ID=StrandBias,Description="Number of non-reference reads in either forward or reverse strand below threshold %d">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Mutation present only in tumor">
##INFO=<ID=NORMAL,Number=0,Type=Flag,Description="Mutation present only in normal">
##INFO=<ID=SHARED,Number=0,Type=Flag,Description="Mutation present in both tumor and normal">
##INFO=<ID=NONE,Number=0,Type=Flag,Description="Mutation not supported by data">
##INFO=<ID=FETS,Number=1,Type=Float,Description="Phred-scaled probability of Fisher exact test of ref and alt allele counts in tumor and normal">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type. Possible values are snv, ins, del and complex">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="Variant length in base pairs">
##INFO=<ID=KMERSIZE,Number=1,Type=Integer,Description="K-mer length used to assemble the locus">
##INFO=<ID=SB,Number=1,Type=Float,Description="Phred-scaled probability of Fisher exact test of fwd and rev tumor read counts in ref and alt alleles. i.e. Strand Bias">
##INFO=<ID=TMR_REF_QUAL,Number=1,Type=Float,Description="Average base qualities for the reference allele in tumor sample">
##INFO=<ID=TMR_ALT_QUAL,Number=1,Type=Float,Description="Average base qualities for the alternate allele in tumor sample">
##INFO=<ID=NML_REF_QUAL,Number=1,Type=Float,Description="Average base qualities for the reference allele in normal sample">
##INFO=<ID=NML_ALT_QUAL,Number=1,Type=Float,Description="Average base qualities for the alternate allele in normal sample">
##INFO=<ID=MS,Number=1,Type=String,Description="Description of the microsatellite length and motif, if found (format: LENGTH:MOTIF)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Number of reads supporting ref and alt alleles at the site">
##FORMAT=<ID=SR,Number=2,Type=Integer,Description="Number of reads in the forward and reverse strands supporting the reference allele">
##FORMAT=<ID=SA,Number=2,Type=Integer,Description="Number of reads in the forward and reverse strands supporting the alternate allele">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
)raw", absl::FormatTime(absl::RFC3339_sec, absl::Now(), absl::LocalTimeZone()), LONG_VERSION, // NOLINT
       p.commandLine, p.referencePath, GetContigHeaderLines(p.referencePath), p.minSTRFisher, p.minFisher,
       p.minNmlCov, p.maxNmlCov, p.minTmrCov, p.maxTmrCov, p.minBaseQual, p.minBaseQual, p.minBaseQual, p.minBaseQual,
       p.minTmrVAF, p.maxNmlVAF, p.minTmrAltCnt, p.maxNmlAltCnt, p.minStrandCnt);

  static constexpr auto tenxTemplate = R"raw(##FILTER=<ID=MultiHP,Description="Reads supporting alternate allele found in multiple haplotypes">
##INFO=<ID=HPS,Number=1,Type=Float,Description="Phred-scaled probability of Fisher exact test of total haplotype counts (HP1 + HP2) in normal and tumor">
##INFO=<ID=HPSN,Number=1,Type=Float,Description="Phred-scaled probability of haplotype counts (HP1, HP2) in ref and alt alleles for normal">
##INFO=<ID=HPST,Number=1,Type=Float,Description="Phred-scaled probability of haplotype counts (HP1, HP2) in ref and alt alleles for tumor">
##FORMAT=<ID=HPR,Number=3,Type=Integer,Description="Number of reads supporting reference allele in haplotype 1, 2 and 0 respectively (0 = Reads with unassigned haplotype)">
##FORMAT=<ID=HPA,Number=3,Type=Integer,Description="Number of reads supporting alternate allele in haplotype 1, 2 and 0 respectively (0 = Reads with unassigned haplotype)">)raw";
  // clang-format on

  const auto chromLine =
      absl::StrFormat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", absl::StrJoin(sample_names, "\t"));

  return p.tenxMode ? absl::StrFormat("%s%s%s", stdResult, tenxTemplate, chromLine)
                    : absl::StrFormat("%s%s", stdResult, chromLine);
}

auto VariantStore::TryAddVariants(absl::Span<const Variant> variants) -> bool {
  const auto isLocked = spinLock.TryLock();
  if (!isLocked) return false;

  UnsafeAddVariantBatch(variants);
  spinLock.Unlock();
  return true;
}

void VariantStore::ForceAddVariants(absl::Span<const Variant> variants) {
  std::lock_guard<utils::SpinLock> guard(spinLock);
  UnsafeAddVariantBatch(variants);
}

auto VariantStore::FlushWindow(const RefWindow& w, std::ostream& outStream, const ContigIDs& ctg_ids) -> bool {
  std::lock_guard<utils::SpinLock> guard(spinLock);

  using vDBPair = std::pair<u64, Variant>;
  std::vector<VariantID> variantIDsToFlush;
  variantIDsToFlush.reserve(1024);
  std::for_each(data.cbegin(), data.cend(), [&variantIDsToFlush, &w, &ctg_ids](const vDBPair& p) {
    if (IsVariantInOrBefore(p.second, w, ctg_ids)) variantIDsToFlush.emplace_back(p.first);
  });

  return Flush(absl::MakeConstSpan(variantIDsToFlush), outStream, ctg_ids);
}

auto VariantStore::FlushAll(std::ostream& outStream, const ContigIDs& ctg_ids) -> bool {
  std::lock_guard<utils::SpinLock> guard(spinLock);
  if (data.empty()) return false;

  std::vector<u64> variantIDsToFlush;
  variantIDsToFlush.reserve(data.size());
  for (const auto& p : data) variantIDsToFlush.emplace_back(p.first);
  return Flush(absl::MakeConstSpan(variantIDsToFlush), outStream, ctg_ids);
}

auto VariantStore::Flush(absl::Span<const VariantID> ids, std::ostream& outStream, const ContigIDs& ctg_ids) -> bool {
  std::vector<Variant> variantsToFlush;
  variantsToFlush.reserve(ids.size());

  for (const auto variantId : ids) {
    const auto handle = data.extract(variantId);
    if (handle.empty()) continue;
    variantsToFlush.emplace_back(handle.mapped());
  }

  std::sort(variantsToFlush.begin(), variantsToFlush.end(),
            [&ctg_ids](const Variant& v1, const Variant& v2) -> bool { return IsVariant1LessThan2(v1, v2, ctg_ids); });

  std::for_each(variantsToFlush.cbegin(), variantsToFlush.cend(),
                [&outStream, this](const Variant& v) { outStream << v.MakeVcfLine(*params); });

  if (!variantsToFlush.empty() && data.load_factor() < 0.8) {
    // swap data with temp to force releasing hash table memory
    absl::flat_hash_map<VariantID, Variant> temp(data);
    data.swap(temp);
  }

  return !variantsToFlush.empty();
}

auto VariantStore::IsVariant1LessThan2(const Variant& v1, const Variant& v2,
                                       const absl::flat_hash_map<std::string, i64>& ctg_ids) -> bool {
  if (v1.ChromName != v2.ChromName) return ctg_ids.at(v1.ChromName) < ctg_ids.at(v2.ChromName);
  if (v1.Position != v2.Position) return v1.Position < v2.Position;
  if (v1.RefAllele != v2.RefAllele) return v1.RefAllele < v2.RefAllele;
  return v1.AltAllele < v2.AltAllele;
}

auto VariantStore::IsVariantInOrBefore(const Variant& v, const RefWindow& w, const ContigIDs& ctg_ids) -> bool {
  if (v.ChromName != w.GetChromName()) return ctg_ids.at(v.ChromName) < ctg_ids.at(w.GetChromName());
  return v.Position <= (static_cast<usize>(w.GetEndPos0()) + 1);
}

void VariantStore::UnsafeAddVariantBatch(absl::Span<const Variant> variants) {
  for (const auto& variant : variants) {
    const auto vID = variant.ID();
    auto itr = data.find(vID);
    if (itr == data.end()) {
      data.emplace(vID, variant);
      continue;
    }

    const auto oldTotalCov = itr->second.TumorCov.TotalCov() + itr->second.NormalCov.TotalCov();
    const auto newTotalCov = variant.TumorCov.TotalCov() + variant.NormalCov.TotalCov();
    if (oldTotalCov < newTotalCov) {
      itr->second.KmerSize = variant.KmerSize;
      itr->second.TumorCov = variant.TumorCov;
      itr->second.NormalCov = variant.NormalCov;
    }
  }
}
}  // namespace lancet2
