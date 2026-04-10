#include "lancet/hts/extractor.h"

#include "lancet/hts/alignment.h"
#include "lancet/hts/reference.h"

#include "benchmark/benchmark.h"
#include "lancet_benchmark_config.h"

#include <vector>

namespace {

void ExtractorCramCoreQname(benchmark::State& state) {
  using lancet::hts::Alignment;
  using lancet::hts::Extractor;
  using lancet::hts::Reference;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    Reference const ref(Hg38Reference);
    Extractor extractor(TumorCram, ref, Alignment::Fields::CORE_QNAME);
    extractor.SetNumThreads(state.threads());
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
  }
}

void ExtractorCramCigarSeqQual(benchmark::State& state) {
  using lancet::hts::Alignment;
  using lancet::hts::Extractor;
  using lancet::hts::Reference;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    Reference const ref(Hg38Reference);
    Extractor extractor(TumorCram, ref, Alignment::Fields::CIGAR_SEQ_QUAL);
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
  }
}

void ExtractorCramAuxRgaux(benchmark::State& state) {
  using lancet::hts::Alignment;
  using lancet::hts::Extractor;
  using lancet::hts::Reference;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    Reference const ref(Hg38Reference);
    Extractor extractor(TumorCram, ref, Alignment::Fields::AUX_RGAUX,
                        {"RG", "MC", "NM", "SA", "XS", "MD", "AS"});
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
  }
}

void ExtractorBamCoreQname(benchmark::State& state) {
  using lancet::hts::Alignment;
  using lancet::hts::Extractor;
  using lancet::hts::Reference;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    Reference const ref(Hg38Reference);
    Extractor extractor(TumorBam, ref, Alignment::Fields::CORE_QNAME);
    extractor.SetNumThreads(state.threads());
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
    ;
  }
}

void ExtractorBamCigarSeqQual(benchmark::State& state) {
  using lancet::hts::Alignment;
  using lancet::hts::Extractor;
  using lancet::hts::Reference;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    Reference const ref(Hg38Reference);
    Extractor extractor(TumorBam, ref, Alignment::Fields::CIGAR_SEQ_QUAL);
    extractor.SetNumThreads(state.threads());
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
  }
}

void ExtractorBamAuxRgaux(benchmark::State& state) {
  using lancet::hts::Alignment;
  using lancet::hts::Extractor;
  using lancet::hts::Reference;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    Reference const ref(Hg38Reference);
    Extractor extractor(TumorBam, ref, Alignment::Fields::AUX_RGAUX,
                        {"RG", "MC", "NM", "SA", "XS", "MD", "AS"});
    extractor.SetNumThreads(state.threads());
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
  }
}

}  // namespace

// NOLINTBEGIN(cert-err58-cpp, cppcoreguidelines-owning-memory, readability-identifier-length, misc-use-anonymous-namespace)
BENCHMARK(ExtractorCramCoreQname)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);
BENCHMARK(ExtractorCramCigarSeqQual)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);
BENCHMARK(ExtractorCramAuxRgaux)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);

BENCHMARK(ExtractorBamCoreQname)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);
BENCHMARK(ExtractorBamCigarSeqQual)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);
BENCHMARK(ExtractorBamAuxRgaux)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);
// NOLINTEND(cert-err58-cpp, cppcoreguidelines-owning-memory, readability-identifier-length, misc-use-anonymous-namespace)
