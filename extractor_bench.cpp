#include <algorithm>
#include <type_traits>

#include "benchmark/benchmark.h"
#include "lancet/hts/extractor.h"
#include "lancet_benchmark_config.h"

namespace {

void ExtractorCramCoreQname(benchmark::State& state) {
  using namespace lancet::hts;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    const Reference ref(Hg38Reference);
    Extractor extractor(TumorCram, ref, Alignment::Fields::CORE_QNAME);
    extractor.SetNumThreads(static_cast<int>(state.threads()));
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
  }
}

void ExtractorCramCigarSeqQual(benchmark::State& state) {
  using namespace lancet::hts;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    const Reference ref(Hg38Reference);
    Extractor extractor(TumorCram, ref, Alignment::Fields::CIGAR_SEQ_QUAL);
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
  }
}

void ExtractorCramAuxRgaux(benchmark::State& state) {
  using namespace lancet::hts;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    const Reference ref(Hg38Reference);
    Extractor extractor(TumorCram, ref, Alignment::Fields::AUX_RGAUX, {"RG", "MC", "NM", "SA", "XS", "MD", "AS"});
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
  }
}

void ExtractorBamCoreQname(benchmark::State& state) {
  using namespace lancet::hts;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    const Reference ref(Hg38Reference);
    Extractor extractor(TumorBam, ref, Alignment::Fields::CORE_QNAME);
    extractor.SetNumThreads(static_cast<int>(state.threads()));
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
    ;
  }
}

void ExtractorBamCigarSeqQual(benchmark::State& state) {
  using namespace lancet::hts;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    const Reference ref(Hg38Reference);
    Extractor extractor(TumorBam, ref, Alignment::Fields::CIGAR_SEQ_QUAL);
    extractor.SetNumThreads(static_cast<int>(state.threads()));
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
  }
}

void ExtractorBamAuxRgaux(benchmark::State& state) {
  using namespace lancet::hts;
  // NOLINTNEXTLINE(readability-identifier-length)
  for ([[maybe_unused]] auto _ : state) {
    const Reference ref(Hg38Reference);
    Extractor extractor(TumorBam, ref, Alignment::Fields::AUX_RGAUX, {"RG", "MC", "NM", "SA", "XS", "MD", "AS"});
    extractor.SetNumThreads(static_cast<int>(state.threads()));
    std::vector<Alignment> results(extractor.begin(), extractor.end());
    benchmark::DoNotOptimize(results);
  }
}

}  // namespace

// NOLINTBEGIN
BENCHMARK(ExtractorCramCoreQname)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);
BENCHMARK(ExtractorCramCigarSeqQual)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);
BENCHMARK(ExtractorCramAuxRgaux)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);

BENCHMARK(ExtractorBamCoreQname)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);
BENCHMARK(ExtractorBamCigarSeqQual)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);
BENCHMARK(ExtractorBamAuxRgaux)->Unit(benchmark::kMillisecond)->DenseThreadRange(1, 8);
// NOLINTEND
