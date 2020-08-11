#include "lancet/cli.h"

#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <thread>

#include "CLI/CLI.hpp"
#include "absl/debugging/failure_signal_handler.h"
#include "absl/debugging/symbolize.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "generated/lancet_version.h"
#include "lancet/cli_params.h"
#include "lancet/run_pipeline.h"
#include "spdlog/sinks/stdout_color_sinks-inl.h"
#include "spdlog/spdlog.h"

#ifdef SPDLOG_ACTIVE_LEVEL
#undef SPDLOG_ACTIVE_LEVEL
#endif

#ifndef NDEBUG
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE  // NOLINT
#else
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO  // NOLINT
#endif

namespace lancet {
auto PipelineSubcmd(CLI::App* app, std::shared_ptr<CliParams> params) -> void;

auto RunCli(int argc, char** argv) noexcept -> int {
  absl::InitializeSymbolizer(argv[0]);  // NOLINT
  absl::FailureSignalHandlerOptions fshOpts{};
  absl::InstallFailureSignalHandler(fshOpts);

  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  auto stderrLogger = spdlog::stderr_color_mt("stderr", spdlog::color_mode::automatic);
  spdlog::set_default_logger(stderrLogger);
  spdlog::set_pattern("%^%Y-%m-%dT%H:%M:%S%z | [%L] | %v%$");

  constexpr auto appDesc = "Microassembly based somatic variant caller";
  CLI::App app(appDesc, "Lancet");
  app.require_subcommand(1);

  const auto pipelineParams = std::make_shared<CliParams>();
  PipelineSubcmd(&app, pipelineParams);

  static const auto printVersion = [](std::size_t count) -> void {
    if (count <= 0) return;
    std::cout << absl::StreamFormat("Lancet %s\n", lancet::LONG_VERSION);
    std::exit(EXIT_SUCCESS);
  };

  static const auto printHelp = [&app](std::size_t count) -> void {
    if (count <= 0) return;
    std::cerr << app.help(app.get_name(), CLI::AppFormatMode::Normal);
    std::exit(EXIT_SUCCESS);
  };

  pipelineParams->commandLine = std::string(argv[0]);  // NOLINT
  for (auto idx = 1; idx < argc; idx++) {
    absl::StrAppend(&pipelineParams->commandLine, " ", argv[idx]);  // NOLINT
  }

  app.set_help_flag();
  app.failure_message(CLI::FailureMessage::help);
  constexpr std::size_t HELP_FORMATTER_WIDTH = 65;
  app.get_formatter()->column_width(HELP_FORMATTER_WIDTH);
  app.add_flag_function("-v,--version", printVersion, "Print version information")->group("Flags");
  app.add_flag_function("-h,--help", printHelp, "Print this help message and exit")->group("Flags");

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& err) {
    return app.exit(err);
  } catch (std::exception& err) {
    return app.exit(CLI::ParseError(err.what(), CLI::ExitCodes::ValidationError));
  }

  return EXIT_SUCCESS;
}

auto PipelineSubcmd(CLI::App* app, std::shared_ptr<CliParams> params) -> void {  // NOLINT
  auto* subcmd = app->add_subcommand("pipeline", "Run Lancet variant calling pipeline");

  // Required
  subcmd->add_option("-t,--tumor", params->tumorPath, "Path to tumor BAM/CRAM file")
      ->required(true)
      ->group("Required")
      ->check(CLI::ExistingFile);

  subcmd->add_option("-n,--normal", params->normalPath, "Path to normal BAM/CRAM file")
      ->required(true)
      ->group("Required")
      ->check(CLI::ExistingFile);

  subcmd->add_option("-r,--reference", params->referencePath, "Path to reference FASTA file")
      ->required(true)
      ->group("Required")
      ->check(CLI::ExistingFile);

  subcmd->add_option("-o,--out-prefix", params->outVcfPrefix, "Output prefix to compressed VCF - ${prefix}.vcf.gz")
      ->required(true)
      ->group("Required");

  // Regions
  subcmd->add_option("--region", params->inRegions, "One or more regions to process (samtools region format)")
      ->group("Regions")
      ->type_name("REF:[:START[-END]]");

  subcmd->add_option("-b,--bed-file", params->bedFilePath, "Path to BED file with regions to process")
      ->group("Regions")
      ->check(CLI::ExistingFile);

  subcmd->add_option("-P,--padding", params->regionPadLength, "Padding to apply for all input regions", true)
      ->group("Regions");

  subcmd->add_option("-w,--window-size", params->windowLength, "Window size for each microassemby task", true)
      ->group("Regions");

  subcmd->add_option("--pct-overlap", params->pctOverlap, "Percent overlap between consecutive windows", true)
      ->group("Regions")
      ->check(CLI::Range(std::uint32_t(5), std::uint32_t(95)));

  // Parameters
  const auto maxNumThreads = static_cast<std::uint32_t>(std::thread::hardware_concurrency());
  subcmd->add_option("-T,--num-threads", params->numWorkerThreads, "Number of additional worker threads", true)
      ->group("Parameters")
      ->check(CLI::Range(std::uint32_t(1), maxNumThreads));

  subcmd->add_option("-k,--min-kmer-length", params->minKmerSize, "Min. kmer length for graph nodes", true)
      ->group("Parameters")
      ->check(CLI::Range(std::uint32_t(11), std::uint32_t(99)));

  subcmd->add_option("-K,--max-kmer-length", params->maxKmerSize, "Max. kmer length for graph nodes", true)
      ->group("Parameters")
      ->check(CLI::Range(std::uint32_t(13), std::uint32_t(101)));

  subcmd->add_option("--min-trim-qual", params->trimBelowQual, "Min. base quality to trim 5' and 3' read bases", true)
      ->group("Parameters")
      ->check(CLI::Range(std::uint32_t(0), std::uint32_t(30)));

  subcmd->add_option("-q,--min-base-qual", params->minBaseQual, "Min. base quality to consider for SNV calling", true)
      ->group("Parameters")
      ->check(CLI::Range(std::uint32_t(0), std::uint32_t(30)));

  subcmd->add_option("-Q,--min-mapping-qual", params->minReadMappingQual, "Min. mapping quality to use a read", true)
      ->group("Parameters");

  subcmd->add_option("--max-rpt-mismatch", params->maxRptMismatch, "Max. mismatches to detect approx. repeats", true)
      ->group("Parameters");

  subcmd->add_option("--max-tip-length", params->minGraphTipLength, "Max. allowed tip length in the graph", true)
      ->group("Parameters");

  subcmd->add_option("--graph-traversal-limit", params->graphTraversalLimit, "Max. BFS/DFS graph traversal limit", true)
      ->group("Parameters");

  subcmd->add_option("--max-indel-length", params->maxIndelLength, "Max. limit on the indel length to detect", true)
      ->group("Parameters");

  subcmd->add_option("--min-anchor-cov", params->minAnchorCov, "Min. coverage for anchor nodes (source & sink)", true)
      ->group("Parameters");

  subcmd->add_option("--min-node-cov", params->minNodeCov, "Min. coverage for all nodes in the graph", true)
      ->group("Parameters");

  subcmd->add_option("--min-cov-ratio", params->minCovRatio, "Min. node by window coverage for all nodes", true)
      ->group("Parameters");

  subcmd->add_option("--max-window-cov", params->maxWindowCov, "Max. average window coverage before downsampling", true)
      ->group("Parameters");

  subcmd->add_option("--min-as-xs-diff", params->minReadAsXsDiff, "Min. diff. between AS & XS scores (BWA-mem)", true)
      ->group("Parameters");

  //  STR parameters
  subcmd->add_option("--max-str-unit-len", params->maxSTRUnitLength, "Max. unit length for an STR motif", true)
      ->group("STR parameters");

  subcmd->add_option("--min-str-units", params->minSTRUnits, "Min. number of STR units required to report", true)
      ->group("STR parameters");

  subcmd->add_option("--min-str-len", params->minSTRLen, "Min. required STR length (in bp) to report", true)
      ->group("STR parameters");

  subcmd->add_option("--max-str-dist", params->maxSTRDist, "Max. distance (in bp) of variant from the STR motif", true)
      ->group("STR parameters");

  // Filters
  subcmd->add_option("-c,--max-nml-alt-cnt", params->maxNmlAltCnt, "Max. ALT allele count in normal sample", true)
      ->group("Filters");

  subcmd->add_option("-C,--min-tmr-alt-cnt", params->minTmrAltCnt, "Min. ALT allele count in tumor sample", true)
      ->group("Filters");

  subcmd->add_option("-v,--max-nml-vaf", params->maxNmlVAF, "Max. variant allele frequency in normal sample", true)
      ->group("Filters");

  subcmd->add_option("-V,--min-tmr-vaf", params->minTmrVAF, "Min. variant allele frequency in tumor sample", true)
      ->group("Filters");

  subcmd->add_option("--min-nml-cov", params->minNmlCov, "Min. variant coverage in the normal sample", true)
      ->group("Filters");

  subcmd->add_option("--min-tmr-cov", params->minTmrCov, "Min. variant coverage in the tumor sample", true)
      ->group("Filters");

  subcmd->add_option("--max-nml-cov", params->maxNmlCov, "Max. variant coverage in the normal sample", true)
      ->group("Filters");

  subcmd->add_option("--max-tmr-cov", params->maxTmrCov, "Max. variant coverage in the tumor sample", true)
      ->group("Filters");

  subcmd->add_option("--min-fisher", params->minFisher, "Min. phred scaled score for somatic variants", true)
      ->group("Filters");

  subcmd->add_option("--min-str-fisher", params->minSTRFisher, "Min. phred scaled score for STR variants", true)
      ->group("Filters");

  subcmd->add_option("--min-strand-cnt", params->minStrandCnt, "Min. per strand contribution for a variant", true)
      ->group("Filters");

  // Feature flags
  subcmd->add_flag("--tenx-mode", params->tenxMode, "Run Lancet in 10X Linked Reads mode")->group("Flags");
  subcmd->add_flag("--active-region-off", params->activeRegionOff, "Turn off active region detection")->group("Flags");
  subcmd->add_flag("--kmer-recovery-on", params->kmerRecoveryOn, "Turn on experimental kmer recovery")->group("Flags");
  subcmd->add_flag("--xa-filter", params->skipMultipleHits, "Skip reads with XA tag (BWA-mem only)")->group("Flags");
  subcmd->add_flag("--skip-secondary", params->skipSecondary, "Skip secondary read alignments")->group("Flags");
  subcmd->add_flag("--extract-pairs", params->extractReadPairs, "Extract read pairs for each window")->group("Flags");
  subcmd->add_flag("--no-contig-check", params->noCtgCheck, "Skip checking for same contigs in BAM/CRAMs and reference")
      ->group("Flags");
  subcmd->add_flag("--use-overlap-reads", params->useOverlapReads, "Use reads overlapping windows to build graph")
      ->group("Flags");

  // Optional
  subcmd->add_option("--graphs-dir", params->outGraphsDir, "Output path to dump serialized graphs for the run", true)
      ->check(CLI::NonexistentPath)
      ->group("Optional");

  // clang-format off
  // http://patorjk.com/software/taag/#p=display&f=Big%20Money-nw&t=Lancet
  static constexpr auto logo = R"raw(
$$\                                               $$\
$$ |                                              $$ |
$$ |      $$$$$$\  $$$$$$$\   $$$$$$$\  $$$$$$\ $$$$$$\
$$ |      \____$$\ $$  __$$\ $$  _____|$$  __$$\\_$$  _|
$$ |      $$$$$$$ |$$ |  $$ |$$ /      $$$$$$$$ | $$ |
$$ |     $$  __$$ |$$ |  $$ |$$ |      $$   ____| $$ |$$\
$$$$$$$$\\$$$$$$$ |$$ |  $$ |\$$$$$$$\ \$$$$$$$\  \$$$$  |
\________|\_______|\__|  \__| \_______| \_______|  \____/

)raw";
  // clang-format on

  subcmd->callback([params]() -> void {
    const auto isTty = static_cast<bool>(isatty(fileno(stderr)));
    if (isTty) std::cerr << logo;
    RunPipeline(params);
  });
}
}  // namespace lancet
