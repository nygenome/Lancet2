#include "lancet/cli/cli_interface.h"

#include <unistd.h>

#include <cstdlib>
#include <exception>
#include <limits>
#include <string>
#include <thread>

#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "lancet/base/logging.h"
#include "lancet/base/types.h"
#include "lancet/base/version.h"
#include "lancet/cli/pipeline_runner.h"
#include "spdlog/fmt/fmt.h"
#include "spdlog/fmt/ostr.h"

namespace {

inline auto MakeCmdLine(const int argc, const char** argv) -> std::string {
  std::string result;
  static constexpr usize LINUX_MAX_CMDLINE_LENGTH = 2097152;
  result.reserve(LINUX_MAX_CMDLINE_LENGTH);
  absl::StrAppend(&result, argv[0]);           // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  for (auto idx = 1; idx < argc; ++idx) {
    absl::StrAppend(&result, " ", argv[idx]);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  }
  result.shrink_to_fit();
  return result;
}

}  // namespace

// http://patorjk.com/software/taag/#p=display&f=Big%20Money-nw&t=Lancet
// clang-format off
  static constexpr auto FIGLET_LANCET_LOGO = R"raw(
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

static constexpr auto APP_NAME_FMT_STR = "Lancet, {}\nMicrosssembly based somatic variant caller\n";

namespace lancet::cli {

namespace detail {

struct ThreadCountValidator : public CLI::Validator {
  ThreadCountValidator() {
    name_ = "VALID_THREAD_COUNT";
    func_ = [](const std::string& val) {
      i64 num = 0;
      absl::SimpleAtoi(val, &num);
      static const auto max_threads = std::thread::hardware_concurrency();
      return num < 1             ? std::string("thread count cannot be less than 1")
             : num > max_threads ? std::string("thread count cannot be greater than physical cores")
                                 : std::string();
    };
  }
};

}  // namespace detail

CliInterface::CliInterface()
    : mCliApp(fmt::format(APP_NAME_FMT_STR, LancetFullVersion())), mParamsPtr(std::make_shared<CliParams>()) {
  mCliApp.require_subcommand(1);
  PipelineSubcmd(&mCliApp, mParamsPtr);

  static const auto version_printer = [](const usize count) -> void {
    if (count > 0) {
      fmt::print(std::cout, "Lancet {}\n", LancetFullVersion());
      std::exit(EXIT_SUCCESS);
    }
  };

  static const auto help_printer = [this](const usize count) -> void {
    if (count > 0) {
      fmt::print(std::cerr, "{}", this->mCliApp.help(this->mCliApp.get_name(), CLI::AppFormatMode::Normal));
      std::exit(EXIT_SUCCESS);
    }
  };

  mCliApp.set_help_flag();
  mCliApp.failure_message(CLI::FailureMessage::help);
  static constexpr usize DEFAULT_TERMINAL_FORMATTER_WIDTH = 70;
  mCliApp.get_formatter()->column_width(DEFAULT_TERMINAL_FORMATTER_WIDTH);
  mCliApp.add_flag_function("-v,--version", version_printer, "Print Lancet version information")->group("Flags");
  mCliApp.add_flag_function("-h,--help", help_printer, "Print this help message and exit")->group("Flags");
}

auto CliInterface::RunMain(const int argc, const char** argv) -> int {
  try {
    mParamsPtr->mFullCmdLine = MakeCmdLine(argc, argv);
    mCliApp.parse(argc, argv);
  } catch (const CLI::ParseError& err) {
    return mCliApp.exit(err);
  }

  return EXIT_SUCCESS;
}

void CliInterface::PipelineSubcmd(CLI::App* app, std::shared_ptr<CliParams>& params) {
  auto* subcmd = app->add_subcommand("pipeline", "Run the Lancet variant calling pipeline");
  subcmd->option_defaults()->always_capture_default();

  auto& wb_prms = params->mWindowBuilder;
  auto& vb_prms = params->mVariantBuilder;
  auto& rc_prms = vb_prms.mRdCollParams;
  auto& grph_prms = vb_prms.mGraphParams;
  auto& fltr_prms = vb_prms.mVariantParams;

  // Datasets
  subcmd->add_option("-n,--normal", params->mNormalPaths, "Path to one (or) more normal BAM/CRAM file(s)")
      ->required(true)
      ->group("Datasets");
  subcmd->add_option("-t,--tumor", params->mTumorPaths, "Path to one (or) more tumor BAM/CRAM file(s)")
      ->required(false)
      ->group("Datasets");

  // Required
  subcmd->add_option("-r,--reference", rc_prms.mRefPath, "Path to the reference FASTA file")
      ->required(true)
      ->group("Required")
      ->check(CLI::ExistingFile);
  subcmd->add_option("-o,--out-vcfgz", params->mOutVcfGz, "Output path to the compressed VCF file")
      ->required(true)
      ->group("Required")
      ->check(CLI::ExistingFile | CLI::NonexistentPath);

  // Regions
  subcmd->add_option("-R,--region", params->mInRegions, "One (or) more regions (1-based both inclusive)")
      ->group("Regions")
      ->type_name("REF:[:START[-END]]");
  subcmd->add_option("-b,--bed-file", params->mBedFile, "Path to BED file with regions to process")
      ->group("Regions")
      ->check(CLI::ExistingFile);
  subcmd->add_option("--padding", wb_prms.mRegionPadding, "Padding for both sides of all input regions")
      ->group("Regions")
      ->check(CLI::Range(u32(0), core::WindowBuilder::MAX_ALLOWED_REGION_PAD));
  subcmd->add_option("-p,--pct-overlap", wb_prms.mPercentOverlap, "Percent overlap between consecutive windows")
      ->group("Regions")
      ->check(CLI::Range(core::WindowBuilder::MIN_ALLOWED_PCT_OVERLAP, core::WindowBuilder::MAX_ALLOWED_PCT_OVERLAP));
  subcmd->add_option("-w,--window-size", params->mWindowBuilder.mWindowLength, "Window size for variant calling tasks")
      ->group("Regions")
      ->check(CLI::Range(core::WindowBuilder::MIN_ALLOWED_WINDOW_LEN, core::WindowBuilder::MAX_ALLOWED_WINDOW_LEN));

  // Parameters
  static const detail::ThreadCountValidator THREAD_COUNT_VALIDATOR;
  subcmd->add_option("-T,--num-threads", params->mNumWorkerThreads, "Number of additional async worker threads")
      ->group("Parameters")
      ->check(THREAD_COUNT_VALIDATOR);
  subcmd->add_option("-k,--min-kmer", vb_prms.mGraphParams.mMinKmerLen, "Min. kmer length to try for graph nodes")
      ->group("Parameters")
      ->check(CLI::Range(cbdg::Graph::DEFAULT_MIN_KMER_LEN, cbdg::Graph::MAX_ALLOWED_KMER_LEN - 2));
  subcmd->add_option("-K,--max-kmer", vb_prms.mGraphParams.mMaxKmerLen, "Max. kmer length to try for graph nodes")
      ->group("Parameters")
      ->check(CLI::Range(cbdg::Graph::DEFAULT_MIN_KMER_LEN + 2, cbdg::Graph::MAX_ALLOWED_KMER_LEN));
  subcmd->add_option("--min-anchor-cov", grph_prms.mMinRefAnchorCov, "Min. coverage for anchor nodes (source/sink)")
      ->group("Parameters");
  subcmd->add_option("--min-node-cov", grph_prms.mMinNodeCoverage, "Min. coverage needed for all nodes in the graph")
      ->group("Parameters");
  subcmd->add_option("--min-cov-ratio", grph_prms.mMinNodeCovRatio, "Min. node to window coverage for all nodes")
      ->group("Parameters");
  subcmd->add_option("--max-win-cov", rc_prms.mMaxWinCov, "Max. combined window coverage for downsampling")
      ->group("Parameters");

  // Filters
  subcmd->add_option("--min-nml-cov", fltr_prms.mMinNmlCov, "Min. normal coverage")->group("Filters");
  subcmd->add_option("--min-tmr-cov", fltr_prms.mMinTmrCov, "Min. tumor coverage")->group("Filters");
  subcmd->add_option("--max-nml-alt-cnt", fltr_prms.mMaxNmlAltCnt, "Max. ALT allele count in normal")->group("Filters");
  subcmd->add_option("--min-tmr-alt-cnt", fltr_prms.mMinTmrAltCnt, "Min. ALT allele count in tumor")->group("Filters");
  subcmd->add_option("--max-nml-vaf", fltr_prms.mMaxNmlVaf, "Max. ALT frequency in normal")
      ->group("Filters")
      ->check(CLI::Range(0.0, 1.0));
  subcmd->add_option("--min-tmr-vaf", fltr_prms.mMinTmrVaf, "Min. ALT frequency in tumor")
      ->group("Filters")
      ->check(CLI::Range(0.0, 1.0));
  subcmd->add_option("--min-phred-score", fltr_prms.mMinPhredScore, "Min. phred score for PASS variants")
      ->group("Filters")
      ->check(CLI::Range(core::VariantBuilder::MIN_PHRED_SCORE, core::VariantBuilder::MAX_PHRED_SCORE));

  // Feature flags
  subcmd->add_flag("--verbose", params->mEnableVerboseLogging, "Turn on verbose logging")->group("Flags");
#ifndef LANCET_DEVELOP_MODE
  subcmd->add_flag("--cpu-profile", params->mEnableCpuProfiling, "Turn on CPU profiling")->group("Flags");
#endif
  subcmd->add_flag("--extract-pairs", rc_prms.mExtractReadPairs, "Extract all useful read pairs")->group("Flags");
  subcmd->add_flag("--no-active-region", vb_prms.mSkipActiveRegion, "Force assemble all windows")->group("Flags");
  subcmd->add_flag("--no-contig-check", rc_prms.mNoCtgCheck, "Skip contig check with reference")->group("Flags");
  subcmd->add_flag("--no-filter-reads", rc_prms.mNoFilterRds, "Use all sample reads skipping filters")->group("Flags");

  // Optional
  subcmd->add_option("--runtime-stats", params->mRunStats, "Output text file with per window runtime & status")
      ->check(CLI::NonexistentPath | CLI::ExistingFile)
      ->group("Optional");
  subcmd->add_option("--graphs-dir", vb_prms.mOutGraphsDir, "Output directory to write per window graphs")
      ->check(CLI::NonexistentPath | CLI::ExistingDirectory)
      ->group("Optional");

  subcmd->callback([params]() {
    // NOLINTBEGIN(readability-braces-around-statements)
    if (static_cast<bool>(isatty(fileno(stderr)))) fmt::print(std::cerr, FIGLET_LANCET_LOGO);
    if (params->mEnableVerboseLogging) SetLancetLoggerLevel(spdlog::level::trace);
    // NOLINTEND(readability-braces-around-statements)

    LOG_INFO("Starting Lancet {}", LancetFullVersion())
    LOG_INFO("Lancet Command line: {}", params->mFullCmdLine)

    PipelineRunner pipeline_runner(params);
    pipeline_runner.Run();
  });
}

}  // namespace lancet::cli