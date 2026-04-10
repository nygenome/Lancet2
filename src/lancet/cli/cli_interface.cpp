#include "lancet/cli/cli_interface.h"

#include "CLI/CLI.hpp"
#include "absl/strings/match.h"
#include "absl/strings/str_cat.h"

#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>

extern "C" {
#include "htslib/hfile.h"
}

#include "lancet/base/logging.h"
#include "lancet/base/types.h"
#include "lancet/base/version.h"
#include "lancet/cbdg/graph.h"
#include "lancet/cli/cli_params.h"
#include "lancet/cli/pipeline_runner.h"
#include "lancet/core/window_builder.h"

#include "spdlog/common.h"
#include "spdlog/fmt/bundled/core.h"
#include "spdlog/fmt/bundled/format.h"
#include "spdlog/fmt/bundled/ostream.h"

namespace {

inline auto MakeCmdLine(int const argc, char const** argv) -> std::string {
  std::string result;
  static constexpr usize LINUX_MAX_CMDLINE_LENGTH = 2'097'152;
  result.reserve(LINUX_MAX_CMDLINE_LENGTH);
  absl::StrAppend(&result, argv[0]);
  for (auto idx = 1; idx < argc; ++idx) {
    absl::StrAppend(&result, " ", argv[idx]);
  }
  result.shrink_to_fit();
  return result;
}

class CliExistingUriOrFile : public CLI::Validator {
 public:
  CliExistingUriOrFile() {
    name_ = "EXISTING_URI_OR_FILE";
    func_ = [](std::string const& str) -> std::string {
      // By natively injecting HTSlib remote read pings here, we decouple the strict
      // local filepath assertions baked inherently into CLI::ExistingFile checking logic.
      // This allows Lancet2 pipeline execution to gracefully accept and validate
      // network streams (S3/GCS paths) uniformly alongside traditional local system binaries.
      if (absl::StartsWith(str, "gs://") ||
          absl::StartsWith(str, "s3://") ||
          absl::StartsWith(str, "http://") ||
          absl::StartsWith(str, "https://") ||
          absl::StartsWith(str, "ftp://") ||
          absl::StartsWith(str, "ftps://")) {
        hFILE* fptr = hopen(str.c_str(), "r");
        if (fptr == nullptr) {
          return fmt::format("Could not open existing cloud/web resource: {}", str);
        }
        if (hclose(fptr) < 0) {
          return fmt::format("Failed to properly close existing cloud/web resource connection: {}",
                             str);
        }
        return "";
      }
      return CLI::ExistingFile(str);
    };
  }
};

class CliNonexistentUriOrPath : public CLI::Validator {
 public:
  CliNonexistentUriOrPath() {
    name_ = "NONEXISTENT_URI_OR_PATH";
    func_ = [](std::string const& str) -> std::string {
      // Skips native CLI11 filesystem bounds checking specifically for cloud output schemas
      // ensuring that Lancet2 doesn't immediately abort if `--out-vcfgz` points completely
      // off-disk.
      if (absl::StartsWith(str, "gs://") ||
          absl::StartsWith(str, "s3://") ||
          absl::StartsWith(str, "http://") ||
          absl::StartsWith(str, "https://") ||
          absl::StartsWith(str, "ftp://") ||
          absl::StartsWith(str, "ftps://")) {
        return "";
      }
      return CLI::NonexistentPath(str);
    };
  }
};

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

CliInterface::CliInterface()
    : mCliApp(fmt::format(APP_NAME_FMT_STR, LancetFullVersion())),
      mParamsPtr(std::make_shared<CliParams>()) {
  // Allow 0 or 1 subcommands so that bare `Lancet2`, `Lancet2 -h`, and `Lancet2 -v`
  // can all be parsed without CLI11 throwing a RequiredError (exit code 106).
  mCliApp.require_subcommand(0, 1);
  PipelineSubcmd(&mCliApp, mParamsPtr);

  static constexpr usize DEFAULT_TERMINAL_FORMATTER_WIDTH = 100;
  mCliApp.get_formatter()->column_width(DEFAULT_TERMINAL_FORMATTER_WIDTH);
  mCliApp.failure_message(CLI::FailureMessage::help);

  // Use CLI11's native version/help flag handlers. These throw CallForVersion/CallForHelp
  // (both subclasses of Success with exit code 0) which are caught cleanly in RunMain.
  mCliApp.set_version_flag("-v,--version", LancetFullVersion(), "Print Lancet version information");
  mCliApp.set_help_flag("-h,--help", "Print this help message and exit");
}

auto CliInterface::RunMain(int const argc, char const** argv) -> int {
  try {
    mParamsPtr->mFullCmdLine = MakeCmdLine(argc, argv);
    mCliApp.parse(argc, argv);
  } catch (CLI::ParseError const& err) {
    return mCliApp.exit(err);
  }

  // If no subcommand was provided (bare `Lancet2`), print help and exit successfully.
  if (mCliApp.get_subcommands().empty()) {
    fmt::print(std::cout, "{}", mCliApp.help());
    return EXIT_SUCCESS;
  }

  return EXIT_SUCCESS;
}

void CliInterface::PipelineSubcmd(CLI::App* app, std::shared_ptr<CliParams>& params) {
  auto* subcmd = app->add_subcommand("pipeline", "Run the Lancet variant calling pipeline");
  subcmd->option_defaults()->always_capture_default();

  // Print help and exit 0 when pipeline subcommand is invoked with no arguments.
  subcmd->preparse_callback([subcmd](std::size_t remaining) -> void {
    if (remaining == 0) throw CLI::CallForHelp();
  });

  auto& wb_prms = params->mWindowBuilder;
  auto& vb_prms = params->mVariantBuilder;
  auto& rc_prms = vb_prms.mRdCollParams;
  auto& grph_prms = vb_prms.mGraphParams;

  static constexpr f64 MIN_TUMOR_VS_NORMAL_VAF_ODDS = 0.0;
  static constexpr f64 MAX_TUMOR_VS_NORMAL_VAF_ODDS = 255.0;
  static int const MAX_NUM_THREADS = static_cast<int>(std::thread::hardware_concurrency());

  // Datasets
  subcmd
      ->add_option("-n,--normal", rc_prms.mNormalPaths,
                   "Path to one (or) more normal BAM/CRAM file(s)")
      ->required(true)
      ->group("Datasets");
  subcmd
      ->add_option("-t,--tumor", rc_prms.mTumorPaths,
                   "Path to one (or) more tumor BAM/CRAM file(s)")
      ->required(false)
      ->group("Datasets");

  // Required
  subcmd->add_option("-r,--reference", rc_prms.mRefPath, "Path to the reference FASTA file")
      ->required(true)
      ->group("Required")
      ->check(CliExistingUriOrFile{});
  subcmd->add_option("-o,--out-vcfgz", params->mOutVcfGz, "Output path to the compressed VCF file")
      ->required(true)
      ->group("Required")
      ->check(CliExistingUriOrFile{} | CliNonexistentUriOrPath{});

  // Regions
  subcmd
      ->add_option("-R,--region", params->mInRegions,
                   "One (or) more regions (1-based both inclusive)")
      ->group("Regions")
      ->type_name("REF:[:START[-END]]");
  subcmd->add_option("-b,--bed-file", params->mBedFile, "Path to BED file with regions to process")
      ->group("Regions")
      ->check(CliExistingUriOrFile{});
  subcmd
      ->add_option("-P,--padding", wb_prms.mRegionPadding,
                   "Padding for both sides of all input regions")
      ->group("Regions")
      ->check(CLI::Range(static_cast<u32>(0), core::WindowBuilder::MAX_ALLOWED_REGION_PAD));
  subcmd
      ->add_option("-p,--pct-overlap", wb_prms.mPercentOverlap,
                   "Percent overlap between consecutive windows")
      ->group("Regions")
      ->check(CLI::Range(core::WindowBuilder::MIN_ALLOWED_PCT_OVERLAP,
                         core::WindowBuilder::MAX_ALLOWED_PCT_OVERLAP));
  subcmd
      ->add_option("-w,--window-size", params->mWindowBuilder.mWindowLength,
                   "Window size for variant calling tasks")
      ->group("Regions")
      ->check(CLI::Range(core::WindowBuilder::MIN_ALLOWED_WINDOW_LEN,
                         core::WindowBuilder::MAX_ALLOWED_WINDOW_LEN));

  // Parameters
  subcmd
      ->add_option("-T,--num-threads", params->mNumWorkerThreads,
                   "Number of additional async worker threads")
      ->group("Parameters")
      ->check(CLI::Range(0, MAX_NUM_THREADS));
  subcmd
      ->add_option("-k,--min-kmer", vb_prms.mGraphParams.mMinKmerLen,
                   "Min. kmer length to try for graph nodes")
      ->group("Parameters")
      ->check(CLI::Range(cbdg::Graph::DEFAULT_MIN_KMER_LEN, cbdg::Graph::MAX_ALLOWED_KMER_LEN - 2));
  subcmd
      ->add_option("-K,--max-kmer", vb_prms.mGraphParams.mMaxKmerLen,
                   "Max. kmer length to try for graph nodes")
      ->group("Parameters")
      ->check(CLI::Range(cbdg::Graph::DEFAULT_MIN_KMER_LEN + 2, cbdg::Graph::MAX_ALLOWED_KMER_LEN));
  subcmd
      ->add_option("-s,--kmer-step", vb_prms.mGraphParams.mKmerStepLen,
                   "Kmer step length to try for graph nodes")
      ->group("Parameters")
      ->check(CLI::IsMember({2, 4, 6, 8, 10}));
  subcmd
      ->add_option("--min-anchor-cov", grph_prms.mMinAnchorCov,
                   "Min. coverage for anchor nodes (source/sink)")
      ->group("Parameters")
      ->check(CLI::Range(static_cast<u32>(1), std::numeric_limits<u32>::max()));
  subcmd
      ->add_option("--min-node-cov", grph_prms.mMinNodeCov,
                   "Min. coverage for nodes kept in the graph")
      ->group("Parameters")
      ->check(CLI::Range(static_cast<u32>(0), std::numeric_limits<u32>::max()));
  subcmd
      ->add_option("--max-sample-cov", rc_prms.mMaxSampleCov,
                   "Max. per sample coverage before downsampling")
      ->group("Parameters")
      ->check(CLI::Range(static_cast<u32>(0), std::numeric_limits<u32>::max()));

  // Feature flags
  subcmd->add_flag("--verbose", params->mEnableVerboseLogging, "Turn on verbose logging")
      ->group("Flags");
  subcmd->add_flag("--extract-pairs", rc_prms.mExtractPairs, "Extract all useful read pairs")
      ->group("Flags");
  subcmd->add_flag("--no-active-region", vb_prms.mSkipActiveRegion, "Force assemble all windows")
      ->group("Flags");
  subcmd->add_flag("--no-contig-check", rc_prms.mNoCtgCheck, "Skip contig check with reference")
      ->group("Flags");

  // Optional
  subcmd
      ->add_option("--graphs-dir", vb_prms.mOutGraphsDir,
                   "Output directory to write per window graphs")
      ->check(CLI::NonexistentPath | CLI::ExistingDirectory)
      ->group("Optional");
  subcmd
      ->add_flag("--enable-graph-complexity-features", vb_prms.mEnableGraphComplexity,
                 "Emit GRAPH_CX INFO tag with per-variant graph complexity metrics")
      ->group("Optional");
  subcmd
      ->add_flag("--enable-sequence-complexity-features", vb_prms.mEnableSequenceComplexity,
                 "Emit ULTRA/MICRO/MACRO_*_CX INFO tags with multi-scale "
                 "sequence complexity metrics (HRun, entropy, TR motifs, LongdustQ)")
      ->group("Optional");
  subcmd
      ->add_option("--genome-gc-bias", vb_prms.mGcFraction,
                   "Global genome GC fraction for LongdustQ score correction. "
                   "Default: 0.41 (human genome-wide average). "
                   "Set to 0.5 to disable GC correction (uniform model).")
      ->group("Optional")
      ->check(CLI::Range(0.0, 1.0));

  subcmd->callback([params]() -> void {
    if (static_cast<bool>(isatty(fileno(stderr)))) fmt::print(std::cerr, FIGLET_LANCET_LOGO);
    if (params->mEnableVerboseLogging) SetLancetLoggerLevel(spdlog::level::trace);

    LOG_INFO("Starting Lancet {}", LancetFullVersion())
    PipelineRunner pipeline_runner(params);
    pipeline_runner.Run();
  });
}

}  // namespace lancet::cli
