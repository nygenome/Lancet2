#include "lancet/cli/cli_interface.h"

#include "lancet/base/logging.h"
#include "lancet/base/types.h"
#include "lancet/base/version.h"
#include "lancet/cbdg/dot_plan.h"
#include "lancet/cbdg/graph_params.h"
#include "lancet/cli/cli_params.h"
#include "lancet/cli/pipeline_runner.h"
#include "lancet/core/window_builder.h"
#include "lancet/hts/uri_utils.h"

#include "CLI/CLI.hpp"
#include "absl/strings/str_cat.h"
#include "spdlog/common.h"
#include "spdlog/fmt/bundled/format.h"
#include "spdlog/fmt/bundled/ostream.h"

// POSIX header — not part of the C/C++ standard library
extern "C" {
#include "unistd.h"
}

#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <cstdio>
#include <cstdlib>

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
      if (lancet::hts::IsCloudUri(str)) return lancet::hts::ValidateCloudAccess(str, "r");
      return CLI::ExistingFile(str);
    };
  }
};

class CliNonexistentUriOrPath : public CLI::Validator {
 public:
  CliNonexistentUriOrPath() {
    name_ = "NONEXISTENT_URI_OR_PATH";
    func_ = [](std::string const& str) -> std::string {
      if (lancet::hts::IsCloudUri(str)) return std::string{};
      return CLI::NonexistentPath(str);
    };
  }
};

// Custom validator: enforce a `.tar.gz` suffix. The path itself can be
// any location (existing or not); the parent dir is created lazily by
// PipelineRunner::SetupGraphOutputArchive.
class CliTarGzSuffixValidator : public CLI::Validator {
 private:
  static constexpr std::string_view EXPECTED_SUFFIX = ".tar.gz";
  static constexpr std::string_view ERR_MSG_FMT =
      "--out-graphs-tgz path must end in `.tar.gz` (got: {})";

 public:
  CliTarGzSuffixValidator() {
    name_ = "TAR_GZ_SUFFIX";
    func_ = [](std::string const& str) -> std::string {
      if (!str.ends_with(EXPECTED_SUFFIX)) return fmt::format(ERR_MSG_FMT, str);
      return std::string{};
    };
  }
};

// ============================================================================
// CLI option groups
// ============================================================================
constexpr auto GRP_DATASETS = "Datasets";
constexpr auto GRP_REQUIRED = "Required";
constexpr auto GRP_REGIONS = "Regions";
constexpr auto GRP_PARAMETERS = "Parameters";
constexpr auto GRP_FLAGS = "Flags";
constexpr auto GRP_OPTIONAL = "Optional";

/// Register a CLI option with standard required/group boilerplate.
/// Returns CLI::Option* for optional ->check() or ->type_name() chaining.
template <typename T>
auto AddOpt(CLI::App* sub, std::string_view flags, T& target, std::string_view desc,
            std::string_view group, bool required = false) -> CLI::Option* {
  auto* opt = sub->add_option(std::string(flags), target, std::string(desc));
  opt->required(required);
  opt->group(std::string(group));
  return opt;
}

/// Register a CLI boolean flag with standard group boilerplate.
template <typename T>
auto AddFlag(CLI::App* sub, std::string_view flags, T& target, std::string_view desc,
             std::string_view group) -> CLI::Option* {
  auto* opt = sub->add_flag(std::string(flags), target, std::string(desc));
  opt->group(std::string(group));
  return opt;
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

static constexpr auto APP_NAME_FMT_STR = "Lancet, {}\nMicroassembly based variant caller\n";

namespace lancet::cli {

CliInterface::CliInterface()
    : mCliApp(fmt::format(APP_NAME_FMT_STR, lancet::base::LancetFullVersion())),
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
  mCliApp.set_version_flag("-v,--version", lancet::base::LancetFullVersion(),
                           "Print Lancet version information");
  mCliApp.set_help_flag("-h,--help", "Print this help message and exit");
}

auto CliInterface::RunMain(int const argc, char const** argv) -> int {
  try {
    mParamsPtr->mFullCmdLine = MakeCmdLine(argc, argv);
    mCliApp.parse(argc, argv);
  } catch (CLI::ParseError const& err) {
    // Delegate to the parsed subcommand so that its help is shown, not the root app's.
    auto const& subs = mCliApp.get_subcommands();
    if (!subs.empty()) return subs.front()->exit(err);
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
  auto* sub = app->add_subcommand("pipeline", "Run the Lancet variant calling pipeline");
  sub->option_defaults()->always_capture_default();
  sub->failure_message(CLI::FailureMessage::help);

  // Print help and exit 0 when pipeline subcommand is invoked with no arguments.
  sub->preparse_callback([sub](std::size_t remaining) -> void {
    if (remaining == 0) throw CLI::CallForHelp();
  });

  auto& win_params = params->mWindowBuilder;
  auto& var_params = params->mVariantBuilder;
  auto& rc_params = var_params.mRdCollParams;
  auto& graph_params = var_params.mGraphParams;

  static int const MAX_THREADS = static_cast<int>(std::thread::hardware_concurrency());

  // ============================================================================
  // Datasets
  // ============================================================================
  AddOpt(sub, "-n,--normal", rc_params.mCtrlPaths,
         "Path to one (or) more control (normal) BAM/CRAM file(s)", GRP_DATASETS, true);
  AddOpt(sub, "-t,--tumor", rc_params.mCasePaths,
         "Path to one (or) more case (tumor) BAM/CRAM file(s)", GRP_DATASETS);
  AddOpt(sub, "-s,--sample", rc_params.mSampleSpecs,
         "Sample input as <path>:<role> (roles: control, case)", GRP_DATASETS);

  // ============================================================================
  // Required
  // ============================================================================
  AddOpt(sub, "-r,--reference", rc_params.mRefPath, "Path to the reference FASTA file",
         GRP_REQUIRED, true)
      ->check(CliExistingUriOrFile{});
  AddOpt(sub, "-o,--out-vcfgz", params->mOutVcfGz, "Output path to the compressed VCF file",
         GRP_REQUIRED, true)
      ->check(CliExistingUriOrFile{} | CliNonexistentUriOrPath{});

  // ============================================================================
  // Regions
  // ============================================================================
  AddOpt(sub, "-R,--region", params->mInRegions, "One (or) more regions (1-based both inclusive)",
         GRP_REGIONS)
      ->type_name("REF:[:START[-END]]");
  AddOpt(sub, "-b,--bed-file", params->mBedFile, "Path to BED file with regions to process",
         GRP_REGIONS)
      ->check(CliExistingUriOrFile{});
  AddOpt(sub, "-P,--padding", win_params.mRegionPadding,
         "Padding for both sides of all input regions", GRP_REGIONS)
      ->check(CLI::Range(u32{0}, core::WindowBuilder::MAX_ALLOWED_REGION_PAD));
  AddOpt(sub, "-p,--pct-overlap", win_params.mPercentOverlap,
         "Percent overlap between consecutive windows", GRP_REGIONS)
      ->check(CLI::Range(core::WindowBuilder::MIN_ALLOWED_PCT_OVERLAP,
                         core::WindowBuilder::MAX_ALLOWED_PCT_OVERLAP));
  AddOpt(sub, "-w,--window-size", win_params.mWindowLength, "Window size for variant calling tasks",
         GRP_REGIONS)
      ->check(CLI::Range(core::WindowBuilder::MIN_ALLOWED_WINDOW_LEN,
                         core::WindowBuilder::MAX_ALLOWED_WINDOW_LEN));

  // ============================================================================
  // Parameters
  // ============================================================================
  AddOpt(sub, "-T,--num-threads", params->mNumWorkerThreads,
         "Number of additional async worker threads", GRP_PARAMETERS)
      ->check(CLI::Range(0, MAX_THREADS));
  AddOpt(sub, "-k,--min-kmer", graph_params.mMinKmerLen, "Min. kmer length to try for graph nodes",
         GRP_PARAMETERS)
      ->check(CLI::Range(cbdg::DEFAULT_MIN_KMER_LEN, cbdg::MAX_ALLOWED_KMER_LEN - 2));
  AddOpt(sub, "-K,--max-kmer", graph_params.mMaxKmerLen, "Max. kmer length to try for graph nodes",
         GRP_PARAMETERS)
      ->check(CLI::Range(cbdg::DEFAULT_MIN_KMER_LEN + 2, cbdg::MAX_ALLOWED_KMER_LEN));
  AddOpt(sub, "--kmer-step", graph_params.mKmerStepLen, "Kmer step length to try for graph nodes",
         GRP_PARAMETERS)
      ->check(CLI::IsMember({2, 4, 6, 8, 10}));
  AddOpt(sub, "--min-anchor-cov", graph_params.mMinAnchorCov,
         "Min. coverage for anchor nodes (source/sink)", GRP_PARAMETERS)
      ->check(CLI::Range(u32{1}, std::numeric_limits<u32>::max()));
  AddOpt(sub, "--min-node-cov", graph_params.mMinNodeCov,
         "Min. coverage for nodes kept in the graph", GRP_PARAMETERS)
      ->check(CLI::Range(u32{0}, std::numeric_limits<u32>::max()));
  AddOpt(sub, "--max-sample-cov", rc_params.mMaxSampleCov,
         "Max. per sample coverage before downsampling", GRP_PARAMETERS)
      ->check(CLI::Range(u32{0}, std::numeric_limits<u32>::max()));

  // ============================================================================
  // Flags
  // ============================================================================
  AddFlag(sub, "--verbose", params->mEnableVerboseLogging, "Turn on verbose logging", GRP_FLAGS);
  AddFlag(sub, "--extract-pairs", rc_params.mExtractPairs, "Extract all useful read pairs",
          GRP_FLAGS);
  AddFlag(sub, "--no-active-region", var_params.mSkipActiveRegion, "Force assemble all windows",
          GRP_FLAGS);
  AddFlag(sub, "--no-contig-check", rc_params.mNoCtgCheck, "Skip contig check with reference",
          GRP_FLAGS);

  // ============================================================================
  // Optional
  // ============================================================================
  static auto const SNAPSHOT_MAP = std::map<std::string, cbdg::GraphSnapshotMode>{
      {"final", cbdg::GraphSnapshotMode::FINAL}, {"verbose", cbdg::GraphSnapshotMode::VERBOSE}};

  AddOpt(sub, "--out-graphs-tgz", var_params.mOutGraphsTgz,
         "Output path for the tar.gz archive of per-window assembly graphs.", GRP_OPTIONAL)
      ->check(CliTarGzSuffixValidator{});
  AddOpt(sub, "--graph-snapshots", graph_params.mSnapshotMode,
         "Control the verbosity of per-window assembly graph snapshots.", GRP_OPTIONAL)
      ->transform(CLI::CheckedTransformer(SNAPSHOT_MAP, CLI::ignore_case));
  AddOpt(sub, "--genome-gc-bias", var_params.mGcFraction,
         "Global genome GC fraction for LongdustQ score correction. Default 0.41", GRP_OPTIONAL)
      ->check(CLI::Range(0.0, 1.0));

  auto* probe_variants_opt =
      AddOpt(sub, "--probe-variants", var_params.mProbeVariantsPath,
             "Path to missed_variants.txt for k-mer probe diagnostics", GRP_OPTIONAL)
          ->check(CLI::ExistingFile);

  auto* probe_results_opt = AddOpt(sub, "--probe-results", var_params.mProbeResultsPath,
                                   "Output path for k-mer probe results TSV", GRP_OPTIONAL);

  // Bidirectional: both must be provided together or neither.
  probe_variants_opt->needs(probe_results_opt);
  probe_results_opt->needs(probe_variants_opt);

  // ============================================================================
  // Subcommand callback
  // ============================================================================
  sub->callback([params]() -> void {
    if (static_cast<bool>(isatty(fileno(stderr)))) fmt::print(std::cerr, FIGLET_LANCET_LOGO);
    if (params->mEnableVerboseLogging) SetLancetLoggerLevel(spdlog::level::trace);

    LOG_INFO("Starting Lancet {}", lancet::base::LancetFullVersion())
    PipelineRunner pipeline_runner(params);
    pipeline_runner.Run();
  });
}

}  // namespace lancet::cli
