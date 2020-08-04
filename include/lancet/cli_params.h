#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace lancet {
constexpr double DEFAULT_MIN_NODE_COV_RATIO = 0.01F;
constexpr double DEFAULT_MAX_WINDOW_COV = 1000.0F;
constexpr double DEFAULT_MIN_TUMOR_VAF = 0.04F;
constexpr double DEFAULT_MAX_NORMAL_VAF = 0.0F;
constexpr double DEFAULT_MIN_PHRED_FISHER = 5.0F;
constexpr double DEFAULT_MIN_PHRED_FISHER_STRS = 25.0F;

constexpr std::uint32_t DEFAULT_NUM_WORKER_THREADS = 1;
constexpr std::uint32_t DEFAULT_REGION_PAD_LENGTH = 250;
constexpr std::uint32_t DEFAULT_WINDOW_LENGTH = 600;
constexpr std::uint32_t DEFAULT_PCT_WINDOW_OVERLAP = 84;
constexpr std::uint32_t DEFAULT_MIN_KMER_SIZE = 11;
constexpr std::uint32_t DEFAULT_MAX_KMER_SIZE = 101;
constexpr std::uint32_t DEFAULT_TRIM_BELOW_QUAL = 10;
constexpr std::uint32_t DEFAULT_MAX_GRAPH_TIP_LENGTH = 11;
constexpr std::uint32_t DEFAULT_MIN_REF_ANCHOR_COV = 5;
constexpr std::uint32_t DEFAULT_MIN_NODE_COV = 1;
constexpr std::uint32_t DEFAULT_GRAPH_TRAVERSAL_LIMIT = 1e6;
constexpr std::uint32_t DEFAULT_MAX_INDEL_LENGTH = 500;
constexpr std::uint32_t DEFAULT_MAX_REPEAT_MISMATCH = 2;
constexpr std::uint32_t DEFAULT_MIN_BASE_QUAL = 17;
constexpr std::uint32_t DEFAULT_MIN_READ_MAPPING_QUAL = 15;
constexpr std::uint32_t DEFAULT_MIN_STRAND_CONTRIB = 1;
constexpr std::uint32_t DEFAULT_MIN_TUMOR_ALT_CNT = 3;
constexpr std::uint32_t DEFAULT_MAX_NORMAL_ALT_CNT = 0;
constexpr std::uint32_t DEFAULT_MIN_TUMOR_COV = 4;
constexpr std::uint32_t DEFAULT_MIN_NORMAL_COV = 10;
constexpr std::uint32_t DEFAULT_MAX_TUMOR_COV = 1000;
constexpr std::uint32_t DEFAULT_MAX_NORMAL_COV = 1000;
constexpr std::uint32_t DEFAULT_MAX_STR_UNIT_LENGTH = 4;
constexpr std::uint32_t DEFAULT_MIN_STR_UNITS_TO_REPORT = 3;
constexpr std::uint32_t DEFAULT_MIN_STR_LENGTH_TO_REPORT = 7;
constexpr std::uint32_t DEFAULT_MAX_DIST_FROM_STR = 1;
constexpr std::uint32_t DEFAULT_MIN_READ_AS_XS_DIFF = 5;

class CliParams {
 public:
  CliParams() = default;

  [[nodiscard]] auto ValidateParams() -> bool;

  std::vector<std::string> inRegions;  // NOLINT
  std::string bedFilePath;             // NOLINT
  std::string outGraphsDir;            // NOLINT
  std::string referencePath;           // NOLINT
  std::string tumorPath;               // NOLINT
  std::string normalPath;              // NOLINT
  std::string outVcfPath;              // NOLINT
  std::string commandLine;             // NOLINT

  double minCovRatio = DEFAULT_MIN_NODE_COV_RATIO;      // NOLINT
  double maxWindowCov = DEFAULT_MAX_WINDOW_COV;         // NOLINT
  double minTmrVAF = DEFAULT_MIN_TUMOR_VAF;             // NOLINT
  double maxNmlVAF = DEFAULT_MAX_NORMAL_VAF;            // NOLINT
  double minFisher = DEFAULT_MIN_PHRED_FISHER;          // NOLINT
  double minSTRFisher = DEFAULT_MIN_PHRED_FISHER_STRS;  // NOLINT

  std::uint32_t numWorkerThreads = DEFAULT_NUM_WORKER_THREADS;        // NOLINT
  std::uint32_t regionPadLength = DEFAULT_REGION_PAD_LENGTH;          // NOLINT
  std::uint32_t windowLength = DEFAULT_WINDOW_LENGTH;                 // NOLINT
  std::uint32_t pctOverlap = DEFAULT_PCT_WINDOW_OVERLAP;              // NOLINT
  std::uint32_t minKmerSize = DEFAULT_MIN_KMER_SIZE;                  // NOLINT
  std::uint32_t maxKmerSize = DEFAULT_MAX_KMER_SIZE;                  // NOLINT
  std::uint32_t trimBelowQual = DEFAULT_TRIM_BELOW_QUAL;              // NOLINT
  std::uint32_t minGraphTipLength = DEFAULT_MAX_GRAPH_TIP_LENGTH;     // NOLINT
  std::uint32_t minAnchorCov = DEFAULT_MIN_REF_ANCHOR_COV;            // NOLINT
  std::uint32_t minNodeCov = DEFAULT_MIN_NODE_COV;                    // NOLINT
  std::uint32_t graphTraversalLimit = DEFAULT_GRAPH_TRAVERSAL_LIMIT;  // NOLINT
  std::uint32_t maxIndelLength = DEFAULT_MAX_INDEL_LENGTH;            // NOLINT
  std::uint32_t maxRptMismatch = DEFAULT_MAX_REPEAT_MISMATCH;         // NOLINT
  std::uint32_t minBaseQual = DEFAULT_MIN_BASE_QUAL;                  // NOLINT
  std::uint32_t minReadMappingQual = DEFAULT_MIN_READ_MAPPING_QUAL;   // NOLINT
  std::uint32_t minStrandCnt = DEFAULT_MIN_STRAND_CONTRIB;            // NOLINT
  std::uint32_t minTmrAltCnt = DEFAULT_MIN_TUMOR_ALT_CNT;             // NOLINT
  std::uint32_t maxNmlAltCnt = DEFAULT_MAX_NORMAL_ALT_CNT;            // NOLINT
  std::uint32_t minTmrCov = DEFAULT_MIN_TUMOR_COV;                    // NOLINT
  std::uint32_t minNmlCov = DEFAULT_MIN_NORMAL_COV;                   // NOLINT
  std::uint32_t maxTmrCov = DEFAULT_MAX_TUMOR_COV;                    // NOLINT
  std::uint32_t maxNmlCov = DEFAULT_MAX_NORMAL_COV;                   // NOLINT
  std::uint32_t maxSTRUnitLength = DEFAULT_MAX_STR_UNIT_LENGTH;       // NOLINT
  std::uint32_t minSTRUnits = DEFAULT_MIN_STR_UNITS_TO_REPORT;        // NOLINT
  std::uint32_t minSTRLen = DEFAULT_MIN_STR_LENGTH_TO_REPORT;         // NOLINT
  std::uint32_t maxSTRDist = DEFAULT_MAX_DIST_FROM_STR;               // NOLINT
  std::uint32_t minReadAsXsDiff = DEFAULT_MIN_READ_AS_XS_DIFF;        // NOLINT

  bool activeRegionOff = false;   // NOLINT
  bool kmerRecoveryOn = false;    // NOLINT
  bool skipMultipleHits = false;  // NOLINT
  bool skipSecondary = false;     // NOLINT
  bool tenxMode = false;          // NOLINT
  bool extractReadPairs = false;  // NOLINT
  bool noCtgCheck = false;        // NOLINT
  bool skipTruncSeq = false;      // NOLINT
};
}  // namespace lancet
