#pragma once

#include <string>
#include <vector>

#include "lancet2/sized_ints.h"

namespace lancet2 {
constexpr double DEFAULT_MIN_NODE_COV_RATIO = 0.01;
constexpr double DEFAULT_MAX_WINDOW_COV = 1000.0;
constexpr double DEFAULT_MIN_TUMOR_VAF = 0.01;
constexpr double DEFAULT_MAX_NORMAL_VAF = 0.0;
constexpr double DEFAULT_MIN_PHRED_FISHER = 5.0;
constexpr double DEFAULT_MIN_PHRED_FISHER_STRS = 25.0;

constexpr u32 DEFAULT_NUM_WORKER_THREADS = 1;
constexpr u32 DEFAULT_REGION_PAD_LENGTH = 250;
constexpr u32 DEFAULT_WINDOW_LENGTH = 600;
constexpr u32 DEFAULT_PCT_WINDOW_OVERLAP = 50;
constexpr u32 DEFAULT_MIN_KMER_SIZE = 11;
constexpr u32 DEFAULT_MAX_KMER_SIZE = 101;
constexpr u32 DEFAULT_TRIM_BELOW_QUAL = 10;
constexpr u32 DEFAULT_MAX_GRAPH_TIP_LENGTH = 11;
constexpr u32 DEFAULT_MIN_REF_ANCHOR_COV = 5;
constexpr u32 DEFAULT_MIN_NODE_COV = 1;
constexpr u32 DEFAULT_GRAPH_TRAVERSAL_LIMIT = 1e5;
constexpr u32 DEFAULT_MAX_INDEL_LENGTH = 500;
constexpr u32 DEFAULT_MAX_REPEAT_MISMATCH = 2;
constexpr u32 DEFAULT_MIN_BASE_QUAL = 17;
constexpr u32 DEFAULT_MIN_READ_MAPPING_QUAL = 15;
constexpr u32 DEFAULT_MIN_STRAND_CONTRIB = 1;
constexpr u32 DEFAULT_MIN_TUMOR_ALT_CNT = 3;
constexpr u32 DEFAULT_MAX_NORMAL_ALT_CNT = 0;
constexpr u32 DEFAULT_MIN_TUMOR_COV = 4;
constexpr u32 DEFAULT_MIN_NORMAL_COV = 10;
constexpr u32 DEFAULT_MAX_TUMOR_COV = 1000;
constexpr u32 DEFAULT_MAX_NORMAL_COV = 1000;
constexpr u32 DEFAULT_MAX_STR_UNIT_LENGTH = 4;
constexpr u32 DEFAULT_MIN_STR_UNITS_TO_REPORT = 3;
constexpr u32 DEFAULT_MIN_STR_LENGTH_TO_REPORT = 7;
constexpr u32 DEFAULT_MAX_DIST_FROM_STR = 1;
constexpr u32 DEFAULT_MIN_READ_AS_XS_DIFF = 5;

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
  std::string outPrefix;               // NOLINT
  std::string commandLine;             // NOLINT

  double minCovRatio = DEFAULT_MIN_NODE_COV_RATIO;      // NOLINT
  double maxWindowCov = DEFAULT_MAX_WINDOW_COV;         // NOLINT
  double minTmrVAF = DEFAULT_MIN_TUMOR_VAF;             // NOLINT
  double maxNmlVAF = DEFAULT_MAX_NORMAL_VAF;            // NOLINT
  double minFisher = DEFAULT_MIN_PHRED_FISHER;          // NOLINT
  double minSTRFisher = DEFAULT_MIN_PHRED_FISHER_STRS;  // NOLINT

  u32 numWorkerThreads = DEFAULT_NUM_WORKER_THREADS;        // NOLINT
  u32 regionPadLength = DEFAULT_REGION_PAD_LENGTH;          // NOLINT
  u32 windowLength = DEFAULT_WINDOW_LENGTH;                 // NOLINT
  u32 pctOverlap = DEFAULT_PCT_WINDOW_OVERLAP;              // NOLINT
  u32 minKmerSize = DEFAULT_MIN_KMER_SIZE;                  // NOLINT
  u32 maxKmerSize = DEFAULT_MAX_KMER_SIZE;                  // NOLINT
  u32 trimBelowQual = DEFAULT_TRIM_BELOW_QUAL;              // NOLINT
  u32 minGraphTipLength = DEFAULT_MAX_GRAPH_TIP_LENGTH;     // NOLINT
  u32 minAnchorCov = DEFAULT_MIN_REF_ANCHOR_COV;            // NOLINT
  u32 minNodeCov = DEFAULT_MIN_NODE_COV;                    // NOLINT
  u32 graphTraversalLimit = DEFAULT_GRAPH_TRAVERSAL_LIMIT;  // NOLINT
  u32 maxIndelLength = DEFAULT_MAX_INDEL_LENGTH;            // NOLINT
  u32 maxRptMismatch = DEFAULT_MAX_REPEAT_MISMATCH;         // NOLINT
  u32 minBaseQual = DEFAULT_MIN_BASE_QUAL;                  // NOLINT
  u32 minReadMappingQual = DEFAULT_MIN_READ_MAPPING_QUAL;   // NOLINT
  u32 minStrandCnt = DEFAULT_MIN_STRAND_CONTRIB;            // NOLINT
  u32 minTmrAltCnt = DEFAULT_MIN_TUMOR_ALT_CNT;             // NOLINT
  u32 maxNmlAltCnt = DEFAULT_MAX_NORMAL_ALT_CNT;            // NOLINT
  u32 minTmrCov = DEFAULT_MIN_TUMOR_COV;                    // NOLINT
  u32 minNmlCov = DEFAULT_MIN_NORMAL_COV;                   // NOLINT
  u32 maxTmrCov = DEFAULT_MAX_TUMOR_COV;                    // NOLINT
  u32 maxNmlCov = DEFAULT_MAX_NORMAL_COV;                   // NOLINT
  u32 maxSTRUnitLength = DEFAULT_MAX_STR_UNIT_LENGTH;       // NOLINT
  u32 minSTRUnits = DEFAULT_MIN_STR_UNITS_TO_REPORT;        // NOLINT
  u32 minSTRLen = DEFAULT_MIN_STR_LENGTH_TO_REPORT;         // NOLINT
  u32 maxSTRDist = DEFAULT_MAX_DIST_FROM_STR;               // NOLINT
  u32 minReadAsXsDiff = DEFAULT_MIN_READ_AS_XS_DIFF;        // NOLINT

  bool verboseLogging = false;     // NOLINT
  bool activeRegionOff = false;    // NOLINT
  bool kmerRecoveryOn = false;     // NOLINT
  bool skipMultipleHits = false;   // NOLINT
  bool skipSecondary = false;      // NOLINT
  bool tenxMode = false;           // NOLINT
  bool extractReadPairs = false;   // NOLINT
  bool noCtgCheck = false;         // NOLINT
  bool useContainedReads = false;  // NOLINT
};
}  // namespace lancet2
