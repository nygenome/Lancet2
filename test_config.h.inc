#ifndef LANCET_TEST_CONFIG_H_INC
#define LANCET_TEST_CONFIG_H_INC

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <utility>

#include "spdlog/fmt/fmt.h"

static constexpr auto FULL_DATA_DIR = "@LANCET_FULL_DATA_DIR@";
static constexpr auto TEST_DATA_DIR = "@LANCET_TEST_DATA_DIR@";

static constexpr auto GRCH38_REF_NAME = "GRCh38_full_analysis_set_plus_decoy_hla.fa.gz";
static constexpr auto TUMOR_CRAM_NAME = "chr4_with_pairs.HCC1395_SAMN10102573_SRR7890893.cram";
static constexpr auto NORMAL_CRAM_NAME = "chr4_with_pairs.HCC1395BL_SAMN10102574_SRR7890943.cram";
static constexpr auto TUMOR_BAM_NAME = "chr4_with_pairs.HCC1395_SAMN10102573_SRR7890893.bam";
static constexpr auto NORMAL_BAM_NAME = "chr4_with_pairs.HCC1395BL_SAMN10102574_SRR7890943.bam";

[[nodiscard]] static inline auto MakePath(std::string_view prefix, std::string_view suffix) -> std::filesystem::path {
  return fmt::format("{}/{}", prefix, suffix);
}

static constexpr auto GRCH38_NAMES_AUTOSOMES_XY = std::array<const char*, 25>{
    "chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9",  "chr10", "chr11", "chr12",
    "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",  "chrY"};

static constexpr auto GRCH38_LENGTHS_AUTOSOMES_XY = std::array<std::uint64_t, 25>{
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
    138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
    83257441,  80373285,  58617616,  64444167,  46709983,  50818468,  156040895, 57227415};

#endif  // LANCET_TEST_CONFIG_H_INC
