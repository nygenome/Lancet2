#include "lancet/align.h"

#include <string_view>

#include "catch2/catch.hpp"
#include "spdlog/spdlog.h"

TEST_CASE("performs alignment using custom aligner", "align.h") {
  constexpr std::string_view ref =
      "AGCCACCCTCCAAACCACAGGGGTCATTCTGGAAACCAAGGCATATCTGTGAACAGGCAGGGTGTTCAGCATCAGAGAGGGGAGCAGGGGCTTCTCCGTGAGCCAGGAGG"
      "CAGGAAGAGACCTGCCCGAAGCTCAGCAGGACAGGGCCAGGCACAGCCGGGTAGCTGCTGGAAACACAGACAGGACCCGCTACTGGCCCAGCTCCCGATATGGCAGACCA"
      "CGATCTCATCCACTGCTCACAGCAGCCCCGTGAGGCAAGCATCATGACCCAAGGTGTACATGGGGAAACTGAGGTTCAGAGGTTTCTCTGTATGAGAAACAAAGGCTCAG"
      "AGAGGTTAGGCAGCCTGCCCAAGGACACACAGCATCAGCAGGTACTCCTGACTGGGGAGCACTGAGTTCACCTCTCCTTCGACAGACAGGGGAGCTTCCCCCCACTCCCT"
      "ACCTGGAGCCCATAGACCTTTGATAATGCCAGGGCTGCCTCTGAGTCTCAACTCACCCCTTGATTTCTTGCACCCCACACTGTACACTAGGCTTCCTGTTTTACACCTAT"
      "CATCTCATCTCATTTTTACAATGGCCCCATC";

  constexpr std::string_view qry =
      "GTTTGGAGGGTGGCTGCTGTATCTCTATTTTCTGTATGAGAAACAAAGGCTCAGAGAGGTTAGGCAGCCTGCCCAAGGACACACAGCATCAGCAGGTACTCCTGACTGGG"
      "GAGCACTGAGTTCACCTCTCCTTCGACAGACAGGGGAGCTTCCCCCCACTCCCTACCTGGAGCCCATAGACCTTTGATAATGCCAGGGCTGCCTCTGAGTCTCAACTCAC"
      "CCCTTGATTTCTTGCACCCCACACTGTACACTACGCTTCCTGTTTTACACCTATCATCTCATCTCATTTTTACAATGGCCCCATC";

  const auto result = lancet::Align(ref, qry);
  spdlog::info("first: {}, len: {}, qry: {}", static_cast<int>(result.ref[0]), result.ref.length(),
               result.qry.length());
  spdlog::info("cust ref: {}", result.ref);
  spdlog::info("cust qry: {}", result.qry);
  CHECK(result.ref.length() == result.qry.length());
}