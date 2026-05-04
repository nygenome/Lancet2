# ═══════════════════════════════════════════════════════════════════════════════
# Lancet2 — Test Suite (Catch2 v3.14.0)
#
# Test sources grouped by module layer (matching the library dependency order):
#   base → hts → cbdg → caller → core
# ═══════════════════════════════════════════════════════════════════════════════
set(LANCET_FULL_DATA_DIR "${PROJECT_SOURCE_DIR}/data")
set(LANCET_TEST_CONFIG_H "${CMAKE_BINARY_DIR}/generated/lancet_test_config.h")
configure_file(test_config.h.inc ${LANCET_TEST_CONFIG_H} @ONLY)

# ── Main test executable ──────────────────────────────────────────────────────
add_executable(TestLancet2
		# Layer 1: base — math, repeat, sequence complexity, gzip + tar.gz I/O
		base/repeat_test.cpp
		base/longdust_scorer_test.cpp
		base/sequence_complexity_test.cpp
		base/mann_whitney_test.cpp
		base/gzip_ostream_test.cpp
		base/tar_gz_writer_test.cpp
		# Layer 2: hts — CIGAR, reference, alignment, extractor
		hts/cigar_utils_test.cpp
		hts/reference_test.cpp
		hts/alignment_test.cpp
		hts/extractor_test.cpp
		# Layer 3: cbdg — k-mer, graph, complexity, sample mask, dot renderer
		cbdg/kmer_test.cpp
		cbdg/sample_mask_test.cpp
		cbdg/graph_complexity_test.cpp
		cbdg/graph_test.cpp
		cbdg/dot_renderer_test.cpp
		# Layer 4: caller — variant set, support metrics, VCF output
		caller/variant_set_test.cpp
		caller/variant_support_metrics_test.cpp
		caller/variant_call_test.cpp
		# Layer 5: core — per-worker shard merge after compute phase
		core/tar_gz_shard_merger_test.cpp
		# External: longdust C sources for cross-validation
		${longdust_SOURCE_DIR}/longdust.c
		${longdust_SOURCE_DIR}/kalloc.c)
target_include_directories(TestLancet2 PRIVATE "${CMAKE_BINARY_DIR}/generated" "${CMAKE_SOURCE_DIR}"
		"${PROJECT_SOURCE_DIR}/tests" "${longdust_SOURCE_DIR}")
target_link_libraries(TestLancet2 PRIVATE Catch2 absl::strings lancet_cli zlibstatic)
# LANCET_DEBUG_MODE enables LANCET_ASSERT (see src/lancet/base/assert.h). The macro
# is normally tied to Debug builds via CMakeLists.txt, but the project's pre-commit
# gate runs the Release-tree test binary (reusing the same build that lint-check /
# iwyu-check produced). Defining the macro on the test target unconditionally keeps
# assertion coverage regardless of CMAKE_BUILD_TYPE.
target_compile_definitions(TestLancet2 PRIVATE LANCET_DEBUG_MODE)
set_target_properties(TestLancet2 PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
set_property(TARGET TestLancet2 PROPERTY $<$<ENABLE_LTO>:INTERPROCEDURAL_OPTIMIZATION TRUE>)

# ── Standalone Longdust tools (calibration utilities) ─────────────────────────
# Score computation tool: reads BAM regions, outputs scored TSV
add_executable(ScoreRegionLongdust tools/score_region_longdust.cpp)
target_include_directories(ScoreRegionLongdust PRIVATE "${CMAKE_BINARY_DIR}/generated" "${CMAKE_SOURCE_DIR}")
target_link_libraries(ScoreRegionLongdust
		PRIVATE lancet_hts concurrentqueue absl::btree absl::fixed_array absl::strings absl::synchronization)
set_target_properties(ScoreRegionLongdust PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

# Score analysis tool: reads scored TSV, generates aggregate statistics
add_executable(AnalyzeLongdustScores tools/analyze_longdust_scores.cpp)
target_include_directories(AnalyzeLongdustScores PRIVATE "${CMAKE_BINARY_DIR}/generated" "${CMAKE_SOURCE_DIR}")
target_link_libraries(AnalyzeLongdustScores PRIVATE lancet_hts absl::btree absl::flat_hash_set absl::strings)
set_target_properties(AnalyzeLongdustScores PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
