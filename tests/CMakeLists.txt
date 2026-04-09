set(LANCET_FULL_DATA_DIR "${PROJECT_SOURCE_DIR}/data")
set(LANCET_TEST_DATA_DIR "${PROJECT_SOURCE_DIR}/tests/data")
set(LANCET_TEST_CONFIG_H "${CMAKE_BINARY_DIR}/generated/lancet_test_config.h")
configure_file(test_config.h.inc ${LANCET_TEST_CONFIG_H} @ONLY)

# Fetch longdust C sources for cross-validation of LongdustQScorer.
# Longdust is a C library (no CMakeLists.txt), so we use FetchContent_Populate
# to download the sources and compile just the files we need.
FetchContent_Declare(longdust
		GIT_REPOSITORY https://github.com/lh3/longdust
		GIT_TAG v1.4
		SYSTEM)
FetchContent_GetProperties(longdust)
if (NOT longdust_POPULATED)
	FetchContent_Populate(longdust)
endif ()

add_executable(TestLancet2
		base/repeat_test.cpp base/longdust_scorer_test.cpp base/sequence_complexity_test.cpp
		hts/cigar_utils_test.cpp hts/reference_test.cpp hts/extractor_test.cpp hts/alignment_test.cpp
		cbdg/kmer_test.cpp cbdg/graph_test.cpp cbdg/graph_complexity_test.cpp caller/variant_set_test.cpp
		${longdust_SOURCE_DIR}/longdust.c ${longdust_SOURCE_DIR}/kalloc.c)
target_include_directories(TestLancet2 PRIVATE
		"${CMAKE_BINARY_DIR}/generated" "${CMAKE_SOURCE_DIR}"
		"${PROJECT_SOURCE_DIR}/tests" "${longdust_SOURCE_DIR}")
target_link_libraries(TestLancet2 PRIVATE Catch2 absl::strings lancet_cli)
set_target_properties(TestLancet2 PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
set_property(TARGET TestLancet2 PROPERTY $<$<ENABLE_LTO>:INTERPROCEDURAL_OPTIMIZATION TRUE>)

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
	target_compile_definitions(TestLancet2 PRIVATE -DLANCET_DEVELOP_MODE)
endif ()

# Standalone Longdust score computation tool (for full calibration scoring)
add_executable(ScoreRegionLongdust base/score_region_longdust.cpp)
target_include_directories(ScoreRegionLongdust PRIVATE
		"${CMAKE_BINARY_DIR}/generated" "${CMAKE_SOURCE_DIR}")
target_link_libraries(ScoreRegionLongdust PRIVATE lancet_hts concurrentqueue absl::btree absl::fixed_array absl::strings absl::synchronization)
set_target_properties(ScoreRegionLongdust PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

# Standalone Longdust score analysis tool (reads scored TSV, generates statistics)
add_executable(AnalyzeLongdustScores base/analyze_longdust_scores.cpp)
target_include_directories(AnalyzeLongdustScores PRIVATE
		"${CMAKE_BINARY_DIR}/generated" "${CMAKE_SOURCE_DIR}")
target_link_libraries(AnalyzeLongdustScores PRIVATE lancet_hts absl::btree absl::flat_hash_set absl::strings)
set_target_properties(AnalyzeLongdustScores PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
