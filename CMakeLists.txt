set(LANCET_FULL_DATA_DIR "${PROJECT_SOURCE_DIR}/data")
set(LANCET_TEST_DATA_DIR "${PROJECT_SOURCE_DIR}/tests/data")
set(LANCET_TEST_CONFIG_H "${CMAKE_BINARY_DIR}/generated/lancet_test_config.h")
configure_file(test_config.h.inc ${LANCET_TEST_CONFIG_H} @ONLY)

add_executable(TestLancet2 base/repeat.cpp hts/reference_test.cpp
		hts/extractor_test.cpp hts/alignment_test.cpp cbdg/kmer_test.cpp)
target_include_directories(TestLancet2 PRIVATE "${CMAKE_BINARY_DIR}/generated" "${CMAKE_SOURCE_DIR}")
target_link_libraries(TestLancet2 PRIVATE Catch2 absl::strings lancet_cli)
set_target_properties(TestLancet2 PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
set_property(TARGET TestLancet2 PROPERTY $<$<ENABLE_LTO>:INTERPROCEDURAL_OPTIMIZATION TRUE>)

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
	target_compile_definitions(TestLancet2 PRIVATE -DLANCET_DEVELOP_MODE)
endif ()
