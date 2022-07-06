set(LANCET2_TEST_DATA_DIR "${PROJECT_SOURCE_DIR}/tests/data")
configure_file(test_config.h.inc "${CMAKE_BINARY_DIR}/generated/test_config.h")

add_executable(lancet2_tests "${CMAKE_BINARY_DIR}/generated/test_config.h" lancet2_test.cc
        fractional_sampler_test.cc online_stats_test.cc)
target_include_directories(lancet2_tests PRIVATE ${CMAKE_BINARY_DIR})
target_link_libraries(lancet2_tests PRIVATE lancet2_core Catch2 spdlog::spdlog)
set_target_properties(lancet2_tests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON
        INTERPROCEDURAL_OPTIMIZATION TRUE POSITION_INDEPENDENT_CODE 1)
