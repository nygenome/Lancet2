set(LANCET_FULL_DATA_DIR "${PROJECT_SOURCE_DIR}/data")
set(LANCET_BENCHMARK_CONFIG_H "${CMAKE_BINARY_DIR}/generated/lancet_benchmark_config.h")
configure_file(benchmark_config.h.inc ${LANCET_BENCHMARK_CONFIG_H} @ONLY)

add_executable(BenchmarkLancet2 main.cpp extractor_bench.cpp hamming_bench.cpp)
target_include_directories(BenchmarkLancet2 PRIVATE "${CMAKE_BINARY_DIR}/generated" "${CMAKE_SOURCE_DIR}")
target_link_libraries(BenchmarkLancet2 PRIVATE mimalloc-static benchmark lancet_cli)
set_target_properties(BenchmarkLancet2 PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
set_property(TARGET BenchmarkLancet2 PROPERTY $<$<ENABLE_LTO>:INTERPROCEDURAL_OPTIMIZATION TRUE>)
