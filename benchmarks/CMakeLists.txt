add_executable(lancet2_benchmarks main.cc align_bench.cc)
target_link_libraries(lancet2_benchmarks PRIVATE lancet2_core benchmark)
set_target_properties(lancet2_benchmarks PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON
        INTERPROCEDURAL_OPTIMIZATION TRUE POSITION_INDEPENDENT_CODE 1)
