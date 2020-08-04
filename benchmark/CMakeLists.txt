add_executable(lancet_benchmark main.cpp)
target_set_warnings(lancet_benchmark ENABLE ALL AS_ERROR ALL DISABLE Annoying)
target_link_libraries(lancet_benchmark PRIVATE lancet_core benchmark)
set_target_properties(lancet_benchmark PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON
        INTERPROCEDURAL_OPTIMIZATION TRUE POSITION_INDEPENDENT_CODE 1)
