set(LANCET_TEST_DATA_DIR "${PROJECT_SOURCE_DIR}/tests/data")
configure_file(test_config.h.in "${CMAKE_BINARY_DIR}/generated/test_config.h")

add_executable(lancet_test "${CMAKE_BINARY_DIR}/generated/test_config.h"
        lancet_test.cpp align_test.cpp)

target_include_directories(lancet_test PRIVATE ${CMAKE_BINARY_DIR})
target_set_warnings(lancet_test ENABLE ALL AS_ERROR ALL DISABLE Annoying)
target_link_libraries(lancet_test PRIVATE lancet_core Catch2 spdlog::spdlog)
set_target_properties(lancet_test PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON
        INTERPROCEDURAL_OPTIMIZATION TRUE POSITION_INDEPENDENT_CODE 1)
