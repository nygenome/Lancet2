set(LANCET_CAN_RUN_COVERAGE "OFF" CACHE INTERNAL "")
find_program(GCOV_EXE NAMES gcov)

if(NOT GCOV_EXE)
    message(STATUS "gcov not found in PATH. Skipping code coverage analysis.")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    message(STATUS "Adding coverage flags for ${CMAKE_CXX_COMPILER_ID}-${CMAKE_BUILD_TYPE} build")
    set(COVERAGE_COMPILER_FLAGS "-g -O1 -fprofile-arcs -ftest-coverage" CACHE INTERNAL "")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(LANCET_CAN_RUN_COVERAGE "ON" CACHE INTERNAL "")
endif()

mark_as_advanced(FORCE LANCET_CAN_RUN_COVERAGE GCOV_EXE LCOV_EXE COVERAGE_COMPILER_FLAGS)
