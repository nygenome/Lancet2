include(ExternalProject)
include(FetchContent)
include(ProcessorCount)
ProcessorCount(NumCores)
find_program(MAKE_EXE NAMES gmake nmake make)

# Overwrite the message() function to prevent subproject deps to write messages
function(message)
    if (NOT MESSAGE_QUIET)
        _message(${ARGN})
    endif ()
endfunction()

set(LIBDEFLATE_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/libdeflate")
set(LIBDEFLATE "${LIBDEFLATE_ROOT_DIR}/libdeflate.a")
ExternalProject_Add(libdeflate
        URL https://github.com/ebiggers/libdeflate/archive/refs/tags/v1.17.tar.gz
        URL_MD5 88cecda43f00fcffb62fd2a398f7554e
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${LIBDEFLATE_ROOT_DIR}
        BUILD_IN_SOURCE 1 INSTALL_COMMAND "" CONFIGURE_COMMAND "" BUILD_BYPRODUCTS ${LIBDEFLATE}
        BUILD_COMMAND ${MAKE_EXE} -j${NumCores}
        LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON USES_TERMINAL_DOWNLOAD OFF
        USER_TERMINAL_CONFIGURE OFF USES_TERMINAL_BUILD OFF)

set(HTSLIB_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/htslib")
set(HTSLIB_INCLUDE_DIR "${HTSLIB_ROOT_DIR}")
set(LIBHTS "${HTSLIB_ROOT_DIR}/libhts.a")
ExternalProject_Add(htslib
        URL https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
        URL_MD5 ff47f4b09e5202cebc3f80e7e02c7728
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${HTSLIB_ROOT_DIR}
        BUILD_IN_SOURCE 1 INSTALL_COMMAND "" BUILD_COMMAND ${MAKE_EXE} -j${NumCores} lib-static
        CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/scripts/configure_htslib.sh ${HTSLIB_ROOT_DIR} ${CMAKE_C_COMPILER}
        BUILD_BYPRODUCTS ${LIBHTS} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
        USES_TERMINAL_DOWNLOAD OFF USER_TERMINAL_CONFIGURE OFF USES_TERMINAL_BUILD OFF)
add_dependencies(htslib libdeflate)

set(GPERFTOOLS_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/gperftools")
set(LIBPROFILER "${GPERFTOOLS_ROOT_DIR}/lib/libprofiler${CMAKE_SHARED_LIBRARY_SUFFIX}")
ExternalProject_Add(gperftools
        URL https://github.com/gperftools/gperftools/releases/download/gperftools-2.10/gperftools-2.10.tar.gz
        URL_MD5 62bf6c76ba855ed580de5e139bd2a483
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${GPERFTOOLS_ROOT_DIR}
        BUILD_IN_SOURCE 1 INSTALL_COMMAND ${MAKE_EXE} install BUILD_COMMAND ${MAKE_EXE} -j${NumCores}
        CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/scripts/configure_gperftools.sh ${GPERFTOOLS_ROOT_DIR} ${CMAKE_C_COMPILER} ${CMAKE_CXX_COMPILER}
        BUILD_BYPRODUCTS ${LIBPROFILER} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
        USES_TERMINAL_DOWNLOAD OFF USER_TERMINAL_CONFIGURE OFF USES_TERMINAL_BUILD OFF)

set(MESSAGE_QUIET ON)

FetchContent_Declare(concurrentqueue GIT_REPOSITORY https://github.com/cameron314/concurrentqueue.git GIT_TAG v1.0.3)
FetchContent_GetProperties(concurrentqueue)
if (NOT concurrentqueue_POPULATED)
    FetchContent_Populate(concurrentqueue)
    add_library(concurrentqueue INTERFACE)
    target_include_directories(concurrentqueue SYSTEM INTERFACE "${concurrentqueue_SOURCE_DIR}")
endif ()

FetchContent_Declare(spdlog GIT_REPOSITORY https://github.com/gabime/spdlog.git GIT_TAG v1.11.0)
FetchContent_MakeAvailable(spdlog)

FetchContent_Declare(cli11 GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git GIT_TAG v2.3.2)
FetchContent_MakeAvailable(cli11)

set(MI_BUILD_TESTS OFF)
FetchContent_Declare(mimalloc GIT_REPOSITORY https://github.com/microsoft/mimalloc.git GIT_TAG v2.0.9)
FetchContent_MakeAvailable(mimalloc)

FetchContent_Declare(abseil GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git GIT_TAG 20230125.1)
FetchContent_GetProperties(abseil)
if (NOT abseil_POPULATED)
    set(BUILD_TESTING OFF)
    FetchContent_Populate(abseil)
    add_subdirectory(${abseil_SOURCE_DIR} ${abseil_BINARY_DIR})
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${abseil_SOURCE_DIR}/absl/copts)
    include(${abseil_SOURCE_DIR}/absl/copts/AbseilConfigureCopts.cmake)
endif ()

if (LANCET2_TESTS)
    FetchContent_Declare(catch GIT_REPOSITORY https://github.com/catchorg/Catch2.git GIT_TAG v3.3.2)
    FetchContent_MakeAvailable(catch)
endif ()

if (LANCET2_BENCHMARKS)
    set(BENCHMARK_ENABLE_TESTING OFF)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
    set(BENCHMARK_ENABLE_ASSEMBLY_TESTS OFF)
    set(BENCHMARK_ENABLE_INSTALL OFF)
    FetchContent_Declare(benchmark GIT_REPOSITORY https://github.com/google/benchmark.git GIT_TAG v1.7.1)
    FetchContent_MakeAvailable(benchmark)
endif ()

unset(MESSAGE_QUIET)
