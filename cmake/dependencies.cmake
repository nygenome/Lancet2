include(ExternalProject)
include(FetchContent)
mark_as_advanced(FORCE FETCHCONTENT_BASE_DIR FETCHCONTENT_QUIET
        FETCHCONTENT_FULLY_DISCONNECTED FETCHCONTENT_UPDATES_DISCONNECTED)

include(ProcessorCount)
ProcessorCount(NumCores)

find_program(MAKE_EXE NAMES gmake nmake make)
mark_as_advanced(FORCE MAKE_EXE)

# Overwrite the message() function to prevent subproject deps to write messages
function(message)
    if (NOT MESSAGE_QUIET)
        _message(${ARGN})
    endif()
endfunction()

set(HTSLIB_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/htslib")
set(HTSLIB_INCLUDE_DIR "${HTSLIB_ROOT_DIR}")
set(LIBHTS "${HTSLIB_ROOT_DIR}/libhts.a")

ExternalProject_Add(htslib
        URL https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2
        URL_MD5 c488c7a79283e5252c04cae335ee80e8
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${HTSLIB_ROOT_DIR}
        BUILD_IN_SOURCE 1 INSTALL_COMMAND "" BUILD_COMMAND ${MAKE_EXE} -j${NumCores} lib-static
        CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/scripts/configure_htslib.sh ${HTSLIB_ROOT_DIR} ${CMAKE_C_COMPILER}
        BUILD_BYPRODUCTS ${LIBHTS} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
        USES_TERMINAL_DOWNLOAD OFF USER_TERMINAL_CONFIGURE OFF USES_TERMINAL_BUILD OFF)

set(GPERFTOOLS_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/gperftools")
set(LIBPROFILER "${GPERFTOOLS_ROOT_DIR}/lib/libprofiler${CMAKE_SHARED_LIBRARY_SUFFIX}")

ExternalProject_Add(gperftools
        URL https://github.com/gperftools/gperftools/releases/download/gperftools-2.8/gperftools-2.8.tar.gz
        URL_MD5 ea61ad1e886f53e4069873ce52cae7df
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${GPERFTOOLS_ROOT_DIR}
        BUILD_IN_SOURCE 1 INSTALL_COMMAND ${MAKE_EXE} install BUILD_COMMAND ${MAKE_EXE} -j${NumCores}
        CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/scripts/configure_gperftools.sh ${GPERFTOOLS_ROOT_DIR} ${CMAKE_C_COMPILER} ${CMAKE_CXX_COMPILER}
        BUILD_BYPRODUCTS ${LIBPROFILER} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON
        USES_TERMINAL_DOWNLOAD OFF USER_TERMINAL_CONFIGURE OFF USES_TERMINAL_BUILD OFF)

set(MESSAGE_QUIET ON)

#FetchContent_Declare(gzipxx GIT_REPOSITORY https://github.com/mapbox/gzip-hpp.git GIT_TAG v0.1.0)
#FetchContent_GetProperties(gzipxx)
#if(NOT gzipxx_POPULATED)
#    FetchContent_Populate(gzipxx)
#    add_library(gzipxx INTERFACE)
#    target_include_directories(gzipxx SYSTEM INTERFACE ${gzipxx_SOURCE_DIR}/include)
#endif()

#set(FLATBUFFERS_BUILD_TESTS OFF)
#set(FLATBUFFERS_MAX_PARSING_DEPTH 16)
#FetchContent_Declare(flatbuffers GIT_REPOSITORY https://github.com/google/flatbuffers.git GIT_TAG v1.12.0)
#FetchContent_MakeAvailable(flatbuffers)

FetchContent_Declare(concurrentqueue GIT_REPOSITORY https://github.com/cameron314/concurrentqueue.git GIT_TAG v1.0.1)
FetchContent_GetProperties(concurrentqueue)
if(NOT concurrentqueue_POPULATED)
    FetchContent_Populate(concurrentqueue)
    add_library(concurrentqueue INTERFACE)
    target_include_directories(concurrentqueue SYSTEM INTERFACE "${concurrentqueue_SOURCE_DIR}")
endif()

FetchContent_Declare(spdlog GIT_REPOSITORY https://github.com/gabime/spdlog.git GIT_TAG v1.8.1)
FetchContent_MakeAvailable(spdlog)

FetchContent_Declare(cli11 GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git GIT_TAG v1.9.1)
FetchContent_MakeAvailable(cli11)

set(MI_BUILD_TESTS OFF)
FetchContent_Declare(mimalloc GIT_REPOSITORY https://github.com/microsoft/mimalloc.git GIT_TAG v1.6.7)
FetchContent_MakeAvailable(mimalloc)

# Add abseil for general utilities. Update as often as possible.
# See https://abseil.io/about/philosophy#upgrade-support
FetchContent_Declare(abseil GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
                            GIT_TAG e63a5a61045e516b7b3dbca090e2b9ff1d057e46)
FetchContent_GetProperties(abseil)
if(NOT abseil_POPULATED)
    set(BUILD_TESTING OFF)
    FetchContent_Populate(abseil)
    add_subdirectory(${abseil_SOURCE_DIR} ${abseil_BINARY_DIR})
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${abseil_SOURCE_DIR}/absl/copts)
    include(${abseil_SOURCE_DIR}/absl/copts/AbseilConfigureCopts.cmake)
endif()

if(LANCET_UNIT_TESTS)
    FetchContent_Declare(catch GIT_REPOSITORY https://github.com/catchorg/Catch2.git GIT_TAG v2.13.2)
    FetchContent_MakeAvailable(catch)
endif()


if(LANCET_BENCHMARKS)
    set(BENCHMARK_ENABLE_TESTING OFF)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
    set(BENCHMARK_ENABLE_ASSEMBLY_TESTS OFF)
    set(BENCHMARK_ENABLE_INSTALL OFF)
    FetchContent_Declare(benchmark GIT_REPOSITORY https://github.com/google/benchmark.git GIT_TAG v1.5.2)
    FetchContent_MakeAvailable(benchmark)
endif()

unset(MESSAGE_QUIET)
