message(STATUS "Current CMake version: ${CMAKE_VERSION}")
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(lancet2 LANGUAGES C CXX
        HOMEPAGE_URL https://github.com/nygenome/Lancet2
        DESCRIPTION "Micro assembly based somatic variant caller")

option(LANCET2_TESTS "Build lancet2 tests" "OFF")
option(LANCET2_BENCHMARKS "Build lancet2 benchmarks" "OFF")
option(LANCET2_PROFILER "Build with gperftools CPU profiler" "OFF")
set(LANCET2_SANITIZER "OFF" CACHE STRING "Build with an LLVM sanitizer")
set_property(CACHE LANCET2_SANITIZER PROPERTY STRINGS "OFF" "Address" "Undefined" "Thread" "Memory")

# https://gitlab.kitware.com/cmake/cmake/-/issues/16716
set(PROJECT_VERSION "2.0.0.alpha1")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In source builds are not recommended. Please run cmake in a separate build directory")
endif ()

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to Release as none was specified.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED)
find_package(OpenSSL REQUIRED)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(defaults)
include(version)
include(dependencies)
include(coverage)
include(sanitizers)

# Add lancet2 core library target
add_library(lancet2_core STATIC
        include/lancet2/timer.h
        include/lancet2/sized_ints.h
        include/lancet2/assert_macro.h
        include/lancet2/log_macros.h
        include/lancet2/fractional_sampler.h
        include/lancet2/spinlock.h
        include/lancet2/utils.h src/utils.cc

        include/lancet2/contig_info.h
        include/lancet2/genomic_region.h
        include/lancet2/cigar.h src/cigar.cc
        include/lancet2/fasta_reader.h src/fasta_reader.cc
        include/lancet2/hts_alignment.h src/hts_alignment.cc
        include/lancet2/hts_reader.h src/hts_reader.cc

        include/lancet2/ref_window.h
        include/lancet2/read_info.h
        include/lancet2/window_builder.h src/window_builder.cc
        include/lancet2/cli_params.h src/cli_params.cc
        include/lancet2/core_enums.h src/core_enums.cc
        include/lancet2/read_extractor.h src/read_extractor.cc
        include/lancet2/fisher_exact.h src/fisher_exact.cc
        include/lancet2/variant.h src/variant.cc
        include/lancet2/variant_store.h src/variant_store.cc
        include/lancet2/tandem_repeat.h src/tandem_repeat.cc
        include/lancet2/canonical_kmers.h src/canonical_kmers.cc
        include/lancet2/dot_serializer.h src/dot_serializer.cc
        include/lancet2/align.h src/align.cc
        include/lancet2/edmond_karp.h src/edmond_karp.cc
        include/lancet2/transcript.h src/transcript.cc
        include/lancet2/variant_hpcov.h
        include/lancet2/path.h src/path.cc
        include/lancet2/path_builder.h src/path_builder.cc
        include/lancet2/edge.h src/edge.cc
        include/lancet2/merge_node_info.h src/merge_node_info.cc
        include/lancet2/barcode_set.h src/barcode_set.cc
        include/lancet2/node_cov.h src/node_cov.cc
        include/lancet2/base_label.h src/base_label.cc
        include/lancet2/node_label.h src/node_label.cc
        include/lancet2/online_stats.h
        include/lancet2/node_qual.h src/node_qual.cc
        include/lancet2/node_neighbour.h
        include/lancet2/kmer.h src/kmer.cc
        include/lancet2/node.h src/node.cc
        include/lancet2/graph.h src/graph.cc
        include/lancet2/graph_builder.h src/graph_builder.cc
        include/lancet2/micro_assembler.h src/micro_assembler.cc
        include/lancet2/run_pipeline.h src/run_pipeline.cc
        include/lancet2/cli.h src/cli.cc)

add_dependencies(lancet2_core htslib)
target_include_directories(lancet2_core PUBLIC include PRIVATE ${CMAKE_BINARY_DIR} ${HTSLIB_INCLUDE_DIR})
target_link_libraries(lancet2_core
        PUBLIC concurrentqueue absl::status absl::flat_hash_map absl::flat_hash_set
        absl::str_format absl::status absl::statusor absl::fixed_array absl::inlined_vector
        absl::span absl::time absl::random_random absl::strings spdlog::spdlog
        PRIVATE CLI11 absl::failure_signal_handler ${LIBHTS}
        INTERFACE ZLIB::ZLIB BZip2::BZip2 CURL::libcurl OpenSSL::Crypto ${LIBLZMA_LIBRARIES} ${LIBDEFLATE})
set_target_properties(lancet2_core PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON POSITION_INDEPENDENT_CODE 1)

# Add main lancet2 executable target
set(LANCET_MAIN_CPP "${CMAKE_BINARY_DIR}/generated/lancet_main.cc")
configure_file(src/main.cc.inc ${LANCET_MAIN_CPP} @ONLY)
add_executable(lancet2 ${LANCET_MAIN_CPP})
target_link_libraries(lancet2 PRIVATE mimalloc lancet2_core)
set_target_properties(lancet2 PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON POSITION_INDEPENDENT_CODE 1)

if (LANCET2_PROFILER)
    add_dependencies(lancet2 gperftools)
    target_link_libraries(lancet2 PRIVATE ${LIBPROFILER})
endif ()

if (LANCET2_TESTS)
    add_subdirectory(tests)
endif ()

if (LANCET2_BENCHMARKS)
    add_subdirectory(benchmarks)
endif ()
