message(STATUS "Current CMake version: ${CMAKE_VERSION}")
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

option(LANCET_WARNINGS "Set this to OFF to disable all warnings" "ON")
option(LANCET_UNIT_TESTS "Build snifty unit tests" "OFF")
option(LANCET_BENCHMARKS "Build snifty benchmarks" "OFF")
set(LANCET_SANITIZER "OFF" CACHE STRING "Build with an LLVM sanitizer")
set_property(CACHE LANCET_SANITIZER PROPERTY STRINGS "OFF" "Address" "Undefined" "Thread" "Memory")

project(lancet VERSION 1.0.0 LANGUAGES C CXX
        HOMEPAGE_URL https://github.com/omicsnut/v2_lancet
        DESCRIPTION "Microassembly based somatic variant caller")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In source builds are not recommended. Please run cmake in a separate build directory")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to Release as none was specified.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

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
include(warnings)
include(coverage)
include(sanitizers)

# Add lancet core library target
add_library(lancet_core STATIC
        include/lancet/logger.h src/logger.cpp
        include/lancet/assert_macro.h include/lancet/timer.h
        include/lancet/fractional_sampler.h
        include/lancet/statusor_internals.h
        include/lancet/statusor.h src/statusor.cpp
        include/lancet/utils.h src/utils.cpp

        include/lancet/contig_info.h
        include/lancet/genomic_region.h
        include/lancet/fasta_reader.h src/fasta_reader.cpp
        include/lancet/cigar.h src/cigar.cpp
        include/lancet/hts_alignment.h src/hts_alignment.cpp
        include/lancet/hts_reader.h src/hts_reader.cpp
        include/lancet/vcf_writer.h src/vcf_writer.cpp

        include/lancet/ref_window.h
        include/lancet/notification.h
        include/lancet/read_info.h
        include/lancet/window_builder.h src/window_builder.cpp
        include/lancet/cli_params.h src/cli_params.cpp
        include/lancet/core_enums.h src/core_enums.cpp
        include/lancet/read_extractor.h src/read_extractor.cpp
        include/lancet/fisher_exact.h src/fisher_exact.cpp
        include/lancet/variant.h src/variant.cpp
        include/lancet/variant_store.h src/variant_store.cpp
        include/lancet/tandem_repeat.h src/tandem_repeat.cpp
        include/lancet/canonical_kmers.h src/canonical_kmers.cpp
        include/lancet/dot_serializer.h src/dot_serializer.cpp
        include/lancet/align.h src/align.cpp
        include/lancet/edmond_karp.h src/edmond_karp.cpp
        include/lancet/transcript.h src/transcript.cpp
        include/lancet/sample_cov.h src/sample_cov.cpp
        include/lancet/variant_hpcov.h
        include/lancet/path.h src/path.cpp
        include/lancet/path_builder.h src/path_builder.cpp
        include/lancet/edge.h src/edge.cpp
        include/lancet/merge_node_info.h src/merge_node_info.cpp
        include/lancet/base_cov.h include/lancet/base_hp.h
        include/lancet/base_hpcov.h src/base_hpcov.cpp
        include/lancet/barcode_set.h src/barcode_set.cpp
        include/lancet/node_cov.h src/node_cov.cpp
        include/lancet/node_hp.h src/node_hp.cpp
        include/lancet/base_label.h src/base_label.cpp
        include/lancet/node_label.h src/node_label.cpp
        include/lancet/cov_stats.h src/cov_stats.cpp
        include/lancet/online_stats.h
        include/lancet/node_qual.h src/node_qual.cpp
        include/lancet/node_neighbour.h
        include/lancet/kmer.h src/kmer.cpp
        include/lancet/node.h src/node.cpp
        include/lancet/graph.h src/graph.cpp
        include/lancet/graph_builder.h src/graph_builder.cpp
        include/lancet/micro_assembler.h src/micro_assembler.cpp
        include/lancet/run_pipeline.h src/run_pipeline.cpp
        include/lancet/cli.h src/cli.cpp)

add_dependencies(lancet_core htslib)
target_set_warnings(lancet_core ENABLE ALL AS_ERROR ALL DISABLE Annoying)
target_include_directories(lancet_core PUBLIC include PRIVATE ${CMAKE_BINARY_DIR} ${HTSLIB_INCLUDE_DIR})

target_link_libraries(lancet_core
        PUBLIC absl::status absl::flat_hash_map absl::flat_hash_set concurrentqueue absl::str_format
               absl::fixed_array absl::span absl::time absl::random_random absl::synchronization
        PRIVATE CLI11 absl::failure_signal_handler absl::strings ${LIBHTS}
        INTERFACE ZLIB::ZLIB BZip2::BZip2 CURL::libcurl OpenSSL::Crypto ${LIBLZMA_LIBRARIES})

set_target_properties(lancet_core PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON POSITION_INDEPENDENT_CODE 1)

# Add main lancet executable target
set(LANCET_MAIN_CPP "${CMAKE_BINARY_DIR}/generated/lancet_main.cpp")
configure_file(main.cpp.in ${LANCET_MAIN_CPP} @ONLY)
add_executable(lancet ${LANCET_MAIN_CPP})
target_set_warnings(lancet ENABLE ALL AS_ERROR ALL DISABLE Annoying)
target_link_libraries(lancet PRIVATE mimalloc lancet_core)
set_target_properties(lancet PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON POSITION_INDEPENDENT_CODE 1)

if(LANCET_UNIT_TESTS)
    add_subdirectory(tests)
endif()

if(LANCET_BENCHMARKS)
    add_subdirectory(benchmark)
endif()
