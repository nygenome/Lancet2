message(STATUS "Current CMake version: ${CMAKE_VERSION}")
cmake_minimum_required(VERSION 3.25 FATAL_ERROR)

project(lancet LANGUAGES C CXX
		HOMEPAGE_URL https://github.com/nygenome/Lancet2
		DESCRIPTION "Micro-assembly based somatic variant caller")

# https://gitlab.kitware.com/cmake/cmake/-/issues/16716
include(cmake/version.cmake)
set(PROJECT_VERSION ${LANCET2_GIT_TAG})
message(STATUS "Lancet v${PROJECT_VERSION}")
option(LANCET_TESTS "Builds the Lancet tests executable" "OFF")
option(LANCET_BENCHMARKS "Builds the Lancet benchmarks executable" "OFF")
option(LANCET_BUILD_STATIC "Build a statically linked Lancet executable" "ON")
option(LANCET_BUILD_ARCH "Build the Lancet executable for given CPU architecture")
set_property(CACHE LANCET_TESTS PROPERTY STRINGS "OFF" "ON")
set_property(CACHE LANCET_BENCHMARKS PROPERTY STRINGS "OFF" "ON")
set_property(CACHE LANCET_BUILD_STATIC PROPERTY STRINGS "OFF" "ON")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
message(STATUS "CMAKE_HOST_SYSTEM_NAME: ${CMAKE_HOST_SYSTEM_NAME}")
message(STATUS "CMAKE_HOST_SYSTEM_PROCESSOR: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
MESSAGE(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
MESSAGE(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

if (NOT CMAKE_HOST_SYSTEM_NAME MATCHES Linux)
	message(WARNING "Running Lancet build on a ${CMAKE_HOST_SYSTEM_NAME} machine")
	message(FATAL_ERROR "Lancet can only be built on systems with Linux Operating System")
endif ()

if (NOT (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES x86_64 OR CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL amd64))
	message(WARNING "Running Lancet build on a ${CMAKE_HOST_SYSTEM_PROCESSOR} machine")
	message(FATAL_ERROR "Lancet can only be built on a 64-bit x86 machine")
endif ()

if (NOT ${CMAKE_C_COMPILER_ID} STREQUAL "GNU" OR NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	message(WARNING "Running Lancet build with ${CMAKE_C_COMPILER} and ${CMAKE_CXX_COMPILER}")
	message(FATAL_ERROR "Lancet can only be built from source with GCC compiler")
endif ()

if (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "12.0")
	message(WARNING "Running Lancet build with GCC version ${CMAKE_CXX_COMPILER_VERSION}")
	message(FATAL_ERROR "Lancet requires GCC version to be atleast 12.0 or greater")
endif ()

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting CMAKE_BUILD_TYPE as Release since none was specified.")
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)

include(cmake/cpu_family.cmake)
include(cmake/defaults.cmake)
include(cmake/dependencies.cmake)

add_library(lancet_base STATIC
		src/lancet/base/version.h src/lancet/base/timer.h
		src/lancet/base/types.h src/lancet/base/assert.h
		src/lancet/base/logging.h src/lancet/base/rev_comp.h
		src/lancet/base/compute_stats.h src/lancet/base/sliding.h
		src/lancet/base/hash.cpp src/lancet/base/hash.h
		src/lancet/base/repeat.cpp src/lancet/base/repeat.h
		src/lancet/base/find_str.cpp src/lancet/base/find_str.h)
target_link_libraries(lancet_base PRIVATE wyhash MurmurHash3 absl::flat_hash_set
		PUBLIC spdlog::spdlog absl::span absl::fixed_array absl::strings absl::time)
target_include_directories(lancet_base PUBLIC "${CMAKE_SOURCE_DIR}/src" "${CMAKE_BINARY_DIR}/generated")
set_target_properties(lancet_base PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

add_library(lancet_hts STATIC
		src/lancet/hts/bgzf_ostream.cpp src/lancet/hts/bgzf_ostream.h
		src/lancet/hts/phred_quality.cpp src/lancet/hts/phred_quality.h
		src/lancet/hts/fisher_exact.cpp src/lancet/hts/fisher_exact.h
		src/lancet/hts/reference.cpp src/lancet/hts/reference.h
		src/lancet/hts/cigar_unit.h src/lancet/hts/aux_tag.h
		src/lancet/hts/alignment.cpp src/lancet/hts/alignment.h
		src/lancet/hts/iterator.cpp src/lancet/hts/iterator.h
		src/lancet/hts/extractor.cpp src/lancet/hts/extractor.h)
add_dependencies(lancet_hts htslib)
set_target_properties(lancet_hts PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
target_include_directories(lancet_hts SYSTEM PUBLIC ${HTSLIB_ROOT_DIR})
target_include_directories(lancet_hts PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(lancet_hts PUBLIC lancet_base absl::statusor absl::fixed_array absl::flat_hash_set ${LIB_HTS}
		PRIVATE absl::flat_hash_map INTERFACE BZip2::BZip2 LibLZMA::LibLZMA zlibstatic libdeflate_static)

add_library(lancet_cbdg STATIC src/lancet/cbdg/label.h
		src/lancet/cbdg/kmer.cpp src/lancet/cbdg/kmer.h
		src/lancet/cbdg/edge.h src/lancet/cbdg/read.h
		src/lancet/cbdg/node.cpp src/lancet/cbdg/node.h
		src/lancet/cbdg/max_flow.cpp src/lancet/cbdg/max_flow.h
		src/lancet/cbdg/graph.cpp src/lancet/cbdg/graph.h)
target_include_directories(lancet_cbdg PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(lancet_cbdg PUBLIC lancet_hts absl::cord absl::inlined_vector absl::flat_hash_map)
set_target_properties(lancet_cbdg PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

add_library(lancet_caller STATIC src/lancet/caller/raw_variant.h
		src/lancet/caller/variant_support.cpp src/lancet/caller/variant_support.h
		src/lancet/caller/variant_call.cpp src/lancet/caller/variant_call.h
		src/lancet/caller/msa_builder.cpp src/lancet/caller/msa_builder.h
		src/lancet/caller/variant_set.cpp src/lancet/caller/variant_set.h
		src/lancet/caller/genotyper.cpp src/lancet/caller/genotyper.h)
add_dependencies(lancet_caller minimap2)
target_include_directories(lancet_caller PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_include_directories(lancet_caller SYSTEM PUBLIC "${MM2_ROOT_DIR}")
target_link_libraries(lancet_caller PUBLIC lancet_cbdg spoa absl::btree ${LIB_MM2} PRIVATE boost_math zlibstatic)
set_target_properties(lancet_caller PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

add_library(lancet_core STATIC src/lancet/core/window.h src/lancet/core/sample_info.h
		src/lancet/core/window_builder.cpp src/lancet/core/window_builder.h
		src/lancet/core/read_collector.cpp src/lancet/core/read_collector.h
		src/lancet/core/variant_store.cpp src/lancet/core/variant_store.h
		src/lancet/core/variant_builder.cpp src/lancet/core/variant_builder.h
		src/lancet/core/async_worker.cpp src/lancet/core/async_worker.h)
target_include_directories(lancet_core PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(lancet_core PUBLIC lancet_caller absl::synchronization concurrentqueue PRIVATE absl::hash)
set_target_properties(lancet_core PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

add_library(lancet_cli STATIC src/lancet/cli/cli_params.h
		src/lancet/cli/eta_timer.cpp src/lancet/cli/eta_timer.h
		src/lancet/cli/pipeline_runner.cpp src/lancet/cli/pipeline_runner.h
		src/lancet/cli/cli_interface.cpp src/lancet/cli/cli_interface.h)
target_include_directories(lancet_cli PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(lancet_cli PUBLIC lancet_core CLI11)
set_target_properties(lancet_cli PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

add_executable(Lancet2 src/lancet/main.cpp)
target_include_directories(Lancet2 PRIVATE "${CMAKE_SOURCE_DIR}")
target_link_libraries(Lancet2 PRIVATE lancet_cli absl::symbolize
		absl::failure_signal_handler absl::cleanup mimalloc-static)
set_target_properties(Lancet2 PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
set_property(TARGET Lancet2 PROPERTY $<$<ENABLE_LTO>:INTERPROCEDURAL_OPTIMIZATION TRUE>)

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
	add_compile_definitions(LANCET_DEVELOP_MODE)
endif ()

# Helpful links for profiling – https://github.com/google/pprof
# https://gperftools.github.io/gperftools/cpuprofile.html
#if (${CMAKE_BUILD_TYPE} MATCHES Release)
#	add_dependencies(lancet_cli gperftools)
#	target_include_directories(lancet_cli SYSTEM PRIVATE ${GPERFTOOLS_INC_DIR})
#	target_link_libraries(lancet_cli PRIVATE ${LIB_PROFILER})
#endif ()

# Helpful links to run tests
# https://github.com/catchorg/Catch2/blob/v3.3.2/docs/Readme.md
if (LANCET_TESTS)
	add_subdirectory(tests)
endif ()

if (LANCET_BENCHMARKS)
	add_subdirectory(benchmarks)
endif ()
