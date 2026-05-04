# ═══════════════════════════════════════════════════════════════════════════════
# Third-Party Dependencies — cmake policies, vendored libraries
#
# All dependencies are pinned to exact versions for reproducible builds.
# FetchContent downloads at configure time; ExternalProject builds at build time.
#
# Core libraries:
#   mimalloc      — high-performance allocator (replaces system malloc)
#   abseil-cpp    — Abseil C++ utilities (containers, strings, hashing)
#   spdlog        — structured logging (bundles fmtlib)
#   CLI11         — command-line argument parsing
#   concurrentqueue — lock-free multi-producer/consumer queue
#
# Compression (required by HTSlib):
#   libdeflate    — fast DEFLATE/gzip compression
#   zlib-ng       — zlib-compatible compression (SIMD-optimized)
#
# Bioinformatics:
#   htslib        — BAM/CRAM/VCF I/O (ExternalProject, builds libhts.a)
#   minimap2      — read-to-haplotype alignment (ExternalProject, builds libminimap2.a)
#   spoa          — SIMD Partial Order Alignment (graph-based MSA)
#
# Testing / Benchmarking / Profiling:
#   Catch2        — unit test framework (amalgamated, only when -DLANCET_TESTS=ON)
#   longdust      — C reference implementation for LongdustQScorer cross-validation
#   benchmark     — Google Benchmark (only when -DLANCET_BENCHMARKS=ON)
#   gperftools    — CPU profiler (only when -DLANCET_PROFILE_MODE=ON)
# ═══════════════════════════════════════════════════════════════════════════════

# ── Policy Defaults for FetchContent Subprojects ──────────────────────────────
# With cmake_minimum_required(VERSION 3.25), all these policies are already NEW
# in Lancet2's own scope. These CMAKE_POLICY_DEFAULT_* variables force NEW
# behavior on FetchContent-downloaded subprojects that may declare a lower
# cmake_minimum_required (e.g., abseil, spoa, zlib-ng).
#
# CMP0012: if() recognizes numbers and boolean constants
# CMP0048: project() supports VERSION keyword
# CMP0063: Honor visibility properties for all target types
# CMP0069: Enforce INTERPROCEDURAL_OPTIMIZATION (LTO)
# CMP0074: find_package uses <PackageName>_ROOT variables
# CMP0077: option() honors normal variables set before it
# CMP0135: FetchContent sets file timestamps to extraction time
# CMP0067: Honor language standard in try_compile checks
set(CMAKE_POLICY_DEFAULT_CMP0012 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0063 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0074 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0135 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0067 NEW)

# ── Vendored Dependencies ────────────────────────────────────────────────────
include(ExternalProject)
include(FetchContent)
include(ProcessorCount)
ProcessorCount(NumCores)
find_program(MAKE_EXE NAMES gmake nmake make REQUIRED)

# ── Third-Party Output Suppression ───────────────────────────────────────────
# When LANCET_QUIET_DEPS is ON (default), third-party STATUS messages are
# suppressed. lancet_dep_status() briefly lifts the suppression to emit
# Lancet2's own progress messages. Pass -DLANCET_QUIET_DEPS=OFF to cmake
# for full unfiltered output when debugging dependency issues.
if(LANCET_QUIET_DEPS)
	set(_LANCET_DEP_LOG_LEVEL WARNING)
else()
	set(_LANCET_DEP_LOG_LEVEL "")
endif()

set(CMAKE_MESSAGE_LOG_LEVEL "${_LANCET_DEP_LOG_LEVEL}")
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ${LANCET_QUIET_DEPS})

macro(lancet_dep_status msg)
	set(CMAKE_MESSAGE_LOG_LEVEL "")
	message(STATUS "${msg}")
	set(CMAKE_MESSAGE_LOG_LEVEL "${_LANCET_DEP_LOG_LEVEL}")
endmacro()

set(MI_SECURE OFF)
set(MI_PADDING OFF)
set(MI_BUILD_STATIC ON)
set(MI_BUILD_SHARED OFF)
set(MI_BUILD_OBJECT OFF)
set(MI_BUILD_TESTS OFF)
set(MI_OVERRIDE ON)
if (CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
	set(MI_OVERRIDE OFF)
endif ()
FetchContent_Declare(mimalloc GIT_REPOSITORY https://github.com/microsoft/mimalloc.git GIT_TAG v3.3.2 SYSTEM)
lancet_dep_status("Configuring mimalloc v3.3.2 ...")
FetchContent_MakeAvailable(mimalloc)

# abseil-cpp deaf349 - Committed on May 1st 2026
lancet_dep_status("Configuring abseil-cpp deaf349 ...")
FetchContent_Declare(abseil GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git GIT_TAG deaf349 SYSTEM)
FetchContent_GetProperties(abseil)
if (NOT abseil_POPULATED)
	set(BUILD_TESTING OFF)
	set(ABSL_PROPAGATE_CXX_STD ON)
	set(ABSL_USE_SYSTEM_INCLUDES ON)
	FetchContent_Populate(abseil)
	add_subdirectory(${abseil_SOURCE_DIR} ${abseil_BINARY_DIR} SYSTEM)
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${abseil_SOURCE_DIR}/absl/copts)
	include(${abseil_SOURCE_DIR}/absl/copts/AbseilConfigureCopts.cmake)
endif ()

lancet_dep_status("Configuring spdlog v1.17.0 ...")
FetchContent_Declare(spdlog GIT_REPOSITORY https://github.com/gabime/spdlog.git GIT_TAG v1.17.0 SYSTEM)
FetchContent_MakeAvailable(spdlog)

lancet_dep_status("Configuring CLI11 v2.6.2 ...")
FetchContent_Declare(cli11 GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git GIT_TAG v2.6.2 SYSTEM)
FetchContent_MakeAvailable(cli11)

lancet_dep_status("Configuring concurrentqueue v1.0.5 ...")
FetchContent_Declare(concurrentqueue GIT_REPOSITORY https://github.com/cameron314/concurrentqueue.git GIT_TAG v1.0.5 SYSTEM)
FetchContent_GetProperties(concurrentqueue)
if (NOT concurrentqueue_POPULATED)
	FetchContent_Populate(concurrentqueue)
	add_library(concurrentqueue INTERFACE)
	target_include_directories(concurrentqueue SYSTEM INTERFACE "${concurrentqueue_SOURCE_DIR}")
endif ()

set(LIBDEFLATE_BUILD_STATIC_LIB ON)
set(LIBDEFLATE_BUILD_SHARED_LIB OFF)
set(LIBDEFLATE_BUILD_GZIP OFF)
set(LIBDEFLATE_BUILD_TESTS OFF)
set(LIBDEFLATE_USE_SHARED_LIB OFF)
lancet_dep_status("Configuring libdeflate v1.25 ...")
FetchContent_Declare(libdeflate GIT_REPOSITORY https://github.com/ebiggers/libdeflate.git GIT_TAG v1.25 SYSTEM)
FetchContent_MakeAvailable(libdeflate)

set(ZLIB_COMPAT ON)
set(ZLIBNG_ENABLE_TESTS OFF)
set(WITH_GTEST OFF)
set(WITH_BENCHMARKS OFF)
set(BUILD_SHARED_LIBS OFF)
set(WITH_NEW_STRATEGIES ON)
set(WITH_OPTIM ON)
set(WITH_NATIVE_INSTRUCTIONS ${LANCET_NATIVE_BUILD})
lancet_dep_status("Configuring zlib-ng 2.3.3 ...")
FetchContent_Declare(zlib-ng GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git GIT_TAG 2.3.3 SYSTEM)
FetchContent_MakeAvailable(zlib-ng)

set(HTSLIB_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/htslib")
set(LIB_HTS "${HTSLIB_ROOT_DIR}/libhts.a")
set(HTSLIB_CONFIG_PARAMS ${HTSLIB_ROOT_DIR} ${CMAKE_C_COMPILER} ${LANCET_ENABLE_CLOUD_IO} ${CMAKE_INSTALL_PREFIX})
lancet_dep_status("Configuring htslib 1.23.1 ...")
ExternalProject_Add(htslib
		URL https://github.com/samtools/htslib/releases/download/1.23.1/htslib-1.23.1.tar.bz2
		URL_MD5 aee2c757fd8c88b9b5b61e8a1eae99de PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
		SOURCE_DIR ${HTSLIB_ROOT_DIR} BUILD_IN_SOURCE 1 INSTALL_COMMAND ""
		BUILD_COMMAND ${MAKE_EXE} -j${NumCores} lib-static BUILD_BYPRODUCTS ${LIB_HTS}
		CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/configure_htslib.sh ${HTSLIB_CONFIG_PARAMS}
		LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON USES_TERMINAL_DOWNLOAD OFF
		USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)
add_dependencies(htslib zlibstatic libdeflate_static)

set(MM2_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/minimap2")
set(LIB_MM2 "${MM2_ROOT_DIR}/libminimap2.a")
set(MM2_BUILD_PARAMS ${MM2_ROOT_DIR} ${CMAKE_C_COMPILER} ${CMAKE_HOST_SYSTEM_PROCESSOR})
lancet_dep_status("Configuring minimap2 2.30 ...")
ExternalProject_Add(minimap2
		URL https://github.com/lh3/minimap2/releases/download/v2.30/minimap2-2.30.tar.bz2
		URL_MD5 e016e3578bf6c763cefe08e7f22f440c PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
		SOURCE_DIR ${MM2_ROOT_DIR} BUILD_IN_SOURCE 1 CONFIGURE_COMMAND "" INSTALL_COMMAND ""
		BUILD_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/build_minimap2.sh ${MM2_BUILD_PARAMS}
		BUILD_BYPRODUCTS ${LIB_MM2} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
		USES_TERMINAL_DOWNLOAD OFF USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)
add_dependencies(minimap2 zlibstatic)

set(spoa_optimize_for_native ${LANCET_NATIVE_BUILD})
lancet_dep_status("Configuring spoa 4.1.5 ...")
FetchContent_Declare(spoa GIT_REPOSITORY https://github.com/rvaser/spoa GIT_TAG 4.1.5 SYSTEM)
FetchContent_MakeAvailable(spoa)

if (LANCET_PROFILE_MODE)
	lancet_dep_status("Configuring gperftools 2.18.1 ...")
	set(GPERFTOOLS_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/gperftools")
	set(GPERFTOOLS_INC_DIR "${GPERFTOOLS_ROOT_DIR}/include")
	set(LIB_PROFILER "${GPERFTOOLS_ROOT_DIR}/lib/libprofiler.a")
	set(GPERFTOOLS_CONFIG_PARAMS ${GPERFTOOLS_ROOT_DIR} ${CMAKE_C_COMPILER} ${CMAKE_CXX_COMPILER})
	ExternalProject_Add(gperftools
		URL https://github.com/gperftools/gperftools/releases/download/gperftools-2.18.1/gperftools-2.18.1.tar.gz
		URL_MD5 129c01f6f5297f0482b33d431b5ec555
		PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${GPERFTOOLS_ROOT_DIR} BUILD_IN_SOURCE 1
		INSTALL_COMMAND ${MAKE_EXE} install BUILD_COMMAND ${MAKE_EXE} -j${NumCores}
		CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/configure_gperftools.sh ${GPERFTOOLS_CONFIG_PARAMS}
		BUILD_BYPRODUCTS ${LIB_PROFILER} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
		USES_TERMINAL_DOWNLOAD OFF USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)
endif ()

if (LANCET_TESTS)
	lancet_dep_status("Configuring Catch2 v3.14.0 ...")
	file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_deps/Catch2")
	set(CATCH_ROOT "${CMAKE_CURRENT_BINARY_DIR}/_deps/Catch2")
	set(CATCH_URL "https://github.com/catchorg/Catch2/releases/download/v3.14.0")
	set(CATCH_MD5c "c7c5431b9bef8a27bfad0a8a45581ba7")
	set(CATCH_MD5h "a0e371008e8fd95bab1ff2fc8ff2058c")
	file(DOWNLOAD "${CATCH_URL}/catch_amalgamated.cpp" "${CATCH_ROOT}/catch_amalgamated.cpp" EXPECTED_MD5 ${CATCH_MD5c})
	file(DOWNLOAD "${CATCH_URL}/catch_amalgamated.hpp" "${CATCH_ROOT}/catch_amalgamated.hpp" EXPECTED_MD5 ${CATCH_MD5h})
	add_library(Catch2 STATIC "${CATCH_ROOT}/catch_amalgamated.cpp" "${CATCH_ROOT}/catch_amalgamated.hpp")
	target_include_directories(Catch2 SYSTEM PUBLIC "${CATCH_ROOT}")

	# Longdust is a C library (no CMakeLists.txt), so FetchContent_Populate
	# downloads the sources; tests/CMakeLists.txt compiles just the files it needs.
	lancet_dep_status("Configuring longdust v1.4 ...")
	FetchContent_Declare(longdust GIT_REPOSITORY https://github.com/lh3/longdust GIT_TAG v1.4 SYSTEM)
	FetchContent_GetProperties(longdust)
	if (NOT longdust_POPULATED)
		FetchContent_Populate(longdust)
	endif ()
endif ()

if (LANCET_BENCHMARKS)
	lancet_dep_status("Configuring Google Benchmark v1.9.5 ...")
	set(BENCHMARK_ENABLE_TESTING OFF)
	set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
	set(BENCHMARK_ENABLE_ASSEMBLY_TESTS OFF)
	set(BENCHMARK_ENABLE_INSTALL OFF)
	set(BENCHMARK_ENABLE_INSTALL OFF)
	set(BENCHMARK_INSTALL_DOCS OFF)
	set(BENCHMARK_ENABLE_LTO OFF)
	# Pre-set to skip benchmark's try_run regex probe. The probe fails in
	# -static builds because the tiny test binary can't resolve pthreads
	# symbols from libc++abi.a at runtime. libc++'s std::regex works fine.
	set(HAVE_STD_REGEX TRUE CACHE BOOL "" FORCE)
	FetchContent_Declare(benchmark GIT_REPOSITORY https://github.com/google/benchmark.git GIT_TAG v1.9.5 SYSTEM)
	FetchContent_MakeAvailable(benchmark)
endif ()

# Restore default log level after all dependencies are configured
set(CMAKE_MESSAGE_LOG_LEVEL "")
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS OFF)
