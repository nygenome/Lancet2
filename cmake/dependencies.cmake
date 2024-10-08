include(ExternalProject)
include(FetchContent)
include(ProcessorCount)
ProcessorCount(NumCores)
find_program(MAKE_EXE NAMES gmake nmake make REQUIRED)

set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

set(MI_SECURE OFF)
set(MI_PADDING OFF)
set(MI_OVERRIDE ON)
set(MI_BUILD_STATIC ON)
set(MI_BUILD_SHARED OFF)
set(MI_BUILD_OBJECT OFF)
set(MI_BUILD_TESTS OFF)
FetchContent_Declare(mimalloc GIT_REPOSITORY https://github.com/microsoft/mimalloc.git GIT_TAG v2.1.7 SYSTEM)
FetchContent_MakeAvailable(mimalloc)

FetchContent_Declare(abseil GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git GIT_TAG 0df5674 SYSTEM)
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

FetchContent_Declare(spdlog GIT_REPOSITORY https://github.com/gabime/spdlog.git GIT_TAG v1.14.1 SYSTEM)
FetchContent_MakeAvailable(spdlog)

FetchContent_Declare(cli11 GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git GIT_TAG v2.4.2 SYSTEM)
FetchContent_MakeAvailable(cli11)

set(CONCURRENTQUEUE_GIT_REPO "https://github.com/cameron314/concurrentqueue.git")
FetchContent_Declare(concurrentqueue GIT_REPOSITORY ${CONCURRENTQUEUE_GIT_REPO} GIT_TAG v1.0.4 SYSTEM)
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
FetchContent_Declare(libdeflate GIT_REPOSITORY https://github.com/ebiggers/libdeflate.git GIT_TAG v1.21 SYSTEM)
FetchContent_MakeAvailable(libdeflate)

set(ZLIB_COMPAT ON)
set(ZLIB_ENABLE_TESTS OFF)
set(ZLIBNG_ENABLE_TESTS OFF)
set(WITH_GTEST OFF)
set(WITH_BENCHMARKS OFF)
set(BUILD_SHARED_LIBS OFF)
set(WITH_NEW_STRATEGIES ON)
if (${LANCET_USER_REQUESTED_ARCH})
	set(WITH_OPTIM OFF)
	set(WITH_NATIVE_INSTRUCTIONS OFF)
else ()
	set(WITH_OPTIM ON)
	set(WITH_NATIVE_INSTRUCTIONS ON)
endif ()
FetchContent_Declare(zlib-ng GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git GIT_TAG 2.1.8 SYSTEM)
FetchContent_MakeAvailable(zlib-ng)

set(HTSLIB_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/htslib")
set(LIB_HTS "${HTSLIB_ROOT_DIR}/libhts.a")
set(HTSLIB_CONFIG_PARAMS ${HTSLIB_ROOT_DIR} ${CMAKE_C_COMPILER})
ExternalProject_Add(htslib
		URL https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
		URL_MD5 e10407ad47233ab762e2aa6eefec9921 PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
		SOURCE_DIR ${HTSLIB_ROOT_DIR} BUILD_IN_SOURCE 1 INSTALL_COMMAND ""
		BUILD_COMMAND ${MAKE_EXE} -j${NumCores} lib-static BUILD_BYPRODUCTS ${LIB_HTS}
		CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/configure_htslib.sh ${HTSLIB_CONFIG_PARAMS}
		LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON USES_TERMINAL_DOWNLOAD OFF
		USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)
add_dependencies(htslib zlibstatic libdeflate_static)

set(MM2_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/minimap2")
set(LIB_MM2 "${MM2_ROOT_DIR}/libminimap2.a")
set(MM2_BUILD_PARAMS ${MM2_ROOT_DIR} ${CMAKE_C_COMPILER})
ExternalProject_Add(minimap2
		URL https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2
		URL_MD5 600bd8d2d9f365d8c5dc6bc38b89ee8c PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
		SOURCE_DIR ${MM2_ROOT_DIR} BUILD_IN_SOURCE 1 CONFIGURE_COMMAND "" INSTALL_COMMAND ""
		BUILD_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/build_minimap2.sh ${MM2_BUILD_PARAMS}
		BUILD_BYPRODUCTS ${LIB_MM2} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
		USES_TERMINAL_DOWNLOAD OFF USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)
add_dependencies(minimap2 zlibstatic)

set(GPERFTOOLS_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/gperftools")
set(GPERFTOOLS_INC_DIR "${GPERFTOOLS_ROOT_DIR}/include")
set(LIB_PROFILER "${GPERFTOOLS_ROOT_DIR}/lib/libprofiler.a")
set(GPERFTOOLS_CONFIG_PARAMS ${GPERFTOOLS_ROOT_DIR} ${CMAKE_C_COMPILER} ${CMAKE_CXX_COMPILER})
ExternalProject_Add(gperftools
		URL https://github.com/gperftools/gperftools/releases/download/gperftools-2.15.90/gperftools-2.15.90.tar.gz
		URL_MD5 fbab1d58619c76857a2a46fc02b9dffa
		PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${GPERFTOOLS_ROOT_DIR} BUILD_IN_SOURCE 1
		INSTALL_COMMAND ${MAKE_EXE} install BUILD_COMMAND ${MAKE_EXE} -j${NumCores}
		CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/configure_gperftools.sh ${GPERFTOOLS_CONFIG_PARAMS}
		BUILD_BYPRODUCTS ${LIB_PROFILER} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
		USES_TERMINAL_DOWNLOAD OFF USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)

set(spoa_optimize_for_native OFF)
FetchContent_Declare(spoa GIT_REPOSITORY https://github.com/rvaser/spoa GIT_TAG 73f7a25 SYSTEM)
FetchContent_MakeAvailable(spoa)

FetchContent_Declare(boost_math GIT_REPOSITORY https://github.com/boostorg/math.git GIT_TAG boost-1.86.0 SYSTEM)
FetchContent_MakeAvailable(boost_math)

if (LANCET_TESTS)
	file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_deps/Catch2")
	set(CATCH_ROOT "${CMAKE_CURRENT_BINARY_DIR}/_deps/Catch2")
	set(CATCH_URL "https://github.com/catchorg/Catch2/releases/download/v3.7.1")
	set(CATCH_MD5c "dc70126409d62593ce53f45830f98065")
	set(CATCH_MD5h "9466b920b07009697e2a152284cfa221")
	file(DOWNLOAD "${CATCH_URL}/catch_amalgamated.cpp" "${CATCH_ROOT}/catch_amalgamated.cpp" EXPECTED_MD5 ${CATCH_MD5c})
	file(DOWNLOAD "${CATCH_URL}/catch_amalgamated.hpp" "${CATCH_ROOT}/catch_amalgamated.hpp" EXPECTED_MD5 ${CATCH_MD5h})
	add_library(Catch2 STATIC "${CATCH_ROOT}/catch_amalgamated.cpp" "${CATCH_ROOT}/catch_amalgamated.hpp")
	target_include_directories(Catch2 SYSTEM PUBLIC "${CATCH_ROOT}")
endif ()

if (LANCET_BENCHMARKS)
	set(BENCHMARK_ENABLE_TESTING OFF)
	set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
	set(BENCHMARK_ENABLE_ASSEMBLY_TESTS OFF)
	set(BENCHMARK_ENABLE_INSTALL OFF)
	set(BENCHMARK_ENABLE_INSTALL OFF)
	set(BENCHMARK_INSTALL_DOCS OFF)
	set(BENCHMARK_ENABLE_LTO OFF)
	FetchContent_Declare(benchmark GIT_REPOSITORY https://github.com/google/benchmark.git GIT_TAG v1.9.0 SYSTEM)
	FetchContent_MakeAvailable(benchmark)
endif ()
