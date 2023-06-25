include(ExternalProject)
include(FetchContent)
include(ProcessorCount)
ProcessorCount(NumCores)
find_program(MAKE_EXE NAMES gmake nmake make REQUIRED)
find_program(AUTOCONF_EXE NAMES autoconf REQUIRED)
find_program(AUTOHEADER_EXE NAMES autoheader REQUIRED)
find_program(AUTORECONF_EXE NAMES autoreconf REQUIRED)

message(STATUS "Setting up library dependencies required for the build")
function(message)
	if (NOT MESSAGE_QUIET)
		_message(${ARGN})
	endif ()
endfunction()

set(MESSAGE_QUIET ON)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

set(MI_SECURE ON)
set(MI_PADDING OFF)
set(MI_OVERRIDE ON)
set(MI_BUILD_STATIC ON)
set(MI_BUILD_SHARED OFF)
set(MI_BUILD_OBJECT OFF)
set(MI_BUILD_TESTS OFF)
if (${CMAKE_BUILD_TYPE} MATCHES Debug)
	set(MI_DEBUG_FULL ON)
endif ()
FetchContent_Declare(mimalloc GIT_REPOSITORY https://github.com/microsoft/mimalloc.git GIT_TAG v2.1.2 SYSTEM)
FetchContent_MakeAvailable(mimalloc)
target_compile_options(mimalloc-static PRIVATE "-Wno-stringop-overflow")

FetchContent_Declare(abseil GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git GIT_TAG 94d77fe SYSTEM)
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

FetchContent_Declare(spdlog GIT_REPOSITORY https://github.com/gabime/spdlog.git GIT_TAG v1.11.0 SYSTEM)
FetchContent_MakeAvailable(spdlog)

set(LIBDEFLATE_BUILD_STATIC_LIB ON)
set(LIBDEFLATE_BUILD_SHARED_LIB OFF)
set(LIBDEFLATE_BUILD_GZIP OFF)
set(LIBDEFLATE_BUILD_TESTS OFF)
set(LIBDEFLATE_USE_SHARED_LIB OFF)
FetchContent_Declare(libdeflate GIT_REPOSITORY https://github.com/ebiggers/libdeflate.git GIT_TAG v1.18 SYSTEM)
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
FetchContent_Declare(zlib-ng GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git GIT_TAG 2.1.2 SYSTEM)
FetchContent_MakeAvailable(zlib-ng)

set(BZIP2_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/bzip2")
set(LIB_BZ2 "${BZIP2_ROOT_DIR}/libbz2.a")
ExternalProject_Add(bzip2
		GIT_REPOSITORY git://sourceware.org/git/bzip2.git GIT_TAG bzip2-1.0.8
		PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${BZIP2_ROOT_DIR}
		BUILD_IN_SOURCE 1 CONFIGURE_COMMAND "" INSTALL_COMMAND ""
		BUILD_COMMAND ${MAKE_EXE} -j${NumCores} libbz2.a "CFLAGS=${CMAKE_C_FLAGS_RELEASE}"
		BUILD_BYPRODUCTS ${LIB_BZ2} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
		USES_TERMINAL_DOWNLOAD OFF USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)

set(LZMA_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/xz")
set(LIB_LZMA "${LZMA_ROOT_DIR}/lib/liblzma.a")
set(LZMA_CONFIG_PARAMS ${LZMA_ROOT_DIR} ${CMAKE_C_COMPILER})
ExternalProject_Add(lzma
		URL https://tukaani.org/xz/xz-5.4.2.tar.gz URL_MD5 4ac4e5da95aa8604a81e32079cb00d42
		PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${LZMA_ROOT_DIR}
		BUILD_IN_SOURCE 1 INSTALL_COMMAND ${MAKE_EXE} install BUILD_COMMAND ${MAKE_EXE}
		CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/configure_lzma.sh ${LZMA_CONFIG_PARAMS}
		BUILD_BYPRODUCTS ${LIB_LZMA} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
		USES_TERMINAL_DOWNLOAD OFF USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)

set(HTSLIB_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/htslib")
set(LIB_HTS "${HTSLIB_ROOT_DIR}/libhts.a")
set(HTSLIB_CONFIG_PARAMS ${HTSLIB_ROOT_DIR} ${CMAKE_C_COMPILER})
ExternalProject_Add(htslib
		URL https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
		URL_MD5 ff47f4b09e5202cebc3f80e7e02c7728 PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
		SOURCE_DIR ${HTSLIB_ROOT_DIR} BUILD_IN_SOURCE 1 INSTALL_COMMAND ""
		BUILD_COMMAND ${MAKE_EXE} -j${NumCores} lib-static BUILD_BYPRODUCTS ${LIB_HTS}
		CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/configure_htslib.sh ${HTSLIB_CONFIG_PARAMS}
		LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON USES_TERMINAL_DOWNLOAD OFF
		USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)
add_dependencies(htslib zlibstatic bzip2 lzma libdeflate_static)

set(MM2_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/minimap2")
set(LIB_MM2 "${MM2_ROOT_DIR}/libminimap2.a")
set(MM2_BUILD_PARAMS ${MM2_ROOT_DIR} ${CMAKE_C_COMPILER})
ExternalProject_Add(minimap2
		URL https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26.tar.bz2
		URL_MD5 b55b69773f07b15ceae9a6b86d907c78 PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
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
		URL https://github.com/gperftools/gperftools/releases/download/gperftools-2.10/gperftools-2.10.tar.gz
		URL_MD5 62bf6c76ba855ed580de5e139bd2a483
		PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${GPERFTOOLS_ROOT_DIR}
		BUILD_IN_SOURCE 1 INSTALL_COMMAND ${MAKE_EXE} install BUILD_COMMAND ${MAKE_EXE} -j${NumCores}
		CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/configure_gperftools.sh ${GPERFTOOLS_CONFIG_PARAMS}
		BUILD_BYPRODUCTS ${LIB_PROFILER} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
		USES_TERMINAL_DOWNLOAD OFF USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)

FetchContent_Declare(wyhash GIT_REPOSITORY https://github.com/wangyi-fudan/wyhash.git GIT_TAG 77e50f2 SYSTEM)
FetchContent_GetProperties(wyhash)
if (NOT wyhash_POPULATED)
	FetchContent_Populate(wyhash)
	add_library(wyhash INTERFACE)
	target_include_directories(wyhash SYSTEM INTERFACE "${wyhash_SOURCE_DIR}")
endif ()

set(spoa_optimize_for_native OFF)
FetchContent_Declare(spoa GIT_REPOSITORY https://github.com/rvaser/spoa GIT_TAG 08957f6 SYSTEM)
FetchContent_MakeAvailable(spoa)

FetchContent_Declare(cli11 GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git GIT_TAG v2.3.2 SYSTEM)
FetchContent_MakeAvailable(cli11)

set(CONCURRENTQUEUE_GIT_REPO "https://github.com/cameron314/concurrentqueue.git")
FetchContent_Declare(concurrentqueue GIT_REPOSITORY ${CONCURRENTQUEUE_GIT_REPO} GIT_TAG v1.0.4 SYSTEM)
FetchContent_GetProperties(concurrentqueue)
if (NOT concurrentqueue_POPULATED)
	FetchContent_Populate(concurrentqueue)
	add_library(concurrentqueue INTERFACE)
	target_include_directories(concurrentqueue SYSTEM INTERFACE "${concurrentqueue_SOURCE_DIR}")
endif ()

set(ROARING_ROOT "${CMAKE_CURRENT_BINARY_DIR}/_deps/roaring")
set(ROARING_URL "https://github.com/RoaringBitmap/CRoaring/releases/download/v1.3.0")
file(MAKE_DIRECTORY "${ROARING_ROOT}")
file(DOWNLOAD "${ROARING_URL}/roaring.c" "${ROARING_ROOT}/roaring.c" EXPECTED_MD5 "a9440ccdb38715aac83fca86dc4809d0")
file(DOWNLOAD "${ROARING_URL}/roaring.h" "${ROARING_ROOT}/roaring.h" EXPECTED_MD5 "fe2b2d24d80c4df5f8dfdf6138564426")
file(DOWNLOAD "${ROARING_URL}/roaring.hh" "${ROARING_ROOT}/roaring.hh" EXPECTED_MD5 "b141f2f268d067e522eed161f480c714")
add_library(RoaringBitmap STATIC "${ROARING_ROOT}/roaring.c" "${ROARING_ROOT}/roaring.h" "${ROARING_ROOT}/roaring.hh")
target_include_directories(RoaringBitmap SYSTEM PUBLIC "${ROARING_ROOT}")

if (LANCET_TESTS)
	file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_deps/Catch2")
	set(CATCH_ROOT "${CMAKE_CURRENT_BINARY_DIR}/_deps/Catch2")
	set(CATCH_URL "https://github.com/catchorg/Catch2/releases/download/v3.3.2")
	set(CATCH_MD5c "39762707976cb30c301b3762c55a58e6")
	set(CATCH_MD5h "e68fa70a3fece7c6de62f224900a9d23")
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
	FetchContent_Declare(benchmark GIT_REPOSITORY https://github.com/google/benchmark.git GIT_TAG v1.8.0 SYSTEM)
	FetchContent_MakeAvailable(benchmark)
endif ()

unset(MESSAGE_QUIET)
