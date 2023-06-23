# Set default policies to ensure same behaviour
if (POLICY CMP0012)
	cmake_policy(SET CMP0012 NEW)
endif ()

if (POLICY CMP0048)
	cmake_policy(SET CMP0048 NEW)
endif ()

if (POLICY CMP0063)
	cmake_policy(SET CMP0063 NEW)
endif ()

if (POLICY CMP0069)
	cmake_policy(SET CMP0069 NEW)
endif ()

if (POLICY CMP0074)
	cmake_policy(SET CMP0074 NEW)
endif ()

if (POLICY CMP0077)
	cmake_policy(SET CMP0077 NEW)
endif ()

if (POLICY CMP0135)
	cmake_policy(SET CMP0135 NEW)
endif ()

set(CMAKE_POLICY_DEFAULT_CMP0012 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0063 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0074 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0135 NEW)

# Use ccache if found to cache previously built object files
find_program(CCACHE_EXE ccache)
if (CCACHE_EXE)
	message(STATUS "Found ccache in PATH. Using ccache to speed up recompilation")
	set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ${CCACHE_EXE})
	set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ${CCACHE_EXE})
endif ()

# Optional IPO. Do not use IPO if it's not supported by compiler.
include(CheckIPOSupported)
check_ipo_supported(RESULT COMPILER_SUPPORTS_IPO LANGUAGES C CXX)
if (COMPILER_SUPPORTS_IPO AND ${CMAKE_BUILD_TYPE} MATCHES Release)
	message(STATUS "Enabling Link Time Optimization because compiler supports it")
	set(ENABLE_LTO ON)
endif ()

if (${LANCET_BUILD_ARCH} STREQUAL "OFF")
	DetectHostCpuArch()
	set(LANCET_USER_REQUESTED_ARCH OFF CACHE INTERNAL "" FORCE)
else ()
	set(LANCET_USER_REQUESTED_ARCH ON CACHE INTERNAL "" FORCE)
endif ()

if (${LANCET_BUILD_STATIC})
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
endif ()

message(STATUS "LANCET_BUILD_ARCH: ${LANCET_BUILD_ARCH}")
set(LANCET_OPT_FLAGS "-Ofast -DNDEBUG -m64 -march=${LANCET_BUILD_ARCH} -mtune=${LANCET_BUILD_ARCH}" CACHE STRING "")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} ${LANCET_OPT_FLAGS} -fno-omit-frame-pointer" CACHE INTERNAL "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${LANCET_OPT_FLAGS} -fno-omit-frame-pointer" CACHE INTERNAL "" FORCE)
