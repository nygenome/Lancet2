# ═══════════════════════════════════════════════════════════════════════════════
# Version Detection
#
# Reads the base version from the VERSION file at the project root (single
# source of truth for CMake, Conda, and Docker). Optionally augments with
# Git branch and short SHA when available. Generates lancet_version.h for
# compile-time version embedding.
# ═══════════════════════════════════════════════════════════════════════════════
set(LANCET2_VERSION_TAG ${PROJECT_VERSION})

find_package(Git QUIET)
if (GIT_FOUND)
	execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			OUTPUT_VARIABLE LANCET2_GIT_BRANCH ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short=10 --verify HEAD WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			OUTPUT_VARIABLE LANCET2_GIT_REVISION ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
endif ()

if (NOT LANCET2_GIT_BRANCH)
	set(LANCET2_GIT_BRANCH "")
endif ()

if (NOT LANCET2_GIT_REVISION)
	set(LANCET2_GIT_REVISION "")
endif ()

set(LANCET2_VERSION_HEADER "${CMAKE_BINARY_DIR}/generated/lancet_version.h")
configure_file("${CMAKE_SOURCE_DIR}/src/lancet/version.h.inc" ${LANCET2_VERSION_HEADER} @ONLY)
