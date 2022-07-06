# Set default policies to ensure same behaviour
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0069 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0077 NEW)

set(CMAKE_POLICY_DEFAULT_CMP0012 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0074 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

# Use ccache if found to cache previously built object files
find_program(CCACHE_EXE ccache)
if(CCACHE_EXE)
    message(STATUS "Found ccache in PATH. Using ccache to speed up recompilation")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ${CCACHE_EXE})
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ${CCACHE_EXE})
endif()

# Add Debug flags
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")
endif()
