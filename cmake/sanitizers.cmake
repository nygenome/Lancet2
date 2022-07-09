function(append value)
    foreach (variable ${ARGN})
        set(${variable} "${${variable}} ${value}" PARENT_SCOPE)
    endforeach (variable)
endfunction()

if (LANCET2_SANITIZER MATCHES "([Mm]emory)" AND APPLE)
    message(FATAL_ERROR "Memory sanitizer is not supported in ${CMAKE_HOST_SYSTEM}")
endif ()

if (LANCET2_SANITIZER)
    message(STATUS "Adding debug symbols for sanitizer builds")
    append("-g -fno-omit-frame-pointer -fsanitize-blacklist=${CMAKE_SOURCE_DIR}/cmake/sanitizer-blacklist.txt" CMAKE_C_FLAGS_DEBUG CMAKE_CXX_FLAGS_DEBUG)

    if (UNIX)
        if (LANCET2_SANITIZER MATCHES "([Aa]ddress)") # ASAN
            message(STATUS "Using Address sanitizer to detect memory corruption")
            append("-fsanitize=address -fsanitize-address-use-after-scope -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-common" CMAKE_C_FLAGS_DEBUG CMAKE_CXX_FLAGS_DEBUG)
            set(ENV{ASAN_OPTIONS} "detect_stack_use_after_return=true halt_on_error=true detect_leaks=1")
        elseif (LANCET2_SANITIZER MATCHES "([Uu]ndefined)") # UBSAN
            message(STATUS "Using Undefined sanitizer to perform runtime checks for undefined behavior")
            append("-fsanitize=undefined -fsanitize=implicit-integer-truncation -fsanitize=implicit-integer-arithmetic-value-change -fsanitize=implicit-conversion -fsanitize=integer" CMAKE_C_FLAGS_DEBUG CMAKE_CXX_FLAGS_DEBUG)
        elseif (LANCET2_SANITIZER MATCHES "([Tt]hread)") # TSAN
            message(STATUS "Using Thread sanitizer to detect data races")
            append("-fsanitize=thread" CMAKE_C_FLAGS_DEBUG CMAKE_CXX_FLAGS_DEBUG)
        elseif (LANCET2_SANITIZER MATCHES "([Mm]emory)") # MSAN
            message(STATUS "Using Memory sanitizer to detect uninitialized memory reads")
            append("-fsanitize=memory -fsanitize-memory-track-origins -fno-optimize-sibling-calls" CMAKE_C_FLAGS_DEBUG CMAKE_CXX_FLAGS_DEBUG)
        else ()
            message(FATAL_ERROR "Unsupported SANITIZER value: " LANCET2_SANITIZER)
        endif ()
    else ()
        message(FATAL_ERROR "LANCET2w_SANITIZER is not supported on ${CMAKE_HOST_SYSTEM}")
    endif ()
endif ()
