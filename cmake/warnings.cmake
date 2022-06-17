function(target_set_warnings)
    if (NOT LANCET_WARNINGS)
        return()
    endif ()

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(WGCC TRUE)
    elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(WCLANG TRUE)
    endif ()

    set(multiValueArgs ENABLE DISABLE AS_ERROR)
    cmake_parse_arguments(this "" "" "${multiValueArgs}" ${ARGN})
    list(FIND this_ENABLE "ALL" enable_all)
    list(FIND this_DISABLE "ALL" disable_all)
    list(FIND this_AS_ERROR "ALL" as_error_all)

    if (NOT ${enable_all} EQUAL -1)
        if (WGCC)
            list(APPEND WarningFlags "-Wall" "-Wextra")
        elseif (WCLANG)
            list(APPEND WarningFlags "-Wall" "-Weverything" "-Wpedantic")
        endif ()
    elseif (NOT ${disable_all} EQUAL -1)
        # Assume all includes coming from system
        set(SystemIncludes TRUE)
        list(APPEND WarningFlags "-w")
    endif ()

    list(FIND this_DISABLE "Annoying" disable_annoying)
    if (NOT ${disable_annoying} EQUAL -1)
        if (WCLANG)
            list(APPEND WarningFlags ${ABSL_DEFAULT_COPTS} "-Wno-newline-eof"
                    "-Wno-documentation-unknown-command" "-Wno-documentation"
                    "-Wno-unused-variable" "-Wno-missing-prototypes" "-Wno-c++2a-compat"
                    "-Wno-signed-enum-bitfield" "-Wno-undefined-func-template" "-Wno-shadow"
                    "-Wno-c++98-compat" "-Wno-c++98-compat-pedantic" "-Wno-padded"
                    "-Wno-old-style-cast" "-Wno-reserved-id-macro" "-Wno-gcc-compat"
                    "-Wno-zero-as-null-pointer-constant" "-Wno-thread-safety-negative"
                    "-Wno-double-promotion" "-Wno-weak-vtables" "-Wno-covered-switch-default"
                    "-Wno-extra-semi-stmt" "-Wno-switch-enum" "-Wno-global-constructors"
                    "-Wno-float-equal" "-Wno-exit-time-destructors" "-Wno-reserved-identifier"
                    "-Wno-disabled-macro-expansion" "-Wno-unused-parameter")
        endif ()

        if (WGCC)
            list(APPEND WarningFlags ${ABSL_DEFAULT_COPTS} "-Wno-redundant-move"
                    "-Wno-unused-parameter")
        endif ()
    endif ()

    if (NOT ${as_error_all} EQUAL -1)
        list(APPEND WarningFlags "-Werror")
    endif ()

    foreach (target IN LISTS this_UNPARSED_ARGUMENTS)
        if (WarningFlags)
            target_compile_options(${target} PRIVATE ${WarningFlags})
        endif ()

        if (WarningDefinitions)
            target_compile_definitions(${target} PRIVATE ${WarningDefinitions})
        endif ()

        if (SystemIncludes)
            set_target_properties(${target} PROPERTIES
                    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                    $<TARGET_PROPERTY:${target},INTERFACE_INCLUDE_DIRECTORIES>)
        endif ()
    endforeach ()
endfunction(target_set_warnings)
