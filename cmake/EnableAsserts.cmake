foreach(compile_flags
        CMAKE_C_FLAGS_DEBUG
        CMAKE_C_FLAGS_RELWITHDEBINFO
        CMAKE_C_FLAGS_MINSIZEREL
        CMAKE_C_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_DEBUG
        CMAKE_CXX_FLAGS_RELWITHDEBINFO
        CMAKE_CXX_FLAGS_MINSIZEREL
        CMAKE_CXX_FLAGS_RELEASE)
    string(REGEX REPLACE "[/-]D *NDEBUG($| )" ""
        ${compile_flags} ${${compile_flags}})
endforeach()
