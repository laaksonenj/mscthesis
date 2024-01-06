function(add_unit_test testname)
    cmake_parse_arguments(
        PARSE_ARGV 1
        ARG
        "" "" "SOURCES;LIBRARIES"
    )

    set(VCPKG_APPLOCAL_DEPS OFF)
    add_executable(${testname} ${ARG_SOURCES})
    target_link_libraries(${testname}
        PRIVATE
            GTest::gtest
            GTest::gtest_main
            ${ARG_LIBRARIES}
    )

    add_test(NAME ${testname} COMMAND ${testname})
    set_tests_properties(${testname} PROPERTIES
        ENVIRONMENT_MODIFICATION PATH=path_list_prepend:${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}$<$<CONFIG:Debug>:/debug>/bin
    )
endfunction()
