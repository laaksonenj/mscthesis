function(add_unit_test testname)
    cmake_parse_arguments(
        PARSE_ARGV 1
        ARG
        "" "" "SOURCES;LIBRARIES"
    )

    set(TEST_EXE_NAME "${testname}_app")
    set(VCPKG_APPLOCAL_DEPS OFF)
    add_executable(${TEST_EXE_NAME} ${ARG_SOURCES})
    target_link_libraries(${TEST_EXE_NAME}
        PRIVATE
            GTest::gtest
            GTest::gtest_main
            ${ARG_LIBRARIES}
    )

    add_test(NAME ${testname} COMMAND ${TEST_EXE_NAME})
    set_tests_properties(${testname} PROPERTIES
        ENVIRONMENT_MODIFICATION PATH=path_list_prepend:$<TARGET_FILE_DIR:GTest::gtest_main>
    )
endfunction()
