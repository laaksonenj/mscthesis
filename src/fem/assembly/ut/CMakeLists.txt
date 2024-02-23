add_library(fem_assembly_ut_refdata_lib INTERFACE)

target_include_directories(fem_assembly_ut_refdata_lib
    INTERFACE
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(fem_assembly_ut_refdata_lib
    INTERFACE
        fem_math_lib
        fem_multiprecision_lib
        GSL::gsl
)

add_unit_test(fem_assembly_test
    SOURCES
        DiracLoadVectorTest.cpp
        LoadVectorTest.cpp
        NeumannLoadVectorTest.cpp
        StiffnessMatrixTest.cpp
    LIBRARIES
        fem_assembly_lib
        fem_assembly_ut_refdata_lib
)

target_compile_definitions(fem_assembly_test
    PRIVATE
        SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)
