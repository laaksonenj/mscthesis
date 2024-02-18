add_unit_test(fem_assembly_test
    SOURCES
        DiracLoadVectorTest.cpp
        StiffnessMatrixTest.cpp
    LIBRARIES
        fem_assembly_lib
)

target_compile_definitions(fem_assembly_test
    PRIVATE
        SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)
