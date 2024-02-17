add_unit_test(fem_basis_test
    SOURCES
        BasisFunctionIndexerTest.cpp
        ShapeFunctionFactoryTest.cpp
        ShapeFunctionIndexerTest.cpp
    LIBRARIES
        fem_basis_lib
)