add_unit_test(fem_math_test
    SOURCES
        AffineMapTest.cpp
        QuadratureTest.cpp
    LIBRARIES
        fem_math_lib
        GSL::gsl
)
