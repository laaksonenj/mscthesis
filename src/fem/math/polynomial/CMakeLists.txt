add_library(fem_math_polynomial_lib STATIC
    Monomial1D.cpp
    Monomial2D.cpp
    Polynomial1D.cpp
    Polynomial2D.cpp
    PolynomialCalculus.cpp
    PolynomialComposition.cpp
    PolynomialStringUtils.cpp
)

target_include_directories(fem_math_polynomial_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(fem_math_polynomial_lib
    PUBLIC
        fem_domain_lib
        fem_multiprecision_lib
    PRIVATE
        Boost::boost
)

add_subdirectory(ut)
