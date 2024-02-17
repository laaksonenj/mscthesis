add_subdirectory(polynomial)

add_library(fem_math_lib STATIC
    AffineMap.cpp
    Quadrature.cpp
)

target_include_directories(fem_math_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(fem_math_lib
    PUBLIC
        fem_domain_lib
        fem_multiprecision_lib
        fem_math_polynomial_lib
    PRIVATE
        GSL::gsl
)

add_subdirectory(ut)
