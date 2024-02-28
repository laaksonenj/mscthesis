add_subdirectory(polynomial)
add_subdirectory(quadrature)

add_library(fem_math_common_lib STATIC
    AffineMap.cpp
)

target_include_directories(fem_math_common_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(fem_math_common_lib
    PUBLIC
        fem_multiprecision_lib
)

add_library(fem_math_lib INTERFACE)

target_link_libraries(fem_math_lib
    INTERFACE
        fem_math_common_lib
        fem_math_polynomial_lib
        fem_math_quadrature_lib
)

add_subdirectory(ut)
