add_library(fem_math_quadrature_lib
    STATIC
        GaussLegendreTable1D.cpp
        GaussLegendreTableQuadrilateral.cpp
        GaussLegendreTableTriangle.cpp
        Quadrature1D.cpp
        Quadrature2D.cpp
)

target_include_directories(fem_math_quadrature_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(fem_math_quadrature_lib
    PUBLIC
        fem_math_common_lib
        fem_multiprecision_lib
        GSL::gsl
)

add_subdirectory(ut)
