add_library(fem_multiprecision_lib STATIC
    Arithmetic.cpp
)

target_include_directories(fem_multiprecision_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(fem_multiprecision_lib
    PUBLIC
        Eigen3::Eigen
        GMP::GMP
)

add_subdirectory(ut)
