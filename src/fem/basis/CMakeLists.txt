add_library(fem_basis_lib STATIC
    BasisFunctionFactory.cpp
    BasisFunctionIndexer.cpp
    ShapeFunctionFactory.cpp
    ShapeFunctionIndexer.cpp
    #TrialFunction.cpp
)

target_include_directories(fem_basis_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(fem_basis_lib
    PUBLIC
        Boost::boost
        fem_domain_lib
        fem_math_lib
    PRIVATE
        OpenMP::OpenMP_CXX
)

add_subdirectory(ut)
