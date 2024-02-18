add_library(fem_assembly_lib STATIC
    DiracLoadVector.cpp
    StiffnessMatrix.cpp
)

target_include_directories(fem_assembly_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(fem_assembly_lib
    PUBLIC
        fem_basis_lib
        fem_domain_lib
        fem_math_lib
        fem_multiprecision_lib
)

add_subdirectory(ut)
