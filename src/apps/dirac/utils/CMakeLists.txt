add_library(apps_common_dirac_utils_lib
    GreensFunction.cpp
    L2Error.cpp
)

target_include_directories(apps_common_dirac_utils_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(apps_common_dirac_utils_lib
    PUBLIC
        fem_basis_lib
        fem_math_lib
        fem_multiprecision_lib
        fem_domain_lib
)
