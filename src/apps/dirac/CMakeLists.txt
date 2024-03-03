add_subdirectory(utils)

add_executable(dirac_pfem
    DiracPfemMain.cpp
)

target_link_libraries(dirac_pfem
    PRIVATE
        apps_common_dirac_utils_lib
        apps_common_lib
        fem_assembly_lib
        fem_basis_lib
        fem_math_lib
        fem_multiprecision_lib
        OpenMP::OpenMP_CXX
)
