add_library(apps_refdata_utils_lib STATIC
    Utils.cpp
)

target_include_directories(apps_refdata_utils_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(apps_refdata_utils_lib
    PUBLIC
        fem_multiprecision_lib
)

add_executable(generate_stiffness_matrix_header
    GenerateStiffnessMatrixHeaderMain.cpp
)

target_link_libraries(generate_stiffness_matrix_header
    PRIVATE
        apps_common_lib
        apps_refdata_utils_lib
        fem_assembly_lib
)

add_executable(generate_dirac_load_header
    GenerateDiracLoadHeaderMain.cpp
)

target_link_libraries(generate_dirac_load_header
    PRIVATE
        apps_common_lib
        apps_refdata_utils_lib
        fem_assembly_lib
)

add_executable(generate_neumann_load_header
    GenerateNeumannLoadHeaderMain.cpp
)

target_link_libraries(generate_neumann_load_header
    PRIVATE
        apps_common_lib
        apps_refdata_utils_lib
        fem_assembly_lib
        fem_assembly_ut_refdata_lib
)
