add_library(apps_common_lib STATIC
    Arguments.cpp
    LinearSolver.cpp
    Timer.cpp
)

target_include_directories(apps_common_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(apps_common_lib
    PUBLIC
        Boost::program_options
    PRIVATE
        fem_basis_lib
        fem_multiprecision_lib
)
