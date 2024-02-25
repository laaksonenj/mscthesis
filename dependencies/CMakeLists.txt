set(Boost_NO_WARN_NEW_VERSIONS 1)
find_package(Boost REQUIRED GLOBAL
    COMPONENTS
        program_options
)
find_package(Eigen3 CONFIG REQUIRED GLOBAL)
find_package(GMP REQUIRED GLOBAL)
find_package(GSL REQUIRED GLOBAL)
find_package(GTest CONFIG REQUIRED GLOBAL)
find_package(OpenMP REQUIRED GLOBAL
    COMPONENTS
        CXX
)
