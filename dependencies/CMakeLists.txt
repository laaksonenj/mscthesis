find_package(GMP REQUIRED GLOBAL)
find_package(Eigen3 CONFIG REQUIRED GLOBAL)

if (BUILD_TESTING)
    find_package(GTest CONFIG REQUIRED GLOBAL)
endif()
