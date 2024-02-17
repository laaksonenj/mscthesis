add_library(fem_domain_lib STATIC
    Side.cpp
    Element.cpp
    Mesh.cpp
    Parallelogram.cpp
    Triangle.cpp
)

target_include_directories(fem_domain_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(fem_domain_lib
    PUBLIC
        fem_math_lib
        fem_multiprecision_lib
)

add_subdirectory(ut)
