cmake_minimum_required(VERSION 3.24)

project(mscthesis)

option(ENABLE_RUNTIME_ASSERTS "" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

if(ENABLE_RUNTIME_ASSERTS)
    include(EnableAsserts)
endif()
include(CTest)
include(UtUtils)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/MP /bigobj)
endif()

add_subdirectory(dependencies)

add_subdirectory(src)
