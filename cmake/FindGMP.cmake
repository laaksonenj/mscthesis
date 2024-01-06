find_package(PkgConfig REQUIRED)

pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)
pkg_check_modules(gmpxx REQUIRED IMPORTED_TARGET gmpxx)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
    DEFAULT_MSG
    gmp_FOUND
    gmpxx_FOUND
)

if(GMP_FOUND)
    add_library(GMP INTERFACE IMPORTED GLOBAL)
    target_link_libraries(GMP INTERFACE PkgConfig::gmp PkgConfig::gmpxx)
    add_library(GMP::GMP ALIAS GMP)
endif()
