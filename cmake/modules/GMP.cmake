find_package(PkgConfig REQUIRED)

pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)
pkg_check_modules(gmpxx REQUIRED IMPORTED_TARGET gmpxx)

add_library(gmp INTERFACE IMPORTED GLOBAL)
target_link_libraries(gmp INTERFACE PkgConfig::gmp PkgConfig::gmpxx)
add_library(gmp::gmp ALIAS gmp)
