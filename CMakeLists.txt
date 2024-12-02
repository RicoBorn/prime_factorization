cmake_minimum_required(VERSION 3.15)

# Project name and version
project(prime_factorization VERSION 1.0 LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find GMP library
find_package(PkgConfig REQUIRED)
pkg_check_modules(GMPXX REQUIRED IMPORTED_TARGET gmpxx)


# Gather all source and header files
set(SOURCES
        main.cpp
        elliptic_curve.cpp
        helper.cpp
        trial_division.cpp
        test_numbers.cpp
)

set(HEADERS
        elliptic_curve.h
        helper.h
        trial_division.h
        test_numbers.h
)

# Add the executable
add_executable(prime_factorization ${SOURCES} ${HEADERS})

# Link the GMP and GMPXX libraries
target_link_libraries(prime_factorization PkgConfig::GMPXX)

# Additional compiler flags for warnings (optional)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(prime_factorization PRIVATE -Wall -Wextra -pedantic)
endif()

# Enable IDE project files generation (useful for CLion)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT prime_factorization)