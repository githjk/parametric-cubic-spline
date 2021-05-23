cmake_minimum_required(VERSION 3.5)
project(parametric_cubic_spline)

set(${PROJECT_NAME}_MAJOR_VERSION 0)
set(${PROJECT_NAME}_MINOR_VERSION 1)
set(${PROJECT_NAME}_PATCH_VERSION 0)
set(${PROJECT_NAME}_VERSION
  ${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${${PROJECT_NAME}_PATCH_VERSION})

# ------------------------------------------------------------------------------
# Compiler Options
# ------------------------------------------------------------------------------
# Default to C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Build
# ------------------------------------------------------------------------------
include_directories(include)

add_library(${PROJECT_NAME} INTERFACE)

# ------------------------------------------------------------------------------
# Test
# ------------------------------------------------------------------------------
enable_testing()
add_subdirectory(test)

# ------------------------------------------------------------------------------
# Install
# ------------------------------------------------------------------------------
