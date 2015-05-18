#
#  Copyright 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
cmake_minimum_required(VERSION 2.8.12)
project(ConstExprBug)

# Ensure that we can find the FindCXXFeatures.cmake module to determine the
# compiler flags needed to enable C++11
list(APPEND CMAKE_MODULE_PATH 
  "${CMAKE_CURRENT_SOURCE_DIR}/../../../external/cmake/")
include(FindCXXFeatures)
message(STATUS "CXX11_COMPILER_FLAGS=${CXX11_COMPILER_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11_COMPILER_FLAGS}")

option(TEST_CONSTEXPR "Test the MSVC bug via the use of a constexpr?" ON)
if(TEST_CONSTEXPR)
  add_definitions(-DUSE_CONSTEXPR)
endif()

option(BUILD_SHARED_LIBS "Build shared libraries?" ON)
file(GLOB CONSTEXPR_CPP RELATIVE ${PROJECT_SOURCE_DIR} "*.cpp")

add_library(ConstExpr ${CONSTEXPR_CPP})
