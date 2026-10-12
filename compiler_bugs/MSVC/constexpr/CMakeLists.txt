#
#  Copyright 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
cmake_minimum_required(VERSION 3.5)
project(ConstExprBug CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(TEST_CONSTEXPR "Test the MSVC bug via the use of a constexpr?" ON)
if(TEST_CONSTEXPR)
  add_definitions(-DUSE_CONSTEXPR)
endif()

option(BUILD_SHARED_LIBS "Build shared libraries?" ON)
file(GLOB CONSTEXPR_CPP RELATIVE ${PROJECT_SOURCE_DIR} "*.cpp")

add_library(ConstExpr ${CONSTEXPR_CPP})
