#
#  Copyright (c) 2009-2010, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#   - Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#
#   - Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#   - Neither the name of the owner nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.
#
cmake_minimum_required(VERSION 2.6)
project(Elemental)
set(Elemental_VERSION_MAJOR 0)
set(Elemental_VERSION_MINOR 5)

option(BUILD_HYBRID "Also build a hybrid OpenMP/MPI version of Elemental" OFF)
option(BUILD_DEBUG "Also build a version with many consistency checks" OFF)
option(BUILD_TESTS "Also build a collection of test executables" ON)
option(BUILD_PMRRR "Build and link against the eigensolver PMRRR" ON)

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER ${MPI_COMPILER})

if(BUILD_HYBRID)
  find_package(OpenMP)
  if(NOT OpenMP_FOUND)
    message("Hybrid build disabled because OpenMP support was not found")
  endif(NOT OpenMP_FOUND)
endif(BUILD_HYBRID)

# A simple extension of find_library from a single lib to a list of libs
#
# Due to the find_library routine caching the return status for fixed variable 
# names within a loop, we must change it for each library.
macro(find_library_list FOUND_LIST SEARCH_LIST PATH PATH_SUFFIX)
  set(${FOUND_LIST})
  foreach(LIB ${${SEARCH_LIST}})
    find_library(FOUND_${LIB} ${LIB} PATHS ${PATH} PATH_SUFFIXES ${PATH_SUFFIX})
    if(FOUND_${LIB})
      list(APPEND ${FOUND_LIST} ${FOUND_${LIB}})
    endif(FOUND_${LIB})
  endforeach(LIB)
  # If FOUND_LIST is not the same length as SEARCH_LIST, then empty it so that
  # it is easy to check whether or not all of the libraries were found.
  list(LENGTH ${SEARCH_LIST} NUM_SEARCHED)
  list(LENGTH ${FOUND_LIST} NUM_FOUND)
  if(NUM_FOUND LESS NUM_SEARCHED)
    set(${FOUND_LIST})
  endif(NUM_FOUND LESS NUM_SEARCHED)
endmacro(find_library_list)

# Query the size of a void pointer in order to determine whether or not this is
# a 32-bit system
if(${CMAKE_SIZEOF_VOID_P} MATCHES 4)
  set(32_BIT_SYSTEM TRUE)
else(${CMAKE_SIZEOF_VOID_P} MATCHES 4)
  set(32_BIT_SYSTEM FALSE)    
endif(${CMAKE_SIZEOF_VOID_P} MATCHES 4)

# Since the CMake find_blas for some reason requires a Fortran compiler, we
# must build our own.
if(APPLE)
  set(BLAS_LAPACK_LIBS "-framework vecLib")
  set(FOUND_BLAS "vecLib")
  set(FOUND_LAPACK "vecLib")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Look for MKL first
  if(32_BIT_SYSTEM)
    set(MKL_SUBDIR lib/32)
  else(32_BIT_SYSTEM)
    set(MKL_SUBDIR lib/em64t)
  endif(32_BIT_SYSTEM)
  set(MKL_LIBS mkl_intel_lp64 mkl_intel_thread mkl_core guide iomp5)
  find_library_list(BLAS_LAPACK_LIBS MKL_LIBS "" ${MKL_SUBDIR})
  if(BLAS_LAPACK_LIBS)
    find_library(FOUND_PTHREAD pthread)
    if(FOUND_PTHREAD)
      set(FOUND_BLAS "MKL")    
      set(FOUND_LAPACK "MKL")
    else(FOUND_PTHREAD)
      message("Found MKL but could not find a libpthread")  
    endif(FOUND_PTHREAD)
  else(BLAS_LAPACK_LIBS)
    # Look for default BLAS and LAPACK
    set(REFERENCE_LIBS lapack blas)
    find_library_list(BLAS_LAPACK_LIBS REFERENCE_LIBS "" "")
    if(BLAS_LAPACK_LIBS)
      set(FOUND_BLAS "Default")
      set(FOUND_LAPACK "Default")
    else(BLAS_LAPACK_LIBS)
      message(FATAL_ERROR "Could not find BLAS/LAPACK libraries")
    endif(BLAS_LAPACK_LIBS)
  endif(BLAS_LAPACK_LIBS)
endif(APPLE)

# Attempt to detect the BLAS and LAPACK underscore conventions
include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LAPACK_LIBS})
check_function_exists(daxpy HAVE_DAXPY)
if(NOT HAVE_DAXPY)
  check_function_exists(_daxpy HAVE_DAXPY_PRE)
  if(HAVE_DAXPY_PRE)
    list(APPEND CMAKE_CXX_FLAGS "-DBLAS_PRE")
  else(HAVE_DAXPY_PRE)
    check_function_exists(daxpy_ HAVE_DAXPY_POST)
    if(HAVE_DAXPY_POST)
      list(APPEND CMAKE_CXX_FLAGS "-DBLAS_POST")
    else(HAVE_DAXPY_POST)
      check_function_exists(_daxpy_ HAVE_DAXPY_PRE_POST)
      if(HAVE_DAXPY_PRE_POST)
        list(APPEND CMAKE_CXX_FLAGS "-DBLAS_PRE")
        list(APPEND CMAKE_CXX_FLAGS "-DBLAS_POST")
      else(HAVE_DAXPY_PRE_POST)
        message(FATAL_ERROR "Could not determine BLAS format.")
      endif(HAVE_DAXPY_PRE_POST)
    endif(HAVE_DAXPY_POST)
  endif(HAVE_DAXPY_PRE)
endif(NOT HAVE_DAXPY)
check_function_exists(dpotrf HAVE_DPOTRF)
if(NOT HAVE_DPOTRF)
  check_function_exists(_dpotrf HAVE_DPOTRF_PRE)
  if(HAVE_DPOTRF_PRE)
    list(APPEND CMAKE_CXX_FLAGS "-DLAPACK_PRE")
  else(HAVE_DPOTRF_PRE)
    check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
    if(HAVE_DPOTRF_POST)
      list(APPEND CMAKE_CXX_FLAGS "-DLAPACK_POST")
    else(HAVE_DPOTRF_POST)
      check_function_exists(_dpotrf_ HAVE_DPOTRF_PRE_POST)
      if(HAVE_DPOTRF_PRE_POST)
        list(APPEND CMAKE_CXX_FLAGS "-DLAPACK_PRE")
        list(APPEND CMAKE_CXX_FLAGS "-DLAPACK_POST")
      else(HAVE_DPOTRF_PRE_POST)
        message(FATAL_ERROR "Could not determine LAPACK format.")
      endif(HAVE_DPOTRF_PRE_POST)
    endif(HAVE_DPOTRF_POST)
  endif(HAVE_DPOTRF_PRE)
endif(NOT HAVE_DPOTRF)

# Grab all of the .cpp and .hpp Elemental files
file(GLOB_RECURSE ELEMENTAL_CPP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "src/*.cpp")
file(GLOB_RECURSE ELEMENTAL_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "include/*.hpp")
set(ELEMENTAL_SRC "${ELEMENTAL_CPP};${ELEMENTAL_HEADERS}")

# Configure for PMRRR
if(BUILD_PMRRR)
  find_library(FOUND_PMRRR pmrrr)
  if(NOT FOUND_PMRRR)
    message("Could not find PMRRR; Elemental will not have eigensolvers.")
  endif(NOT FOUND_PMRRR)
endif(BUILD_PMRRR)

# Create the lists of test executables
if(BUILD_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/test)
  set(CORE_TESTS DistMatrix Matrix)
  set(BLAS_TESTS Gemm Hemm Her2k Herk Symm Symv Syr2k Syrk Trmm Trsm)
  set(LAPACK_TESTS Chol Hegst LU QR Tridiag Trinv UT)
endif(BUILD_TESTS)

# Set up for the build
include_directories("${PROJECT_SOURCE_DIR}/include")
set(ELEMENTAL_REQUIRED_LIBS ${BLAS_LAPACK_LIBS})
if(BUILD_PMRRR AND FOUND_PMRRR)
  list(APPEND ELEMENTAL_REQUIRED_LIBS FOUND_PMRRR)
  list(APPEND LAPACK_TESTS GeneralizedHermitianEig)
  list(APPEND LAPACK_TESTS HermitianEig)
else(BUILD_PMRRR AND FOUND_PMRRR)
  list(APPEND CMAKE_CXX_FLAGS "-DWITHOUT_PMRRR")
endif(BUILD_PMRRR AND FOUND_PMRRR)

# Build the pure-debug version (if necessary)
if(BUILD_DEBUG)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/pure/debug)
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/pure/release)
  add_library(pure-debug STATIC ${ELEMENTAL_SRC})
  set_target_properties(pure-debug PROPERTIES OUTPUT_NAME elemental)
  if(BUILD_TESTS)
    # Build the core tests
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/pure/debug/core)
    foreach(TEST ${CORE_TESTS})
      add_executable(pure-debug-core-${TEST} ${TEST_DIR}/core/${TEST}.cpp)
      target_link_libraries(pure-debug-core-${TEST} 
                            pure-debug ${ELEMENTAL_REQUIRED_LIBS})
      set_target_properties(pure-debug-core-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST})
    endforeach(TEST)
    # Build the blas tests
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/pure/debug/blas)
    foreach(TEST ${CORE_TESTS})
      add_executable(pure-debug-blas-${TEST} ${TEST_DIR}/blas/${TEST}.cpp)
      target_link_libraries(pure-debug-blas-${TEST} 
                            pure-debug ${ELEMENTAL_REQUIRED_LIBS})
      set_target_properties(pure-debug-blas-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST})
    endforeach(TEST)
    # Build the lapack tests
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/pure/debug/lapack)
    foreach(TEST ${CORE_TESTS})
      add_executable(pure-debug-lapack-${TEST} ${TEST_DIR}/lapack/${TEST}.cpp)
      target_link_libraries(pure-debug-lapack-${TEST} 
                            pure-debug ${ELEMENTAL_REQUIRED_LIBS})
      set_target_properties(pure-debug-lapack-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST})
    endforeach(TEST)
  endif(BUILD_TESTS)
endif(BUILD_DEBUG)

# Build the pure-release version
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/pure/release)
list(APPEND CMAKE_CXX_FLAGS "-DRELEASE")
add_library(pure-release STATIC ${ELEMENTAL_SRC})
set_target_properties(pure-release PROPERTIES OUTPUT_NAME elemental)
list(REMOVE_ITEM CMAKE_CXX_FLAGS "-DRELEASE")
if(BUILD_TESTS)
  # Build the core tests
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/pure/release/core) 
  foreach(TEST ${CORE_TESTS})
    add_executable(pure-release-core-${TEST} ${TEST_DIR}/core/${TEST}.cpp)
    target_link_libraries(pure-release-core-${TEST} 
                          pure-release ${ELEMENTAL_REQUIRED_LIBS})
    set_target_properties(pure-release-core-${TEST} 
                          PROPERTIES OUTPUT_NAME ${TEST})
  endforeach(TEST)
  # Build the blas tests
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/pure/release/blas) 
  foreach(TEST ${BLAS_TESTS})
    add_executable(pure-release-blas-${TEST} ${TEST_DIR}/blas/${TEST}.cpp) 
    target_link_libraries(pure-release-blas-${TEST} 
                          pure-release ${ELEMENTAL_REQUIRED_LIBS})
    set_target_properties(pure-release-blas-${TEST} 
                          PROPERTIES OUTPUT_NAME ${TEST})
  endforeach(TEST)
  # Build the lapack tests
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/pure/release/lapack)
  foreach(TEST ${LAPACK_TESTS})
    add_executable(pure-release-lapack-${TEST} ${TEST_DIR}/lapack/${TEST}.cpp)
    target_link_libraries(pure-release-lapack-${TEST} 
                          pure-release ${ELEMENTAL_REQUIRED_LIBS})
    set_target_properties(pure-release-lapack-${TEST} 
                          PROPERTIES OUTPUT_NAME ${TEST})
  endforeach(TEST)
endif(BUILD_TESTS)

if(BUILD_HYBRID AND OpenMP_FOUND)
  list(APPEND CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})
  # Build the hybrid-debug version (if necessary)
  if(BUILD_DEBUG)
    set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/hybrid/debug)
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/hybrid/debug)
    add_library(hybrid-debug STATIC ${ELEMENTAL_SRC})
    set_target_properties(hybrid-debug PROPERTIES OUTPUT_NAME elemental)
    if(BUILD_TESTS)
      # Build the core tests
      set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/hybrid/debug/core) 
      foreach(TEST ${CORE_TESTS})
        add_executable(hybrid-debug-core-${TEST} ${TEST_DIR}/core/${TEST}.cpp)
        target_link_libraries(hybrid-debug-core-${TEST} 
                              hybrid-debug ${ELEMENTAL_REQUIRED_LIBS})
        set_target_properties(hybrid-debug-core-${TEST} 
                              PROPERTIES OUTPUT_NAME ${TEST})
      endforeach(TEST)
      # Build the blas tests
      set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/hybrid/debug/blas) 
      foreach(TEST ${BLAS_TESTS})
        add_executable(hybrid-debug-blas-${TEST} ${TEST_DIR}/blas/${TEST}.cpp) 
        target_link_libraries(hybrid-debug-blas-${TEST} 
                              hybrid-debug ${ELEMENTAL_REQUIRED_LIBS})
        set_target_properties(hybrid-debug-blas-${TEST} 
                              PROPERTIES OUTPUT_NAME ${TEST})
      endforeach(TEST)
      # Build the lapack tests
      set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/hybrid/debug/lapack)
      foreach(TEST ${LAPACK_TESTS})
        add_executable(hybrid-debug-lapack-${TEST} 
                       ${TEST_DIR}/lapack/${TEST}.cpp)
        target_link_libraries(hybrid-debug-lapack-${TEST} 
                              hybrid-debug ${ELEMENTAL_REQUIRED_LIBS})
        set_target_properties(hybrid-debug-lapack-${TEST} 
                              PROPERTIES OUTPUT_NAME ${TEST})
      endforeach(TEST)
    endif(BUILD_TESTS)
  endif(BUILD_DEBUG)

  # Build the hybrid-release version
  set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/hybrid/release)
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/hybrid/release)
  list(APPEND CMAKE_CXX_FLAGS "-DRELEASE")
  add_library(hybrid-release STATIC ${ELEMENTAL_SRC})
  set_target_properties(hybrid-release PROPERTIES OUTPUT_NAME elemental)
  list(REMOVE_ITEM CMAKE_CXX_FLAGS "-DRELEASE")
  if(BUILD_TESTS)
    # Build the core tests
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/hybrid/release/core) 
    foreach(TEST ${CORE_TESTS})
      add_executable(hybrid-release-core-${TEST} ${TEST_DIR}/core/${TEST}.cpp)
      target_link_libraries(hybrid-release-core-${TEST} 
                            hybrid-release ${ELEMENTAL_REQUIRED_LIBS})
      set_target_properties(hybrid-release-core-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST})
    endforeach(TEST)
    # Build the blas tests
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/hybrid/release/blas) 
    foreach(TEST ${BLAS_TESTS})
      add_executable(hybrid-release-blas-${TEST} ${TEST_DIR}/blas/${TEST}.cpp) 
      target_link_libraries(hybrid-release-blas-${TEST} 
                            hybrid-release ${ELEMENTAL_REQUIRED_LIBS})
      set_target_properties(hybrid-release-blas-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST})
    endforeach(TEST)
    # Build the lapack tests
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/hybrid/release/lapack)
    foreach(TEST ${LAPACK_TESTS})
      add_executable(hybrid-release-lapack-${TEST} 
                     ${TEST_DIR}/lapack/${TEST}.cpp)
      target_link_libraries(hybrid-release-lapack-${TEST} 
                            hybrid-release ${ELEMENTAL_REQUIRED_LIBS})
      set_target_properties(hybrid-release-lapack-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST})
    endforeach(TEST)
  endif(BUILD_TESTS)
endif(BUILD_HYBRID AND OpenMP_FOUND)

