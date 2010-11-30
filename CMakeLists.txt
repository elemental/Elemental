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

if(CMAKE_BUILD_TYPE)
  message("Elemental does not use CMAKE_BUILD_TYPE.") 
endif(CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "" CACHE INTERNAL "Removing GUI option" FORCE)

option(BUILD_PURE    "Build a pure MPI version of Elemental"              ON)
option(BUILD_HYBRID  "Build a hybrid OpenMP/MPI version of Elemental"     OFF)
option(BUILD_DEBUG   "Build a version with many consistency checks"       OFF)
option(BUILD_RELEASE "Build the release (optimized) version of Elemental" ON)
option(BUILD_TESTS   "Build a collection of test executables"             ON)
option(LINK_PMRRR    "Link against the eigensolver PMRRR"                 ON)

option(WITHOUT_COMPLEX "Avoid all usage of complex numbers" OFF)
if(NOT WITHOUT_COMPLEX)
  option(AVOID_COMPLEX_MPI "Avoid complex MPI routines for robustness" ON)
endif(NOT WITHOUT_COMPLEX)
option(TIMING "If defined, some routines will accumulate statistics"        ON)
option(CACHE_WARNINGS "Warns when using cache-unfriendly routines"          OFF)
option(UNALIGNED_WARNINGS "Warn when performing unaligned redistributions"  OFF)
option(VECTOR_WARNINGS "Warn when vector redistribution chances are missed" OFF)
option(ENABLE_ALL_DISTRIBUTED_DOT "Build every parallel dot product."       OFF)

if(BUILD_PURE AND BUILD_DEBUG)
  set(CXX_PURE_DEBUG_FLAGS "-O2 -g" CACHE STRING
      "CXX optimization/debugging flags for pure MPI debug build")
endif(BUILD_PURE AND BUILD_DEBUG)
if(BUILD_PURE AND BUILD_RELEASE)
  set(CXX_PURE_RELEASE_FLAGS "-O3" CACHE STRING
      "CXX optimization flags for pure MPI release build")
endif(BUILD_PURE AND BUILD_RELEASE)
if(BUILD_HYBRID AND BUILD_DEBUG)
  set(CXX_HYBRID_DEBUG_FLAGS "-O2 -g" CACHE STRING
      "CXX optimization/debug flags for hybrid OpenMP/MPI debug build")
endif(BUILD_HYBRID AND BUILD_DEBUG)
if(BUILD_HYBRID AND BUILD_RELEASE)
  set(CXX_HYBRID_RELEASE_FLAGS "-O3" CACHE STRING
      "CXX optimization flags for hybrid OpenMP/MPI release build")
endif(BUILD_HYBRID AND BUILD_RELEASE)

if(NOT MPI_CXX_COMPILER)
  find_package(MPI)
  set(MPI_CXX_DESC "Full path to MPI CXX compiler. If blank, CMake will search for it.")
  if(NOT MPI_FOUND)
    set(MPI_CXX_COMPILER "" CACHE STRING ${MPI_CXX_MSG})
    message(FATAL_ERROR "Could not find MPI CXX compiler. Please manually set MPI_CXX_COMPILER")
  endif(NOT MPI_FOUND)
  # Regrettably, 'find_package(MPI)' does not qualify the compiler as CXX
  set(MPI_CXX_COMPILER ${MPI_COMPILER} CACHE STRING ${MPI_CXX_DESC})
endif(NOT MPI_CXX_COMPILER)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

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
if(NOT BLAS_LAPACK_LIBS)
  set(BLAS_LAPACK_DESC "BLAS/LAPACK link flags")
  if(APPLE)
    set(BLAS_LAPACK_LIBS "-framework vecLib" CACHE STRING ${BLAS_LAPACK_DESC})
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
      find_library(PTHREAD_LIB pthread DOC "Pthreads library")
      if(PTHREAD_LIB)
        set(FOUND_BLAS "MKL")    
        set(FOUND_LAPACK "MKL")
        set(BLAS_LAPACK_LIBS ${BLAS_LAPACK_LIBS} 
            CACHE STRING ${BLAS_LAPACK_DESC})
      else(PTHREAD_LIB)
        message("Found MKL but could not find pthread library")  
        set(BLAS_LAPACK_LIBS "")
      endif(PTHREAD_LIB)
    endif(BLAS_LAPACK_LIBS)
    if(NOT BLAS_LAPACK_LIBS)
      # Look for default BLAS and LAPACK
      set(REFERENCE_LIBS lapack blas)
      find_library_list(BLAS_LAPACK_LIBS REFERENCE_LIBS "" "")
      if(BLAS_LAPACK_LIBS)
        set(FOUND_BLAS "Default")
        set(FOUND_LAPACK "Default")
        set(BLAS_LAPACK_LIBS ${BLAS_LAPACK_LIBS}
            CACHE STRING ${BLAS_LAPACK_DESC})
      else(BLAS_LAPACK_LIBS)
        set(BLAS_LAPACK_LIBS "" CACHE STRING ${BLAS_LAPACK_DESC})
        message(FATAL_ERROR 
          "Could not find BLAS/LAPACK libs. Please manually specify BLAS_LAPACK_LIBS.")
      endif(BLAS_LAPACK_LIBS)
    endif(BLAS_LAPACK_LIBS)
  endif(APPLE)
endif(NOT BLAS_LAPACK_LIBS)

# Begin building CXX_BASE_FLAGS
set(CXX_BASE_FLAGS)
if(WITHOUT_COMPLEX)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DWITHOUT_COMPLEX")
elseif(AVOID_COMPLEX_MPI)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DAVOID_COMPLEX_MPI")
endif(WITHOUT_COMPLEX)
if(TIMING)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DTIMING")
endif(TIMING)
if(CACHE_WARNINGS)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DCACHE_WARNINGS")
endif(CACHE_WARNINGS)
if(UNALIGNED_WARNINGS)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DUNALIGNED_WARNINGS")
endif(UNALIGNED_WARNINGS)
if(VECTOR_WARNINGS)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DVECTOR_WARNINGS")
endif(VECTOR_WARNINGS)
if(ENABLE_ALL_DISTRIBUTED_DOT)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DENABLE_ALL_DISTRIBUTED_DOT")
endif(ENABLE_ALL_DISTRIBUTED_DOT)

# Attempt to detect the BLAS and LAPACK underscore conventions
include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LAPACK_LIBS})
check_function_exists(daxpy HAVE_DAXPY)
if(NOT HAVE_DAXPY)
  check_function_exists(_daxpy HAVE_DAXPY_PRE)
  if(HAVE_DAXPY_PRE)
    set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DBLAS_PRE")
  else(HAVE_DAXPY_PRE)
    check_function_exists(daxpy_ HAVE_DAXPY_POST)
    if(HAVE_DAXPY_POST)
      set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DBLAS_POST")
    else(HAVE_DAXPY_POST)
      check_function_exists(_daxpy_ HAVE_DAXPY_PRE_POST)
      if(HAVE_DAXPY_PRE_POST)
        set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DBLAS_PRE -DBLAS_POST")
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
    set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DLAPACK_PRE")
  else(HAVE_DPOTRF_PRE)
    check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
    if(HAVE_DPOTRF_POST)
      set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DLAPACK_PRE")
    else(HAVE_DPOTRF_POST)
      check_function_exists(_dpotrf_ HAVE_DPOTRF_PRE_POST)
      if(HAVE_DPOTRF_PRE_POST)
        set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DLAPACK_PRE -DLAPACK_POST")
      else(HAVE_DPOTRF_PRE_POST)
        message(FATAL_ERROR "Could not determine LAPACK format.")
      endif(HAVE_DPOTRF_PRE_POST)
    endif(HAVE_DPOTRF_POST)
  endif(HAVE_DPOTRF_PRE)
endif(NOT HAVE_DPOTRF)

# Search for PMRRR
if(LINK_PMRRR)
  find_library(PMRRR_LIB pmrrr 
      DOC "Parallel Multiple Relatively Robust Representations (PMRRR) library")
  if(NOT PMRRR_LIB)
    set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DWITHOUT_PMRRR")
    message("Could not find PMRRR; Elemental will not have eigensolvers.")
  endif(NOT PMRRR_LIB)
endif(LINK_PMRRR)

# Grab all of the .cpp and .hpp Elemental files
file(GLOB_RECURSE ELEMENTAL_CPP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "src/*.cpp")
file(GLOB_RECURSE ELEMENTAL_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "include/*.hpp")
if(NOT LINK_PMRRR OR NOT PMRRR_LIB)
  # It might be cleaner to force these to be empty through the preprocessor
  list(REMOVE_ITEM ELEMENTAL_CPP 
       "src/lapack/HermitianEig/HermitianEig.cpp")
  list(REMOVE_ITEM ELEMENTAL_CPP 
       "src/lapack/GeneralizedHermitianEig/GeneralizedHermitianEig.cpp")
endif(NOT LINK_PMRRR OR NOT PMRRR_LIB)
set(ELEMENTAL_SRC "${ELEMENTAL_CPP};${ELEMENTAL_HEADERS}")

# Create the lists of test executables
if(BUILD_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/test)
  set(CORE_TESTS DistMatrix Matrix)
  set(BLAS_TESTS Gemm Hemm Her2k Herk Symm Symv Syr2k Syrk Trmm Trsm)
  set(LAPACK_TESTS Chol Hegst LU QR Tridiag Trinv UT)
  if(LINK_PMRRR AND PMRRR_LIB)
    list(APPEND LAPACK_TESTS GeneralizedHermitianEig)
    list(APPEND LAPACK_TESTS HermitianEig)
  endif(LINK_PMRRR AND PMRRR_LIB)
endif(BUILD_TESTS)

# Set up the include directories and the required libraries
include_directories("${PROJECT_SOURCE_DIR}/include")
set(ELEMENTAL_REQUIRED_LIBS ${BLAS_LAPACK_LIBS})
if(LINK_PMRRR AND PMRRR_LIB)
  list(APPEND ELEMENTAL_REQUIRED_LIBS ${PMRRR_LIB})
endif(LINK_PMRRR AND PMRRR_LIB)

# Build the pure-debug version (if necessary)
if(BUILD_PURE AND BUILD_DEBUG)
  set(CMAKE_CXX_FLAGS "${CXX_BASE_FLAGS} ${CXX_PURE_DEBUG_FLAGS}")
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
    foreach(TEST ${BLASE_TESTS})
      add_executable(pure-debug-blas-${TEST} ${TEST_DIR}/blas/${TEST}.cpp)
      target_link_libraries(pure-debug-blas-${TEST} 
                            pure-debug ${ELEMENTAL_REQUIRED_LIBS})
      set_target_properties(pure-debug-blas-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST})
    endforeach(TEST)
    # Build the lapack tests
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/pure/debug/lapack)
    foreach(TEST ${LAPACK_TESTS})
      add_executable(pure-debug-lapack-${TEST} ${TEST_DIR}/lapack/${TEST}.cpp)
      target_link_libraries(pure-debug-lapack-${TEST} 
                            pure-debug ${ELEMENTAL_REQUIRED_LIBS})
      set_target_properties(pure-debug-lapack-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST})
    endforeach(TEST)
  endif(BUILD_TESTS)
endif(BUILD_PURE AND BUILD_DEBUG)

# Build the pure-release version (if necessary)
if(BUILD_PURE AND BUILD_RELEASE)
  set(CMAKE_CXX_FLAGS "${CXX_BASE_FLAGS} ${CXX_PURE_RELEASE_FLAGS} -DRELEASE")
  set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/pure/release)
  add_library(pure-release STATIC ${ELEMENTAL_SRC})
  set_target_properties(pure-release PROPERTIES OUTPUT_NAME elemental)
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
endif(BUILD_PURE AND BUILD_RELEASE)

if(BUILD_HYBRID AND OpenMP_FOUND)
  # Build the hybrid-debug version (if necessary)
  if(BUILD_DEBUG)
    set(CMAKE_CXX_FLAGS "${CXX_BASE_FLAGS} ${CXX_HYBRID_DEBUG_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
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

  # Build the hybrid-release version (if necessary)
  if(BUILD_RELEASE)
    set(CMAKE_CXX_FLAGS "${CXX_BASE_FLAGS} ${CXX_HYBRID_RELEASE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DRELEASE")
    set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/hybrid/release)
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin/hybrid/release)
    add_library(hybrid-release STATIC ${ELEMENTAL_SRC})
    set_target_properties(hybrid-release PROPERTIES OUTPUT_NAME elemental)
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
  endif(BUILD_RELEASE)
endif(BUILD_HYBRID AND OpenMP_FOUND)

