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
set(Elemental_VERSION_MINOR 52)

if(CMAKE_BUILD_TYPE)
  message("Elemental does not use CMAKE_BUILD_TYPE.") 
endif(CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "" CACHE INTERNAL "Removing GUI option" FORCE)

option(BUILD_PURE    "Build a pure MPI version of Elemental"              ON)
option(BUILD_HYBRID  "Build a hybrid OpenMP/MPI version of Elemental"     ON)
option(BUILD_DEBUG   "Build a version with many consistency checks"       OFF)
option(BUILD_RELEASE "Build the release (optimized) version of Elemental" ON)
option(BUILD_TESTS   "Build a collection of test executables"             ON)
option(BUILD_PMRRR   "Build and link against eigensolver PMRRR"           ON)

option(WITHOUT_COMPLEX "Avoid all usage of complex numbers" OFF)
if(NOT WITHOUT_COMPLEX)
  option(AVOID_COMPLEX_MPI "Avoid complex MPI routines for robustness" ON)
endif(NOT WITHOUT_COMPLEX)
option(TIMING "If defined, some routines will accumulate statistics"        ON)
option(CACHE_WARNINGS "Warns when using cache-unfriendly routines"          OFF)
option(UNALIGNED_WARNINGS "Warn when performing unaligned redistributions"  OFF)
option(VECTOR_WARNINGS "Warn when vector redistribution chances are missed" OFF)
option(ENABLE_ALL_DISTRIBUTED_DOT "Build every parallel dot product."       OFF)
MARK_AS_ADVANCED(AVOID_COMPLEX_MPI TIMING CACHE_WARNINGS UNALIGNED_WARNINGS 
                 VECTOR_WARNINGS ENABLE_ALL_DISTRIBUTED_DOT)

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

if(MPI_CXX_COMPILER AND NOT EXISTS ${MPI_CXX_COMPILER})
  message(FATAL_ERROR 
          "Invalid MPI_CXX_COMPILER: ${MPI_CXX_COMPILER} does not exist.")
endif(MPI_CXX_COMPILER AND NOT EXISTS ${MPI_CXX_COMPILER})
if(NOT MPI_CXX_COMPILER)
  find_package(MPI)
  set(MPI_CXX_DESC 
      "Full path to MPI CXX compiler. If blank, CMake will search for it.")
  if(NOT MPI_FOUND)
    # Try a manual search for the common MPI names
    find_program(MPI_COMPILER mpicxx mpic++ mpixlcxx_r mpixlcxx mpixlC mpixlC_r 
                 DOC "Full path to MPI CXX compiler")
    if(NOT MPI_COMPILER)
      set(MPI_CXX_COMPILER "" CACHE STRING ${MPI_CXX_MSG})
      message(FATAL_ERROR 
        "Could not find MPI CXX compiler. Please manually set MPI_CXX_COMPILER")
    endif(NOT MPI_COMPILER)
  endif(NOT MPI_FOUND)
  # Regrettably, 'find_package(MPI)' does not qualify the compiler as CXX
  set(MPI_CXX_COMPILER ${MPI_COMPILER} CACHE STRING ${MPI_CXX_DESC})
endif(NOT MPI_CXX_COMPILER)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

if(BUILD_HYBRID)
  if(MANUAL_OPENMP_CXX_FLAGS)
    set(HAVE_OPENMP TRUE)
    set(OpenMP_CXX_FLAGS ${MANUAL_OPENMP_CXX_FLAGS})
  else(MANUAL_OPENMP_CXX_FLAGS)
    find_package(OpenMP)
    if(OPENMP_FOUND)
      set(HAVE_OPENMP TRUE)
    else(OPENMP_FOUND)
      set(MANUAL_OPENMP_CXX_FLAGS "" CACHE STRING "OpenMP CXX FLAGS")
      message("Hybrid build disabled because OpenMP support was not verified.")
      message("You can alternatively manually specify MANUAL_OPENMP_CXX_FLAGS.")
    endif(OPENMP_FOUND)
  endif(MANUAL_OPENMP_CXX_FLAGS)
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
    # On some machines, all of the above FOUND_${LIB} variables show up in 
    # the CMake wizard's cache. This is annoying, so we should explicitly 
    # remove them by forcing them to be internal.
    set(FOUND_${LIB} CACHE INTERNAL "")
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
if(BUILD_PURE OR BUILD_PMRRR)
  if(PURE_BLAS_LAPACK_LIBS)
    set(NEED_PURE_BLAS_LAPACK FALSE)
  else(PURE_BLAS_LAPACK_LIBS)
    set(NEED_PURE_BLAS_LAPACK TRUE)
  endif(PURE_BLAS_LAPACK_LIBS)
else(BUILD_PURE OR BUILD_PMRRR)
  set(NEED_PURE_BLAS_LAPACK FALSE)
endif(BUILD_PURE OR BUILD_PMRRR)
if(BUILD_HYBRID AND HAVE_OPENMP AND NOT HYBRID_BLAS_LAPACK_LIBS)
  set(NEED_HYBRID_BLAS_LAPACK TRUE)
else(BUILD_HYBRID AND HAVE_OPENMP AND NOT HYBRID_BLAS_LAPACK_LIBS)
  set(NEED_HYBRID_BLAS_LAPACK FALSE)
endif(BUILD_HYBRID AND HAVE_OPENMP AND NOT HYBRID_BLAS_LAPACK_LIBS)
if(NEED_PURE_BLAS_LAPACK OR NEED_HYBRID_BLAS_LAPACK)
  set(PURE_BLAS_LAPACK_DESC "Serial BLAS/LAPACK link flags")
  set(HYBRID_BLAS_LAPACK_DESC "Threaded BLAS/LAPACK link flags")
  if(APPLE)
    if(NEED_PURE_BLAS_LAPACK)
      set(PURE_BLAS_LAPACK_LIBS "-framework vecLib" 
          CACHE STRING ${PURE_BLAS_LAPACK_DESC})
    endif(NEED_PURE_BLAS_LAPACK)
    if(NEED_HYBRID_BLAS_LAPACK)
      set(HYBRID_BLAS_LAPACK_LIBS "-framework vecLib" 
          CACHE STRING ${HYBRID_BLAS_LAPACK_DESC})
    endif(NEED_HYBRID_BLAS_LAPACK)
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(FOUND_BLAS_LAPACK_LIBS FALSE)
    # Look for MKL first
    if(32_BIT_SYSTEM)
      set(MKL_SUBDIR lib/32)
    else(32_BIT_SYSTEM)
      set(MKL_SUBDIR lib/em64t)
    endif(32_BIT_SYSTEM)
    set(MKL_LIBS mkl_intel_lp64 mkl_intel_thread mkl_core guide iomp5)
    find_library_list(MKL_BLAS_LAPACK_LIBS MKL_LIBS "" ${MKL_SUBDIR})
    if(MKL_BLAS_LAPACK_LIBS)
      find_library(PTHREAD_LIB pthread DOC "Pthreads library")
      if(PTHREAD_LIB)
        if(NEED_PURE_BLAS_LAPACK)
          set(PURE_BLAS_LAPACK_LIBS ${MKL_BLAS_LAPACK_LIBS} 
              CACHE STRING ${PURE_BLAS_LAPACK_DESC})
          message(STATUS "Using MKL for serial BLAS/LAPACK.")
        endif(NEED_PURE_BLAS_LAPACK)
        if(NEED_HYBRID_BLAS_LAPACK)
          set(HYBRID_BLAS_LAPACK_LIBS ${MKL_BLAS_LAPACK_LIBS} 
              CACHE STRING ${HYBRID_BLAS_LAPACK_DESC})
          message(STATUS "Using MKL for threaded BLAS/LAPACK.")
        endif(NEED_HYBRID_BLAS_LAPACK)
        set(FOUND_BLAS_LAPACK_LIBS TRUE)
      else(PTHREAD_LIB)
        message(STATUS "WARNING: Found MKL but could not find pthread lib")  
        if(NEED_PURE_BLAS_LAPACK)
          set(PURE_BLAS_LAPACK_LIBS "" 
              CACHE STRING ${PURE_BLAS_LAPACK_DESC})
        endif(NEED_PURE_BLAS_LAPACK)
        if(NEED_HYBRID_BLAS_LAPACK)
          set(HYBRID_BLAS_LAPACK_LIBS "" 
              CACHE STRING ${HYBRID_BLAS_LAPACK_DESC})
        endif(NEED_HYBRID_BLAS_LAPACK)
      endif(PTHREAD_LIB)
    endif(MKL_BLAS_LAPACK_LIBS)
    if(NOT FOUND_BLAS_LAPACK_LIBS)
      # Look for default BLAS and LAPACK
      set(REFERENCE_LIBS lapack blas)
      find_library_list(REF_BLAS_LAPACK_LIBS REFERENCE_LIBS "" "")
      if(REF_BLAS_LAPACK_LIBS)
        if(NEED_PURE_BLAS_LAPACK)
          set(PURE_BLAS_LAPACK_LIBS ${REF_BLAS_LAPACK_LIBS}
              CACHE STRING ${PURE_BLAS_LAPACK_DESC})
          message(STATUS "WARNING: Using reference serial BLAS/LAPACK.")
        endif(NEED_PURE_BLAS_LAPACK)
        if(NEED_HYBRID_BLAS_LAPACK)
          set(HYBRID_BLAS_LAPACK_LIBS ${REF_BLAS_LAPACK_LIBS}
              CACHE STRING ${HYBRID_BLAS_LAPACK_DESC})
          message(STATUS "WARNING: Using reference threaded BLAS/LAPACK.")
        endif(NEED_HYBRID_BLAS_LAPACK)
      else(REF_BLAS_LAPACK_LIBS)
        if(NEED_PURE_BLAS_LAPACK)
          set(PURE_BLAS_LAPACK_LIBS "" 
              CACHE STRING ${PURE_BLAS_LAPACK_DESC})
        endif(NEED_PURE_BLAS_LAPACK)
        if(NEED_HYBRID_BLAS_LAPACK)
          set(HYBRID_BLAS_LAPACK_LIBS "" 
              CACHE STRING ${HYBRID_BLAS_LAPACK_DESC})
        endif(NEED_HYBRID_BLAS_LAPACK)
        message(FATAL_ERROR 
          "Could not find BLAS/LAPACK libs. Please manually specify BLAS_LAPACK_LIBS.")
      endif(REF_BLAS_LAPACK_LIBS)
    endif(NOT FOUND_BLAS_LAPACK_LIBS)
  endif(APPLE)
endif(NEED_PURE_BLAS_LAPACK OR NEED_HYBRID_BLAS_LAPACK)

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
if(BUILD_PURE)
  set(CMAKE_REQUIRED_LIBRARIES ${PURE_BLAS_LAPACK_LIBS})
  check_function_exists(daxpy HAVE_DAXPY)
  if(HAVE_DAXPY)
    set(PURE_BLAS_DEFS "")
  else(HAVE_DAXPY)
    check_function_exists(_daxpy HAVE_DAXPY_PRE)
    if(HAVE_DAXPY_PRE)
      set(PURE_BLAS_DEFS "-DBLAS_PRE")
    else(HAVE_DAXPY_PRE)
      check_function_exists(daxpy_ HAVE_DAXPY_POST)
      if(HAVE_DAXPY_POST)
        set(PURE_BLAS_DEFS "-DBLAS_POST")
      else(HAVE_DAXPY_POST)
        check_function_exists(_daxpy_ HAVE_DAXPY_PRE_POST)
        if(HAVE_DAXPY_PRE_POST)
          set(PURE_BLAS_DEFS "-DBLAS_PRE -DBLAS_POST")
        else(HAVE_DAXPY_PRE_POST)
          message(FATAL_ERROR "Could not determine BLAS format.")
        endif(HAVE_DAXPY_PRE_POST)
      endif(HAVE_DAXPY_POST)
    endif(HAVE_DAXPY_PRE)
  endif(HAVE_DAXPY)
  check_function_exists(dpotrf HAVE_DPOTRF)
  if(HAVE_DPOTRF)
    set(PURE_LAPACK_DEFS "")
  else(HAVE_DPOTRF)
    check_function_exists(_dpotrf HAVE_DPOTRF_PRE)
    if(HAVE_DPOTRF_PRE)
      set(PURE_LAPACK_DEFS "-DLAPACK_PRE")
    else(HAVE_DPOTRF_PRE)
      check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
      if(HAVE_DPOTRF_POST)
        set(PURE_LAPACK_DEFS "-DLAPACK_POST")
      else(HAVE_DPOTRF_POST)
        check_function_exists(_dpotrf_ HAVE_DPOTRF_PRE_POST)
        if(HAVE_DPOTRF_PRE_POST)
          set(PURE_LAPACK_DEFS "-DLAPACK_PRE -DLAPACK_POST")
        else(HAVE_DPOTRF_PRE_POST)
          message(FATAL_ERROR "Could not determine LAPACK format.")
        endif(HAVE_DPOTRF_PRE_POST)
      endif(HAVE_DPOTRF_POST)
    endif(HAVE_DPOTRF_PRE)
  endif(HAVE_DPOTRF)
endif(BUILD_PURE)
if(BUILD_HYBRID AND HAVE_OPENMP)
  set(CMAKE_REQUIRED_LIBRARIES ${HYBRID_BLAS_LAPACK_LIBS})
  check_function_exists(daxpy HAVE_DAXPY)
  if(HAVE_DAXPY)
    set(HYBRID_BLAS_DEFS "")
  else(HAVE_DAXPY)
    check_function_exists(_daxpy HAVE_DAXPY_PRE)
    if(HAVE_DAXPY_PRE)
      set(HYBRID_BLAS_DEFS "-DBLAS_PRE")
    else(HAVE_DAXPY_PRE)
      check_function_exists(daxpy_ HAVE_DAXPY_POST)
      if(HAVE_DAXPY_POST)
        set(HYBRID_BLAS_DEFS "-DBLAS_POST")
      else(HAVE_DAXPY_POST)
        check_function_exists(_daxpy_ HAVE_DAXPY_PRE_POST)
        if(HAVE_DAXPY_PRE_POST)
          set(HYBRID_BLAS_DEFS "-DBLAS_PRE -DBLAS_POST")
        else(HAVE_DAXPY_PRE_POST)
          message(FATAL_ERROR "Could not determine BLAS format.")
        endif(HAVE_DAXPY_PRE_POST)
      endif(HAVE_DAXPY_POST)
    endif(HAVE_DAXPY_PRE)
  endif(HAVE_DAXPY)
  check_function_exists(dpotrf HAVE_DPOTRF)
  if(HAVE_DPOTRF)
    set(HYBRID_LAPACK_DEFS "")
  else(HAVE_DPOTRF)
    check_function_exists(_dpotrf HAVE_DPOTRF_PRE)
    if(HAVE_DPOTRF_PRE)
      set(HYBRID_LAPACK_DEFS "-DLAPACK_PRE")
    else(HAVE_DPOTRF_PRE)
      check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
      if(HAVE_DPOTRF_POST)
        set(HYBRID_LAPACK_DEFS "-DLAPACK_POST")
      else(HAVE_DPOTRF_POST)
        check_function_exists(_dpotrf_ HAVE_DPOTRF_PRE_POST)
        if(HAVE_DPOTRF_PRE_POST)
          set(HYBRID_LAPACK_DEFS "-DLAPACK_PRE -DLAPACK_POST")
        else(HAVE_DPOTRF_PRE_POST)
          message(FATAL_ERROR "Could not determine LAPACK format.")
        endif(HAVE_DPOTRF_PRE_POST)
      endif(HAVE_DPOTRF_POST)
    endif(HAVE_DPOTRF_PRE)
  endif(HAVE_DPOTRF)
endif(BUILD_HYBRID AND HAVE_OPENMP)

if(BUILD_PMRRR)
  add_subdirectory(contrib/pmrrr)
endif(BUILD_PMRRR)

# Grab all of the .cpp and .hpp Elemental files
file(GLOB_RECURSE ELEMENTAL_CPP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "src/*.cpp")
file(GLOB_RECURSE ELEMENTAL_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "include/*.hpp")
set(ELEMENTAL_SRC "${ELEMENTAL_CPP};${ELEMENTAL_HEADERS}")

# Create the lists of test executables
if(BUILD_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/test)
  set(CORE_TESTS DistMatrix Matrix)
  set(BLAS_TESTS Gemm Hemm Her2k Herk Symm Symv Syr2k Syrk Trmm Trsm)
  set(LAPACK_TESTS_WITHOUT_PMRRR Chol Hegst LU QR Tridiag Trinv UT)

  set(LAPACK_TESTS_WITH_PMRRR ${LAPACK_TESTS_WITHOUT_PMRRR})
  list(APPEND LAPACK_TESTS_WITH_PMRRR GeneralizedHermitianEig)
  list(APPEND LAPACK_TESTS_WITH_PMRRR HermitianEig)
endif(BUILD_TESTS)

# Make sure the Elemental headers can be found
include_directories("${PROJECT_SOURCE_DIR}/include")

# Unroll the BLAS/LAPACK libs into a string
set(PURE_BLAS_LAPACK_FLAGS "")
foreach(LIB ${PURE_BLAS_LAPACK_LIBS})
  set(PURE_BLAS_LAPACK_FLAGS "${PURE_BLAS_LAPACK_FLAGS} ${LIB}")
endforeach(LIB)
set(HYBRID_BLAS_LAPACK_FLAGS "")
foreach(LIB ${HYBRID_BLAS_LAPACK_LIBS})
  set(HYBRID_BLAS_LAPACK_FLAGS "${HYBRID_BLAS_LAPACK_FLAGS} ${LIB}")
endforeach(LIB)

# Build the pure-debug version (if necessary)
if(BUILD_PURE AND BUILD_DEBUG)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/lib/pure/debug")
  if(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
    set(CMAKE_CXX_FLAGS 
        "${CXX_BASE_FLAGS} ${PURE_BLAS_DEFS} ${PURE_LAPACK_DEFS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_PURE_DEBUG_FLAGS}")
    add_library(pure-debug STATIC ${ELEMENTAL_SRC})
    target_link_libraries(pure-debug pmrrr-pure)
  else(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
    set(CMAKE_CXX_FLAGS 
        "${CXX_BASE_FLAGS} ${PURE_BLAS_DEFS} ${PURE_LAPACK_DEFS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_PURE_DEBUG_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITHOUT_PMRRR")
    add_library(pure-debug STATIC ${ELEMENTAL_SRC})
  endif(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
  set_target_properties(pure-debug 
                        PROPERTIES OUTPUT_NAME elemental
                        RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
  if(BUILD_TESTS)
    # Build the core tests
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/pure/debug/core")
    foreach(TEST ${CORE_TESTS})
      add_executable(pure-debug-core-${TEST} ${TEST_DIR}/core/${TEST}.cpp)
      target_link_libraries(pure-debug-core-${TEST} pure-debug)
      set_target_properties(pure-debug-core-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST}
                            RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                            LINK_FLAGS "${PURE_BLAS_LAPACK_FLAGS}")
    endforeach(TEST)
    # Build the blas tests
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/pure/debug/blas")
    foreach(TEST ${BLAS_TESTS})
      add_executable(pure-debug-blas-${TEST} ${TEST_DIR}/blas/${TEST}.cpp)
      target_link_libraries(pure-debug-blas-${TEST} pure-debug)
      set_target_properties(pure-debug-blas-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST}
                            RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                            LINK_FLAGS "${PURE_BLAS_LAPACK_FLAGS}")
    endforeach(TEST)
    # Build the lapack tests
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/pure/debug/lapack")
    if(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
      set(LAPACK_TESTS ${LAPACK_TESTS_WITH_PMRRR})
    else(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
      set(LAPACK_TESTS ${LAPACK_TESTS_WITHOUT_PMRRR})
    endif(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
    foreach(TEST ${LAPACK_TESTS})
      add_executable(pure-debug-lapack-${TEST} ${TEST_DIR}/lapack/${TEST}.cpp)
      target_link_libraries(pure-debug-lapack-${TEST} pure-debug)
      set_target_properties(pure-debug-lapack-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST}
                            RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                            LINK_FLAGS "${PURE_BLAS_LAPACK_FLAGS}")
    endforeach(TEST)
  endif(BUILD_TESTS)
endif(BUILD_PURE AND BUILD_DEBUG)

# Build the pure-release version (if necessary)
if(BUILD_PURE AND BUILD_RELEASE)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/lib/pure/release")
  if(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
    set(CMAKE_CXX_FLAGS 
        "${CXX_BASE_FLAGS} ${PURE_BLAS_DEFS} ${PURE_LAPACK_DEFS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_PURE_RELEASE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRELEASE")
    add_library(pure-release STATIC ${ELEMENTAL_SRC})
    target_link_libraries(pure-release pmrrr-pure)
  else(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
    set(CMAKE_CXX_FLAGS 
        "${CXX_BASE_FLAGS} ${PURE_BLAS_DEFS} ${PURE_LAPACK_DEFS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_PURE_RELEASE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRELEASE -DWITHOUT_PMRRR")
    add_library(pure-release STATIC ${ELEMENTAL_SRC})
  endif(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
  set_target_properties(pure-release
                        PROPERTIES OUTPUT_NAME elemental
                        RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
  if(BUILD_TESTS)
    # Build the core tests
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/pure/release/core")
    foreach(TEST ${CORE_TESTS})
      add_executable(pure-release-core-${TEST} ${TEST_DIR}/core/${TEST}.cpp)
      target_link_libraries(pure-release-core-${TEST} pure-release)
      set_target_properties(pure-release-core-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST}
                            RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                            LINK_FLAGS "${PURE_BLAS_LAPACK_FLAGS}")
    endforeach(TEST)
    # Build the blas tests
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/pure/release/blas")
    foreach(TEST ${BLAS_TESTS})
      add_executable(pure-release-blas-${TEST} ${TEST_DIR}/blas/${TEST}.cpp) 
      target_link_libraries(pure-release-blas-${TEST} pure-release)
      set_target_properties(pure-release-blas-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST}
                            RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                            LINK_FLAGS "${PURE_BLAS_LAPACK_FLAGS}")
    endforeach(TEST)
    # Build the lapack tests
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/pure/release/lapack")
    if(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
      set(LAPACK_TESTS ${LAPACK_TESTS_WITH_PMRRR})
    else(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
      set(LAPACK_TESTS ${LAPACK_TESTS_WITHOUT_PMRRR})
    endif(BUILD_PMRRR AND NOT FAILED_PURE_PMRRR)
    foreach(TEST ${LAPACK_TESTS})
      add_executable(pure-release-lapack-${TEST} ${TEST_DIR}/lapack/${TEST}.cpp)
      target_link_libraries(pure-release-lapack-${TEST} pure-release)
      set_target_properties(pure-release-lapack-${TEST} 
                            PROPERTIES OUTPUT_NAME ${TEST}
                            RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                            LINK_FLAGS "${PURE_BLAS_LAPACK_FLAGS}")
    endforeach(TEST)
  endif(BUILD_TESTS)
endif(BUILD_PURE AND BUILD_RELEASE)

if(BUILD_HYBRID AND HAVE_OPENMP)
  # Build the hybrid-debug version (if necessary)
  if(BUILD_DEBUG)
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/lib/hybrid/debug")
    if(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
      set(CMAKE_CXX_FLAGS 
          "${CXX_BASE_FLAGS} ${HYBRID_BLAS_DEFS} ${HYBRID_LAPACK_DEFS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_HYBRID_DEBUG_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      add_library(hybrid-debug STATIC ${ELEMENTAL_SRC})
      target_link_libraries(hybrid-debug pmrrr-hybrid)
    else(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
      set(CMAKE_CXX_FLAGS 
          "${CXX_BASE_FLAGS} ${HYBRID_BLAS_DEFS} ${HYBRID_LAPACK_DEFS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_HYBRID_DEBUG_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITHOUT_PMRRR")
      add_library(hybrid-debug STATIC ${ELEMENTAL_SRC})
    endif(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
    set_target_properties(hybrid-debug
                          PROPERTIES OUTPUT_NAME elemental
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(BUILD_TESTS)
      # Build the core tests
      set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/hybrid/debug/core")
      foreach(TEST ${CORE_TESTS})
        add_executable(hybrid-debug-core-${TEST} ${TEST_DIR}/core/${TEST}.cpp)
        target_link_libraries(hybrid-debug-core-${TEST} hybrid-debug)
        set_target_properties(hybrid-debug-core-${TEST} 
                              PROPERTIES OUTPUT_NAME ${TEST}
                              RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                              LINK_FLAGS "${HYBRID_BLAS_LAPACK_FLAGS}")
      endforeach(TEST)
      # Build the blas tests
      set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/hybrid/debug/blas")
      foreach(TEST ${BLAS_TESTS})
        add_executable(hybrid-debug-blas-${TEST} ${TEST_DIR}/blas/${TEST}.cpp) 
        target_link_libraries(hybrid-debug-blas-${TEST} hybrid-debug)
        set_target_properties(hybrid-debug-blas-${TEST} 
                              PROPERTIES OUTPUT_NAME ${TEST}
                              RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                              LINK_FLAGS "${HYBRID_BLAS_LAPACK_FLAGS}")
      endforeach(TEST)
      # Build the lapack tests
      set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/hybrid/debug/lapack")
      if(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
        set(LAPACK_TESTS ${LAPACK_TESTS_WITH_PMRRR})
      else(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
        set(LAPACK_TESTS ${LAPACK_TESTS_WITHOUT_PMRRR})
      endif(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
      foreach(TEST ${LAPACK_TESTS})
        add_executable(hybrid-debug-lapack-${TEST} 
                       ${TEST_DIR}/lapack/${TEST}.cpp)
        target_link_libraries(hybrid-debug-lapack-${TEST} hybrid-debug)
        set_target_properties(hybrid-debug-lapack-${TEST} 
                              PROPERTIES OUTPUT_NAME ${TEST}
                              RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                              LINK_FLAGS "${HYBRID_BLAS_LAPACK_FLAGS}")
      endforeach(TEST)
    endif(BUILD_TESTS)
  endif(BUILD_DEBUG)

  # Build the hybrid-release version (if necessary)
  if(BUILD_RELEASE)
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/lib/hybrid/release")
    if(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
      set(CMAKE_CXX_FLAGS 
          "${CXX_BASE_FLAGS} ${HYBRID_BLAS_DEFS} ${HYBRID_LAPACK_DEFS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_HYBRID_RELEASE_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRELEASE")
      add_library(hybrid-release STATIC ${ELEMENTAL_SRC})
      target_link_libraries(hybrid-release pmrrr-hybrid)
    else(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
      set(CMAKE_CXX_FLAGS 
          "${CXX_BASE_FLAGS} ${HYBRID_BLAS_DEFS} ${HYBRID_LAPACK_DEFS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_HYBRID_RELEASE_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRELEASE -DWITHOUT_PMRRR")
      add_library(hybrid-release STATIC ${ELEMENTAL_SRC})
    endif(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
    set_target_properties(hybrid-release
                          PROPERTIES OUTPUT_NAME elemental
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(BUILD_TESTS)
      # Build the core tests
      set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/hybrid/release/core")
      foreach(TEST ${CORE_TESTS})
        add_executable(hybrid-release-core-${TEST} ${TEST_DIR}/core/${TEST}.cpp)
        target_link_libraries(hybrid-release-core-${TEST} hybrid-release)
        set_target_properties(hybrid-release-core-${TEST} 
                              PROPERTIES OUTPUT_NAME ${TEST}
                              RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                              LINK_FLAGS "${HYBRID_BLAS_LAPACK_FLAGS}")
      endforeach(TEST)
      # Build the blas tests
      set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/hybrid/release/blas")
      foreach(TEST ${BLAS_TESTS})
        add_executable(hybrid-release-blas-${TEST} ${TEST_DIR}/blas/${TEST}.cpp) 
        target_link_libraries(hybrid-release-blas-${TEST} hybrid-release)
        set_target_properties(hybrid-release-blas-${TEST} 
                              PROPERTIES OUTPUT_NAME ${TEST}
                              RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                              LINK_FLAGS "${HYBRID_BLAS_LAPACK_FLAGS}")
      endforeach(TEST)
      # Build the lapack tests
      set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/hybrid/release/lapack")
      if(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
        set(LAPACK_TESTS ${LAPACK_TESTS_WITH_PMRRR})
      else(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
        set(LAPACK_TESTS ${LAPACK_TESTS_WITHOUT_PMRRR})
      endif(BUILD_PMRRR AND NOT FAILED_HYBRID_PMRRR)
      foreach(TEST ${LAPACK_TESTS})
        add_executable(hybrid-release-lapack-${TEST} 
                       ${TEST_DIR}/lapack/${TEST}.cpp)
        target_link_libraries(hybrid-release-lapack-${TEST} hybrid-release)
        set_target_properties(hybrid-release-lapack-${TEST} 
                              PROPERTIES OUTPUT_NAME ${TEST}
                              RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                              LINK_FLAGS "${HYBRID_BLAS_LAPACK_FLAGS}")
      endforeach(TEST)
    endif(BUILD_TESTS)
  endif(BUILD_RELEASE)
endif(BUILD_HYBRID AND HAVE_OPENMP)

