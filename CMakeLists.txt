#
#  Copyright (c) 2009-2012, Jack Poulson
#                     2012, Jed Brown (fixing installation issues)
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
cmake_minimum_required(VERSION 2.8.5)

project(Elemental)
set(Elemental_VERSION_MAJOR 0)
set(Elemental_VERSION_MINOR 75)

################################################################################
# Common build options                                                         #
################################################################################

# By default, Elemental builds static libraries, but shared libraries are 
# also supported through this option.
option(SHARED_LIBRARIES "Switch to shared libraries" OFF)

# Whether or not to build a collection of simple example drivers.
option(BUILD_EXAMPLES "A collection of simple examples" ON)

# Whether or not to build a collection of performance and correctness tests
option(BUILD_TESTS "A collection of performance and correctness tests" OFF)

# Whether or not to build current 'experimental' code, such as:
# 1) prototype C and F90 interfaces
# 2) Martin's 3D Gemm code
option(BUILD_EXPERIMENTAL "Build experimental code" OFF)

# Petschow and Bientinesi's Parallel Multiple Relatively Robust 
# Representations (PMRRR) is used at the heart of all of Elemental's current 
# eigensolvers. Since it makes use of pthreads, Windows builds cannot use 
# PMRRR.
if(WIN32)
  option(BUILD_PMRRR "Build and link against eigensolver PMRRR" OFF)
else(WIN32)
  option(BUILD_PMRRR "Build and link against eigensolver PMRRR" ON)
endif(WIN32)

################################################################################
# Advanced build options                                                       #
################################################################################

# Since MPI libraries frequently mishandle complex data, Elemental defaults
# to using routines designed for real datatypes whenever it does not result
# in a significant performance loss.
option(AVOID_COMPLEX_MPI "Avoid potentially buggy complex MPI routines" ON)
mark_as_advanced(AVOID_COMPLEX_MPI)

# At one point, a bug was found in IBM's C++ compiler for Blue Gene/P, 
# where OpenMP statements of the form a[i] += alpha b[i], with complex data,
# would segfault and/or return incorrect results
option(AVOID_OMP_FMA "Avoid a bug in the IBM compilers." ON)
mark_as_advanced(AVOID_OMP_FMA)

# Due to a subtle flaw in the Blue Gene/P extensions for MPICH2, treating 
# floating-point data as a collection of byte-sized objects results in a 
# better algorithm being chosen for MPI_Allgather. This should not effect
# performance on most machines.
option(USE_BYTE_ALLGATHERS "Avoid BG/P allgather performance bug." ON)
mark_as_advanced(USE_BYTE_ALLGATHERS)

# If MPI_Reduce_scatter_block doesn't exist, perform it by composing 
# MPI_Allreduce and std::memcpy rather than MPI_Reduce and MPI_Scatter
option(REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
       "AllReduce based block MPI_Reduce_scatter" OFF)
mark_as_advanced(REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE)

# This can easily be performed, but is likely not a good idea
option(POOL_MEMORY "Make Memory class accumulate memory until destruction" OFF)
mark_as_advanced(POOL_MEMORY)

################################################################################
# Elemental-development build options                                          #
################################################################################

# Print a warning any time a redistribution is performed which unpacks a 
# large amount of data with a non-unit stride
option(CACHE_WARNINGS "Warns when using cache-unfriendly routines" OFF)
mark_as_advanced(CACHE_WARNINGS)

# Print a warning when an improperly aligned redistribution is performed, 
# i.e., if an unnecessary permutation communication stage must take place
option(UNALIGNED_WARNINGS "Warn when performing unaligned redistributions" OFF)
mark_as_advanced(UNALIGNED_WARNINGS)

# Print a warning if an opportunity was missed to implement a redistribution
# approach specifically for vectors (instead of matrices)
option(VECTOR_WARNINGS "Warn when vector redistribution chances are missed" OFF)
mark_as_advanced(CACHE_WARNINGS UNALIGNED_WARNINGS VECTOR_WARNINGS)

################################################################################
# Significant command-line variable definitions                                #
################################################################################

# "CMAKE_BUILD_TYPE"
#   Elemental requires it to be one of the following choices:
#     1. "PureDebug": Vanilla MPI build meant for development and debugging
#     2. "PureRelease": Vanilla MPI build meant for production runs
#     3. "HybridDebug": MPI+OpenMP build meant for development and debugging
#     4. "HybridRelease": MPI+OpenMP build meant for production runs
#   If "CMAKE_BUILD_TYPE" is not defined, the default is "PureRelease".

# "CXX_PURE_DEBUG_FLAGS"/"CXX_PURE_RELEASE_FLAGS"/
# "CXX_HYBRID_DEBUG_FLAGS"/"CXX_HYBRID_RELEASE_FLAGS"
#   Override the default optimization (and debugging/symbol) flags

# "OpenMP_CXX_FLAGS"
#   Overrides the default compile flags for adding OpenMP support to CXX code

# TODO: Add list of MPI variables here, such as "MPI_CXX_COMPILE_FLAGS"

# TODO: Add list of MKL-specific optional variables, "INTEL_ROOT" and "MKL_ROOT"

# TODO: Add discussion of "MATH_LIBS"

# TODO: Add discussion of "REFERENCE_ROOT"

# TODO: Check to see if there are any others...

################################################################################
# The build logic starts here                                                  #
################################################################################

# Elemental must be built "out-of-source", so we start by ensuring that the
# source and build directories are different.
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source build attempted; please clean the CMake cache and then switch to an out-of-source build, e.g., rm CMakeCache.txt && rm -Rf CMakeFiles/ && mkdir build/ && cd build/ && cmake ..")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# Set the library type from the 'SHARED_LIBRARIES' option
if(SHARED_LIBRARIES)
  set(LIBRARY_TYPE SHARED)
else(SHARED_LIBRARIES)
  set(LIBRARY_TYPE STATIC)
endif(SHARED_LIBRARIES)

# Ensure that we have a valid build type
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Build mode not specified, defaulting to PureRelease build.")
  set(CMAKE_BUILD_TYPE PureRelease)
else(NOT CMAKE_BUILD_TYPE)
  if(NOT (CMAKE_BUILD_TYPE MATCHES HybridDebug) AND
     NOT (CMAKE_BUILD_TYPE MATCHES HybridRelease) AND 
     NOT (CMAKE_BUILD_TYPE MATCHES PureDebug) AND
     NOT (CMAKE_BUILD_TYPE MATCHES PureRelease))
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be [Hybrid,Pure][Debug,Release]")
  endif(NOT (CMAKE_BUILD_TYPE MATCHES HybridDebug) AND
        NOT (CMAKE_BUILD_TYPE MATCHES HybridRelease) AND 
        NOT (CMAKE_BUILD_TYPE MATCHES PureDebug) AND
        NOT (CMAKE_BUILD_TYPE MATCHES PureRelease))
endif(NOT CMAKE_BUILD_TYPE)

# Extract a few booleans from the build type
if(CMAKE_BUILD_TYPE STREQUAL "HybridDebug")
  set(HYBRID TRUE)
  set(DEBUG TRUE)
  set(HYBRID_DEBUG TRUE)
endif(CMAKE_BUILD_TYPE STREQUAL "HybridDebug")
if(CMAKE_BUILD_TYPE STREQUAL "HybridRelease")
  set(HYBRID TRUE)
  set(RELEASE TRUE)
  set(HYBRID_RELEASE TRUE)
endif(CMAKE_BUILD_TYPE STREQUAL "HybridRelease")
if(CMAKE_BUILD_TYPE STREQUAL "PureDebug")
  set(PURE TRUE)
  set(DEBUG TRUE)
  set(PURE_DEBUG TRUE)
endif(CMAKE_BUILD_TYPE STREQUAL "PureDebug")
if(CMAKE_BUILD_TYPE STREQUAL "PureRelease")
  set(PURE TRUE)
  set(RELEASE TRUE)
  set(PURE_RELEASE TRUE)   
endif(CMAKE_BUILD_TYPE STREQUAL "PureRelease")

# Set the basic compile flags from the build type
if(NOT WIN32)
  if(HYBRID_DEBUG AND NOT CXX_HYBRID_DEBUG_FLAGS)
    set(CXX_HYBRID_DEBUG_FLAGS "-O2 -g" CACHE STRING
        "CXX optimization/debug flags for hybrid OpenMP/MPI debug build")
  endif(HYBRID_DEBUG AND NOT CXX_HYBRID_DEBUG_FLAGS)
  if(HYBRID_RELEASE AND NOT CXX_HYBRID_RELEASE_FLAGS)
    set(CXX_HYBRID_RELEASE_FLAGS "-O3" CACHE STRING
        "CXX optimization flags for hybrid OpenMP/MPI release build")
  endif(HYBRID_RELEASE AND NOT CXX_HYBRID_RELEASE_FLAGS)
  if(PURE_DEBUG AND NOT CXX_PURE_DEBUG_FLAGS)
    set(CXX_PURE_DEBUG_FLAGS "-O2 -g" CACHE STRING
        "CXX optimization/debugging flags for pure MPI debug build")
  endif(PURE_DEBUG AND NOT CXX_PURE_DEBUG_FLAGS)
  if(PURE_RELEASE AND NOT CXX_PURE_RELEASE_FLAGS)
    set(CXX_PURE_RELEASE_FLAGS "-O3" CACHE STRING
        "CXX optimization flags for pure MPI release build")
  endif(PURE_RELEASE AND NOT CXX_PURE_RELEASE_FLAGS)
endif(NOT WIN32)

# Go ahead and check for Fortran, but keep in mind that CMake's 'OPTIONAL' 
# argument for enable_language is still completely broken as of 2.8.8
include(./cmake/language_support_v2.cmake)
workaround_9220(Fortran FORTRAN_WORKS)
if(FORTRAN_WORKS)
  enable_language(Fortran)  
  set(HAVE_F90_INTERFACE FALSE)
  if(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    include(FortranCInterface)
    FortranCInterface_VERIFY(CXX)
    if(FortranCInterface_VERIFIED_CXX)
      set(HAVE_F90_INTERFACE TRUE)
      FortranCInterface_HEADER(
        ${CMAKE_CURRENT_BINARY_DIR}/include/elemental/FCMangle.h 
        MACRO_NAMESPACE "FC_")
      install(FILES ${PROJECT_BINARY_DIR}/include/elemental/FCMangle.h
              DESTINATION include/elemental/)
    endif(FortranCInterface_VERIFIED_CXX)
  endif(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
endif(FORTRAN_WORKS)

# Make sure that MPI is available
find_package(MPI)
if(NOT MPI_CXX_FOUND)
  message(FATAL_ERROR "MPI C++ compiler was not found and is required")
endif(NOT MPI_CXX_FOUND)
include_directories(${MPI_CXX_INCLUDE_PATH})

# Search for OpenMP if necessary
if(HYBRID)
  if(OpenMP_CXX_FLAGS)
    set(HAVE_OPENMP TRUE)
    set(OpenMP_CXX_FLAGS ${OpenMP_CXX_FLAGS})
    message(STATUS "Using user-defined OpenMP_CXX_FLAGS=${OpenMP_CXX_FLAGS}")
  else(OpenMP_CXX_FLAGS)
    find_package(OpenMP)
    if(OPENMP_FOUND)
      set(HAVE_OPENMP TRUE)
    else(OPENMP_FOUND)
      set(OpenMP_CXX_FLAGS "" CACHE STRING "OpenMP CXX FLAGS")
      message(
        FATAL_ERROR
        "Hybrid build failed because OpenMP support was not detected. Please try specifying OpenMP_CXX_FLAGS."
      )
    endif(OPENMP_FOUND)
  endif(OpenMP_CXX_FLAGS)
endif(HYBRID)

# Append MPI and OpenMP compile flags as necessary
set(CMAKE_CXX_FLAGS_HYBRIDDEBUG 
    "${CXX_HYBRID_DEBUG_FLAGS} ${OpenMP_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_HYBRIDRELEASE 
    "${CXX_HYBRID_RELEASE_FLAGS} ${OpenMP_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_PUREDEBUG 
    "${CXX_PURE_DEBUG_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_PURERELEASE 
    "${CXX_PURE_RELEASE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")

# Query the size of a void pointer in order to determine whether or not this is
# a 32-bit system.
# Unfortunately, CMAKE_SIZEOF_VOID_P seems to be incorrect on some Windows 
# platforms, so I am switching to directly checking the size of void*
include(CheckTypeSize)
if(NOT SIZEOF_VOID_P)
  check_type_size("void*" SIZEOF_VOID_P)
endif(NOT SIZEOF_VOID_P)
set(32_BIT_SYSTEM FALSE)
if(${SIZEOF_VOID_P} EQUAL 4)
  set(32_BIT_SYSTEM TRUE)
endif(${SIZEOF_VOID_P} EQUAL 4)

# Since the CMake find_blas for some reason requires a Fortran compiler, we
# must build our own.
if(PURE)
  set(MATH_DESC "Unthreaded BLAS/LAPACK link flags")
else(PURE)
  set(MATH_DESC "Threaded BLAS/LAPACK link flags")
endif(PURE)
if(MATH_LIBS)
  set(MISSING_MATH FALSE)
  message(STATUS "Using user-defined MATH_LIBS=${MATH_LIBS}")
else(MATH_LIBS)
  if(APPLE)
    set(MATH_LIBS "-framework vecLib" CACHE STRING ${MATH_DESC})
    message(STATUS "Using Apple vecLib framework.")
    set(MISSING_MATH FALSE)
  else(APPLE)
    set(MISSING_MATH TRUE)
  endif(APPLE)
endif(MATH_LIBS)
if(MISSING_MATH)
  # Look for MKL first (no promises that this will work...)
  if(32_BIT_SYSTEM)
    set(INTEL_ARCH_SUBDIRS 32/lib lib/32)
    set(MKL_ARCH_SUBDIRS 32/lib lib/32)
  else(32_BIT_SYSTEM)
    set(INTEL_ARCH_SUBDIRS em64t/lib intel64/lib lib/em64t lib/intel64)
    set(MKL_ARCH_SUBDIRS em64t/lib intel64/lib lib/em64t lib/intel64 
        mkl/em64t/lib mkl/intel64/lib mkl/lib/em64t mkl/lib/intel64)
  endif(32_BIT_SYSTEM)
  if(INTEL_ROOT)
    message(STATUS "Will also search INTEL_ROOT=${INTEL_ROOT} for MKL")
    message(STATUS "INTEL_ARCH_SUBDIRS=${INTEL_ARCH_SUBDIRS}")
  endif(INTEL_ROOT)
  if(MKL_ROOT)
    message(STATUS "Will also search MKL_ROOT=${MKL_ROOT} for MKL")
    message(STATUS "MKL_ARCH_SUBDIRS=${INTEL_ARCH_SUBDIRS}")
  endif(MKL_ROOT)
  # There is a cyclic dependency between MKL_CORE and 
  # MKL_SEQUENTIAL/MKL_INTEL_THREAD with some versions of MKL
  if(NOT WIN32)
    find_library(PTHREAD_LIB pthread DOC "Pthreads library")
  endif(NOT WIN32)
  find_library(
    CMATH_LIB
    NAMES m libmmd.lib
    PATHS ${INTEL_ROOT} 
    PATH_SUFFIXES ${INTEL_ARCH_SUBDIRS}
    DOC "C math library"
  )
  find_library(
    IOMP5_LIB
    NAMES libiomp5.a libiomp5md.lib
    PATHS ${INTEL_ROOT} 
    PATH_SUFFIXES ${INTEL_ARCH_SUBDIRS}
  )
  find_library(
    GUIDE_LIB
    NAMES libguide.a libguide.lib
    PATHS ${INTEL_ROOT}
    PATH_SUFFIXES ${INTEL_ARCH_SUBDIRS}
  )
  find_library(
    MKL_CORE_LIB
    NAMES libmkl_core.a mkl_core.lib 
    PATHS ${MKL_ROOT}
    PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
  )
  find_library(
    MKL_INTEL_LIB
    NAMES libmkl_intel_lp64.a mkl_intel_lp64.lib 
    PATHS ${MKL_ROOT}
    PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
  )
  if(PURE)
    find_library(
      MKL_SEQUENTIAL_LIB
      NAMES libmkl_sequential.a mkl_sequential.lib
      PATHS ${MKL_ROOT}
      PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
    )
    set(MKL_REQUIRED MKL_INTEL MKL_SEQUENTIAL MKL_CORE MKL_SEQUENTIAL CMATH
                     GUIDE IOMP5)
  else(PURE)
    find_library(
      MKL_INTEL_THREAD_LIB
      NAMES libmkl_intel_thread.a mkl_intel_thread.lib
      PATHS ${MKL_ROOT}
      PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
    )
    set(MKL_REQUIRED MKL_INTEL MKL_INTEL_THREAD MKL_CORE MKL_INTEL_THREAD 
                     CMATH GUIDE IOMP5)
  endif(PURE)
  if(NOT WIN32)
    list(APPEND MKL_REQUIRED PTHREAD)
  endif(NOT WIN32)
  set(MKL_FOUND TRUE)
  set(MATH_LIBS "")
  foreach(NAME ${MKL_REQUIRED})
    if(${NAME}_LIB)
      message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
      list(APPEND MATH_LIBS ${${NAME}_LIB})
    else(${NAME}_LIB)
      message(STATUS "Could not find ${NAME}_LIB") 
      if(NOT "${NAME}" STREQUAL "GUIDE")
        set(MATH_LIBS "")
        set(MKL_FOUND FALSE)
      endif(NOT "${NAME}" STREQUAL "GUIDE")
    endif(${NAME}_LIB)
  endforeach(NAME)
  if(MKL_FOUND)
    message(STATUS "Using MKL math libraries.")
    message(STATUS "MATH_LIBS=${MATH_LIBS}")
    set(MISSING_MATH FALSE)
  endif(MKL_FOUND)
  if(MISSING_MATH)
    # Look for default BLAS and LAPACK
    if(REFERENCE_ROOT)
      message(STATUS 
        "Will search REFERENCE_ROOT=${REFERENCE_ROOT} for math libs")
    endif(REFERENCE_ROOT)
    set(REFERENCE_REQUIRED LAPACK BLAS)
    find_library(BLAS_LIB NAMES blas PATHS ${REFERENCE_ROOT})
    find_library(LAPACK_LIB NAMES lapack reflapack PATHS ${REFERENCE_ROOT}) 
    set(REFERENCE_FOUND TRUE)
    foreach(NAME ${REFERENCE_REQUIRED})
      if(${NAME}_LIB)
        message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
        list(APPEND MATH_LIBS ${${NAME}_LIB})
      else(${NAME}_LIB)
        message(STATUS "Could not find ${NAME}_LIB")
        set(MATH_LIBS "")
        set(REFERENCE_FOUND FALSE)
      endif(${NAME}_LIB)
    endforeach(NAME)
    if(REFERENCE_FOUND)
      message(WARNING "Using reference BLAS/LAPACK.")
    else(REFERENCE_FOUND) 
      message(FATAL_ERROR 
        "Could not find BLAS/LAPACK libraries. Please specify MATH_LIBS")
    endif(REFERENCE_FOUND)
  endif(MISSING_MATH)
endif(MISSING_MATH)

# Attempt to detect the BLAS and LAPACK underscore conventions. 
# We currently only handle whether or not there is an underscore appended.
include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
check_function_exists(daxpy HAVE_DAXPY)
if(HAVE_DAXPY)
  set(BLAS_POST FALSE)
  set(BLAS_DEFS "")
else(HAVE_DAXPY)
  check_function_exists(daxpy_ HAVE_DAXPY_POST)
  if(HAVE_DAXPY_POST)
    set(BLAS_POST TRUE)
    set(BLAS_DEFS "-DBLAS_POST")
  else(HAVE_DAXPY_POST)
    message(FATAL_ERROR "Could not determine BLAS format.")
  endif(HAVE_DAXPY_POST)
endif(HAVE_DAXPY)
check_function_exists(dpotrf HAVE_DPOTRF)
if(HAVE_DPOTRF)
  set(LAPACK_POST FALSE)
  set(LAPACK_DEFS "")
else(HAVE_DPOTRF)
  check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
  if(HAVE_DPOTRF_POST)
    set(LAPACK_POST TRUE)
    set(LAPACK_DEFS "-DLAPACK_POST")
  else(HAVE_DPOTRF_POST)
    message(FATAL_ERROR "Could not determine LAPACK format.")
  endif(HAVE_DPOTRF_POST)
endif(HAVE_DPOTRF)

# Attempt to detect whether or not we have MKL
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
check_function_exists(mkl_set_num_threads MATH_HAS_MKL)

# Check whether or not the FLAME bidiagonal QR routines are available
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
check_function_exists(FLA_Bsvd_v_opd_var1 HAVE_FLA_BSVD)

# Look for MPI_Reduce_scatter_block (and MPI_Reduce_scatter as sanity check)
set(CMAKE_REQUIRED_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_C_LINK_FLAGS}")
set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
check_function_exists(MPI_Reduce_scatter HAVE_MPI_REDUCE_SCATTER)
check_function_exists(MPI_Reduce_scatter_block HAVE_MPI_REDUCE_SCATTER_BLOCK)
if(NOT HAVE_MPI_REDUCE_SCATTER)
  message(FATAL_ERROR "Could not find MPI_Reduce_scatter")
endif(NOT HAVE_MPI_REDUCE_SCATTER)

# Check for nonblocking collectives
check_function_exists(MPI_Iallgather HAVE_MPI3_NONBLOCKING_COLLECTIVES)
check_function_exists(MPIX_Iallgather HAVE_MPIX_NONBLOCKING_COLLECTIVES)

# Look for restrict support
include(CheckCXXSourceCompiles)
set(RESTRICT_CODE
    "int main(void)
     {
         int* RESTRICT a;
         return 0;
     }
    ")
set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict__")
check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict__)
if(HAVE___restrict__)
  set(RESTRICT "__restrict__")
  message(STATUS "Using __restrict__ keyword.")
else(HAVE___restrict__)
  set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict")
  check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict)
  if(HAVE___restrict)
    set(RESTRICT "__restrict")
    message(STATUS "Using __restrict keyword.")
  else(HAVE___restrict)
    set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=restrict")
    check_cxx_source_compiles("${RESTRICT_CODE}" HAVE_restrict)
    if(HAVE_restrict)
      set(RESTRICT "restrict")
      message(STATUS "Using restrict keyword.")
    else(HAVE_restrict)
      set(RESTRICT "")
      message(STATUS "Could not find a restrict keyword.")
    endif(HAVE_restrict)
  endif(HAVE___restrict)
endif(HAVE___restrict__)

# Add the Parallel Linear Congruential Generator (PLCG) project
add_subdirectory(external/plcg)

# Add the Parallel Multiple Relatively Robust Representations (PMRRR) project
# if necessary
if(BUILD_PMRRR)
  add_subdirectory(external/pmrrr)
endif(BUILD_PMRRR)
if(BUILD_PMRRR AND NOT FAILED_PMRRR)
  set(WITHOUT_PMRRR FALSE)
else(BUILD_PMRRR AND NOT FAILED_PMRRR)
  set(WITHOUT_PMRRR TRUE)
endif(BUILD_PMRRR AND NOT FAILED_PMRRR)

# Create the Elemental configuration header
configure_file( 
  ${PROJECT_SOURCE_DIR}/include/elemental/config.h.cmake
  ${PROJECT_BINARY_DIR}/include/elemental/config.h 
) 
install(FILES ${PROJECT_BINARY_DIR}/include/elemental/config.h 
        DESTINATION include/elemental)

# Grab all of the .c, .cpp, .h, and .hpp Elemental files
file(GLOB_RECURSE ELEMENTAL_CPP RELATIVE ${PROJECT_SOURCE_DIR} 
     "src/*.c" "src/*.cpp")
file(GLOB_RECURSE ELEMENTAL_HEADERS RELATIVE ${PROJECT_SOURCE_DIR} 
     "include/*.h" "include/*.hpp")
set(ELEMENTAL_SRC "${ELEMENTAL_CPP};${ELEMENTAL_HEADERS}")

# Create a dummy library in order to be able to force the math (and potentially
# MPI libraries) to be linked last
add_library(cmake-dummy-lib ${LIBRARY_TYPE} cmake/CMakeDummyFunction.cpp)
target_link_libraries(cmake-dummy-lib ${MATH_LIBS} ${MPI_CXX_LIBRARIES})
install(TARGETS cmake-dummy-lib DESTINATION lib)

# The main library
if(BUILD_PMRRR AND NOT FAILED_PMRRR)
  add_library(elemental ${LIBRARY_TYPE} ${ELEMENTAL_SRC})
  target_link_libraries(elemental pmrrr plcg cmake-dummy-lib)
  if(MPI_CXX_LINK_FLAGS AND MPI_C_LINK_FLAGS)
    set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS} ${MPI_C_LINK_FLAGS}")
  elseif(MPI_CXX_LINK_FLAGS)
    set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  elseif(MPI_C_LINK_FLAGS)
    set(MPI_LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif(MPI_CXX_LINK_FLAGS AND MPI_C_LINK_FLAGS)
else(BUILD_PMRRR AND NOT FAILED_PMRRR)
  add_library(elemental ${LIBRARY_TYPE} ${ELEMENTAL_SRC})
  target_link_libraries(elemental plcg cmake-dummy-lib)
  set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif(BUILD_PMRRR AND NOT FAILED_PMRRR)
install(TARGETS elemental DESTINATION lib)

# Define the header-file preparation rules
set(PREPARED_HEADERS)
foreach(HEADER ${ELEMENTAL_HEADERS})
  add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/${HEADER}
            ${PROJECT_BINARY_DIR}/${HEADER}
    DEPENDS "${PROJECT_SOURCE_DIR}/${HEADER}"
  )
  list(APPEND PREPARED_HEADERS ${PROJECT_BINARY_DIR}/${HEADER})

  get_filename_component(HEADER_PATH ${HEADER} PATH)
  install(FILES ${PROJECT_BINARY_DIR}/${HEADER} DESTINATION ${HEADER_PATH})
endforeach(HEADER)
add_custom_target(prepare_elemental_headers DEPENDS ${PREPARED_HEADERS})
add_dependencies(elemental prepare_elemental_headers)

# Make sure the Elemental headers can be found
include_directories("${PROJECT_BINARY_DIR}/include")

# Build the test drivers if necessary
if(BUILD_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/tests)

  set(CORE_TESTS AxpyInterface Complex DifferentGrids DistMatrix Matrix)
  set(BLAS_TESTS Gemm Hemm Her2k Herk Symm Symv Syr2k Syrk Trmm Trsm)
  set(LAPACK_TESTS ApplyPackedReflectors Cholesky Hegst HermitianTridiag LDL 
                   LU LQ QR TriangularInverse)

  if(BUILD_PMRRR AND NOT FAILED_PMRRR)
    list(APPEND LAPACK_TESTS HermitianGenDefiniteEig)
    list(APPEND LAPACK_TESTS HermitianEig)
  endif(BUILD_PMRRR AND NOT FAILED_PMRRR)

  # Build the core tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/core")
  foreach(TEST ${CORE_TESTS})
    add_executable(tests-core-${TEST} ${TEST_DIR}/core/${TEST}.cpp)
    target_link_libraries(tests-core-${TEST} elemental)
    set_target_properties(tests-core-${TEST} PROPERTIES
                          OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(tests-core-${TEST} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS tests-core-${TEST} DESTINATION bin/tests/core)
  endforeach(TEST)
  # Build the BLAS-like tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/blas-like")
  foreach(TEST ${BLAS_TESTS})
    add_executable(tests-blas-${TEST} ${TEST_DIR}/blas-like/${TEST}.cpp)
    target_link_libraries(tests-blas-${TEST} elemental)
    set_target_properties(tests-blas-${TEST} PROPERTIES
                          OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(tests-blas-${TEST} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS tests-blas-${TEST} DESTINATION bin/tests/blas-like)
  endforeach(TEST)
  # Build the LAPACK-like tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/lapack-like")
  foreach(TEST ${LAPACK_TESTS})
    add_executable(tests-lapack-${TEST} ${TEST_DIR}/lapack-like/${TEST}.cpp)
    target_link_libraries(tests-lapack-${TEST} elemental)
    set_target_properties(tests-lapack-${TEST} PROPERTIES
                          OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(tests-lapack-${TEST} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS tests-lapack-${TEST} DESTINATION bin/tests/lapack-like)
  endforeach(TEST)
endif(BUILD_TESTS)

# Build the example drivers if necessary
if(BUILD_EXAMPLES)
  set(EXAMPLE_DIR ${PROJECT_SOURCE_DIR}/examples)

  set(CORE_EXAMPLES Constructors)
  set(BLAS_EXAMPLES Gemv)
  set(LAPACK_EXAMPLES GaussianElimination HermitianQDWH HouseholderSolve 
                      HPDInverse LDL LDLInverse Polar Pseudoinverse QDWH 
                      SequentialQR QR SequentialSVD SimpleSVD SVD)
  set(SPECIAL_EXAMPLES Cauchy CauchyLike Circulant Diagonal DiscreteFourier 
                       Hankel HermitianUniformSpectrum Hilbert Identity 
                       NormalUniformSpectrum Ones OneTwoOne Toeplitz Uniform
                       Walsh Wilkinson Zeros)

  if(BUILD_PMRRR AND NOT FAILED_PMRRR)
    list(APPEND LAPACK_EXAMPLES HermitianEig 
                                HermitianEigFromSequential
                                HermitianSVD
                                ComplexHermitianFunction 
                                RealHermitianFunction
                                HermitianPseudoinverse
                                HPSDCholesky
                                HPSDSquareRoot)
    list(APPEND LAPACK_EXAMPLES SymmetricEig RealSymmetricFunction)
  endif(BUILD_PMRRR AND NOT FAILED_PMRRR)

  # Build the core tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples/core")
  foreach(EXAMPLE ${CORE_EXAMPLES})
    add_executable(examples-core-${EXAMPLE} ${EXAMPLE_DIR}/core/${EXAMPLE}.cpp)
    target_link_libraries(examples-core-${EXAMPLE} elemental)
    set_target_properties(examples-core-${EXAMPLE} PROPERTIES 
                          OUTPUT_NAME ${EXAMPLE} 
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(examples-core-${EXAMPLE} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS examples-core-${EXAMPLE} DESTINATION bin/examples/core)
  endforeach(EXAMPLE)

  # Build the BLAS-like tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples/blas-like")
  foreach(EXAMPLE ${BLAS_EXAMPLES})
    add_executable(examples-blas-${EXAMPLE} 
                   ${EXAMPLE_DIR}/blas-like/${EXAMPLE}.cpp)
    target_link_libraries(examples-blas-${EXAMPLE} elemental)
    set_target_properties(examples-blas-${EXAMPLE} PROPERTIES
                          OUTPUT_NAME ${EXAMPLE}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(examples-blas-${EXAMPLE} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS examples-blas-${EXAMPLE} 
            DESTINATION bin/examples/blas-like)
  endforeach(EXAMPLE)

  # Build the LAPACK-like tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples/lapack-like")
  foreach(EXAMPLE ${LAPACK_EXAMPLES})
    add_executable(examples-lapack-${EXAMPLE} 
                   ${EXAMPLE_DIR}/lapack-like/${EXAMPLE}.cpp)
    target_link_libraries(examples-lapack-${EXAMPLE} elemental)
    set_target_properties(examples-lapack-${EXAMPLE} PROPERTIES
                          OUTPUT_NAME ${EXAMPLE}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(examples-lapack-${EXAMPLE} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS examples-lapack-${EXAMPLE} 
            DESTINATION bin/examples/lapack-like)
  endforeach(EXAMPLE)

  # Build the special matrix test(s)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples/special_matrices")
  foreach(EXAMPLE ${SPECIAL_EXAMPLES})
    add_executable(examples-special-${EXAMPLE}
                   ${EXAMPLE_DIR}/special_matrices/${EXAMPLE}.cpp)
    target_link_libraries(examples-special-${EXAMPLE} elemental)
    set_target_properties(examples-special-${EXAMPLE} PROPERTIES
                          OUTPUT_NAME ${EXAMPLE}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(examples-special-${EXAMPLE} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS examples-special-${EXAMPLE}
            DESTINATION bin/examples/special_matrices)
  endforeach(EXAMPLE)
endif(BUILD_EXAMPLES)

# Build experimental drivers that do not require an eigensolver
if(BUILD_EXPERIMENTAL)
  set(EXPERIMENTAL_DIR ${PROJECT_SOURCE_DIR}/experimental)

  set(G3D_EXPERS G3DGemm)

  # Build the G3D example(s)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/experimental/g3d")
  foreach(EXPER ${G3D_EXPERS})
    add_executable(experimental-g3d-${EXPER} 
                   ${EXPERIMENTAL_DIR}/g3d/${EXPER}.cpp)
    target_link_libraries(experimental-g3d-${EXPER} elemental)
    set_target_properties(experimental-g3d-${EXPER} PROPERTIES 
                          OUTPUT_NAME ${EXPER} 
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(experimental-g3d-${EXPER} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS experimental-g3d-${EXPER} DESTINATION bin/experimental/g3d)
  endforeach(EXPER)
endif(BUILD_EXPERIMENTAL)

# Build experimental drivers that DO require an eigensolver
if(BUILD_EXPERIMENTAL AND BUILD_PMRRR AND NOT FAILED_PMRRR)
  if(NOT MPI_C_FOUND)
    message(FATAL_ERROR "Could not build C interface without C MPI")
  endif(NOT MPI_C_FOUND)
  include_directories(${MPI_C_INCLUDE_PATH})
  set(CMAKE_C_FLAGS_HYBRIDDEBUG "${MPI_C_COMPILE_FLAGS}")
  set(CMAKE_C_FLAGS_HYBRIDRELEASE "${MPI_C_COMPILE_FLAGS}")
  set(CMAKE_C_FLAGS_PUREDEBUG "${MPI_C_COMPILE_FLAGS}")
  set(CMAKE_C_FLAGS_PURERELEASE "${MPI_C_COMPILE_FLAGS}")

  set(EXPERIMENTAL_DIR ${PROJECT_SOURCE_DIR}/experimental)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/experimental/c_interface")

  # Build the experimental C interface
  file(GLOB_RECURSE EXPERIMENTAL_C_CPP RELATIVE ${PROJECT_SOURCE_DIR} 
       "experimental/c_interface/c_interface.cpp")
  add_library(experimental-c-interface ${LIBRARY_TYPE} ${EXPERIMENTAL_C_CPP})

  add_executable(
      experimental-c-interface-test
      ${EXPERIMENTAL_DIR}/c_interface/test.c
  )
  target_link_libraries(
      experimental-c-interface-test experimental-c-interface elemental
  )
  set_target_properties(
      experimental-c-interface-test PROPERTIES OUTPUT_NAME c_interface_test
                                    RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
  )
  if(MPI_LINK_FLAGS)
    set_target_properties(
        experimental-c-interface-test PROPERTIES LINK_FLAGS ${MPI_LINK_FLAGS}
    )
  endif(MPI_LINK_FLAGS)
  install(TARGETS experimental-c-interface-test 
          DESTINATION bin/experimental/c_interface)  
  add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/include/elemental.h
    COMMAND ${CMAKE_COMMAND} -E copy 
            ${PROJECT_SOURCE_DIR}/experimental/c_interface/elemental.h
            ${PROJECT_BINARY_DIR}/include/elemental.h
    DEPENDS "${PROJECT_SOURCE_DIR}/experimental/c_interface/elemental.h"
  )
  install(FILES ${PROJECT_BINARY_DIR}/include/elemental.h 
          DESTINATION include/elemental.h)

  if(HAVE_F90_INTERFACE)
    if(NOT MPI_Fortran_FOUND)
      message(FATAL_ERROR "Could not build F90 interface without Fortran MPI")
    endif(NOT MPI_Fortran_FOUND)

    # We cannot simply use check_function_exists, since MPI_Comm_f2c is 
    # a macro in some MPI implementations
    set(MPI_COMM_F2C_CODE
    "#include \"mpi.h\"
     int main( int argc, char* argv[] )
     {
         MPI_Init( &argc, &argv );
         MPI_Fint fint;
         MPI_Comm comm = MPI_Comm_f2c( fint );
         MPI_Finalize();
         return 0;
     }
    ")
    include(CheckCSourceCompiles)
    check_c_source_compiles("${MPI_COMM_F2C_CODE}" HAVE_MPI_COMM_F2C)
    if(NOT HAVE_MPI_COMM_F2C)  
      message(WARNING 
        "Could not find MPI_Comm_f2c, but it is needed for the F90 interface")
    endif(NOT HAVE_MPI_COMM_F2C)

    include_directories(${MPI_Fortran_INCLUDE_PATH})
    set(CMAKE_Fortran_FLAGS_HYBRIDDEBUG "${MPI_Fortran_COMPILE_FLAGS}")
    set(CMAKE_Fortran_FLAGS_HYBRIDRELEASE "${MPI_Fortran_COMPILE_FLAGS}")
    set(CMAKE_Fortran_FLAGS_PUREDEBUG "${MPI_Fortran_COMPILE_FLAGS}")
    set(CMAKE_Fortran_FLAGS_PURERELEASE "${MPI_Fortran_COMPILE_FLAGS}")

    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/experimental/f90_interface")

    # Build the experimental F90 interface
    file(GLOB_RECURSE EXPERIMENTAL_F90_CPP RELATIVE ${PROJECT_SOURCE_DIR} 
         "experimental/f90_interface/f90_interface.cpp")
    add_library(experimental-f90-interface ${LIBRARY_TYPE} ${EXPERIMENTAL_F90_CPP})

    add_executable(
        experimental-f90-interface-test
        ${EXPERIMENTAL_DIR}/f90_interface/test.f90
    )
    target_link_libraries(
        experimental-f90-interface-test experimental-f90-interface elemental
    )
    set_target_properties(
        experimental-f90-interface-test 
        PROPERTIES OUTPUT_NAME f90_interface_test
        RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
    )
    if(MPI_LINK_FLAGS)
      set_target_properties(
          experimental-f90-interface-test 
          PROPERTIES LINK_FLAGS ${MPI_LINK_FLAGS}
      )
    endif(MPI_LINK_FLAGS)
    install(TARGETS experimental-f90-interface-test 
            DESTINATION bin/experimental/f90_interface)  
  endif(HAVE_F90_INTERFACE)
endif(BUILD_EXPERIMENTAL AND BUILD_PMRRR AND NOT FAILED_PMRRR)

################################################################################
# Uncomment if including Elemental as a subproject in another build system     #
################################################################################
#set(LIBRARY_TYPE ${LIBRARY_TYPE} PARENT_SCOPE)
#set(MPI_C_COMPILER ${MPI_C_COMPILER} PARENT_SCOPE)
#set(MPI_C_INCLUDE_PATH ${MPI_C_INCLUDE_PATH} PARENT_SCOPE)
#set(MPI_C_COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS} PARENT_SCOPE)
#set(MPI_C_LINK_FLAGS ${MPI_C_LINK_FLAGS} PARENT_SCOPE)
#set(MPI_C_LIBRARIES ${MPI_C_LIBRARIES} PARENT_SCOPE)
#set(MPI_CXX_COMPILER ${MPI_CXX_COMPILER} PARENT_SCOPE)
#set(MPI_CXX_INCLUDE_PATH ${MPI_CXX_INCLUDE_PATH} PARENT_SCOPE)
#set(MPI_CXX_COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS} PARENT_SCOPE)
#set(MPI_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS} PARENT_SCOPE)
#set(MPI_CXX_LIBRARIES ${MPI_CXX_LIBRARIES} PARENT_SCOPE)
#set(MATH_LIBS ${MATH_LIBS} PARENT_SCOPE)
#set(RESTRICT ${RESTRICT} PARENT_SCOPE)
#set(RELEASE ${RELEASE} PARENT_SCOPE)
#set(BLAS_POST ${BLAS_POST} PARENT_SCOPE)
#set(LAPACK_POST ${LAPACK_POST} PARENT_SCOPE)
#set(HAVE_F90_INTERFACE ${HAVE_F90_INTERFACE} PARENT_SCOPE)
#set(WITHOUT_PMRRR ${WITHOUT_PMRRR} PARENT_SCOPE)
#set(AVOID_COMPLEX_MPI ${AVOID_COMPLEX_MPI} PARENT_SCOPE)
#set(HAVE_REDUCE_SCATTER_BLOCK ${HAVE_REDUCE_SCATTER_BLOCK} PARENT_SCOPE)
#set(REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE ${REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE} PARENT_SCOPE)
#set(USE_BYTE_ALLGATHERS ${USE_BYTE_ALLGATHERS} PARENT_SCOPE)
