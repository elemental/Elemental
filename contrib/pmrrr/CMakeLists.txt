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

if(CMAKE_BUILD_TYPE)
  message(STATUS "PMRRR does not use CMAKE_BUILD_TYPE.")
endif(CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "" CACHE INTERNAL "Removing GUI option" FORCE)

option(HAVE_SPINLOCKS "Enable if pthread lib supports spinlocks" OFF)
MARK_AS_ADVANCED(HAVE_SPINLOCKS)

set(C_OPT_FLAGS "-O3" CACHE STRING "C optimization flags")

# Search the BLAS/LAPACK libs for  to test whether it is a new enough version
# for our purposes. If the LAPACK library is version 3.1 or later, then 
# dstemr will support subset computations. However, dlamch only became 
# threadsafe at version 3.3. If Elemental only requires a pure MPI build, 
# then we don't need to worry about threadsafe routines and can settle for 
# 3.1 or later.
#
# Based upon the netlib changesets, we can determine if the LAPACK version 
# is >= 3.1 by searching for claqr0, and likewise searching for cuncsd to test
# if the version is >= 3.3.
include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES ${PURE_BLAS_LAPACK_LIBS})
if(HAVE_DPOTRF)
  check_function_exists(claqr0 LAPACK_3_1)
  check_function_exists(cuncsd LAPACK_3_3)
elseif(HAVE_DPOTRF_PRE)
  check_function_exists(_claqr0 LAPACK_3_1)
  check_function_exists(_cuncsd LAPACK_3_3)
elseif(HAVE_DPOTRF_POST)
  check_function_exists(claqr0_ LAPACK_3_1)
  check_function_exists(cuncsd_ LAPACK_3_3)
elseif(HAVE_DPOTRF_PRE_POST)
  check_function_exists(_claqr0_ LAPACK_3_1)
  check_function_exists(_cuncsd_ LAPACK_3_3)
endif(HAVE_DPOTRF)

if(BUILD_PURE AND NOT LAPACK_3_1)
  set(MISSING_LAPACK_3_1 TRUE)
  message(STATUS "Pure Elemental+PMRRR requires LAPACK >= 3.1.")
else(BUILD_PURE AND NOT LAPACK_3_1)
  set(MISSING_LAPACK_3_1 FALSE)
endif(BUILD_PURE AND NOT LAPACK_3_1)

if(BUILD_HYBRID AND HAVE_OPENMP AND NOT LAPACK_3_3)
  set(MISSING_LAPACK_3_3 TRUE)
  message(STATUS "Hybrid Elemental+PMRRR requires LAPACK >= 3.3.")
elseif(BUILD_HYBRID AND HAVE_OPENMP AND NOT LAPACK_3_3)
  set(MISSING_LAPACK_3_3 FALSE)
endif(BUILD_HYBRID AND HAVE_OPENMP AND NOT LAPACK_3_3)

if(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
  # For some reason enable_language( ... OPTIONAL) causes an error, so we 
  # should instead manually probe for a Fortran compiler first
  #enable_language(Fortran OPTIONAL)
  find_program(FORTRAN_COMPILER gfortran g77 ifort xlf_r xlf xlf90_r xlf90
                                xlf95_r xlf95 xlf2003_r xlf2003)
  if(FORTRAN_COMPILER)
    enable_language(Fortran)

    # Now determine the way C should interface with our Fortran addons. 
    get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
    set(LAPACK_ADDONS_LINK_LIBS "")
    if(Fortran_COMPILER_NAME STREQUAL "gfortran")
      # gfortran should append an underscore by default
      set(F_FROM_C_DEFS "-DBLAS_POST -DLAPACK_POST")
      if(APPLE)
        find_library(GFORTRAN_LIB gfortran)
        if(GFORTRAN_LIB)
          set(LAPACK_ADDONS_LINK_LIBS ${GFORTRAN_LIB})
        else(GFORTRAN_LIB)
          message(FATAL_ERROR "Could not find gfortran library")
        endif(GFORTRAN_LIB)
      endif(APPLE)
    else(Fortran_COMPILER_NAME STREQUAL "gfortran")
      # The list of other compiler options should grow quickly. It is best
      # to default to assuming an appended underscore.
      set(F_FROM_C_DEFS "-DBLAS_POST -DLAPACK_POST")
    endif(Fortran_COMPILER_NAME STREQUAL "gfortran")
    set(F_FROM_C_DEFS "${F_FROM_C_DEFS} -DPREPEND_X")
    
    file(GLOB_RECURSE LAPACK_ADDON_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
         "*.f")
    add_library(lapack-addons STATIC ${LAPACK_ADDON_SRC})
    target_link_libraries(lapack-addons ${LAPACK_ADDONS_LINK_LIBS})
  else(FORTRAN_COMPILER)
    if(MISSING_LAPACK_3_1)
      set(FAILED_PURE_PMRRR TRUE)
      set(FAILED_HYBRID_PMRRR TRUE)
      set(FAILED_PURE_PMRRR TRUE PARENT_SCOPE)
      set(FAILED_HYBRID_PMRRR TRUE PARENT_SCOPE)
      message(STATUS "Fortran compiler not found; we cannot build PMRRR")
    else(MISSING_LAPACK_3_1)
      set(FAILED_HYBRID_PMRRR TRUE)
      set(FAILED_HYBRID_PMRRR TRUE PARENT_SCOPE)
      message(STATUS "Fortran compiler not found; we cannot build hybrid PMRRR")
      # Salvage the pure PMRRR build without a Fortran compiler
      set(MISSING_LAPACK_3_3 FALSE)
    endif(MISSING_LAPACK_3_1)
  endif(FORTRAN_COMPILER)
endif(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)

if(FAILED_PURE_PMRRR)
  # Do nothing, we are falling back to building Elemental without PMRRR
else(FAILED_PURE_PMRRR)
  # Add in the preprocessor definitions 
  if(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
    set(C_BASE_FLAGS "${F_FROM_C_DEFS}")
  else(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
    set(C_BASE_FLAGS "${PURE_BLAS_DEFS} ${PURE_LAPACK_DEFS}")
  endif(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
  if(NOT HAVE_SPINLOCKS)
    set(C_BASE_FLAGS "${C_BASE_FLAGS} -DNOSPINLOCKS")
  endif(NOT HAVE_SPINLOCKS)

  # Grab all of the .c and .h files
  file(GLOB_RECURSE PMRRR_C RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.c")
  file(GLOB_RECURSE PMRRR_H RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.h")
  set(PMRRR_SRC "${PMRRR_C};${PMRRR_H}")

  # Search for the MPI C compiler
  if(MPI_C_COMPILER AND NOT EXISTS ${MPI_C_COMPILER})
    message(SEND_ERROR
            "Invalid MPI_C_COMPILER: ${MPI_C_COMPILER} does not exist.")
  endif(MPI_C_COMPILER AND NOT EXISTS ${MPI_C_COMPILER})
  if(NOT MPI_C_COMPILER)
    find_program(MPI_C_COMPILER mpicc mpixlc mpixlc_r 
      DOC "Full path to MPI C compiler. If blank, CMake will search for it.")
    if(NOT MPI_C_COMPILER)
      message(SEND_ERROR 
        "Could not find MPI C compiler. Please manually set MPI_C_COMPILER.")
    endif(NOT MPI_C_COMPILER)
  endif(NOT MPI_C_COMPILER)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  
  # Findpthreads.cmake does not seem to be standard
  #
  #find_package(pthreads REQUIRED)
  #set(C_BASE_FLAGS "${C_BASE_FLAGS} -I${PTHREADS_INCLUDE_DIR}")
  #foreach(DEF ${PTHREADS_DEFINITIONS})
  #  set(C_BASE_FLAGS "${C_BASE_FLAGS} -D${DEF}")
  #endforeach(DEF)
  set(C_PTHREADS_FLAGS "-pthread" CACHE STRING 
      "C compile flag for pthreads support.")

  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
  set(CMAKE_C_FLAGS "${C_BASE_FLAGS} ${C_OPT_FLAGS} ${C_PTHREADS_FLAG}")
  set(CMAKE_Fortran_FLAGS "${FORTRAN_OPT_FLAGS}")
  if(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
    add_library(pmrrr-pure STATIC ${PMRRR_SRC})
    target_link_libraries(pmrrr-pure lapack-addons)
    if(BUILD_HYBRID AND HAVE_OPENMP AND NOT FAILED_HYBRID_PMRRR)
      add_library(pmrrr-hybrid STATIC ${PMRRR_SRC})
      target_link_libraries(pmrrr-hybrid lapack-addons)
    endif(BUILD_HYBRID AND HAVE_OPENMP AND NOT FAILED_HYBRID_PMRRR)
  else(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
    add_library(pmrrr-pure STATIC ${PMRRR_SRC})
    if(BUILD_HYBRID AND HAVE_OPENMP AND NOT FAILED_HYBRID_PMRRR)
      add_library(pmrrr-hybrid STATIC ${PMRRR_SRC})
    endif(BUILD_HYBRID AND HAVE_OPENMP AND NOT FAILED_HYBRID_PMRRR)
  endif(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
endif(FAILED_PURE_PMRRR)

