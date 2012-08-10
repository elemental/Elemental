#
#  Copyright (c) 2009-2012, Jack Poulson
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

option(HAVE_SPINLOCKS "Enable if pthread lib supports spinlocks" OFF)
MARK_AS_ADVANCED(HAVE_SPINLOCKS)

if(HYBRID_DEBUG)
  set(
    C_HYBRID_DEBUG_FLAGS "-O2 -g" CACHE STRING
    "C optimization/debug flags for hybrid OpenMP/MPI debug build"
  )
endif()
if(HYBRID_RELEASE)
  set(
    C_HYBRID_RELEASE_FLAGS "-O3" CACHE STRING
    "C optimization flags for hybrid OpenMP/MPI release build"
  )
endif()
if(PURE_DEBUG)
  set(
    C_PURE_DEBUG_FLAGS "-O2 -g" CACHE STRING
    "C optimization/debugging flags for pure MPI debug build"
  )
endif()
if(PURE_RELEASE)
  set(
    C_PURE_RELEASE_FLAGS "-O3" CACHE STRING
    "C optimization flags for pure MPI release build"
  )
endif()

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
#
# Unfortunately, these checks are not sufficient for MKL.
include(CheckFunctionExists)

set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
if(HAVE_DPOTRF)
  check_function_exists(claqr0 LAPACK_3_1)
  check_function_exists(cuncsd LAPACK_3_3)
elseif(HAVE_DPOTRF_POST)
  check_function_exists(claqr0_ LAPACK_3_1)
  check_function_exists(cuncsd_ LAPACK_3_3)
endif()

set(CMAKE_REQUIRED_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_C_LINK_FLAGS}")
set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
check_function_exists(MPI_Query_thread HAVE_MPI_QUERY_THREAD)
if(NOT HAVE_MPI_QUERY_THREAD)
  message(WARNING "Do not have MPI_Query_thread, so cannot build PMRRR")
  set(FAILED_PMRRR TRUE PARENT_SCOPE)
  return()
endif()

if(PURE)
  if(MATH_HAS_MKL)
    set(MISSING_LAPACK TRUE)
    message(STATUS "Pure Elemental+PMRRR cannot use LAPACK's MRRR.")
  elseif(NOT LAPACK_3_1)
    set(MISSING_LAPACK TRUE)
    message(STATUS "Pure Elemental+PMRRR requires LAPACK >= 3.1.")
  else()
    set(MISSING_LAPACK FALSE)
  endif()
else()
  # We need at least 3.3
  if(MATH_HAS_MKL)
    set(MISSING_LAPACK TRUE)
    message(STATUS "Hybrid PMRRR cannot use MKL's LAPACK routines.")
  elseif(NOT LAPACK_3_3)
    set(MISSING_LAPACK TRUE)
    message(STATUS "Hybrid Elemental+PMRRR requires LAPACK >= 3.3.")
  else()
    set(MISSING_LAPACK FALSE)
  endif()
endif()

if(MISSING_LAPACK)
  if(CMAKE_Fortran_COMPILER_WORKS)
    if(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
      message(
        WARNING 
        "${CMAKE_Fortran_COMPILER} does not support F90, so PMRRR cannot be built"
      )
      set(FAILED_PMRRR TRUE PARENT_SCOPE)
      return()
    endif()

    if(HYBRID_DEBUG)
      set(
        Fortran_HYBRID_DEBUG_FLAGS "-O2 -g" CACHE STRING
        "Fortran optimization/debug flags for hybrid threaded/MPI debug build"
      )
    endif()
    if(HYBRID_RELEASE)
      set(
        Fortran_HYBRID_RELEASE_FLAGS "-O3" CACHE STRING
        "Fortran optimization flags for hybrid threaded/MPI release build"
      )
    endif()
    if(PURE_DEBUG)
      set(
        Fortran_PURE_DEBUG_FLAGS "-O2 -g" CACHE STRING
        "Fortran optimization/debugging flags for pure MPI debug build"
      )
    endif()
    if(PURE_RELEASE)
      set(
        Fortran_PURE_RELEASE_FLAGS "-O3" CACHE STRING
        "Fortran optimization flags for pure MPI release build"
      )
    endif()

    set(CMAKE_Fortran_FLAGS_HYBRIDDEBUG ${Fortran_HYBRID_DEBUG_FLAGS})
    set(CMAKE_Fortran_FLAGS_HYBRIDRELEASE ${Fortran_HYBRID_RELEASE_FLAGS})
    set(CMAKE_Fortran_FLAGS_PUREDEBUG ${Fortran_PURE_DEBUG_FLAGS})
    set(CMAKE_Fortran_FLAGS_PURERELEASE ${Fortran_PURE_RELEASE_FLAGS})

    # Now determine the way C should interface with our Fortran addons. 
    get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
    set(LAPACK_ADDONS_LIBS "")
    if(Fortran_COMPILER_NAME STREQUAL "gfortran")
      # This library should be in your LD_LIBRARY_PATH. CMake does not appear
      # to find it in cases where it is a shared library with extensions, e.g.,
      # if /usr/lib/libgfortran.so.1 exists, but /usr/lib/libgfortran.so 
      # does not exist. 
      if(GFORTRAN_LIB)
        set(LAPACK_ADDONS_LIBS ${GFORTRAN_LIB})
        message(STATUS "Using user-defined GFORTRAN_LIB=${GFORTRAN_LIB}")
      else()
        string(
          REGEX REPLACE "bin/gfortran" "lib" GFORTRAN_LIB_DIR_HINT
          "${CMAKE_Fortran_COMPILER}"
        )
        message(
          STATUS 
          "Will also search ${GFORTRAN_LIB_DIR_HINT} for gfortran lib."
        )
        find_library(
          GFORTRAN_LIB 
          NAMES gfortran libgfortran.a libgfortran.so libgfortran.so.1 
                libgfortran.so.2 libgfortran.so.3
          HINTS ${GFORTRAN_LIB_DIR_HINT}
        )
        if(GFORTRAN_LIB)
          message(STATUS "Found gfortran lib at ${GFORTRAN_LIB}")
        else()
          message(
            FATAL_ERROR 
            "Could not find gfortran library. Please specify it by defining GFORTRAN_LIB."
          )
        endif()
      endif()
      set(LAPACK_ADDONS_LIBS ${GFORTRAN_LIB})
    elseif(Fortran_COMPILER_NAME STREQUAL "ifort")
      if(IFCORE_LIB)
        message(STATUS "Using user-defined IFCORE_LIB=${IFCORE_LIB}")
        set(LAPACK_ADDONS_LIBS ${IFCORE_LIB})
      else()
        string(
          REGEX REPLACE "/bin/([^/]*)/?ifort" "/lib/\\1" 
          IFCORE_LIB_DIR_HINT "${CMAKE_Fortran_COMPILER}"
        )
        message(STATUS "Will also search ${IFCORE_LIB_DIR_HINT} for ifcore lib")
        find_library(IFCORE_LIB ifcore HINTS ${IFCORE_LIB_DIR_HINT})
        if(IFCORE_LIB)
          set(LAPACK_ADDONS_LIBS ${IFCORE_LIB})
          message("Found ifcore lib at ${IFCORE_LIB}")
        else()
          message(FATAL_ERROR "Could not find ifcore library")
        endif()
      endif()
    endif()
   
    file(
      GLOB_RECURSE LAPACK_ADDON_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.f"
    )
    add_library(lapack-addons ${LIBRARY_TYPE} ${LAPACK_ADDON_SRC})
    target_link_libraries(lapack-addons ${LAPACK_ADDONS_LIBS})
  else()
    set(FAILED_PMRRR TRUE PARENT_SCOPE)
    message(WARNING "F90 compiler not found; cannot build PMRRR using the supplied MATH_LIBS")
    return()
  endif()
endif()
set(LAPACK_ADDONS_LIBS ${LAPACK_ADDONS_LIBS} PARENT_SCOPE)

# Add in the preprocessor definitions and fortran name-mangling
if(MISSING_LAPACK)
  set(C_BASE_FLAGS "-DCUSTOM_BLAS_LAPACK")
  include(FortranCInterface)
  FortranCInterface_VERIFY()
  FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/include/FCMangle.h 
                           MACRO_NAMESPACE "FC_")
else()
  set(C_BASE_FLAGS "${BLAS_DEFS} ${LAPACK_DEFS}")
endif()
if(NOT HAVE_SPINLOCKS)
  set(C_BASE_FLAGS "${C_BASE_FLAGS} -DNOSPINLOCKS")
endif()

# Ensure that an MPI C compiler was found
if(NOT MPI_C_FOUND)
  message(WARNING "No MPI C compiler was found, so PMRRR cannot be built")
  set(FAILED_PMRRR TRUE PARENT_SCOPE)
  return()
endif()
include_directories(${MPI_C_INCLUDE_PATH})
    
# Set the various build flags
set(CMAKE_C_FLAGS_HYBRIDDEBUG "${C_BASE_FLAGS} ${C_HYBRID_DEBUG_FLAGS}")
set(CMAKE_C_FLAGS_HYBRIDRELEASE "${C_BASE_FLAGS} ${C_HYBRID_RELEASE_FLAGS}")
set(CMAKE_C_FLAGS_PUREDEBUG "${C_BASE_FLAGS} ${C_PURE_DEBUG_FLAGS}")
set(CMAKE_C_FLAGS_PURERELEASE "${C_BASE_FLAGS} ${C_PURE_RELEASE_FLAGS}")

# Look for a PTHREADS library unless one was already specified
if(NOT CMAKE_THREAD_LIBS_INIT)
  set(CMAKE_THREAD_PREFER_PTHREAD ON)
  find_package(Threads)
  if(NOT CMAKE_USE_PTHREADS_INIT)
    message(WARNING "Could not find a pthreads library, cannot build PMRRR")
    set(FAILED_PMRRR TRUE PARENT_SCOPE)
    return()
  endif()
endif()

# Grab all of the .c and .h files
file(GLOB_RECURSE PMRRR_C RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.c")
file(GLOB_RECURSE PMRRR_H RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.h")
set(PMRRR_SRC "${PMRRR_C};${PMRRR_H}")

# Define the main library and its link libraries
set(MISSING_LAPACK ${MISSING_LAPACK} PARENT_SCOPE)
if(MISSING_LAPACK)
  add_library(pmrrr ${LIBRARY_TYPE} ${PMRRR_SRC})
  target_link_libraries(
    pmrrr lapack-addons ${MPI_C_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}
  )
  install(TARGETS pmrrr lapack-addons DESTINATION lib)
else()
  add_library(pmrrr ${LIBRARY_TYPE} ${PMRRR_SRC})
  target_link_libraries(pmrrr ${MPI_C_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  install(TARGETS pmrrr DESTINATION lib)
endif()
set(CMAKE_THREAD_LIBS_INIT ${CMAKE_THREAD_LIBS_INIT} PARENT_SCOPE)

# Define the header-preparation rules
set(PREPARED_HEADERS)
foreach(HEADER ${PMRRR_H})
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}
            ${CMAKE_CURRENT_BINARY_DIR}/${HEADER}
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}"
  )
  list(APPEND PREPARED_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${HEADER})

  get_filename_component(HEADER_PATH ${HEADER} PATH)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${HEADER} 
          DESTINATION ${HEADER_PATH})
endforeach()
add_custom_target(prepare_pmrrr_headers DEPENDS ${PREPARED_HEADERS})
add_dependencies(pmrrr prepare_pmrrr_headers)

# Make sure that the PMRRR headers can be found
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")

