#
#  Copyright 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
include(ExternalProject)
include(ElCheckFunctionExists)
include(ElLibraryName)

# NOTE: Elemental has not yet added support for BLIS with 64-bit integers

find_package(OpenMP)

if(CMAKE_COMPILER_IS_GNUCC)
  if(NOT CMAKE_THREAD_LIBS_INIT)
    set(CMAKE_THREAD_PREFER_PTHREAD ON)
    find_package(Threads REQUIRED)
    if(NOT CMAKE_USE_PTHREADS_INIT)
      message(FATAL_ERROR "Could not find a pthreads library")
    endif()
  endif()
  if(NOT STD_MATH_LIB)
    find_library(STD_MATH_LIB m)
    if(NOT STD_MATH_LIB)
      message(FATAL_ERROR "Could not find standard math library")
    endif()
  endif()
  set(GNU_ADDONS ${CMAKE_THREAD_LIBS_INIT} ${STD_MATH_LIB})
else()
  set(GNU_ADDONS)
endif()

if(NOT EL_FORCE_BLIS_LAPACK_BUILD AND NOT EL_USE_64BIT_BLAS_INTS)
  message(STATUS "Searching for previously installed BLIS+LAPACK")
  # NOTE: The following tests will assume that liblapack is NOT sufficient
  #       by itself and must be supported with a valid BLAS library 
  find_library(BLIS NAMES blis PATHS ${MATH_PATHS})
  find_library(LAPACK NAMES lapack PATHS ${MATH_PATHS})
  if(BLIS AND LAPACK)
    set(CMAKE_REQUIRED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_REQUIRED_LIBRARIES ${LAPACK} ${BLIS} ${GNU_ADDONS})
    El_check_function_exists(dgemm   EL_HAVE_DGEMM)
    El_check_function_exists(dgemm_  EL_HAVE_DGEMM_POST)
    El_check_function_exists(dstegr  EL_HAVE_DSTEGR)
    El_check_function_exists(dstegr_ EL_HAVE_DSTEGR_POST)
    if(EL_HAVE_DGEMM AND (EL_HAVE_DSTEGR OR EL_HAVE_DSTEGR_POST))
      set(USE_FOUND_BLIS_LAPACK TRUE)
      set(EL_BLAS_SUFFIX)
      if(EL_HAVE_DSTEGR)
        set(EL_LAPACK_SUFFIX)
      else()
        set(EL_LAPACK_SUFFIX _)
      endif()
    elseif(EL_HAVE_DGEMM_POST AND (EL_HAVE_DSTEGR OR EL_HAVE_DSTEGR_POST))
      set(USE_FOUND_BLIS_LAPACK TRUE)
      set(EL_BLAS_SUFFIX _)
      if(EL_HAVE_DSTEGR)
        set(EL_LAPACK_SUFFIX)
      else()
        set(EL_LAPACK_SUFFIX _)
      endif()
    endif()
    unset(CMAKE_REQUIRED_LINKER_FLAGS)
    unset(CMAKE_REQUIRED_LIBRARIES)
  endif()
endif()

if(USE_FOUND_BLIS_LAPACK)
  set(BLIS_LAPACK_LIBS ${LAPACK} ${BLIS} ${GNU_ADDONS})
  set(BLIS_LAPACK_LIBS_AT_CONFIG ${BLIS_LAPACK_LIBS})
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${BLIS_LAPACK_LIBS})
  set(EL_BUILT_BLIS_LAPACK FALSE)
  set(EL_HAVE_BLIS_LAPACK TRUE)
elseif(EL_HAVE_F90_INTERFACE AND
    NOT EL_USE_64BIT_BLAS_INTS AND
    NOT MSVC AND
    NOT EL_CONFIG_TIME_BLIS_LAPACK)
  if(NOT DEFINED LAPACK_URL)
    set(LAPACK_URL https://github.com/scibuilder/lapack)
  endif()
  message(STATUS "Will download LAPACK from ${LAPACK_URL}")

  set(LAPACK_SOURCE_DIR ${PROJECT_BINARY_DIR}/download/blis_lapack/source)
  set(LAPACK_BINARY_DIR ${PROJECT_BINARY_DIR}/download/blis_lapack/build)

  # Provide a way for Elemental to pass down the BLIS architecture
  if(BLIS_ARCH)
    set(BLIS_ARCH_COMMAND -D BLIS_ARCH=${BLIS_ARCH})
  else()
    set(BLIS_ARCH_COMMAND)
  endif()

  # Set up a target for building and installing OpenBLAS
  ExternalProject_Add(project_blis_lapack
    PREFIX ${CMAKE_INSTALL_PREFIX}
    GIT_REPOSITORY ${LAPACK_URL}
    STAMP_DIR  ${LAPACK_BINARY_DIR}/stamp
    SOURCE_DIR ${LAPACK_SOURCE_DIR}
    BINARY_DIR ${LAPACK_BINARY_DIR}
    TMP_DIR    ${LAPACK_BINARY_DIR}/tmp
    UPDATE_COMMAND ""
    CMAKE_ARGS 
      -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -D CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -D CMAKE_C_FLAGS=${CMAKE_C_FLAGS}
      -D CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
      -D USE_OPTIMIZED_BLAS=ON
      -D FORCE_BLIS_BUILD=ON ${BLIS_ARCH_COMMAND}
      -D CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
      -D CMAKE_MACOSX_RPATH=${CMAKE_MACOSX_RPATH}
      -D CMAKE_SKIP_RPATH=${CMAKE_SKIP_RPATH}
      -D CMAKE_SKIP_BUILD_RPATH=${CMAKE_SKIP_BUILD_RPATH}
      -D CMAKE_BUILD_WITH_INSTALL_RPATH=${CMAKE_BUILD_WITH_INSTALL_RPATH}
      -D CMAKE_INSTALL_RPATH_USE_LINK_PATH=${CMAKE_INSTALL_RPATH_USE_LINK_PATH} 
      -D CMAKE_INSTALL_RPATH=${CMAKE_INSTALL_RPATH}
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
  )

  # Extract the source and install directories
  ExternalProject_Get_Property(project_blis_lapack source_dir install_dir)

  # Add a target for libblis 
  add_library(libblis ${LIBRARY_TYPE} IMPORTED)
  El_library_name(blis_name blis)
  set(BLIS_LIB ${install_dir}/lib/${blis_name})
  set_property(TARGET libblis PROPERTY IMPORTED_LOCATION ${BLIS_LIB})

  # Add a target for liblapack
  add_library(liblapack ${LIBRARY_TYPE} IMPORTED)
  El_library_name(lapack_name lapack)
  set(LAPACK_LIB ${install_dir}/lib/${lapack_name})
  set_property(TARGET liblapack PROPERTY IMPORTED_LOCATION ${LAPACK_LIB})

  set(BLIS_LAPACK_LIBS ${LAPACK_LIB} ${BLIS_LIB} ${GNU_ADDONS})
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${BLIS_LAPACK_LIBS})
  set(EL_BUILT_BLIS_LAPACK TRUE)
  set(EL_HAVE_BLIS_LAPACK TRUE)
elseif(EL_HAVE_F90_INTERFACE AND NOT EL_USE_64BIT_BLAS_INTS)
  set(EL_PLAN_BLIS_LAPACK TRUE)
  set(EL_BUILT_BLIS_LAPACK FALSE)
  set(EL_HAVE_BLIS_LAPACK FALSE)
else()
  set(EL_BUILT_BLIS_LAPACK FALSE)
  set(EL_HAVE_BLIS_LAPACK FALSE)
endif()
