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

if(EL_USE_64BIT_BLAS_INTS)
  set(OPENBLAS_SUFFIX 64)
else()
  set(OPENBLAS_SUFFIX)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(GFORTRAN_PATHS /usr/lib
                     /usr/lib64
                     /usr/local/lib
                     /usr/local/lib64
                     /usr/lib/openmpi/lib
                     /usr/lib/gcc/x86_64-linux-gnu/4.8
                     /usr/lib/gcc/x86_64-linux-gnu/4.9
                     /lib/x86_64-linux-gnu
                     /usr/lib/x86_64-linux-gnu
                     /usr/lib/openblas-base
                     /usr/lib64/openblas-base
                     /usr/local/lib/gcc/4.8
                     /usr/local/lib/gcc/4.9)
  if(NOT GFORTRAN_LIB)
    find_library(GFORTRAN_LIB NAMES gfortran PATHS ${GFORTRAN_PATHS})
    if(NOT GFORTRAN_LIB)
      if(APPLE)
        message(FATAL_ERROR "Could not find gfortran library; please consider setting the GFORTRAN_LIB variable. If you installed gfortran via homebrew, please see the issue filed at https://github.com/Homebrew/homebrew/issues/8539")
      else()
        message(FATAL_ERROR "Could not find gfortran library; please consider setting the GFORTRAN_LIB variable.")
      endif()
    endif()
  endif()
  if(NOT CMAKE_THREAD_LIBS_INIT)
    set(CMAKE_THREAD_PREFER_PTHREAD ON)
    find_package(Threads)
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
  set(GNU_ADDONS ${GFORTRAN_LIB} ${CMAKE_THREAD_LIBS_INIT} ${STD_MATH_LIB})
else()
  set(GNU_ADDONS)
endif()

if(NOT EL_FORCE_OPENBLAS_BUILD)
  message(STATUS "Searching for previously installed OpenBLAS+LAPACK")
  find_library(OpenBLAS NAMES openblas${OPENBLAS_SUFFIX} PATHS ${MATH_PATHS})
  if(OpenBLAS)
    set(CMAKE_REQUIRED_LIBRARIES ${OpenBLAS} ${GNU_ADDONS})
    El_check_function_exists(dgemm${OPENBLAS_SUFFIX}   EL_HAVE_DGEMM_OPENBLAS)
    El_check_function_exists(dgemm${OPENBLAS_SUFFIX}_  EL_HAVE_DGEMM_POST_OPENBLAS)
    El_check_function_exists(dsytrd${OPENBLAS_SUFFIX}  EL_HAVE_DSYTRD_OPENBLAS)
    El_check_function_exists(dsytrd${OPENBLAS_SUFFIX}_ EL_HAVE_DSYTRD_POST_OPENBLAS)
    if(EL_HAVE_DGEMM_OPENBLAS OR EL_HAVE_DGEMM_POST_OPENBLAS)
      set(EL_HAVE_OPENBLAS_BLAS TRUE)
      if(EL_HYBRID)
        message("Found OpenBLAS as ${OpenBLAS}; you may want to experiment with the environment variable OPENBLAS_NUM_THREADS")
      else()
        message("Found OpenBLAS as ${OpenBLAS}; you may want to set the environment variable OPENBLAS_NUM_THREADS=1")
      endif()
    else()
      message(WARNING "OpenBLAS was found as ${OpenBLAS}, but BLAS support was not detected")
    endif()
    if(EL_HAVE_DSYTRD_OPENBLAS OR EL_HAVE_DSYTRD_POST_OPENBLAS)
      set(EL_HAVE_OPENBLAS_LAPACK TRUE)
    else()
      message(WARNING "OpenBLAS was found as ${OpenBLAS}, but LAPACK support was not detected")
    endif()
    if(EL_HAVE_OPENBLAS_BLAS AND EL_HAVE_OPENBLAS_LAPACK)
      if(EL_HAVE_DGEMM_OPENBLAS)
        set(EL_BLAS_SUFFIX ${OPENBLAS_SUFFIX})
      else()
        set(EL_BLAS_SUFFIX ${OPENBLAS_SUFFIX}_)
      endif()
      if(EL_HAVE_DSYTRD_OPENBLAS)
        set(EL_LAPACK_SUFFIX ${OPENBLAS_SUFFIX})
      else()
        set(EL_LAPACK_SUFFIX ${OPENBLAS_SUFFIX}_)
      endif()
    endif()
    unset(CMAKE_REQUIRED_LIBRARIES)
  endif()
endif()

if(EL_HAVE_OPENBLAS_BLAS AND EL_HAVE_OPENBLAS_LAPACK)
  set(OPENBLAS_LIBS ${OpenBLAS} ${GNU_ADDONS})
  set(OPENBLAS_LIBS_AT_CONFIG ${OPENBLAS_LIBS})
  set(EL_HAVE_OPENBLAS TRUE)
  set(EL_BUILT_OPENBLAS FALSE) 
  if(EL_HAVE_DGEMM_OPENBLAS
  message(STATUS "Using OpenBLAS+LAPACK found at ${OpenBLAS}")
elseif(FORTRAN_WORKS AND NOT MSVC AND NOT EL_CONFIG_TIME_OPENBLAS)
  if(NOT DEFINED OPENBLAS_URL)
    set(OPENBLAS_URL https://github.com/xianyi/OpenBLAS.git)
  endif()
  message(STATUS "Will pull OpenBLAS from ${OPENBLAS_URL}")

  set(OPENBLAS_SOURCE_DIR ${PROJECT_BINARY_DIR}/download/OpenBLAS/source)
  set(OPENBLAS_BINARY_DIR ${PROJECT_BINARY_DIR}/download/OpenBLAS/build)

  if(APPLE)
    if(NOT OPENBLAS_ARCH_COMMAND)
      # This is a hack but is a good default for modern Mac's
      set(OPENBLAS_ARCH_COMMAND TARGET=SANDYBRIDGE)
    endif()
  else()
    if(NOT OPENBLAS_ARCH_COMMAND)
      set(OPENBLAS_ARCH_COMMAND)
    endif()
  endif()

  if(NOT OPENBLAS_THREAD_COMMAND)
    if(EL_HYBRID)
      set(OPENBLAS_THREAD_COMMAND USE_OPENMP=1)
    else()
      set(OPENBLAS_THREAD_COMMAND USE_THREAD=0)
    endif()
  endif()

  if(EL_USE_64BIT_BLAS_INTS)
    set(OPENBLAS_INTERFACE_COMMAND INTERFACE64=1 SYMBOLSUFFIX=64)
  else()
    set(OPENBLAS_INTERFACE_COMMAND INTERFACE64=0)
  endif()

  ExternalProject_Add(project_openblas
    PREFIX ${CMAKE_INSTALL_PREFIX}
    GIT_REPOSITORY ${OPENBLAS_URL}
    GIT_TAG "v0.2.15"
    STAMP_DIR ${OPENBLAS_BINARY_DIR}/stamp
    BUILD_IN_SOURCE 1
    SOURCE_DIR ${OPENBLAS_SOURCE_DIR}
    TMP_DIR    ${OPENBLAS_BINARY_DIR}/tmp
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    CONFIGURE_COMMAND ""
    UPDATE_COMMAND "" 
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} \
      CC=${CMAKE_C_COMPILER} \
      FC=${CMAKE_Fortran_COMPILER} \
      ${OPENBLAS_THREAD_COMMAND} \
      ${OPENBLAS_ARCH_COMMAND} \
      ${OPENBLAS_INTERFACE_COMMAND} \
      libs netlib shared
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install PREFIX=<INSTALL_DIR>
  )

  # Extract the installation directory
  ExternalProject_Get_Property(project_openblas install_dir)

  # Add a target for libopenblas(64) (either shared or static)
  add_library(libopenblas ${LIBRARY_TYPE} IMPORTED)
  El_library_name(openblas_name openblas${OPENBLAS_SUFFIX})
  set(OPENBLAS_LIB ${install_dir}/lib/${openblas_name})
  set_property(TARGET libopenblas PROPERTY IMPORTED_LOCATION ${OPENBLAS_LIB})

  set(OPENBLAS_LIBS ${OPENBLAS_LIB} ${GNU_ADDONS})
  set(EL_HAVE_OPENBLAS TRUE)
  set(EL_BUILT_OPENBLAS TRUE)
elseif(FORTRAN_WORKS AND NOT MSVC)
  set(EL_PLAN_OPENBLAS TRUE)
  set(EL_HAVE_OPENBLAS FALSE)
  set(EL_BUILT_OPENBLAS FALSE)
else()
  set(EL_HAVE_OPENBLAS FALSE)
  set(EL_BUILT_OPENBLAS FALSE)
endif()
