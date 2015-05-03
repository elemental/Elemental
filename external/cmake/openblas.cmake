#
#  Copyright 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
include(ExternalProject)

if(NOT EL_BUILD_OPENBLAS)
  find_library(OpenBLAS NAMES openblas PATHS ${MATH_PATHS})
  if(OpenBLAS)
    if(CMAKE_COMPILER_IS_GNUCC)
      set(CMAKE_REQUIRED_LIBRARIES ${OpenBLAS} gfortran pthread m)
    else()
      set(CMAKE_REQUIRED_LIBRARIES ${OpenBLAS})
    endif()
    check_function_exists(dgemm   EL_HAVE_DGEMM_OPENBLAS)
    check_function_exists(dgemm_  EL_HAVE_DGEMM_POST_OPENBLAS)
    check_function_exists(dsytrd  EL_HAVE_DSYTRD_OPENBLAS)
    check_function_exists(dsytrd_ EL_HAVE_DSYTRD_POST_OPENBLAS)
    if(EL_HAVE_DGEMM_OPENBLAS OR EL_HAVE_DGEMM_POST_OPENBLAS)
      set(EL_HAVE_OPENBLAS_BLAS TRUE)
    else()
      message(WARNING "OpenBLAS was found as ${OpenBLAS}, but BLAS support was not detected")
    endif()
    if(EL_HAVE_DSYTRD_OPENBLAS OR EL_HAVE_DSYTRD_POST_OPENBLAS)
      set(EL_HAVE_OPENBLAS_LAPACK TRUE)
    else()
      message(WARNING "OpenBLAS was found as ${OpenBLAS}, but LAPACK support was not detected")
    endif()
  endif()
endif()

if(EL_HAVE_OPENBLAS_BLAS AND EL_HAVE_OPENBLAS_LAPACK)
  set(OPENBLAS_LIBS ${OpenBLAS})
  if(CMAKE_COMPILER_IS_GNUCC)
    set(OPENBLAS_LIBS ${OPENBLAS_LIBS} gfortran pthread m)
  endif()
  set(EL_HAVE_OPENBLAS TRUE)
  set(EL_BUILT_OPENBLAS FALSE) 
  message(STATUS "Using OpenBLAS+LAPACK found at ${OpenBLAS}")
elseif(FORTRAN_WORKS)
  if(MSVC)
    set(EL_HAVE_OPENBLAS FALSE)
    set(EL_BUILT_OPENBLAS FALSE)
  else()
    # TODO: Check out from GitHub
    if(NOT DEFINED OPENBLAS_URL)
      set(OPENBLAS_URL https://github.com/xianyi/OpenBLAS.git)
    endif()
    message(STATUS "Will pull OpenBLAS from ${OPENBLAS_URL}")

    set(OPENBLAS_SOURCE_DIR ${PROJECT_BINARY_DIR}/download/OpenBLAS/source)
    set(OPENBLAS_BINARY_DIR ${PROJECT_BINARY_DIR}/download/OpenBLAS/build)

    ExternalProject_Add(project_openblas
      PREFIX ${CMAKE_INSTALL_PREFIX}
      GIT_REPOSITORY ${OPENBLAS_URL}
      GIT_TAG "v0.2.14"
      STAMP_DIR ${OPENBLAS_BINARY_DIR}/stamp
      BUILD_IN_SOURCE 1
      SOURCE_DIR ${OPENBLAS_SOURCE_DIR}
      TMP_DIR    ${OPENBLAS_BINARY_DIR}/tmp
      INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
      CONFIGURE_COMMAND ""
      UPDATE_COMMAND "" 
      BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} CC=${CMAKE_C_COMPILER} FC=${CMAKE_Fortran_COMPILER} libs netlib shared
      INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install PREFIX=<INSTALL_DIR>
    )
    add_dependencies(External project_openblas)

    # Extract the installation directory
    ExternalProject_Get_Property(project_openblas install_dir)

    # Add a target for libopenblas (either shared or static)
    if(BUILD_SHARED_LIBS)
      add_library(libopenblas SHARED IMPORTED)
      set_property(TARGET libopenblas PROPERTY IMPORTED_LOCATION ${install_dir}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}openblas${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
      add_library(libopenblas STATIC IMPORTED)
      set_property(TARGET libopenblas PROPERTY IMPORTED_LOCATION ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif() 

    set(OPENBLAS_LIBS libopenblas)
    if(CMAKE_COMPILER_IS_GNUCC)
      set(OPENBLAS_LIBS ${OPENBLAS_LIBS} gfortran pthread m)
    endif()
    set(EL_HAVE_OPENBLAS TRUE)
    set(EL_BUILT_OPENBLAS TRUE)
  endif()
else()
  set(EL_HAVE_OPENBLAS FALSE)
  set(EL_BUILT_OPENBLAS FALSE)
endif()
