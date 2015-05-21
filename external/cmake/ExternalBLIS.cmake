#
#  Copyright 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
include(ExternalProject)
include(ElCheckFunctionExists)

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

if(NOT EL_BUILD_BLIS)
  find_library(BLIS NAMES blis PATHS ${MATH_PATHS})
  if(BLIS)
    set(CMAKE_REQUIRED_LINKER_FLAGS)
    set(CMAKE_REQUIRED_LIBRARIES ${BLIS} ${GNU_ADDONS})
    El_check_function_exists(dgemm_ EL_HAVE_DGEMM_BLIS)
    set(CMAKE_REQUIRED_LINKER_FLAGS ${OpenMP_C_FLAGS})
    El_check_function_exists(dgemm_ EL_HAVE_DGEMM_OPENMP_BLIS)
    if(EL_HAVE_DGEMM_POST_BLIS OR EL_HAVE_DGEMM_OPENMP_BLIS)
      set(EL_HAVE_BLIS_BLAS TRUE)
    else()
      message(WARNING "BLIS was found as ${BLIS}, but BLAS support was not detected")
    endif()
    set(CMAKE_REQUIRED_LIBRARIES)
  endif()
endif()

if(EL_HAVE_BLIS_BLAS)
  set(BLIS_LIBS ${BLIS} ${GNU_ADDONS})
  if(EL_HAVE_DGEMM_POST_BLIS)
    set(BLIS_LINK_FLAGS)
  else()
    set(BLIS_LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
  set(EL_HAVE_BLIS TRUE)
  set(EL_BUILT_BLIS FALSE) 
  message(STATUS "Using BLIS found at ${BLIS}")
elseif(NOT MSVC)
  if(NOT DEFINED BLIS_URL)
    set(BLIS_URL https://github.com/flame/blis.git)
  endif()
  message(STATUS "Will pull BLIS from ${BLIS_URL}")

  set(BLIS_SOURCE_DIR ${PROJECT_BINARY_DIR}/download/blis/source)
  set(BLIS_BINARY_DIR ${PROJECT_BINARY_DIR}/download/blis/build)

  if(NOT BLIS_ARCH)
    if(APPLE)
      # This is a hack but is a good default since Mac detection can be poor
      set(BLIS_ARCH sandybridge)
    else()
      set(BLIS_ARCH auto)
    endif()
  endif()

  ExternalProject_Add(project_blis
    PREFIX ${CMAKE_INSTALL_PREFIX}
    GIT_REPOSITORY ${BLIS_URL}
    STAMP_DIR ${BLIS_BINARY_DIR}/stamp
    BUILD_IN_SOURCE 1
    SOURCE_DIR ${BLIS_SOURCE_DIR}
    TMP_DIR    ${BLIS_BINARY_DIR}/tmp
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    CONFIGURE_COMMAND ""
    UPDATE_COMMAND "" 
    BUILD_COMMAND ${BLIS_SOURCE_DIR}/configure -p <INSTALL_DIR> ${BLIS_ARCH}
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
  )

  # Extract the installation directory
  ExternalProject_Get_Property(project_blis install_dir)

  # Add a target for libblis 
  add_library(libblis STATIC IMPORTED)
  set(BLIS_LIB ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}blis${CMAKE_STATIC_LIBRARY_SUFFIX})
  set_property(TARGET libblis PROPERTY IMPORTED_LOCATION ${BLIS_LIB})

  set(BLIS_LIBS ${BLIS_LIB} ${GNU_ADDONS})
  # If we used the 'auto' configuration, we don't know if BLIS used OpenMP or
  # not, and so we might as well add it as a link flag
  set(BLIS_LINK_FLAGS ${OpenMP_CXX_FLAGS})
  set(EL_HAVE_BLIS TRUE)
  set(EL_BUILT_BLIS TRUE)
else()
  set(EL_HAVE_BLIS FALSE)
  set(EL_BUILT_BLIS FALSE)
endif()
