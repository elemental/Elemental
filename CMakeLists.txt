#
#  Copyright 2009-2016, Jack Poulson
#  All rights reserved.
#
#  Copyright 2012, Jed Brown
#  All rights reserved.
#
#  Copyright 2015, Jake Bolewski and Jiahao Chen
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
if(APPLE)
  # RPATH support for OS X was introduced in this version
  cmake_minimum_required(VERSION 2.8.12)
else()
  # It is likely/possible that this version could be decreased, but this version
  # was chosen to support the default version of CMake in CentOS 7
  cmake_minimum_required(VERSION 2.8.11)
endif()

project(Elemental C CXX)
set(EL_VERSION_MAJOR 0)
set(EL_VERSION_MINOR 86-dev)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# Configuration options
# =====================

# Basic options
# -------------

# "CMAKE_BUILD_TYPE"
#   Elemental requires it to be one of the following two choices:
#     1. Debug: Meant for development and debugging
#     2. Release: Meant for production runs
#   If either MinSizeRel or RelWithDebInfo are specified, then Elemental falls
#   back to Release mode.
#
#   If "CMAKE_BUILD_TYPE" is not defined, the default is Release mode.

# "CXX_FLAGS"
#   Optimization and debugging/symbol flags

# "OpenMP_CXX_FLAGS"
#   Overrides the default compile flags for adding OpenMP support to CXX code

# TODO: Add list of MPI variables here, such as "MPI_CXX_COMPILE_FLAGS"

# TODO: Add discussion of "MATH_LIBS" and "REFERENCE_ROOT"

get_property(SUPPORTS_SHARED GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
if(SUPPORTS_SHARED)
  option(BUILD_SHARED_LIBS "Build shared libraries?" ON)
else()
  option(BUILD_SHARED_LIBS "Build shared libraries?" OFF)
endif()
if(BUILD_SHARED_LIBS AND NOT SUPPORTS_SHARED)
  message("Elemental was configured to build shared libraries, but they are not supported on the target platform: BUILD_SHARED_LIBS will be manually set to OFF")
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries?")
endif()

include(GNUInstallDirs)
if(NOT CMAKE_INSTALL_INCLUDEDIR)
  set(CMAKE_INSTALL_INCLUDEDIR "include")
endif()
if(NOT CMAKE_INSTALL_LIBDIR) 
  set(CMAKE_INSTALL_LIBDIR "lib") 
endif()
if(NOT CMAKE_INSTALL_BINDIR)
  set(CMAKE_INSTALL_BINDIR "bin")
endif()
if(NOT CMAKE_INSTALL_DATADIR)
  set(CMAKE_INSTALL_DATADIR "share") 
endif()
if(NOT CMAKE_INSTALL_DOCDIR)
  set(CMAKE_INSTALL_DOCDIR "share/doc/elemental") 
endif()
if(NOT CMAKE_INSTALL_SYSCONFDIR)
  set(CMAKE_INSTALL_SYSCONFDIR "etc") 
endif()

# Whether or not to attempt to use OpenMP within hot-spots in Elemental
# (NOTE: This option is not actively maintained)
option(EL_HYBRID "Make use of OpenMP within MPI packing/unpacking" OFF)

option(EL_C_INTERFACE "Build C interface" ON)

if(BUILD_SHARED_LIBS AND EL_C_INTERFACE)
  # Define PYTHON_SITE_PACKAGES if you want to install the python package
  # somewhere other than the default system-wide location (e.g., within your
  # home directory). Alternatively, if you do NOT define PYTHON_SITE_PACKAGES 
  # and instead define INSTALL_PYTHON_INTO_USER_SITE=ON, Elemental will
  # install its Python interface into site.USER_SITE in your home directory.
  option(INSTALL_PYTHON_PACKAGE "Install Python interface? (requires Python 2.x)" OFF)

  set(Python_ADDITIONAL_VERSIONS 2.7)
  include(FindPythonInterp) #Check for Python version

  #Only Python 2.x supported
  if (NOT PYTHON_VERSION_MAJOR EQUAL 2)
    set(INSTALL_PYTHON_PACKAGE FALSE)
  endif()
endif()

# Whether or not to use Qt5 for visualization if it is found.
# Setting Qt5Widgets_DIR, Qt5_DIR, or Qt5_LIBDIR overrides this variable
# and is often necessary anyway.
option(EL_USE_QT5 "Attempt to use Qt5?" OFF)

option(EL_EXAMPLES "Build simple examples?" OFF)
option(EL_TESTS "Build performance and correctness tests?" OFF)
option(EL_EXPERIMENTAL "Build experimental code" OFF)

# Attempt to use 64-bit integers?
option(EL_USE_64BIT_INTS "Use 64-bit integers for El indexing" OFF)
option(EL_USE_64BIT_BLAS_INTS "Use 64-bit integers for BLAS/LAPACK" OFF)
if(EL_USE_64BIT_BLAS_INTS)
  set(EL_USE_64BIT_INTS ON)
endif()

option(EL_DISABLE_MKL "Do not use MKL even if installed?" OFF)
option(EL_DISABLE_MKL_CSRMV "Avoid MKL's CSR mat-vec?" ON)

option(EL_DISABLE_BLIS_LAPACK "Avoid BLIS+LAPACK" OFF)
option(EL_PREFER_BLIS_LAPACK "Prefer BLIS+LAPACK over alternatives?" OFF)
option(EL_FORCE_BLIS_LAPACK_BUILD 
  "Avoid search for BLIS and instead buiild?" OFF)

option(EL_DISABLE_OPENBLAS "Avoid OpenBLAS" OFF)
option(EL_PREFER_OPENBLAS "Prefer OpenBLAS over alternatives?" OFF)
option(EL_FORCE_OPENBLAS_BUILD
  "Avoid search for OpenBLAS and instead build?" OFF)

# Ugh. Elemental currently supports ScaLAPACK building its own custom
# BLAS and LAPACK if Elemental is using a custom symbol suffix for its own.
if(EL_USE_64BIT_BLAS_INTS AND
   (NOT CUSTOM_BLAS_SUFFIX OR NOT CUSTOM_LAPACK_SUFFIX))
  option(EL_DISABLE_SCALAPACK "Avoid ScaLAPACK" ON)
else()
  option(EL_DISABLE_SCALAPACK "Avoid ScaLAPACK" OFF)
endif()
option(EL_FORCE_SCALAPACK_BUILD
  "Avoid search for ScaLAPACK and instead build?" OFF)

# By default, a search is performed for ParMETIS via find_package, and, if this
# search fails, ParMETIS is downloaded from the official location:
#    http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
# If EL_FORCE_PARMETIS_BUILD is turned on, then no initial call to
# find_package is performed.
#
# Should you want to manually specify a ParMETIS installation, you can set the
# variables PARMETIS_INCLUDE_DIRS and PARMETIS_LIBRARIES
#
option(EL_FORCE_PARMETIS_BUILD "Force build of ParMETIS?" OFF)
# ParMETIS cannot be commercially used without a license and so users interested
# in using Elemental in commercial applications should disable ParMETIS
option(EL_DISABLE_PARMETIS "Disable ParMETIS?" OFF)

# If ParMETIS is disabled, a search is performed for METIS via find_package, 
# and, if this search fails, METIS is cloned from the git repo 
#     https://github.com/poulson/Metis.git.
# If EL_FORCE_METIS_BUILD is turned on, then no initial call to find_package is 
# performed.
#
# Should you want to manually specify a METIS installation, you can set the
# variables METIS_INCLUDE_DIRS and METIS_LIBRARIES
option(EL_FORCE_METIS_BUILD "Force a build of METIS?" OFF)

# Advanced options
# ----------------

if(APPLE)
  option(HAVE_PACKAGEMAKER "Have Mac OS X PackageMaker?" OFF)
  mark_as_advanced(HAVE_PACKAGEMAKER)
endif()

# Whether or not to have the Memory class zero initialize what it allocates.
# If valgrind was detected and is running, this will be forced anyway.
option(EL_ZERO_INIT "Initialize buffers to zero by default?" OFF)
mark_as_advanced(EL_ZERO_INIT)

option(EL_DISABLE_VALGRIND "Prevent Elemental from looking for valgrind?" OFF)
mark_as_advanced(EL_DISABLE_VALGRIND)

option(EL_USE_CUSTOM_ALLTOALLV "Avoid MPI_Alltoallv for performance reasons" ON)
mark_as_advanced(EL_USE_CUSTOM_ALLTOALLV)

# Since it is surprisingly common for MPI libraries to have bugs in their
# support for complex data, the following option forces Elemental to cast
# all possible MPI communications in terms of twice as many real units of data.
option(EL_AVOID_COMPLEX_MPI "Avoid potentially buggy complex MPI routines" ON)
mark_as_advanced(EL_AVOID_COMPLEX_MPI)

# At one point, a bug was found in IBM's C++ compiler for Blue Gene/P,
# where OpenMP statements of the form a[i] += alpha b[i], with complex data,
# would segfault and/or return incorrect results
option(EL_AVOID_OMP_FMA "Avoid a bug in the IBM compilers." OFF)
mark_as_advanced(EL_AVOID_OMP_FMA)

# Due to a subtle flaw in the Blue Gene/P extensions for MPICH2, treating
# floating-point data as a collection of byte-sized objects results in a
# better algorithm being chosen for MPI_Allgather. This should not effect
# performance on most machines.
option(EL_USE_BYTE_ALLGATHERS "Avoid BG/P allgather performance bug." ON)
mark_as_advanced(EL_USE_BYTE_ALLGATHERS)

# If MPI_Reduce_scatter_block doesn't exist, perform it by composing
# MPI_Allreduce and std::memcpy rather than MPI_Reduce and MPI_Scatter
option(EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
       "AllReduce based block MPI_Reduce_scatter" OFF)
mark_as_advanced(EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE)

# Print a warning any time a redistribution is performed which unpacks a
# large amount of data with a non-unit stride
option(EL_CACHE_WARNINGS "Warns when using cache-unfriendly routines" OFF)
mark_as_advanced(EL_CACHE_WARNINGS)

# Print a warning when an improperly aligned redistribution is performed,
# i.e., if an unnecessary permutation communication stage must take place
option(EL_UNALIGNED_WARNINGS
       "Warn when performing unaligned redistributions" OFF)
mark_as_advanced(EL_UNALIGNED_WARNINGS)

# Print a warning if an opportunity was missed to implement a redistribution
# approach specifically for vectors (instead of matrices)
option(EL_VECTOR_WARNINGS
       "Warn when vector redistribution chances are missed" OFF)
mark_as_advanced(EL_VECTOR_WARNINGS)

# Build logic
# ===========

# Add the header file include directories first (so they have priority)
# ---------------------------------------------------------------------
message(STATUS 
  "Appending ${PROJECT_SOURCE_DIR}/include for Elemental's source includes")
include_directories("${PROJECT_SOURCE_DIR}/include")
message(STATUS 
  "Appending ${PROJECT_BINARY_DIR}/include for Elemental's binary includes")
include_directories("${PROJECT_BINARY_DIR}/include")

if(BUILD_SHARED_LIBS)
  set(LIBRARY_TYPE SHARED)
else()
  set(LIBRARY_TYPE STATIC)
endif()

# Elemental must be built "out-of-source", so we start by ensuring that the
# source and build directories are different.
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source build attempted; please clean the CMake cache and then switch to an out-of-source build, e.g.,\nrm CMakeCache.txt && rm -Rf CMakeFiles/\nmkdir build/ && cd build/ && cmake ..")
endif()

# Get the Git revision
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Ensure that the build type is set to either Release or Debug
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  # This option is okay as-is
  set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # This option is okay as-is
  set(CMAKE_BUILD_TYPE Debug)
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  message(WARNING "RelWithDebInfo not supported; switching to Release")
  set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
  message(WARNING "MinSizeRel not supported; switching to Release")
  set(CMAKE_BUILD_TYPE Release)
else()
  message(WARNING "Build mode not specified, defaulting to Release build.")
  set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(EL_RELEASE TRUE)
else()
  set(EL_RELEASE FALSE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} UPPER_BUILD_TYPE)

# Set the basic compile flags from the build type
if(NOT WIN32)
  set(LANGUAGES CXX C Fortran)
  foreach(LANG ${LANGUAGES})
    if(NOT ${LANG}_FLAGS)
      if(EL_RELEASE)
        set(${LANG}_FLAGS "-O3" CACHE STRING
          "${LANG} optimization flags for ${UPPER_BUILD_TYPE} mode")
      else()
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
          set(${LANG}_FLAGS "-O2 -g -Wall -Wunused-variable -Wunused-but-set-variable -Wunused-local-typedefs" CACHE STRING "${LANG} optimization/debug flags for ${UPPER_BUILD_TYPE} mode")
        else()
          set(${LANG}_FLAGS "-O2 -g" CACHE STRING
            "${LANG} optimization/debug flags for ${UPPER_BUILD_TYPE} mode")
        endif()
      endif()
    endif()
  endforeach()
endif()

# Disable the requirement for safe exception handling in MSVC builds (otherwise there are conflicts with OpenBLAS)
if(MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SAFESEH:NO")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /SAFESEH:NO")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /SAFESEH:NO")
  set(CMAKE_EXE_LINKER_FLAGS_${UPPER_BUILD_TYPE} "${CMAKE_EXE_LINKER_FLAGS_${UPPER_BUILD_TYPE}} /SAFESEH:NO")
  set(CMAKE_SHARED_LINKER_FLAGS_${UPPER_BUILD_TYPE} "${CMAKE_SHARED_LINKER_FLAGS_${UPPER_BUILD_TYPE}} /SAFESEH:NO")
  set(CMAKE_MODULE_LINKER_FLAGS_${UPPER_BUILD_TYPE} "${CMAKE_MODULE_LINKER_FLAGS_${UPPER_BUILD_TYPE}} /SAFESEH:NO")
endif()

# Handle RPATHs for Mac
# ---------------------
if(APPLE)
  set(CMAKE_MACOSX_RPATH TRUE)

  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  # Add automatically determined parts of RPATH which point to directories 
  # outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES 
    "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif()
endif()

# Detect Fortran
# --------------
include(detect/Fortran)

# Detect C++11
# ------------
include(detect/CXX)
set(CXX_FLAGS "${CXX_FLAGS} ${CXX11_COMPILER_FLAGS}")

# Detect MPI
# ----------
include(detect/MPI)
message(STATUS "Appending ${MPI_C_INCLUDE_PATH} for MPI headers")
include_directories(${MPI_C_INCLUDE_PATH})
set(CXX_FLAGS "${CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")

# Detect OpenMP
# -------------
include(detect/OpenMP)
if(EL_HYBRID)
  set(CXX_FLAGS "${CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Detect Qt5
# ----------
include(detect/Qt5)
if(EL_HAVE_QT5)
  set(CXX_FLAGS "${CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
endif()

# Allow valgrind support if possible (if running valgrind, explicitly zero init)
# ------------------------------------------------------------------------------
if(NOT EL_DISABLE_VALGRIND)
  include(FindValgrind)
  if(VALGRIND_FOUND)
    message(STATUS "Appending ${VALGRIND_INCLUDE_DIR} for Valgrind headers")
    include_directories(${VALGRIND_INCLUDE_DIR})
    set(EL_HAVE_VALGRIND TRUE)
  endif()
endif()

# Import several external math libraries (e.g., BLAS/LAPACK/METIS/ParMETIS)
# -------------------------------------------------------------------------
include(external_projects/ElMath)

# Now that testing is done, set the CXX build and link flags
# ----------------------------------------------------------
set(CMAKE_CXX_FLAGS_${UPPER_BUILD_TYPE} ${CXX_FLAGS})
set(EL_LINK_FLAGS)
if(MPI_CXX_LINK_FLAGS)
  set(EL_LINK_FLAGS ${MPI_CXX_LINK_FLAGS})
endif()
if(EL_BUILT_BLIS_LAPACK)
  if(EL_LINK_FLAGS)
    set(EL_LINK_FLAGS "${EL_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")
  else()
    set(EL_LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endif()

# Elemental's mod's of Parallel Multiple Relatively Robust Representations
# ------------------------------------------------------------------------
add_subdirectory(external/pmrrr)
if(EL_BUILT_SCALAPACK)
  add_dependencies(pmrrr project_scalapack)
else()
  if(EL_BUILT_BLIS_LAPACK)
    add_dependencies(pmrrr project_blis_lapack)
  endif()
  if(EL_BUILT_OPENBLAS)
    add_dependencies(pmrrr project_openblas)
  endif()
endif()

# Elemental's modifications of (a subset of) SuiteSparse
# ------------------------------------------------------
add_subdirectory(external/suite_sparse)
include_directories(external/suite_sparse/include)

# Create a file which can be included in Makefile's.
# This is meant to be analogous to PETSc's 'conf/petscvariables' file
include(CreateElVars)

# The main library
# ================

# Get the combined list of source files for Elemental
# ---------------------------------------------------
file(GLOB_RECURSE EL_C_CPP_SOURCE RELATIVE "${PROJECT_SOURCE_DIR}"
  "src/*.c" "src/*.cpp" "include/*.h" "include/*.hpp" )
if(EL_HAVE_QT5)
  set(EL_HEADERS_PREMOC
    "include/El/io/DisplayWindow-premoc.hpp;include/El/io/ComplexDisplayWindow-premoc.hpp")
  qt_wrap_cpp(El EL_MOC_SRC ${EL_HEADERS_PREMOC})
  message(STATUS "Appending ${Qt5Widgets_INCLUDE_DIRS} for Qt5 headers")
  include_directories(${Qt5Widgets_INCLUDE_DIRS})
  add_definitions(${Qt5Widgets_DEFINITIONS})

  # Qt5Widgets_DIR = Qt5_LIBDIR/cmake/Qt5Widgets
  get_filename_component(Qt5_CMAKEDIR ${Qt5Widgets_DIR} PATH)
  get_filename_component(Qt5_LIBDIR ${Qt5_CMAKEDIR} PATH)

  set(EL_C_CPP_FILES "${EL_C_CPP_SOURCE};${EL_MOC_SRC}")
else()
  set(EL_C_CPP_FILES "${EL_C_CPP_SOURCE}")
endif()

# Handle the header preparation and installation
# ----------------------------------------------
configure_file("${PROJECT_SOURCE_DIR}/cmake/configure_files/config.h.in"
               "${PROJECT_BINARY_DIR}/include/El/config.h")
install(FILES "${PROJECT_BINARY_DIR}/include/El/config.h"
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/El)
add_custom_target(El_config DEPENDS "${PROJECT_BINARY_DIR}/include/El/config.h")
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/"
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} 
  	FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

# Define the main library for Elemental
# -------------------------------------
add_library(El ${EL_C_CPP_FILES})
add_dependencies(El El_config)
if(EL_BUILT_SCALAPACK)
  add_dependencies(El project_scalapack)
else()
  if(EL_BUILT_BLIS_LAPACK)
    add_dependencies(El project_blis_lapack)
  endif()
  if(EL_BUILT_OPENBLAS)
    add_dependencies(El project_openblas)
  endif()
endif()
if(EL_BUILT_METIS)
  add_dependencies(El project_metis)
endif()
if(EL_BUILT_PARMETIS)
  add_dependencies(El project_parmetis)
endif()
set(LINK_LIBS pmrrr ElSuiteSparse
  ${EXTERNAL_LIBS} ${MATH_LIBS} ${MPI_CXX_LIBRARIES})
if(EL_HAVE_QT5)
  set(LINK_LIBS ${LINK_LIBS} ${Qt5Widgets_LIBRARIES})
endif()
target_link_libraries(El ${LINK_LIBS})
if(EL_LINK_FLAGS)
  set_target_properties(El PROPERTIES LINK_FLAGS ${EL_LINK_FLAGS})
endif()
install(TARGETS El EXPORT ElementalTargets DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Decide if/how to install the Python interface
# ---------------------------------------------
if(BUILD_SHARED_LIBS)
  # Create the environment file for Python which hardcodes the installation dir
  configure_file("${PROJECT_SOURCE_DIR}/cmake/configure_files/environment.py.in"
                 "${PROJECT_BINARY_DIR}/python/core/environment.py")
  if(INSTALL_PYTHON_PACKAGE)
    if(NOT PYTHON_SITE_PACKAGES)
      if(INSTALL_PYTHON_INTO_USER_SITE)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
          "import site; print site.USER_SITE"
          OUTPUT_VARIABLE PYTHON_SITE_PACKAGES 
          OUTPUT_STRIP_TRAILING_WHITESPACE)
      else()
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
          "from distutils.sysconfig import get_python_lib; print get_python_lib()"
          OUTPUT_VARIABLE PYTHON_SITE_PACKAGES 
          OUTPUT_STRIP_TRAILING_WHITESPACE)
      endif()  
      set(PYTHON_SITE_PACKAGES "${PYTHON_SITE_PACKAGES}"
        CACHE PATH "Destination for Python site package directory")
    endif()
    install(DIRECTORY "${PROJECT_SOURCE_DIR}/python/"
      DESTINATION "${PYTHON_SITE_PACKAGES}/El"
      FILES_MATCHING PATTERN "*.py")
    install(FILES "${PROJECT_BINARY_DIR}/python/core/environment.py"
      DESTINATION "${PYTHON_SITE_PACKAGES}/El/core")
  else()
    install(DIRECTORY "${PROJECT_SOURCE_DIR}/python/"
      DESTINATION python/El
      FILES_MATCHING PATTERN "*.py")
    install(FILES "${PROJECT_BINARY_DIR}/python/core/environment.py"
      DESTINATION python/El/core)
  endif()
endif()
if(NOT NO_BINARY_SUBDIRECTORIES)
   set(EXPERIMENTAL_INSTALL_DIR experimental/g3d)
endif()
# Experimental
# ============
# NOTE: This is insanely outdated
if(EL_EXPERIMENTAL)
  set(EXPERIMENTAL_DIR "${PROJECT_SOURCE_DIR}/experimental")
  # Build the G3D example(s)
  set(G3D_EXPERS G3DGemm)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/experimental/g3d")
  foreach(EXPER ${G3D_EXPERS})
    set(DRIVER "${EXPERIMENTAL_DIR}/g3d/${EXPER}.cpp")
    add_executable(experimental-g3d-${EXPER} "${DRIVER}")
    set_source_files_properties("${DRIVER}" PROPERTIES
      OBJECT_DEPENDS "${PREPARED_HEADERS}")
    target_link_libraries(experimental-g3d-${EXPER} El)
    set_target_properties(experimental-g3d-${EXPER} PROPERTIES
      OUTPUT_NAME ${EXPER}
      SUFFIX "${CMAKE_EXECUTABLE_SUFFIX_CXX}"
      RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(EL_LINK_FLAGS)
      set_target_properties(experimental-g3d-${EXPER} PROPERTIES
        LINK_FLAGS ${EL_LINK_FLAGS})
    endif()
    install(TARGETS experimental-g3d-${EXPER} DESTINATION ${CMAKE_INSTALL_BINDIR}/${EXPERIMENTAL_INSTALL_DIR})
  endforeach()
endif()

# Install data
# ============
# TODO(poulson): Also install into the build directory? How to have drivers
# automatically find the right path?
install(DIRECTORY "${PROJECT_SOURCE_DIR}/data/"
	DESTINATION ${CMAKE_INSTALL_DATADIR}/elemental/ FILES_MATCHING PATTERN "*.*")

# Install license and author information
# ======================================
install(FILES "${PROJECT_SOURCE_DIR}/AUTHORS"
              "${PROJECT_SOURCE_DIR}/LICENSE"
              "${PROJECT_SOURCE_DIR}/PUBLICATIONS.bib"
              "${PROJECT_SOURCE_DIR}/README.md"
              "${PROJECT_SOURCE_DIR}/REFERENCES.bib"
	      DESTINATION ${CMAKE_INSTALL_DOCDIR})

# Tests and examples
# ==================
include(CTest)

# Sandbox test
# ------------
install(FILES "${PROJECT_SOURCE_DIR}/sandbox/Makefile" DESTINATION ${CMAKE_INSTALL_DOCDIR}/sandbox)
install(FILES "${PROJECT_SOURCE_DIR}/sandbox/test.cpp" DESTINATION ${CMAKE_INSTALL_DOCDIR}/sandbox)
add_executable(sandbox-test "${PROJECT_SOURCE_DIR}/sandbox/test.cpp")
set_source_files_properties("${PROJECT_SOURCE_DIR}/sandbox/test.cpp"
  PROPERTIES OBJECT_DEPENDS "${PREPARED_HEADERS}")
target_link_libraries(sandbox-test El)
set_target_properties(sandbox-test PROPERTIES
  OUTPUT_NAME sandbox-test
  SUFFIX "${CMAKE_EXECUTABLE_SUFFIX_CXX}"
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
if(EL_LINK_FLAGS)
  set_target_properties(sandbox-test PROPERTIES LINK_FLAGS ${EL_LINK_FLAGS})
endif()
add_test(NAME sandbox/test WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/bin"
  COMMAND sandbox-test)

# Test drivers
# ------------
if(EL_TESTS)
  set(TEST_DIR "${PROJECT_SOURCE_DIR}/tests")
  set(TEST_TYPES core blas_like lapack_like optimization)
  foreach(TYPE ${TEST_TYPES})
    file(GLOB_RECURSE ${TYPE}_TESTS
      RELATIVE "${PROJECT_SOURCE_DIR}/tests/${TYPE}/" "tests/${TYPE}/*.cpp")

    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/${TYPE}")
    foreach(TEST ${${TYPE}_TESTS})
      set(DRIVER "${TEST_DIR}/${TYPE}/${TEST}")
      get_filename_component(TESTNAME ${TEST} NAME_WE)
      add_executable(tests-${TYPE}-${TESTNAME} "${DRIVER}")
      set_source_files_properties("${DRIVER}" PROPERTIES
        OBJECT_DEPENDS "${PREPARED_HEADERS}")
      target_link_libraries(tests-${TYPE}-${TESTNAME} El)
      set_target_properties(tests-${TYPE}-${TESTNAME} PROPERTIES
        OUTPUT_NAME ${TESTNAME}
        SUFFIX "${CMAKE_EXECUTABLE_SUFFIX_CXX}"
        RUNTIME_OUTPUT_DIRECTORY "${OUTPUT_DIR}")
      if(EL_LINK_FLAGS)
        set_target_properties(tests-${TYPE}-${TESTNAME} PROPERTIES
          LINK_FLAGS ${EL_LINK_FLAGS})
      endif()
      if(NOT NO_BINARY_SUBDIRECTORIES)
      	set(TEST_INSTALL_DIR test/${TYPE})
      endif()
      install(TARGETS tests-${TYPE}-${TESTNAME} DESTINATION ${CMAKE_INSTALL_BINDIR}/${TEST_INSTALL_DIR})
      if(NOT TESTNAME STREQUAL "SparseLDLRange") #Skip tests that can time out
        add_test(NAME Tests/${TYPE}/${TESTNAME} 
          WORKING_DIRECTORY "${TEST_DIR}" COMMAND tests-${TYPE}-${TESTNAME})
      endif()
    endforeach()
  endforeach()
endif()

# Examples
# --------
if(EL_EXAMPLES)
  set(EXAMPLE_DIR "${PROJECT_SOURCE_DIR}/examples")
  set(EXAMPLE_TYPES
      optimization core blas_like interface io lapack_like matrices
      number_theory)
  foreach(TYPE ${EXAMPLE_TYPES})
    if(EL_C_INTERFACE)
      file(GLOB_RECURSE ${TYPE}_EXAMPLES RELATIVE
           "${PROJECT_SOURCE_DIR}/examples/${TYPE}/"
           "examples/${TYPE}/*.cpp"
           "examples/${TYPE}/*.c" )
    else()
      file(GLOB_RECURSE ${TYPE}_EXAMPLES RELATIVE
           "${PROJECT_SOURCE_DIR}/examples/${TYPE}/"
           "examples/${TYPE}/*.cpp")
    endif()
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples/${TYPE}")
    foreach(EXAMPLE ${${TYPE}_EXAMPLES})
      set(DRIVER "${EXAMPLE_DIR}/${TYPE}/${EXAMPLE}")
      get_filename_component(EXNAME ${EXAMPLE} NAME_WE)
      add_executable(examples-${TYPE}-${EXNAME} "${DRIVER}")
      set_source_files_properties("${DRIVER}" PROPERTIES
        OBJECT_DEPENDS "${PREPARED_HEADERS}")
      target_link_libraries(examples-${TYPE}-${EXNAME} El)
      set_target_properties(examples-${TYPE}-${EXNAME} PROPERTIES
        OUTPUT_NAME ${EXNAME}
        SUFFIX "${CMAKE_EXECUTABLE_SUFFIX_CXX}"
        RUNTIME_OUTPUT_DIRECTORY "${OUTPUT_DIR}")
      if(EL_LINK_FLAGS)
        set_target_properties(examples-${TYPE}-${EXNAME} PROPERTIES
          LINK_FLAGS ${EL_LINK_FLAGS})
      endif()
      if(NOT NO_BINARY_SUBDIRECTORIES)
      	set(EXAMPLE_INSTALL_DIR examples/${TYPE})
      endif()
      install(TARGETS examples-${TYPE}-${EXNAME}
	      DESTINATION ${CMAKE_INSTALL_BINDIR}/${EXAMPLE_INSTALL_DIR})
      add_test(NAME Examples/${TYPE}/${EXNAME} WORKING_DIRECTORY "${OUTPUT_DIR}"
        COMMAND examples-${TYPE}-${EXNAME})
    endforeach()
  endforeach()
endif()

# CPack
# =====
# While PackageMaker is deprecated, productbuild is not yet supported by CPack.
# Furthermore, it is possible to download the "Late July 2012" version of
# PackageMaker (which curiously has a release date of August 7, 2012) from 
# Apple, but it appears to be incompatible with Yosemite.
#
# Lastly, the current CPack configuration does not yet allow the Python
# interface to be installed within a different directory than the main library
# (e.g., into the user-level or system-wide site-packages directory).
set(CPACK_PACKAGE_NAME "Elemental")
set(CPACK_PACKAGE_VENDOR "libelemental.org")
set(CPACK_PACKAGE_VERSION_MAJOR ${EL_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${EL_VERSION_MINOR})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY 
  "Distributed-memory linear algebra and optimization")
set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/cmake/El.png")
include(CPack)

# ProjectConfig
# =====
# Create ElementalConfig.cmake for find_package.
# Add all targets to the build-tree export set
export(TARGETS El pmrrr ElSuiteSparse
  FILE "${PROJECT_BINARY_DIR}/ElementalTargets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE Elemental)

# Create the ElementalConfig.cmake and ElementalConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${PROJECT_BINARY_DIR}"
  "${PROJECT_BINARY_DIR}/include")
# ... for the build tree
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/configure_files/ElementalConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/ElementalConfig.cmake" @ONLY)
# ... for both
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/configure_files/ElementalConfigVersion.cmake.in"
  "${PROJECT_BINARY_DIR}/ElementalConfigVersion.cmake" @ONLY)

# Install the ElementalConfig.cmake and ElementalConfigVersion.cmake
set(INSTALL_CMAKE_DIR CMake)
install(FILES
  "${PROJECT_BINARY_DIR}/ElementalConfig.cmake"
  "${PROJECT_BINARY_DIR}/ElementalConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_SYSCONFDIR}/elemental/${INSTALL_CMAKE_DIR}" COMPONENT dev)

install(EXPORT ElementalTargets
	DESTINATION "${CMAKE_INSTALL_SYSCONFDIR}/elemental/${INSTALL_CMAKE_DIR}" COMPONENT dev)
