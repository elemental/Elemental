#
#  Copyright (c) 2009-2013, Jack Poulson
#                     2012, Jed Brown 
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
cmake_minimum_required(VERSION 2.8.8)
include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)
include(CheckFunctionExists)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
include(language_support_v2)

project(Elemental)
set(Elemental_VERSION_MAJOR 0)
set(Elemental_VERSION_MINOR 83-dev)

################################################################################
# Common build options                                                         #
################################################################################

# By default, Elemental builds static libraries, but shared libraries are 
# also supported through this option. If the SWIG interface is to be enabled,
# then Elemental switches to dynamic libraries in order to avoid massive 
# duplication.
option(SHARED_LIBRARIES "Switch to shared libraries" OFF)

# Whether or not to build a collection of simple example drivers.
option(ELEM_EXAMPLES "A collection of simple examples" OFF)

# Whether or not to build a collection of performance and correctness tests
option(ELEM_TESTS "A collection of performance and correctness tests" OFF)

# Whether or not to build current 'experimental' code, such as:
# 1) prototype C and F90 interfaces
# 2) Martin's 3D Gemm code
option(ELEM_EXPERIMENTAL "Build experimental code" OFF)

# Attempt to use 64-bit integers
option(USE_64BIT_INTS "Use 64-bit integers where possible" OFF)

# Petschow and Bientinesi's Parallel Multiple Relatively Robust 
# Representations (PMRRR) is used at the heart of all of Elemental's current 
# eigensolvers. Since it makes use of pthreads, Windows builds cannot use 
# PMRRR.
if(WIN32)
  option(USE_PMRRR "Build and link against eigensolver PMRRR" OFF)
else()
  option(USE_PMRRR "Build and link against eigensolver PMRRR" ON)
endif()

# Whether or not to use Qt5 for visualization if it is found (experimental)
# NOTE: Setting Qt5Widgets_DIR, Qt5_DIR, or Qt5_LIBDIR overrides this variable
option(USE_QT5 "Attempt to use Qt5?" OFF)

# Whether or not to attempt to build a Python interface via SWIG (experimental)
option(USE_SWIG "Attempt to build Python wrappers?" OFF)

################################################################################
# Advanced build options                                                       #
################################################################################

# Since it is surprisingly common for MPI libraries to have bugs in their 
# support for complex data, the following option forces Elemental to cast 
# all possible MPI communications in terms of twice as many real units of data.
option(AVOID_COMPLEX_MPI "Avoid potentially buggy complex MPI routines" OFF)
mark_as_advanced(AVOID_COMPLEX_MPI)

# At one point, a bug was found in IBM's C++ compiler for Blue Gene/P, 
# where OpenMP statements of the form a[i] += alpha b[i], with complex data,
# would segfault and/or return incorrect results
option(AVOID_OMP_FMA "Avoid a bug in the IBM compilers." OFF)
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

# "CXX_FLAGS"
#   Optimization and debugging/symbol flags

# "OpenMP_CXX_FLAGS"
#   Overrides the default compile flags for adding OpenMP support to CXX code

# TODO: Add list of MPI variables here, such as "MPI_CXX_COMPILE_FLAGS"

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
endif()

# Get the Git revision
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Extract a few booleans from the build type and default to PureRelease
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Build mode not specified, defaulting to PureRelease build.")
  set(CMAKE_BUILD_TYPE PureRelease)
endif()
if(CMAKE_BUILD_TYPE STREQUAL "HybridDebug")
  set(HYBRID TRUE)
  set(DEBUG TRUE)
elseif(CMAKE_BUILD_TYPE STREQUAL "HybridRelease")
  set(HYBRID TRUE)
  set(RELEASE TRUE)
elseif(CMAKE_BUILD_TYPE STREQUAL "PureDebug")
  set(PURE TRUE)
  set(DEBUG TRUE)
elseif(CMAKE_BUILD_TYPE STREQUAL "PureRelease")
  set(PURE TRUE)
  set(RELEASE TRUE)
else()
  message(FATAL_ERROR "CMAKE_BUILD_TYPE must be [Hybrid,Pure][Debug,Release]")
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} UPPER_BUILD_TYPE)

# Set the basic compile flags from the build type
if(NOT WIN32) # I forget why we need this guard
  set(LANGUAGES CXX C Fortran)
  foreach(LANG ${LANGUAGES})
    if(NOT ${LANG}_FLAGS_${UPPER_BUILD_TYPE})
      if(DEBUG)
        set(${LANG}_FLAGS_${UPPER_BUILD_TYPE} "-O2 -g" CACHE STRING
          "${LANG} optimization/debug flags for ${UPPER_BUILD_TYPE} mode")
      else()
        set(${LANG}_FLAGS_${UPPER_BUILD_TYPE} "-O3" CACHE STRING
          "${LANG} optimization flags for ${UPPER_BUILD_TYPE} mode")
      endif()
    endif()
  endforeach()
endif()

# Go ahead and check for Fortran, but keep in mind that CMake's 'OPTIONAL' 
# argument for enable_language is still completely broken as of 2.8.8
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
    endif()
  endif()
endif()

#
# Check C++ support
#
include(FindCXXFeatures)
message(STATUS "CXX11_COMPILER_FLAGS=${CXX11_COMPILER_FLAGS}")
set(EXTRA_FLAGS "${CXX11_COMPILER_FLAGS}")
set(CMAKE_REQUIRED_FLAGS ${CXX11_COMPILER_FLAGS})
set(RANDOM_CODE
    "#include <random>
     int main()
     {
         std::random_device rd;
         std::mt19937 mt(rd());
         std::normal_distribution<double> dist(0,1);
         const double x = dist(mt);
     }")
set(ALIAS_CODE
    "#include <complex>
     template<typename Real>
     using Complex = std::complex<Real>;
     int main()
     {
         Complex<double> a;
     }")
check_cxx_source_compiles("${RANDOM_CODE}" HAVE_NORMAL_DIST)
check_cxx_source_compiles("${ALIAS_CODE}" HAVE_TEMPLATE_ALIAS)
if(NOT HAVE_NORMAL_DIST OR NOT HAVE_TEMPLATE_ALIAS)
  message(FATAL_ERROR "C++11 random number generation and template aliasing not found. You may want to make sure that your compiler is sufficiently up-to-date (e.g., g++ >= 4.7)")
endif()
# Check for restrict
set(RESTRICT_CODE "int main() { int* RESTRICT a; return 0; }")
set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict__")
check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict__)
set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict")
check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict)
set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=restrict")
check_cxx_source_compiles("${RESTRICT_CODE}" HAVE_restrict)
if(HAVE___restrict__)
  set(RESTRICT "__restrict__")
  message(STATUS "Using __restrict__ keyword.")
elseif(HAVE___restrict)
  set(RESTRICT "__restrict")
  message(STATUS "Using __restrict keyword.")
elseif(HAVE_restrict)
  set(RESTRICT "restrict")
  message(STATUS "Using restrict keyword.")
else()
  set(RESTRICT "")
  message(STATUS "Could not find a restrict keyword.")
endif()

#
# Check for MPI support
#
find_package(MPI)
if(NOT MPI_CXX_FOUND)
  message(FATAL_ERROR "MPI C++ compiler was not found and is required")
endif()
include_directories(${MPI_CXX_INCLUDE_PATH})
set(EXTRA_FLAGS "${EXTRA_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
# Look for MPI_Reduce_scatter_block (and MPI_Reduce_scatter as sanity check)
set(CMAKE_REQUIRED_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_C_LINK_FLAGS}")
set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
check_function_exists(MPI_Reduce_scatter HAVE_MPI_REDUCE_SCATTER)
if(NOT HAVE_MPI_REDUCE_SCATTER)
  message(FATAL_ERROR "Could not find MPI_Reduce_scatter")
endif()
check_function_exists(MPI_Reduce_scatter_block HAVE_MPI_REDUCE_SCATTER_BLOCK)
check_function_exists(MPI_Iallgather HAVE_MPI3_NONBLOCKING_COLLECTIVES)
check_function_exists(MPIX_Iallgather HAVE_MPIX_NONBLOCKING_COLLECTIVES)
check_function_exists(MPI_Init_thread HAVE_MPI_INIT_THREAD)
check_function_exists(MPI_Query_thread HAVE_MPI_QUERY_THREAD)
check_function_exists(MPI_Comm_set_errhandler HAVE_MPI_COMM_SET_ERRHANDLER)
check_function_exists(MPI_Type_create_struct HAVE_MPI_TYPE_CREATE_STRUCT)
if(NOT HAVE_MPI_TYPE_CREATE_STRUCT)
  message(FATAL_ERROR "Could not find MPI_Type_create_struct")
endif()
# Check for MPI_IN_PLACE (essentially MPI2 support)
set(MPI_IN_PLACE_CODE
    "#include \"mpi.h\"
     int main( int argc, char* argv[] ) 
     {
         MPI_Init( &argc, &argv );
         float a;
         MPI_Allreduce
         ( MPI_IN_PLACE, &a, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
         MPI_Finalize(); 
         return 0; 
     }")
set(MPI_LONG_LONG_CODE
    "#include \"mpi.h\"
     int main( int argc, char* argv[] )
     {
         MPI_Init( &argc, &argv );
         MPI_Datatype lli = MPI_LONG_LONG_INT;
         MPI_Datatype llu = MPI_UNSIGNED_LONG_LONG;
         MPI_Finalize();
         return 0;
     }")
check_c_source_compiles("${MPI_IN_PLACE_CODE}" HAVE_MPI_IN_PLACE)
check_c_source_compiles("${MPI_LONG_LONG_CODE}" HAVE_MPI_LONG_LONG)
if(USE_64BIT_INTS AND NOT HAVE_MPI_LONG_LONG)
  message(FATAL_ERROR 
    "Did not detect MPI_LONG_LONG_INT and MPI_UNSIGNED_LONG_LONG")
endif()

#
# Check for OpenMP support if necessary
#
if(HYBRID)
  if(OpenMP_CXX_FLAGS)
    set(HAVE_OPENMP TRUE)
    set(OpenMP_CXX_FLAGS ${OpenMP_CXX_FLAGS})
    message(STATUS "Using OpenMP_CXX_FLAGS=${OpenMP_CXX_FLAGS}")
  else()
    find_package(OpenMP)
    if(OPENMP_FOUND)
      set(HAVE_OPENMP TRUE)
    else()
      set(OpenMP_CXX_FLAGS "" CACHE STRING "OpenMP CXX FLAGS")
      message(FATAL_ERROR
        "Hybrid build failed because OpenMP support not detected. Please specify OpenMP_CXX_FLAGS.")
    endif()
  endif()
  set(EXTRA_FLAGS "${OpenMP_CXX_FLAGS} ${EXTRA_FLAGS}")
endif()
# See if we have 'collapse' support, since Clang is currently setting _OPENMP
# greater than 200805 despite not supporting OpenMP 3.0 (i.e., collapse)
set(CMAKE_REQUIRED_FLAGS ${OpenMP_CXX_FLAGS})
set(OMP_COLLAPSE_CODE
    "#include <omp.h>
     int main( int argc, char* argv[] ) 
     {
         int k[100];
     #pragma omp collapse(2)
         for( int i=0; i<10; ++i ) 
             for( int j=0; j<10; ++j )
                 k[i+j*10] = i+j; 
         return 0; 
     }")
check_cxx_source_compiles("${OMP_COLLAPSE_CODE}" HAVE_OMP_COLLAPSE)

#
# Check for Qt5 support
#
set(HAVE_QT5 FALSE)
if(Qt5_DIR)
  set(Qt5_LIBDIR "${Qt5_DIR}/lib")
endif()
if(Qt5_LIBDIR)
  set(Qt5Widgets_DIR "${Qt5_LIBDIR}/cmake/Qt5Widgets")
endif()
if(USE_QT5 OR Qt5Widgets_DIR)
  # Search for Qt5
  find_package(Qt5Widgets)
  if(Qt5Widgets_FOUND)
    set(ELEM_HEADERS_PREMOC 
        "include/elemental/io/display_window_decl.hpp;include/elemental/io/complex_display_window_decl.hpp")
    qt_wrap_cpp(elemental ELEM_MOC_SRC ${ELEM_HEADERS_PREMOC})
    include_directories(${Qt5Widgets_INCLUDE_DIRS})
    add_definitions(${Qt5Widgets_DEFINITIONS})
    set(EXTRA_FLAGS "${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS} ${EXTRA_FLAGS}")

    # Qt5Widgets_DIR = Qt5_LIBDIR/cmake/Qt5Widgets
    get_filename_component(Qt5_CMAKEDIR ${Qt5Widgets_DIR} PATH)
    get_filename_component(Qt5_LIBDIR ${Qt5_CMAKEDIR} PATH)

    set(HAVE_QT5 TRUE)
    message(STATUS "Found Qt5")
  else()
    message(STATUS "Did NOT find Qt5")
  endif()
endif()

#
# Check for SWIG support
#
set(HAVE_SWIG FALSE)
if(USE_SWIG)
  # Search for SWIG, Python, and NumPy
  find_package(SWIG)
  find_package(PythonLibs)
  set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
  find_package(NumPy)
  if(SWIG_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)
    if(${NUMPY_VERSION_MAJOR} GREATER 0 AND ${NUMPY_VERSION_MINOR} GREATER 6)
      include_directories(${NUMPY_INCLUDE_DIRS})
      set(HAVE_SWIG TRUE)
    else()
      message(STATUS "NumPy version must be >= 1.7")
    endif()
  else()
    message(STATUS "Did not find necessary SWIG, Python, and NumPy libs")
  endif()
endif()

# Set the library type from the 'SHARED_LIBRARIES' option
if(SHARED_LIBRARIES OR HAVE_SWIG)
  if(NOT SHARED_LIBRARIES)
    message(WARNING "Switched to a dynamic build since SWIG is enabled")
  endif()
  set(LIBRARY_TYPE SHARED)
else()
  set(LIBRARY_TYPE STATIC)
endif()

#
# Check for BLAS and LAPACK support
#
if(PURE)
  set(MATH_DESC "Unthreaded BLAS/LAPACK link flags")
else()
  set(MATH_DESC "Threaded BLAS/LAPACK link flags")
endif()
if(MATH_LIBS)
  message(STATUS "Using user-defined MATH_LIBS=${MATH_LIBS}")
elseif(APPLE)
  set(MATH_LIBS "-framework vecLib" CACHE STRING ${MATH_DESC})
  message(STATUS "Using Apple vecLib framework.")
else()
  # Look for default BLAS and LAPACK
  if(REFERENCE_ROOT)
    message(STATUS "Searching REFERENCE_ROOT=${REFERENCE_ROOT} for math libs")
  endif()
  set(REFERENCE_REQUIRED LAPACK BLAS)
  find_library(BLAS_LIB NAMES blas PATHS ${REFERENCE_ROOT})
  find_library(LAPACK_LIB NAMES lapack reflapack PATHS ${REFERENCE_ROOT}) 
  set(REFERENCE_FOUND TRUE)
  foreach(NAME ${REFERENCE_REQUIRED})
    if(${NAME}_LIB)
      message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
      list(APPEND MATH_LIBS ${${NAME}_LIB})
    else()
      message(STATUS "Could not find ${NAME}_LIB")
      set(MATH_LIBS "")
      set(REFERENCE_FOUND FALSE)
    endif()
  endforeach()
  if(REFERENCE_FOUND)
    message(WARNING "Using reference BLAS/LAPACK; performance will be poor")
  else()
    message(FATAL_ERROR "Could not find BLAS/LAPACK. Please specify MATH_LIBS")
  endif()
endif()
# Check the BLAS and LAPACK underscore conventions
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
check_function_exists(daxpy  HAVE_DAXPY)
check_function_exists(daxpy_ HAVE_DAXPY_POST)
check_function_exists(dpotrf  HAVE_DPOTRF)
check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
if(HAVE_DAXPY)
  set(BLAS_POST FALSE)
  set(BLAS_DEFS "")
elseif(HAVE_DAXPY_POST)
  set(BLAS_POST TRUE)
  set(BLAS_DEFS "-DBLAS_POST")
else()
  message(FATAL_ERROR "Could not determine BLAS format.")
endif()
if(HAVE_DPOTRF)
  set(LAPACK_POST FALSE)
  set(LAPACK_DEFS "")
elseif(HAVE_DPOTRF_POST)
  set(LAPACK_POST TRUE)
  set(LAPACK_DEFS "-DLAPACK_POST")
else()
  message(FATAL_ERROR "Could not determine LAPACK format.")
endif()
# Ensure that we have a relatively new version of LAPACK
if(HAVE_DPOTRF)
  check_function_exists(dsyevr HAVE_DSYEVR)
  if(NOT HAVE_DSYEVR)
    message(FATAL_ERROR "LAPACK is missing dsyevr")
  endif()
else()
  check_function_exists(dsyevr_ HAVE_DSYEVR_POST)
  if(NOT HAVE_DSYEVR_POST)
    message(FATAL_ERROR "LAPACK is missing dsyevr_")
  endif()
endif()

#
# Check for libFLAME support
#
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
check_function_exists(FLA_Bsvd_v_opd_var1 HAVE_FLA_BSVD)

# Now append the extra flags
set(CMAKE_CXX_FLAGS_${UPPER_BUILD_TYPE} "${CXX_FLAGS} ${EXTRA_FLAGS}")

# Add the Parallel Multiple Relatively Robust Representations (PMRRR) project
# if necessary
set(HAVE_PMRRR FALSE)
if(USE_PMRRR)
  add_subdirectory(external/pmrrr)
endif()

# Create the Elemental configuration header
configure_file(${PROJECT_SOURCE_DIR}/cmake/config.h.cmake
               ${PROJECT_BINARY_DIR}/include/elemental/config.h)
install(FILES ${PROJECT_BINARY_DIR}/include/elemental/config.h 
        DESTINATION include/elemental)

# Create a file which can be included in Makefile's.
# This is meant to be analogous to PETSc's 'conf/petscvariables' file
set(MPI_CXX_INCSTRING)
foreach(INC_PATH ${MPI_CXX_INCLUDE_PATH})
  set(MPI_CXX_INCSTRING "${MPI_CXX_INCSTRING} -I${INC_PATH}")  
endforeach()
set(MATH_LIBSTRING)
foreach(LIB ${MATH_LIBS})
  set(MATH_LIBSTRING "${MATH_LIBSTRING} ${LIB}")
endforeach()
set(MPI_CXX_LIBSTRING)
foreach(LIB ${MPI_CXX_LIBRARIES})
  set(MPI_CXX_LIBSTRING "${MPI_CXX_LIBSTRING} ${LIB}")
endforeach()
set(LAPACK_ADDONSTRING)
foreach(LIB ${LAPACK_ADDONS_LIBS})
  set(LAPACK_ADDONSTRING "${LAPACK_ADDONSTRING} ${LIB}")
endforeach()
set(QT5_DEFSTRING)
foreach(DEF ${Qt5Widgets_DEFINITIONS})
  set(QT5_DEFSTRING "${QT5_DEFSTRING} ${DEF}")
endforeach()
set(QT5_INCSTRING)
foreach(INC ${Qt5Widgets_INCLUDE_DIRS})
  set(QT5_INCSTRING "${QT5_INCSTRING} -I${INC}")
endforeach()
set(QT5_COMPILESTRING)
foreach(FLAG ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})
  set(QT5_COMPILESTRING "${QT5_COMPILE_STRING} ${FLAG}")
endforeach()
set(QT5_LIBSTRING "-L${Qt5_LIBDIR} -lQt5Widgets -lQt5Gui -lQt5Core")
configure_file(${PROJECT_SOURCE_DIR}/cmake/elemvariables.cmake
               ${PROJECT_BINARY_DIR}/conf/elemvariables @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/conf/elemvariables DESTINATION conf)

# Grab all of the .c, .cpp, .h, and .hpp Elemental files
file(GLOB_RECURSE ELEM_CPP RELATIVE ${PROJECT_SOURCE_DIR} "src/*.c" "src/*.cpp")
file(GLOB_RECURSE ELEM_HEADERS RELATIVE ${PROJECT_SOURCE_DIR} 
  "include/*.h" "include/*.hpp")
if(ELEM_EXPERIMENTAL AND HAVE_PMRRR)
  if(NOT MPI_C_FOUND)
    message(FATAL_ERROR "Could not build C interface without C MPI")
  endif()
  include_directories(${MPI_C_INCLUDE_PATH})
  set(CMAKE_C_FLAGS_${UPPER_BUILD_TYPE} "${MPI_C_COMPILE_FLAGS}")
  set(ELEM_CPP "${ELEM_CPP};experimental/c/interface.cpp")
  set(ELEM_HEADERS "${ELEM_HEADERS};experimental/c/elemental.h")

  if(HAVE_F90_INTERFACE)
    if(NOT MPI_Fortran_FOUND)
      message(FATAL_ERROR "Could not build F90 interface without Fortran MPI")
    endif()
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
     }")
    check_c_source_compiles("${MPI_COMM_F2C_CODE}" HAVE_MPI_COMM_F2C)
    if(NOT HAVE_MPI_COMM_F2C)
      message(FATAL_ERROR
        "Could not find MPI_Comm_f2c, but it is needed for the F90 interface")
    endif()
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    set(CMAKE_Fortran_FLAGS_${UPPER_BUILD_TYPE} "${MPI_Fortran_COMPILE_FLAGS}")
    set(ELEM_CPP "${ELEM_CPP};experimental/f90/interface.cpp")
  endif()
endif()
set(ELEM_SRC "${ELEM_CPP};${ELEM_HEADERS};${ELEM_MOC_SRC}")

# The main library
set(LINK_LIBS ${MATH_LIBS} ${MPI_CXX_LIBRARIES})
if(HAVE_QT5)
  set(LINK_LIBS ${LINK_LIBS} ${Qt5Widgets_LIBRARIES})  
endif()
if(HAVE_PMRRR)
  add_library(elemental ${LIBRARY_TYPE} ${ELEM_SRC})
  target_link_libraries(elemental pmrrr ${LINK_LIBS})
  if(MPI_CXX_LINK_FLAGS AND MPI_C_LINK_FLAGS)
    set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS} ${MPI_C_LINK_FLAGS}")
  elseif(MPI_CXX_LINK_FLAGS)
    set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  elseif(MPI_C_LINK_FLAGS)
    set(MPI_LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()
else()
  add_library(elemental ${LIBRARY_TYPE} ${ELEM_SRC})
  target_link_libraries(elemental ${LINK_LIBS})
  if(MPI_CXX_LINK_FLAGS)
    set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif()
endif()
install(TARGETS elemental DESTINATION lib)

# Define the header-file preparation rules
set(PREPARED_HEADERS)
foreach(HEADER ${ELEM_HEADERS})
  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/${HEADER}
            ${PROJECT_BINARY_DIR}/${HEADER}
    DEPENDS "${PROJECT_SOURCE_DIR}/${HEADER}")
  list(APPEND PREPARED_HEADERS ${PROJECT_BINARY_DIR}/${HEADER})

  get_filename_component(HEADER_PATH ${HEADER} PATH)
  install(FILES ${PROJECT_BINARY_DIR}/${HEADER} DESTINATION ${HEADER_PATH})
endforeach()
add_custom_target(prepare_elemental_headers DEPENDS ${PREPARED_HEADERS})
add_dependencies(elemental prepare_elemental_headers)

# Make sure the Elemental headers can be found
include_directories("${PROJECT_BINARY_DIR}/include")

# Build the SWIG wrappers if necessary
if(HAVE_SWIG)
  include(${SWIG_USE_FILE})
  include_directories(${PYTHON_INCLUDE_PATH})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  set(CMAKE_SWIG_FLAGS "-fvirtual")
  # NOTE: There seems to be a bug in the SWIG module where the necessary
  #       subdirectories (i.e., swig/) are not created automatically
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/swig)
  set_source_files_properties(swig/elem.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(swig/mpi.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(swig/view.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(swig/blas.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(swig/lapack.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(swig/control.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(swig/convex.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(swig/matrices.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(swig/io.i PROPERTIES CPLUSPLUS ON)
  swig_add_module(elem python swig/elem.i)
  swig_add_module(elem_mpi python swig/mpi.i)
  swig_add_module(elem_view python swig/view.i)
  swig_add_module(elem_blas python swig/blas.i)
  swig_add_module(elem_lapack python swig/lapack.i)
  swig_add_module(elem_control python swig/control.i)
  swig_add_module(elem_convex python swig/convex.i)
  swig_add_module(elem_matrices python swig/matrices.i)
  swig_add_module(elem_io python swig/io.i)
  swig_link_libraries(elem ${PYTHON_LIBRARIES} elemental)
  swig_link_libraries(elem_mpi ${PYTHON_LIBRARIES} elemental)
  swig_link_libraries(elem_view ${PYTHON_LIBRARIES} elemental)
  swig_link_libraries(elem_blas ${PYTHON_LIBRARIES} elemental)
  swig_link_libraries(elem_lapack ${PYTHON_LIBRARIES} elemental)
  swig_link_libraries(elem_control ${PYTHON_LIBRARIES} elemental)
  swig_link_libraries(elem_convex ${PYTHON_LIBRARIES} elemental)
  swig_link_libraries(elem_matrices ${PYTHON_LIBRARIES} elemental)
  swig_link_libraries(elem_io ${PYTHON_LIBRARIES} elemental)
endif()

# Build the test drivers if necessary
if(ELEM_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/tests)
  set(TEST_TYPES core blas-like lapack-like convex)

  set(core_TESTS AxpyInterface DifferentGrids DistMatrix Matrix Version)
  set(blas-like_TESTS 
    Gemm Hemm Her2k Herk Symm Symv Syr2k Syrk Trmm Trsm Trsv TwoSidedTrmm
    TwoSidedTrsm)
  set(lapack-like_TESTS 
    ApplyPackedReflectors Cholesky CholeskyQR HermitianTridiag LDL LU LQ QR RQ
    SequentialLU TriangularInverse TSQR)
  if(HAVE_PMRRR)
    list(APPEND lapack-like_TESTS HermitianEig HermitianGenDefiniteEig)
  endif()
  set(convex_TESTS TSSVT)

  foreach(TYPE ${TEST_TYPES})
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/${TYPE}")
    foreach(TEST ${${TYPE}_TESTS})
      add_executable(tests-${TYPE}-${TEST} ${TEST_DIR}/${TYPE}/${TEST}.cpp)
      target_link_libraries(tests-${TYPE}-${TEST} elemental)
      set_target_properties(tests-${TYPE}-${TEST} PROPERTIES
        OUTPUT_NAME ${TEST} RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
      if(MPI_LINK_FLAGS)
        set_target_properties(tests-${TYPE}-${TEST} PROPERTIES
          LINK_FLAGS ${MPI_LINK_FLAGS})
      endif()
      install(TARGETS tests-${TYPE}-${TEST} DESTINATION bin/tests/${TYPE})
    endforeach()
  endforeach()
endif()

# Build the example drivers if necessary
if(ELEM_EXAMPLES)
  set(EXAMPLE_DIR ${PROJECT_SOURCE_DIR}/examples)
  set(EXAMPLE_TYPES convex core blas-like lapack-like matrices)

  set(convex_EXAMPLES LogDetDivergence RPCA)
  set(core_EXAMPLES Constructors LogicError)
  set(blas-like_EXAMPLES Cannon Gemm Gemv)
  set(lapack-like_EXAMPLES 
    BunchKaufman BusingerGolub GaussianElimination HermitianQDWH HermitianSDC 
    HPDInverse ID KyFanAndSchatten LDL LDLInverse LeastSquares Polar 
    Pseudoinverse QDWH RealSchur Schur SequentialBunchKaufman SequentialQR 
    SequentialSVD Skeleton QR Sign SimpleSVD SVD)
  set(matrices_EXAMPLES 
    Cauchy CauchyLike Circulant Diagonal Egorov Fourier Hankel Helmholtz1D
    Helmholtz2D Helmholtz3D HermitianUniformSpectrum Hilbert Identity Kahan 
    Legendre LehmerParterRis NormalUniformSpectrum Ones OneTwoOne PSFW 
    RiemannRedhefferGCD Toeplitz Uniform Walsh Wilkinson Zeros)
  if(HAVE_PMRRR)
    list(APPEND lapack-like_EXAMPLES
      HermitianEig HermitianEigFromSequential HermitianSVD SkewHermitianEig
      ComplexHermitianFunction RealHermitianFunction HermitianPseudoinverse
      HPSDCholesky HPSDSquareRoot RealSymmetricFunction)
  endif()

  foreach(TYPE ${EXAMPLE_TYPES})
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples/${TYPE}")
    foreach(EXAMPLE ${${TYPE}_EXAMPLES})
      add_executable(examples-${TYPE}-${EXAMPLE} 
        ${EXAMPLE_DIR}/${TYPE}/${EXAMPLE}.cpp)
      target_link_libraries(examples-${TYPE}-${EXAMPLE} elemental)
      set_target_properties(examples-${TYPE}-${EXAMPLE} PROPERTIES 
        OUTPUT_NAME ${EXAMPLE} RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
      if(MPI_LINK_FLAGS)
        set_target_properties(examples-${TYPE}-${EXAMPLE} PROPERTIES
          LINK_FLAGS ${MPI_LINK_FLAGS})
      endif()
      install(TARGETS examples-${TYPE}-${EXAMPLE} 
        DESTINATION bin/examples/${TYPE})
    endforeach()
  endforeach()
endif()

# Build experimental drivers that do not require an eigensolver
if(ELEM_EXPERIMENTAL)
  set(EXPERIMENTAL_DIR ${PROJECT_SOURCE_DIR}/experimental)
  set(G3D_EXPERS G3DGemm)

  # Build the G3D example(s)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/experimental/g3d")
  foreach(EXPER ${G3D_EXPERS})
    add_executable(experimental-g3d-${EXPER} 
      ${EXPERIMENTAL_DIR}/g3d/${EXPER}.cpp)
    target_link_libraries(experimental-g3d-${EXPER} elemental)
    set_target_properties(experimental-g3d-${EXPER} PROPERTIES 
      OUTPUT_NAME ${EXPER} RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(experimental-g3d-${EXPER} PROPERTIES
        LINK_FLAGS ${MPI_LINK_FLAGS})
    endif()
    install(TARGETS experimental-g3d-${EXPER} DESTINATION bin/experimental/g3d)
  endforeach()
endif()

# Build experimental drivers that DO require an eigensolver
if(ELEM_EXPERIMENTAL AND HAVE_PMRRR)
  set(EXPERIMENTAL_DIR ${PROJECT_SOURCE_DIR}/experimental)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/experimental/c")
  add_executable(experimental-c-test ${EXPERIMENTAL_DIR}/c/test.c)
  target_link_libraries(experimental-c-test elemental)
  set_target_properties(experimental-c-test 
    PROPERTIES OUTPUT_NAME c_test RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
  if(MPI_LINK_FLAGS)
    set_target_properties(experimental-c-test 
      PROPERTIES LINK_FLAGS ${MPI_LINK_FLAGS})
  endif()
  install(TARGETS experimental-c-test DESTINATION bin/experimental/c)  

  if(HAVE_F90_INTERFACE)
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/experimental/f90")
    add_executable(experimental-f90-test ${EXPERIMENTAL_DIR}/f90/test.f90) 
    target_link_libraries(experimental-f90-test elemental)
    set_target_properties(experimental-f90-test PROPERTIES
      LINKER_LANGUAGE Fortran
      OUTPUT_NAME f90_test 
      RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(experimental-f90-test 
        PROPERTIES LINK_FLAGS ${MPI_LINK_FLAGS})
    endif()
    install(TARGETS experimental-f90-test DESTINATION bin/experimental/f90)  
  endif()
endif()

# If Elemental is a subproject, then pass some variables to the parent
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(./cmake/elemental_sub.cmake)
endif()
