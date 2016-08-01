#
#  Copyright 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
include(ElCheckFunctionExists)
include(CheckCXXSourceCompiles)

# Default locations (currently, Linux-centric) for searching for math libs
# ========================================================================
# NOTE: The variables BLAS_ROOT, LAPACK_ROOT, MATH_ROOT, BLAS_PATH, LAPACK_PATH,
#       and MATH_PATH are all optional and to ease the process
set(MATH_PATHS /usr/lib
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
               ${BLAS_ROOT}/lib
               ${BLAS_ROOT}/lib64
               ${LAPACK_ROOT}/lib
               ${LAPACK_ROOT}/lib64
               ${MATH_ROOT}/lib
               ${MATH_ROOT}/lib64
               ${BLAS_PATH}
               ${LAPACK_PATH}
               ${MATH_PATH})

# Check for BLAS and LAPACK support
# =================================
# The default priorities of BLAS/LAPACK libraries are:
#   1. MKL
#   2. Accelerate/vecLib
#   3. OpenBLAS
#   4. BLIS
# Since MKL and Accelerate/vecLib cannot be built, their support is detected
# immediately before calling out to the nested CMake ScaLAPACK build.
# There is currently a bug when ScaLAPACK is built on Macs.

if(EL_BLAS_SUFFIX AND NOT EL_BLAS_SUFFIX STREQUAL _)
  set(CUSTOM_BLAS_SUFFIX TRUE)
endif()
if(EL_LAPACK_SUFFIX AND NOT EL_LAPACK_SUFFIX STREQUAL _)
  set(CUSTOM_LAPACK_SUFFIX TRUE)
endif()

# Test for pre-built libraries
# ----------------------------
if(MATH_LIBS)
  set(MATH_LIBS_AT_CONFIG ${MATH_LIBS})
  message(STATUS "Will attempt to extend user-defined MATH_LIBS=${MATH_LIBS}")
endif()

# Check for MKL
# ^^^^^^^^^^^^^
if(MATH_LIBS)
  set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
  El_check_function_exists(mkl_dcsrmv EL_HAVE_MKL_DCSRMV)
  if(EL_HAVE_MKL_DCSRMV)
    set(EL_HAVE_MKL TRUE)
    message(STATUS "Using Intel MKL via ${MATH_LIBS}")
    El_check_function_exists(dgemmt EL_HAVE_DGEMMT)
    El_check_function_exists(dgemmt_ EL_HAVE_DGEMMT_POST)
    if(EL_HAVE_DGEMMT OR EL_HAVE_DGEMMT_POST)
      # The potential underscore suffix will be taken care of elsewhere
      set(EL_HAVE_MKL_GEMMT TRUE)
    endif()
  endif()
  unset(CMAKE_REQUIRED_LIBRARIES)
# NOTE:
# There is a bug in MKL such that, if the following line is used,
# there is sometimes an error of the form
#
#   Intel MKL FATAL ERROR: Cannot load libmkl_avx2.so or libmkl_def.so.
#
# when running the executable. Should this occur, manually specify
#
#   MATH_LIBS="-L/path/to/mkl/libs -lmkl_rt"
#
# e.g.,
#
#   MATH_LIBS="-L/opt/intel/mkl/lib/intel64 -lmkl_rt"
#       
# Due to this error, as well as the fact that there is no simple support for
# the ILP64 interface (which is needed when EL_USE_64BIT_BLAS_INTS is enabled),
# the "-mkl=" detection is now disabled (by the 'AND FALSE' clause).
elseif(NOT EL_DISABLE_MKL AND
       NOT EL_USE_64BIT_BLAS_INTS AND
       NOT EL_PREFER_OPENBLAS AND
       NOT EL_PREFER_APPLE_MATH AND
       NOT EL_PREFER_BLIS_LAPACK AND
       FALSE)
  if(EL_HYBRID)
    set(MKL_LIBS "-mkl=parallel")
    message(STATUS "Attempting to link MKL using ${MKL_LIBS}")
    set(CMAKE_REQUIRED_FLAGS ${MKL_LIBS})
    El_check_function_exists(dpotrf  EL_HAVE_DPOTRF_MKL)
    El_check_function_exists(dpotrf_ EL_HAVE_DPOTRF_POST_MKL)
    El_check_function_exists(mkl_dcsrmv EL_HAVE_MKL_DCSRMV)
    if((EL_HAVE_DPOTRF_MKL OR EL_HAVE_DPOTRF_POST_MKL) AND EL_HAVE_MKL_DCSRMV)
      set(EL_FOUND_MKL TRUE)
      El_check_function_exists(dgemmt EL_HAVE_DGEMMT)
      El_check_function_exists(dgemmt_ EL_HAVE_DGEMMT_POST)
      if(EL_HAVE_DGEMMT OR EL_HAVE_DGEMMT_POST)
        # The potential underscore suffix will be taken care of elsewhere
        set(EL_HAVE_MKL_GEMMT TRUE)
      endif()
    endif()
    unset(CMAKE_REQUIRED_FLAGS)
  else()
    set(MKL_LIBS "-mkl=cluster")
    message(STATUS "Attempting to link MKL using ${MKL_LIBS}")
    set(CMAKE_REQUIRED_FLAGS "${MKL_LIBS} ${MPI_C_COMPILE_FLAGS}")
    set(CMAKE_REQUIRED_LINKER_FLAGS
      "${MPI_C_LINK_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
    set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
    El_check_function_exists(dpotrf  EL_HAVE_DPOTRF_MKL_CLUSTER)
    El_check_function_exists(dpotrf_ EL_HAVE_DPOTRF_POST_MKL_CLUSTER)
    El_check_function_exists(mkl_dcsrmv EL_HAVE_MKL_DCSRMV_CLUSTER)
    if((EL_HAVE_DPOTRF_MKL_CLUSTER OR EL_HAVE_DPOTRF_POST_MKL_CLUSTER) 
       AND EL_HAVE_MKL_DCSRMV_CLUSTER)
      set(EL_FOUND_MKL TRUE)
      El_check_function_exists(dgemmt EL_HAVE_DGEMMT)
      El_check_function_exists(dgemmt_ EL_HAVE_DGEMMT_POST)
      if(EL_HAVE_DGEMMT OR EL_HAVE_DGEMMT_POST)
        # The potential underscore suffix will be taken care of elsewhere
        set(EL_HAVE_MKL_GEMMT TRUE)
      endif()
    endif()
    unset(CMAKE_REQUIRED_FLAGS)
    unset(CMAKE_REQUIRED_LINKER_FLAGS)
    unset(CMAKE_REQUIRED_INCLUDES)
    unset(CMAKE_REQUIRED_LIBRARIES)

    if(NOT EL_FOUND_MKL)
      set(MKL_LIBS "-mkl=sequential")
      message(STATUS "Attempting to link MKL using ${MKL_LIBS}")
      set(CMAKE_REQUIRED_FLAGS ${MKL_LIBS})
      El_check_function_exists(dpotrf  EL_HAVE_DPOTRF_MKL_SEQ)
      El_check_function_exists(dpotrf_ EL_HAVE_DPOTRF_POST_MKL_SEQ)
      El_check_function_exists(mkl_dcsrmv EL_HAVE_MKL_DCSRMV_SEQ)
      if((EL_HAVE_DPOTRF_MKL_SEQ OR EL_HAVE_DPOTRF_POST_MKL_SEQ) 
         AND EL_HAVE_MKL_DCSRMV_SEQ)
        set(EL_FOUND_MKL TRUE)
        El_check_function_exists(dgemmt EL_HAVE_DGEMMT)
        El_check_function_exists(dgemmt_ EL_HAVE_DGEMMT_POST)
        if(EL_HAVE_DGEMMT OR EL_HAVE_DGEMMT_POST)
          # The potential underscore suffix will be taken care of elsewhere
          set(EL_HAVE_MKL_GEMMT TRUE)
        endif()
      endif()
      unset(CMAKE_REQUIRED_FLAGS)
    endif()
  endif()
  if(EL_FOUND_MKL)
    set(MATH_LIBS_AT_CONFIG ${MKL_LIBS}) 
    set(EL_HAVE_MKL TRUE)
    message(STATUS "Using Intel MKL via ${MKL_LIBS}")
    message(WARNING "If you receive an MKL Fatal Error when running your application, this is due to a bug in MKL that can be avoided by reconfiguring with MATH_LIBS containing an explicit link line for MKL, e.g., MATH_LIBS=\"-L/opt/mkl/lib/intel64 -lmkl_rt\"")
  endif()
endif()

if(APPLE AND NOT EL_USE_64BIT_BLAS_INTS)
  # Check for Accelerate
  # ^^^^^^^^^^^^^^^^^^^^
  message(STATUS "Testing for LAPACK support via Accelerate")
  set(ACCELERATE_LIBS "-framework Accelerate")
  set(CMAKE_REQUIRED_LIBRARIES ${ACCELERATE_LIBS})
  El_check_function_exists(dpotrf  EL_HAVE_ACCELERATE_NO_UNDER)
  El_check_function_exists(dpotrf_ EL_HAVE_ACCELERATE_UNDER)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(EL_HAVE_ACCELERATE_NO_UNDER OR EL_HAVE_ACCELERATE_UNDER)
    set(EL_HAVE_ACCELERATE TRUE)
  endif()
  if(EL_HAVE_ACCELERATE AND NOT MATH_LIBS_AT_CONFIG 
                        AND NOT EL_DISABLE_APPLE_MATH
                        AND NOT EL_PREFER_OPENBLAS
                        AND NOT EL_PREFER_BLIS_LAPACK)
    set(MATH_LIBS_AT_CONFIG ${ACCELERATE_LIBS})
    message(STATUS "Using Apple Accelerate framework")
  endif()

  # Check for vecLib
  # ^^^^^^^^^^^^^^^^
  message(STATUS "Testing for LAPACK support via vecLib")
  set(VECLIB_LIBS "-framework vecLib")
  set(CMAKE_REQUIRED_LIBRARIES ${VECLIB_LIBS})
  El_check_function_exists(dpotrf  EL_HAVE_VECLIB_NO_UNDER)
  El_check_function_exists(dpotrf_ EL_HAVE_VECLIB_UNDER)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(EL_HAVE_VECLIB_NO_UNDER OR EL_HAVE_VECLIB_UNDER)
    set(EL_HAVE_VECLIB TRUE)
  endif()
  if(EL_HAVE_VECLIB AND NOT MATH_LIBS_AT_CONFIG
                    AND NOT EL_DISABLE_APPLE_MATH
                    AND NOT EL_PREFER_OPENBLAS
                    AND NOT EL_PREFER_BLIS_LAPACK)
    set(MATH_LIBS_AT_CONFIG ${VECLIB_LIBS})
    message(STATUS "Using Apple vecLib framework")
  endif()
endif()

# Check for reference BLAS/LAPACK
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# TODO: Enable support for 64-bit integer version
if(NOT EL_USE_64BIT_BLAS_INTS)
  find_package(LAPACK)
  if(LAPACK_FOUND AND NOT MATH_LIBS_AT_CONFIG)
    if(MSVC OR NOT FORTRAN_WORKS OR 
      (EL_DISABLE_OPENBLAS AND EL_DISABLE_BLIS_LAPACK))
      set(MATH_LIBS_AT_CONFIG "${LAPACK_LINKER_FLAGS};${LAPACK_LIBRARIES}")
    endif()
  endif()
endif()

if(NOT EL_DISABLE_SCALAPACK)
  # Attempt to build ScaLAPACK
  # --------------------------
  include(external_projects/ElMath/ScaLAPACK)
  if(EL_HAVE_SCALAPACK)
    if(CUSTOM_BLAS_SUFFIX OR CUSTOM_LAPACK_SUFFIX)
      set(MATH_LIBS ${MATH_LIBS} ${SCALAPACK_LIBS})
      set(MATH_LIBS_AT_CONFIG ${MATH_LIBS_AT_CONFIG} ${SCALAPACK_LIBS_AT_CONFIG})
    else()
      set(MATH_LIBS ${SCALAPACK_LIBS})
      set(MATH_LIBS_AT_CONFIG ${SCALAPACK_LIBS_AT_CONFIG})
    endif()
  endif()
endif()
if(MATH_LIBS_AT_CONFIG AND NOT MATH_LIBS)
  # We could not extend MATH_LIBS_AT_CONFIG with ScaLAPACK, so use the
  # specified BLAS/LAPACK libraries
  set(MATH_LIBS ${MATH_LIBS_AT_CONFIG})
endif()

# Decide on the BLAS/LAPACK libraries
# -----------------------------------
if(NOT MATH_LIBS AND NOT EL_DISABLE_OPENBLAS
                 AND NOT EL_PREFER_BLIS_LAPACK)
  # Attempt to find/build OpenBLAS
  # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  include(external_projects/ElMath/OpenBLAS)
  if(EL_HAVE_OPENBLAS)
    set(MATH_LIBS ${OPENBLAS_LIBS})
    set(MATH_LIBS_AT_CONFIG ${OPENBLAS_LIBS_AT_CONFIG})
    message("Will use OpenBLAS+LAPACK via MATH_LIBS=${MATH_LIBS}")
  endif()
endif()

if(NOT MATH_LIBS AND NOT EL_DISABLE_BLIS_LAPACK)
  # Attempt to find/build BLIS+LAPACK
  # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  include(external_projects/ElMath/BLIS_LAPACK)
  if(EL_HAVE_BLIS_LAPACK)
    set(MATH_LIBS ${BLIS_LAPACK_LIBS})
    set(MATH_LIBS_AT_CONFIG ${BLIS_LAPACK_LIBS_AT_CONFIG})
    message("Will use BLIS+LAPACK via MATH_LIBS=${MATH_LIBS}")
  endif()
endif()

if(NOT MATH_LIBS)
  # Custom search for reference implementations of BLAS+LAPACK
  # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  message(STATUS "Custom search for reference BLAS and LAPACK")
  set(REFERENCE_REQUIRED LAPACK BLAS)
  find_library(BLAS_LIB NAMES blas PATHS ${MATH_PATHS})
  find_library(LAPACK_LIB NAMES lapack reflapack PATHS ${MATH_PATHS})
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
  set(MATH_LIBS_AT_CONFIG ${MATH_LIBS})
  if(REFERENCE_FOUND)
    message(WARNING "Using reference BLAS/LAPACK; performance will be poor")
  else()
    message(FATAL_ERROR "Could not find BLAS/LAPACK. Please specify MATH_LIBS")
  endif()
endif()

if(EL_BLAS_SUFFIX)
  set(EL_HAVE_BLAS_SUFFIX TRUE)
endif()
if(EL_LAPACK_SUFFIX)
  set(EL_HAVE_LAPACK_SUFFIX TRUE)
endif()

# Check/predict the BLAS and LAPACK underscore conventions
# ========================================================
if(NOT EL_BUILT_BLIS_LAPACK AND NOT EL_BUILT_OPENBLAS)
  set(CMAKE_REQUIRED_FLAGS "${MPI_C_COMPILE_FLAGS}")
  set(CMAKE_REQUIRED_LINKER_FLAGS
    "${MPI_C_LINK_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
  set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
  set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS_AT_CONFIG} ${MPI_C_LIBRARIES})
  # Check BLAS
  # ----------
  # NOTE: MATH_LIBS may involve MPI functionality (e.g., ScaLAPACK) and so
  #       MPI flags should be added for the detection
  if(EL_BLAS_SUFFIX)
    El_check_function_exists(daxpy${EL_BLAS_SUFFIX} EL_HAVE_DAXPY_SUFFIX)
    if(NOT EL_HAVE_DAXPY_SUFFIX)
      message(FATAL_ERROR "daxpy${EL_BLAS_SUFFIX} was not detected")
    endif()
  else()
    El_check_function_exists(daxpy  EL_HAVE_DAXPY)
    El_check_function_exists(daxpy_ EL_HAVE_DAXPY_POST)
    if(EL_HAVE_DAXPY)
      set(EL_HAVE_BLAS_SUFFIX FALSE)
    elseif(EL_HAVE_DAXPY_POST)
      set(EL_HAVE_BLAS_SUFFIX TRUE)
      set(EL_BLAS_SUFFIX _)
    else()
      message(FATAL_ERROR "Could not determine BLAS format.")
    endif()
  endif()
  # Check LAPACK
  # ------------
  if(EL_LAPACK_SUFFIX)
    El_check_function_exists(dpotrf${EL_BLAS_SUFFIX} EL_HAVE_DPOTRF_SUFFIX)
    if(NOT EL_HAVE_DPOTRF_SUFFIX)
      message(FATAL_ERROR "dpotrf${EL_LAPACK_SUFFIX} was not detected")
    endif()
  else()
    El_check_function_exists(dpotrf  EL_HAVE_DPOTRF)
    El_check_function_exists(dpotrf_ EL_HAVE_DPOTRF_POST)
    if(EL_HAVE_DPOTRF)
      set(EL_HAVE_LAPACK_SUFFIX FALSE)
    elseif(EL_HAVE_DPOTRF_POST)
      set(EL_HAVE_LAPACK_SUFFIX TRUE)
      set(EL_LAPACK_SUFFIX _)
    else()
      message(FATAL_ERROR "Could not determine LAPACK format.")
    endif()
  endif()
  # Ensure that we have a relatively new version of LAPACK
  El_check_function_exists(dsyevr${EL_LAPACK_SUFFIX} EL_HAVE_DSYEVR)
  if(NOT EL_HAVE_DSYEVR)
    message(FATAL_ERROR "LAPACK is missing dsyevr${EL_LAPACK_SUFFIX}")
  endif()

  # Check for libFLAME support
  # ==========================
  # TODO: Make this an external project
  El_check_function_exists(FLA_Bsvd_v_opd_var1 EL_HAVE_FLA_BSVD)

  # Clean up the requirements since they cause problems in other Find packages,
  # such as FindThreads
  unset(CMAKE_REQUIRED_FLAGS)
  unset(CMAKE_REQUIRED_LINKER_FLAGS)
  unset(CMAKE_REQUIRED_INCLUDES)
  unset(CMAKE_REQUIRED_LIBRARIES)
endif()

# Check for quad-precision support
# ================================
if(NOT EL_DISABLE_QUAD)
  find_library(QUADMATH_LIB NAMES quadmath PATHS ${MATH_PATHS})
  if(QUADMATH_LIB)
    set(CMAKE_REQUIRED_LIBRARIES ${QUADMATH_LIB})
    set(QUADMATH_CODE
      "#include <complex>
       #include <iostream>
       #include <quadmath.h>
       int main( int argc, char* argv[] )
       {
           __float128 a = 2.0q;

           char aStr[128];
           quadmath_snprintf( aStr, sizeof(aStr), \"%Q\", a );
           std::cout << aStr << std::endl;

           __complex128 y;
           std::complex<__float128> z;

           return 0;    
       }")
    check_cxx_source_compiles("${QUADMATH_CODE}" EL_HAVE_QUADMATH)
    if(EL_HAVE_QUADMATH)
      set(EL_HAVE_QUAD TRUE)
      list(APPEND MATH_LIBS ${QUADMATH_LIB})
      list(APPEND MATH_LIBS_AT_CONFIG ${QUADMATH_LIB})
    else()
      message(WARNING "Found libquadmath but could not use it in C++")
    endif()
    unset(CMAKE_REQUIRED_LIBRARIES)
  endif()
endif()

# Check for QD
# ============
if(NOT EL_DISABLE_QD)
  find_package(QD)
  if(QD_FOUND)
    set(CMAKE_REQUIRED_LIBRARIES ${QD_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES ${QD_INCLUDES})
    set(QD_CODE
      "#include <iostream>
       #include <qd/qd_real.h>
       int main( int argc, char* argv[] )
       {
           double a1=1., a2=2., b1=3., b2=4.;
           dd_real a(a1,a2), b(b1,b2);
           dd_real c = a*b;
           std::cout << \"c=\" << c << std::endl;
           qd_real d(a1,a2,b1,b2);
           std::cout << \"d=\" << d << std::endl;
       }")
    check_cxx_source_compiles("${QD_CODE}" EL_HAVE_QD)
    if(NOT EL_HAVE_QD)
      message(WARNING "Found QD but could not successfully compile with it")
    endif()
    unset(CMAKE_REQUIRED_LIBRARIES)
    unset(CMAKE_REQUIRED_INCLUDES)
  endif()
  if(EL_HAVE_QD)
    set(EL_HAVE_QD TRUE) # Switch from '1' to 'TRUE' for Make
    list(APPEND MATH_LIBS ${QD_LIBRARIES})
    list(APPEND MATH_LIBS_AT_CONFIG ${QD_LIBRARIES})
    message(STATUS "Including ${QD_INCLUDES} to add support for QD")
    include_directories(${QD_INCLUDES})
  endif()
endif()

# Check for GMP, MPFR, *and* MPC support
# ======================================
if(NOT EL_HAVE_MPI_LONG_LONG AND NOT EL_DISABLE_MPFR)
  message("Disabling MPFR since MPI_LONG_LONG was not detected")
endif()
if(EL_HAVE_MPI_LONG_LONG AND NOT EL_DISABLE_MPFR)
  find_package(GMP 6.0.0)
  if(GMP_FOUND)
    set(CMAKE_REQUIRED_LIBRARIES ${GMP_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES ${GMP_INCLUDES})
    set(GMP_CODE
      "#include <gmp.h>
       int main( int argc, char* argv[] )
       {
           gmp_randstate_t randState;
           gmp_randinit_default( randState );
           const long seed = 1024;
           gmp_randseed_ui( randState, seed );
           return 0;
       }")
    check_cxx_source_compiles("${GMP_CODE}" EL_HAVE_GMP)
    if(NOT EL_HAVE_GMP)
      message(WARNING "Found GMP but could not successfully compile with it")
    endif()
    unset(CMAKE_REQUIRED_LIBRARIES)
    unset(CMAKE_REQUIRED_INCLUDES)
  endif()

  if(EL_HAVE_GMP)
    # TODO: See if this requirement could be lowered
    find_package(MPFR 3.1.0)
    if(MPFR_FOUND)
      set(CMAKE_REQUIRED_LIBRARIES ${MPFR_LIBRARIES} ${GMP_LIBRARIES})
      set(CMAKE_REQUIRED_INCLUDES ${MPFR_INCLUDES} ${GMP_INCLUDES})
      set(MPFR_CODE
        "#include <mpfr.h>
         int main( int argc, char* argv[] )
         {
             mpfr_t a;
             mpfr_prec_t prec = 256;
             mpfr_init2( a, prec );
  
             /* Also test that GMP links */
             gmp_randstate_t randState;
             gmp_randinit_default( randState );
             const long seed = 1024;
             gmp_randseed_ui( randState, seed );
             
             return 0;
         }")
      check_cxx_source_compiles("${MPFR_CODE}" EL_HAVE_MPFR)
      if(NOT EL_HAVE_MPFR)
        message(WARNING "Found MPFR but could not successfully compile with it")
      endif()
      unset(CMAKE_REQUIRED_LIBRARIES)
      unset(CMAKE_REQUIRED_INCLUDES)
    endif()
  endif()

  if(EL_HAVE_GMP AND EL_HAVE_MPFR)
    find_package(MPC 1.0.0)
    if(MPC_FOUND) 
      set(CMAKE_REQUIRED_LIBRARIES
        ${MPC_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES})
      set(CMAKE_REQUIRED_INCLUDES
        ${MPC_INCLUDES} ${MPFR_INCLUDES} ${GMP_INCLUDES})
      set(MPC_CODE
        "#include <mpc.h>
         int main( int argc, char* argv[] )
         {
             mpc_t a;
             mpfr_prec_t prec = 256;
             mpc_init2( a, prec );
              
             /* Also test that GMP links */
             gmp_randstate_t randState;
             gmp_randinit_default( randState );
             const long seed = 1024;
             gmp_randseed_ui( randState, seed );
              
             return 0;
         }")
      check_cxx_source_compiles("${MPC_CODE}" EL_HAVE_MPC)
      if(EL_HAVE_MPC)
        set(EL_HAVE_MPC TRUE) # Switch from '1' to 'TRUE' for Make
        list(APPEND MATH_LIBS
          ${MPC_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES})
        list(APPEND MATH_LIBS_AT_CONFIG
          ${MPC_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES})
        message(STATUS "Including ${GMP_INCLUDES}, ${MPFR_INCLUDES}, and ${MPC_INCLUDES} to add support for GMP, MPFR, and MPC")
        include_directories(${GMP_INCLUDES} ${MPFR_INCLUDES} ${MPC_INCLUDES})
      else()
        message(WARNING "Found MPC but could not successfully compile with it")
      endif()
      unset(CMAKE_REQUIRED_LIBRARIES)
      unset(CMAKE_REQUIRED_INCLUDES)
    endif()
  endif()
endif()

if(EL_DISABLE_PARMETIS)
  include(external_projects/ElMath/METIS)
else()
  include(external_projects/ElMath/ParMETIS)
endif()
if(NOT EL_HAVE_METIS)
  message(FATAL_ERROR "METIS support is required for Elemental but existing support was not detected and downloading was prevented")
endif()
