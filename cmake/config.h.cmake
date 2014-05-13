/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CONFIG_H
#define EL_CONFIG_H

/* Build type and version information */
#define EL_GIT_SHA1 "@GIT_SHA1@"
#define EL_VERSION_MAJOR "@EL_VERSION_MAJOR@"
#define EL_VERSION_MINOR "@EL_VERSION_MINOR@"
#define EL_CMAKE_BUILD_TYPE "@CMAKE_BUILD_TYPE@"
#cmakedefine EL_RELEASE

/* C compiler info */
#define EL_CMAKE_C_COMPILER    "@CMAKE_C_COMPILER@"
#define EL_MPI_C_COMPILER      "@MPI_C_COMPILER@"
#define EL_MPI_C_INCLUDE_PATH  "@MPI_C_INCLUDE_PATH@"
#define EL_MPI_C_COMPILE_FLAGS "@MPI_C_COMPILE_FLAGS@"
#define EL_MPI_C_LINK_FLAGS    "@MPI_C_LINK_FLAGS@"
#define EL_MPI_C_LIBRARIES     "@MPI_C_LIBRARIES@"

/* C++ compiler info */
#define EL_CMAKE_CXX_COMPILER    "@CMAKE_CXX_COMPILER@"
#define EL_CXX_FLAGS             "@CXX_FLAGS@"
#define EL_MPI_CXX_COMPILER      "@MPI_CXX_COMPILER@"
#define EL_MPI_CXX_INCLUDE_PATH  "@MPI_CXX_INCLUDE_PATH@"
#define EL_MPI_CXX_COMPILE_FLAGS "@MPI_CXX_COMPILE_FLAGS@"
#define EL_MPI_CXX_LINK_FLAGS    "@MPI_CXX_LINK_FLAGS@"
#define EL_MPI_CXX_LIBRARIES     "@MPI_CXX_LIBRARIES@"

/* Math libraries */
#define EL_MATH_LIBS "@MATH_LIBS@"
#cmakedefine EL_BLAS_POST
#cmakedefine EL_LAPACK_POST
#cmakedefine EL_HAVE_SCALAPACK
#cmakedefine EL_SCALAPACK_POST
#cmakedefine EL_HAVE_FLA_BSVD
#define EL_FORT_LOGICAL @EL_FORT_LOGICAL@
#define EL_FORT_TRUE    @EL_FORT_TRUE@
#define EL_FORT_FALSE   @EL_FORT_FALSE@

/* Basic configuration options */
#define EL_RESTRICT @RESTRICT@
#cmakedefine EL_HAVE_OPENMP
#cmakedefine EL_HAVE_OMP_COLLAPSE
#cmakedefine EL_HAVE_QT5
#cmakedefine EL_HAVE_F90_INTERFACE
#cmakedefine EL_AVOID_COMPLEX_MPI
#cmakedefine EL_HAVE_CXX11RANDOM
#cmakedefine EL_HAVE_STEADYCLOCK
#cmakedefine EL_HAVE_NOEXCEPT
#cmakedefine EL_HAVE_MPI_REDUCE_SCATTER_BLOCK
#cmakedefine EL_HAVE_MPI_IN_PLACE
#cmakedefine EL_HAVE_MPI_LONG_LONG
#cmakedefine EL_HAVE_MPI_COMM_SET_ERRHANDLER
#cmakedefine EL_HAVE_MPI_INIT_THREAD
#cmakedefine EL_HAVE_MPI_QUERY_THREAD
#cmakedefine EL_HAVE_MPI3_NONBLOCKING_COLLECTIVES
#cmakedefine EL_HAVE_MPIX_NONBLOCKING_COLLECTIVES
#cmakedefine EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
#cmakedefine EL_USE_BYTE_ALLGATHERS
#cmakedefine EL_USE_64BIT_INTS

/* Advanced configuration options */
#cmakedefine EL_CACHE_WARNINGS
#cmakedefine EL_UNALIGNED_WARNINGS
#cmakedefine EL_VECTOR_WARNINGS
#cmakedefine EL_POOL_MEMORY
#cmakedefine EL_AVOID_OMP_FMA

#endif /* EL_CONFIG_H */
