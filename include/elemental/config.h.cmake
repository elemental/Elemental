/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELEMENTAL_CONFIG_H
#define ELEMENTAL_CONFIG_H 1

/* Build type and version information */
#define CMAKE_BUILD_TYPE @CMAKE_BUILD_TYPE@
#define Elemental_VERSION_MAJOR @Elemental_VERSION_MAJOR@
#define Elemental_VERSION_MINOR @Elemental_VERSION_MINOR@

/* C compiler info */
#define CMAKE_C_COMPILER    @CMAKE_C_COMPILER@
#define MPI_C_COMPILER      @MPI_C_COMPILER@
#define MPI_C_INCLUDE_PATH  @MPI_C_INCLUDE_PATH@
#define MPI_C_COMPILE_FLAGS @MPI_C_COMPILE_FLAGS@
#define MPI_C_LINK_FLAGS    @MPI_C_LINK_FLAGS@
#define MPI_C_LIBRARIES     @MPI_C_LIBRARIES@

/* C++ compiler info */
#define CMAKE_CXX_COMPILER    @CMAKE_CXX_COMPILER@
#define MPI_CXX_COMPILER      @MPI_CXX_COMPILER@
#define MPI_CXX_INCLUDE_PATH  @MPI_CXX_INCLUDE_PATH@
#define MPI_CXX_COMPILE_FLAGS @MPI_CXX_COMPILE_FLAGS@
#define MPI_CXX_LINK_FLAGS    @MPI_CXX_LINK_FLAGS@
#define MPI_CXX_LIBRARIES     @MPI_CXX_LIBRARIES@

/* Math libraries */
#define MATH_LIBS @MATH_LIBS@

/* Basic configuration options */
#define RESTRICT @RESTRICT@
#cmakedefine RELEASE
#cmakedefine BLAS_POST
#cmakedefine LAPACK_POST
#cmakedefine WITHOUT_PMRRR
#cmakedefine DISABLE_SCALAR_WRAPPER
#cmakedefine AVOID_COMPLEX_MPI
#cmakedefine HAVE_REDUCE_SCATTER_BLOCK
#cmakedefine HAVE_MPI3_NONBLOCKING_COLLECTIVES
#cmakedefine HAVE_MPIX_NONBLOCKING_COLLECTIVES
#cmakedefine REDUCE_SCATTER_BLOCK_VIA_ALL_REDUCE
#cmakedefine USE_BYTE_ALLGATHERS

/* Advanced configuration options */
#cmakedefine CACHE_WARNINGS
#cmakedefine UNALIGNED_WARNINGS
#cmakedefine VECTOR_WARNINGS
#cmakedefine POOL_MEMORY
#cmakedefine AVOID_OMP_FMA

#endif // ELEMENTAL_CONFIG_H
