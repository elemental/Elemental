/*
   Copyright (c) 2009-2011, Jack Poulson
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

/* Basic variables */
#define CMAKE_BUILD_TYPE @CMAKE_BUILD_TYPE@
#define Elemental_FLAME_VERSION_MAJOR @Elemental_FLAME_VERSION_MAJOR@
#define Elemental_FLAME_VERSION_MINOR @Elemental_FLAME_VERSION_MINOR@
#cmakedefine RELEASE
#cmakedefine BLAS_POST
#cmakedefine LAPACK_POST
#cmakedefine WITHOUT_PMRRR
#cmakedefine WITHOUT_COMPLEX
#cmakedefine AVOID_COMPLEX_MPI

/* Advanced variables */
#cmakedefine CACHE_WARNINGS
#cmakedefine UNALIGNED_WARNINGS
#cmakedefine VECTOR_WARNINGS
#cmakedefine ENABLE_ALL_DISTRIBUTED_DOT
#cmakedefine POOL_MEMORY
#cmakedefine AVOID_OMP_FMA

#endif // ELEMENTAL_CONFIG_H
