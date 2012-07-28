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

#include "mpi.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

// If defined, the _OPENMP macro contains the date of the specification
#ifdef _OPENMP
# include <omp.h>
# if _OPENMP >= 200805
#  define COLLAPSE(N) collapse(N)
# else
#  define COLLAPSE(N) 
# endif
#endif

#if defined(BLAS_POST)
#define BLAS(name) name ## _
#else
#define BLAS(name) name
#endif

#if defined(LAPACK_POST)
#define LAPACK(name) name ## _
#else
#define LAPACK(name) name
#endif

// Declare the intertwined core parts of our library
#include "elemental/core/timer_decl.hpp"
#include "elemental/core/memory_decl.hpp"
#include "elemental/core/complex_decl.hpp"
#include "elemental/core/types_decl.hpp"
#include "elemental/core/matrix_decl.hpp"
#include "elemental/core/imports/mpi.hpp"
#include "elemental/core/grid_decl.hpp"
#include "elemental/core/dist_matrix_decl.hpp"
#include "elemental/core/environment_decl.hpp"
#include "elemental/core/indexing_decl.hpp"

#include "elemental/core/imports/blas.hpp"
#include "elemental/core/imports/lapack.hpp"
#include "elemental/core/imports/plcg.hpp"
#ifndef WITHOUT_PMRRR
  #include "elemental/core/imports/pmrrr.hpp"
#endif

// Implement the intertwined parts of the library
#include "elemental/core/timer_impl.hpp"
#include "elemental/core/memory_impl.hpp"
#include "elemental/core/complex_impl.hpp"
#include "elemental/core/types_impl.hpp"
#include "elemental/core/matrix_impl.hpp"
#include "elemental/core/grid_impl.hpp"
#include "elemental/core/dist_matrix_impl.hpp"
#include "elemental/core/environment_impl.hpp"
#include "elemental/core/indexing_impl.hpp"

// Declare and implement the decoupled parts of the core of the library
// (perhaps these should be moved into their own directory?)
#include "elemental/core/partition_decl.hpp"
#include "elemental/core/partition_impl.hpp"
#include "elemental/core/repartition_decl.hpp"
#include "elemental/core/repartition_impl.hpp"
#include "elemental/core/slide_partition_decl.hpp"
#include "elemental/core/slide_partition_impl.hpp"
#include "elemental/core/random_decl.hpp"
#include "elemental/core/random_impl.hpp"
#include "elemental/core/axpy_interface_decl.hpp"
#include "elemental/core/axpy_interface_impl.hpp"
