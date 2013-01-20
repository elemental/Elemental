/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_HPP
#define CORE_HPP

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
#ifdef HAVE_OPENMP
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

// TODO: Think of how to better decouple the following components

// Declare the intertwined core parts of our library
#include "elemental/core/timer_decl.hpp"
#include "elemental/core/memory_decl.hpp"
#include "elemental/core/complex_decl.hpp"
#include "elemental/core/types_decl.hpp"
#include "elemental/core/matrix_forward_decl.hpp"
#include "elemental/core/dist_matrix_forward_decl.hpp"
#include "elemental/core/view_decl.hpp"
#include "elemental/core/matrix_decl.hpp"
#include "elemental/core/imports/mpi.hpp"
#include "elemental/core/grid_decl.hpp"
#include "elemental/core/dist_matrix_decl.hpp"
#include "elemental/core/environment_decl.hpp"
#include "elemental/core/indexing_decl.hpp"

#include "elemental/core/imports/blas.hpp"
#include "elemental/core/imports/lapack.hpp"
#include "elemental/core/imports/flame.hpp"
#include "elemental/core/imports/plcg.hpp"
#ifndef WITHOUT_PMRRR
  #include "elemental/core/imports/pmrrr.hpp"
#endif
#include "elemental/core/imports/choice.hpp"
#include "elemental/core/imports/mpi_choice.hpp"

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
#include "elemental/core/view_impl.hpp"
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

#include "elemental/core/ReduceComm.hpp"

#endif // ifndef CORE_HPP
