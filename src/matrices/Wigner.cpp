/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
void Wigner( Matrix<F>& A, Int n, F mean, Base<F> stddev )
{
    DEBUG_ONLY(CallStackEntry cse("Wigner"))
    Gaussian( A, n, n, mean, stddev );
    MakeHermitian( LOWER, A );
}

template<typename F>
void Wigner( AbstractDistMatrix<F>& A, Int n, F mean, Base<F> stddev )
{
    DEBUG_ONLY(CallStackEntry cse("Wigner"))
    Gaussian( A, n, n, mean, stddev );
    MakeHermitian( LOWER, A );
}

#define PROTO(F) \
  template void Wigner( Matrix<F>& A, Int n, F mean, Base<F> stddev ); \
  template void Wigner \
  ( AbstractDistMatrix<F>& A, Int n, F mean, Base<F> stddev );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
