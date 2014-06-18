/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Wigner( Matrix<T>& A, Int n, T mean, Base<T> stddev )
{
    DEBUG_ONLY(CallStackEntry cse("Wigner"))
    Gaussian( A, n, n, mean, stddev );
    MakeHermitian( LOWER, A );
}

template<typename T>
void Wigner( AbstractDistMatrix<T>& A, Int n, T mean, Base<T> stddev )
{
    DEBUG_ONLY(CallStackEntry cse("Wigner"))
    Gaussian( A, n, n, mean, stddev );
    MakeHermitian( LOWER, A );
}

#define PROTO(T) \
  template void Wigner( Matrix<T>& A, Int n, T mean, Base<T> stddev ); \
  template void Wigner \
  ( AbstractDistMatrix<T>& A, Int n, T mean, Base<T> stddev );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
