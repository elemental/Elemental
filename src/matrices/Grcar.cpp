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
void Grcar( Matrix<T>& A, Int n, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("Grcar"))
    if( k < 0 )
        LogicError("Number of superdiagonals of ones must be non-negative");
    Zeros( A, n, n );
    if( n > 1 )
        SetDiagonal( A, -1, -1 );
    for( Int j=0; j<Min(n,k+1); ++j )
        SetDiagonal( A, 1, j );
}

template<typename T>
void Grcar( AbstractDistMatrix<T>& A, Int n, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("Grcar"))
    if( k < 0 )
        LogicError("Number of superdiagonals of ones must be non-negative");
    Zeros( A, n, n );
    if( n > 1 )
        SetDiagonal( A, -1, -1 );
    for( Int j=0; j<Min(n,k+1); ++j )
        SetDiagonal( A, 1, j );
}

template<typename T>
void Grcar( AbstractBlockDistMatrix<T>& A, Int n, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("Grcar"))
    if( k < 0 )
        LogicError("Number of superdiagonals of ones must be non-negative");
    Zeros( A, n, n );
    if( n > 1 )
        SetDiagonal( A, -1, -1 );
    for( Int j=0; j<Min(n,k+1); ++j )
        SetDiagonal( A, 1, j );
}

#define PROTO(T) \
  template void Grcar( Matrix<T>& A, Int n, Int k ); \
  template void Grcar( AbstractDistMatrix<T>& A, Int n, Int k ); \
  template void Grcar( AbstractBlockDistMatrix<T>& A, Int n, Int k );

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
