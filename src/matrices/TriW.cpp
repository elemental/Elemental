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
void TriW( Matrix<T>& A, Int m, Int n, T alpha, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("TriW"))
    if( k < 0 )
        LogicError("Number of superdiagonals must be non-negative");
    Zeros( A, n, n );
    SetDiagonal( A, 1 );
    for( Int j=0; j<Min(n-1,k); ++j ) 
        SetDiagonal( A, alpha, j+1 );
}

template<typename T>
void TriW( AbstractDistMatrix<T>& A, Int m, Int n, T alpha, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("TriW"))
    if( k < 0 )
        LogicError("Number of superdiagonals must be non-negative");
    Zeros( A, n, n );
    SetDiagonal( A, 1 );
    for( Int j=0; j<Min(n-1,k); ++j )
        SetDiagonal( A, alpha, j+1 );
}

template<typename T>
void TriW( AbstractBlockDistMatrix<T>& A, Int m, Int n, T alpha, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("TriW"))
    if( k < 0 )
        LogicError("Number of superdiagonals must be non-negative");
    Zeros( A, n, n );
    SetDiagonal( A, 1 );
    for( Int j=0; j<Min(n-1,k); ++j )
        SetDiagonal( A, alpha, j+1 );
}

#define PROTO(T) \
  template void TriW \
  ( Matrix<T>& A, Int m, Int n, T alpha, Int k ); \
  template void TriW \
  ( AbstractDistMatrix<T>& A, Int m, Int n, T alpha, Int k ); \
  template void TriW \
  ( AbstractBlockDistMatrix<T>& A, Int m, Int n, T alpha, Int k );

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
