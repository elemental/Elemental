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
void OneTwoOne( Matrix<T>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("OneTwoOne"))
    Zeros( A, n, n );
    SetDiagonal( A, T(1), -1 );
    SetDiagonal( A, T(2),  0 );
    SetDiagonal( A, T(1),  1 );
}

template<typename T>
void OneTwoOne( AbstractDistMatrix<T>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("OneTwoOne"))
    Zeros( A, n, n );
    SetDiagonal( A, T(1), -1 );
    SetDiagonal( A, T(2),  0 );
    SetDiagonal( A, T(1),  1 );
}

template<typename T>
void OneTwoOne( AbstractBlockDistMatrix<T>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("OneTwoOne"))
    Zeros( A, n, n );
    SetDiagonal( A, T(1), -1 );
    SetDiagonal( A, T(2),  0 );
    SetDiagonal( A, T(1),  1 );
}

#define PROTO(T) \
  template void OneTwoOne( Matrix<T>& A, Int n ); \
  template void OneTwoOne( AbstractDistMatrix<T>& A, Int n ); \
  template void OneTwoOne( AbstractBlockDistMatrix<T>& A, Int n ); 

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
