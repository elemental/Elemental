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
void MakeIdentity( Matrix<T>& I )
{
    DEBUG_ONLY(CallStackEntry cse("MakeIdentity"))
    Zero( I );
    SetDiagonal( I, T(1) );
}

template<typename T>
void MakeIdentity( AbstractDistMatrix<T>& I )
{
    DEBUG_ONLY(CallStackEntry cse("MakeIdentity"))
    Zero( I );
    SetDiagonal( I, T(1) );
}

template<typename T>
void MakeIdentity( AbstractBlockDistMatrix<T>& I )
{
    DEBUG_ONLY(CallStackEntry cse("MakeIdentity"))
    Zero( I );
    SetDiagonal( I, T(1) );
}

template<typename T>
void Identity( Matrix<T>& I, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Identity"))
    I.Resize( m, n );
    MakeIdentity( I );
}

template<typename T>
void Identity( AbstractDistMatrix<T>& I, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Identity"))
    I.Resize( m, n );
    MakeIdentity( I );
}

template<typename T>
void Identity( AbstractBlockDistMatrix<T>& I, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Identity"))
    I.Resize( m, n );
    MakeIdentity( I );
}

#define PROTO(T) \
  template void MakeIdentity( Matrix<T>& I ); \
  template void MakeIdentity( AbstractDistMatrix<T>& I ); \
  template void MakeIdentity( AbstractBlockDistMatrix<T>& I ); \
  template void Identity( Matrix<T>& I, Int m, Int n ); \
  template void Identity( AbstractDistMatrix<T>& I, Int m, Int n ); \
  template void Identity( AbstractBlockDistMatrix<T>& I, Int m, Int n );

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
