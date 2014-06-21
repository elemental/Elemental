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
void MinIJ( Matrix<T>& M, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("MinIJ"))
    M.Resize( n, n );
    IndexDependentFill( M, []( Int i, Int j ) { return T(Min(i+1,j+1)); } );
}

template<typename T>
void MinIJ( AbstractDistMatrix<T>& M, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("MinIJ"))
    M.Resize( n, n );
    IndexDependentFill( M, []( Int i, Int j ) { return T(Min(i+1,j+1)); } );
}

template<typename T>
void MinIJ( AbstractBlockDistMatrix<T>& M, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("MinIJ"))
    M.Resize( n, n );
    IndexDependentFill( M, []( Int i, Int j ) { return T(Min(i+1,j+1)); } );
}

#define PROTO(T) \
  template void MinIJ( Matrix<T>& M, Int n ); \
  template void MinIJ( AbstractDistMatrix<T>& M, Int n ); \
  template void MinIJ( AbstractBlockDistMatrix<T>& M, Int n );

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
