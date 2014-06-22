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
void Redheffer( Matrix<T>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Redheffer"))
    R.Resize( n, n );
    IndexDependentFill
    ( R, []( Int i, Int j ) 
         { if( j == 0 || ((j+1)%(i+1))==0 ) { return T(1); }
           else                             { return T(0); } } );
}

template<typename T>
void Redheffer( AbstractDistMatrix<T>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Redheffer"))
    R.Resize( n, n );
    IndexDependentFill
    ( R, []( Int i, Int j ) 
         { if( j == 0 || ((j+1)%(i+1))==0 ) { return T(1); }
           else                             { return T(0); } } );
}

template<typename T>
void Redheffer( AbstractBlockDistMatrix<T>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Redheffer"))
    R.Resize( n, n );
    IndexDependentFill
    ( R, []( Int i, Int j ) 
         { if( j == 0 || ((j+1)%(i+1))==0 ) { return T(1); }
           else                             { return T(0); } } );
}

#define PROTO(T) \
  template void Redheffer( Matrix<T>& R, Int n ); \
  template void Redheffer( AbstractDistMatrix<T>& R, Int n ); \
  template void Redheffer( AbstractBlockDistMatrix<T>& R, Int n );

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
