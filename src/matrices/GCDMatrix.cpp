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
void GCDMatrix( Matrix<T>& G, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GCDMatrix"))
    G.Resize( m, n );
    IndexDependentFill( G, []( Int i, Int j ) { return T(GCD(i+1,j+1)); } );
}

template<typename T>
void GCDMatrix( AbstractDistMatrix<T>& G, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GCDMatrix"))
    G.Resize( m, n );
    IndexDependentFill( G, []( Int i, Int j ) { return T(GCD(i+1,j+1)); } );
}

template<typename T>
void GCDMatrix( AbstractBlockDistMatrix<T>& G, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GCDMatrix"))
    G.Resize( m, n );
    IndexDependentFill( G, []( Int i, Int j ) { return T(GCD(i+1,j+1)); } );
}

#define PROTO(T) \
  template void GCDMatrix( Matrix<T>& G, Int m, Int n ); \
  template void GCDMatrix( AbstractDistMatrix<T>& G, Int m, Int n ); \
  template void GCDMatrix( AbstractBlockDistMatrix<T>& G, Int m, Int n );

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
