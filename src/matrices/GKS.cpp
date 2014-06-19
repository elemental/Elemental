/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// The Golub Klema Stewart matrix is upper-triangular with 1/sqrt(j) on its 
// j'th diagonal entry and -1/sqrt(j) elsewhere in the upper triangle.
// 
// It was originally introduced as an example of where greedy RRQR fails.

namespace El {

template<typename F>
void GKS( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GKS"))
    A.Resize( n, n );
    IndexDependentFill
    ( A, []( Int i, Int j ) 
         { if( i < j )       { return -F(1)/Sqrt(F(j)); }
           else if( i == j ) { return  F(1)/Sqrt(F(j)); }
           else              { return  F(0);            } } );
}

template<typename F>
void GKS( AbstractDistMatrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GKS"))
    A.Resize( n, n );
    IndexDependentFill
    ( A, []( Int i, Int j )
         { if( i < j )       { return -F(1)/Sqrt(F(j)); }
           else if( i == j ) { return  F(1)/Sqrt(F(j)); }
           else              { return  F(0);            } } );
}

template<typename F>
void GKS( AbstractBlockDistMatrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GKS"))
    A.Resize( n, n );
    IndexDependentFill
    ( A, []( Int i, Int j )
         { if( i < j )       { return -F(1)/Sqrt(F(j)); }
           else if( i == j ) { return  F(1)/Sqrt(F(j)); }
           else              { return  F(0);            } } );
}

#define PROTO(F) \
  template void GKS( Matrix<F>& A, Int n ); \
  template void GKS( AbstractDistMatrix<F>& A, Int n ); \
  template void GKS( AbstractBlockDistMatrix<F>& A, Int n );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
