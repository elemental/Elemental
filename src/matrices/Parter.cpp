/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F> 
void Parter( Matrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Parter"))
    P.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    IndexDependentFill
    ( P, [=]( Int i, Int j ) { return F(1)/(F(i)-F(j)+oneHalf); } );
}

template<typename F>
void Parter( AbstractDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Parter"))
    P.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    IndexDependentFill
    ( P, [=]( Int i, Int j ) { return F(1)/(F(i)-F(j)+oneHalf); } );
}

template<typename F>
void Parter( AbstractBlockDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Parter"))
    P.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    IndexDependentFill
    ( P, [=]( Int i, Int j ) { return F(1)/(F(i)-F(j)+oneHalf); } );
}

#define PROTO(F) \
  template void Parter( Matrix<F>& P, Int n ); \
  template void Parter( AbstractDistMatrix<F>& P, Int n ); \
  template void Parter( AbstractBlockDistMatrix<F>& P, Int n );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
