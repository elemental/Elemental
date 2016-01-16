/*
   Copyright (c) 2009-2016, Jack Poulson
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
    DEBUG_ONLY(CSE cse("Parter"))
    P.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    auto parterFill = [=]( Int i, Int j ) { return F(1)/(F(i)-F(j)+oneHalf); };
    IndexDependentFill( P, function<F(Int,Int)>(parterFill) );
}

template<typename F>
void Parter( AbstractDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CSE cse("Parter"))
    P.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    auto parterFill = [=]( Int i, Int j ) { return F(1)/(F(i)-F(j)+oneHalf); };
    IndexDependentFill( P, function<F(Int,Int)>(parterFill) );
}

#define PROTO(F) \
  template void Parter( Matrix<F>& P, Int n ); \
  template void Parter( AbstractDistMatrix<F>& P, Int n );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
