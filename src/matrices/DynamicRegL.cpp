/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F> 
void DynamicRegL( Matrix<F>& L, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("DynamicRegL"))
    if( n < 3 )
      LogicError("n must be at least three");
    Zeros( L, n, n );
    FillDiagonal( L,  F(1),        0 );
    FillDiagonal( L,  F(11)/F(2), -1 );
    FillDiagonal( L, -F(1)/F(2),  -2 ); 
}

template<typename F>
void DynamicRegL( AbstractDistMatrix<F>& L, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("DynamicRegL"))
    if( n < 3 )
      LogicError("n must be at least three");
    Zeros( L, n, n );
    FillDiagonal( L,  F(1),        0 );
    FillDiagonal( L,  F(11)/F(2), -1 );
    FillDiagonal( L, -F(1)/F(2),  -2 ); 
}

#define EL_NO_INT_PROTO
#define PROTO(F) \
  template void DynamicRegL( Matrix<F>& L, Int n ); \
  template void DynamicRegL( AbstractDistMatrix<F>& L, Int n );

#include "El/macros/Instantiate.h"

} // namespace El
