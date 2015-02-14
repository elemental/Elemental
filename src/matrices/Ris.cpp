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
void Ris( Matrix<F>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ris"))
    R.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    auto risFill = [=]( Int i, Int j ) { return oneHalf/(F(n-i-j)-oneHalf); };
    IndexDependentFill( R, function<F(Int,Int)>(risFill) );
}

template<typename F>
void Ris( AbstractDistMatrix<F>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ris"))
    R.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    auto risFill = [=]( Int i, Int j ) { return oneHalf/(F(n-i-j)-oneHalf); };
    IndexDependentFill( R, function<F(Int,Int)>(risFill) );
}

template<typename F>
void Ris( AbstractBlockDistMatrix<F>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ris"))
    R.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    auto risFill = [=]( Int i, Int j ) { return oneHalf/(F(n-i-j)-oneHalf); };
    IndexDependentFill( R, function<F(Int,Int)>(risFill) );
}

#define PROTO(F) \
  template void Ris( Matrix<F>& R, Int n ); \
  template void Ris( AbstractDistMatrix<F>& R, Int n ); \
  template void Ris( AbstractBlockDistMatrix<F>& R, Int n );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
