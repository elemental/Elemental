/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

void ExplicitPermutation( const Matrix<Int>& p, Matrix<Int>& P )
{
    DEBUG_ONLY(
        CallStackEntry cse("ExplicitPermutation");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
    )
    const Int n = p.Height();
    Zeros( P, n, n );
    for( Int i=0; i<n; ++i )
        P.Set( i, p.Get(i,0), 1 );
}

void ExplicitPermutation
( const AbstractDistMatrix<Int>& p, AbstractDistMatrix<Int>& P )
{
    DEBUG_ONLY(
        CallStackEntry cse("ExplicitPermutation");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
    )
    const Int n = p.Height();
    Zeros( P, n, n );
    for( Int i=0; i<n; ++i )
        P.Set( i, p.Get(i,0), 1 );
}

} // namespace El
