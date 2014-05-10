/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_EXPLICITPERMUTATION_HPP
#define ELEM_LAPACK_EXPLICITPERMUTATION_HPP

#include ELEM_ZEROS_INC

namespace elem {

inline void
ExplicitPermutation( const Matrix<Int>& perm, Matrix<Int>& P )
{
    DEBUG_ONLY(
        CallStackEntry cse("ExplicitPermutation");
        if( perm.Width() != 1 )
            LogicError("perm must be a column vector");
    )
    const Int n = perm.Height();
    Zeros( P, n, n );
    for( Int i=0; i<n; ++i )
        P.Set( i, perm.Get(i,0), 1 );
}

template<Dist UPerm,Dist U,Dist V>
inline void
ExplicitPermutation
( const DistMatrix<Int,UPerm,STAR>& perm, DistMatrix<Int,U,V>& P )
{
    DEBUG_ONLY(
        CallStackEntry cse("ExplicitPermutation");
        if( perm.Width() != 1 )
            LogicError("perm must be a column vector");
    )
    const Int n = perm.Height();
    Zeros( P, n, n );
    for( Int i=0; i<n; ++i )
        P.Set( i, perm.Get(i,0), 1 );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_EXPLICITPERMUTATION_HPP
