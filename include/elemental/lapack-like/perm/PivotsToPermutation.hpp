/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_PIVOTSTOPERMUTATION_HPP
#define ELEM_LAPACK_PIVOTSTOPERMUTATION_HPP

#include ELEM_SWAP_INC
#include ELEM_MAXNORM_INC

namespace elem {

inline void
PivotsToPermutation
( const Matrix<Int>& pivots, Matrix<Int>& perm, Int offset=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("PivotsToPermutation");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
    )

    const Int n = pivots.Height();
    if( n == 0 )
    {
        perm.Resize( 0, 1 );
        return;
    }
 
    // Initialize to the identity permutation
    const Int range = MaxNorm( pivots ) + 1 - offset;
    perm.Resize( range, 1 );
    for( Int i=0; i<range; ++i )
        perm.Set( i, 0, i );

    // Track the location of the nonzero column in each row of the permutation
    // NOTE: Assuming that we have enough memory, it would be faster to perform
    //       this procedure sequentially
    for( Int i=0; i<n; ++i )
    {
        const Int j = pivots.Get(i,0)-offset;
        RowSwap( perm, i, j );
    }
}

template<Dist U,Dist UPerm>
inline void
PivotsToPermutation
( const DistMatrix<Int,U,STAR>& pivots, DistMatrix<Int,UPerm,STAR>& perm, 
  Int offset=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("PivotsToInversePermutation");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
    )

    const Int n = pivots.Height();
    if( n == 0 )
    {
        perm.Resize( 0, 1 );
        return;
    }

    // Initialize to the identity permutation
    const Int range = MaxNorm( pivots ) + 1 - offset;
    perm.SetGrid( pivots.Grid() );
    perm.AlignWith( pivots );
    perm.Resize( range, 1 );
    for( Int iLoc=0; iLoc<perm.LocalHeight(); ++iLoc )
        perm.SetLocal( iLoc, 0, perm.GlobalRow(iLoc) );

    // Track the location of the nonzero column in each row of the permutation
    // NOTE: Assuming that we have enough memory, it would be faster to perform
    //       this procedure sequentially
    for( Int i=0; i<n; ++i )
    {
        const Int j = pivots.Get(i,0)-offset;
        RowSwap( perm, i, j );
    }
}

inline void
PivotsToInversePermutation
( const Matrix<Int>& pivots, Matrix<Int>& invPerm, Int offset=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("PivotsToInversePermutation");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
    )

    const Int n = pivots.Height();
    if( n == 0 )
    {
        invPerm.Resize( 0, 1 );
        return;
    }
 
    // Initialize to the identity permutation
    const Int range = MaxNorm( pivots ) + 1 - offset;
    invPerm.Resize( range, 1 );
    for( Int i=0; i<range; ++i )
        invPerm.Set( i, 0, i );

    // Track the location of the nonzero in each column of the permutation
    // NOTE: Assuming that we have enough memory, it would be faster to perform
    //       this procedure sequentially
    for( Int i=n-1; i>=0; --i )
    {
        const Int j = pivots.Get(i,0)-offset;
        RowSwap( invPerm, i, j );
    }
}

template<Dist U,Dist UPerm>
inline void
PivotsToInversePermutation
( const DistMatrix<Int,U,STAR>& pivots, DistMatrix<Int,UPerm,STAR>& invPerm, 
  Int offset=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("PivotsToInversePermutation");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
    )

    const Int n = pivots.Height();
    if( n == 0 )
    {
        invPerm.Resize( 0, 1 );
        return;
    }

    // Initialize to the identity permutation
    const Int range = MaxNorm( pivots ) + 1 - offset;
    invPerm.SetGrid( pivots.Grid() );
    invPerm.AlignWith( pivots );
    invPerm.Resize( range, 1 );
    for( Int iLoc=0; iLoc<invPerm.LocalHeight(); ++iLoc )
        invPerm.SetLocal( iLoc, 0, invPerm.GlobalRow(iLoc) );

    // Track the location of the nonzero in each column of the permutation
    // NOTE: Assuming that we have enough memory, it would be faster to perform
    //       this procedure sequentially
    for( Int i=n-1; i>=0; --i )
    {
        const Int j = pivots.Get(i,0)-offset;
        RowSwap( invPerm, i, j );
    }
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_PIVOTSTOPERMUTATION_HPP
