/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PERMUTATIONPARITY_HPP
#define ELEM_PERMUTATIONPARITY_HPP

#include ELEM_INVERTPERMUTATION_INC

namespace elem {

// Walking through the process of LU factorization with partial pivoting
// for a permutation matrix, which never requires a Schur-complement update,
// yields an algorithm for expressing the inverse of a permutation in terms of 
// a sequence of transpositions in linear time. Note that performing the swaps 
// requires access to the inverse permutation, which can be formed in linear 
// time.

inline bool
PermutationParity( const Matrix<Int>& origPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("PermutationParity");
        if( origPerm.Width() != 1 )
            LogicError("permutation must be a column vector");
    )

    Matrix<Int> perm( origPerm );

    Matrix<Int> invPerm;
    InvertPermutation( perm, invPerm );

    bool isOdd = false;
    const Int n = perm.Height();
    for( Int k=0; k<n; ++k )
    {
        const Int permVal = perm.Get(k,0);
        if( permVal != k )
        {
            isOdd = !isOdd;
            const Int invPermVal = invPerm.Get(k,0);
            // We only need to perform half of the swaps
            //      perm[k] <-> perm[invPerm[k]]
            //   invPerm[k] <-> invPerm[perk[k]] 
            // since we will not need to access perm[k] and invPerm[k] again.
            perm.Set( invPermVal, 0, permVal );
            invPerm.Set( permVal, 0, invPermVal );
        }
    }
    return isOdd;
}

template<Dist UPerm>
inline bool
PermutationParity( const DistMatrix<Int,UPerm,STAR>& origPerm ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("PermutationParity");
        if( origPerm.Width() != 1 )
            LogicError("permutation must be a column vector");
    )

    DistMatrix<Int,UPerm,STAR> perm( origPerm );

    DistMatrix<Int,UPerm,STAR> invPerm( origPerm.Grid() );
    InvertPermutation( perm, invPerm );

    bool isOdd = false;
    const Int n = perm.Height();
    for( Int k=0; k<n; ++k )
    {
        const Int permVal = perm.Get(k,0);
        if( permVal != k )
        {
            isOdd = !isOdd;
            const Int invPermVal = invPerm.Get(k,0);
            // We only need to perform half of the swaps
            //      perm[k] <-> perm[invPerm[k]]
            //   invPerm[k] <-> invPerm[perk[k]] 
            // since we will not need to access perm[k] and invPerm[k] again.
            perm.Set( invPermVal, 0, permVal );
            invPerm.Set( permVal, 0, invPermVal );
        }
    }
    return isOdd;
}

} // namespace elem

#endif // ifndef ELEM_PERMUTATIONPARITY_HPP
