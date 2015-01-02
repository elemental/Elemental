/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Walking through the process of LU factorization with partial pivoting
// for a permutation matrix, which never requires a Schur-complement update,
// yields an algorithm for expressing the inverse of a permutation in terms of 
// a sequence of transpositions in linear time. Note that performing the swaps 
// requires access to the inverse permutation, which can be formed in linear 
// time.

bool PermutationParity( const Matrix<Int>& origPerm )
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

bool PermutationParity( const AbstractDistMatrix<Int>& pPre ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("PermutationParity");
        if( pPre.Width() != 1 )
            LogicError("permutation must be a column vector");
    )

    DistMatrix<Int,VC,STAR> p( pPre ), pInv( pPre.Grid() );
    InvertPermutation( p, pInv );

    bool isOdd = false;
    const Int n = p.Height();
    for( Int k=0; k<n; ++k )
    {
        const Int permVal = p.Get(k,0);
        if( permVal != k )
        {
            isOdd = !isOdd;
            const Int invPermVal = pInv.Get(k,0);
            // We only need to perform half of the swaps
            //      p[k] <-> p[pInv[k]]
            //   pInv[k] <-> pInv[p[k]] 
            // since we will not need to access p[k] and pInv[k] again.
            p.Set( invPermVal, 0, permVal );
            pInv.Set( permVal, 0, invPermVal );
        }
    }
    return isOdd;
}

} // namespace El
