/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_QR_COL_SWAP_HPP
#define EL_QR_COL_SWAP_HPP

namespace El {
namespace qr {

// Update an explicit QR factorization with a neighboring column swap
// (j,j+1).
//
// For each column swap P,
//
//     A = Q R -> A P = Q R P,
//
// and, which introduces a small bulge underneath a single diagonal entry of
// R which can be zeroed with a Givens rotation applied to R from the left
// (and so the inverse Givens rotation must be applied to the columns of Q).
//

template<typename F>
void NeighborColSwap
(       Matrix<F>& Q,
        Matrix<F>& R,
  Int j )
{
    DEBUG_ONLY(CSE cse("qr::NeighborColSwap"))
    typedef Base<F> Real;
    const Int m = Q.Height();

    ColSwap( R, j, j+1 );
    if( j < m-1 )
    {
        Real c; F s;
        lapack::Givens( R.Get(j,j), R.Get(j+1,j), &c, &s );

        auto RBR = R( IR(j,END), IR(j,END) );
        RotateRows( c, s, RBR, 0, 1 );

        auto col1 = Q( ALL, IR(j)   );
        auto col2 = Q( ALL, IR(j+1) );
        RotateCols( c, Conj(s), Q, j, j+1 );
    }
}

// Update an explicit QR factorization with disjoint neighboring column-swaps
// denoted by an integer vector containing each column j to be swapped with
// column j+1.
//
//
// For each column swap P,
//
//     A = Q R -> A P = Q R P,
//
// and, which introduces a small bulge underneath a single diagonal entry of
// R which can be zeroed with a Givens rotation applied to R from the left
// (and so the inverse Givens rotation must be applied to the columns of Q).
//

template<typename F>
void DisjointNeighborColSwaps
(       Matrix<F>& Q,
        Matrix<F>& R,
  const Matrix<Int>& swapInds )
{
    DEBUG_ONLY(CSE cse("qr::DisjointNeighborColSwaps"))
    const Int numSwaps = swapInds.Width();
    for( Int swap=0; swap<numSwaps; ++swap )
        NeighborColSwap( Q, R, swapInds.Get(swap,0) );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_COL_SWAP_HPP
