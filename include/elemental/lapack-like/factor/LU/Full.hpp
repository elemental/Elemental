/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LU_FULL_HPP
#define ELEM_LU_FULL_HPP

#include ELEM_MAXABS_INC
#include ELEM_SCALE_INC
#include ELEM_SWAP_INC
#include ELEM_GERU_INC

namespace elem {
namespace lu {

template<typename F>
inline void
Full( Matrix<F>& A, Matrix<Int>& pPerm, Matrix<Int>& qPerm )
{
    DEBUG_ONLY(CallStackEntry cse("lu::Full"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);

    // Initialize the permutations P and Q
    pPerm.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
        pPerm.Set( i, 0, i );
    qPerm.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
        qPerm.Set( j, 0, j );

    for( Int k=0; k<minDim; ++k )
    {
        // Find the index and value of the pivot candidate
        auto ABR = ViewRange( A, k, k, m, n );
        auto pivot = MaxAbs( ABR );
        const Int iPiv = pivot.indices[0] + k;
        const Int jPiv = pivot.indices[1] + k;

        RowSwap( A,     k, iPiv );
        RowSwap( pPerm, k, iPiv );

        ColSwap( A,     k, jPiv );
        RowSwap( qPerm, k, jPiv );

        // Now we can perform the update of the current panel
        const F alpha11 = A.Get(k,k);
        auto a21 = ViewRange( A, k+1, k,   m,   k+1 );
        auto a12 = ViewRange( A, k,   k+1, k+1, n   );
        auto A22 = ViewRange( A, k+1, k+1, m,   n   );
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha11;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
    }
}

template<typename F,Dist UPerm>
inline void
Full
( DistMatrix<F>& A, 
  DistMatrix<Int,UPerm,STAR>& pPerm, 
  DistMatrix<Int,UPerm,STAR>& qPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::Full");
        if( A.Grid() != pPerm.Grid() || pPerm.Grid() != qPerm.Grid() )
            LogicError("Matrices must be distributed over the same grid");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);

    // Initialize the permutations P and Q
    qPerm.AlignWith( pPerm );
    pPerm.Resize( m, 1 );
    qPerm.Resize( n, 1 );
    for( Int iLoc=0; iLoc<pPerm.LocalHeight(); ++iLoc )
        pPerm.SetLocal( iLoc, 0, pPerm.GlobalRow(iLoc) );
    for( Int jLoc=0; jLoc<qPerm.LocalHeight(); ++jLoc )
        qPerm.SetLocal( jLoc, 0, qPerm.GlobalRow(jLoc) );

    for( Int k=0; k<minDim; ++k )
    {
        // Find the index and value of the pivot candidate
        auto ABR = ViewRange( A, k, k, m, n );
        auto pivot = MaxAbs( ABR );
        const Int iPiv = pivot.indices[0] + k;
        const Int jPiv = pivot.indices[1] + k;

        RowSwap( A,     iPiv, k );
        RowSwap( pPerm, iPiv, k );

        ColSwap( A,     jPiv, k );
        RowSwap( qPerm, jPiv, k );

        // Now we can perform the update of the current panel
        const F alpha11 = A.Get(k,k);
        auto a21 = ViewRange( A, k+1, k,   m,   k+1 );
        auto a12 = ViewRange( A, k,   k+1, k+1, n   );
        auto A22 = ViewRange( A, k+1, k+1, m,   n   );
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha11;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
    }
}

} // namespace lu
} // namespace elem

#endif // ifndef ELEM_LU_FULL_HPP
