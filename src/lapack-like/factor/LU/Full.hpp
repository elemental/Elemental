/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LU_FULL_HPP
#define EL_LU_FULL_HPP

namespace El {
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
        const IndexRange ind1( k, k+1 );
        const IndexRange indB( k, m );
        const IndexRange indR( k, n );
        const IndexRange ind2Vert( k+1, m );
        const IndexRange ind2Horz( k+1, n );

        // Find the index and value of the pivot candidate
        auto ABR = View( A, indB, indR );
        auto pivot = MaxAbs( ABR );
        const Int iPiv = pivot.indices[0] + k;
        const Int jPiv = pivot.indices[1] + k;

        RowSwap( A,     k, iPiv );
        RowSwap( pPerm, k, iPiv );

        ColSwap( A,     k, jPiv );
        RowSwap( qPerm, k, jPiv );

        // Now we can perform the update of the current panel
        const F alpha11 = A.Get(k,k);
        auto a21 = View( A, ind2Vert, ind1     );
        auto a12 = View( A, ind1,     ind2Horz );
        auto A22 = View( A, ind2Vert, ind2Horz );
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha11;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
    }
}

template<typename F>
inline void
Full
( AbstractDistMatrix<F>& APre, 
  AbstractDistMatrix<Int>& pPermPre, AbstractDistMatrix<Int>& qPermPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::Full");
        AssertSameGrids( APre, pPermPre, qPermPre );
    )
    const Int m = APre.Height();
    const Int n = APre.Width();
    const Int minDim = Min(m,n);
    const Grid& g = APre.Grid();

    pPermPre.Resize( m, 1 );
    qPermPre.Resize( n, 1 );

    DistMatrix<F> A(g);
    DistMatrix<Int,VC,STAR> pPerm(g), qPerm(g);
    Copy( APre, A, READ_WRITE_PROXY );
    Copy( pPermPre, pPerm, WRITE_PROXY );
    Copy( qPermPre, qPerm, WRITE_PROXY );

    // Initialize the permutations P and Q
    for( Int iLoc=0; iLoc<pPerm.LocalHeight(); ++iLoc )
        pPerm.SetLocal( iLoc, 0, pPerm.GlobalRow(iLoc) );
    for( Int jLoc=0; jLoc<qPerm.LocalHeight(); ++jLoc )
        qPerm.SetLocal( jLoc, 0, qPerm.GlobalRow(jLoc) );

    for( Int k=0; k<minDim; ++k )
    {
        const IndexRange ind1( k, k+1 );
        const IndexRange ind2Vert( k+1, m );
        const IndexRange ind2Horz( k+1, n );
        const IndexRange indB( k, m );
        const IndexRange indR( k, n );

        // Find the index and value of the pivot candidate
        auto ABR = View( A, indB, indR );
        auto pivot = MaxAbs( ABR );
        const Int iPiv = pivot.indices[0] + k;
        const Int jPiv = pivot.indices[1] + k;

        RowSwap( A,     iPiv, k );
        RowSwap( pPerm, iPiv, k );

        ColSwap( A,     jPiv, k );
        RowSwap( qPerm, jPiv, k );

        // Now we can perform the update of the current panel
        const F alpha11 = A.Get(k,k);
        auto a21 = View( A, ind2Vert, ind1     );
        auto a12 = View( A, ind1,     ind2Horz );
        auto A22 = View( A, ind2Vert, ind2Horz );
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha11;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
    }
    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
    Copy( pPerm, pPermPre, RESTORE_WRITE_PROXY );
    Copy( qPerm, qPermPre, RESTORE_WRITE_PROXY );
}

} // namespace lu
} // namespace El

#endif // ifndef EL_LU_FULL_HPP
