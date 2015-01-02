/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LDL_PIVOTED_BLOCKED_HPP
#define EL_LDL_PIVOTED_BLOCKED_HPP

namespace El {
namespace ldl {
namespace pivot {

template<typename F>
inline void
Blocked
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& p, bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A, Base<F> gamma=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::pivot::Blocked");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();
    if( n == 0 )
    {
        dSub.Resize( 0, 1 );
        p.Resize( 0, 1 );
        return;
    }
    Zeros( dSub, n-1, 1 );

    // Initialize the permutation to the identity
    p.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        p.Set( i, 0, i );

    Matrix<F> XB1, YB1;
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        const Range<Int> indB( k, n ), indBSub( k, n-1 );
        auto dSubB = dSub( indBSub, IR(0,1) );
        auto pB = p( indB, IR(0,1) );
        Panel( A, dSubB, pB, XB1, YB1, nbProp, k, conjugate, pivotType, gamma );
        const Int nb = XB1.Width();

        // Update the bottom-right panel
        const Range<Int> ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 = A( ind2, ind2 );
        auto X21 = XB1( ind2Pan, ind1Pan );
        auto Y21 = YB1( ind2Pan, ind1Pan );
        Trrk( LOWER, NORMAL, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );

        k += nb;
    }
}

template<typename F>
inline void
Blocked
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& dSubPre,
  AbstractDistMatrix<Int>& pPre, bool conjugate=false, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A, Base<F> gamma=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::pivot::Blocked");
        AssertSameGrids( APre, dSubPre, pPre );
        if( APre.Height() != APre.Width() )
            LogicError("A must be square");
    )
    const Int n = APre.Height();
    pPre.Resize( n, 1 );
    if( n == 0 )
    {
        dSubPre.Resize( 0, 1 );
        return;
    }
    dSubPre.Resize( n-1, 1 );

    auto APtr    = ReadWriteProxy<F,MC,MR>( &APre );  auto& A    = *APtr;
    auto dSubPtr = WriteProxy<F,MC,STAR>( &dSubPre ); auto& dSub = *dSubPtr;
    auto pPtr    = WriteProxy<Int,VC,STAR>( &pPre );  auto& p    = *pPtr;

    Zero( dSub );

    // Initialize the permutation to the identity
    if( p.IsLocalCol(0) )
        for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
            p.SetLocal( iLoc, 0, p.GlobalRow(iLoc) );

    const Grid& g = APre.Grid();
    DistMatrix<F,MC,STAR> XB1(g);
    DistMatrix<F,MR,STAR> YB1(g);
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        const Range<Int> indB( k, n ), indBSub( k, n-1 );
        auto dSubB = dSub( indBSub, IR(0,1) );
        auto pB = p( indB, IR(0,1) );
        Panel( A, dSubB, pB, XB1, YB1, nbProp, k, conjugate, pivotType, gamma );
        const Int nb = XB1.Width();

        // Update the bottom-right panel
        const Range<Int> ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 = A( ind2, ind2 );
        auto X21 = XB1( ind2Pan, ind1Pan );
        auto Y21 = YB1( ind2Pan, ind1Pan );
        LocalTrrk( LOWER, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );

        k += nb;
    }
}

} // namespace pivot
} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PIVOTED_BLOCKED_HPP
