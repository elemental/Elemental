/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BIDIAG_L_HPP
#define EL_BIDIAG_L_HPP

#include "./LUnb.hpp"
#include "./LPan.hpp"

namespace El {
namespace bidiag {

// NOTE: Very little is changed versus the upper case. Perhaps they should be
//       combined.
template<typename F>
inline void L( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ )
{
    DEBUG_ONLY(CallStackEntry cse("bidiag::L"))
    const Int m = A.Height();
    const Int n = A.Width();
    DEBUG_ONLY(
        if( m > n )
            LogicError("A must be at least as wide as it is tall");
        // Are these requirements necessary?!?
        if( tP.Viewing() || tQ.Viewing() )
            LogicError("tP and tQ must not be views");
    )
    const Int tPHeight = m;
    const Int tQHeight = Max(m-1,0);
    tP.Resize( tPHeight, 1 );
    tQ.Resize( tQHeight, 1 );

    Matrix<F> X, Y;

    const Int bsize = Blocksize();
    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind1( k, k+nb ),
                         indB( k, m ), indR( k, n ),
                         ind2Vert( k+nb, m ), ind2Horz( k+nb, n );

        auto A22 = A( ind2Vert, ind2Horz );
        auto ABR = A( indB,     indR     );

        auto tP1 = tP( ind1, IR(0,1) );

        if( A22.Height() > 0 )
        {
            auto A12 = A( ind1,     ind2Horz );
            auto A21 = A( ind2Vert, ind1     );

            auto tQ1 = tQ( ind1, IR(0,1) );
            X.Resize( m-k, nb  );
            Y.Resize( nb,  n-k );
            bidiag::LPan( ABR, tP1, tQ1, X, Y );

            auto X21 = X( IR(nb,m-k), IR(0,nb)   );
            auto Y12 = Y( IR(0,nb),   IR(nb,n-k) );

            // Set top-right entry of A21 to 1
            const F epsilon = A21.Get(0,nb-1);
            A21.Set(0,nb-1,F(1));

            Gemm( NORMAL, NORMAL, F(-1), A21, Y12, F(1), A22 );
            Conjugate( A12 );
            Gemm( NORMAL, NORMAL, F(-1), X21, A12, F(1), A22 );
            Conjugate( A12 );

            // Put back top-right entry of A21
            A21.Set(0,nb-1,epsilon);
        }
        else
        {
            auto tQ1 = tQ( IR(k,k+nb-1), IR(0,1) );
            bidiag::LUnb( ABR, tP1, tQ1 );
        }
    }
}

// NOTE: Very little is different from the upper case. Perhaps they should
//       be combined.
template<typename F> 
inline void
L
( AbstractDistMatrix<F>& APre, 
  AbstractDistMatrix<F>& tPPre, AbstractDistMatrix<F>& tQPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("bidiag::U");
        AssertSameGrids( APre, tPPre, tQPre );
    )

    auto APtr  = ReadWriteProxy<F,MC,MR>( &APre );  auto& A  = *APtr;
    auto tPPtr = WriteProxy<F,STAR,STAR>( &tPPre ); auto& tP = *tPPtr;
    auto tQPtr = WriteProxy<F,STAR,STAR>( &tQPre ); auto& tQ = *tQPtr;

    const Int m = A.Height();
    const Int n = A.Width();
    DEBUG_ONLY(
        if( m > n )
            LogicError("A must be at least as wide as it is tall");
        // Are these requirements necessary?!?
        if( tP.Viewing() || tQ.Viewing() )
            LogicError("tP and tQ must not be views");
    )
    const Grid& g = A.Grid();
    const Int tPHeight = m;
    const Int tQHeight = Max(m-1,0);
    tP.Resize( tPHeight, 1 );
    tQ.Resize( tQHeight, 1 );

    DistMatrix<F> X(g), Y(g);
    DistMatrix<F,MC,STAR> X21_MC_STAR(g);
    DistMatrix<F,MR,STAR> Y12Adj_MR_STAR(g);

    DistMatrix<F,MC,  STAR> AB1_MC_STAR(g);
    DistMatrix<F,STAR,MR  > A1R_STAR_MR(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind1( k, k+nb ),
                         indB( k, m ), indR( k, n ),
                         ind2Vert( k+nb, m ), ind2Horz( k+nb, n );

        auto A22 = A( ind2Vert, ind2Horz );
        auto ABR = A( indB,     indR     );

        auto tP1 = tP( ind1, IR(0,1) );

        if( A22.Height() > 0 )
        {
            X.AlignWith( ABR );
            Y.AlignWith( ABR );
            X.Resize( m-k, nb  );
            Y.Resize( nb,  n-k );

            AB1_MC_STAR.AlignWith( ABR );
            A1R_STAR_MR.AlignWith( ABR );
            AB1_MC_STAR.Resize( m-k, nb  );
            A1R_STAR_MR.Resize( nb,  n-k );

            auto tQ1 = tQ( ind1, IR(0,1) );
            bidiag::LPan( ABR, tP1, tQ1, X, Y, AB1_MC_STAR, A1R_STAR_MR );

            auto X21 = X( IR(nb,m-k), IR(0,nb)   );
            auto Y12 = Y( IR(0,nb),   IR(nb,n-k) );
            X21_MC_STAR.AlignWith( A22 );
            Y12Adj_MR_STAR.AlignWith( A22 );
            X21_MC_STAR = X21;
            Adjoint( Y12, Y12Adj_MR_STAR );

            auto A21_MC_STAR = AB1_MC_STAR( IR(nb,m-k), IR(0,nb)   );
            auto A12_STAR_MR = A1R_STAR_MR( IR(0,nb),   IR(nb,n-k) );

            LocalGemm
            ( NORMAL, ADJOINT, F(-1), A21_MC_STAR, Y12Adj_MR_STAR, F(1), A22 );
            Conjugate( A12_STAR_MR );
            LocalGemm
            ( NORMAL, NORMAL, F(-1), X21_MC_STAR, A12_STAR_MR, F(1), A22 );
        }
        else
        {
            auto tQ1 = tQ( IR(k,k+nb-1), IR(0,1) );
            bidiag::LUnb( ABR, tP1, tQ1 );
        }
    }
}

} // namespace bidiag
} // namespace El

#endif // ifndef EL_LAPACK_CONDENSE_BIDIAG_L_HPP
