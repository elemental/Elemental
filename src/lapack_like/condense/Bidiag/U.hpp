/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BIDIAG_U_HPP
#define EL_BIDIAG_U_HPP

#include "./UUnb.hpp"
#include "./UPan.hpp"

namespace El {
namespace bidiag {

template<typename F>
inline void U( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ )
{
    DEBUG_ONLY(CallStackEntry cse("bidiag::U"))
    const Int m = A.Height();
    const Int n = A.Width();
    DEBUG_ONLY(
        if( m < n ) 
            LogicError("A must be at least as tall as it is wide");
        // Are these requirements necessary?!?
        if( tP.Viewing() || tQ.Viewing() )
            LogicError("tP and tQ must not be views");
    )
    const Int tPHeight = Max(n-1,0);
    const Int tQHeight = n;
    tP.Resize( tPHeight, 1 );
    tQ.Resize( tQHeight, 1 );

    Matrix<F> X, Y;

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k, k+nb ),
                         indB( k, m    ), indR( k, n ),
                         ind2Vert( k+nb, m ), ind2Horz( k+nb, n );

        auto ABR = A( indB,     indR     );
        auto A22 = A( ind2Vert, ind2Horz );

        auto tQ1 = tQ( ind1, IR(0,1) );

        if( A22.Width() > 0 )
        {
            auto tP1 = tP( ind1, IR(0,1) );
            X.Resize( m-k, nb  );
            Y.Resize( nb,  n-k );
            bidiag::UPan( ABR, tP1, tQ1, X, Y );

            auto A12 = A( ind1,     ind2Horz );
            auto A21 = A( ind2Vert, ind1     );
            auto X21 = X( IR(nb,m-k), IR(0,nb)   );
            auto Y12 = Y( IR(0,nb),   IR(nb,n-k) );

            // Set bottom-left entry of A12 to 1
            const F epsilon = A12.Get(nb-1,0);
            A12.Set(nb-1,0,F(1));

            Gemm( NORMAL, NORMAL, F(-1), A21, Y12, F(1), A22 );
            Conjugate( A12 );
            Gemm( NORMAL, NORMAL, F(-1), X21, A12, F(1), A22 );
            Conjugate( A12 );

            // Put back bottom-left entry of A12
            A12.Set(nb-1,0,epsilon);
        }
        else
        {
            auto tP1 = tP( IR(k,k+nb-1), IR(0,1) );
            bidiag::UUnb( ABR, tP1, tQ1 );
        }
    }
}

template<typename F> 
inline void
U
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
        if( m < n ) 
            LogicError("A must be at least as tall as it is wide");
        // Are these requirements necessary?!?
        if( tP.Viewing() || tQ.Viewing() )
            LogicError("tP and tQ must not be views");
    )
    const Grid& g = A.Grid();
    const Int tPHeight = Max(n-1,0);
    const Int tQHeight = n;
    tP.Resize( tPHeight, 1 );
    tQ.Resize( tQHeight, 1 );

    DistMatrix<F> X(g), Y(g);
    DistMatrix<F,MC,STAR> X21_MC_STAR(g);
    DistMatrix<F,MR,STAR> Y12Adj_MR_STAR(g);

    DistMatrix<F,MC,  STAR> AB1_MC_STAR(g);
    DistMatrix<F,STAR,MR  > A1R_STAR_MR(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k, k+nb ),
                         indB( k, m    ), indR( k, n ),
                         ind2Vert( k+nb, m ), ind2Horz( k+nb, n );

        auto A22 = A( ind2Vert, ind2Horz );
        auto ABR = A( indB,     indR     );

        auto tQ1 = tQ( ind1, IR(0,1) );

        if( A22.Width() > 0 )
        {
            X.AlignWith( ABR );
            Y.AlignWith( ABR );
            X.Resize( m-k, nb  );
            Y.Resize( nb,  n-k );

            AB1_MC_STAR.AlignWith( ABR );
            A1R_STAR_MR.AlignWith( ABR );
            AB1_MC_STAR.Resize( m-k, nb  );
            A1R_STAR_MR.Resize( nb,  n-k );

            auto tP1 = tP( ind1, IR(0,1) );
            bidiag::UPan( ABR, tP1, tQ1, X, Y, AB1_MC_STAR, A1R_STAR_MR );

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
            auto tP1 = tP( IR(k,k+nb-1), IR(0,1) );
            bidiag::UUnb( ABR, tP1, tQ1 );
        }
    }
}

} // namespace bidiag
} // namespace El

#endif // ifndef EL_BIDIAG_U_HPP
