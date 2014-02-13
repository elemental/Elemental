/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BIDIAG_L_HPP
#define ELEM_BIDIAG_L_HPP

#include ELEM_CONJUGATE_INC
#include ELEM_GEMV_INC
#include ELEM_GER_INC
#include ELEM_GEMM_INC

#include ELEM_REFLECTOR_INC

#include "./LUnb.hpp"
#include "./LPan.hpp"

namespace elem {
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
        auto A22 = ViewRange( A, k+nb, k+nb, m, n );
        auto ABR = ViewRange( A, k,    k,    m, n );
        if( A22.Height() > 0 )
        {
            auto tP1 = View( tP, k, 0, nb, 1 );
            auto tQ1 = View( tQ, k, 0, nb, 1 );
            X.Resize( m-k, nb  );
            Y.Resize( nb,  n-k );
            bidiag::LPan( ABR, tP1, tQ1, X, Y );

            auto A12 = ViewRange( A, k,    k+nb, k+nb, n    );
            auto A21 = ViewRange( A, k+nb, k,    m,    k+nb );
            auto X21 = ViewRange( X, nb, 0,  m-k, nb  );
            auto Y12 = ViewRange( Y, 0,  nb, nb,  n-k );

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
            auto tP1 = View( tP, k, 0, nb,   1 );
            auto tQ1 = View( tQ, k, 0, nb-1, 1 );
            bidiag::LUnb( ABR, tP1, tQ1 );
        }
    }
}

// NOTE: Very little is different from the upper case. Perhaps they should
//       be combined.
template<typename F> 
inline void
L( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& tP, DistMatrix<F,STAR,STAR>& tQ )
{
    DEBUG_ONLY(
        CallStackEntry cse("bidiag::U");
        if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() )
            LogicError
            ("{A,tP,tQ} must be distributed over the same grid");
    )
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

        auto A11 = ViewRange( A, k,    k,    k+nb, k+nb );
        auto A12 = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A21 = ViewRange( A, k+nb, k,    m,    k+nb );
        auto A22 = ViewRange( A, k+nb, k+nb, m,    n    );
        auto ABR = ViewRange( A, k,    k,    m,    n    );

        if( A22.Height() > 0 )
        {
            X.AlignWith( A11 );
            Y.AlignWith( A11 );
            X.Resize( m-k, nb  );
            Y.Resize( nb,  n-k );

            AB1_MC_STAR.AlignWith( A11 );
            A1R_STAR_MR.AlignWith( A11 );
            AB1_MC_STAR.Resize( m-k, nb  );
            A1R_STAR_MR.Resize( nb,  n-k );

            auto tP1 = View( tP, k, 0, nb, 1 );
            auto tQ1 = View( tQ, k, 0, nb, 1 );
            bidiag::LPan( ABR, tP1, tQ1, X, Y, AB1_MC_STAR, A1R_STAR_MR );

            auto X21 = ViewRange( X, nb, 0,  m-k, nb  );
            auto Y12 = ViewRange( Y, 0,  nb, nb,  n-k );
            X21_MC_STAR.AlignWith( A21 );
            Y12Adj_MR_STAR.AlignWith( A12 );
            X21_MC_STAR = X21;
            Y12.AdjointColAllGather( Y12Adj_MR_STAR );

            auto A21_MC_STAR = ViewRange( AB1_MC_STAR, nb, 0,  m-k, nb  );
            auto A12_STAR_MR = ViewRange( A1R_STAR_MR, 0,  nb, nb,  n-k );

            LocalGemm
            ( NORMAL, ADJOINT, F(-1), A21_MC_STAR, Y12Adj_MR_STAR, F(1), A22 );
            Conjugate( A12_STAR_MR );
            LocalGemm
            ( NORMAL, NORMAL, F(-1), X21_MC_STAR, A12_STAR_MR, F(1), A22 );
        }
        else
        {
            auto tP1 = View( tP, k, 0, nb,   1 );
            auto tQ1 = View( tQ, k, 0, nb-1, 1 );
            bidiag::LUnb( ABR, tP1, tQ1 );
        }
    }
}

} // namespace bidiag
} // namespace elem

#endif // ifndef ELEM_LAPACK_CONDENSE_BIDIAG_L_HPP
