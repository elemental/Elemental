/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_BIDIAG_U_HPP
#define ELEM_LAPACK_BIDIAG_U_HPP

#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/Bidiag/PanelU.hpp"
#include "elemental/lapack-like/Bidiag/UUnb.hpp"
#include "elemental/lapack-like/Reflector.hpp"

namespace elem {
namespace bidiag {

template<typename F>
inline void U( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ )
{
    DEBUG_ONLY(CallStackEntry cse("bidiag::U"))
    // TODO: Sequential blocked implementation
    UUnb( A, tP, tQ );
}

template<typename F> 
inline void
U( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& tP, DistMatrix<F,STAR,STAR>& tQ )
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
        if( m < n ) 
            LogicError("A must be at least as tall as it is wide");
        // Are these requirements necessary?!?
        if( tP.Viewing() || tQ.Viewing() )
            LogicError("tP and tQ must not be views");
    )
    const Grid& g = A.Grid();
    const Int tPHeight = Max(n-1,0);
    const Int tQHeight = n;
    DistMatrix<F,MD,STAR> tPDiag(g), tQDiag(g);
    tPDiag.AlignWithDiagonal( A, 1 );
    tQDiag.AlignWithDiagonal( A, 0 );
    tPDiag.Resize( tPHeight, 1 );
    tQDiag.Resize( tQHeight, 1 );

    DistMatrix<F> X(g), Y(g);
    DistMatrix<F,MC,  STAR> X21_MC_STAR(g), AColPan_MC_STAR(g);
    DistMatrix<F,MR,  STAR> Y21_MR_STAR(g);
    DistMatrix<F,STAR,MR  > ARowPan_STAR_MR(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        
        auto A11 = ViewRange( A, k,    k,    k+nb, k+nb );
        auto A12 = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A21 = ViewRange( A, k+nb, k,    m,    k+nb );
        auto A22 = ViewRange( A, k+nb, k+nb, m,    n    );
        auto ABR = ViewRange( A, k,    k,    m,    n    );

        if( A22.Width() > 0 )
        {
            X.AlignWith( A11 );
            Y.AlignWith( A11 );
            X.Resize( m-k, nb );
            Y.Resize( n-k, nb );

            AColPan_MC_STAR.AlignWith( A11 );
            ARowPan_STAR_MR.AlignWith( A11 );
            AColPan_MC_STAR.Resize( m-k, nb  );
            ARowPan_STAR_MR.Resize( nb,  n-k );

            auto tP1 = View( tPDiag, k, 0, nb, 1 );
            auto tQ1 = View( tQDiag, k, 0, nb, 1 );
            bidiag::PanelU
            ( ABR, tP1, tQ1, X, Y, AColPan_MC_STAR, ARowPan_STAR_MR );

            auto X21 = ViewRange( X, nb, 0, m-k, nb );
            auto Y21 = ViewRange( Y, nb, 0, n-k, nb );
            X21_MC_STAR.AlignWith( A21 );
            Y21_MR_STAR.AlignWith( A12 );
            X21_MC_STAR = X21;
            Y21_MR_STAR = Y21;

            auto A21_MC_STAR = ViewRange( AColPan_MC_STAR, nb, 0, m-k, nb );
            auto A12_STAR_MR = ViewRange( ARowPan_STAR_MR, 0, nb, nb, n-k );
            LocalGemm
            ( NORMAL, ADJOINT, F(-1), A21_MC_STAR, Y21_MR_STAR, F(1), A22 );
            LocalGemm
            ( NORMAL, NORMAL, F(-1), X21_MC_STAR, A12_STAR_MR, F(1), A22 );
        }
        else
        {
            auto tP1 = View( tPDiag, k, 0, nb-1, 1 );
            auto tQ1 = View( tQDiag, k, 0, nb,   1 );
            bidiag::UUnb( ABR, tP1, tQ1 );
        }
    }

    // Redistribute from matrix-diagonal form to fully replicated
    tP = tPDiag;
    tQ = tQDiag;
}

} // namespace bidiag
} // namespace elem

#endif // ifndef ELEM_LAPACK_BIDIAG_U_HPP
