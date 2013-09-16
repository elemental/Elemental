/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_BIDIAG_L_HPP
#define ELEM_LAPACK_BIDIAG_L_HPP

#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/Bidiag/LUnb.hpp"
#include "elemental/lapack-like/Bidiag/PanelL.hpp"
#include "elemental/lapack-like/Reflector.hpp"

namespace elem {
namespace bidiag {

template<typename F>
inline void L( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ )
{
#ifndef RELEASE
    CallStackEntry entry("bidiag::L");
    if( A.Height() > A.Width() )
        LogicError("A must be at least as wide as it is tall");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    const Int tPHeight = m;
    const Int tQHeight = Max(m-1,0);
    tP.ResizeTo( tPHeight, 1 );
    tQ.ResizeTo( tQHeight, 1 );

    // Views
    Matrix<F> alpha21T, a21B;

    // Temporaries
    Matrix<F> x12Adj, w21;

    for( Int k=0; k<m; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto A22     = ViewRange( A, k+1, k+1, m,   n   );
        auto a1R     = ViewRange( A, k,   k,   k+1, n   );
        auto A2R     = ViewRange( A, k+1, k,   m,   n   );

        // Due to deficiencies in the BLAS ?gemv routines, this section is 
        // easier if we temporarily conjugate a1R = | alpha11, a12 |
        Conjugate( a1R );

        // Find tauP, v, and epsilonP such that
        //     I - conj(tauP) | 1 | | 1, v^H | | alpha11 | = | epsilonP | 
        //                    | v |            |   a12^T |   |    0     |
        const F tauP = Reflector( alpha11, a12 );
        const F epsilonP = alpha11.Get(0,0);
        tP.Set(k,0,tauP);

        // Set a1R^T = | 1 | and form w21 := A2R a1R^T = A2R | 1 |
        //             | v |                                 | v |
        alpha11.Set(0,0,F(1));
        Zeros( w21, a21.Height(), 1 );
        Gemv( NORMAL, F(1), A2R, a1R, F(0), w21 );

        // A2R := A2R - tauP w21 conj(a1R)
        //      = A2R - tauP A2R a1R^T conj(a1R)
        //      = A22 (I - tauP a1R^T conj(a1R))
        //      = A22 conj(I - conj(tauP) a1R^H a1R)
        // which compensates for the fact that the reflector was generated
        // on the conjugated a1R.
        Ger( -tauP, w21, a1R, A2R );

        // Put epsilonP back instead of the temporary value, 1
        alpha11.Set(0,0,epsilonP);

        // Undo the temporary conjugation
        Conjugate( a1R );

        if( A22.Height() != 0 )
        {
            // Expose the subvector we seek to zero, a21B
            PartitionDown
            ( a21, alpha21T,
                   a21B, 1 );

            // Find tauQ, u, and epsilonQ such that
            //     I - conj(tauQ) | 1 | | 1, u^H | | alpha21T | = | epsilonQ |
            //                    | u |            |   a21B   | = |    0     |
            const F tauQ = Reflector( alpha21T, a21B );
            const F epsilonQ = alpha21T.Get(0,0);
            tQ.Set(k,0,tauQ);

            // Set a21 = | 1 | and form x12^H = (a21^H A22)^H = A22^H a21
            //           | u |
            alpha21T.Set(0,0,F(1));
            Zeros( x12Adj, a12.Width(), 1 );
            Gemv( ADJOINT, F(1), A22, a21, F(0), x12Adj );

            // A22 := A22 - conj(tauQ) a21 x12 
            //      = A22 - conj(tauQ) a21 a21^H A22
            //      = (I - conj(tauQ) a21 a21^H) A22
            Ger( -Conj(tauQ), a21, x12Adj, A22 );

            // Put epsilonQ back instead of the temporary value, 1
            alpha21T.Set(0,0,epsilonQ);
        }
    }
}

template<typename F> 
inline void
L( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& tP, DistMatrix<F,STAR,STAR>& tQ )
{
#ifndef RELEASE
    CallStackEntry entry("bidiag::L");
    if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() )
        LogicError("{A,tP,tQ} must be distributed over the same grid");
    if( A.Height() > A.Width() )
        LogicError("A must be at least as wide as it is tall");
    // Are these requirements necessary?!?
    if( tP.Viewing() || tQ.Viewing() )
        LogicError("tP and tQ must not be views");
#endif
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int tPHeight = Max(m-1,0);
    const Int tQHeight = m;
    DistMatrix<F,MD,STAR> tPDiag(g), tQDiag(g);
    tPDiag.AlignWithDiagonal( A, -1 );
    tQDiag.AlignWithDiagonal( A, 0 );
    tPDiag.ResizeTo( tPHeight, 1 );
    tQDiag.ResizeTo( tQHeight, 1 );

    DistMatrix<F> X(g), Y(g);
    DistMatrix<F,MC,  STAR> X21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> Y21_MR_STAR(g);
    DistMatrix<F,MC,  STAR> AColPan_MC_STAR(g);
    DistMatrix<F,STAR,MR  > ARowPan_STAR_MR(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);
        const Int nbtP = Min(bsize,m-1-k);

        auto A11 = ViewRange( A,      k,    k,    k+nb, k+nb );
        auto A12 = ViewRange( A,      k,    k+nb, k+nb, n    );
        auto A21 = ViewRange( A,      k+nb, k,    m,    k+nb );
        auto A22 = ViewRange( A,      k+nb, k+nb, m,    n    );
        auto ABR = ViewRange( A,      k,    k,    m,    n    );
        auto tP1 = ViewRange( tPDiag, k,    0,    nbtP, 1    );
        auto tQ1 = ViewRange( tQDiag, k,    0,    nb,   1    );
        
        if( A22.Height() > 0 )
        {
            X.AlignWith( A11 );
            Y.AlignWith( A11 );
            X.ResizeTo( ABR.Height(), nb );
            Y.ResizeTo( ABR.Width(), nb );

            AColPan_MC_STAR.AlignWith( A11 );
            ARowPan_STAR_MR.AlignWith( A11 );
            AColPan_MC_STAR.ResizeTo( ABR.Height(), nb );
            ARowPan_STAR_MR.ResizeTo( nb, ABR.Width() );

            bidiag::PanelL
            ( ABR, tP1, tQ1, X, Y, AColPan_MC_STAR, ARowPan_STAR_MR );

            auto X21 = LockedViewRange( X, nb, 0, ABR.Height(), nb );
            auto Y21 = LockedViewRange( Y, nb, 0, ABR.Width(), nb );
            X21_MC_STAR.AlignWith( A21 );
            Y21_MR_STAR.AlignWith( A12 );
            X21_MC_STAR = X21;
            Y21_MR_STAR = Y21;

            auto A21_MC_STAR = 
                LockedViewRange( AColPan_MC_STAR, nb, 0, ABR.Height(), nb );
            auto A12_STAR_MR = 
                LockedViewRange( ARowPan_STAR_MR, 0, nb, nb, ABR.Width() );
            LocalGemm
            ( NORMAL, ADJOINT, F(-1), A21_MC_STAR, Y21_MR_STAR, F(1), A22 );
            LocalGemm
            ( NORMAL, NORMAL, F(-1), X21_MC_STAR, A12_STAR_MR, F(1), A22 );
        }
        else
        {
            bidiag::LUnb( ABR, tP1, tQ1 );
        }
    }

    // Redistribute from matrix-diagonal form to fully replicated
    tP = tPDiag;
    tQ = tQDiag;
}

} // namespace bidiag
} // namespace elem

#endif // ifndef ELEM_LAPACK_BIDIAG_L_HPP
