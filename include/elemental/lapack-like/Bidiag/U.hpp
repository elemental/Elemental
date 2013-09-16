/*
   Copyright (c) 2009-2013, Jack Poulson
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
#ifndef RELEASE
    CallStackEntry entry("bidiag::U");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    const Int tPHeight = Max(n-1,0);
    const Int tQHeight = n;
#ifndef RELEASE
    if( m < n )
        LogicError("A must be at least as tall as it is wide");
#endif
    tP.ResizeTo( tPHeight, 1 );
    tQ.ResizeTo( tQHeight, 1 );

    Matrix<F> x12Adj, w21;

    for( Int k=0; k<n; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto A22     = ViewRange( A, k+1, k+1, m,   n   );
        auto aB1     = ViewRange( A, k,   k,   m,   k+1 );
        auto AB2     = ViewRange( A, k,   k+1, m,   n   );

        // Find tauQ, u, and epsilonQ such that
        //     I - conj(tauQ) | 1 | | 1, u^H | | alpha11 | = | epsilonQ |
        //                    | u |            |    a21  |   |    0     |
        const F tauQ = Reflector( alpha11, a21 );
        const F epsilonQ = alpha11.Get(0,0);
        tQ.Set(k,0,tauQ );

        // Set aB1 = | 1 | and form x12^H := (aB1^H AB2)^H = AB2^H aB1
        //           | u |
        alpha11.Set(0,0,F(1));
        Zeros( x12Adj, a12.Width(), 1 );
        Gemv( ADJOINT, F(1), AB2, aB1, F(0), x12Adj );

        // Update AB2 := AB2 - conj(tauQ) aB1 x12
        //             = AB2 - conj(tauQ) aB1 aB1^H AB2 
        //             = (I - conj(tauQ) aB1 aB1^H) AB2
        Ger( -Conj(tauQ), aB1, x12Adj, AB2 );
        // Put epsilonQ back instead of the temporary value, 1
        alpha11.Set(0,0,epsilonQ);

        if( A22.Width() != 0 )
        {
            // Due to the deficiencies in the BLAS ?gemv routines, this section
            // is easier if we temporarily conjugate a12
            Conjugate( a12 );

            // Expose the subvector we seek to zero, a12R
            Matrix<F> alpha12L, a12R;
            PartitionRight( a12, alpha12L, a12R, 1 );

            // Find tauP, v, and epsilonP such that
            //     I - conj(tauP) | 1 | | 1, v^H | | alpha12L | = | epsilonP |
            //                    | v |            |  a12R^T  |   |    0     |
            const F tauP = Reflector( alpha12L, a12R );
            const F epsilonP = alpha12L.Get(0,0);
            tP.Set(k,0,tauP);

            // Set a12^T = | 1 | and form w21 := A22 a12^T = A22 | 1 |
            //             | v |                                 | v |
            alpha12L.Set(0,0,F(1));
            Zeros( w21, a21.Height(), 1 );
            Gemv( NORMAL, F(1), A22, a12, F(0), w21 );

            // A22 := A22 - tauP w21 conj(a12)
            //      = A22 - tauP A22 a12^T conj(a12)
            //      = A22 (I - tauP a12^T conj(a12))
            //      = A22 conj(I - conj(tauP) a12^H a12)
            // which compensates for the fact that the reflector was generated
            // on the conjugated a12.
            Ger( -tauP, w21, a12, A22 );

            // Put epsilonP back instead of the temporary value, 1
            alpha12L.Set(0,0,epsilonP);

            // Undue the temporary conjugation
            Conjugate( a12 );
        }
    }
}

template<typename F> 
inline void
U( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& tP, DistMatrix<F,STAR,STAR>& tQ )
{
#ifndef RELEASE
    CallStackEntry entry("bidiag::U");
    if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() )
        LogicError
        ("{A,tP,tQ} must be distributed over the same grid");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
#ifndef RELEASE
    if( m < n ) 
        LogicError("A must be at least as tall as it is wide");
    // Are these requirements necessary?!?
    if( tP.Viewing() || tQ.Viewing() )
        LogicError("tP and tQ must not be views");
#endif
    const Grid& g = A.Grid();
    const Int tPHeight = Max(n-1,0);
    const Int tQHeight = n;
    DistMatrix<F,MD,STAR> tPDiag(g), tQDiag(g);
    tPDiag.AlignWithDiagonal( A, 1 );
    tQDiag.AlignWithDiagonal( A, 0 );
    tPDiag.ResizeTo( tPHeight, 1 );
    tQDiag.ResizeTo( tQHeight, 1 );

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
            X.ResizeTo( m-k, nb );
            Y.ResizeTo( n-k, nb );

            AColPan_MC_STAR.AlignWith( A11 );
            ARowPan_STAR_MR.AlignWith( A11 );
            AColPan_MC_STAR.ResizeTo( m-k, nb  );
            ARowPan_STAR_MR.ResizeTo( nb,  n-k );

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
