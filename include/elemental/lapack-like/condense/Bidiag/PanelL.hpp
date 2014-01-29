/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BIDIAG_PANELL_HPP
#define ELEM_BIDIAG_PANELL_HPP

#include ELEM_ADJOINT_INC
#include ELEM_AXPY_INC
#include ELEM_CONJUGATE_INC
#include ELEM_SCALE_INC
#include ELEM_GEMV_INC

#include ELEM_REFLECTOR_INC

namespace elem {
namespace bidiag {

template<typename F>
inline void
PanelL( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ, Matrix<F>& X, Matrix<F>& Y )
{
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int nX = X.Width();
    DEBUG_ONLY(
        CallStackEntry cse("bidiag::PanelL");
        if( tP.Height() != nX || tP.Width() != 1 )
            LogicError("tP was not the right size");
        if( tQ.Height() != nX || tQ.Width() != 1 )
            LogicError("tQ was not the right size");
        if( mA > nA )
            LogicError("A must be at least as wide as it is tall"); 
        if( mA != X.Height() )
            LogicError("A and X must have the same height");
        if( nA != Y.Width() )
            LogicError("A and Y must have the same width");
        if( Y.Height() != nX )
            LogicError("X is the wrong height");
        if( Y.Width() < nX )
            LogicError("Y must be a row panel");
    )
    typedef Base<F> Real;

    Matrix<F> zT1, z01, z21;
    
    Matrix<Real> d, e;
    d.Resize( nX, 1 );
    e.Resize( nX, 1 );

    for( Int k=0; k<nX; ++k )
    {
        auto aT1      = ViewRange( A, 0,   k,   k+1, k+1 );
        auto AT2      = ViewRange( A, 0,   k+1, k+1, nA  );
        auto A0R      = ViewRange( A, 0,   k,   k,   nA  );
        auto a10      = ViewRange( A, k,   0,   k+1, k   );
        auto alpha11  = ViewRange( A, k,   k,   k+1, k+1 );
        auto a1R      = ViewRange( A, k,   k,   k+1, nA  );
        auto a12      = ViewRange( A, k,   k+1, k+1, nA  );
        auto a21      = ViewRange( A, k+1, k,   mA,  k+1 );
        auto alpha21T = ViewRange( A, k+1, k,   k+2, k+1 );
        auto a21B     = ViewRange( A, k+2, k,   mA,  k+1 );
        auto A20      = ViewRange( A, k+1, 0,   mA,  k   );
        auto A2R      = ViewRange( A, k+1, k,   mA,  nA  );
        auto A22      = ViewRange( A, k+1, k+1, mA,  nA  );

        auto x10 = ViewRange( X, k,   0,   k+1, k   );
        auto XB0 = ViewRange( X, k,   0,   mA,  k   );
        auto X20 = ViewRange( X, k+1, 0,   mA,  k   );
        auto X2L = ViewRange( X, k+1, 0,   mA,  k+1 );
        auto x21 = ViewRange( X, k+1, k,   mA,  k+1 );

        auto y01 = ViewRange( Y, 0, k,   k,   k+1 );
        auto Y0R = ViewRange( Y, 0, k,   k,   nA  );
        auto Y02 = ViewRange( Y, 0, k+1, k,   nA  );
        auto y12 = ViewRange( Y, k, k+1, k+1, nA  );

        // Apply all previous reflectors to a1R:    
        //   a1R := a1R - a10 Y0R         - x10 conj(A0R)
        //        = a1R - (Y0R^T a10^T)^T - (A0R^H x10^T)^T
        Gemv( TRANSPOSE, F(-1), Y0R, a10, F(1), a1R );
        Gemv( ADJOINT,   F(-1), A0R, x10, F(1), a1R );

        // Find tauP and v such that
        //  |alpha11 a12| /I - tauP |1  | |1 conj(v)|\ = |delta 0|
        //                \         |v^T|            /
        const F tauP = RightReflector( alpha11, a12 );
        tP.Set(k,0,tauP);

        // Temporarily set a1R = | 1 v |
        d.Set(k,0,alpha11.GetRealPart(0,0));
        alpha11.Set(0,0,F(1));

        // Form half of the right-reflector using an implicitly-updated A2R:
        // x21 := tauP (A2R - A20 Y0R - X20 conj(A0R)) a1R^T
        //      = tauP (A2R a1R^T - A20 (Y0R a1R^T) - X20 (conj(A0R) a1R^T))
        // -----------------------------------------------------------------
        // x21 := A2R a1R^T
        Zeros( x21, A2R.Height(), 1 );
        Gemv( NORMAL, F(1), A2R, a1R, F(0), x21 );
        // x21 := x21 - A20 (Y0R a1R^T) 
        Gemv( NORMAL, F(1),  Y0R, a1R, z01 );
        Gemv( NORMAL, F(-1), A20, z01, F(1), x21 );
        // x21 := x21 - X20 (conj(A0R) a1R^T)
        //      = x21 - X20 conj(A0R a1R^H)
        Conjugate( a1R );
        Gemv( NORMAL, F(1),  A0R, a1R, z01 );
        Conjugate( a1R );
        Conjugate( z01 );
        Gemv( NORMAL, F(-1), X20, z01, F(1), x21 );
        // x21 := tauP x21
        Scale( tauP, x21 );

        // Apply all previous reflectors to a21:
        //   a21 := a21 - A20 y01 - X2L conj(aT1)
        Gemv( NORMAL, F(-1), A20, y01, F(1), a21 );
        Conjugate( aT1 );
        Gemv( NORMAL, F(-1), X2L, aT1, F(1), a21 );
        Conjugate( aT1 );

        // Find tauQ and u such that
        //  / I - tauQ | 1 | | 1, u^H | \ | alpha21T | = | epsilon |
        //  \          | u |            / |     a21B |   |    0    |
        const F tauQ = LeftReflector( alpha21T, a21B );
        tQ.Set(k,0,tauQ);

        // Temporarily set a21 = | 1 |
        //                       | u |
        e.Set(k,0,alpha21T.GetRealPart(0,0));
        alpha21T.Set(0,0,F(1));

        // Form half of the left-reflector using an implicitly-updated A22:
        // y12 := tauQ a21^H ( A22 - A20 Y02 - X2L conj(AT2) )
        //      = tauQ ( a21^H A22 - (a21^H A20) Y02 - (a21^H X2L) conj(AT2) )
        //      = tauQ ( A22^H a21 - Y02^H (A20^H a21) - AT2^T (X2L^H a21) )^H
        // -------------------------------------------------------------------
        // z21 := A22^H a21
        Gemv( ADJOINT, F(1), A22, a21, z21 );
        // z21 := z21 - Y02^H (A20^H a21)
        Gemv( ADJOINT, F(1),  A20, a21, z01 );
        Gemv( ADJOINT, F(-1), Y02, z01, F(1), z21 );
        // z21 := z21 - AT2^T (X2L^H a21)
        Gemv( ADJOINT, F(1),  X2L, a21, zT1 );
        Gemv( TRANSPOSE, F(-1), AT2, zT1, F(1), z21 );
        // y12 := tauQ z21^H
        Adjoint( z21, y12 );
        Scale( tauQ, y12 );
    }

    // Put back d and e
    auto ATL = View( A, 0, 0, nX, nX );
    auto ATLExpanded = View( A, 0, 0, nX+1, nX );
    ATL.SetRealPartOfDiagonal( d, 0 );
    ATLExpanded.SetRealPartOfDiagonal( e, -1 );
}

template<typename F>
inline void
PanelL
( DistMatrix<F>& A,
  DistMatrix<F,MD,STAR>& tP,
  DistMatrix<F,MD,STAR>& tQ,
  DistMatrix<F>& X,
  DistMatrix<F>& Y,
  DistMatrix<F,MC,STAR>& AColPan_MC_STAR,
  DistMatrix<F,STAR,MR>& ARowPan_STAR_MR )
{
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int nX = X.Width();
    DEBUG_ONLY(
        CallStackEntry cse("bidiag::PanelL");
    )
    LogicError("This routine not yet written");
}

} // namespace bidiag
} // namespace elem

#endif // ifndef ELEM_BIDIAG_PANELL_HPP
