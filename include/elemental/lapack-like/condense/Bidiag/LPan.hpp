/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BIDIAG_LPAN_HPP
#define ELEM_BIDIAG_LPAN_HPP

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
LPan( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ, Matrix<F>& X, Matrix<F>& Y )
{
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int nX = X.Width();
    DEBUG_ONLY(
        CallStackEntry cse("bidiag::LPan");
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
        auto a01      = ViewRange( A, 0,   k,   k,   k+1 );
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
        //        = a21 - A20 y01 - (X20 conj(a01) + x21*1)
        Gemv( NORMAL, F(-1), A20, y01, F(1), a21 );
        Conjugate( a01 );
        Gemv( NORMAL, F(-1), X20, a01, F(1), a21 );
        Conjugate( a01 );
        Axpy( F(-1), x21, a21 );

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
LPan
( DistMatrix<F>& A,
  DistMatrix<F,STAR,STAR>& tP,
  DistMatrix<F,STAR,STAR>& tQ,
  DistMatrix<F>& X,
  DistMatrix<F>& Y,
  DistMatrix<F,MC,  STAR>& AL_MC_STAR,
  DistMatrix<F,STAR,MR  >& AT_STAR_MR )
{
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int nX = X.Width();
    DEBUG_ONLY(
        CallStackEntry cse("bidiag::LPan");
        if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() || 
            tQ.Grid() != X.Grid() || X.Grid() != Y.Grid() ||
            Y.Grid() != AL_MC_STAR.Grid() || 
            Y.Grid() != AT_STAR_MR.Grid() )
            LogicError("Grids must match");
        if( A.ColAlign() != X.ColAlign() ||
            A.RowAlign() != X.RowAlign() )
            LogicError("A and X must be aligned");
        if( A.ColAlign() != Y.ColAlign() ||
            A.RowAlign() != Y.RowAlign() )
            LogicError("A and Y must be aligned");
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
    const Grid& g = A.Grid();

    DistMatrix<F,MC,  STAR> z01_MC_STAR(g),
                            zT1_MC_STAR(g),
                            z21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> a01_MR_STAR(g),
                            y01_MR_STAR(g),
                            z01_MR_STAR(g),
                            zT1_MR_STAR(g),
                            z21_MR_STAR(g); 
    DistMatrix<F,STAR,MC  > a10_STAR_MC(g),
                            x10_STAR_MC(g);
    DistMatrix<F,STAR,MR  > z1R_STAR_MR(g);
    
    DistMatrix<Real,MD,STAR> d(g), e(g);
    d.SetRoot( A.DiagonalRoot( 0) );
    e.SetRoot( A.DiagonalRoot(-1) );
    d.AlignCols( A.DiagonalAlign( 0) );
    e.AlignCols( A.DiagonalAlign(-1) );
    d.Resize( nX, 1 );
    e.Resize( nX, 1 );

    for( Int k=0; k<nX; ++k )
    {
        auto a01      = ViewRange( A, 0,   k,   k,   k+1 );
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

        auto a1R_STAR_MR = ViewRange( AT_STAR_MR, k,   k, k+1, nA  );
        auto a21_MC_STAR = ViewRange( AL_MC_STAR, k+1, k, mA,  k+1 ); 

        auto x10 = ViewRange( X, k,   0,   k+1, k   );
        auto XB0 = ViewRange( X, k,   0,   mA,  k   );
        auto X20 = ViewRange( X, k+1, 0,   mA,  k   );
        auto X2L = ViewRange( X, k+1, 0,   mA,  k+1 );
        auto x21 = ViewRange( X, k+1, k,   mA,  k+1 );

        auto y01 = ViewRange( Y, 0, k,   k,   k+1 );
        auto Y0R = ViewRange( Y, 0, k,   k,   nA  );
        auto Y02 = ViewRange( Y, 0, k+1, k,   nA  );
        auto y12 = ViewRange( Y, k, k+1, k+1, nA  );

        auto delta1   = View( d,  k, 0, 1, 1 );
        auto epsilon1 = View( e,  k, 0, 1, 1 );

        // Apply all previous reflectors to a1R:    
        //   a1R := a1R - a10 Y0R         - x10 conj(A0R)
        //        = a1R - (Y0R^T a10^T)^T - (A0R^H x10^T)^T
        if( k > 0 )
        {
            a10_STAR_MC.AlignWith( Y0R );
            x10_STAR_MC.AlignWith( A0R );
            a10_STAR_MC = a10;
            x10_STAR_MC = x10;

            z1R_STAR_MR.AlignWith( a1R );
            Zeros( z1R_STAR_MR, 1, a1R.Width() );
            // z1R[* ,MR] := (Y0R^T[MR,MC] a10^T[MC,* ])^T +
            //               (A0R^H[MR,MC] x10^T[MC,* ])^T
            LocalGemv( TRANSPOSE, F(1), Y0R, a10_STAR_MC, F(0), z1R_STAR_MR );
            LocalGemv( ADJOINT,   F(1), A0R, x10_STAR_MC, F(1), z1R_STAR_MR ); 
            // Sum the partial contributions and subtract from a1R
            a1R.ColSumScatterUpdate( F(-1), z1R_STAR_MR );
        }

        // Find tauP and v such that
        //  |alpha11 a12| /I - tauP |1  | |1 conj(v)|\ = |delta 0|
        //                \         |v^T|            /
        const F tauP = RightReflector( alpha11, a12 );
        tP.Set(k,0,tauP);

        // Temporarily set a1R = | 1 v |
        if( alpha11.IsLocal(0,0) )
        {
            delta1.SetLocal(0,0,alpha11.GetLocalRealPart(0,0));
            alpha11.SetLocal(0,0,F(1));
        }

        // Form half of the right-reflector using an implicitly-updated A2R:
        // x21 := tauP (A2R - A20 Y0R - X20 conj(A0R)) a1R^T
        //      = tauP (A2R a1R^T - A20 (Y0R a1R^T) - X20 (conj(A0R) a1R^T))
        // -----------------------------------------------------------------
        a1R_STAR_MR = a1R;

        // z21[MC,* ] := A2R[MC,MR] a1R^T[MR,* ]
        z21_MC_STAR.AlignWith( A2R );
        Zeros( z21_MC_STAR, A2R.Height(), 1 );
        LocalGemv( NORMAL, F(1), A2R, a1R_STAR_MR, F(0), z21_MC_STAR );

        // z01[MR,* ] := (Y01 a1R^T)[MR,* ]
        z01_MC_STAR.AlignWith( Y0R );
        Zeros( z01_MC_STAR, Y0R.Height(), 1 );
        LocalGemv( NORMAL, F(1), Y0R, a1R_STAR_MR, F(0), z01_MC_STAR );
        z01_MC_STAR.SumOver( Y0R.RowComm() );
        z01_MR_STAR.AlignWith( A20 );
        z01_MR_STAR = z01_MC_STAR;
        // z21[MC,* ] -= A20[MC,MR] z01[MR,* ]
        LocalGemv( NORMAL, F(-1), A20, z01_MR_STAR, F(1), z21_MC_STAR );

        // z01[MR,* ] := conj(A0R a1R^H)[MR,* ]
        //             = (conj(A0R) a1R^T)[MR,* ]
        z01_MC_STAR.AlignWith( A0R );
        Zeros( z01_MC_STAR, A0R.Height(), 1 );
        Conjugate( a1R_STAR_MR );
        LocalGemv( NORMAL, F(1), A0R, a1R_STAR_MR, F(0), z01_MC_STAR );
        z01_MC_STAR.SumOver( A0R.RowComm() );
        Conjugate( a1R_STAR_MR );
        Conjugate( z01_MC_STAR );
        z01_MR_STAR.AlignWith( X20 );
        z01_MR_STAR = z01_MC_STAR;
        // z21[MC,* ] -= X20[MC,MR] z01[MR,* ] 
        LocalGemv( NORMAL, F(-1), X20, z01_MR_STAR, F(1), z21_MC_STAR );

        // Finally perform the row summation and then scale by tauP
        x21.RowSumScatterFrom( z21_MC_STAR );
        Scale( tauP, x21 );

        // Apply all previous reflectors to a21:
        //   a21 := a21 - A20 y01 - X2L conj(aT1)
        //        = a21 - A20 y01 - (X20 conj(a01) + x21*1)
        // ------------------------------------------------
        // a21 := a21 - A20 y01 (do not sum over rows yet)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        y01_MR_STAR.AlignWith( A20 );
        y01_MR_STAR = y01;
        // z21[MC,* ] := A20[MC,MR] y01[MR,* ]
        z21_MC_STAR.AlignWith( A20 );
        Zeros( z21_MC_STAR, A20.Height(), 1 );
        LocalGemv( NORMAL, F(1), A20, y01_MR_STAR, F(0), z21_MC_STAR );

        // a21 := a21 - X20 conj(a01) (and incorporate last update into sum)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        a01_MR_STAR.AlignWith( X20 );
        a01_MR_STAR = a01;
        // z21[MC,* ] += X20[MC,MR] conj(a01)[MR,* ]
        Conjugate( a01_MR_STAR );
        LocalGemv( NORMAL, F(1), X20, a01_MR_STAR, F(1), z21_MC_STAR );
        Conjugate( a01_MR_STAR );
        // Sum the partial contributions from the past two updates
        a21.RowSumScatterUpdate( F(-1), z21_MC_STAR );

        // a21 := a21 - x21
        // ^^^^^^^^^^^^^^^^
        Axpy( F(-1), x21, a21 );

        // Find tauQ and u such that
        //  / I - tauQ | 1 | | 1, u^H | \ | alpha21T | = | epsilon |
        //  \          | u |            / |     a21B |   |    0    |
        const F tauQ = LeftReflector( alpha21T, a21B );
        tQ.Set(k,0,tauQ);

        // Temporarily set a21 = | 1 |
        //                       | u |
        if( alpha21T.IsLocal(0,0) )
        {
            epsilon1.SetLocal(0,0,alpha21T.GetLocalRealPart(0,0));
            alpha21T.SetLocal(0,0,F(1));
        }

        // Form half of the left-reflector using an implicitly-updated A22:
        // y12 := tauQ a21^H ( A22 - A20 Y02 - X2L conj(AT2) )
        //      = tauQ ( a21^H A22 - (a21^H A20) Y02 - (a21^H X2L) conj(AT2) )
        //      = tauQ ( A22^H a21 - Y02^H (A20^H a21) - AT2^T (X2L^H a21) )^H
        // -------------------------------------------------------------------
        a21_MC_STAR = a21;

        // z21[MR,* ] := A22^H[MR,MC] a21[MC,* ]
        z21_MR_STAR.AlignWith( A22 );
        Zeros( z21_MR_STAR, A22.Width(), 1 );
        LocalGemv( ADJOINT, F(1), A22, a21_MC_STAR, F(0), z21_MR_STAR );

        // z01[MC,* ] := (A20^H a21)[MC,* ]
        z01_MR_STAR.AlignWith( A20 );
        Zeros( z01_MR_STAR, A20.Width(), 1 );
        LocalGemv( ADJOINT, F(1), A20, a21_MC_STAR, F(0), z01_MR_STAR );
        z01_MR_STAR.SumOver( A20.ColComm() );
        z01_MC_STAR.AlignWith( Y02 );
        z01_MC_STAR = z01_MR_STAR;
        // z21[MR,* ] -= Y02^H[MR,MC] z01[MC,* ]
        LocalGemv( ADJOINT, F(-1), Y02, z01_MC_STAR, F(1), z21_MR_STAR );

        // z01[MR,* ] := X2L^H[MR,MC] a21[MC,* ]
        zT1_MR_STAR.AlignWith( X2L );
        Zeros( zT1_MR_STAR, X2L.Width(), 1 );
        LocalGemv( ADJOINT, F(1), X2L, a21_MC_STAR, F(0), zT1_MR_STAR );
        zT1_MR_STAR.SumOver( X2L.ColComm() );
        zT1_MC_STAR.AlignWith( AT2 );
        zT1_MC_STAR = zT1_MR_STAR; 
        // z21[MR,* ] -= AT2^T[MR,MC] (X2L^H a21)[MC,* ]
        LocalGemv( TRANSPOSE, F(-1), AT2, zT1_MC_STAR, F(1), z21_MR_STAR );

        // Finally perform the column summation and then scale by tauQ
        y12.AdjointColSumScatterFrom( z21_MR_STAR );
        Scale( tauQ, y12 );
    }

    // Put back d and e
    auto ATL = View( A, 0, 0, nX, nX );
    auto ATLExpanded = View( A, 0, 0, nX+1, nX );
    ATL.SetRealPartOfDiagonal( d, 0 );
    ATLExpanded.SetRealPartOfDiagonal( e, -1 );
}

} // namespace bidiag
} // namespace elem

#endif // ifndef ELEM_BIDIAG_LPAN_HPP
