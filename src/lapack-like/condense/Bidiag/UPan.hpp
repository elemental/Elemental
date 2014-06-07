/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BIDIAG_UPAN_HPP
#define EL_BIDIAG_UPAN_HPP

namespace El {
namespace bidiag {

template<typename F> 
inline void
UPan( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ, Matrix<F>& X, Matrix<F>& Y )
{
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int nX = X.Width();
    DEBUG_ONLY(
        CallStackEntry cse("bidiag::UPan");
        if( tP.Height() != nX || tP.Width() != 1 )
            LogicError("tP was not the right size");
        if( tQ.Height() != nX || tQ.Width() != 1 )
            LogicError("tQ was not the right size");
        if( mA < nA )
            LogicError("A must be at least as tall as it is wide");
        if( mA != X.Height() )
            LogicError("A and X must have the same height");
        if( nA != Y.Width() )
            LogicError("A and Y must have the same width");
        if( X.Height() < nX )
            LogicError("X must be a column panel");
        if( Y.Height() != nX )
            LogicError("Y is the wrong height");
    )
    typedef Base<F> Real;

    Matrix<F> zT1, z01, z21;

    Matrix<Real> d, e;
    d.Resize( nX, 1 );
    e.Resize( nX, 1 );

    for( Int k=0; k<nX; ++k )
    {
        auto a01      = ViewRange( A, 0,   k,   k,   k+1 );
        auto A02      = ViewRange( A, 0,   k+1, k,   nA  );
        auto a10      = ViewRange( A, k,   0,   k+1, k   );
        auto alpha11  = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12      = ViewRange( A, k,   k+1, k+1, nA  );
        auto alpha12L = ViewRange( A, k,   k+1, k+1, k+2 );
        auto a12R     = ViewRange( A, k,   k+2, k+1, nA  );
        auto a21      = ViewRange( A, k+1, k,   mA,  k+1 );
        auto A22      = ViewRange( A, k+1, k+1, mA,  nA  );
        auto AB0      = ViewRange( A, k,   0,   mA,  k   );
        auto aB1      = ViewRange( A, k,   k,   mA,  k+1 );
        auto AB2      = ViewRange( A, k,   k+1, mA,  nA  );
        auto A2L      = ViewRange( A, k+1, 0,   mA,  k+1 );

        auto x10 = ViewRange( X, k,   0,   k+1, k   );
        auto X20 = ViewRange( X, k+1, 0,   mA,  k   );
        auto x21 = ViewRange( X, k+1, k,   mA,  k+1 );
        auto XB0 = ViewRange( X, k,   0,   mA,  k   );

        auto y01 = ViewRange( Y, 0, k,   k,   k+1 );
        auto Y02 = ViewRange( Y, 0, k+1, k,   nA  );
        auto y12 = ViewRange( Y, k, k+1, k+1, nA  );
        auto YT2 = ViewRange( Y, 0, k+1, k+1, nA  );

        // Apply all previous reflectors to aB1:
        //   aB1 := aB1 - AB0 y01 - XB0 conj(a01)
        Gemv( NORMAL, F(-1), AB0, y01, F(1), aB1 );
        Conjugate( a01 );
        Gemv( NORMAL, F(-1), XB0, a01, F(1), aB1 );
        Conjugate( a01 );

        // Find tauQ and u such that
        //  / I - tauQ | 1 | | 1, u^H | \ | alpha11 | = | delta |
        //  \          | u |            / |     a21 |   |    0  |
        const F tauQ = LeftReflector( alpha11, a21 );
        tQ.Set(k,0,tauQ);

        // Temporarily set aB1 = | 1 |
        //                       | u |
        d.Set(k,0,alpha11.GetRealPart(0,0));
        alpha11.Set(0,0,F(1));

        // Form half of the left-reflector using an implicitly-updated AB2:
        // y12 := tauQ aB1^H ( AB2 - AB0 Y02 - XB0 conj(A02) )
        //      = tauQ ( aB1^H AB2 - (aB1^H AB0) Y02 - (aB1^H XB0) conj(A02) )
        //      = tauQ ( AB2^H aB1 - Y02^H (AB0^H aB1) - A02^T (XB0^H aB1) )^H
        // -------------------------------------------------------------------
        // z21 := AB2^H aB1
        Gemv( ADJOINT, F(1), AB2, aB1, z21 );
        // z21 := z21 - Y02^H (AB0^H aB1)
        Gemv( ADJOINT, F(1),  AB0, aB1, z01 );
        Gemv( ADJOINT, F(-1), Y02, z01, F(1), z21 );
        // z21 := z21 - A02^T (XB0^H aB1)
        Gemv( ADJOINT, F(1),  XB0, aB1, z01 );
        Gemv( TRANSPOSE, F(-1), A02, z01, F(1), z21 );
        // y12 := tauQ z21^H
        Adjoint( z21, y12 );
        Scale( tauQ, y12 ); 

        // Apply all previous reflectors to a12:
        // a12 := a12 - a1L yT2           - x10 conj(A02)
        //      = a12 - (a10 Y02 + 1*y12) - x10 conj(A02)
        Gemv( TRANSPOSE, F(-1), Y02, a10, F(1), a12 );
        Axpy( F(-1), y12, a12 );
        Gemv( ADJOINT, F(-1), A02, x10, F(1), a12 ); 

        // Find tauP and v such that
        //  |alpha12L a12R| /I - tauP |1  | |1 conj(v)|\ = |epsilon 0|
        //                  \         |v^T|            /
        const F tauP = RightReflector( alpha12L, a12R );
        tP.Set(k,0,tauP);

        // Temporarily set a12 = | 1 v |
        e.Set(k,0,alpha12L.GetRealPart(0,0));
        alpha12L.Set(0,0,F(1));

        // Form half of the right-reflector using an implicitly-updated A22:
        // x21 := tauP (A22 - A2L YT2 - X20 conj(A02)) a12^T
        //      = tauP (A22 a12^T - A2L (YT2 a12^T) - X20 (conj(A02) a12^T))
        // -----------------------------------------------------------------
        // x21 := A22 a12^T
        Zeros( x21, A22.Height(), 1 );
        Gemv( NORMAL, F(1), A22, a12, F(0), x21 );
        // x21 := x21 - A2L (YT2 a12^T) 
        Gemv( NORMAL, F(1),  YT2, a12, zT1 );
        Gemv( NORMAL, F(-1), A2L, zT1, F(1), x21 );
        // x21 := x21 - X20 (conj(A02) a12^T)
        //      = x21 - X20 conj(A02 a12^H)
        Conjugate( a12 );
        Gemv( NORMAL, F(1),  A02, a12, z01 );
        Conjugate( a12 );
        Conjugate( z01 );
        Gemv( NORMAL, F(-1), X20, z01, F(1), x21 );
        // x21 := tauP x21
        Scale( tauP, x21 );
    }

    // Put back d and e
    auto ATL = View( A, 0, 0, nX, nX );
    auto ATLExpanded = View( A, 0, 0, nX, nX+1 );
    ATL.SetRealPartOfDiagonal( d, 0 );
    ATLExpanded.SetRealPartOfDiagonal( e, 1 );
}

template<typename F> 
inline void
UPan
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
        CallStackEntry cse("bidiag::UPan");
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
        if( mA < nA )
            LogicError("A must be at least as tall as it is wide");
        if( mA != X.Height() )
            LogicError("A and X must have the same height");
        if( nA != Y.Width() )
            LogicError("A and Y must have the same width");
        if( X.Height() < nX )
            LogicError("X must be a column panel");
        if( Y.Height() != nX )
            LogicError("Y is the wrong height");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();

    DistMatrix<F,MC,  STAR> z01_MC_STAR(g),
                            zT1_MC_STAR(g),
                            z21_MC_STAR(g),
                            zB1_MC_STAR(g);
    DistMatrix<F,MR,  STAR> z01_MR_STAR(g),
                            zT1_MR_STAR(g),
                            z21_MR_STAR(g),
                            y01_MR_STAR(g),
                            a01_MR_STAR(g);
    DistMatrix<F,STAR,MC  > x10_STAR_MC(g),
                            a10_STAR_MC(g);

    DistMatrix<Real,MD,STAR> d(g), e(g);
    d.SetRoot( A.DiagonalRoot(0) );
    e.SetRoot( A.DiagonalRoot(1) );
    d.AlignCols( A.DiagonalAlign(0) );
    e.AlignCols( A.DiagonalAlign(1) ); 
    d.Resize( nX, 1 );
    e.Resize( nX, 1 );

    for( Int k=0; k<nX; ++k )
    {
        auto a01      = ViewRange( A, 0,   k,   k,   k+1 );
        auto A02      = ViewRange( A, 0,   k+1, k,   nA  );
        auto a10      = ViewRange( A, k,   0,   k+1, k   );
        auto alpha11  = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12      = ViewRange( A, k,   k+1, k+1, nA  );
        auto alpha12L = ViewRange( A, k,   k+1, k+1, k+2 );
        auto a12R     = ViewRange( A, k,   k+2, k+1, nA  );
        auto a21      = ViewRange( A, k+1, k,   mA,  k+1 );
        auto A22      = ViewRange( A, k+1, k+1, mA,  nA  );
        auto AB0      = ViewRange( A, k,   0,   mA,  k   );
        auto aB1      = ViewRange( A, k,   k,   mA,  k+1 );
        auto AB2      = ViewRange( A, k,   k+1, mA,  nA  );
        auto A2L      = ViewRange( A, k+1, 0,   mA,  k+1 );

        auto a12_STAR_MR = ViewRange( AT_STAR_MR, k, k+1, k+1, nA  );
        auto aB1_MC_STAR = ViewRange( AL_MC_STAR, k, k,   mA,  k+1 );

        auto x10 = ViewRange( X, k,   0,   k+1, k   );
        auto X20 = ViewRange( X, k+1, 0,   mA,  k   );
        auto x21 = ViewRange( X, k+1, k,   mA,  k+1 );
        auto XB0 = ViewRange( X, k,   0,   mA,  k   );

        auto y01 = ViewRange( Y, 0, k,   k,   k+1 );
        auto Y02 = ViewRange( Y, 0, k+1, k,   nA  );
        auto y12 = ViewRange( Y, k, k+1, k+1, nA  );
        auto YT2 = ViewRange( Y, 0, k+1, k+1, nA  );

        auto delta1   = View( d,  k, 0, 1, 1 );
        auto epsilon1 = View( e,  k, 0, 1, 1 );

        // Apply all previous reflectors to aB1:
        //   aB1 := aB1 - AB0 y01 - XB0 conj(a01)
        if( k > 0 )
        {
            y01_MR_STAR.AlignWith( AB0 );
            a01_MR_STAR.AlignWith( AB0 );
            y01_MR_STAR = y01;
            a01_MR_STAR = a01;

            zB1_MC_STAR.AlignWith( aB1 );
            Zeros( zB1_MC_STAR, aB1.Height(), 1 );
            // zB1[MC,* ] := AB0[MC,MR] y01[MR,* ] + XB0[MC,MR] conj(a01[MR,* ])
            LocalGemv( NORMAL, F(1), AB0, y01_MR_STAR, F(0), zB1_MC_STAR );
            Conjugate( a01_MR_STAR );
            LocalGemv( NORMAL, F(1), XB0, a01_MR_STAR, F(1), zB1_MC_STAR );
            // Sum the partial contributions and subtract from aB1
            aB1.RowSumScatterUpdate( F(-1), zB1_MC_STAR );
        }

        // Find tauQ and u such that
        //  / I - tauQ | 1 | | 1, u^H | \ | alpha11 | = | delta |
        //  \          | u |            / |     a21 |   |    0  |
        const F tauQ = LeftReflector( alpha11, a21 );
        tQ.Set(k,0,tauQ);

        // Temporarily set aB1 = | 1 |
        //                       | u |
        if( delta1.IsLocal(0,0) )
        {
            delta1.SetLocal(0,0,alpha11.GetLocalRealPart(0,0));
            alpha11.SetLocal(0,0,F(1));
        }

        // Form half of the left-reflector using an implicitly-updated AB2:
        // y12 := tauQ aB1^H ( AB2 - AB0 Y02 - XB0 conj(A02) )
        //      = tauQ ( aB1^H AB2 - (aB1^H AB0) Y02 - (aB1^H XB0) conj(A02) )
        //      = tauQ ( AB2^H aB1 - Y02^H (AB0^H aB1) - A02^T (XB0^H aB1) )^H
        // -------------------------------------------------------------------
        aB1_MC_STAR = aB1;

        // z21[MR,* ] := AB2^H[MR,MC] aB1[MC,* ]
        z21_MR_STAR.AlignWith( AB2 );
        Zeros( z21_MR_STAR, A22.Width(), 1 );
        LocalGemv( ADJOINT, F(1), AB2, aB1_MC_STAR, F(0), z21_MR_STAR );

        // z01[MC,* ] := (AB0^H aB1)[MC,* ]
        z01_MR_STAR.AlignWith( AB0 );
        Zeros( z01_MR_STAR, AB0.Width(), 1 );
        LocalGemv( ADJOINT, F(1), AB0, aB1_MC_STAR, F(0), z01_MR_STAR );
        z01_MR_STAR.SumOver( AB0.ColComm() );
        z01_MC_STAR.AlignWith( Y02 );
        z01_MC_STAR = z01_MR_STAR;
        // z21[MR,* ] -= Y02^H[MR,MC] z01[MC,* ] 
        LocalGemv( ADJOINT, F(-1), Y02, z01_MC_STAR, F(1), z21_MR_STAR );

        // z01[MR,* ] := XB0^H[MR,MC] aB1[MC,* ]
        z01_MR_STAR.AlignWith( XB0 );
        Zeros( z01_MR_STAR, XB0.Width(), 1 );
        LocalGemv( ADJOINT, F(1), XB0, aB1_MC_STAR, F(0), z01_MR_STAR );
        z01_MR_STAR.SumOver( XB0.ColComm() );
        z01_MC_STAR.AlignWith( A02 );
        z01_MC_STAR = z01_MR_STAR;
        // z21[MR,* ] -= A02^T[MR,MC] (XB0^H aB1)[MC,* ]
        LocalGemv( TRANSPOSE, F(-1), A02, z01_MC_STAR, F(1), z21_MR_STAR );

        // Finally perform the column summation and then scale by tauQ
        y12.AdjointColSumScatterFrom( z21_MR_STAR );
        Scale( tauQ, y12 );

        // Apply all previous reflectors to a12:
        // a12 := a12 - a1L yT2           - x10 conj(A02)
        //      = a12 - (a10 Y02 + 1*y12) - x10 conj(A02)
        // ----------------------------------------------
        // a12 := a12 - a10 Y02 (do not sum over columns yet)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        a10_STAR_MC.AlignWith( Y02 );
        a10_STAR_MC = a10;
        // z21[MR,* ] := Y02^T[MR,MC] a10^T[MC,* ]
        z21_MR_STAR.AlignWith( Y02 );
        Zeros( z21_MR_STAR, Y02.Width(), 1 );
        LocalGemv( TRANSPOSE, F(1), Y02, a10_STAR_MC, F(0), z21_MR_STAR );

        // a12 := a12 - x10 conj(A02) (and incorporate last update into sum)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        x10_STAR_MC.AlignWith( A02 );
        x10_STAR_MC = x10;
        // z21[MR,* ] := A02^H[MR,MC] x10^T[MC,* ]
        LocalGemv( ADJOINT, F(1), A02, x10_STAR_MC, F(1), z21_MR_STAR );
        // Sum the partial contributions from the past two updates
        a12.TransposeColSumScatterUpdate( F(-1), z21_MR_STAR );

        // a12 := a12 - y12
        // ^^^^^^^^^^^^^^^^
        Axpy( F(-1), y12, a12 );

        // Find tauP and v such that
        //  |alpha12L a12R| /I - tauP |1  | |1, conj(v)|\ = |epsilon 0|
        //                  \         |v^T|             /
        const F tauP = RightReflector( alpha12L, a12R );
        tP.Set(k,0,tauP);

        // Temporarily set a12 = | 1 v |
        if( epsilon1.IsLocal(0,0) )
        {
            epsilon1.SetLocal(0,0,alpha12L.GetLocalRealPart(0,0));
            alpha12L.SetLocal(0,0,F(1));
        }

        // Form half of the right-reflector using an implicitly-updated A22:
        // x21 := tauP (A22 - A2L YT2 - X20 conj(A02)) a12^T
        //      = tauP (A22 a12^T - A2L (YT2 a12^T) - X20 (conj(A02) a12^T))
        // -----------------------------------------------------------------
        a12_STAR_MR = a12;

        // z21[MC,* ] := A22[MC,MR] a12^T[MR,* ]
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        z21_MC_STAR.AlignWith( A22 );
        Zeros( z21_MC_STAR, A22.Height(), 1 );
        LocalGemv( NORMAL, F(1), A22, a12_STAR_MR, F(0), z21_MC_STAR );
       
        // z21[MC,* ] -= A2L[MC,MR] (YT2 a12^T)[MR,* ]
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // zT1[MC,* ] := YT2[MC,MR] a12^T[MR,* ]
        zT1_MC_STAR.AlignWith( YT2 );
        Zeros( zT1_MC_STAR, YT2.Height(), 1 );
        LocalGemv( NORMAL, F(1), YT2, a12_STAR_MR, F(0), zT1_MC_STAR );
        zT1_MC_STAR.SumOver( YT2.RowComm() );
        // Redistribute and perform local Gemv 
        zT1_MR_STAR.AlignWith( A2L );
        zT1_MR_STAR = zT1_MC_STAR;
        LocalGemv( NORMAL, F(-1), A2L, zT1_MR_STAR, F(1), z21_MC_STAR );

        // z21[MC,* ] -= X20[MC,MR] (conj(A02) a12^T)[MR,* ]
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // z01[MC,* ] := conj(A02[MC,MR] a12^H[MR,* ])
        z01_MC_STAR.AlignWith( A02 ); 
        Zeros( z01_MC_STAR, A02.Height(), 1 );
        Conjugate( a12_STAR_MR );
        LocalGemv( NORMAL, F(1), A02, a12_STAR_MR, F(0), z01_MC_STAR );
        z01_MC_STAR.SumOver( A02.RowComm() );
        Conjugate( a12_STAR_MR );
        Conjugate( z01_MC_STAR );
        // Redistribute and perform local Gemv
        z01_MR_STAR.AlignWith( X20 );
        z01_MR_STAR = z01_MC_STAR;
        LocalGemv( NORMAL, F(-1), X20, z01_MR_STAR, F(1), z21_MC_STAR );
        // Sum the various contributions within process rows
        x21.RowSumScatterFrom( z21_MC_STAR );
        Scale( tauP, x21 );
    }

    // Put back d and e
    auto ATL = View( A, 0, 0, nX, nX );
    auto ATLExpanded = View( A, 0, 0, nX, nX+1 );
    ATL.SetRealPartOfDiagonal( d, 0 );
    ATLExpanded.SetRealPartOfDiagonal( e, 1 );
}

} // namespace bidiag
} // namespace El

#endif // ifndef EL_BIDIAG_UPAN_HPP
