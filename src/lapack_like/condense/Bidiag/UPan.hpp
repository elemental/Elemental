/*
   Copyright (c) 2009-2015, Jack Poulson
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
    const Int nX = X.Width();
    DEBUG_ONLY(
      CSE cse("bidiag::UPan");
      if( tP.Height() != nX || tP.Width() != 1 )
          LogicError("tP was not the right size");
      if( tQ.Height() != nX || tQ.Width() != 1 )
          LogicError("tQ was not the right size");
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( A.Height() != X.Height() )
          LogicError("A and X must have the same height");
      if( A.Width() != Y.Width() )
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
        const Range<Int> ind0( 0, k ), ind1( k ), ind2( k+1, END ),
                         indT( 0, k+1 ), indL( 0, k+1 ),
                         indB( k, END );

        auto a01      = A( ind0, ind1 );
        auto A02      = A( ind0, ind2 );
        auto a10      = A( ind1, ind0 );
        auto alpha11  = A( ind1, ind1 );
        auto a12      = A( ind1, ind2 );
        auto a21      = A( ind2, ind1 );
        auto A22      = A( ind2, ind2 );
        auto AB0      = A( indB, ind0 );
        auto aB1      = A( indB, ind1 );
        auto AB2      = A( indB, ind2 );
        auto A2L      = A( ind2, indL );

        auto alpha12L = A( IR(k), IR(k+1)     );
        auto a12R     = A( IR(k), IR(k+2,END) );

        auto x10 = X( ind1, ind0 );
        auto X20 = X( ind2, ind0 );
        auto x21 = X( ind2, ind1 );
        auto XB0 = X( indB, ind0 );

        auto y01 = Y( ind0, ind1 );
        auto Y02 = Y( ind0, ind2 );
        auto y12 = Y( ind1, ind2 );
        auto YT2 = Y( indT, ind2 );

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
        y12 *= tauQ;

        // Apply all previous reflectors to a12:
        // a12 := a12 - a1L yT2           - x10 conj(A02)
        //      = a12 - (a10 Y02 + 1*y12) - x10 conj(A02)
        Gemv( TRANSPOSE, F(-1), Y02, a10, F(1), a12 );
        a12 -= y12;
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
        x21 *= tauP;
    }

    // Put back d and e
    auto ATL = A( IR(0,nX), IR(0,nX) );
    auto ATLExpanded = A( IR(0,nX), IR(0,nX+1) );
    SetRealPartOfDiagonal( ATL, d, 0 );
    SetRealPartOfDiagonal( ATLExpanded, e, 1 );
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
    const Int nX = X.Width();
    DEBUG_ONLY(
      CSE cse("bidiag::UPan");
      AssertSameGrids( A, tP, tQ, X, Y, AL_MC_STAR, AT_STAR_MR );
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
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( A.Height() != X.Height() )
          LogicError("A and X must have the same height");
      if( A.Width() != Y.Width() )
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
        const Range<Int> ind0( 0, k ), ind1( k ), ind2( k+1, END ),
                         indT( 0, k+1 ), indL( 0, k+1 ), indB( k, END );

        auto a01      = A( ind0, ind1 );
        auto A02      = A( ind0, ind2 );
        auto a10      = A( ind1, ind0 );
        auto alpha11  = A( ind1, ind1 );
        auto a12      = A( ind1, ind2 );
        auto a21      = A( ind2, ind1 );
        auto A22      = A( ind2, ind2 );
        auto AB0      = A( indB, ind0 );
        auto aB1      = A( indB, ind1 );
        auto AB2      = A( indB, ind2 );
        auto A2L      = A( ind2, indL );

        auto alpha12L = A( ind1, IR(k+1)     );
        auto a12R     = A( ind1, IR(k+2,END) );

        auto x10 = X( ind1, ind0 );
        auto X20 = X( ind2, ind0 );
        auto x21 = X( ind2, ind1 );
        auto XB0 = X( indB, ind0 );

        auto y01 = Y( ind0, ind1 );
        auto Y02 = Y( ind0, ind2 );
        auto y12 = Y( ind1, ind2 );
        auto YT2 = Y( indT, ind2 );

        auto a12_STAR_MR = AT_STAR_MR( ind1, ind2 );
        auto aB1_MC_STAR = AL_MC_STAR( indB, ind1 );

        auto delta1   = d( ind1, ALL );
        auto epsilon1 = e( ind1, ALL );

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
            AxpyContract( F(-1), zB1_MC_STAR, aB1 );
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
        El::AllReduce( z01_MR_STAR, AB0.ColComm() );
        z01_MC_STAR.AlignWith( Y02 );
        z01_MC_STAR = z01_MR_STAR;
        // z21[MR,* ] -= Y02^H[MR,MC] z01[MC,* ] 
        LocalGemv( ADJOINT, F(-1), Y02, z01_MC_STAR, F(1), z21_MR_STAR );

        // z01[MR,* ] := XB0^H[MR,MC] aB1[MC,* ]
        z01_MR_STAR.AlignWith( XB0 );
        Zeros( z01_MR_STAR, XB0.Width(), 1 );
        LocalGemv( ADJOINT, F(1), XB0, aB1_MC_STAR, F(0), z01_MR_STAR );
        El::AllReduce( z01_MR_STAR, XB0.ColComm() );
        z01_MC_STAR.AlignWith( A02 );
        z01_MC_STAR = z01_MR_STAR;
        // z21[MR,* ] -= A02^T[MR,MC] (XB0^H aB1)[MC,* ]
        LocalGemv( TRANSPOSE, F(-1), A02, z01_MC_STAR, F(1), z21_MR_STAR );

        // Finally perform the column summation and then scale by tauQ
        AdjointContract( z21_MR_STAR, y12 );
        y12 *= tauQ;

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
        TransposeAxpyContract( F(-1), z21_MR_STAR, a12 );

        // a12 := a12 - y12
        // ^^^^^^^^^^^^^^^^
        a12 -= y12;

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
        El::AllReduce( zT1_MC_STAR, YT2.RowComm() );
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
        El::AllReduce( z01_MC_STAR, A02.RowComm() );
        Conjugate( a12_STAR_MR );
        Conjugate( z01_MC_STAR );
        // Redistribute and perform local Gemv
        z01_MR_STAR.AlignWith( X20 );
        z01_MR_STAR = z01_MC_STAR;
        LocalGemv( NORMAL, F(-1), X20, z01_MR_STAR, F(1), z21_MC_STAR );
        // Sum the various contributions within process rows
        Contract( z21_MC_STAR, x21 );
        x21 *= tauP;
    }

    // Put back d and e
    auto ATL = A( IR(0,nX), IR(0,nX) );
    auto ATLExpanded = A( IR(0,nX), IR(0,nX+1) );
    SetRealPartOfDiagonal( ATL, d, 0 );
    SetRealPartOfDiagonal( ATLExpanded, e, 1 );
}

} // namespace bidiag
} // namespace El

#endif // ifndef EL_BIDIAG_UPAN_HPP
