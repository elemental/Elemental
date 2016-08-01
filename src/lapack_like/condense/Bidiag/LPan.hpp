/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BIDIAG_LPAN_HPP
#define EL_BIDIAG_LPAN_HPP

namespace El {
namespace bidiag {

template<typename F>
void
LPan
( Matrix<F>& A,
  Matrix<F>& phaseP,
  Matrix<F>& phaseQ,
  Matrix<F>& X,
  Matrix<F>& Y )
{
    DEBUG_CSE
    const Int nX = X.Width();
    DEBUG_ONLY(
      if( phaseP.Height() != nX || phaseP.Width() != 1 )
          LogicError("phaseP was not the right size");
      if( phaseQ.Height() != nX || phaseQ.Width() != 1 )
          LogicError("phaseQ was not the right size");
      if( A.Height() > A.Width() )
          LogicError("A must be at least as wide as it is tall"); 
      if( A.Height() != X.Height() )
          LogicError("A and X must have the same height");
      if( A.Width() != Y.Width() )
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
        const Range<Int> ind0( 0, k ), ind1( k ), ind2( k+1, END ),
                         indT( 0, k+1 ), indL( 0, k+1 ),
                         indB( k, END ), indR( k, END );

        auto a01      = A( ind0, ind1 );
        auto AT2      = A( indT, ind2 );
        auto A0R      = A( ind0, indR );
        auto a10      = A( ind1, ind0 );
        auto alpha11  = A( ind1, ind1 );
        auto a1R      = A( ind1, indR );
        auto a12      = A( ind1, ind2 );
        auto a21      = A( ind2, ind1 );
        auto A20      = A( ind2, ind0 );
        auto A2R      = A( ind2, indR );
        auto A22      = A( ind2, ind2 );

        auto alpha21T = A( IR(k+1),     ind1 );
        auto a21B     = A( IR(k+2,END), ind1 );

        auto x10 = X( ind1, ind0 );
        auto XB0 = X( indB, ind0 );
        auto X20 = X( ind2, ind0 );
        auto X2L = X( ind2, indL );
        auto x21 = X( ind2, ind1 );

        auto y01 = Y( ind0, ind1 );
        auto Y0R = Y( ind0, indR );
        auto Y02 = Y( ind0, ind2 );
        auto y12 = Y( ind1, ind2 );

        // Apply all previous reflectors to a1R:    
        //   a1R := a1R - a10 Y0R         - x10 conj(A0R)
        //        = a1R - (Y0R^T a10^T)^T - (A0R^H x10^T)^T
        Gemv( TRANSPOSE, F(-1), Y0R, a10, F(1), a1R );
        Gemv( ADJOINT,   F(-1), A0R, x10, F(1), a1R );

        // Find tauP and v such that
        //  |alpha11 a12| /I - tauP |1  | |1 conj(v)|\ = |delta 0|
        //                \         |v^T|            /
        const F tauP = RightReflector( alpha11, a12 );
        phaseP(k) = tauP;

        // Temporarily set a1R = | 1 v |
        d(k) = RealPart(alpha11(0));
        alpha11(0) = F(1);

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
        x21 *= tauP;

        // Apply all previous reflectors to a21:
        //   a21 := a21 - A20 y01 - X2L conj(aT1)
        //        = a21 - A20 y01 - (X20 conj(a01) + x21*1)
        Gemv( NORMAL, F(-1), A20, y01, F(1), a21 );
        Conjugate( a01 );
        Gemv( NORMAL, F(-1), X20, a01, F(1), a21 );
        Conjugate( a01 );
        a12 -= x21;

        // Find tauQ and u such that
        //  / I - tauQ | 1 | | 1, u^H | \ | alpha21T | = | epsilon |
        //  \          | u |            / |     a21B |   |    0    |
        const F tauQ = LeftReflector( alpha21T, a21B );
        phaseQ(k) = tauQ;

        // Temporarily set a21 = | 1 |
        //                       | u |
        e(k) = RealPart(alpha21T(0));
        alpha21T(0) = F(1);

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
        y12 *= tauQ;
    }

    // Put back d and e
    auto ATL = A( IR(0,nX), IR(0,nX) );
    auto ATLExpanded = A( IR(0,nX+1), IR(0,nX) );
    SetRealPartOfDiagonal( ATL, d, 0 );
    SetRealPartOfDiagonal( ATLExpanded, e, -1 );
}

template<typename F>
void
LPan
( DistMatrix<F>& A,
  DistMatrix<F,STAR,STAR>& phaseP,
  DistMatrix<F,STAR,STAR>& phaseQ,
  DistMatrix<F>& X,
  DistMatrix<F>& Y,
  DistMatrix<F,MC,  STAR>& AL_MC_STAR,
  DistMatrix<F,STAR,MR  >& AT_STAR_MR )
{
    DEBUG_CSE
    const Int nX = X.Width();
    DEBUG_ONLY(
      AssertSameGrids( A, phaseP, phaseQ, X, Y, AL_MC_STAR, AT_STAR_MR );
      if( A.ColAlign() != X.ColAlign() ||
          A.RowAlign() != X.RowAlign() )
          LogicError("A and X must be aligned");
      if( A.ColAlign() != Y.ColAlign() ||
          A.RowAlign() != Y.RowAlign() )
          LogicError("A and Y must be aligned");
      if( phaseP.Height() != nX || phaseP.Width() != 1 )
          LogicError("phaseP was not the right size");
      if( phaseQ.Height() != nX || phaseQ.Width() != 1 )
          LogicError("phaseQ was not the right size");
      if( A.Height() > A.Width() )
          LogicError("A must be at least as wide as it is tall"); 
      if( A.Height() != X.Height() )
          LogicError("A and X must have the same height");
      if( A.Width() != Y.Width() )
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
        const Range<Int> ind0( 0, k ), ind1( k ), ind2( k+1, END ),
                         indT( 0, k+1 ), indL( 0, k+1 ),
                         indB( k, END ), indR( k, END );

        auto a01      = A( ind0, ind1 );
        auto AT2      = A( indT, ind2 );
        auto A0R      = A( ind0, indR );
        auto a10      = A( ind1, ind0 );
        auto alpha11  = A( ind1, ind1 );
        auto a1R      = A( ind1, indR );
        auto a12      = A( ind1, ind2 );
        auto a21      = A( ind2, ind1 );
        auto A20      = A( ind2, ind0 );
        auto A2R      = A( ind2, indR );
        auto A22      = A( ind2, ind2 );

        auto alpha21T = A( IR(k+1),     ind1 );
        auto a21B     = A( IR(k+2,END), ind1 );

        auto x10 = X( ind1, ind0 );
        auto XB0 = X( indB, ind0 );
        auto X20 = X( ind2, ind0 );
        auto X2L = X( ind2, indL );
        auto x21 = X( ind2, ind1 );

        auto y01 = Y( ind0, ind1 );
        auto Y0R = Y( ind0, indR );
        auto Y02 = Y( ind0, ind2 );
        auto y12 = Y( ind1, ind2 );

        auto a1R_STAR_MR = AT_STAR_MR( ind1, indR );
        auto a21_MC_STAR = AL_MC_STAR( ind2, ind1 ); 

        auto delta1   = d( ind1, ALL );
        auto epsilon1 = e( ind1, ALL );

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
            AxpyContract( F(-1), z1R_STAR_MR, a1R );
        }

        // Find tauP and v such that
        //  |alpha11 a12| /I - tauP |1  | |1 conj(v)|\ = |delta 0|
        //                \         |v^T|            /
        const F tauP = RightReflector( alpha11, a12 );
        phaseP.Set(k,0,tauP);

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
        El::AllReduce( z01_MC_STAR, Y0R.RowComm() );
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
        El::AllReduce( z01_MC_STAR, A0R.RowComm() );
        Conjugate( a1R_STAR_MR );
        Conjugate( z01_MC_STAR );
        z01_MR_STAR.AlignWith( X20 );
        z01_MR_STAR = z01_MC_STAR;
        // z21[MC,* ] -= X20[MC,MR] z01[MR,* ] 
        LocalGemv( NORMAL, F(-1), X20, z01_MR_STAR, F(1), z21_MC_STAR );

        // Finally perform the row summation and then scale by tauP
        Contract( z21_MC_STAR, x21 );
        x21 *= tauP;

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
        AxpyContract( F(-1), z21_MC_STAR, a21 );

        // a21 := a21 - x21
        // ^^^^^^^^^^^^^^^^
        a21 -= x21;

        // Find tauQ and u such that
        //  / I - tauQ | 1 | | 1, u^H | \ | alpha21T | = | epsilon |
        //  \          | u |            / |     a21B |   |    0    |
        const F tauQ = LeftReflector( alpha21T, a21B );
        phaseQ.Set(k,0,tauQ);

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
        El::AllReduce( z01_MR_STAR, A20.ColComm() );
        z01_MC_STAR.AlignWith( Y02 );
        z01_MC_STAR = z01_MR_STAR;
        // z21[MR,* ] -= Y02^H[MR,MC] z01[MC,* ]
        LocalGemv( ADJOINT, F(-1), Y02, z01_MC_STAR, F(1), z21_MR_STAR );

        // z01[MR,* ] := X2L^H[MR,MC] a21[MC,* ]
        zT1_MR_STAR.AlignWith( X2L );
        Zeros( zT1_MR_STAR, X2L.Width(), 1 );
        LocalGemv( ADJOINT, F(1), X2L, a21_MC_STAR, F(0), zT1_MR_STAR );
        El::AllReduce( zT1_MR_STAR, X2L.ColComm() );
        zT1_MC_STAR.AlignWith( AT2 );
        zT1_MC_STAR = zT1_MR_STAR; 
        // z21[MR,* ] -= AT2^T[MR,MC] (X2L^H a21)[MC,* ]
        LocalGemv( TRANSPOSE, F(-1), AT2, zT1_MC_STAR, F(1), z21_MR_STAR );

        // Finally perform the column summation and then scale by tauQ
        AdjointContract( z21_MR_STAR, y12 );
        y12 *= tauQ;
    }

    // Put back d and e
    auto ATL = A( IR(0,nX), IR(0,nX) );
    auto ATLExpanded = A( IR(0,nX+1), IR(0,nX) );
    SetRealPartOfDiagonal( ATL, d, 0 );
    SetRealPartOfDiagonal( ATLExpanded, e, -1 );
}

} // namespace bidiag
} // namespace El

#endif // ifndef EL_BIDIAG_LPAN_HPP
