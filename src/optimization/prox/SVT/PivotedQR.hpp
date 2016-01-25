/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SVT_PIVOTEDQR_HPP
#define EL_SVT_PIVOTEDQR_HPP

namespace El {
namespace svt {

// Preprocess with numSteps iterations of pivoted QR factorization

template<typename F>
Int PivotedQR( Matrix<F>& A, Base<F> tau, Int numSteps, bool relative )
{
    DEBUG_ONLY(
      CSE cse("svt::PivotedQR");
      if( numSteps > Min(A.Height(),A.Width()) )
          LogicError("number of steps is too large");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<F> ACopy( A ), t;
    Matrix<Real> d;
    Permutation Omega;
    QRCtrl<Base<F>> qrCtrl;
    qrCtrl.boundRank = true;
    qrCtrl.maxRank = numSteps;
    QR( ACopy, t, d, Omega, qrCtrl );
    auto ACopyUpper = ACopy( IR(0,numSteps), IR(0,n) );

    Matrix<F> R( ACopyUpper );
    MakeTrapezoidal( UPPER, R );

    Matrix<F> U, V;
    Matrix<Real> s;
    SVDCtrl<Real> svdCtrl;
    svdCtrl.approach = PRODUCT_SVD;
    svdCtrl.tol = tau;
    svdCtrl.relative = relative;
    SVD( R, U, s, V, svdCtrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Omega.PermuteRows( V );
    Matrix<F> RThresh;
    Gemm( NORMAL, ADJOINT, F(1), U, V, RThresh );

    ACopy.Resize( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, ACopy, t );
    DiagonalScale( RIGHT, NORMAL, d, ACopy );
    Gemm( NORMAL, NORMAL, F(1), ACopy, RThresh, F(0), A );

    return ZeroNorm( s );
}

template<typename F>
Int PivotedQR
( ElementalMatrix<F>& APre, Base<F> tau, Int numSteps, bool relative )
{
    DEBUG_ONLY(
      CSE cse("svt::PivotedQR");
      if( numSteps > Min(APre.Height(),APre.Width()) )
          LogicError("number of steps is too large");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<F> ACopy( A );
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Real,MD,STAR> d(g);
    DistPermutation Omega(g);
    QRCtrl<Base<F>> qrCtrl;
    qrCtrl.boundRank = true;
    qrCtrl.maxRank = numSteps;
    QR( ACopy, t, d, Omega, qrCtrl );
    auto ACopyUpper = ACopy( IR(0,numSteps), IR(0,n) );

    DistMatrix<F> R( ACopyUpper );
    MakeTrapezoidal( UPPER, R );

    DistMatrix<F> U(g), V(g);
    DistMatrix<Real,VR,STAR> s(g);
    SVDCtrl<Real> svdCtrl;
    svdCtrl.approach = PRODUCT_SVD;
    svdCtrl.tol = tau;
    svdCtrl.relative = relative;
    SVD( R, U, s, V, svdCtrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Omega.PermuteRows( V );
    DistMatrix<F> RThresh(g);
    Gemm( NORMAL, ADJOINT, F(1), U, V, RThresh );

    ACopy.Resize( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, ACopy, t );
    DiagonalScale( RIGHT, NORMAL, d, ACopy );
    Gemm( NORMAL, NORMAL, F(1), ACopy, RThresh, F(0), A );

    return ZeroNorm( s );
}

} // namespace svt
} // namespace El

#endif // ifndef EL_SVT_PIVOTEDQR_HPP
