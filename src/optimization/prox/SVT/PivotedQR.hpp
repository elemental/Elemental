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

template<typename Field>
Int PivotedQR
( Matrix<Field>& A, const Base<Field>& tau, Int numSteps, bool relative )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( numSteps > Min(A.Height(),A.Width()) )
          LogicError("number of steps is too large");
    )
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Field> ACopy( A ), t;
    Matrix<Real> d;
    Permutation Omega;
    QRCtrl<Real> qrCtrl;
    qrCtrl.boundRank = true;
    qrCtrl.maxRank = numSteps;
    QR( ACopy, t, d, Omega, qrCtrl );
    auto ACopyUpper = ACopy( IR(0,numSteps), IR(0,n) );

    Matrix<Field> R( ACopyUpper );
    MakeTrapezoidal( UPPER, R );

    Matrix<Field> U, V;
    Matrix<Real> s;
    SVDCtrl<Real> svdCtrl;
    svdCtrl.bidiagSVDCtrl.approach = PRODUCT_SVD;
    // TODO(poulson): A more aggressive tolerance? tol = tau is too string
    SVD( R, U, s, V, svdCtrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Omega.PermuteRows( V );
    Matrix<Field> RThresh;
    Gemm( NORMAL, ADJOINT, Field(1), U, V, RThresh );

    ACopy.Resize( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, ACopy, t );
    DiagonalScale( RIGHT, NORMAL, d, ACopy );
    Gemm( NORMAL, NORMAL, Field(1), ACopy, RThresh, Field(0), A );

    return ZeroNorm( s );
}

template<typename Field>
Int PivotedQR
( AbstractDistMatrix<Field>& APre, const Base<Field>& tau, Int numSteps,
  bool relative )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( numSteps > Min(APre.Height(),APre.Width()) )
          LogicError("number of steps is too large");
    )
    typedef Base<Field> Real;

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<Field> ACopy( A );
    DistMatrix<Field,MD,STAR> t(g);
    DistMatrix<Real,MD,STAR> d(g);
    DistPermutation Omega(g);
    QRCtrl<Real> qrCtrl;
    qrCtrl.boundRank = true;
    qrCtrl.maxRank = numSteps;
    QR( ACopy, t, d, Omega, qrCtrl );
    auto ACopyUpper = ACopy( IR(0,numSteps), IR(0,n) );

    DistMatrix<Field> R( ACopyUpper );
    MakeTrapezoidal( UPPER, R );

    DistMatrix<Field> U(g), V(g);
    DistMatrix<Real,VR,STAR> s(g);
    SVDCtrl<Real> svdCtrl;
    svdCtrl.bidiagSVDCtrl.approach = PRODUCT_SVD;
    // TODO(poulson): A more aggressive tolerance? tol = tau is too string
    SVD( R, U, s, V, svdCtrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Omega.PermuteRows( V );
    DistMatrix<Field> RThresh(g);
    Gemm( NORMAL, ADJOINT, Field(1), U, V, RThresh );

    ACopy.Resize( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, ACopy, t );
    DiagonalScale( RIGHT, NORMAL, d, ACopy );
    Gemm( NORMAL, NORMAL, Field(1), ACopy, RThresh, Field(0), A );

    return ZeroNorm( s );
}

} // namespace svt
} // namespace El

#endif // ifndef EL_SVT_PIVOTEDQR_HPP
