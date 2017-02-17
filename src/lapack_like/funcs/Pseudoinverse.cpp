/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

// Replace A with its pseudoinverse

template<typename Field>
void Pseudoinverse( Matrix<Field>& A, Base<Field> tolerance )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Real eps = limits::Epsilon<Real>();

    // Get the SVD of A
    Matrix<Real> s;
    Matrix<Field> U, V;
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    ctrl.bidiagSVDCtrl.approach = COMPACT_SVD;
    // TODO(poulson): Let the user change these defaults
    ctrl.bidiagSVDCtrl.tolType = RELATIVE_TO_MAX_SING_VAL_TOL;
    ctrl.bidiagSVDCtrl.tol =
      ( tolerance == Real(0) ? Max(m,n)*eps : tolerance );
    SVD( A, U, s, V, ctrl );

    // Scale U with the inverted (nonzero) singular values, U := U / Sigma
    DiagonalSolve( RIGHT, NORMAL, s, U );

    // Form pinvA = (U Sigma V^H)^H = V (U Sigma)^H
    Gemm( NORMAL, ADJOINT, Field(1), V, U, A );
}

template<typename Field>
void HermitianPseudoinverse
( UpperOrLower uplo, Matrix<Field>& A, Base<Field> tolerance )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    // Get the EVD of A
    // TODO(poulson): Use a relative eigenvalue lower bound
    Matrix<Real> w;
    Matrix<Field> Z;
    HermitianEig( uplo, A, w, Z );

    if( tolerance == Real(0) )
    {
        // Set the tolerance equal to n ||A||_2 eps
        const Int n = Z.Height();
        const Real eps = limits::Epsilon<Real>();
        const Real twoNorm = MaxNorm( w );
        tolerance = n*twoNorm*eps;
    }
    // Invert above the tolerance
    auto omegaMap =
      [=]( const Real& omega )
      { return ( omega < tolerance ? Real(0) : 1/omega ); };
    EntrywiseMap( w, MakeFunction(omegaMap) );

    // Form the pseudoinverse
    HermitianFromEVD( uplo, A, w, Z );
}

template<typename Field>
void Pseudoinverse( AbstractDistMatrix<Field>& APre, Base<Field> tolerance )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();

    const Real eps = limits::Epsilon<Real>();

    // Get the SVD of A
    DistMatrix<Real,VR,STAR> s(g);
    DistMatrix<Field> U(g), V(g);
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    ctrl.bidiagSVDCtrl.approach = COMPACT_SVD;
    // TODO(poulson): Let the user change these defaults
    ctrl.bidiagSVDCtrl.tolType = RELATIVE_TO_MAX_SING_VAL_TOL;
    ctrl.bidiagSVDCtrl.tol =
      ( tolerance == Real(0) ? Max(m,n)*eps : tolerance );
    SVD( A, U, s, V, ctrl );

    // Scale U with the inverted (nonzero) singular values, U := U / Sigma
    DiagonalSolve( RIGHT, NORMAL, s, U );

    // Form pinvA = (U Sigma V^H)^H = V (U Sigma)^H
    Gemm( NORMAL, ADJOINT, Field(1), V, U, A );
}

template<typename Field>
void HermitianPseudoinverse
( UpperOrLower uplo, AbstractDistMatrix<Field>& APre, Base<Field> tolerance )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    const Grid& g = A.Grid();

    // Get the EVD of A
    // TODO(poulson): Use a relative eigenvalue lower-bound
    DistMatrix<Real,VR,STAR> w(g);
    DistMatrix<Field> Z(g);
    HermitianEig( uplo, A, w, Z );

    if( tolerance == Real(0) )
    {
        // Set the tolerance equal to n ||A||_2 eps
        const Int n = Z.Height();
        const Real eps = limits::Epsilon<Real>();
        const Real twoNorm = MaxNorm( w );
        tolerance = n*twoNorm*eps;
    }
    // Invert above the tolerance
    auto omegaMap =
      [=]( const Real& omega )
      { return ( omega < tolerance ? Real(0) : 1/omega ); };
    EntrywiseMap( w, MakeFunction(omegaMap) );

    // Form the pseudoinverse
    HermitianFromEVD( uplo, A, w, Z );
}

#define PROTO(Field) \
  template void Pseudoinverse( Matrix<Field>& A, Base<Field> tolerance ); \
  template void Pseudoinverse \
  ( AbstractDistMatrix<Field>& A, Base<Field> tolerance ); \
  template void HermitianPseudoinverse \
  ( UpperOrLower uplo, Matrix<Field>& A, Base<Field> tolerance ); \
  template void HermitianPseudoinverse \
  ( UpperOrLower uplo, AbstractDistMatrix<Field>& A, Base<Field> tolerance );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
