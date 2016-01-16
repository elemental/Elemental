/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Replace A with its pseudoinverse

// TODO: Avoid unnecessary work for backtransformation with zero singular
//       values

template<typename F>
void Pseudoinverse( Matrix<F>& A, Base<F> tolerance )
{
    DEBUG_ONLY(CSE cse("Pseudoinverse"))
    typedef Base<F> Real;

    // Get the SVD of A
    Matrix<Real> s;
    Matrix<F> U, V;
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    SVD( A, U, s, V );

    if( tolerance == Real(0) )
    {
        // Set the tolerance equal to k ||A||_2 eps
        const Int k = Max( U.Height(), V.Height() );
        const Real eps = limits::Epsilon<Real>();
        const Real twoNorm = MaxNorm( s );
        tolerance = k*twoNorm*eps;
    }
    // Invert above the tolerance
    auto sigmaMap = 
      [=]( Real sigma ) { return ( sigma < tolerance ? Real(0) : 1/sigma ); };
    EntrywiseMap( s, function<Real(Real)>(sigmaMap) );

    // Scale U with the singular values, U := U Sigma
    DiagonalScale( RIGHT, NORMAL, s, U );

    // Form pinvA = (U Sigma V^H)^H = V (U Sigma)^H
    Gemm( NORMAL, ADJOINT, F(1), V, U, A );
}

template<typename F>
void HermitianPseudoinverse
( UpperOrLower uplo, Matrix<F>& A, Base<F> tolerance )
{
    DEBUG_ONLY(CSE cse("HermitianPseudoinverse"))
    typedef Base<F> Real;

    // Get the EVD of A
    Matrix<Real> w;
    Matrix<F> Z;
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
      [=]( Real omega ) { return ( omega < tolerance ? Real(0) : 1/omega ); };
    EntrywiseMap( w, function<Real(Real)>(omegaMap) );

    // Form the pseudoinverse
    HermitianFromEVD( uplo, A, w, Z );
}

template<typename F>
void Pseudoinverse( ElementalMatrix<F>& APre, Base<F> tolerance )
{
    DEBUG_ONLY(CSE cse("Pseudoinverse"))
    typedef Base<F> Real;

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    const Grid& g = A.Grid();

    // Get the SVD of A
    DistMatrix<Real,VR,STAR> s(g);
    DistMatrix<F> U(g), V(g);
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    SVD( A, U, s, V, ctrl );

    if( tolerance == Real(0) )
    {
        // Set the tolerance equal to k ||A||_2 eps
        const Int k = Max( U.Height(), V.Height() );
        const Real eps = limits::Epsilon<Real>();
        const Real twoNorm = MaxNorm( s );
        tolerance = k*twoNorm*eps;
    }
    // Invert above the tolerance
    auto sigmaMap = 
      [=]( Real sigma ) { return ( sigma < tolerance ? Real(0) : 1/sigma ); };
    EntrywiseMap( s, function<Real(Real)>(sigmaMap) );

    // Scale U with the singular values, U := U Sigma
    DiagonalScale( RIGHT, NORMAL, s, U );

    // Form pinvA = (U Sigma V^H)^H = V (U Sigma)^H
    Gemm( NORMAL, ADJOINT, F(1), V, U, A );
}

template<typename F>
void HermitianPseudoinverse
( UpperOrLower uplo, ElementalMatrix<F>& APre, Base<F> tolerance )
{
    DEBUG_ONLY(CSE cse("HermitianPseudoinverse"))
    typedef Base<F> Real;

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    const Grid& g = A.Grid();

    // Get the EVD of A
    DistMatrix<Real,VR,STAR> w(g);
    DistMatrix<F> Z(g);
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
      [=]( Real omega ) { return ( omega < tolerance ? Real(0) : 1/omega ); };
    EntrywiseMap( w, function<Real(Real)>(omegaMap) );

    // Form the pseudoinverse
    HermitianFromEVD( uplo, A, w, Z );
}

#define PROTO(F) \
  template void Pseudoinverse( Matrix<F>& A, Base<F> tolerance ); \
  template void Pseudoinverse( ElementalMatrix<F>& A, Base<F> tolerance ); \
  template void HermitianPseudoinverse \
  ( UpperOrLower uplo, Matrix<F>& A, Base<F> tolerance ); \
  template void HermitianPseudoinverse \
  ( UpperOrLower uplo, ElementalMatrix<F>& A, Base<F> tolerance );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
