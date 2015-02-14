/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FUNCS_HPP
#define EL_FUNCS_HPP

#include "El/lapack_like/factor.hpp"
#include "El/lapack_like/spectral.hpp"

namespace El {

namespace SignScalingNS {
enum SignScaling {
    SIGN_SCALE_NONE,
    SIGN_SCALE_DET,
    SIGN_SCALE_FROB
};
}
using namespace SignScalingNS;

template<typename Real>
struct SignCtrl {
    Int maxIts;
    Real tol;
    Real power;
    SignScaling scaling;
    bool progress;

    SignCtrl()
    : maxIts(100), tol(0), power(1), scaling(SIGN_SCALE_FROB), progress(false)
    { }
};

template<typename Real>
struct SquareRootCtrl {
    Int maxIts;
    Real tol;
    Real power;
    bool progress;

    SquareRootCtrl()
    : maxIts(100), tol(0), power(1), progress(false)
    { }
};

// Hermitian function
// ==================
template<typename F>
void HermitianFunction
( UpperOrLower uplo, Matrix<F>& A, function<Base<F>(Base<F>)> func );
template<typename F>
void HermitianFunction
( UpperOrLower uplo, AbstractDistMatrix<F>& A, 
  function<Base<F>(Base<F>)> func );

template<typename Real>
void HermitianFunction
( UpperOrLower uplo, Matrix<Complex<Real>>& A,
  function<Complex<Real>(Real)> func );
template<typename Real>
void HermitianFunction
( UpperOrLower uplo, AbstractDistMatrix<Complex<Real>>& A,
  function<Complex<Real>(Real)> func );

// Inverse
// =======
template<typename F>
void Inverse( Matrix<F>& A );
template<typename F>
void Inverse( AbstractDistMatrix<F>& A );
template<typename F>
void LocalInverse( DistMatrix<F,STAR,STAR>& A );
namespace inverse {
template<typename F>
void AfterLUPartialPiv( Matrix<F>& A, const Matrix<Int>& p );
template<typename F>
void AfterLUPartialPiv
( AbstractDistMatrix<F>& A, const AbstractDistMatrix<Int>& p );
} // namespace inverse

template<typename F>
void HPDInverse( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void HPDInverse( UpperOrLower uplo, AbstractDistMatrix<F>& A );
template<typename F>
void LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A );

template<typename F>
void HermitianInverse
( UpperOrLower uplo, Matrix<F>& A, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void HermitianInverse
( UpperOrLower uplo, AbstractDistMatrix<F>& A, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void LocalHermitianInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );

template<typename F>
void SymmetricInverse
( UpperOrLower uplo, Matrix<F>& A, bool conjugate=false, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void SymmetricInverse
( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool conjugate=false, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void LocalSymmetricInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, bool conjugate=false,
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );

template<typename F>
void TriangularInverse( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A );
template<typename F>
void TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, AbstractDistMatrix<F>& A  );
template<typename F>
void LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A );

// Pseudoinverse
// =============
template<typename F>
void Pseudoinverse( Matrix<F>& A, Base<F> tolerance=0 );
template<typename F>
void Pseudoinverse( AbstractDistMatrix<F>& A, Base<F> tolerance=0 );

template<typename F>
void HermitianPseudoinverse
( UpperOrLower uplo, Matrix<F>& A, Base<F> tolerance=0 );
template<typename F>
void HermitianPseudoinverse
( UpperOrLower uplo, AbstractDistMatrix<F>& A, Base<F> tolerance=0 );

// Sign
// ====
template<typename F>
void Sign
( Matrix<F>& A, const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );
template<typename F>
void Sign
( AbstractDistMatrix<F>& A, const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );

template<typename F>
void Sign
( Matrix<F>& A, Matrix<F>& N, 
  const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );
template<typename F>
void Sign
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& N, 
  const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );

template<typename F>
void HermitianSign
( UpperOrLower uplo, Matrix<F>& A, 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HermitianSign
( UpperOrLower uplo, AbstractDistMatrix<F>& A, 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

template<typename F>
void HermitianSign
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& N, 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HermitianSign
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& N, 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

// SquareRoot
// ==========
template<typename F>
void SquareRoot
( Matrix<F>& A,
  const SquareRootCtrl<Base<F>> ctrl=SquareRootCtrl<Base<F>>() );
template<typename F>
void SquareRoot
( AbstractDistMatrix<F>& A,
  const SquareRootCtrl<Base<F>> ctrl=SquareRootCtrl<Base<F>>() );

template<typename F>
void HPSDSquareRoot
( UpperOrLower uplo, Matrix<F>& A,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HPSDSquareRoot
( UpperOrLower uplo, AbstractDistMatrix<F>& A,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

} // namespace El

#endif // ifndef EL_FUNCS_HPP
