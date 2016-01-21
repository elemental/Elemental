/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
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
struct SignCtrl 
{
    Int maxIts=100;
    Real tol=0;
    Real power=1;
    SignScaling scaling=SIGN_SCALE_FROB;
    bool progress=false;
};

template<typename Real>
struct SquareRootCtrl 
{
    Int maxIts=100;
    Real tol=0;
    Real power=1;
    bool progress=false;
};

// Hermitian function
// ==================
template<typename F>
void HermitianFunction
( UpperOrLower uplo, Matrix<F>& A, function<Base<F>(Base<F>)> func );
template<typename F>
void HermitianFunction
( UpperOrLower uplo, ElementalMatrix<F>& A, 
  function<Base<F>(Base<F>)> func );

template<typename Real>
void HermitianFunction
( UpperOrLower uplo, Matrix<Complex<Real>>& A,
  function<Complex<Real>(Real)> func );
template<typename Real>
void HermitianFunction
( UpperOrLower uplo, ElementalMatrix<Complex<Real>>& A,
  function<Complex<Real>(Real)> func );

// Inverse
// =======
template<typename F>
void Inverse( Matrix<F>& A );
template<typename F>
void Inverse( ElementalMatrix<F>& A );
template<typename F>
void LocalInverse( DistMatrix<F,STAR,STAR>& A );
namespace inverse {
template<typename F>
void AfterLUPartialPiv
(       Matrix<F>& A,
  const Permutation& P );
template<typename F>
void AfterLUPartialPiv
(       ElementalMatrix<F>& A,
  const DistPermutation& P );
} // namespace inverse

template<typename F>
void HPDInverse( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void HPDInverse( UpperOrLower uplo, ElementalMatrix<F>& A );
template<typename F>
void LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A );

template<typename F>
void HermitianInverse
( UpperOrLower uplo, Matrix<F>& A, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void HermitianInverse
( UpperOrLower uplo, ElementalMatrix<F>& A, 
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
( UpperOrLower uplo, ElementalMatrix<F>& A, bool conjugate=false, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void LocalSymmetricInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, bool conjugate=false,
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );

template<typename F>
void TriangularInverse( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A );
template<typename F>
void TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, ElementalMatrix<F>& A  );
template<typename F>
void LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A );

// Pseudoinverse
// =============
template<typename F>
void Pseudoinverse( Matrix<F>& A, Base<F> tolerance=0 );
template<typename F>
void Pseudoinverse( ElementalMatrix<F>& A, Base<F> tolerance=0 );

template<typename F>
void HermitianPseudoinverse
( UpperOrLower uplo, Matrix<F>& A, Base<F> tolerance=0 );
template<typename F>
void HermitianPseudoinverse
( UpperOrLower uplo, ElementalMatrix<F>& A, Base<F> tolerance=0 );

// Sign
// ====
template<typename F>
void Sign
( Matrix<F>& A, const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );
template<typename F>
void Sign
( ElementalMatrix<F>& A, const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );

template<typename F>
void Sign
( Matrix<F>& A, Matrix<F>& N, 
  const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );
template<typename F>
void Sign
( ElementalMatrix<F>& A, ElementalMatrix<F>& N, 
  const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );

template<typename F>
void HermitianSign
( UpperOrLower uplo, Matrix<F>& A, 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HermitianSign
( UpperOrLower uplo, ElementalMatrix<F>& A, 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

template<typename F>
void HermitianSign
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& N, 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HermitianSign
( UpperOrLower uplo, ElementalMatrix<F>& A, ElementalMatrix<F>& N, 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

// SquareRoot
// ==========
template<typename F>
void SquareRoot
( Matrix<F>& A,
  const SquareRootCtrl<Base<F>> ctrl=SquareRootCtrl<Base<F>>() );
template<typename F>
void SquareRoot
( ElementalMatrix<F>& A,
  const SquareRootCtrl<Base<F>> ctrl=SquareRootCtrl<Base<F>>() );

template<typename F>
void HPSDSquareRoot
( UpperOrLower uplo, Matrix<F>& A,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HPSDSquareRoot
( UpperOrLower uplo, ElementalMatrix<F>& A,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

} // namespace El

#endif // ifndef EL_FUNCS_HPP
