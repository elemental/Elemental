/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_FUNCS_HPP
#define EL_FUNCS_HPP

#include <El/lapack_like/factor.hpp>
#include <El/lapack_like/spectral.hpp>

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
    Real tol=Real(0);
    Real power=Real(1);
    SignScaling scaling=SIGN_SCALE_FROB;
    bool progress=false;
};

template<typename Real>
struct SquareRootCtrl
{
    Int maxIts=100;
    Real tol=Real(0);
    Real power=Real(1);
    bool progress=false;
};

// Hermitian function
// ==================
template<typename Field>
void HermitianFunction
( UpperOrLower uplo, Matrix<Field>& A,
  function<Base<Field>(const Base<Field>&)> func );
template<typename Field>
void HermitianFunction
( UpperOrLower uplo, AbstractDistMatrix<Field>& A,
  function<Base<Field>(const Base<Field>&)> func );

template<typename Real>
void HermitianFunction
( UpperOrLower uplo, Matrix<Complex<Real>>& A,
  function<Complex<Real>(const Real&)> func );
template<typename Real>
void HermitianFunction
( UpperOrLower uplo, AbstractDistMatrix<Complex<Real>>& A,
  function<Complex<Real>(const Real&)> func );

// Inverse
// =======
template<typename Field>
void Inverse( Matrix<Field>& A );
template<typename Field>
void Inverse( AbstractDistMatrix<Field>& A );
template<typename Field>
void LocalInverse( DistMatrix<Field,STAR,STAR>& A );
namespace inverse {
template<typename Field>
void AfterLUPartialPiv
(       Matrix<Field>& A,
  const Permutation& P );
template<typename Field>
void AfterLUPartialPiv
(       AbstractDistMatrix<Field>& A,
  const DistPermutation& P );
} // namespace inverse

template<typename Field>
void HPDInverse( UpperOrLower uplo, Matrix<Field>& A );
template<typename Field>
void HPDInverse( UpperOrLower uplo, AbstractDistMatrix<Field>& A );
template<typename Field>
void LocalHPDInverse( UpperOrLower uplo, DistMatrix<Field,STAR,STAR>& A );

template<typename Field>
void HermitianInverse
( UpperOrLower uplo, Matrix<Field>& A,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );
template<typename Field>
void HermitianInverse
( UpperOrLower uplo, AbstractDistMatrix<Field>& A,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );
template<typename Field>
void LocalHermitianInverse
( UpperOrLower uplo, DistMatrix<Field,STAR,STAR>& A,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );

template<typename Field>
void SymmetricInverse
( UpperOrLower uplo, Matrix<Field>& A, bool conjugate=false,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );
template<typename Field>
void SymmetricInverse
( UpperOrLower uplo, AbstractDistMatrix<Field>& A, bool conjugate=false,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );
template<typename Field>
void LocalSymmetricInverse
( UpperOrLower uplo, DistMatrix<Field,STAR,STAR>& A, bool conjugate=false,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );

template<typename Field>
void TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<Field>& A );
template<typename Field>
void TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, AbstractDistMatrix<Field>& A  );
template<typename Field>
void LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<Field,STAR,STAR>& A );

// Pseudoinverse
// =============
template<typename Field>
void Pseudoinverse( Matrix<Field>& A, Base<Field> tolerance=0 );
template<typename Field>
void Pseudoinverse( AbstractDistMatrix<Field>& A, Base<Field> tolerance=0 );

template<typename Field>
void HermitianPseudoinverse
( UpperOrLower uplo, Matrix<Field>& A, Base<Field> tolerance=0 );
template<typename Field>
void HermitianPseudoinverse
( UpperOrLower uplo, AbstractDistMatrix<Field>& A, Base<Field> tolerance=0 );

// Sign
// ====
template<typename Field>
void Sign
( Matrix<Field>& A, const SignCtrl<Base<Field>> ctrl=SignCtrl<Base<Field>>() );
template<typename Field>
void Sign
( AbstractDistMatrix<Field>& A,
  const SignCtrl<Base<Field>> ctrl=SignCtrl<Base<Field>>() );

template<typename Field>
void Sign
( Matrix<Field>& A, Matrix<Field>& N,
  const SignCtrl<Base<Field>> ctrl=SignCtrl<Base<Field>>() );
template<typename Field>
void Sign
( AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& N,
  const SignCtrl<Base<Field>> ctrl=SignCtrl<Base<Field>>() );

template<typename Field>
void HermitianSign
( UpperOrLower uplo, Matrix<Field>& A,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );
template<typename Field>
void HermitianSign
( UpperOrLower uplo, AbstractDistMatrix<Field>& A,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );

template<typename Field>
void HermitianSign
( UpperOrLower uplo, Matrix<Field>& A, Matrix<Field>& N,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );
template<typename Field>
void HermitianSign
( UpperOrLower uplo, AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& N,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );

// SquareRoot
// ==========
template<typename Field>
void SquareRoot
( Matrix<Field>& A,
  const SquareRootCtrl<Base<Field>> ctrl=SquareRootCtrl<Base<Field>>() );
template<typename Field>
void SquareRoot
( AbstractDistMatrix<Field>& A,
  const SquareRootCtrl<Base<Field>> ctrl=SquareRootCtrl<Base<Field>>() );

template<typename Field>
void HPSDSquareRoot
( UpperOrLower uplo, Matrix<Field>& A,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );
template<typename Field>
void HPSDSquareRoot
( UpperOrLower uplo, AbstractDistMatrix<Field>& A,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );

} // namespace El

#endif // ifndef EL_FUNCS_HPP
