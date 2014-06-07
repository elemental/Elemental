/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FUNCS_HPP
#define EL_FUNCS_HPP

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

} // namespace El

#include "El/lapack-like/factor.hpp"
#include "El/lapack-like/decomp.hpp"

#include "./funcs/HermitianFunction.hpp"

namespace El {

// Inverse
// =======
template<typename F>
void Inverse( Matrix<F>& A );
template<typename F>
void Inverse( DistMatrix<F>& A );
template<typename F>
void LocalInverse( DistMatrix<F,STAR,STAR>& A );
namespace inverse {
template<typename F>
void AfterLUPartialPiv
( Matrix<F>& A, const Matrix<Int>& pPerm );
template<typename F>
void AfterLUPartialPiv
( DistMatrix<F>& A, const DistMatrix<Int,VC,STAR>& pPerm );
} // namespace inverse

template<typename F>
void HPDInverse( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void HPDInverse( UpperOrLower uplo, DistMatrix<F>& A );
template<typename F>
void LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A );

template<typename F>
void HermitianInverse
( UpperOrLower uplo, Matrix<F>& A, LDLPivotType pivotType=BUNCH_KAUFMAN_A );
template<typename F>
void HermitianInverse
( UpperOrLower uplo, DistMatrix<F>& A, LDLPivotType pivotType=BUNCH_KAUFMAN_A );
template<typename F>
void LocalHermitianInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A );

template<typename F>
void SymmetricInverse
( UpperOrLower uplo, Matrix<F>& A, bool conjugate=false, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A );
template<typename F>
void SymmetricInverse
( UpperOrLower uplo, DistMatrix<F>& A, bool conjugate=false, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A );
template<typename F>
void LocalSymmetricInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A );

template<typename F>
void TriangularInverse( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A );
template<typename F>
void TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A  );
template<typename F>
void LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A );

// Pseudoinverse
// =============
template<typename F>
void Pseudoinverse( Matrix<F>& A, Base<F> tolerance=0 );
template<typename F>
void Pseudoinverse( DistMatrix<F>& A, Base<F> tolerance=0 );

template<typename F>
void HermitianPseudoinverse
( UpperOrLower uplo, Matrix<F>& A, Base<F> tolerance=0 );
template<typename F>
void HermitianPseudoinverse
( UpperOrLower uplo, DistMatrix<F>& A, Base<F> tolerance=0 );

// Sign
// ====
template<typename F>
void Sign
( Matrix<F>& A, const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );
template<typename F>
void Sign
( DistMatrix<F>& A, const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );

template<typename F>
void Sign
( Matrix<F>& A, Matrix<F>& N, 
  const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );
template<typename F>
void Sign
( DistMatrix<F>& A, DistMatrix<F>& N, 
  const SignCtrl<Base<F>> ctrl=SignCtrl<Base<F>>() );

template<typename F>
void HermitianSign
( UpperOrLower uplo, Matrix<F>& A, 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianSign
( UpperOrLower uplo, DistMatrix<F>& A, 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );

template<typename F>
void HermitianSign
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& N, 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianSign
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& N, 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );

// SquareRoot
// ==========
template<typename F>
void SquareRoot
( Matrix<F>& A,
  const SquareRootCtrl<Base<F>> ctrl=SquareRootCtrl<Base<F>>() );
template<typename F>
void SquareRoot
( DistMatrix<F>& A,
  const SquareRootCtrl<Base<F>> ctrl=SquareRootCtrl<Base<F>>() );

template<typename F>
void HPSDSquareRoot
( UpperOrLower uplo, Matrix<F>& A,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HPSDSquareRoot
( UpperOrLower uplo, DistMatrix<F>& A,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );

} // namespace El

#endif // ifndef EL_FUNCS_HPP
