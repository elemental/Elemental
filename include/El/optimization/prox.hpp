/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_PROX_HPP
#define EL_OPTIMIZATION_PROX_HPP

namespace El {

// Proximal maps
// =============

// Clipping
// --------
template<typename Real>
void LowerClip( Matrix<Real>& X, Real lowerBound=0 );
template<typename Real>
void LowerClip( AbstractDistMatrix<Real>& X, Real lowerBound=0 );
template<typename Real>
void LowerClip( DistMultiVec<Real>& X, Real lowerBound=0 );

template<typename Real>
void UpperClip( Matrix<Real>& X, Real upperBound=0 );
template<typename Real>
void UpperClip( AbstractDistMatrix<Real>& X, Real upperBound=0 );
template<typename Real>
void UpperClip( DistMultiVec<Real>& X, Real upperBound=0 );

template<typename Real>
void Clip
( Matrix<Real>& X, Real lowerBound=0, Real upperBound=1 );
template<typename Real>
void Clip
( AbstractDistMatrix<Real>& X, Real lowerBound=0, Real upperBound=1 );
template<typename Real>
void Clip
( DistMultiVec<Real>& X, Real lowerBound=0, Real upperBound=1 );

// Frobenius-norm proximal map
// ---------------------------
// The Frobenius norm prox returns the solution to
//     arg min || A ||_F + rho/2 || A - A0 ||_F^2
//        A
// where A0 in the input matrix.
template<typename F>
void FrobeniusProx( Matrix<F>& A, Base<F> rho );
template<typename F>
void FrobeniusProx( AbstractDistMatrix<F>& A, Base<F> rho );

// Hinge-loss proximal map
// -----------------------
// TODO: Description
template<typename Real>
void HingeLossProx( Matrix<Real>& A, Real rho );
template<typename Real>
void HingeLossProx( AbstractDistMatrix<Real>& A, Real rho );

// Logistic proximal map
// ---------------------
// The logistic proximal map returns the solution to
//    arg min sum_{i,j}[ log(1+exp(-A_{i,j})) ] + rho/2 || A - A0 ||_F^2
//       A
// where A0 is the input matrix.
template<typename Real>
void LogisticProx( Matrix<Real>& A, Real rho, Int numIts=5 );
template<typename Real>
void LogisticProx( AbstractDistMatrix<Real>& A, Real rho, Int numIts=5 );

// Singular-value soft thresholding
// --------------------------------
template<typename F>
Int SVT( Matrix<F>& A, Base<F> rho, bool relative=false );
template<typename F>
Int SVT( ElementalMatrix<F>& A, Base<F> rho, bool relative=false );
template<typename F>
Int SVT( Matrix<F>& A, Base<F> rho, Int relaxedRank, bool relative=false );
template<typename F>
Int SVT
( ElementalMatrix<F>& A, Base<F> rho, Int relaxedRank, bool relative=false );
template<typename F,Dist U>
Int SVT( DistMatrix<F,U,STAR>& A, Base<F> rho, bool relative=false );

namespace svt {

// TODO: Add SVT control structure

template<typename F>
Int Cross( Matrix<F>& A, Base<F> rho, bool relative=false );
template<typename F>
Int Cross( ElementalMatrix<F>& A, Base<F> rho, bool relative=false );
template<typename F>
Int Cross( DistMatrix<F,VC,STAR>& A, Base<F> rho, bool relative=false );

template<typename F>
Int Normal( Matrix<F>& A, Base<F> rho, bool relative=false );
template<typename F>
Int Normal( ElementalMatrix<F>& A, Base<F> rho, bool relative=false );

template<typename F>
Int PivotedQR
( Matrix<F>& A, Base<F> rho, Int numSteps, bool relative=false );
template<typename F>
Int PivotedQR
( ElementalMatrix<F>& A, Base<F> rho, Int numSteps, bool relative=false );

template<typename F>
Int TSQR( ElementalMatrix<F>& A, Base<F> rho, bool relative=false );

} // namespace svt

// Soft-thresholding
// -----------------
// Returns the solution to
//     arg min || vec(A) ||_1 + rho/2 || A - A0 ||_F^2
//        A 
// where A0 is the input matrix.
template<typename F>
F SoftThreshold( F alpha, Base<F> rho );

template<typename F>
void SoftThreshold
( Matrix<F>& A, Base<F> rho, bool relative=false );
template<typename F>
void SoftThreshold
( AbstractDistMatrix<F>& A, Base<F> rho, bool relative=false );

} // namespace El

#endif // ifndef EL_OPTIMIZATION_PROX_HPP
