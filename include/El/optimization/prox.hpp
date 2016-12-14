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
void LowerClip( Matrix<Real>& X, const Real& lowerBound=0 );
template<typename Real>
void LowerClip( AbstractDistMatrix<Real>& X, const Real& lowerBound=0 );
template<typename Real>
void LowerClip( DistMultiVec<Real>& X, const Real& lowerBound=0 );

template<typename Real>
void UpperClip( Matrix<Real>& X, const Real& upperBound=0 );
template<typename Real>
void UpperClip( AbstractDistMatrix<Real>& X, const Real& upperBound=0 );
template<typename Real>
void UpperClip( DistMultiVec<Real>& X, const Real& upperBound=0 );

template<typename Real>
void Clip
( Matrix<Real>& X,
  const Real& lowerBound=0,
  const Real& upperBound=1 );
template<typename Real>
void Clip
( AbstractDistMatrix<Real>& X,
  const Real& lowerBound=0,
  const Real& upperBound=1 );
template<typename Real>
void Clip
( DistMultiVec<Real>& X,
  const Real& lowerBound=0,
  const Real& upperBound=1 );

// Frobenius-norm proximal map
// ---------------------------
// The Frobenius norm prox returns the solution to
//     arg min || A ||_F + rho/2 || A - A0 ||_F^2
//        A
// where A0 in the input matrix.
template<typename Field>
void FrobeniusProx( Matrix<Field>& A, const Base<Field>& rho );
template<typename Field>
void FrobeniusProx( AbstractDistMatrix<Field>& A, const Base<Field>& rho );

// Hinge-loss proximal map
// -----------------------
// TODO(poulson): Description
template<typename Real>
void HingeLossProx( Matrix<Real>& A, const Real& rho );
template<typename Real>
void HingeLossProx( AbstractDistMatrix<Real>& A, const Real& rho );

// Logistic proximal map
// ---------------------
// The logistic proximal map returns the solution to
//    arg min sum_{i,j}[ log(1+exp(-A_{i,j})) ] + rho/2 || A - A0 ||_F^2
//       A
// where A0 is the input matrix.
template<typename Real>
void LogisticProx
( Matrix<Real>& A, const Real& rho, Int numIts=5 );
template<typename Real>
void LogisticProx
( AbstractDistMatrix<Real>& A, const Real& rho, Int numIts=5 );

// Singular-value soft thresholding
// --------------------------------
template<typename Field>
Int SVT
( Matrix<Field>& A,
  const Base<Field>& rho,
  bool relative=false );
template<typename Field>
Int SVT
( AbstractDistMatrix<Field>& A,
  const Base<Field>& rho,
  bool relative=false );
template<typename Field>
Int SVT
( Matrix<Field>& A,
  const Base<Field>& rho,
  Int relaxedRank,
  bool relative=false );
template<typename Field>
Int SVT
( AbstractDistMatrix<Field>& A,
  const Base<Field>& rho,
  Int relaxedRank,
  bool relative=false );
template<typename Field,Dist U>
Int SVT
( DistMatrix<Field,U,STAR>& A,
  const Base<Field>& rho,
  bool relative=false );

namespace svt {

// TODO(poulson): Add SVT control structure

template<typename Field>
Int Cross
( Matrix<Field>& A,
  const Base<Field>& rho,
  bool relative=false );
template<typename Field>
Int Cross
( AbstractDistMatrix<Field>& A,
  const Base<Field>& rho,
  bool relative=false );
template<typename Field>
Int Cross
( DistMatrix<Field,VC,STAR>& A,
  const Base<Field>& rho,
  bool relative=false );

template<typename Field>
Int Normal
( Matrix<Field>& A,
  const Base<Field>& rho,
  bool relative=false );
template<typename Field>
Int Normal
( AbstractDistMatrix<Field>& A,
  const Base<Field>& rho,
  bool relative=false );

template<typename Field>
Int PivotedQR
( Matrix<Field>& A,
  const Base<Field>& rho,
  Int numSteps,
  bool relative=false );
template<typename Field>
Int PivotedQR
( AbstractDistMatrix<Field>& A,
  const Base<Field>& rho,
  Int numSteps,
  bool relative=false );

template<typename Field>
Int TSQR
( AbstractDistMatrix<Field>& A,
  const Base<Field>& rho,
  bool relative=false );

} // namespace svt

// Soft-thresholding
// -----------------
// Returns the solution to
//     arg min || vec(A) ||_1 + rho/2 || A - A0 ||_F^2
//        A
// where A0 is the input matrix.
template<typename Field>
Field SoftThreshold( const Field& alpha, const Base<Field>& rho );

template<typename Field>
void SoftThreshold
( Matrix<Field>& A, const Base<Field>& rho, bool relative=false );
template<typename Field>
void SoftThreshold
( AbstractDistMatrix<Field>& A, const Base<Field>& rho, bool relative=false );

} // namespace El

#endif // ifndef EL_OPTIMIZATION_PROX_HPP
