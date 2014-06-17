/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS1_HPP
#define EL_BLAS1_HPP

namespace El {

// Adjoint
// =======
template<typename T>
void Adjoint( const Matrix<T>& A, Matrix<T>& B );
template<typename T>
void Adjoint( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );
template<typename T>
void Adjoint
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

// Axpy
// ====
template<typename T,typename S>
void Axpy( S alpha, const Matrix<T>& X, Matrix<T>& Y );
template<typename T,typename S,Dist U,Dist V,Dist W,Dist Z>
void Axpy( S alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,W,Z>& Y );
template<typename T,typename S>
void Axpy( S alpha, const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y );

// AxpyTriangle
// ============
template<typename T,typename S>
void AxpyTriangle
( UpperOrLower uplo, S alpha, const Matrix<T>& X, Matrix<T>& Y );
template<typename T,typename S,Dist U,Dist V>
void AxpyTriangle
( UpperOrLower uplo, S alpha,
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y );
template<typename T,typename S>
void AxpyTriangle
( UpperOrLower uplo, S alpha,
  const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y );

// Conjugate
// =========
template<typename Real>
void Conjugate( Matrix<Real>& A );
template<typename Real>
void Conjugate( Matrix<Complex<Real>>& A );

template<typename T>
void Conjugate( const Matrix<T>& A, Matrix<T>& B );

template<typename T>
void Conjugate( AbstractDistMatrix<T>& A );

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );
template<typename T>
void Conjugate( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );

// Copy
// ====
template<typename T>
void Copy( const Matrix<T>& A, Matrix<T>& B );

template<typename Real>
void Copy( const Matrix<Real>& A, Matrix<Complex<Real>>& B );

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

template<typename Real,Dist U,Dist V,Dist W,Dist Z>
void Copy( const DistMatrix<Real,U,V>& A, DistMatrix<Complex<Real>,W,Z>& B );

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Copy( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,W,Z>& B );

template<typename Real,Dist U,Dist V,Dist W,Dist Z>
void Copy
( const BlockDistMatrix<Real,U,V>& A, BlockDistMatrix<Complex<Real>,W,Z>& B );

template<typename T>
void Copy( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );

// DiagonalScale
// =============
template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<TDiag>& d, Matrix<T>& X );

template<typename TDiag,typename T,Dist U,Dist V,Dist W,Dist Z>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMatrix<TDiag,U,V>& d, DistMatrix<T,W,Z>& X );

template<typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X );

template<typename Real>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<Real>& d, AbstractDistMatrix<Complex<Real>>& X );

// DiagonalScaleTrapezoid
// ======================
template<typename TDiag,typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const Matrix<TDiag>& d, Matrix<T>& A, Int offset=0 );

template<typename TDiag,typename T,Dist U,Dist V,Dist W,Dist Z>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const DistMatrix<TDiag,U,V>& d, DistMatrix<T,W,Z>& A, Int offset=0 );

template<typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& A, Int offset=0 );

template<typename Real>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<Real>& d,
        AbstractDistMatrix<Complex<Real>>& A, Int offset=0 );

// DiagonalSolve
// =============
template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<FDiag>& d, Matrix<F>& X, bool checkIfSingular=true );

template<typename FDiag,typename F,Dist U,Dist V,Dist W,Dist Z>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMatrix<FDiag,U,V>& d, DistMatrix<F,W,Z>& X,
  bool checkIfSingular=true );

template<typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& d, AbstractDistMatrix<F>& X,
  bool checkIfSingular=true );

template<typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& d, AbstractDistMatrix<Complex<F>>& X,
  bool checkIfSingular=true );

// Dot
// ===
template<typename F>
F Dot( const Matrix<F>& A, const Matrix<F>& B );
template<typename F>
F Dot( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B );

// Dotu
// ====
template<typename F>
F Dotu( const Matrix<F>& A, const Matrix<F>& B );
template<typename F>
F Dotu( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B );

// Fill
// ====
template<typename T>
void Fill( Matrix<T>& A, T alpha );
template<typename T>
void Fill( AbstractDistMatrix<T>& A, T alpha );
template<typename T>
void Fill( AbstractBlockDistMatrix<T>& A, T alpha );

// Hadamard
// ========
template<typename T>
void Hadamard( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C );
template<typename T,Dist U,Dist V>
void Hadamard
( const DistMatrix<T,U,V>& A,
  const DistMatrix<T,U,V>& B,
        DistMatrix<T,U,V>& C );
template<typename T>
void Hadamard
( const AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& C );

// HilbertSchmidt
// ==============
template<typename F>
F HilbertSchmidt( const Matrix<F>& A, const Matrix<F>& B );
template<typename F>
F HilbertSchmidt
( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& C );

// MakeHermitian
// =============
template<typename T>
void MakeHermitian( UpperOrLower uplo, Matrix<T>& A );
template<typename T>
void MakeHermitian( UpperOrLower uplo, AbstractDistMatrix<T>& A );

// MakeReal
// ========
template<typename Real>
void MakeReal( Matrix<Real>& A );
template<typename Real>
void MakeReal( Matrix<Complex<Real>>& A );
template<typename T>
void MakeReal( AbstractDistMatrix<T>& A );

// MakeSymmetric
// =============
template<typename T>
void MakeSymmetric( UpperOrLower uplo, Matrix<T>& A, bool conjugate=false );
template<typename T,Dist U,Dist V>
void MakeSymmetric
( UpperOrLower uplo, DistMatrix<T,U,V>& A, bool conjugate=false );
template<typename F>
void MakeSymmetric
( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool conjugate=false );

// MakeTrapezoidal
// ===============
template<typename T>
void MakeTrapezoidal( UpperOrLower uplo, Matrix<T>& A, Int offset=0 );
template<typename T>
void MakeTrapezoidal
( UpperOrLower uplo, AbstractDistMatrix<T>& A, Int offset=0 );
template<typename T>
void MakeTrapezoidal
( UpperOrLower uplo, AbstractBlockDistMatrix<T>& A, Int offset=0 );

// MakeTriangular
// ==============
template<typename T>
void MakeTriangular( UpperOrLower uplo, Matrix<T>& A );
template<typename T>
void MakeTriangular( UpperOrLower uplo, AbstractDistMatrix<T>& A );
template<typename T>
void MakeTriangular( UpperOrLower uplo, AbstractBlockDistMatrix<T>& A );

// Max
// ===
template<typename Real>
ValueInt<Real> VectorMax( const Matrix<Real>& x );
template<typename Real>
ValueInt<Real> VectorMax( const AbstractDistMatrix<Real>& x );

template<typename Real>
ValueIntPair<Real> Max( const Matrix<Real>& A );
template<typename Real>
ValueIntPair<Real> Max( const AbstractDistMatrix<Real>& A );

template<typename Real>
ValueIntPair<Real> SymmetricMax( UpperOrLower uplo, const Matrix<Real>& A );
template<typename Real>
ValueIntPair<Real>
SymmetricMax( UpperOrLower uplo, const AbstractDistMatrix<Real>& A );

template<typename Real>
ValueInt<Real> DiagonalMax( const Matrix<Real>& A );
template<typename Real,Dist U,Dist V>
ValueInt<Real> DiagonalMax( const DistMatrix<Real,U,V>& A );

// MaxAbs
// ======
template<typename F>
ValueInt<Base<F>> VectorMaxAbs( const Matrix<F>& x );
template<typename F>
ValueInt<Base<F>> VectorMaxAbs( const AbstractDistMatrix<F>& x );

template<typename F>
ValueIntPair<Base<F>> MaxAbs( const Matrix<F>& A );
template<typename F>
ValueIntPair<Base<F>> MaxAbs( const AbstractDistMatrix<F>& A );

template<typename F>
ValueIntPair<Base<F>> SymmetricMaxAbs( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
ValueIntPair<Base<F>> SymmetricMaxAbs
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
ValueInt<Base<F>> DiagonalMaxAbs( const Matrix<F>& A );
template<typename F,Dist U,Dist V>
ValueInt<Base<F>> DiagonalMaxAbs( const DistMatrix<F,U,V>& A );

// Min
// ===
template<typename Real>
ValueInt<Real> VectorMin( const Matrix<Real>& x );
template<typename Real>
ValueInt<Real> VectorMin( const AbstractDistMatrix<Real>& x );

template<typename Real>
ValueIntPair<Real> Min( const Matrix<Real>& A );
template<typename Real>
ValueIntPair<Real> Min( const AbstractDistMatrix<Real>& A );

template<typename Real>
ValueIntPair<Real> SymmetricMin( UpperOrLower uplo, const Matrix<Real>& A );
template<typename Real>
ValueIntPair<Real>
SymmetricMin( UpperOrLower uplo, const AbstractDistMatrix<Real>& A );

template<typename Real>
ValueInt<Real> DiagonalMin( const Matrix<Real>& A );
template<typename Real,Dist U,Dist V>
ValueInt<Real> DiagonalMin( const DistMatrix<Real,U,V>& A );

// MinAbs
// ======
template<typename F>
ValueInt<Base<F>> VectorMinAbs( const Matrix<F>& x );
template<typename F>
ValueInt<Base<F>> VectorMinAbs( const AbstractDistMatrix<F>& x );

template<typename F>
ValueIntPair<Base<F>> MinAbs( const Matrix<F>& A );
template<typename F>
ValueIntPair<Base<F>> MinAbs( const AbstractDistMatrix<F>& A );

template<typename F>
ValueIntPair<Base<F>> SymmetricMinAbs( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
ValueIntPair<Base<F>>
SymmetricMinAbs( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
ValueInt<Base<F>> DiagonalMinAbs( const Matrix<F>& A );
template<typename F,Dist U,Dist V>
ValueInt<Base<F>> DiagonalMinAbs( const DistMatrix<F,U,V>& A );

// Nrm2
// ====
template<typename F>
Base<F> Nrm2( const Matrix<F>& x );
template<typename F>
Base<F> Nrm2( const AbstractDistMatrix<F>& x );

// QuasiDiagonalScale
// ==================
template<typename F,typename FMain>
void QuasiDiagonalScale
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<FMain>& d, const Matrix<F>& dSub,
  Matrix<F>& X, bool conjugated=false );
template<typename F,typename FMain,Dist U1,Dist V1,Dist U2,Dist V2>
void QuasiDiagonalScale
( LeftOrRight side, UpperOrLower uplo,
  const DistMatrix<FMain,U1,V1>& d, const DistMatrix<F,U1,V1>& dSub,
  DistMatrix<F,U2,V2>& X, bool conjugated=false );

template<typename F,typename FMain,Dist U,Dist V>
void LeftQuasiDiagonalScale
( UpperOrLower uplo,
  const DistMatrix<FMain,U,STAR>& d,
  const DistMatrix<FMain,U,STAR>& dPrev,
  const DistMatrix<FMain,U,STAR>& dNext,
  const DistMatrix<F,    U,STAR>& dSub,
  const DistMatrix<F,    U,STAR>& dSubPrev,
  const DistMatrix<F,    U,STAR>& dSubNext,
        DistMatrix<F,U,V>& X,
  const DistMatrix<F,U,V>& XPrev,
  const DistMatrix<F,U,V>& XNext,
  bool conjugated=false );

template<typename F,typename FMain,Dist U,Dist V>
void RightQuasiDiagonalScale
( UpperOrLower uplo,
  const DistMatrix<FMain,V,STAR>& d,
  const DistMatrix<FMain,V,STAR>& dPrev,
  const DistMatrix<FMain,V,STAR>& dNext,
  const DistMatrix<F,    V,STAR>& dSub,
  const DistMatrix<F,    V,STAR>& dSubPrev,
  const DistMatrix<F,    V,STAR>& dSubNext,
        DistMatrix<F,U,V>& X,
  const DistMatrix<F,U,V>& XPrev,
  const DistMatrix<F,U,V>& XNext,
  bool conjugated=false );

// QuasiDiagonalSolve
// ==================
template<typename F,typename FMain>
void
QuasiDiagonalSolve
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<FMain>& d, const Matrix<F>& dSub,
  Matrix<F>& X, bool conjugated=false );
template<typename F,typename FMain,Dist U1,Dist V1,Dist U2,Dist V2>
void
QuasiDiagonalSolve
( LeftOrRight side, UpperOrLower uplo,
  const DistMatrix<FMain,U1,V1>& d, const DistMatrix<F,U1,V1>& dSub,
  DistMatrix<F,U2,V2>& X, bool conjugated=false );

template<typename F,typename FMain,Dist U,Dist V>
void
LeftQuasiDiagonalSolve
( UpperOrLower uplo,
  const DistMatrix<FMain,U,STAR>& d,
  const DistMatrix<FMain,U,STAR>& dPrev,
  const DistMatrix<FMain,U,STAR>& dNext,
  const DistMatrix<F,    U,STAR>& dSub,
  const DistMatrix<F,    U,STAR>& dSubPrev,
  const DistMatrix<F,    U,STAR>& dSubNext,
        DistMatrix<F,U,V>& X,
  const DistMatrix<F,U,V>& XPrev,
  const DistMatrix<F,U,V>& XNext,
  bool conjugated=false );

template<typename F,typename FMain,Dist U,Dist V>
void
RightQuasiDiagonalSolve
( UpperOrLower uplo,
  const DistMatrix<FMain,V,STAR>& d,
  const DistMatrix<FMain,V,STAR>& dPrev,
  const DistMatrix<FMain,V,STAR>& dNext,
  const DistMatrix<F,    V,STAR>& dSub,
  const DistMatrix<F,    V,STAR>& dSubPrev,
  const DistMatrix<F,    V,STAR>& dSubNext,
        DistMatrix<F,U,V>& X,
  const DistMatrix<F,U,V>& XPrev,
  const DistMatrix<F,U,V>& XNext,
  bool conjugated=false );

// Scale
// =====
template<typename T,typename S>
void Scale( S alpha, Matrix<T>& A );
template<typename T,typename S>
void Scale( S alpha, AbstractDistMatrix<T>& A );
template<typename T,typename S>
void Scale( S alpha, AbstractBlockDistMatrix<T>& A );

template<typename Real,typename S>
void Scale( S alpha, Matrix<Real>& AReal, Matrix<Real>& AImag );
template<typename Real,typename S>
void Scale
( S alpha, AbstractDistMatrix<Real>& AReal, AbstractDistMatrix<Real>& AImag );
template<typename Real,typename S>
void Scale
( S alpha, AbstractBlockDistMatrix<Real>& AReal,
           AbstractBlockDistMatrix<Real>& AImag );

// ScaleTrapezoid
// ==============
template<typename T,typename S>
void ScaleTrapezoid( S alpha, UpperOrLower uplo, Matrix<T>& A, Int offset=0 );
template<typename T,typename S>
void ScaleTrapezoid
( S alpha, UpperOrLower uplo, AbstractDistMatrix<T>& A, Int offset=0 );

// SetDiagonal
// ===========
template<typename T,typename S>
void SetDiagonal( Matrix<T>& A, S alpha, Int offset=0 );
template<typename T,typename S>
void SetDiagonal( AbstractDistMatrix<T>& A, S alpha, Int offset=0 );
template<typename T,typename S>
void SetDiagonal( AbstractBlockDistMatrix<T>& A, S alpha, Int offset=0 );

// Swap
// ====
template<typename F>
void Swap( Orientation orientation, Matrix<F>& X, Matrix<F>& Y );
template<typename F,Dist U1,Dist V1,Dist U2,Dist V2>
void Swap
( Orientation orientation, DistMatrix<F,U1,V1>& X, DistMatrix<F,U2,V2>& Y );

template<typename F>
void RowSwap( Matrix<F>& A, Int to, Int from );
template<typename F,Dist U,Dist V>
void RowSwap( DistMatrix<F,U,V>& A, Int to, Int from );

template<typename F>
void ColSwap( Matrix<F>& A, Int to, Int from );
template<typename F,Dist U,Dist V>
void ColSwap( DistMatrix<F,U,V>& A, Int to, Int from );

template<typename F>
void SymmetricSwap
( UpperOrLower uplo, Matrix<F>& A, Int to, Int from, bool conjugate=false );
template<typename F,Dist U,Dist V>
void SymmetricSwap
( UpperOrLower uplo, DistMatrix<F,U,V>& A, Int to, Int from,
  bool conjugate=false );

template<typename F>
void HermitianSwap( UpperOrLower uplo, Matrix<F>& A, Int to, Int from );
template<typename F,Dist U,Dist V>
void HermitianSwap( UpperOrLower uplo, DistMatrix<F,U,V>& A, Int to, Int from );

// Symmetric2x2Inv
// ===============
template<typename F>
void Symmetric2x2Inv( UpperOrLower uplo, Matrix<F>& D, bool conjugate=false );

// Symmetric2x2Scale
// =================
template<typename F>
void
Symmetric2x2Scale
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, Matrix<F>& A, bool conjugate=false );
template<typename F>
void
Symmetric2x2Scale
( LeftOrRight side, UpperOrLower uplo,
  const DistMatrix<F,STAR,STAR>& D, AbstractDistMatrix<F>& A,
  bool conjugate=false );

template<typename F>
void
FirstHalfOfSymmetric2x2Scale
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, Matrix<F>& a1, const Matrix<F>& a2,
  bool conjugate=false );

template<typename F>
void
SecondHalfOfSymmetric2x2Scale
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, const Matrix<F>& a1, Matrix<F>& a2,
  bool conjugate=false );

// Symmetric2x2Solve
// =================
template<typename F>
void
Symmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, Matrix<F>& A, bool conjugate=false );
template<typename F>
void
Symmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const DistMatrix<F,STAR,STAR>& D, AbstractDistMatrix<F>& A,
  bool conjugate=false );

template<typename F>
void
FirstHalfOfSymmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, Matrix<F>& a1, const Matrix<F>& a2,
  bool conjugate=false );

template<typename F>
void
SecondHalfOfSymmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, const Matrix<F>& a1, Matrix<F>& a2,
  bool conjugate=false );

// Transpose
// =========
template<typename T>
void Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate=false );
template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Transpose
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B, bool conjugate=false );

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Transpose
( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,W,Z>& B,
  bool conjugate=false );

template<typename T>
void Transpose
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B,
  bool conjugate=false );
template<typename T>
void Transpose
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B,
  bool conjugate=false );

// UpdateDiagonal
// ==============
template<typename T,typename S>
void UpdateDiagonal( Matrix<T>& A, S alpha, Int offset=0 );
template<typename T,typename S>
void UpdateDiagonal( AbstractDistMatrix<T>& A, S alpha, Int offset=0 );

// Zero
// ====
template<typename T>
void Zero( Matrix<T>& A );
template<typename T>
void Zero( AbstractDistMatrix<T>& A );
template<typename T>
void Zero( AbstractBlockDistMatrix<T>& A );

} // namespace El

#include "./level1/EntrywiseMap.hpp"

#endif // ifndef EL_BLAS1_HPP
