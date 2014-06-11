/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PROPS_HPP
#define EL_PROPS_HPP

namespace El {

enum PseudospecNorm {
  PS_TWO_NORM,
  PS_ONE_NORM
  /* For now, handle the infinity norm by using the adjoint matrix */
};

// Configurations for how often and what format numerical (num) and image (img)
// snapshots of the pseudospectral estimates should be saved
struct SnapshotCtrl
{
    Int realSize, imagSize;

    Int imgSaveFreq, numSaveFreq, imgDispFreq;
    Int imgSaveCount, numSaveCount, imgDispCount;
    std::string imgBase, numBase;
    FileFormat imgFormat, numFormat;
    bool itCounts;

    SnapshotCtrl()
    : realSize(0), imagSize(0),
      imgSaveFreq(-1), numSaveFreq(-1), imgDispFreq(-1),
      imgSaveCount(0), numSaveCount(0), imgDispCount(0),
      imgBase("ps"), numBase("ps"), imgFormat(PNG), numFormat(ASCII_MATLAB),
      itCounts(true)
    { }

    void ResetCounts()
    {
        imgSaveCount = 0;
        numSaveCount = 0;
        imgDispCount = 0;
    }
    void Iterate()
    {
        ++imgSaveCount;
        ++numSaveCount;
        ++imgDispCount;
    }
};

template<typename Real>
struct PseudospecCtrl
{
    PseudospecNorm norm;
    Int blockWidth; // block width for block 1-norm estimator

    // Preprocessing configuration
    bool schur; // begin with reduction to Schur form?
    bool forceComplexSchur;
    bool forceComplexPs;
    SchurCtrl<Real> schurCtrl;

    // Convergence and deflation criteria
    Int maxIts;
    Real tol;
    bool deflate;

    // (Implicitly Restarted) Arnoldi/Lanczos. If basisSize > 1, then
    // there is implicit restarting
    bool arnoldi;
    Int basisSize;
    bool reorthog; // only matters for IRL, which isn't currently used

    // Whether or not to print progress information at each iteration
    bool progress;

    SnapshotCtrl snapCtrl;

    PseudospecCtrl()
    : norm(PS_TWO_NORM), blockWidth(10),
      schur(true), forceComplexSchur(false), forceComplexPs(false), schurCtrl(),
      maxIts(200), tol(1e-6), deflate(true),
      arnoldi(true), basisSize(10), reorthog(true),
      progress(false), snapCtrl()
    { }
};

// Condition number
// ================
template<typename F>
Base<F> Condition( const Matrix<F>& A, NormType type=TWO_NORM );
template<typename F,Dist U,Dist V>
Base<F> Condition( const DistMatrix<F,U,V>& A, NormType type=TWO_NORM );

template<typename F>
Base<F> FrobeniusCondition( const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> FrobeniusCondition( const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> InfinityCondition( const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> InfinityCondition( const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> MaxCondition( const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> MaxCondition( const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> OneCondition( const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> OneCondition( const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> TwoCondition( const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> TwoCondition( const DistMatrix<F,U,V>& A );

// Determinant
// ===========
template<typename F>
SafeProduct<F> SafeDeterminant( const Matrix<F>& A );
template<typename F>
SafeProduct<F> SafeDeterminant( const DistMatrix<F>& A );
template<typename F>
SafeProduct<F> SafeDeterminant( Matrix<F>& A, bool canOverwrite=false );
template<typename F>
SafeProduct<F> SafeDeterminant( DistMatrix<F>& A, bool canOverwrite=false );

template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, const DistMatrix<F>& A );
template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false );
template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false );

template<typename F>
F Determinant( const Matrix<F>& A );
template<typename F>
F Determinant( const DistMatrix<F>& A );
template<typename F>
F Determinant( Matrix<F>& A, bool canOverwrite=false );
template<typename F>
F Determinant( DistMatrix<F>& A, bool canOverwrite=false );

template<typename F>
Base<F> HPDDeterminant( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A );
template<typename F>
Base<F> HPDDeterminant
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false );
template<typename F>
Base<F> HPDDeterminant
( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false );

namespace hpd_det {

template<typename F>
SafeProduct<Base<F>> AfterCholesky( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
SafeProduct<Base<F>> AfterCholesky( UpperOrLower uplo, const DistMatrix<F>& A );

} // namespace hpd_det

namespace det {

template<typename F>
SafeProduct<F> AfterLUPartialPiv
( const Matrix<F>& A, const Matrix<Int>& pPerm );
template<typename F,Dist UPerm>
SafeProduct<F> AfterLUPartialPiv
( const DistMatrix<F>& A, const DistMatrix<Int,UPerm,STAR>& pPerm );

} // namespace det

// Inertia
// =======
template<typename F>
InertiaType Inertia
( UpperOrLower uplo, Matrix<F>& A, LDLPivotType pivotType=BUNCH_PARLETT );
template<typename F>
InertiaType Inertia
( UpperOrLower uplo, DistMatrix<F>& A, LDLPivotType pivotType=BUNCH_PARLETT );

// Norm
// ====
template<typename F>
Base<F> Norm( const Matrix<F>& A, NormType type=FROBENIUS_NORM );
template<typename F,Dist U,Dist V>
Base<F> Norm( const DistMatrix<F,U,V>& A, NormType type=FROBENIUS_NORM );

template<typename F>
Base<F> SymmetricNorm
( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM );
template<typename F,Dist U,Dist V>
Base<F> SymmetricNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, 
  NormType type=FROBENIUS_NORM );

template<typename F>
Base<F> HermitianNorm
( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM );
template<typename F,Dist U,Dist V>
Base<F> HermitianNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, 
  NormType type=FROBENIUS_NORM );

// Entrywise norm
// --------------
template<typename F>
Base<F> EntrywiseNorm( const Matrix<F>& A, Base<F> p );
template<typename F>
Base<F> EntrywiseNorm( const AbstractDistMatrix<F>& A, Base<F> p );

template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p );
template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p );

template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p );
template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p );

// Entrywise one-norm
// ------------------
template<typename F>
Base<F> EntrywiseOneNorm( const Matrix<F>& A );
template<typename F>
Base<F> EntrywiseOneNorm( const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> HermitianEntrywiseOneNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianEntrywiseOneNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> SymmetricEntrywiseOneNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricEntrywiseOneNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Frobenius norm
// --------------
template<typename F>
Base<F> FrobeniusNorm( const Matrix<F>& A );
template<typename F>
Base<F> FrobeniusNorm( const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> HermitianFrobeniusNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Infinity norm
// -------------
template<typename F>
Base<F> InfinityNorm( const Matrix<F>& A );
template<typename F>
Base<F> InfinityNorm( const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> HermitianInfinityNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianInfinityNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> SymmetricInfinityNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricInfinityNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Ky-Fan norms
// ------------
template<typename F>
Base<F> KyFanNorm( const Matrix<F>& A, Int k );
template<typename F,Dist U,Dist V>
Base<F> KyFanNorm( const DistMatrix<F,U,V>& A, Int k );

template<typename F>
Base<F> HermitianKyFanNorm
( UpperOrLower uplo, const Matrix<F>& A, Int k );
template<typename F,Dist U,Dist V>
Base<F> HermitianKyFanNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Int k );

template<typename F>
Base<F> SymmetricKyFanNorm
( UpperOrLower uplo, const Matrix<F>& A, Int k );
template<typename F,Dist U,Dist V>
Base<F> SymmetricKyFanNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Int k );

// Max norm
// --------
template<typename F>
Base<F> MaxNorm( const Matrix<F>& A );
template<typename F>
Base<F> MaxNorm( const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> HermitianMaxNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianMaxNorm( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> SymmetricMaxNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricMaxNorm( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Nuclear norm
// ------------
template<typename F>
Base<F> NuclearNorm( const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> NuclearNorm( const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> HermitianNuclearNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> HermitianNuclearNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> SymmetricNuclearNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> SymmetricNuclearNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A );

// One norm
// --------
template<typename F>
Base<F> OneNorm( const Matrix<F>& A );
template<typename F>
Base<F> OneNorm( const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> HermitianOneNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianOneNorm( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> SymmetricOneNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricOneNorm( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Schatten norm
// -------------
template<typename F>
Base<F> SchattenNorm( const Matrix<F>& A, Base<F> p );
template<typename F,Dist U,Dist V>
Base<F> SchattenNorm( const DistMatrix<F,U,V>& A, Base<F> p );

template<typename F>
Base<F> HermitianSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p );
template<typename F,Dist U,Dist V>
Base<F> HermitianSchattenNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Base<F> p );

template<typename F>
Base<F> SymmetricSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p );
template<typename F,Dist U,Dist V>
Base<F> SymmetricSchattenNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Base<F> p );

// Two-norm estimate
// -----------------
template<typename F>
Base<F> TwoNormEstimate
( const Matrix<F>& A, Base<F> tol=1e-6, Int maxIts=1000 );
template<typename F>
Base<F> TwoNormEstimate
( const DistMatrix<F>& A, Base<F> tol=1e-6, Int maxIts=1000 );

template<typename F>
Base<F> HermitianTwoNormEstimate
( UpperOrLower uplo, const Matrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );
template<typename F>
Base<F> HermitianTwoNormEstimate
( UpperOrLower uplo, const DistMatrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );

template<typename F>
Base<F> SymmetricTwoNormEstimate
( UpperOrLower uplo, const Matrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );
template<typename F>
Base<F> SymmetricTwoNormEstimate
( UpperOrLower uplo, const DistMatrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );

// Two norm
// --------
template<typename F>
Base<F> TwoNorm( const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> TwoNorm( const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> HermitianTwoNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> HermitianTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> SymmetricTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A );

// Zero "norm"
// -----------
template<typename F>
Int ZeroNorm( const Matrix<F>& A );
template<typename F>
Int ZeroNorm( const AbstractDistMatrix<F>& A );

// Pseudospectrum
// ==============
template<typename Real>
Matrix<Int> TriangularPseudospectrum
( const Matrix<Complex<Real>>& U, const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> TriangularPseudospectrum
( const DistMatrix<Complex<Real>>& U,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> TriangularPseudospectrum
( const Matrix<Complex<Real>>& U, const Matrix<Complex<Real>>& Q,
  const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> TriangularPseudospectrum
( const DistMatrix<Complex<Real>>& U,
  const DistMatrix<Complex<Real>>& Q,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> TriangularPseudospectrum
( const Matrix<Real>& U, const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> TriangularPseudospectrum
( const DistMatrix<Real>& U, const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> TriangularPseudospectrum
( const Matrix<Real>& U, const Matrix<Real>& Q,
  const Matrix<Complex<Real>>& shifts, Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> TriangularPseudospectrum
( const DistMatrix<Real>& U,
  const DistMatrix<Real>& Q,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U,
  const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U, const Matrix<Real>& Q,
  const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U,
  const DistMatrix<Real>& Q,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<Complex<Real>>& H, const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> HessenbergPseudospectrum
( const DistMatrix<Complex<Real>>& H,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<Complex<Real>>& H,
  const Matrix<Complex<Real>>& Q,
  const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> HessenbergPseudospectrum
( const DistMatrix<Complex<Real>>& H,
  const DistMatrix<Complex<Real>>& Q,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<Real>& H, const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> HessenbergPseudospectrum
( const DistMatrix<Real>& H, const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<Real>& H,
  const Matrix<Real>& Q,
  const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> HessenbergPseudospectrum
( const DistMatrix<Real>& H,
  const DistMatrix<Real>& Q,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> Pseudospectrum
( const Matrix<Real>& A, const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> Pseudospectrum
( const DistMatrix<Real>& A, const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> Pseudospectrum
( const Matrix<Complex<Real>>& A, const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> Pseudospectrum
( const DistMatrix<Complex<Real>>& A,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

// Treat each pixel as being located a cell center and tesselate a box with
// said square cells

template<typename F>
Matrix<Int> TriangularPseudospectrum
( const Matrix<F>& U, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularPseudospectrum
( const DistMatrix<F>& U, DistMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularPseudospectrum
( const Matrix<F>& U, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularPseudospectrum
( const DistMatrix<F>& U, const DistMatrix<F>& Q,
  DistMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename Real>
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U,
  Matrix<Real>& invNormMap,
  Complex<Real> center, Real realWidth, Real imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U,
  DistMatrix<Real>& invNormMap,
  Complex<Real> center, Real realWidth, Real imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U,
  const Matrix<Real>& Q,
  Matrix<Real>& invNormMap,
  Complex<Real> center, Real realWidth, Real imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U,
  const DistMatrix<Real>& Q,
  DistMatrix<Real>& invNormMap,
  Complex<Real> center, Real realWidth, Real imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename F>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<F>& H, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergPseudospectrum
( const DistMatrix<F>& H, DistMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<F>& H, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergPseudospectrum
( const DistMatrix<F>& H,
  const DistMatrix<F>& Q,
  DistMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> Pseudospectrum
( const Matrix<F>& A, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> Pseudospectrum
( const DistMatrix<F>& A, DistMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularPseudospectrum
( const Matrix<F>& U, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularPseudospectrum
( const DistMatrix<F>& U, DistMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularPseudospectrum
( const Matrix<F>& U, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularPseudospectrum
( const DistMatrix<F>& U, const DistMatrix<F>& Q, 
  DistMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename Real>
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U,
  Matrix<Real>& invNormMap, Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U,
  DistMatrix<Real>& invNormMap,
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U, const Matrix<Real>& Q,
  Matrix<Real>& invNormMap,
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U,
  const DistMatrix<Real>& Q,
  DistMatrix<Real>& invNormMap,
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename F>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<F>& H, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergPseudospectrum
( const DistMatrix<F>& H, DistMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<F>& H, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergPseudospectrum
( const DistMatrix<F>& H,
  const DistMatrix<F>& Q, DistMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename Real>
Matrix<Int> Pseudospectrum
( const Matrix<Complex<Real>>& A, Matrix<Real>& invNormMap,
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> Pseudospectrum
( const DistMatrix<Complex<Real>>& A, DistMatrix<Real>& invNormMap,
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> Pseudospectrum
( const Matrix<Real>& A, Matrix<Real>& invNormMap,
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> Pseudospectrum
( const DistMatrix<Real>& A, DistMatrix<Real>& invNormMap,
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

// Trace
// =====
template<typename F>
F Trace( const Matrix<F>& A );
template<typename F>
F Trace( const DistMatrix<F>& A );

} // namespace El

#endif // ifndef EL_PROPS_HPP
