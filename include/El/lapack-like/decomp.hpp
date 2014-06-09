/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DECOMP_HPP
#define EL_DECOMP_HPP

namespace El {

namespace HermitianGenDefiniteEigTypeNS {
enum HermitianGenDefiniteEigType
{
    AXBX=1,
    ABX=2,
    BAX=3
};
}
using namespace HermitianGenDefiniteEigTypeNS;

template<typename Real>
struct HermitianSdcCtrl {
    Int cutoff;
    Int maxInnerIts;
    Int maxOuterIts;
    Real tol;
    Real spreadFactor;
    bool progress;

    HermitianSdcCtrl()
    : cutoff(256), maxInnerIts(2), maxOuterIts(10),
      tol(0), spreadFactor(1e-6),
      progress(false)
    { }
};

template<typename Real>
struct HermitianEigCtrl
{
    HermitianTridiagCtrl tridiagCtrl;
    HermitianSdcCtrl<Real> sdcCtrl;
    bool useSdc;

    HermitianEigCtrl()
    : tridiagCtrl(), sdcCtrl(), useSdc(false)
    { }
};

struct PolarCtrl {
    bool qdwh;
    bool colPiv;
    Int maxIts;
    mutable Int numIts;

    PolarCtrl() : qdwh(false), colPiv(false), maxIts(20), numIts(0) { }
};

struct HessQrCtrl {
    bool aed;
    Int blockHeight, blockWidth;

    HessQrCtrl() 
    : aed(false), 
      blockHeight(DefaultBlockHeight()), blockWidth(DefaultBlockWidth()) 
    { }
};

template<typename Real>
struct SdcCtrl {
    Int cutoff;
    Int maxInnerIts;
    Int maxOuterIts;
    Real tol;
    Real spreadFactor;
    bool random;
    bool progress;

    SignCtrl<Real> signCtrl;

    SdcCtrl()
    : cutoff(256), maxInnerIts(2), maxOuterIts(10),
      tol(0), spreadFactor(1e-6),
      random(true), progress(false), signCtrl()
    { }
};

template<typename Real>
struct SchurCtrl {
    bool useSdc;
    HessQrCtrl qrCtrl;
    SdcCtrl<Real> sdcCtrl;    

    SchurCtrl() : useSdc(false), qrCtrl(), sdcCtrl() { }
};

// Hermitian eigenvalue solvers
// ============================

// Compute the eigenvalues of a Hermitian matrix
// ---------------------------------------------
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,
  DistMatrix<Base<F>,STAR,STAR>& w, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );

// Compute the full eigenvalue decomposition of a Hermitian matrix
// ---------------------------------------------------------------
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,
  DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& paddedZ,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );

// Compute the eigenvalues of a Hermitian matrix within a selected range
// ---------------------------------------------------------------------
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w,
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w,
  Base<F> lowerBound, Base<F> upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A,
  DistMatrix<Base<F>,STAR,STAR>& w,
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A,
  DistMatrix<Base<F>,STAR,STAR>& w,
  Base<F> lowerBound, Base<F> upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,
  Base<F> lowerBound, Base<F> upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );

// Compute a selected set of eigenpairs of a Hermitian matrix
// ----------------------------------------------------------
template<typename F>
void HermitianEig
( UpperOrLower uplo,
  Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z,
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo,
  Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z,
  Base<F> lowerBound, Base<F> upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A,
  DistMatrix<Base<F>,STAR,STAR>& w,
  DistMatrix<F,STAR,STAR>& Z,
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A,
  DistMatrix<Base<F>,STAR,STAR>& w,
  DistMatrix<F,STAR,STAR>& Z,
  Base<F> lowerBound, Base<F> upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& paddedZ,
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& paddedZ,
  Base<F> lowerBound, Base<F> upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );

// Hermitian generalized definite eigenvalue solvers
// =================================================
// Compute the full set of eigenvalues
// -----------------------------------
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
// Compute the full set of eigenpairs
// ----------------------------------
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, Matrix<F>& X,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& X,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
// Compute the eigenvalues within a specified index range
// ------------------------------------------------------
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
// Compute the eigenpairs within a specified index range
// -----------------------------------------------------
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, Matrix<F>& X,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& X,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
// Compute the eigenvalues lying in a particular interval
// ------------------------------------------------------
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w,
  Base<F> a, Base<F> b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<Base<F>,VR,STAR>& w,
  Base<F> a, Base<F> b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
// Compute the eigenpairs with eigenvalues in a particular interval
// ----------------------------------------------------------------
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, Matrix<F>& X,
  Base<F> a, Base<F> b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& X,
  Base<F> a, Base<F> b, SortType sort, const HermitianEigCtrl<Base<F>> ctrl );

// Hermitian tridiagonal eigenvalue solvers
// ========================================
// Compute the full set of eigenvalues
// -----------------------------------
template<typename F>
void HermitianTridiagEig
( Matrix<Base<F>>& d, Matrix<F>& e, Matrix<Base<F>>& w, 
  SortType sort=ASCENDING );
template<typename F,Dist U,Dist V,Dist X>
void HermitianTridiagEig
( const DistMatrix<Base<F>,U,V   >& d,
  const DistMatrix<F,      U,V   >& e,
        DistMatrix<Base<F>,X,STAR>& w, SortType sort=ASCENDING );
// Compute the eigenvalues within the specified index range
// --------------------------------------------------------
template<typename F>
void HermitianTridiagEig
( Matrix<Base<F>>& d, Matrix<F>& e, Matrix<Base<F>>& w,
  Int il, Int iu, SortType sort=ASCENDING );
template<typename F,Dist U,Dist V,Dist X>
void HermitianTridiagEig
( const DistMatrix<Base<F>,U,V   >& d,
  const DistMatrix<F,      U,V   >& e,
        DistMatrix<Base<F>,X,STAR>& w, Int il, Int iu, 
  SortType sort=ASCENDING );
// Compute the eigenvalues within a given interval
// -----------------------------------------------
template<typename F>
void HermitianTridiagEig
( Matrix<Base<F>>& d, Matrix<F>& e, Matrix<Base<F>>& w,
  Base<F> vl, Base<F> vu, SortType sort=ASCENDING );
template<typename F,Dist U,Dist V,Dist X>
void HermitianTridiagEig
( const DistMatrix<Base<F>,U,V   >& d,
  const DistMatrix<F,      U,V   >& e,
        DistMatrix<Base<F>,X,STAR>& w, Base<F> vl, Base<F> vu, 
  SortType sort=ASCENDING );
// Compute the full set of eigenpairs
// ----------------------------------
template<typename F>
void HermitianTridiagEig
( Matrix<Base<F>>& d, Matrix<F>& e, Matrix<Base<F>>& w, Matrix<F>& Z,
  SortType sort=ASCENDING );
template<typename F,Dist U,Dist V,Dist X>
void HermitianTridiagEig
( const DistMatrix<Base<F>,U,   V   >& d,
  const DistMatrix<F,      U,   V   >& e,
        DistMatrix<Base<F>,X,   STAR>& w,
        DistMatrix<F,      STAR,X   >& Z, SortType sort=ASCENDING );
// Compute the eigenpairs within the specified index range
// -------------------------------------------------------
template<typename F>
void HermitianTridiagEig
( Matrix<Base<F>>& d, Matrix<F>& e, Matrix<Base<F>>& w,
  Matrix<F>& Z, Int il, Int iu, SortType sort=ASCENDING );
template<typename F,Dist U,Dist V,Dist X>
void HermitianTridiagEig
( const DistMatrix<Base<F>,U,   V   >& d,
  const DistMatrix<F,      U,   V   >& e,
        DistMatrix<Base<F>,X,   STAR>& w,
        DistMatrix<F,      STAR,X   >& Z, Int il, Int iu, 
  SortType sort=ASCENDING );
// Compute the eigenpairs with eigenvalues within a specified interval
// -------------------------------------------------------------------
template<typename F>
void HermitianTridiagEig
( Matrix<Base<F>>& d, Matrix<F>& e, Matrix<Base<F>>& w,
  Matrix<F>& Z, Base<F> vl, Base<F> vu, SortType sort=ASCENDING );
template<typename F,Dist U,Dist V,Dist X>
void HermitianTridiagEig
( const DistMatrix<Base<F>,U,   V   >& d,
  const DistMatrix<F,      U,   V   >& e,
        DistMatrix<Base<F>,X,   STAR>& w,
        DistMatrix<F,      STAR,X   >& Z,
  Base<F> vl, Base<F> vu, SortType sort=ASCENDING );

template<typename Real,Dist U,Dist V>
Int HermitianTridiagEigEstimate
( const DistMatrix<Real,U,V>& d,
  const DistMatrix<Real,U,V>& e,
        mpi::Comm wColComm, Real vl, Real vu );
// Z is assumed to be sufficiently large and properly aligned
template<typename Real,Dist U,Dist V,Dist X>
void HermitianTridiagEigPostEstimate
( const DistMatrix<Real,U,   V   >& d,
  const DistMatrix<Real,U,   V   >& e,
        DistMatrix<Real,X,   STAR>& w,
        DistMatrix<Real,STAR,X   >& Z, 
  Real vl, Real vu, SortType sort=ASCENDING );

namespace herm_eig {

template<typename F>
void Sort( Matrix<Base<F>>& w, Matrix<F>& Z, SortType sort=ASCENDING );
template<typename F,Dist U1,Dist V1,Dist U2,Dist V2>
void Sort
( DistMatrix<Base<F>,U1,V1>& w, DistMatrix<F,U2,V2>& Z,
  SortType sort=ASCENDING );

} // namespace herm_eig

// Polar decomposition
// ===================
template<typename F>
void Polar( Matrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void Polar( DistMatrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
void Polar
( Matrix<F>& A, Matrix<F>& P, const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void Polar
( DistMatrix<F>& A, DistMatrix<F>& P, const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
void HermitianPolar
( UpperOrLower uplo, Matrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void HermitianPolar
( UpperOrLower uplo, DistMatrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
void HermitianPolar
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& P, 
  const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void HermitianPolar
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& P, 
  const PolarCtrl& ctrl=PolarCtrl() ); 

// Schur decomposition
// ===================
template<typename F>
void Schur
( Matrix<F>& A, Matrix<Complex<Base<F>>>& w,
  bool fullTriangle=true, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( DistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w,
  bool fullTriangle=true, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( BlockDistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w,
  bool fullTriangle=true, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );

template<typename F>
void Schur
( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, Matrix<F>& Q,
  bool fullTriangle=true, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( DistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, DistMatrix<F>& Q,
  bool fullTriangle=true, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( BlockDistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, 
  BlockDistMatrix<F>& Q,
  bool fullTriangle=true, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );

namespace schur {

template<typename Real>
void CheckRealSchur( const Matrix<Real>& U, bool standardForm=false );
template<typename Real>
void CheckRealSchur( const DistMatrix<Real>& U, bool standardForm=false );

// NOTE: These will always throw an error
template<typename Real>
void CheckRealSchur
( const Matrix<Complex<Real>>& U, bool standardForm=false );
template<typename Real>
void CheckRealSchur
( const DistMatrix<Complex<Real>>& U, bool standardForm=false );

template<typename F>
void QuasiTriangEig
( const Matrix<F>& dMain, const Matrix<F>& dSub, const Matrix<F>& dSup,
  Matrix<Complex<Base<F>>>& w );

template<typename F>
void QuasiTriangEig( const Matrix<F>& U, Matrix<Complex<Base<F>>>& w );
template<typename F,Dist colDist,Dist rowDist>
void QuasiTriangEig
( const DistMatrix<F>& U, DistMatrix<Complex<Base<F>>,colDist,rowDist>& w );

template<typename F>
Matrix<Complex<Base<F>>> QuasiTriangEig( const Matrix<F>& U );
template<typename F>
DistMatrix<Complex<Base<F>>,VR,STAR> QuasiTriangEig( const DistMatrix<F>& U );

template<typename Real>
void RealToComplex
( const Matrix<Real>& UQuasi, Matrix<Complex<Real>>& U );
template<typename Real>
void RealToComplex
( const DistMatrix<Real>& UQuasi, DistMatrix<Complex<Real>>& U );

} // namespace schur

// Skew-Hermitian eigenvalue solvers
// =================================
// Compute the full set of eigenvalues
// -----------------------------------
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, Matrix<Real>& wImag,
  SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G,
  Matrix<Real>& wImag, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

// Compute the eigenvalues within the specified index range
// --------------------------------------------------------
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, Matrix<Real>& wImag,
  Int a, Int b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, Int a, Int b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G,
  Matrix<Real>& wImag, Int a, Int b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag, Int a, Int b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

// Compute the eigenvalues which lie within the specified imaginary interval
// -------------------------------------------------------------------------
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, Matrix<Real>& wImag,
  Real a, Real b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, Real a, Real b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G,
  Matrix<Real>& wImag, Real a, Real b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag, Real a, Real b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

// Compute the full set of eigenpairs
// ----------------------------------
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G,
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G,
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

// Compute the eigenpairs within the specified index range
// -------------------------------------------------------
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G,
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  Int a, Int b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  Int a, Int b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G,
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  Int a, Int b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  Int a, Int b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

// Compute eigenpairs with eigenvalues within a specified imaginary interval
// -------------------------------------------------------------------------
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G,
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  Real a, Real b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  Real a, Real b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G,
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  Real a, Real b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );
template<typename Real>
void SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  Real a, Real b, SortType sort=ASCENDING,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() );

// Singular Value Decomposition
// ============================

// Compute the singular values
// ---------------------------
template<typename F>
void SVD( Matrix<F>& A, Matrix<Base<F>>& s );
template<typename F>
void SVD
( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, double heightRatio=1.2 );

template<typename F>
void HermitianSVD( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& s );
template<typename F>
void HermitianSVD
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s );

// Compute the full SVD
// --------------------
template<typename F>
void SVD( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, bool useQR=false );
template<typename F>
void SVD
( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& V,
  double heightRatio=1.5 );

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& U, Matrix<F>& V );
template<typename F>
void HermitianSVD
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& U, DistMatrix<F>& V );

namespace svd {

template<typename F>
void Thresholded
( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V,
  Base<F> tol=0, bool relative=false );
template<typename F>
void Thresholded
( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& V,
  Base<F> tol=0, bool relative=false );

template<typename F>
void TallThresholded
( DistMatrix<F,VC,STAR>& A,
  DistMatrix<Base<F>,STAR,STAR>& s,
  DistMatrix<F,STAR,STAR>& V,
  Base<F> tol=0, bool relative=false );
// NOTE: [* ,VR] WideThresholded would produce U with different distribution
//       than A. It makes more sense to overwrite A with V'.

} // namespace svd

} // namespace El

#endif // ifndef EL_DECOMP_HPP
