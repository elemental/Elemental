/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPECTRAL_HPP
#define EL_SPECTRAL_HPP

namespace El {

namespace PencilNS {
enum Pencil
{
    AXBX=1,
    ABX=2,
    BAX=3
};
}
using namespace PencilNS;

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
struct HermitianEigSubset
{
    bool indexSubset;
    Int lowerIndex, upperIndex;
 
    bool rangeSubset;
    Real lowerBound, upperBound;

    HermitianEigSubset() 
    : indexSubset(false), lowerIndex(0), upperIndex(0), 
      rangeSubset(false), lowerBound(0), upperBound(0) { }
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
    bool distAed;
    Int blockHeight, blockWidth;

    HessQrCtrl() 
    : distAed(false), 
      blockHeight(DefaultBlockHeight()), blockWidth(DefaultBlockWidth()) 
    { }
};

// Forward declaration
template<typename Real> struct SignCtrl;

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
// Compute eigenvalues
// -------------------
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,
  DistMatrix<Base<F>,STAR,STAR>& w, SortType sort=ASCENDING,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& w,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );

// Compute eigenpairs
// ------------------
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,
  DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& w, 
  AbstractDistMatrix<F>& Z, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );

// Hermitian generalized definite eigenvalue solvers
// =================================================
// Compute eigenvalues
// -------------------
template<typename F>
void HermitianGenDefEig
( Pencil pencil, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianGenDefEig
( Pencil pencil, UpperOrLower uplo,
  AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B,
  AbstractDistMatrix<Base<F>>& w, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
// Compute eigenpairs
// ------------------
template<typename F>
void HermitianGenDefEig
( Pencil pencil, UpperOrLower uplo,
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, Matrix<F>& X,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void HermitianGenDefEig
( Pencil pencil, UpperOrLower uplo,
  AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B,
  AbstractDistMatrix<Base<F>>& w, AbstractDistMatrix<F>& X,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() );

// Hermitian tridiagonal eigenvalue solvers
// ========================================
// Compute eigenvalues
// --------------------
template<typename F>
void HermitianTridiagEig
( Matrix<Base<F>>& d, Matrix<F>& e, Matrix<Base<F>>& w, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>() );
template<typename F>
void HermitianTridiagEig
( const AbstractDistMatrix<Base<F>>& d, const AbstractDistMatrix<F>& e,
        AbstractDistMatrix<Base<F>>& w, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>() );
// Compute eigenpairs
// ------------------
template<typename F>
void HermitianTridiagEig
( Matrix<Base<F>>& d, Matrix<F>& e, Matrix<Base<F>>& w, Matrix<F>& Z,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>() );
template<typename F>
void HermitianTridiagEig
( const AbstractDistMatrix<Base<F>>& d, const AbstractDistMatrix<F>& e,
        AbstractDistMatrix<Base<F>>& w,       AbstractDistMatrix<F>& Z, 
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>() );

template<typename Real>
Int HermitianTridiagEigEstimate
( const AbstractDistMatrix<Real>& d, const AbstractDistMatrix<Real>& e,
        mpi::Comm wColComm, Real vl, Real vu );
// Z is assumed to be sufficiently large and properly aligned
template<typename Real>
void HermitianTridiagEigPostEstimate
( const AbstractDistMatrix<Real>& d, const AbstractDistMatrix<Real>& e,
        AbstractDistMatrix<Real>& w,       AbstractDistMatrix<Real>& Z, 
  SortType sort, Real vl, Real vu );

namespace herm_eig {

template<typename F>
void Sort( Matrix<Base<F>>& w, Matrix<F>& Z, SortType sort=ASCENDING );
template<typename F>
void Sort
( AbstractDistMatrix<Base<F>>& w, AbstractDistMatrix<F>& Z,
  SortType sort=ASCENDING );

} // namespace herm_eig

// Polar decomposition
// ===================
template<typename F>
void Polar( Matrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void Polar( AbstractDistMatrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
void Polar
( Matrix<F>& A, Matrix<F>& P, const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void Polar
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& P, 
  const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
void HermitianPolar
( UpperOrLower uplo, Matrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void HermitianPolar
( UpperOrLower uplo, AbstractDistMatrix<F>& A, 
  const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
void HermitianPolar
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& P, 
  const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void HermitianPolar
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& P, 
  const PolarCtrl& ctrl=PolarCtrl() ); 

// Schur decomposition
// ===================
template<typename F>
void Schur
( Matrix<F>& A, Matrix<Complex<Base<F>>>& w,
  bool fullTriangle=false, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w,
  bool fullTriangle=false, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( BlockDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w,
  bool fullTriangle=false, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );

template<typename F>
void Schur
( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, Matrix<F>& Q,
  bool fullTriangle=true, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w, 
  AbstractDistMatrix<F>& Q, bool fullTriangle=true, 
  const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( BlockDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w, 
  BlockDistMatrix<F>& Q, bool fullTriangle=true, 
  const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );

namespace schur {

template<typename Real>
void CheckRealSchur( const Matrix<Real>& U, bool standardForm=false );
template<typename Real>
void CheckRealSchur
( const AbstractDistMatrix<Real>& U, bool standardForm=false );

// NOTE: These will always throw an error
template<typename Real>
void CheckRealSchur
( const Matrix<Complex<Real>>& U, bool standardForm=false );
template<typename Real>
void CheckRealSchur
( const AbstractDistMatrix<Complex<Real>>& U, bool standardForm=false );

template<typename F>
void QuasiTriangEig
( const Matrix<F>& dMain, const Matrix<F>& dSub, const Matrix<F>& dSup,
  Matrix<Complex<Base<F>>>& w );

template<typename F>
void QuasiTriangEig( const Matrix<F>& U, Matrix<Complex<Base<F>>>& w );
template<typename F>
void QuasiTriangEig
( const AbstractDistMatrix<F>& U, AbstractDistMatrix<Complex<Base<F>>>& w );

template<typename F>
Matrix<Complex<Base<F>>> QuasiTriangEig( const Matrix<F>& U );
template<typename F>
DistMatrix<Complex<Base<F>>,VR,STAR> 
QuasiTriangEig( const AbstractDistMatrix<F>& U );

template<typename Real>
void RealToComplex
( const Matrix<Real>& UQuasi, Matrix<Complex<Real>>& U );
template<typename Real>
void RealToComplex
( const AbstractDistMatrix<Real>& UQuasi, 
        AbstractDistMatrix<Complex<Real>>& U );

} // namespace schur

// Skew-Hermitian eigenvalue solvers
// =================================
// Compute the full set of eigenvalues
// -----------------------------------
template<typename F>
void SkewHermitianEig
( UpperOrLower uplo, const Matrix<F>& G, Matrix<Base<F>>& wImag,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>>& ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void SkewHermitianEig
( UpperOrLower uplo, const AbstractDistMatrix<F>& G,
  AbstractDistMatrix<Base<F>>& wImag, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>>& ctrl=HermitianEigCtrl<Base<F>>() );

// Compute eigenpairs
// ------------------
template<typename F>
void SkewHermitianEig
( UpperOrLower uplo, const Matrix<F>& G,
  Matrix<Base<F>>& wImag, Matrix<Complex<Base<F>>>& Z,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>>& ctrl=HermitianEigCtrl<Base<F>>() );
template<typename F>
void SkewHermitianEig
( UpperOrLower uplo, const AbstractDistMatrix<F>& G,
  AbstractDistMatrix<Base<F>>& wImag, AbstractDistMatrix<Complex<Base<F>>>& Z,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Base<F>>& ctrl=HermitianEigCtrl<Base<F>>() );

// Singular Value Decomposition
// ============================

// Compute the singular values
// ---------------------------
template<typename F>
void SVD( Matrix<F>& A, Matrix<Base<F>>& s );
template<typename F>
void SVD
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, 
  double heightRatio=1.2 );

template<typename F>
void HermitianSVD( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& s );
template<typename F>
void HermitianSVD
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s );

// Compute the full SVD
// --------------------
template<typename F>
void SVD( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, bool useQR=false );
template<typename F>
void SVD
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, 
  AbstractDistMatrix<F>& V, double heightRatio=1.5 );

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& U, Matrix<F>& V );
template<typename F>
void HermitianSVD
( UpperOrLower uplo, AbstractDistMatrix<F>& A,
  AbstractDistMatrix<Base<F>>& s, AbstractDistMatrix<F>& U, 
  AbstractDistMatrix<F>& V );

namespace svd {

template<typename F>
void Thresholded
( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V,
  Base<F> tol=0, bool relative=false );
template<typename F>
void Thresholded
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Base<F>>& s, 
  AbstractDistMatrix<F>& V, Base<F> tol=0, bool relative=false );

template<typename F>
void TallThresholded
( DistMatrix<F,VC,STAR>& A, AbstractDistMatrix<Base<F>>& s,
  AbstractDistMatrix<F>& V, Base<F> tol=0, bool relative=false );

// NOTE: [* ,VR] WideThresholded would produce U with different distribution
//       than A. It makes more sense to overwrite A with V'.

} // namespace svd

} // namespace El

#endif // ifndef EL_SPECTRAL_HPP
