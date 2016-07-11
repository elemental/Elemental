/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SPECTRAL_HPP
#define EL_SPECTRAL_HPP

#include <El/lapack_like/condense.hpp>

namespace El {

template<typename Real>
struct HermitianEigSubset
{
    bool indexSubset=false;
    // The valid index range is [lowerIndex,upperIndex]
    Int lowerIndex=0, upperIndex=0;
 
    bool rangeSubset=false;
    // The valid value range is (lowerBound,upperBound]
    Real lowerBound=Real(0), upperBound=Real(0);
};

// Hermitian tridiagonal eigenvalue solvers
// ========================================

namespace herm_tridiag_eig {

struct QRInfo
{
    Int numUnconverged=0;
    Int numIterations=0;
};

struct QRCtrl
{
    Int maxIterPerEig=30;
    bool demandConverged=true;

    bool fullAccuracyTwoByTwo=true;
};

} // namespace herm_tridiag_eig

struct HermitianTridiagEigInfo
{
    herm_tridiag_eig::QRInfo qrInfo;
    // TODO(poulson): MRRR info
};

template<typename Real,typename=EnableIf<IsReal<Real>>>
struct HermitianTridiagEigCtrl
{
    bool wantEigVecs=false;
    bool accumulateEigVecs=false;
    SortType sort=ASCENDING;
    HermitianEigSubset<Real> subset;
    bool progress=false;

    bool useQR=false;
    herm_tridiag_eig::QRCtrl qrCtrl;
    // TODO(poulson): MRRR ctrl
};

// Compute eigenvalues
// --------------------
template<typename F>
HermitianTridiagEigInfo
HermitianTridiagEig
( Matrix<Base<F>>& d,
  Matrix<F>& dSub,
  Matrix<Base<F>>& w, 
  const HermitianTridiagEigCtrl<Base<F>>& ctrl=
        HermitianTridiagEigCtrl<Base<F>>() );
template<typename F>
HermitianTridiagEigInfo
HermitianTridiagEig
( const AbstractDistMatrix<Base<F>>& d,
  const AbstractDistMatrix<F>& dSub,
        AbstractDistMatrix<Base<F>>& w,
  const HermitianTridiagEigCtrl<Base<F>>& ctrl=
        HermitianTridiagEigCtrl<Base<F>>() );
// Compute eigenpairs
// ------------------
template<typename F>
HermitianTridiagEigInfo
HermitianTridiagEig
( Matrix<Base<F>>& d,
  Matrix<F>& dSub,
  Matrix<Base<F>>& w,
  Matrix<F>& Q,
  const HermitianTridiagEigCtrl<Base<F>>& ctrl=
        HermitianTridiagEigCtrl<Base<F>>() );
template<typename F>
HermitianTridiagEigInfo
HermitianTridiagEig
( const AbstractDistMatrix<Base<F>>& d,
  const AbstractDistMatrix<F>& dSub,
        AbstractDistMatrix<Base<F>>& w,
        AbstractDistMatrix<F>& Q, 
  const HermitianTridiagEigCtrl<Base<F>>& ctrl=
        HermitianTridiagEigCtrl<Base<F>>() );

// Hermitian eigenvalue solvers
// ============================
template<typename Real>
struct HermitianSDCCtrl 
{
    Int cutoff=256;
    Int maxInnerIts=2, maxOuterIts=10;
    Real tol=Real(0);
    Real spreadFactor=Real(1e-6);
    bool progress=false;
};

template<typename F>
struct HermitianEigCtrl
{
    HermitianTridiagCtrl<F> tridiagCtrl;
    HermitianTridiagEigCtrl<Base<F>> tridiagEigCtrl;
    HermitianSDCCtrl<Base<F>> sdcCtrl;
    bool useScaLAPACK=false;
    bool useSDC=false;
    bool timeStages=false;
};

struct HermitianEigInfo
{
    HermitianTridiagEigInfo tridiagEigInfo;
    // TODO(poulson): SDC info
};

// Compute eigenvalues
// -------------------
template<typename F>
HermitianEigInfo
HermitianEig
(       UpperOrLower uplo,
        Matrix<F>& A,
        Matrix<Base<F>>& w,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
HermitianEigInfo
HermitianEig
(       UpperOrLower uplo,
        AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& w,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

// Compute eigenpairs
// ------------------
template<typename F>
HermitianEigInfo
HermitianEig
(       UpperOrLower uplo,
        Matrix<F>& A,
        Matrix<Base<F>>& w,
        Matrix<F>& Q,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
HermitianEigInfo
HermitianEig
(       UpperOrLower uplo,
        AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& w, 
        AbstractDistMatrix<F>& Q,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

namespace herm_eig {

template<typename Real,typename=EnableIf<IsReal<Real>>>
void TwoByTwo
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
  Real& lambda0, Real& lambda1,
  bool fullAccuracy=true );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void TwoByTwo
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
  Real& lambda0, Real& lambda1,
  Real& c, Real& s,
  bool fullAccuracy=true );

} // namespace herm_eig

// Skew-Hermitian eigenvalue solvers
// =================================
// Compute the full set of eigenvalues
// -----------------------------------
template<typename F>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const Matrix<F>& G,
        Matrix<Base<F>>& wImag,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl=
        HermitianEigCtrl<Complex<Base<F>>>() );
template<typename F>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const AbstractDistMatrix<F>& G,
        AbstractDistMatrix<Base<F>>& wImag,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl=
        HermitianEigCtrl<Complex<Base<F>>>() );

// Compute eigenpairs
// ------------------
template<typename F>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const Matrix<F>& G,
        Matrix<Base<F>>& wImag,
        Matrix<Complex<Base<F>>>& Q,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl=
        HermitianEigCtrl<Complex<Base<F>>>() );
template<typename F>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const AbstractDistMatrix<F>& G,
        AbstractDistMatrix<Base<F>>& wImag,
        AbstractDistMatrix<Complex<Base<F>>>& Q,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl=
        HermitianEigCtrl<Complex<Base<F>>>() );

// Hermitian generalized definite eigenvalue solvers
// =================================================
// TODO(poulson): Add support for Fix-Heiberger

namespace PencilNS {
enum Pencil
{
  AXBX=1,
  ABX=2,
  BAX=3
};
}
using namespace PencilNS;

// Compute eigenvalues
// -------------------
template<typename F>
HermitianEigInfo
HermitianGenDefEig
(       Pencil pencil,
        UpperOrLower uplo, 
        Matrix<F>& A,
        Matrix<F>& B,
        Matrix<Base<F>>& w,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
HermitianEigInfo
HermitianGenDefEig
(       Pencil pencil,
        UpperOrLower uplo,
        AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& B,
        AbstractDistMatrix<Base<F>>& w,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
// Compute eigenpairs
// ------------------
template<typename F>
HermitianEigInfo
HermitianGenDefEig
(       Pencil pencil,
        UpperOrLower uplo,
        Matrix<F>& A,
        Matrix<F>& B,
        Matrix<Base<F>>& w,
        Matrix<F>& X,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
HermitianEigInfo
HermitianGenDefEig
(       Pencil pencil,
        UpperOrLower uplo,
        AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& B,
        AbstractDistMatrix<Base<F>>& w,
        AbstractDistMatrix<F>& X,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

// Polar decomposition
// ===================
struct QDWHCtrl
{
    bool colPiv=false;
    Int maxIts=20;
};

struct PolarCtrl 
{
    bool qdwh=false;
    QDWHCtrl qdwhCtrl;
};

struct QDWHInfo
{
    Int numIts=0;
    Int numQRIts=0;
    Int numCholIts=0;
};

struct PolarInfo
{
    QDWHInfo qdwhInfo;
};

template<typename F>
PolarInfo Polar( Matrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
PolarInfo Polar( ElementalMatrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
PolarInfo Polar
( Matrix<F>& A,
  Matrix<F>& P,
  const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
PolarInfo Polar
( ElementalMatrix<F>& A,
  ElementalMatrix<F>& P, 
  const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
PolarInfo HermitianPolar
( UpperOrLower uplo,
  Matrix<F>& A,
  const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
PolarInfo HermitianPolar
( UpperOrLower uplo,
  ElementalMatrix<F>& A, 
  const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
PolarInfo HermitianPolar
( UpperOrLower uplo,
  Matrix<F>& A,
  Matrix<F>& P, 
  const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
PolarInfo HermitianPolar
( UpperOrLower uplo,
  ElementalMatrix<F>& A,
  ElementalMatrix<F>& P, 
  const PolarCtrl& ctrl=PolarCtrl() ); 

// Hessenberg Schur decomposition
// ==============================
struct HessenbergSchurInfo
{
    Int numUnconverged=0;
    Int numIterations=0;
};

namespace hess_schur {
namespace aed {

// Cf. LAPACK's IPARMQ for these choices. The primary difference here is that
// we do not use a fixed value (of 256) for windows of size at least 6000.
inline Int NumShifts( Int n, Int winSize )
{
    Int numShifts;
    if( winSize < 30 )
        numShifts = 2; 
    else if( winSize < 60 )
        numShifts = 4; 
    else if( winSize < 150 )
        numShifts = 10;
    else if( winSize < 590 ) 
        numShifts = Max( 10, winSize/Int(Log2(double(winSize))) );
    else if( winSize < 3000 )
        numShifts = 64;
    else if( winSize < 6000 )
        numShifts = 128;
    else
        numShifts = Max( 256, winSize/Int(2*Log2(double(winSize))) );

    numShifts = Min( numShifts, winSize );
    numShifts = Max( 2, numShifts-Mod(numShifts,2) );

    return numShifts;
}

// Cf. LAPACK's IPARMQ for these deflation window sizes
inline Int DeflationSize( Int n, Int winSize, Int numShifts )
{
    Int deflationSize;
    if( winSize <= 500 )
        deflationSize = numShifts;
    else
        deflationSize = (3*numShifts) / 2;

    deflationSize = Min( deflationSize, winSize );
    deflationSize = Min( deflationSize, (n-1)/3 );
    deflationSize = Max( 2, deflationSize-Mod(deflationSize,2) );

    return deflationSize;
}

// Cf. LAPACK's IPARMQ for the choice of skipping a QR sweep if at least
// 14% of the eigenvalues in a window deflated
inline Int SufficientDeflation( Int deflationSize )
{
    const Int nibble = 14;
    return (nibble*deflationSize) / 100;
}

} // namespace aed
} // namespace hess_schur

struct HessenbergSchurCtrl
{
    Int winBeg=0;
    Int winEnd=END;
    bool fullTriangle=true;
    bool wantSchurVecs=false;
    bool demandConverged=true;

    bool useAED=true;
    bool recursiveAED=true;
    bool accumulateReflections=true;

    bool progress=false;

    // Cf. LAPACK's IPARMQ for this choice;
    // note that LAPACK's hard minimum of 12 does not apply to us
    Int minAEDSize = 75;

    function<Int(Int,Int)> numShifts =
      function<Int(Int,Int)>(hess_schur::aed::NumShifts);

    function<Int(Int,Int,Int)> deflationSize =
      function<Int(Int,Int,Int)>(hess_schur::aed::DeflationSize);

    function<Int(Int)> sufficientDeflation =
      function<Int(Int)>(hess_schur::aed::SufficientDeflation);
};

template<typename Real>
HessenbergSchurInfo
HessenbergSchur
( Matrix<Real>& H,
  Matrix<Complex<Real>>& w,
  const HessenbergSchurCtrl& ctrl=HessenbergSchurCtrl() );
template<typename Real>
HessenbergSchurInfo
HessenbergSchur
( Matrix<Real>& H,
  Matrix<Complex<Real>>& w,
  Matrix<Real>& Z,
  const HessenbergSchurCtrl& ctrl=HessenbergSchurCtrl() );

template<typename Real>
HessenbergSchurInfo
HessenbergSchur
( Matrix<Complex<Real>>& H,
  Matrix<Complex<Real>>& w,
  const HessenbergSchurCtrl& ctrl=HessenbergSchurCtrl() );
template<typename Real>
HessenbergSchurInfo
HessenbergSchur
( Matrix<Complex<Real>>& H,
  Matrix<Complex<Real>>& w,
  Matrix<Complex<Real>>& Z,
  const HessenbergSchurCtrl& ctrl=HessenbergSchurCtrl() );

// Schur decomposition
// ===================
// Forward declaration
template<typename Real> struct SignCtrl;

template<typename Real>
struct SDCCtrl 
{
    Int cutoff=256;
    Int maxInnerIts=2, maxOuterIts=10;
    Real tol=0;
    Real spreadFactor=Real(1e-6);
    bool random=true;
    bool progress=false;

    SignCtrl<Real> signCtrl;
};

// TODO: Combine with HessenbergSchurCtrl
struct HessQRCtrl 
{
    bool distAED=false;
    Int blockHeight=DefaultBlockHeight(), blockWidth=DefaultBlockWidth();
};

template<typename Real>
struct SchurCtrl 
{
    bool useSDC=false;
    HessQRCtrl qrCtrl;
    SDCCtrl<Real> sdcCtrl;    
    bool time=false;
};

template<typename F>
void Schur
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  bool fullTriangle=false,
  const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( ElementalMatrix<F>& A,
  ElementalMatrix<Complex<Base<F>>>& w,
  bool fullTriangle=false,
  const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( DistMatrix<F,MC,MR,BLOCK>& A,
  ElementalMatrix<Complex<Base<F>>>& w,
  bool fullTriangle=false,
  const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );

template<typename F>
void Schur
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  Matrix<F>& Q,
  bool fullTriangle=true,
  const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( ElementalMatrix<F>& A,
  ElementalMatrix<Complex<Base<F>>>& w, 
  ElementalMatrix<F>& Q,
  bool fullTriangle=true, 
  const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( DistMatrix<F,MC,MR,BLOCK>& A,
  ElementalMatrix<Complex<Base<F>>>& w, 
  DistMatrix<F,MC,MR,BLOCK>& Q,
  bool fullTriangle=true, 
  const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );

namespace schur {

template<typename Real>
void CheckRealSchur( const Matrix<Real>& U, bool standardForm=false );
template<typename Real>
void CheckRealSchur
( const ElementalMatrix<Real>& U, bool standardForm=false );

// NOTE: These will always throw an error
template<typename Real>
void CheckRealSchur
( const Matrix<Complex<Real>>& U, bool standardForm=false );
template<typename Real>
void CheckRealSchur
( const ElementalMatrix<Complex<Real>>& U, bool standardForm=false );

template<typename F>
void QuasiTriangEig
( const Matrix<F>& dMain,
  const Matrix<F>& dSub,
  const Matrix<F>& dSup,
  Matrix<Complex<Base<F>>>& w );

template<typename F>
void QuasiTriangEig
( const Matrix<F>& U,
        Matrix<Complex<Base<F>>>& w );
template<typename F>
void QuasiTriangEig
( const ElementalMatrix<F>& U,
        ElementalMatrix<Complex<Base<F>>>& w );

template<typename F>
Matrix<Complex<Base<F>>> QuasiTriangEig( const Matrix<F>& U );
template<typename F>
DistMatrix<Complex<Base<F>>,VR,STAR> 
QuasiTriangEig( const ElementalMatrix<F>& U );

template<typename Real>
void RealToComplex
( const Matrix<Real>& UQuasi, Matrix<Complex<Real>>& U );
template<typename Real>
void RealToComplex
( const ElementalMatrix<Real>& UQuasi, 
        ElementalMatrix<Complex<Real>>& U );

template<typename Real>
void RealToComplex
( const Matrix<Real>& UQuasi,
  const Matrix<Real>& QQuasi,
        Matrix<Complex<Real>>& U,
        Matrix<Complex<Real>>& Q );
template<typename Real>
void RealToComplex
( const ElementalMatrix<Real>& UQuasi, 
  const ElementalMatrix<Real>& QQuasi,
        ElementalMatrix<Complex<Real>>& U,
        ElementalMatrix<Complex<Real>>& Q );

template<typename Real,typename=EnableIf<IsReal<Real>>>
void TwoByTwo
( Real& alpha00, Real& alpha01,
  Real& alpha10, Real& alpha11,
  Complex<Real>& lambda0,
  Complex<Real>& lambda1 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void TwoByTwo
( Real& alpha00, Real& alpha01,
  Real& alpha10, Real& alpha11,
  Complex<Real>& lambda0,
  Complex<Real>& lambda1,
  Real& c, Real& s );

template<typename Real>
void TwoByTwo
( Complex<Real>& alpha00, Complex<Real>& alpha01,
  Complex<Real>& alpha10, Complex<Real>& alpha11,
  Complex<Real>& lambda0, Complex<Real>& lambda1 );

} // namespace schur

// Compute eigenvectors of a triangular matrix
// ===========================================
template<typename F>
void TriangEig
(       Matrix<F>& U,
        Matrix<F>& X );
template<typename F>
void TriangEig
( const ElementalMatrix<F>& U,
        ElementalMatrix<F>& X );

// Compute the eigendecomposition of a square matrix
// =================================================
template<typename F>
void Eig
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  Matrix<Complex<Base<F>>>& X );
template<typename F>
void Eig
( ElementalMatrix<F>& A,
  ElementalMatrix<Complex<Base<F>>>& w,
  ElementalMatrix<Complex<Base<F>>>& X );

// Bidiagonal Singular Value Decomposition
// =======================================

// TODO(poulson): Decide if this should be a separate enum, BidiagSVDApproach
enum SVDApproach
{
  // If A is m x n, return A = U S V^H, where U is m x min(m,n) and 
  // V is n x min(m,n).
  THIN_SVD,

  // If A is m x n and rank k, return A = U S V^H, where U is m x k and 
  // V is n x k.
  COMPACT_SVD,

  // If A is m x n, return A = U S V^H, where U is m x m and V is n x n.
  FULL_SVD,

  // If the sufficiently small singular triplets should be thrown away.
  // When thresholded, a cross-product algorithm is used. This is often
  // advantageous since tridiagonal eigensolvers tend to have faster 
  // parallel implementations than bidiagonal SVD's.
  PRODUCT_SVD
};

enum SingularValueToleranceType
{
  ABSOLUTE_SING_VAL_TOL,
  RELATIVE_TO_MAX_SING_VAL_TOL,
  RELATIVE_TO_SELF_SING_VAL_TOL
};

namespace bidiag_svd {

struct QRInfo
{
    Int numUnconverged=0;
    Int numIterations=0;
};

struct QRCtrl
{
    Int maxIterPerVal=6;
    bool demandConverged=true;

    // See the note above MinSingularValueEstimateOfBidiag
    bool looseMinSingValEst=true;

    bool useFLAME=false;
    bool useLAPACK=false;
};

} // namespace bidiag_svd

struct BidiagSVDInfo
{
    bidiag_svd::QRInfo qrInfo;
};

template<typename Real>
struct BidiagSVDCtrl
{
    bool wantU=true, wantV=true;
    bool accumulateU=false, accumulateV=false;
    SVDApproach approach=THIN_SVD;

    SingularValueToleranceType tolType=RELATIVE_TO_MAX_SING_VAL_TOL;
    Real tol=0; // If zero, the default will be chosen

    bool progress=false;

    bidiag_svd::QRCtrl qrCtrl;
};

namespace bidiag_svd {

// For determining a consistent threshold for setting singular values to zero
template<typename Real>
Real APosterioriThreshold
( Int m, Int n,
  const Real& twoNorm,
  const BidiagSVDCtrl<Real>& ctrl );

} // namespace bidiag_svd

template<typename Real,typename=EnableIf<IsReal<Real>>>
BidiagSVDInfo
BidiagSVD
( Matrix<Real>& mainDiag,
  Matrix<Real>& superDiag,
  Matrix<Real>& s,
  const BidiagSVDCtrl<Real>& ctrl=BidiagSVDCtrl<Real>() );
template<typename F>
BidiagSVDInfo
BidiagSVD
( Matrix<Base<F>>& mainDiag,
  Matrix<Base<F>>& superDiag,
  Matrix<F>& U,
  Matrix<Base<F>>& s,
  Matrix<F>& V,
  const BidiagSVDCtrl<Base<F>>& ctrl=BidiagSVDCtrl<Base<F>>() );

// Singular Value Decomposition
// ============================

struct SVDInfo
{
    BidiagSVDInfo bidiagSVDInfo;
};

template<typename Real>
struct SVDCtrl 
{
    bool overwrite=false; // Allow 'A' to be overwritten computing A = U S V^H
    bool time=false;

    // Use LAPACK within sequential SVD?
    bool useLAPACK=false;

    // Use ScaLAPACK within distributed SVD?
    bool useScaLAPACK=false;

    // Chan's algorithm
    // ----------------

    // The minimum height/width ratio before preprocessing with a QR 
    // decomposition when only computing singular values
    double valChanRatio=1.2;

    // The minimum height/width ratio before preprocessing with a QR
    // decomposition when computing a full SVD
    double fullChanRatio=1.5;

    BidiagSVDCtrl<Real> bidiagSVDCtrl;
};

// Compute the singular values
// ---------------------------
template<typename F>
SVDInfo SVD
(       Matrix<F>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl=SVDCtrl<Base<F>>() );
template<typename F>
SVDInfo SVD
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl=SVDCtrl<Base<F>>() );

template<typename F>
SVDInfo SVD
(       AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl=SVDCtrl<Base<F>>() );
template<typename F>
SVDInfo SVD
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl=SVDCtrl<Base<F>>() );

namespace svd {

template<typename F>
SVDInfo TSQR
(       AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s,
  bool overwrite=false );
template<typename F>
SVDInfo TSQR
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s );

template<typename Real,typename=EnableIf<IsReal<Real>>>
void TwoByTwoUpper
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
        Real& sigmaMax,
        Real& sgnMax,
        Real& sigmaMin,
        Real& sgnMin,
        Real& cU,
        Real& sU,
        Real& cV,
        Real& sV );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void TwoByTwoUpperStandard
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
        Real& sigmaMax,
        Real& sigmaMin,
        Real& cU,
        Real& sU,
        Real& cV,
        Real& sV );

template<typename Real,typename=EnableIf<IsReal<Real>>>
void TwoByTwoUpper
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
        Real& sigmaMax,
        Real& sigmaMin );

} // namespace svd

// Compute the full SVD
// --------------------
template<typename F>
SVDInfo SVD
(       Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl=SVDCtrl<Base<F>>() );
template<typename F>
SVDInfo SVD
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl=SVDCtrl<Base<F>>() );

template<typename F>
SVDInfo SVD
(       AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl=SVDCtrl<Base<F>>() );
template<typename F>
SVDInfo SVD
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl=SVDCtrl<Base<F>>() );

namespace svd {

template<typename F>
SVDInfo TSQR
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s,
        AbstractDistMatrix<F>& V );

} // namespace svd

// Hermitian SVD
// =============

// Compute the singular values
// ---------------------------
template<typename F>
void HermitianSVD
( UpperOrLower uplo,
        Matrix<F>& A,
        Matrix<Base<F>>& s,
  bool overwrite=false );
template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  const Matrix<F>& A,
        Matrix<Base<F>>& s );

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
        ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& s,
  bool overwrite=false );
template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  const ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& s );

// Compute the full SVD
// --------------------
template<typename F>
void HermitianSVD
( UpperOrLower uplo,
        Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  bool overwrite=false );
template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V );

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
        ElementalMatrix<F>& A,
        ElementalMatrix<F>& U, 
        ElementalMatrix<Base<F>>& s,
        ElementalMatrix<F>& V,
  bool overwrite=false );
template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  const ElementalMatrix<F>& A,
        ElementalMatrix<F>& U, 
        ElementalMatrix<Base<F>>& s,
        ElementalMatrix<F>& V );

// Image and kernel
// ================
// Return orthonormal bases for the image and/or kernel of a matrix
// TODO: Provide support for various algorithms (e.g., RRQR vs. SVD)

template<typename F>
void ImageAndKernel
( const Matrix<F>& A,
        Matrix<F>& M,
        Matrix<F>& K );
template<typename F>
void ImageAndKernel
( const ElementalMatrix<F>& A,
        ElementalMatrix<F>& M,
        ElementalMatrix<F>& K );

template<typename F>
void Image
( const Matrix<F>& A,
        Matrix<F>& M );
template<typename F>
void Image
( const ElementalMatrix<F>& A,
        ElementalMatrix<F>& M );

template<typename F>
void Kernel
( const Matrix<F>& A,
        Matrix<F>& K );
template<typename F>
void Kernel
( const ElementalMatrix<F>& A,
        ElementalMatrix<F>& K );

// Lanczos
// =======
// Form the Lanczos decomposition
//
//    A V = V T + v (beta e_{k-1})^H,
//
// where A is an (explicitly) Hermitian matrix.

template<typename F>
void Lanczos
( const SparseMatrix<F>& A,
        Matrix<Base<F>>& T,
        Int basisSize=20 );
template<typename F>
void Lanczos
( const DistSparseMatrix<F>& A,
        ElementalMatrix<Base<F>>& T,
        Int basisSize=20 );

template<typename F>
Base<F> LanczosDecomp
( const SparseMatrix<F>& A,
        Matrix<F>& V,
        Matrix<Base<F>>& T,
        Matrix<F>& v,
        Int basisSize=15 );
template<typename F>
Base<F> LanczosDecomp
( const DistSparseMatrix<F>& A,
        DistMultiVec<F>& V,
        ElementalMatrix<Base<F>>& T,
        DistMultiVec<F>& v,
        Int basisSize=15 );

// Product Lanczos
// ===============
// Form the product Lanczos decomposition
//
//    B V = V T + v (beta e_{k-1})^H,
//
// where B is either A^H A or A A^H, depending upon which is larger.

template<typename F>
void ProductLanczos
( const SparseMatrix<F>& A,
        Matrix<Base<F>>& T,
        Int basisSize=20 );
template<typename F>
void ProductLanczos
( const DistSparseMatrix<F>& A,
        ElementalMatrix<Base<F>>& T,
        Int basisSize=20 );

template<typename F>
Base<F> ProductLanczosDecomp
( const SparseMatrix<F>& A,
        Matrix<F>& V,
        Matrix<Base<F>>& T,
        Matrix<F>& v,
        Int basisSize=15 );
template<typename F>
Base<F> ProductLanczosDecomp
( const DistSparseMatrix<F>& A,
        DistMultiVec<F>& V,
        ElementalMatrix<Base<F>>& T,
        DistMultiVec<F>& v,
        Int basisSize=15 );

// Extremal singular value estimates
// =================================
// Form a product Lanczos decomposition and use the square-roots of the 
// Ritz values as estimates of the extremal singular values.
//
// Note that the minimum singular value, which is the first value returned in
// the pair, is likely to be extremely inaccurate for problems with large
// condition numbers (but the dominant singular value is likely to be accurate).

template<typename F>
pair<Base<F>,Base<F>> 
ExtremalSingValEst( const SparseMatrix<F>& A, Int basisSize=20 );
template<typename F>
pair<Base<F>,Base<F>> 
ExtremalSingValEst( const DistSparseMatrix<F>& A, Int basisSize=20 );

template<typename F>
pair<Base<F>,Base<F>> 
HermitianExtremalSingValEst( const SparseMatrix<F>& A, Int basisSize=20 );
template<typename F>
pair<Base<F>,Base<F>> 
HermitianExtremalSingValEst( const DistSparseMatrix<F>& A, Int basisSize=20 );

// Pseudospectra
// =============
enum PseudospecNorm {
  PS_TWO_NORM,
  PS_ONE_NORM
  /* For now, handle the infinity norm by using the adjoint matrix */
};

// Configurations for how often and what format numerical (num) and image (img)
// snapshots of the pseudospectral estimates should be saved
struct SnapshotCtrl
{
    Int realSize=0, imagSize=0;

    Int imgSaveFreq=-1, numSaveFreq=-1, imgDispFreq=-1;
    Int imgSaveCount=0, numSaveCount=0, imgDispCount=0;
    string imgBase="ps", numBase="ps";
    FileFormat imgFormat=PNG, numFormat=ASCII_MATLAB;
    bool itCounts=true;

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
    PseudospecNorm norm=PS_TWO_NORM;
    Int blockWidth=10; // block width for block 1-norm estimator

    // Preprocessing configuration
    bool schur=true; // begin with reduction to Schur form?
    bool forceComplexSchur=false;
    bool forceComplexPs=false;
    SchurCtrl<Real> schurCtrl;

    // Convergence and deflation criteria
    Int maxIts=50;
    Real tol=Real(1e-6);
    bool deflate=true;

    // (Implicitly Restarted) Arnoldi/Lanczos. If basisSize > 1, then
    // there is implicit restarting
    bool arnoldi=true;
    Int basisSize=10;
    bool reorthog=true; // only matters for IRL, which isn't currently used

    // Whether or not to print progress information at each iteration
    bool progress=false;

    SnapshotCtrl snapCtrl;

    mutable Complex<Real> center = Complex<Real>(0);
    mutable Real realWidth=Real(0), imagWidth=Real(0);
};

template<typename Real>
struct SpectralBox
{
    Complex<Real> center;
    Real realWidth, imagWidth;
};

// (Pseudo-)Spectral portrait
// --------------------------
// Treat each pixel as being located a cell center and tesselate a box with
// said square cells
template<typename F>
Matrix<Int> SpectralPortrait
( const Matrix<F>& A,
        Matrix<Base<F>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<F>>& box,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> SpectralPortrait
( const ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<F>>& box,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularSpectralPortrait
( const Matrix<F>& U,
        Matrix<Base<F>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<F>>& box,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularSpectralPortrait
( const ElementalMatrix<F>& U,
        ElementalMatrix<Base<F>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<F>>& box,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularSpectralPortrait
( const Matrix<F>& U,
  const Matrix<F>& Q,
        Matrix<Base<F>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<F>>& box,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularSpectralPortrait
( const ElementalMatrix<F>& U,
  const ElementalMatrix<F>& Q,
        ElementalMatrix<Base<F>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<F>>& box,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename Real>
Matrix<Int> QuasiTriangularSpectralPortrait
( const Matrix<Real>& U,
        Matrix<Real>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Real>& box,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralPortrait
( const ElementalMatrix<Real>& U,
        ElementalMatrix<Real>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Real>& box,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> QuasiTriangularSpectralPortrait
( const Matrix<Real>& U,
  const Matrix<Real>& Q,
        Matrix<Real>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Real>& box,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralPortrait
( const ElementalMatrix<Real>& U,
  const ElementalMatrix<Real>& Q,
        ElementalMatrix<Real>& invNormMap,
        Int realSize,
        Int imagSize, 
        SpectralBox<Real>& box,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename F>
Matrix<Int> HessenbergSpectralPortrait
( const Matrix<F>& H,
        Matrix<Base<F>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<F>>& box,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergSpectralPortrait
( const ElementalMatrix<F>& H,
        ElementalMatrix<Base<F>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<F>>& box,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> HessenbergSpectralPortrait
( const Matrix<F>& H,
  const Matrix<F>& Q,
        Matrix<Base<F>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<F>>& box,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergSpectralPortrait
( const ElementalMatrix<F>& H,
  const ElementalMatrix<F>& Q,
        ElementalMatrix<Base<F>>& invNormMap,
        Int realSize,
        Int imagSize, 
        SpectralBox<Base<F>>& box,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

// (Pseudo-)Spectral window
// ------------------------
// Treat each pixel as being located a cell center and tesselate a box with
// said square cells
template<typename F>
Matrix<Int> SpectralWindow
( const Matrix<F>& A,
        Matrix<Base<F>>& invNormMap,
        Complex<Base<F>> center,
        Base<F> realWidth,
        Base<F> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> SpectralWindow
( const ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& invNormMap,
        Complex<Base<F>> center,
        Base<F> realWidth,
        Base<F> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularSpectralWindow
( const Matrix<F>& U,
        Matrix<Base<F>>& invNormMap,
        Complex<Base<F>> center,
        Base<F> realWidth,
        Base<F> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularSpectralWindow
( const ElementalMatrix<F>& U,
        ElementalMatrix<Base<F>>& invNormMap,
        Complex<Base<F>> center,
        Base<F> realWidth,
        Base<F> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularSpectralWindow
( const Matrix<F>& U,
  const Matrix<F>& Q,
        Matrix<Base<F>>& invNormMap,
        Complex<Base<F>> center,
        Base<F> realWidth,
        Base<F> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularSpectralWindow
( const ElementalMatrix<F>& U,
  const ElementalMatrix<F>& Q,
        ElementalMatrix<Base<F>>& invNormMap,
        Complex<Base<F>> center,
        Base<F> realWidth,
        Base<F> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename Real>
Matrix<Int> QuasiTriangularSpectralWindow
( const Matrix<Real>& U,
        Matrix<Real>& invNormMap,
        Complex<Real> center,
        Real realWidth,
        Real imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralWindow
( const ElementalMatrix<Real>& U,
        ElementalMatrix<Real>& invNormMap,
        Complex<Real> center,
        Real realWidth,
        Real imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> QuasiTriangularSpectralWindow
( const Matrix<Real>& U,
  const Matrix<Real>& Q,
        Matrix<Real>& invNormMap,
        Complex<Real> center,
        Real realWidth,
        Real imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralWindow
( const ElementalMatrix<Real>& U,
  const ElementalMatrix<Real>& Q,
        ElementalMatrix<Real>& invNormMap,
        Complex<Real> center,
        Real realWidth,
        Real imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename F>
Matrix<Int> HessenbergSpectralWindow
( const Matrix<F>& H,
        Matrix<Base<F>>& invNormMap,
        Complex<Base<F>> center,
        Base<F> realWidth,
        Base<F> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergSpectralWindow
( const ElementalMatrix<F>& H,
        ElementalMatrix<Base<F>>& invNormMap,
        Complex<Base<F>> center,
        Base<F> realWidth,
        Base<F> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> HessenbergSpectralWindow
( const Matrix<F>& H,
  const Matrix<F>& Q,
        Matrix<Base<F>>& invNormMap,
        Complex<Base<F>> center,
        Base<F> realWidth,
        Base<F> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergSpectralWindow
( const ElementalMatrix<F>& H,
  const ElementalMatrix<F>& Q,
        ElementalMatrix<Base<F>>& invNormMap,
        Complex<Base<F>> center,
        Base<F> realWidth,
        Base<F> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

// (Pseudo-)Spectral cloud
// -----------------------
template<typename F>
Matrix<Int> SpectralCloud
( const Matrix<F>& A,
  const Matrix<Complex<Base<F>>>& shifts,
        Matrix<Base<F>>& invNorms,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int,VR,STAR> SpectralCloud
( const ElementalMatrix<F>& A,
  const ElementalMatrix<Complex<Base<F>>>& shifts,
        ElementalMatrix<Base<F>>& invNorms,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularSpectralCloud
( const Matrix<F>& U,
  const Matrix<Complex<Base<F>>>& shifts,
        Matrix<Base<F>>& invNorms,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int,VR,STAR> TriangularSpectralCloud
( const ElementalMatrix<F>& U,
  const ElementalMatrix<Complex<Base<F>>>& shifts,
        ElementalMatrix<Base<F>>& invNorms,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularSpectralCloud
( const Matrix<F>& U,
  const Matrix<F>& Q,
  const Matrix<Complex<Base<F>>>& shifts,
        Matrix<Base<F>>& invNorms,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int,VR,STAR> TriangularSpectralCloud
( const ElementalMatrix<F>& U,
  const ElementalMatrix<F>& Q,
  const ElementalMatrix<Complex<Base<F>>>& shifts,
        ElementalMatrix<Base<F>>& invNorms,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename Real>
Matrix<Int> QuasiTriangularSpectralCloud
( const Matrix<Real>& U,
  const Matrix<Complex<Real>>& shifts,
        Matrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> QuasiTriangularSpectralCloud
( const ElementalMatrix<Real>& U,
  const ElementalMatrix<Complex<Real>>& shifts,
        ElementalMatrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> QuasiTriangularSpectralCloud
( const Matrix<Real>& U,
  const Matrix<Real>& Q,
  const Matrix<Complex<Real>>& shifts,
        Matrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> QuasiTriangularSpectralCloud
( const ElementalMatrix<Real>& U,
  const ElementalMatrix<Real>& Q,
  const ElementalMatrix<Complex<Real>>& shifts,
        ElementalMatrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename F>
Matrix<Int> HessenbergSpectralCloud
( const Matrix<F>& H,
  const Matrix<Complex<Base<F>>>& shifts,
        Matrix<Base<F>>& invNorms,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int,VR,STAR> HessenbergSpectralCloud
( const ElementalMatrix<F>& H,
  const ElementalMatrix<Complex<Base<F>>>& shifts,
        ElementalMatrix<Base<F>>& invNorms,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> HessenbergSpectralCloud
( const Matrix<F>& H,
  const Matrix<F>& Q,
  const Matrix<Complex<Base<F>>>& shifts,
        Matrix<Base<F>>& invNorms,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int,VR,STAR> HessenbergSpectralCloud
( const ElementalMatrix<F>& H,
  const ElementalMatrix<F>& Q,
  const ElementalMatrix<Complex<Base<F>>>& shifts,
        ElementalMatrix<Base<F>>& invNorms,
        PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

} // namespace El

#include <El/lapack_like/spectral/Schur.hpp>
#include <El/lapack_like/spectral/HermitianEig.hpp>
#include <El/lapack_like/spectral/SVD.hpp>
#include <El/lapack_like/spectral/Lanczos.hpp>
#include <El/lapack_like/spectral/ProductLanczos.hpp>

#endif // ifndef EL_SPECTRAL_HPP
