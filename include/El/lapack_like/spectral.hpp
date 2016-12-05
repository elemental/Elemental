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

// Cubic Secular
// =============

// Solve for an inner root of the secular equation
//
//   f(x) = rho + z(0) / (d(0)-x) + z(1) / (d(1)-x) + z(2) / (d(2)-x),
//
// where each numerator is positive and d(0) < d(1) < d(2).
//
// Just as in LAPACK's {s,d}laed6 [CITATION], we require that the user pass in
// an accurate evaluation of f(0).
//

struct CubicSecularInfo
{
    Int numIterations = 0;
    bool converged = true;
};

struct CubicSecularCtrl
{
    Int maxIterations = 40; // Cf. LAPACK's {s,d}laed6 for this choice
    FlipOrClip negativeFix = CLIP_NEGATIVES;
};

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
CubicSecularInfo
CubicSecular
( bool initialize,
  bool rightRoot,
  const Real& rho,
  const Matrix<Real>& z,
  const Matrix<Real>& d,
  const Real& originEval,
        Real& root,
  const CubicSecularCtrl& ctrl=CubicSecularCtrl() );

// Secular Eigenvalue Decomposition
// ================================

struct SecularEVDInfo
{
    Int numIterations = 0;
    Int numAlternations = 0;

    Int numCubicIterations = 0;
    Int numCubicFailures = 0;

    Int numDeflations=0;
    Int numCloseDiagonalDeflations=0;
    Int numSmallUpdateDeflations=0;
};

template<typename Real>
struct SecularEVDCtrl
{
    Int maxIterations = 40; // Cf. LAPACK's {s,d}laed4's choice of 30
    // TODO(poulson): Specialize iteration bounds to grows with precision

    Real sufficientDecay = Real(1)/Real(10);
    FlipOrClip negativeFix = CLIP_NEGATIVES;

    // Incorporate the derivative of the secular function times the absolute
    // value of the current relative estimate of the root of the secular
    // equation into the relative error bound? LAPACK incorporates this when
    // solving for eigenvalues but *not* when solving for singular values.
    bool penalizeDerivative = true;

    bool progress = false;

    CubicSecularCtrl cubicCtrl;
};

// Compute a single eigenvalue corresponding to the diagonal plus rank-one
// matrix
//
//     diag(d) + rho z z^T,
//
// where || z ||_2 = 1, with
//
//     d(0) < d(1) < ... < d(n-1)
//
// and rho > 0.
//
// This routine loosely corresponds to LAPACK's {s,d}laed4 [CITATION].
//
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
SecularEVDInfo
SecularEigenvalue
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Real& eigenvalue,
  const SecularEVDCtrl<Real>& ctrl=SecularEVDCtrl<Real>() );
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
SecularEVDInfo
SecularEigenvalue
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Real& eigenvalue,
        Matrix<Real>& dMinusShift,
  const SecularEVDCtrl<Real>& ctrl=SecularEVDCtrl<Real>() );

// Note that this routine requires that d(0) <= d(1) <= ... <= d(n-1) and
// that || z ||_2 = 1.
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
SecularEVDInfo
SecularEVD
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Matrix<Real>& w,
        Matrix<Real>& Q,
  const SecularEVDCtrl<Real>& ctrl=SecularEVDCtrl<Real>() );

// Secular Singular Value Decomposition
// ====================================

struct SecularSVDInfo
{
    Int numIterations = 0;
    Int numAlternations = 0;
    Int numCubicIterations = 0;
    Int numCubicFailures = 0;

    Int numDeflations=0;
    Int numSmallDiagonalDeflations=0;
    Int numCloseDiagonalDeflations=0;
    Int numSmallUpdateDeflations=0;
};

template<typename Real>
struct SecularSVDCtrl
{
    Int maxIterations = 400; // Cf. LAPACK's {s,d}lasd4 for this choice
    // TODO(poulson): Specialize iteration bounds to grows with precision

    Real sufficientDecay = Real(1)/Real(10);
    FlipOrClip negativeFix = CLIP_NEGATIVES;

    // Incorporate the derivative of the secular function times the absolute
    // value of the current relative estimate of the root of the secular
    // equation into the relative error bound? LAPACK incorporates this when
    // solving for eigenvalues but *not* when solving for singular values.
    bool penalizeDerivative = false;

    bool progress = false;

    CubicSecularCtrl cubicCtrl;
};

// Compute a single singular value corresponding to the square-root of an
// eigenvalue of the diagonal plus rank-one matrix
//
//     diag(d)^2 + rho z z^T,
//
// where || z ||_2 = 1, with
//
//     0 <= d(0) < d(1) < ... < d(n-1)
//
// and rho > 0. In the important case where d(0) = 0, we can build a
// matrix
//
//   M = | sqrt(rho)*z(0), sqrt(rho)*z(1), ..., sqrt(rho)*z(n-1) |
//       |                      d(1),                  .         |
//       |                                  .          .         |
//       |                                           d(n-1)      |
//
// which has said singular values.
//
// This routine loosely corresponds to LAPACK's {s,d}lasd4 [CITATION].
//
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
SecularSVDInfo
SecularSingularValue
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Real& singularValue,
  const SecularSVDCtrl<Real>& ctrl=SecularSVDCtrl<Real>() );
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
SecularSVDInfo
SecularSingularValue
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Real& singularValue,
        Matrix<Real>& dMinusShift,
        Matrix<Real>& dPlusShift,
  const SecularSVDCtrl<Real>& ctrl=SecularSVDCtrl<Real>() );

// Note that this routine requires that 0 = d(0) <= d(1) <= ... <= d(n-1) and
// that || z ||_2 = 1.
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
SecularSVDInfo
SecularSVD
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Matrix<Real>& U,
        Matrix<Real>& s,
        Matrix<Real>& V,
  const SecularSVDCtrl<Real>& ctrl=SecularSVDCtrl<Real>() );


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
    // TODO(poulson): Extend similar to bidiag_svd::QRInfo
};

struct QRCtrl
{
    Int maxIterPerEig=30;
    bool demandConverged=true;

    bool fullAccuracyTwoByTwo=true;

    // As noted in
    //
    //   Beresford N. Parlett and Jian Le,
    //   "Forward instability of Tridiagonal QR",
    //   SIAM J. Matrix Anal. & Appl., 14(1), 279--316, 1991.
    //
    // [CITATION], the tridiagonal QR algorithm is not forward stable.
    // Thus, small differences in the floating-point arithmetic on different
    // processes could lead to substantial differences in the computed
    // eigenvectors. While this is a reasonably well-known issue for the
    // Hessenberg Schur decomposition, it seems to occur much less frequently
    // for tridiagonal matrices.
    //
    bool broadcast=true;
};

struct DCInfo
{
    // TODO(poulson): Extend with more information. For example, the maximum
    // number of levels in the tree and/or the information on the QR iteration
    // at the leaves?
    SecularEVDInfo secularInfo;
};

template<typename Real>
struct DCCtrl
{
    SecularEVDCtrl<Real> secularCtrl;

    // Cf. LAPACK's {s,d}laed2 [CITATION] for the choice of Gu/Eisenstat's
    // [CITATION] "tau" as 8.
    Real deflationFudge = Real(8);

    // Stop recursing when the height is at most 'cutoff'
    Int cutoff = 60;

    // Exploit the nonzero structure of Q when composing the secular
    // eigenvectors with the outer singular vectors? This should only be
    // disabled for academic reasons.
    bool exploitStructure = true;
};

// Cf. Section 4 of Gu and Eisenstat's "A Divide-and-Conquer Algorithm for the
// Bidiagonal SVD" [CITATION] and LAPACK's {s,d}lasd2 [CITATION].
//
// We operationalize Gu and Eisenstat's [CITATION] deflation-tracking
// mechanism by initializing the tags for the nonzero structure of the
// columns of the singular vectors:
//
//   0: nonzero in first block
//   1: nonzero in second block
//   2: dense
//   3: deflated
//
// Cf. LAPACK's {s,d}laed2 [CITATION] for this mechanism.
//
enum DCCombinedColumnType {
  COLUMN_NONZERO_IN_FIRST_BLOCK = 0,
  COLUMN_NONZERO_IN_SECOND_BLOCK = 1,
  DENSE_COLUMN = 2,
  DEFLATED_COLUMN = 3
};
const Int NUM_DC_COMBINED_COLUMN_TYPES = 4;

} // namespace herm_tridiag_eig

struct HermitianTridiagEigInfo
{
    herm_tridiag_eig::QRInfo qrInfo;
    herm_tridiag_eig::DCInfo dcInfo;
    // TODO(poulson): MRRR info
};

enum HermitianTridiagEigAlg {
  HERM_TRIDIAG_EIG_QR = 0,
  HERM_TRIDIAG_EIG_DC = 1,
  HERM_TRIDIAG_EIG_MRRR = 2
};

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
struct HermitianTridiagEigCtrl
{
    bool wantEigVecs=true;
    bool accumulateEigVecs=false;
    SortType sort=ASCENDING;
    HermitianEigSubset<Real> subset;
    bool progress=false;

    HermitianTridiagEigAlg alg=HERM_TRIDIAG_EIG_MRRR;
    herm_tridiag_eig::QRCtrl qrCtrl;
    herm_tridiag_eig::DCCtrl<Real> dcCtrl;
    // TODO(poulson): MRRR ctrl
};

// Compute eigenvalues
// --------------------
template<typename Field>
HermitianTridiagEigInfo
HermitianTridiagEig
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& dSub,
        Matrix<Base<Field>>& w,
  const HermitianTridiagEigCtrl<Base<Field>>& ctrl=
        HermitianTridiagEigCtrl<Base<Field>>() );
template<typename Field>
HermitianTridiagEigInfo
HermitianTridiagEig
( const AbstractDistMatrix<Base<Field>>& d,
  const AbstractDistMatrix<Field>& dSub,
        AbstractDistMatrix<Base<Field>>& w,
  const HermitianTridiagEigCtrl<Base<Field>>& ctrl=
        HermitianTridiagEigCtrl<Base<Field>>() );
// Compute eigenpairs
// ------------------
template<typename Field>
HermitianTridiagEigInfo
HermitianTridiagEig
( const Matrix<Base<Field>>& d,
  const Matrix<Field>& dSub,
        Matrix<Base<Field>>& w,
        Matrix<Field>& Q,
  const HermitianTridiagEigCtrl<Base<Field>>& ctrl=
        HermitianTridiagEigCtrl<Base<Field>>() );
template<typename Field>
HermitianTridiagEigInfo
HermitianTridiagEig
( const AbstractDistMatrix<Base<Field>>& d,
  const AbstractDistMatrix<Field>& dSub,
        AbstractDistMatrix<Base<Field>>& w,
        AbstractDistMatrix<Field>& Q,
  const HermitianTridiagEigCtrl<Base<Field>>& ctrl=
        HermitianTridiagEigCtrl<Base<Field>>() );

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

template<typename Field>
struct HermitianEigCtrl
{
    HermitianTridiagCtrl<Field> tridiagCtrl;
    HermitianTridiagEigCtrl<Base<Field>> tridiagEigCtrl;
    HermitianSDCCtrl<Base<Field>> sdcCtrl;
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
template<typename Field>
HermitianEigInfo
HermitianEig
(       UpperOrLower uplo,
        Matrix<Field>& A,
        Matrix<Base<Field>>& w,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );
template<typename Field>
HermitianEigInfo
HermitianEig
(       UpperOrLower uplo,
        AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& w,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );

// Compute eigenpairs
// ------------------
template<typename Field>
HermitianEigInfo
HermitianEig
(       UpperOrLower uplo,
        Matrix<Field>& A,
        Matrix<Base<Field>>& w,
        Matrix<Field>& Q,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );
template<typename Field>
HermitianEigInfo
HermitianEig
(       UpperOrLower uplo,
        AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& w,
        AbstractDistMatrix<Field>& Q,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );

namespace herm_eig {

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void TwoByTwo
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
  Real& lambda0, Real& lambda1,
  bool fullAccuracy=true );
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
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
template<typename Field>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const Matrix<Field>& G,
        Matrix<Base<Field>>& wImag,
  const HermitianEigCtrl<Complex<Base<Field>>>& ctrl=
        HermitianEigCtrl<Complex<Base<Field>>>() );
template<typename Field>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& G,
        AbstractDistMatrix<Base<Field>>& wImag,
  const HermitianEigCtrl<Complex<Base<Field>>>& ctrl=
        HermitianEigCtrl<Complex<Base<Field>>>() );

// Compute eigenpairs
// ------------------
template<typename Field>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const Matrix<Field>& G,
        Matrix<Base<Field>>& wImag,
        Matrix<Complex<Base<Field>>>& Q,
  const HermitianEigCtrl<Complex<Base<Field>>>& ctrl=
        HermitianEigCtrl<Complex<Base<Field>>>() );
template<typename Field>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& G,
        AbstractDistMatrix<Base<Field>>& wImag,
        AbstractDistMatrix<Complex<Base<Field>>>& Q,
  const HermitianEigCtrl<Complex<Base<Field>>>& ctrl=
        HermitianEigCtrl<Complex<Base<Field>>>() );

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
template<typename Field>
HermitianEigInfo
HermitianGenDefEig
(       Pencil pencil,
        UpperOrLower uplo,
        Matrix<Field>& A,
        Matrix<Field>& B,
        Matrix<Base<Field>>& w,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );
template<typename Field>
HermitianEigInfo
HermitianGenDefEig
(       Pencil pencil,
        UpperOrLower uplo,
        AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B,
        AbstractDistMatrix<Base<Field>>& w,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );
// Compute eigenpairs
// ------------------
template<typename Field>
HermitianEigInfo
HermitianGenDefEig
(       Pencil pencil,
        UpperOrLower uplo,
        Matrix<Field>& A,
        Matrix<Field>& B,
        Matrix<Base<Field>>& w,
        Matrix<Field>& X,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );
template<typename Field>
HermitianEigInfo
HermitianGenDefEig
(       Pencil pencil,
        UpperOrLower uplo,
        AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B,
        AbstractDistMatrix<Base<Field>>& w,
        AbstractDistMatrix<Field>& X,
  const HermitianEigCtrl<Field>& ctrl=HermitianEigCtrl<Field>() );

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

template<typename Field>
PolarInfo
Polar( Matrix<Field>& A, const PolarCtrl& ctrl=PolarCtrl() );
template<typename Field>
PolarInfo
Polar( AbstractDistMatrix<Field>& A, const PolarCtrl& ctrl=PolarCtrl() );
template<typename Field>
PolarInfo Polar
( Matrix<Field>& A,
  Matrix<Field>& P,
  const PolarCtrl& ctrl=PolarCtrl() );
template<typename Field>
PolarInfo Polar
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& P,
  const PolarCtrl& ctrl=PolarCtrl() );
template<typename Field>
PolarInfo HermitianPolar
( UpperOrLower uplo,
  Matrix<Field>& A,
  const PolarCtrl& ctrl=PolarCtrl() );
template<typename Field>
PolarInfo HermitianPolar
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& A,
  const PolarCtrl& ctrl=PolarCtrl() );
template<typename Field>
PolarInfo HermitianPolar
( UpperOrLower uplo,
  Matrix<Field>& A,
  Matrix<Field>& P,
  const PolarCtrl& ctrl=PolarCtrl() );
template<typename Field>
PolarInfo HermitianPolar
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& P,
  const PolarCtrl& ctrl=PolarCtrl() );

// Hessenberg Schur decomposition
// ==============================
struct HessenbergSchurInfo
{
    Int numUnconverged=0;
    Int numIterations=0;
};

namespace hess_schur {

namespace multibulge {

inline Int NumBulgesPerBlock( Int blockHeight )
{ return (blockHeight-1) / 6; }

} // namespace multibulge

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

enum HessenbergSchurAlg {
  HESSENBERG_SCHUR_AED=0,
  HESSENBERG_SCHUR_MULTIBULGE=1,
  HESSENBERG_SCHUR_SIMPLE=2
};

struct HessenbergSchurCtrl
{
    Int winBeg=0;
    Int winEnd=END;
    bool fullTriangle=true;
    bool wantSchurVecs=false;
    bool accumulateSchurVecs=false;
    bool demandConverged=true;

    HessenbergSchurAlg alg=HESSENBERG_SCHUR_AED;
    bool recursiveAED=true;
    bool accumulateReflections=true;
    bool sortShifts=true;

    bool progress=false;

    // Cf. LAPACK's IPARMQ for this choice;
    // note that LAPACK's hard minimum of 12 does not apply to us
    Int minMultiBulgeSize = 75;
    Int minDistMultiBulgeSize = 400;

    function<Int(Int,Int)> numShifts =
      function<Int(Int,Int)>(hess_schur::aed::NumShifts);

    function<Int(Int,Int,Int)> deflationSize =
      function<Int(Int,Int,Int)>(hess_schur::aed::DeflationSize);

    function<Int(Int)> sufficientDeflation =
      function<Int(Int)>(hess_schur::aed::SufficientDeflation);

    // For the distributed Hessenberg QR algorithm
    // TODO(poulson): Move this into a substructure?
    bool scalapack=false;
    Int blockHeight=DefaultBlockHeight();
    // A map from the block height to the number of bulges per diagonal block in
    // the distributed multibulge algorithm.
    function<Int(Int)> numBulgesPerBlock =
      function<Int(Int)>(hess_schur::multibulge::NumBulgesPerBlock);
};

template<typename Field>
HessenbergSchurInfo
HessenbergSchur
( Matrix<Field>& H,
  Matrix<Complex<Base<Field>>>& w,
  const HessenbergSchurCtrl& ctrl=HessenbergSchurCtrl() );
template<typename Field>
HessenbergSchurInfo
HessenbergSchur
( Matrix<Field>& H,
  Matrix<Complex<Base<Field>>>& w,
  Matrix<Field>& Z,
  const HessenbergSchurCtrl& ctrl=HessenbergSchurCtrl() );

template<typename Field>
HessenbergSchurInfo
HessenbergSchur
( AbstractDistMatrix<Field>& H,
  AbstractDistMatrix<Complex<Base<Field>>>& w,
  const HessenbergSchurCtrl& ctrl=HessenbergSchurCtrl() );
template<typename Field>
HessenbergSchurInfo
HessenbergSchur
( AbstractDistMatrix<Field>& H,
  AbstractDistMatrix<Complex<Base<Field>>>& w,
  AbstractDistMatrix<Field>& Z,
  const HessenbergSchurCtrl& ctrl=HessenbergSchurCtrl() );

namespace hess_schur {

template<typename Field>
void Sweep
( Matrix<Field>& H,
  Matrix<Complex<Base<Field>>>& shifts,
  Matrix<Field>& Z,
  const HessenbergSchurCtrl& ctrl );
// TODO(poulson): Generalize to AbstractDistMatrix?
template<typename Field>
void Sweep
( DistMatrix<Field,MC,MR,BLOCK>& H,
  DistMatrix<Complex<Base<Field>>,STAR,STAR>& shifts,
  DistMatrix<Field,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl );

} // namespace hess_schur

// Schur decomposition
// ===================
// Forward declaration
template<typename Real> struct SignCtrl;

template<typename Real>
struct SDCCtrl
{
    Int cutoff=256;
    Int maxInnerIts=2, maxOuterIts=10;
    Real tol=Real(0);
    Real spreadFactor=Real(1e-6);
    bool random=true;
    bool progress=false;

    SignCtrl<Real> signCtrl;
};

template<typename Real>
struct SchurCtrl
{
    bool useSDC=false;
    HessenbergSchurCtrl hessSchurCtrl;
    SDCCtrl<Real> sdcCtrl;
    bool time=false;
};

template<typename Field>
void Schur
( Matrix<Field>& A,
  Matrix<Complex<Base<Field>>>& w,
  const SchurCtrl<Base<Field>>& ctrl=SchurCtrl<Base<Field>>() );
template<typename Field>
void Schur
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Complex<Base<Field>>>& w,
  const SchurCtrl<Base<Field>>& ctrl=SchurCtrl<Base<Field>>() );

template<typename Field>
void Schur
( Matrix<Field>& A,
  Matrix<Complex<Base<Field>>>& w,
  Matrix<Field>& Q,
  const SchurCtrl<Base<Field>>& ctrl=SchurCtrl<Base<Field>>() );
template<typename Field>
void Schur
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Complex<Base<Field>>>& w,
  AbstractDistMatrix<Field>& Q,
  const SchurCtrl<Base<Field>>& ctrl=SchurCtrl<Base<Field>>() );

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

template<typename Field>
void QuasiTriangEig
( const Matrix<Field>& dMain,
  const Matrix<Field>& dSub,
  const Matrix<Field>& dSup,
  Matrix<Complex<Base<Field>>>& w );

template<typename Field>
void QuasiTriangEig
( const Matrix<Field>& U,
        Matrix<Complex<Base<Field>>>& w );
template<typename Field>
void QuasiTriangEig
( const AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Complex<Base<Field>>>& w );

template<typename Field>
Matrix<Complex<Base<Field>>> QuasiTriangEig( const Matrix<Field>& U );
template<typename Field>
DistMatrix<Complex<Base<Field>>,VR,STAR>
QuasiTriangEig( const AbstractDistMatrix<Field>& U );

template<typename Real>
void RealToComplex
( const Matrix<Real>& UQuasi, Matrix<Complex<Real>>& U );
template<typename Real>
void RealToComplex
( const AbstractDistMatrix<Real>& UQuasi,
        AbstractDistMatrix<Complex<Real>>& U );

template<typename Real>
void RealToComplex
( const Matrix<Real>& UQuasi,
  const Matrix<Real>& QQuasi,
        Matrix<Complex<Real>>& U,
        Matrix<Complex<Real>>& Q );
template<typename Real>
void RealToComplex
( const AbstractDistMatrix<Real>& UQuasi,
  const AbstractDistMatrix<Real>& QQuasi,
        AbstractDistMatrix<Complex<Real>>& U,
        AbstractDistMatrix<Complex<Real>>& Q );

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void TwoByTwo
( Real& alpha00, Real& alpha01,
  Real& alpha10, Real& alpha11,
  Complex<Real>& lambda0,
  Complex<Real>& lambda1 );
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
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
template<typename Field>
void TriangEig
(       Matrix<Field>& U,
        Matrix<Field>& X );
template<typename Field>
void TriangEig
( const AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Field>& X );

// Compute the eigendecomposition of a square matrix
// =================================================
template<typename Field>
void Eig
( Matrix<Field>& A,
  Matrix<Complex<Base<Field>>>& w,
  Matrix<Complex<Base<Field>>>& X );
template<typename Field>
void Eig
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Complex<Base<Field>>>& w,
  AbstractDistMatrix<Complex<Base<Field>>>& X );

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
    Int numInnerLoops=0;

    Int numZeroShiftForwardIterations=0;
    Int numZeroShiftForwardInnerLoops=0;

    Int numZeroShiftBackwardIterations=0;
    Int numZeroShiftBackwardInnerLoops=0;

    Int numNonzeroShiftForwardIterations=0;
    Int numNonzeroShiftForwardInnerLoops=0;

    Int numNonzeroShiftBackwardIterations=0;
    Int numNonzeroShiftBackwardInnerLoops=0;
};

struct QRCtrl
{
    Int maxIterPerVal=6;
    bool demandConverged=true;

    // See the note above MinSingularValueEstimateOfBidiag
    bool looseMinSingValEst=true;

    bool useFLAME=false;
    bool useLAPACK=false;

    // As noted in
    //
    //   Beresford N. Parlett and Jian Le,
    //   "Forward instability of Tridiagonal QR",
    //   SIAM J. Matrix Anal. & Appl., 14(1), 279--316, 1991.
    //
    // [CITATION], the tridiagonal QR algorithm is not forward stable.
    // Thus, small differences in the floating-point arithmetic on different
    // processes could lead to substantial differences in the computed
    // singular vectors. While this is a reasonably well-known issue for the
    // Hessenberg Schur decomposition, it seems to occur much less frequently
    // for tridiagonal matrices.
    //
    bool broadcast=true;
};

struct DCInfo
{
    // TODO(poulson): Extend with more information. For example, the maximum
    // number of levels in the tree and/or the information on the QR iteration
    // at the leaves?
    SecularSVDInfo secularInfo;
};

template<typename Real>
struct DCCtrl
{
    SecularSVDCtrl<Real> secularCtrl;

    // Cf. LAPACK's {s,d}lasd2 [CITATION] for the choice of Gu/Eisenstat's
    // [CITATION] "tau" as 8.
    Real deflationFudge = Real(8);

    // Stop recursing when the height is at most 'cutoff'
    Int cutoff = 60;

    // Exploit the nonzero structure of U and V when composing the secular
    // singular vectors with the outer singular vectors? This should only be
    // disabled for academic reasons.
    bool exploitStructure = true;
};

// Cf. Section 4 of Gu and Eisenstat's "A Divide-and-Conquer Algorithm for the
// Bidiagonal SVD" [CITATION] and LAPACK's {s,d}lasd2 [CITATION].
//
// We begin with the decomposition
//
// B = | U_0, 0,  0  | |     diag(s_0),              0       | | V_0, 0   |^T,
//     |  0,  1,  0  | | alpha e_{m_0}^T V_0, beta*e_0^T*V_1 | |   0, V_1 |
//     |  0,  0, U_1 | |         0,              diag(s_1)   |
//
// where U_0 is m_0 x m_0, U_1 is m_1 x m_1, V_0 is (m0+1) x (m0+1), and V_1 is
// either m1 x m1 or (m1+1) x (m1+1). Thus, putting m = m_0 + 1 + m_1, B is
// either m x m or m x (m+1). On entry, U and V should be filled with their
// above depictions.
//
// We operationalize Gu and Eisenstat's [CITATION] deflation-tracking
// mechanism by initializing the tags for the nonzero structure of the
// columns of the singular vectors:
//
//   0: nonzero in first block
//   1: nonzero in second block
//   2: dense
//   3: deflated
//
// Cf. LAPACK's {s,d}lasd2 [CITATION] for this mechanism. Note that LAPACK
// currently ignores deflations of the form |d(0)-d(j)| <= deflationTol,
// which results in the first column of U potentially becoming dense. We
// do not ignore such deflations and always mark the first column of U
// as dense for the sake of simplicity.
//
enum DCCombinedColumnType {
  COLUMN_NONZERO_IN_FIRST_BLOCK = 0,
  COLUMN_NONZERO_IN_SECOND_BLOCK = 1,
  DENSE_COLUMN = 2,
  DEFLATED_COLUMN = 3
};
const Int NUM_DC_COMBINED_COLUMN_TYPES = 4;

} // namespace bidiag_svd

struct BidiagSVDInfo
{
    bidiag_svd::QRInfo qrInfo;
    bidiag_svd::DCInfo dcInfo;
};

template<typename Real>
struct BidiagSVDCtrl
{
    bool wantU=true, wantV=true;
    bool accumulateU=false, accumulateV=false;
    SVDApproach approach=THIN_SVD;

    SingularValueToleranceType tolType=RELATIVE_TO_MAX_SING_VAL_TOL;
    Real tol=Real(0); // If zero, the default will be chosen

    bool progress=false;

    bool useQR=false;
    bidiag_svd::QRCtrl qrCtrl;
    bidiag_svd::DCCtrl<Real> dcCtrl;
};

namespace bidiag_svd {

// For determining a consistent threshold for setting singular values to zero
template<typename Real>
Real APosterioriThreshold
( Int m, Int n,
  const Real& twoNorm,
  const BidiagSVDCtrl<Real>& ctrl );

} // namespace bidiag_svd

template<typename Field>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const Matrix<Field>& mainDiag,
  const Matrix<Field>& offDiag,
        Matrix<Base<Field>>& s,
  const BidiagSVDCtrl<Base<Field>>& ctrl=BidiagSVDCtrl<Base<Field>>() );
template<typename Field>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const Matrix<Field>& mainDiag,
  const Matrix<Field>& offDiag,
        Matrix<Field>& U,
        Matrix<Base<Field>>& s,
        Matrix<Field>& V,
  const BidiagSVDCtrl<Base<Field>>& ctrl=BidiagSVDCtrl<Base<Field>>() );
template<typename Field>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const Matrix<Base<Field>>& mainDiag,
  const Matrix<Base<Field>>& offDiag,
        Matrix<Field>& U,
        Matrix<Base<Field>>& s,
        Matrix<Field>& V,
  const BidiagSVDCtrl<Base<Field>>& ctrl=BidiagSVDCtrl<Base<Field>>() );

template<typename Field>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& mainDiag,
  const AbstractDistMatrix<Field>& offDiag,
        AbstractDistMatrix<Base<Field>>& s,
  const BidiagSVDCtrl<Base<Field>>& ctrl=BidiagSVDCtrl<Base<Field>>() );
template<typename Field>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& mainDiag,
  const AbstractDistMatrix<Field>& offDiag,
        AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& s,
        AbstractDistMatrix<Field>& V,
  const BidiagSVDCtrl<Base<Field>>& ctrl=BidiagSVDCtrl<Base<Field>>() );
template<typename Field>
BidiagSVDInfo
BidiagSVD
( UpperOrLower uplo,
  const AbstractDistMatrix<Base<Field>>& mainDiag,
  const AbstractDistMatrix<Base<Field>>& offDiag,
        AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& s,
        AbstractDistMatrix<Field>& V,
  const BidiagSVDCtrl<Base<Field>>& ctrl=BidiagSVDCtrl<Base<Field>>() );

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
template<typename Field>
SVDInfo SVD
(       Matrix<Field>& A,
        Matrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl=SVDCtrl<Base<Field>>() );
template<typename Field>
SVDInfo SVD
( const Matrix<Field>& A,
        Matrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl=SVDCtrl<Base<Field>>() );

template<typename Field>
SVDInfo SVD
(       AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl=SVDCtrl<Base<Field>>() );
template<typename Field>
SVDInfo SVD
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl=SVDCtrl<Base<Field>>() );

namespace svd {

template<typename Field>
SVDInfo TSQR
(       AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& s,
  bool overwrite=false );
template<typename Field>
SVDInfo TSQR
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& s );

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
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
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
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

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void TwoByTwoUpper
( const Real& alpha00,
  const Real& alpha01,
  const Real& alpha11,
        Real& sigmaMax,
        Real& sigmaMin );

} // namespace svd

// Compute the full SVD
// --------------------
template<typename Field>
SVDInfo SVD
(       Matrix<Field>& A,
        Matrix<Field>& U,
        Matrix<Base<Field>>& s,
        Matrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl=SVDCtrl<Base<Field>>() );
template<typename Field>
SVDInfo SVD
( const Matrix<Field>& A,
        Matrix<Field>& U,
        Matrix<Base<Field>>& s,
        Matrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl=SVDCtrl<Base<Field>>() );

template<typename Field>
SVDInfo SVD
(       AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& s,
        AbstractDistMatrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl=SVDCtrl<Base<Field>>() );
template<typename Field>
SVDInfo SVD
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& s,
        AbstractDistMatrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl=SVDCtrl<Base<Field>>() );

namespace svd {

template<typename Field>
SVDInfo TSQR
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& s,
        AbstractDistMatrix<Field>& V );

} // namespace svd

// Hermitian SVD
// =============

// Compute the singular values
// ---------------------------
template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
        Matrix<Field>& A,
        Matrix<Base<Field>>& s,
  bool overwrite=false );
template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  const Matrix<Field>& A,
        Matrix<Base<Field>>& s );

template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
        AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& s,
  bool overwrite=false );
template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& s );

// Compute the full SVD
// --------------------
template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
        Matrix<Field>& A,
        Matrix<Field>& U,
        Matrix<Base<Field>>& s,
        Matrix<Field>& V,
  bool overwrite=false );
template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  const Matrix<Field>& A,
        Matrix<Field>& U,
        Matrix<Base<Field>>& s,
        Matrix<Field>& V );

template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& U,
  AbstractDistMatrix<Base<Field>>& s,
  AbstractDistMatrix<Field>& V,
  bool overwrite=false );
template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& s,
        AbstractDistMatrix<Field>& V );

// Image and kernel
// ================
// Return orthonormal bases for the image and/or kernel of a matrix
// TODO(poulson): Provide support for various algorithms (e.g., RRQR vs. SVD)

template<typename Field>
void ImageAndKernel
( const Matrix<Field>& A,
        Matrix<Field>& M,
        Matrix<Field>& K );
template<typename Field>
void ImageAndKernel
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& M,
        AbstractDistMatrix<Field>& K );

template<typename Field>
void Image
( const Matrix<Field>& A,
        Matrix<Field>& M );
template<typename Field>
void Image
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& M );

template<typename Field>
void Kernel
( const Matrix<Field>& A,
        Matrix<Field>& K );
template<typename Field>
void Kernel
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& K );

// Lanczos
// =======
// Form the Lanczos decomposition
//
//    A V = V T + v (beta e_{k-1})^H,
//
// where A is an (explicitly) Hermitian matrix.

template<typename Field>
void Lanczos
( const SparseMatrix<Field>& A,
        Matrix<Base<Field>>& T,
        Int basisSize=20 );
template<typename Field>
void Lanczos
( const DistSparseMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& T,
        Int basisSize=20 );

template<typename Field>
Base<Field> LanczosDecomp
( const SparseMatrix<Field>& A,
        Matrix<Field>& V,
        Matrix<Base<Field>>& T,
        Matrix<Field>& v,
        Int basisSize=15 );
template<typename Field>
Base<Field> LanczosDecomp
( const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& V,
        AbstractDistMatrix<Base<Field>>& T,
        DistMultiVec<Field>& v,
        Int basisSize=15 );

// Product Lanczos
// ===============
// Form the product Lanczos decomposition
//
//    B V = V T + v (beta e_{k-1})^H,
//
// where B is either A^H A or A A^H, depending upon which is larger.

template<typename Field>
void ProductLanczos
( const SparseMatrix<Field>& A,
        Matrix<Base<Field>>& T,
        Int basisSize=20 );
template<typename Field>
void ProductLanczos
( const DistSparseMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& T,
        Int basisSize=20 );

template<typename Field>
Base<Field> ProductLanczosDecomp
( const SparseMatrix<Field>& A,
        Matrix<Field>& V,
        Matrix<Base<Field>>& T,
        Matrix<Field>& v,
        Int basisSize=15 );
template<typename Field>
Base<Field> ProductLanczosDecomp
( const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& V,
        AbstractDistMatrix<Base<Field>>& T,
        DistMultiVec<Field>& v,
        Int basisSize=15 );

// Extremal singular value estimates
// =================================
// Form a product Lanczos decomposition and use the square-roots of the
// Ritz values as estimates of the extremal singular values.
//
// Note that the minimum singular value, which is the first value returned in
// the pair, is likely to be extremely inaccurate for problems with large
// condition numbers (but the dominant singular value is likely to be accurate).

template<typename Field>
pair<Base<Field>,Base<Field>>
ExtremalSingValEst
( const SparseMatrix<Field>& A, Int basisSize=20 );
template<typename Field>
pair<Base<Field>,Base<Field>>
ExtremalSingValEst
( const DistSparseMatrix<Field>& A, Int basisSize=20 );

template<typename Field>
pair<Base<Field>,Base<Field>>
HermitianExtremalSingValEst
( const SparseMatrix<Field>& A, Int basisSize=20 );
template<typename Field>
pair<Base<Field>,Base<Field>>
HermitianExtremalSingValEst
( const DistSparseMatrix<Field>& A, Int basisSize=20 );

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
template<typename Field>
Matrix<Int> SpectralPortrait
( const Matrix<Field>& A,
        Matrix<Base<Field>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<Field>>& box,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int> SpectralPortrait
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<Field>>& box,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

template<typename Field>
Matrix<Int> TriangularSpectralPortrait
( const Matrix<Field>& U,
        Matrix<Base<Field>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<Field>>& box,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int> TriangularSpectralPortrait
( const AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<Field>>& box,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

template<typename Field>
Matrix<Int> TriangularSpectralPortrait
( const Matrix<Field>& U,
  const Matrix<Field>& Q,
        Matrix<Base<Field>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<Field>>& box,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int> TriangularSpectralPortrait
( const AbstractDistMatrix<Field>& U,
  const AbstractDistMatrix<Field>& Q,
        AbstractDistMatrix<Base<Field>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<Field>>& box,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

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
( const AbstractDistMatrix<Real>& U,
        AbstractDistMatrix<Real>& invNormMap,
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
( const AbstractDistMatrix<Real>& U,
  const AbstractDistMatrix<Real>& Q,
        AbstractDistMatrix<Real>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Real>& box,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Field>
Matrix<Int> HessenbergSpectralPortrait
( const Matrix<Field>& H,
        Matrix<Base<Field>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<Field>>& box,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int> HessenbergSpectralPortrait
( const AbstractDistMatrix<Field>& H,
        AbstractDistMatrix<Base<Field>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<Field>>& box,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

template<typename Field>
Matrix<Int> HessenbergSpectralPortrait
( const Matrix<Field>& H,
  const Matrix<Field>& Q,
        Matrix<Base<Field>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<Field>>& box,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int> HessenbergSpectralPortrait
( const AbstractDistMatrix<Field>& H,
  const AbstractDistMatrix<Field>& Q,
        AbstractDistMatrix<Base<Field>>& invNormMap,
        Int realSize,
        Int imagSize,
        SpectralBox<Base<Field>>& box,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

// (Pseudo-)Spectral window
// ------------------------
// Treat each pixel as being located a cell center and tesselate a box with
// said square cells
template<typename Field>
Matrix<Int> SpectralWindow
( const Matrix<Field>& A,
        Matrix<Base<Field>>& invNormMap,
        Complex<Base<Field>> center,
        Base<Field> realWidth,
        Base<Field> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int> SpectralWindow
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& invNormMap,
        Complex<Base<Field>> center,
        Base<Field> realWidth,
        Base<Field> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

template<typename Field>
Matrix<Int> TriangularSpectralWindow
( const Matrix<Field>& U,
        Matrix<Base<Field>>& invNormMap,
        Complex<Base<Field>> center,
        Base<Field> realWidth,
        Base<Field> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int> TriangularSpectralWindow
( const AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& invNormMap,
        Complex<Base<Field>> center,
        Base<Field> realWidth,
        Base<Field> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

template<typename Field>
Matrix<Int> TriangularSpectralWindow
( const Matrix<Field>& U,
  const Matrix<Field>& Q,
        Matrix<Base<Field>>& invNormMap,
        Complex<Base<Field>> center,
        Base<Field> realWidth,
        Base<Field> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int> TriangularSpectralWindow
( const AbstractDistMatrix<Field>& U,
  const AbstractDistMatrix<Field>& Q,
        AbstractDistMatrix<Base<Field>>& invNormMap,
        Complex<Base<Field>> center,
        Base<Field> realWidth,
        Base<Field> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

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
( const AbstractDistMatrix<Real>& U,
        AbstractDistMatrix<Real>& invNormMap,
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
( const AbstractDistMatrix<Real>& U,
  const AbstractDistMatrix<Real>& Q,
        AbstractDistMatrix<Real>& invNormMap,
        Complex<Real> center,
        Real realWidth,
        Real imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Field>
Matrix<Int> HessenbergSpectralWindow
( const Matrix<Field>& H,
        Matrix<Base<Field>>& invNormMap,
        Complex<Base<Field>> center,
        Base<Field> realWidth,
        Base<Field> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int> HessenbergSpectralWindow
( const AbstractDistMatrix<Field>& H,
        AbstractDistMatrix<Base<Field>>& invNormMap,
        Complex<Base<Field>> center,
        Base<Field> realWidth,
        Base<Field> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

template<typename Field>
Matrix<Int> HessenbergSpectralWindow
( const Matrix<Field>& H,
  const Matrix<Field>& Q,
        Matrix<Base<Field>>& invNormMap,
        Complex<Base<Field>> center,
        Base<Field> realWidth,
        Base<Field> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int> HessenbergSpectralWindow
( const AbstractDistMatrix<Field>& H,
  const AbstractDistMatrix<Field>& Q,
        AbstractDistMatrix<Base<Field>>& invNormMap,
        Complex<Base<Field>> center,
        Base<Field> realWidth,
        Base<Field> imagWidth,
        Int realSize,
        Int imagSize,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

// (Pseudo-)Spectral cloud
// -----------------------
template<typename Field>
Matrix<Int> SpectralCloud
( const Matrix<Field>& A,
  const Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int,VR,STAR> SpectralCloud
( const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Complex<Base<Field>>>& shifts,
        AbstractDistMatrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

template<typename Field>
Matrix<Int> TriangularSpectralCloud
( const Matrix<Field>& U,
  const Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int,VR,STAR> TriangularSpectralCloud
( const AbstractDistMatrix<Field>& U,
  const AbstractDistMatrix<Complex<Base<Field>>>& shifts,
        AbstractDistMatrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

template<typename Field>
Matrix<Int> TriangularSpectralCloud
( const Matrix<Field>& U,
  const Matrix<Field>& Q,
  const Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int,VR,STAR> TriangularSpectralCloud
( const AbstractDistMatrix<Field>& U,
  const AbstractDistMatrix<Field>& Q,
  const AbstractDistMatrix<Complex<Base<Field>>>& shifts,
        AbstractDistMatrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

template<typename Real>
Matrix<Int> QuasiTriangularSpectralCloud
( const Matrix<Real>& U,
  const Matrix<Complex<Real>>& shifts,
        Matrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> QuasiTriangularSpectralCloud
( const AbstractDistMatrix<Real>& U,
  const AbstractDistMatrix<Complex<Real>>& shifts,
        AbstractDistMatrix<Real>& invNorms,
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
( const AbstractDistMatrix<Real>& U,
  const AbstractDistMatrix<Real>& Q,
  const AbstractDistMatrix<Complex<Real>>& shifts,
        AbstractDistMatrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Field>
Matrix<Int> HessenbergSpectralCloud
( const Matrix<Field>& H,
  const Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int,VR,STAR> HessenbergSpectralCloud
( const AbstractDistMatrix<Field>& H,
  const AbstractDistMatrix<Complex<Base<Field>>>& shifts,
        AbstractDistMatrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

template<typename Field>
Matrix<Int> HessenbergSpectralCloud
( const Matrix<Field>& H,
  const Matrix<Field>& Q,
  const Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );
template<typename Field>
DistMatrix<Int,VR,STAR> HessenbergSpectralCloud
( const AbstractDistMatrix<Field>& H,
  const AbstractDistMatrix<Field>& Q,
  const AbstractDistMatrix<Complex<Base<Field>>>& shifts,
        AbstractDistMatrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl=PseudospecCtrl<Base<Field>>() );

} // namespace El

#include <El/lapack_like/spectral/Schur.hpp>
#include <El/lapack_like/spectral/HermitianEig.hpp>
#include <El/lapack_like/spectral/SVD.hpp>
#include <El/lapack_like/spectral/Lanczos.hpp>
#include <El/lapack_like/spectral/ProductLanczos.hpp>

#endif // ifndef EL_SPECTRAL_HPP
