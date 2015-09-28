/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPECTRAL_HPP
#define EL_SPECTRAL_HPP

namespace El {

// Hermitian eigenvalue solvers
// ============================
template<typename Real>
struct HermitianEigSubset
{
    bool indexSubset=false;
    Int lowerIndex=0, upperIndex=0;
 
    bool rangeSubset=false;
    Real lowerBound=0, upperBound=0;
};

template<typename Real>
struct HermitianSDCCtrl 
{
    Int cutoff=256;
    Int maxInnerIts=2, maxOuterIts=10;
    Real tol=0;
    Real spreadFactor=1e-6;
    bool progress=false;
};

template<typename F>
struct HermitianEigCtrl
{
    HermitianTridiagCtrl<F> tridiagCtrl;
    HermitianSDCCtrl<Base<F>> sdcCtrl;
    bool useSDC=false;
    bool timeStages=false;
};

// Compute eigenvalues
// -------------------
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,
  DistMatrix<Base<F>,STAR,STAR>& w, SortType sort=ASCENDING,
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, ElementalMatrix<F>& A, ElementalMatrix<Base<F>>& w,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

// Compute eigenpairs
// ------------------
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,
  DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, ElementalMatrix<F>& A, ElementalMatrix<Base<F>>& w, 
  ElementalMatrix<F>& Z, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

// Hermitian generalized definite eigenvalue solvers
// =================================================
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
void HermitianGenDefEig
( Pencil pencil, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HermitianGenDefEig
( Pencil pencil, UpperOrLower uplo,
  ElementalMatrix<F>& A, ElementalMatrix<F>& B,
  ElementalMatrix<Base<F>>& w, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
// Compute eigenpairs
// ------------------
template<typename F>
void HermitianGenDefEig
( Pencil pencil, UpperOrLower uplo,
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, Matrix<F>& X,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );
template<typename F>
void HermitianGenDefEig
( Pencil pencil, UpperOrLower uplo,
  ElementalMatrix<F>& A, ElementalMatrix<F>& B,
  ElementalMatrix<Base<F>>& w, ElementalMatrix<F>& X,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>> subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<F>& ctrl=HermitianEigCtrl<F>() );

// Hermitian tridiagonal eigenvalue solvers
// ========================================
// Compute eigenvalues
// --------------------
template<typename F>
void HermitianTridiagEig
( Matrix<Base<F>>& d, Matrix<F>& dSub, Matrix<Base<F>>& w, 
  SortType sort=ASCENDING, 
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>() );
template<typename F>
void HermitianTridiagEig
( const ElementalMatrix<Base<F>>& d, const ElementalMatrix<F>& dSub,
        ElementalMatrix<Base<F>>& w, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>() );
// Compute eigenpairs
// ------------------
template<typename F>
void HermitianTridiagEig
( Matrix<Base<F>>& d, Matrix<F>& dSub, Matrix<Base<F>>& w, Matrix<F>& Z,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>() );
template<typename F>
void HermitianTridiagEig
( const ElementalMatrix<Base<F>>& d, const ElementalMatrix<F>& dSub,
        ElementalMatrix<Base<F>>& w,       ElementalMatrix<F>& Z, 
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>() );

template<typename Real>
Int HermitianTridiagEigEstimate
( const ElementalMatrix<Real>& d, const ElementalMatrix<Real>& dSub,
        mpi::Comm wColComm, Real vl, Real vu );
// Z is assumed to be sufficiently large and properly aligned
template<typename Real>
void HermitianTridiagEigPostEstimate
( const ElementalMatrix<Real>& d, const ElementalMatrix<Real>& dSub,
        ElementalMatrix<Real>& w,       ElementalMatrix<Real>& Z, 
  SortType sort, Real vl, Real vu );

namespace herm_eig {

template<typename F>
void Sort( Matrix<Base<F>>& w, Matrix<F>& Z, SortType sort=ASCENDING );
template<typename Real,typename F>
void Sort
( ElementalMatrix<Real>& w, ElementalMatrix<F>& Z,
  SortType sort=ASCENDING );

} // namespace herm_eig

// Polar decomposition
// ===================
struct PolarCtrl 
{
    bool qdwh=false;
    bool colPiv=false;
    Int maxIts=20;
    mutable Int numIts=0;
};

template<typename F>
void Polar( Matrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void Polar( ElementalMatrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
void Polar
( Matrix<F>& A, Matrix<F>& P, const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void Polar
( ElementalMatrix<F>& A, ElementalMatrix<F>& P, 
  const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
void HermitianPolar
( UpperOrLower uplo, Matrix<F>& A, const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void HermitianPolar
( UpperOrLower uplo, ElementalMatrix<F>& A, 
  const PolarCtrl& ctrl=PolarCtrl() ); 
template<typename F>
void HermitianPolar
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& P, 
  const PolarCtrl& ctrl=PolarCtrl() );
template<typename F>
void HermitianPolar
( UpperOrLower uplo, ElementalMatrix<F>& A, ElementalMatrix<F>& P, 
  const PolarCtrl& ctrl=PolarCtrl() ); 

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
    Real spreadFactor=1e-6;
    bool random=true;
    bool progress=false;

    SignCtrl<Real> signCtrl;
};

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
};

template<typename F>
void Schur
( Matrix<F>& A, Matrix<Complex<Base<F>>>& w,
  bool fullTriangle=false, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( ElementalMatrix<F>& A, ElementalMatrix<Complex<Base<F>>>& w,
  bool fullTriangle=false, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( DistMatrix<F,MC,MR,BLOCK>& A,
  ElementalMatrix<Complex<Base<F>>>& w,
  bool fullTriangle=false,
  const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );

template<typename F>
void Schur
( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, Matrix<F>& Q,
  bool fullTriangle=true, const SchurCtrl<Base<F>> ctrl=SchurCtrl<Base<F>>() );
template<typename F>
void Schur
( ElementalMatrix<F>& A, ElementalMatrix<Complex<Base<F>>>& w, 
  ElementalMatrix<F>& Q, bool fullTriangle=true, 
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
( const Matrix<F>& dMain, const Matrix<F>& dSub, const Matrix<F>& dSup,
  Matrix<Complex<Base<F>>>& w );

template<typename F>
void QuasiTriangEig( const Matrix<F>& U, Matrix<Complex<Base<F>>>& w );
template<typename F>
void QuasiTriangEig
( const ElementalMatrix<F>& U, ElementalMatrix<Complex<Base<F>>>& w );

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
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl=
        HermitianEigCtrl<Complex<Base<F>>>() );
template<typename F>
void SkewHermitianEig
( UpperOrLower uplo, const ElementalMatrix<F>& G,
  ElementalMatrix<Base<F>>& wImag, SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl=
        HermitianEigCtrl<Complex<Base<F>>>() );

// Compute eigenpairs
// ------------------
template<typename F>
void SkewHermitianEig
( UpperOrLower uplo, const Matrix<F>& G,
  Matrix<Base<F>>& wImag, Matrix<Complex<Base<F>>>& Z,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl=
        HermitianEigCtrl<Complex<Base<F>>>() );
template<typename F>
void SkewHermitianEig
( UpperOrLower uplo, const ElementalMatrix<F>& G,
  ElementalMatrix<Base<F>>& wImag, ElementalMatrix<Complex<Base<F>>>& Z,
  SortType sort=ASCENDING,
  const HermitianEigSubset<Base<F>>& subset=HermitianEigSubset<Base<F>>(), 
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl=
        HermitianEigCtrl<Complex<Base<F>>>() );

// Singular Value Decomposition
// ============================
template<typename Real>
struct SVDCtrl 
{
    // Bidiagonal SVD options
    // ----------------------

    // Whether or not sequential implementations should use the QR algorithm
    // instead of (Cuppen's) Divide and Conquer when computing singular
    // vectors. When only singular values are requested, a bidiagonal DQDS
    // algorithm is always run.
    bool seqQR=false;

    // Chan's algorithm
    // ----------------

    // The minimum height/width ratio before preprocessing with a QR 
    // decomposition when only computing singular values
    double valChanRatio=1.2;

    // The minimum height/width ratio before preprocessing with a QR
    // decomposition when computing a full SVD
    double fullChanRatio=1.5;

    // Thresholding
    // ------------
    // NOTE: Currently only supported when computing both singular values
    //       and vectors

    // If the sufficiently small singular triplets should be thrown away.
    // When thresholded, a cross-product algorithm is used. This is often
    // advantageous since tridiagonal eigensolvers tend to have faster 
    // parallel implementations than bidiagonal SVD's.
    bool thresholded=false;

    // If the tolerance should be relative to the largest singular value
    bool relative=true;

    // The numerical tolerance for the thresholding. If this value is kept at
    // zero, then a value is automatically chosen based upon the matrix
    Real tol=0; 
};

// Compute the singular values
// ---------------------------
template<typename F>
void SVD( Matrix<F>& A, Matrix<Base<F>>& s );
template<typename F>
void SVD
( ElementalMatrix<F>& A,
  ElementalMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl=SVDCtrl<Base<F>>() );

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  Matrix<F>& A,
  Matrix<Base<F>>& s );
template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  ElementalMatrix<F>& A,
  ElementalMatrix<Base<F>>& s );

// Compute the full SVD
// --------------------
template<typename F>
void SVD
( Matrix<F>& A,
  Matrix<Base<F>>& s,
  Matrix<F>& V, 
  const SVDCtrl<Base<F>>& ctrl=SVDCtrl<Base<F>>() );
template<typename F>
void SVD
( ElementalMatrix<F>& A,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl=SVDCtrl<Base<F>>() );

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  Matrix<F>& A,
  Matrix<Base<F>>& s,
  Matrix<F>& U,
  Matrix<F>& V );
template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  ElementalMatrix<F>& A,
  ElementalMatrix<Base<F>>& s,
  ElementalMatrix<F>& U, 
  ElementalMatrix<F>& V );

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
        Matrix<Base<F>>& T,
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
        Matrix<Base<F>>& T,
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
        Matrix<Base<F>>& T,
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
        Matrix<Base<F>>& T,
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
    Int maxIts=200;
    Real tol=1e-6;
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
    mutable Real realWidth=0, imagWidth=0;
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
( const Matrix<F>& A, Matrix<Base<F>>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Base<F>>& box,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> SpectralPortrait
( const ElementalMatrix<F>& A, ElementalMatrix<Base<F>>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Base<F>>& box,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularSpectralPortrait
( const Matrix<F>& U, Matrix<Base<F>>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Base<F>>& box,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularSpectralPortrait
( const ElementalMatrix<F>& U, ElementalMatrix<Base<F>>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Base<F>>& box,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularSpectralPortrait
( const Matrix<F>& U, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Base<F>>& box,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularSpectralPortrait
( const ElementalMatrix<F>& U, const ElementalMatrix<F>& Q,
  ElementalMatrix<Base<F>>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Base<F>>& box,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename Real>
Matrix<Int> QuasiTriangularSpectralPortrait
( const Matrix<Real>& U, Matrix<Real>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralPortrait
( const ElementalMatrix<Real>& U, ElementalMatrix<Real>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> QuasiTriangularSpectralPortrait
( const Matrix<Real>& U, const Matrix<Real>& Q, Matrix<Real>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralPortrait
( const ElementalMatrix<Real>& U, const ElementalMatrix<Real>& Q,
  ElementalMatrix<Real>& invNormMap, Int realSize, Int imagSize, 
  SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename F>
Matrix<Int> HessenbergSpectralPortrait
( const Matrix<F>& H, Matrix<Base<F>>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Base<F>>& box,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergSpectralPortrait
( const ElementalMatrix<F>& H, ElementalMatrix<Base<F>>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Base<F>>& box,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> HessenbergSpectralPortrait
( const Matrix<F>& H, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap,
  Int realSize, Int imagSize, SpectralBox<Base<F>>& box,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergSpectralPortrait
( const ElementalMatrix<F>& H, const ElementalMatrix<F>& Q,
  ElementalMatrix<Base<F>>& invNormMap, Int realSize, Int imagSize, 
  SpectralBox<Base<F>>& box,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

// (Pseudo-)Spectral window
// ------------------------
// Treat each pixel as being located a cell center and tesselate a box with
// said square cells
template<typename F>
Matrix<Int> SpectralWindow
( const Matrix<F>& A, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> SpectralWindow
( const ElementalMatrix<F>& A, ElementalMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularSpectralWindow
( const Matrix<F>& U, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularSpectralWindow
( const ElementalMatrix<F>& U, ElementalMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> TriangularSpectralWindow
( const Matrix<F>& U, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> TriangularSpectralWindow
( const ElementalMatrix<F>& U, const ElementalMatrix<F>& Q,
  ElementalMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename Real>
Matrix<Int> QuasiTriangularSpectralWindow
( const Matrix<Real>& U,
  Matrix<Real>& invNormMap,
  Complex<Real> center, Real realWidth, Real imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralWindow
( const ElementalMatrix<Real>& U,
  ElementalMatrix<Real>& invNormMap,
  Complex<Real> center, Real realWidth, Real imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename Real>
Matrix<Int> QuasiTriangularSpectralWindow
( const Matrix<Real>& U, const Matrix<Real>& Q,
  Matrix<Real>& invNormMap,
  Complex<Real> center, Real realWidth, Real imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralWindow
( const ElementalMatrix<Real>& U, const ElementalMatrix<Real>& Q,
  ElementalMatrix<Real>& invNormMap,
  Complex<Real> center, Real realWidth, Real imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename F>
Matrix<Int> HessenbergSpectralWindow
( const Matrix<F>& H, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergSpectralWindow
( const ElementalMatrix<F>& H, ElementalMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

template<typename F>
Matrix<Int> HessenbergSpectralWindow
( const Matrix<F>& H, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int> HessenbergSpectralWindow
( const ElementalMatrix<F>& H, const ElementalMatrix<F>& Q,
  ElementalMatrix<Base<F>>& invNormMap,
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

// (Pseudo-)Spectral cloud
// -----------------------
template<typename F>
Matrix<Int> SpectralCloud
( const Matrix<F>& A, const Matrix<Complex<Base<F>>>& shifts,
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
( const Matrix<F>& U, const Matrix<Complex<Base<F>>>& shifts,
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
( const Matrix<F>& U, const Matrix<F>& Q,
  const Matrix<Complex<Base<F>>>& shifts, Matrix<Base<F>>& invNorms,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int,VR,STAR> TriangularSpectralCloud
( const ElementalMatrix<F>& U, const ElementalMatrix<F>& Q,
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
( const Matrix<Real>& U, const Matrix<Real>& Q,
  const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );
template<typename Real>
DistMatrix<Int,VR,STAR> QuasiTriangularSpectralCloud
( const ElementalMatrix<Real>& U, const ElementalMatrix<Real>& Q,
  const ElementalMatrix<Complex<Real>>& shifts,
        ElementalMatrix<Real>& invNorms,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() );

template<typename F>
Matrix<Int> HessenbergSpectralCloud
( const Matrix<F>& H, const Matrix<Complex<Base<F>>>& shifts,
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
( const Matrix<F>& H, const Matrix<F>& Q,
  const Matrix<Complex<Base<F>>>& shifts, Matrix<Base<F>>& invNorms,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );
template<typename F>
DistMatrix<Int,VR,STAR> HessenbergSpectralCloud
( const ElementalMatrix<F>& H, const ElementalMatrix<F>& Q,
  const ElementalMatrix<Complex<Base<F>>>& shifts,
        ElementalMatrix<Base<F>>& invNorms,
  PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

} // namespace El

#endif // ifndef EL_SPECTRAL_HPP
