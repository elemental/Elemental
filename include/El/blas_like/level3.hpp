/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS3_HPP
#define EL_BLAS3_HPP

namespace El {

template<typename T> void SetLocalTrrkBlocksize( Int blocksize );
template<typename T> Int LocalTrrkBlocksize();

template<typename T> void SetLocalTrr2kBlocksize( Int blocksize );
template<typename T> Int LocalTrr2kBlocksize();

// Gemm
// ====
namespace GemmAlgorithmNS {
enum GemmAlgorithm {
  GEMM_DEFAULT,
  GEMM_SUMMA_A,
  GEMM_SUMMA_B,
  GEMM_SUMMA_C,
  GEMM_SUMMA_DOT,
  GEMM_CANNON
};
}
using namespace GemmAlgorithmNS;

template<typename T>
void Gemm
( Orientation orientA, Orientation orientB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

template<typename T>
void Gemm
( Orientation orientA, Orientation orientB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C );

template<typename T>
void Gemm
( Orientation orientA, Orientation orientB,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C, GemmAlgorithm alg=GEMM_DEFAULT );

template<typename T>
void Gemm
( Orientation orientA, Orientation orientB,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
                 AbstractDistMatrix<T>& C, GemmAlgorithm alg=GEMM_DEFAULT );

template<typename T>
void LocalGemm
( Orientation orientA, Orientation orientB,
  T alpha, const AbstractDistMatrix<T>& A,
           const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C );
template<typename T>
void LocalGemm
( Orientation orientA, Orientation orientB,
  T alpha, const AbstractDistMatrix<T>& A,
           const AbstractDistMatrix<T>& B,
                 AbstractDistMatrix<T>& C );

// Hemm
// ====
template<typename T>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

template<typename T>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C );

// Herk
// ====
template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const Matrix<T>& A, Base<T> beta, Matrix<T>& C );
template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const Matrix<T>& A, Matrix<T>& C );

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const AbstractDistMatrix<T>& A, 
  Base<T> beta,        AbstractDistMatrix<T>& C );
template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& C );

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const SparseMatrix<T>& A, 
  Base<T> beta,        SparseMatrix<T>& C );
template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const SparseMatrix<T>& A, 
                       SparseMatrix<T>& C );

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const DistSparseMatrix<T>& A,
  Base<T> beta,        DistSparseMatrix<T>& C );
template<typename T> 
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const DistSparseMatrix<T>& A,
                       DistSparseMatrix<T>& C );

// Her2k
// =====
template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha,       const Matrix<T>& A, const Matrix<T>& B, 
  Base<T> beta,        Matrix<T>& C );

template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C );

template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha,      const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  Base<T> beta,       AbstractDistMatrix<T>& C );

template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
                 AbstractDistMatrix<T>& C );

// Multiply
// ========
// NOTE: The following routine multiplies a sparse matrix by a set of vectors
//       and is obviously not a BLAS routine. However, it is a basic linear
//       algebra routine making use of Elemental's core data structures, and
//       so this is the natural placement
template<typename T>
void Multiply
( Orientation orientation,
  T alpha, const SparseMatrix<T>& A, const Matrix<T>& X,
  T beta,                                  Matrix<T>& Y );

template<typename T>
void Multiply
( Orientation orientation,
  T alpha,
  const DistSparseMatrix<T>& A,
  const DistMultiVec<T>& X,
  T beta,
        DistMultiVec<T>& Y );
template<typename T>
void Multiply
( Orientation orientation,
  T alpha,
  const DistSparseMatrix<T>& A,
  const AbstractDistMatrix<T>& X,
  T beta,
        AbstractDistMatrix<T>& Y );

// MultiShiftQuasiTrsm
// ===================
template<typename F>
void MultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const Matrix<F>& A, const Matrix<F>& shifts, Matrix<F>& B );
template<typename Real>
void MultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  Complex<Real> alpha,
  const Matrix<Real>& A,
  const Matrix<Complex<Real>>& shifts,
        Matrix<Real>& BReal, Matrix<Real>& BImag );
template<typename F>
void MultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& shifts,
  AbstractDistMatrix<F>& B );
template<typename Real>
void MultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  Complex<Real> alpha,
  const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Complex<Real>>& shifts,
        AbstractDistMatrix<Real>& BReal, AbstractDistMatrix<Real>& BImag );

template<typename F>
void LocalMultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const DistMatrix<F,STAR,STAR>& A,
           const AbstractDistMatrix<F>& shifts,
                 AbstractDistMatrix<F>& X );
template<typename Real>
void LocalMultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  Complex<Real> alpha,
  const DistMatrix<Real,STAR,STAR>& A,
  const AbstractDistMatrix<Complex<Real>>& shifts,
        AbstractDistMatrix<Real>& XReal,
        AbstractDistMatrix<Real>& XImag );

// MultiShiftTrsm
// ==============
template<typename F>
void MultiShiftTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, Matrix<F>& U, const Matrix<F>& shifts, Matrix<F>& X );
template<typename F>
void MultiShiftTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const AbstractDistMatrix<F>& U, const AbstractDistMatrix<F>& shifts,
  AbstractDistMatrix<F>& X );

// SafeMultiShiftTrsm
// ==================
template<typename F>
void SafeMultiShiftTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, Matrix<F>& A, const Matrix<F>& shifts,
  Matrix<F>& B, Matrix<F>& scales );
template<typename F>
void SafeMultiShiftTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& shifts,
  AbstractDistMatrix<F>& B, AbstractDistMatrix<F>& scales );
  
// QuasiTrsm
// =========
template<typename F>
void QuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const Matrix<F>& A, Matrix<F>& B,
  bool checkIfSingular=false );
template<typename F>
void QuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B,
  bool checkIfSingular=false );

template<typename F>
void LocalQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const DistMatrix<F,STAR,STAR>& A, AbstractDistMatrix<F>& X,
  bool checkIfSingular=false );

// Symm
// ====
template<typename T>
void Symm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C,
  bool conjugate=false );
template<typename T>
void Symm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C, bool conjugate=false );

namespace symm {

template<typename T>
void LocalAccumulateLL
( Orientation orientation, T alpha,
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,  STAR>& B_MC_STAR,
  const DistMatrix<T,STAR,MR  >& BTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR,
        DistMatrix<T,MR,  STAR>& Z_MR_STAR );

template<typename T>
void LocalAccumulateLU
( Orientation orientation, T alpha,
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,  STAR>& B_MC_STAR,
  const DistMatrix<T,STAR,MR  >& BTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR,
        DistMatrix<T,MR,  STAR>& Z_MR_STAR );

template<typename T>
void LocalAccumulateRL
( Orientation orientation, T alpha,
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC  >& B_STAR_MC,
  const DistMatrix<T,MR,  STAR>& BTrans_MR_STAR,
        DistMatrix<T,MC,  STAR>& ZTrans_MC_STAR,
        DistMatrix<T,MR,  STAR>& ZTrans_MR_STAR );

template<typename T>
void LocalAccumulateRU
( Orientation orientation, T alpha,
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,STAR,MC  >& B_STAR_MC,
  const DistMatrix<T,MR,  STAR>& BTrans_MR_STAR,
        DistMatrix<T,MC,  STAR>& ZTrans_MC_STAR,
        DistMatrix<T,MR,  STAR>& ZTrans_MR_STAR );

} // namespace symm

// Syrk
// ====
template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C,
  bool conjugate=false );
template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, Matrix<T>& C,
  bool conjugate=false );

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, 
  T beta,        AbstractDistMatrix<T>& C, bool conjugate=false );
template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& C,
  bool conjugate=false );

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const SparseMatrix<T>& A, 
  T beta,        SparseMatrix<T>& C, bool conjugate=false );
template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const SparseMatrix<T>& A, 
                 SparseMatrix<T>& C, bool conjugate=false );

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistSparseMatrix<T>& A, 
  T beta,        DistSparseMatrix<T>& C, bool conjugate=false );
template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistSparseMatrix<T>& A, 
                 DistSparseMatrix<T>& C, bool conjugate=false );

// Syr2k
// =====
template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C,
  bool conjugate=false );

template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C,
  bool conjugate=false );

template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C,
  bool conjugate=false );

template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
                 AbstractDistMatrix<T>& C,
  bool conjugate=false );

// Trdtrmm
// =======
template<typename F>
void Trdtrmm( UpperOrLower uplo, Matrix<F>& A, bool conjugate=false );
template<typename F>
void Trdtrmm
( UpperOrLower uplo, Matrix<F>& A, const Matrix<F>& dOff, 
  bool conjugate=false );

template<typename F>
void Trdtrmm
( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool conjugate=false );
template<typename F>
void Trdtrmm
( UpperOrLower uplo,
  AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& dOff, 
  bool conjugate=false );

template<typename F>
void Trdtrmm
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, bool conjugate=false );
template<typename F>
void Trdtrmm
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A, const DistMatrix<F,STAR,STAR>& dOff,
  bool conjugate=false );

// Trmm
// ====
template<typename T>
void Trmm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const Matrix<T>& A, Matrix<T>& B );
template<typename T>
void Trmm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& X );

template<typename T>
void LocalTrmm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const DistMatrix<T,STAR,STAR>& A, AbstractDistMatrix<T>& B );

// Trsm
// ====
namespace TrsmAlgorithmNS {
enum TrsmAlgorithm {
  TRSM_DEFAULT,
  TRSM_LARGE,
  TRSM_MEDIUM,
  TRSM_SMALL
};
}
using namespace TrsmAlgorithmNS;

template<typename F>
void Trsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& B,
  bool checkIfSingular=false );
template<typename F>
void Trsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha,
  const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& B,
  bool checkIfSingular=false, TrsmAlgorithm alg=TRSM_DEFAULT );

template<typename F>
void LocalTrsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha,
  const DistMatrix<F,STAR,STAR>& A,
        AbstractDistMatrix<F>& X,
  bool checkIfSingular=false );

// Trstrm
// ======
template<typename F>
void Trstrm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& X,
  bool checkIfSingular=true );
template<typename F>
void Trstrm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& X,
  bool checkIfSingular=true );
template<typename F>
void Trstrm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,STAR>& A, DistMatrix<F,STAR,STAR>& X,
  bool checkIfSingular=true );

// Trtrmm
// ======
template<typename T>
void Trtrmm( UpperOrLower uplo, Matrix<T>& A, bool conjugate=false );
template<typename T>
void Trtrmm
( UpperOrLower uplo, AbstractDistMatrix<T>& A, bool conjugate=false );
template<typename T>
void Trtrmm
( UpperOrLower uplo, DistMatrix<T,STAR,STAR>& A, bool conjugate=false );

// TwoSidedTrmm
// ============
template<typename T>
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  Matrix<T>& A, const Matrix<T>& B );
template<typename T>
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B );
template<typename T>
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag,
  DistMatrix<T,MC,MR,BLOCK>& A, const DistMatrix<T,MC,MR,BLOCK>& B );
template<typename T>
void LocalTwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag,
  DistMatrix<T,STAR,STAR>& A, const DistMatrix<T,STAR,STAR>& B );

// TwoSidedTrsm
// ============
template<typename F>
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  Matrix<F>& A, const Matrix<F>& B );
template<typename F>
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag,
  AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B );
template<typename F>
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag,
  DistMatrix<F,MC,MR,BLOCK>& A, const DistMatrix<F,MC,MR,BLOCK>& B );
template<typename F>
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag,
  DistMatrix<F,STAR,STAR>& A, const DistMatrix<F,STAR,STAR>& B );

// Trrk
// ====
template<typename T>
void Trrk
( UpperOrLower uplo, 
  Orientation orientA, Orientation orientB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C );
template<typename T>
void Trrk
( UpperOrLower uplo, 
  Orientation orientA, Orientation orientB,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,  STAR>& A,
           const DistMatrix<T,STAR,MR  >& B,
  T beta,        DistMatrix<T,MC,  MR  >& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientB,
  T alpha, const DistMatrix<T,MC,STAR>& A,
           const DistMatrix<T,MR,STAR>& B,
  T beta,        DistMatrix<T>& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientA,
  T alpha, const DistMatrix<T,STAR,MC>& A,
           const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T,MC,  MR>& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientA, Orientation orientB,
  T alpha, const DistMatrix<T,STAR,MC  >& A,
           const DistMatrix<T,MR,  STAR>& B,
  T beta,        DistMatrix<T,MC,  MR  >& C );

// Trr2k
// =====
/*
template<typename T>
void Trr2k
( UpperOrLower uplo, 
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,  const Matrix<T>& C, const Matrix<T>& D,
  Tgamma,        Matrix<T>& E );
*/
template<typename T>
void Trr2k
( UpperOrLower uplo,
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,  const AbstractDistMatrix<T>& C, const AbstractDistMatrix<T>& D,
  T gamma,       AbstractDistMatrix<T>& E );

// The distributions of the oriented matrices must match
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,  const AbstractDistMatrix<T>& C, const AbstractDistMatrix<T>& D,
  T gamma,       AbstractDistMatrix<T>& E );

// Hermitian from EVD
// ==================
// A := Z diag(w) Z^H, where w is real
template<typename F>
void HermitianFromEVD
( UpperOrLower uplo,
        Matrix<F>& A,
  const Matrix<Base<F>>& w,
  const Matrix<F>& Z );
template<typename F>
void HermitianFromEVD
( UpperOrLower uplo,
        AbstractDistMatrix<F>& A,
  const AbstractDistMatrix<Base<F>>& w,
  const AbstractDistMatrix<F>& Z );

// Normal from EVD
// ===============
// A := Z diag(w) Z^H, where w is complex
template<typename Real>
void NormalFromEVD
(       Matrix<Complex<Real>>& A,
  const Matrix<Complex<Real>>& w,
  const Matrix<Complex<Real>>& Z );
template<typename Real>
void NormalFromEVD
(       AbstractDistMatrix<Complex<Real>>& A,
  const AbstractDistMatrix<Complex<Real>>& w,
  const AbstractDistMatrix<Complex<Real>>& Z );

} // namespace El

#endif // ifndef EL_BLAS3_HPP
