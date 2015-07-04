/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS3_HPP
#define EL_BLAS3_HPP

namespace El {

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

template<typename Ring>
void Gemm
( Ring alpha, const OrientedMatrix<Ring>& A, 
              const OrientedMatrix<Ring>& B, 
  Ring beta,                Matrix<Ring>& C );
template<typename Ring>
void Gemm
( Ring alpha, const OrientedMatrix<Ring>& A, 
              const OrientedMatrix<Ring>& B, 
                            Matrix<Ring>& C );

template<typename Ring>
void Gemm
( Ring alpha, const OrientedAbstractDistMatrix<Ring>& A, 
              const OrientedAbstractDistMatrix<Ring>& B,
  Ring beta,                AbstractDistMatrix<Ring>& C, 
  GemmAlgorithm alg=GEMM_DEFAULT );
template<typename Ring>
void Gemm
( Ring alpha, const OrientedAbstractDistMatrix<Ring>& A, 
              const OrientedAbstractDistMatrix<Ring>& B,
                            AbstractDistMatrix<Ring>& C, 
  GemmAlgorithm alg=GEMM_DEFAULT );

template<typename Ring>
void LocalGemm
( Ring alpha, const OrientedAbstractDistMatrix<Ring>& A,
              const OrientedAbstractDistMatrix<Ring>& B,
  Ring beta,                AbstractDistMatrix<Ring>& C );
template<typename Ring>
void LocalGemm
( Ring alpha, const OrientedAbstractDistMatrix<Ring>& A,
              const OrientedAbstractDistMatrix<Ring>& B,
                            AbstractDistMatrix<Ring>& C );

// Deprecated
// ----------
template<typename Ring>
void Gemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const Matrix<Ring>& A, const Matrix<Ring>& B, 
  Ring beta,        Matrix<Ring>& C );
template<typename Ring>
void Gemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B, 
                    Matrix<Ring>& C );
template<typename Ring>
void Gemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,        AbstractDistMatrix<Ring>& C, 
  GemmAlgorithm alg=GEMM_DEFAULT );
template<typename Ring>
void Gemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
                    AbstractDistMatrix<Ring>& C, 
  GemmAlgorithm alg=GEMM_DEFAULT );
template<typename Ring>
void LocalGemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const AbstractDistMatrix<Ring>& A,
              const AbstractDistMatrix<Ring>& B,
  Ring beta,        AbstractDistMatrix<Ring>& C );
template<typename Ring>
void LocalGemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const AbstractDistMatrix<Ring>& A,
              const AbstractDistMatrix<Ring>& B,
                    AbstractDistMatrix<Ring>& C );

// Hemm
// ====
template<typename Ring>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B, 
  Ring beta,        Matrix<Ring>& C );

template<typename Ring>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,        AbstractDistMatrix<Ring>& C );

// Herk
// ====
template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const Matrix<Ring>& A, 
  Base<Ring> beta,        Matrix<Ring>& C );
template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const Matrix<Ring>& A, 
                          Matrix<Ring>& C );

template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const AbstractDistMatrix<Ring>& A, 
  Base<Ring> beta,        AbstractDistMatrix<Ring>& C );
template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const AbstractDistMatrix<Ring>& A, 
                          AbstractDistMatrix<Ring>& C );

template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const SparseMatrix<Ring>& A, 
  Base<Ring> beta,        SparseMatrix<Ring>& C );
template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const SparseMatrix<Ring>& A, 
                          SparseMatrix<Ring>& C );

template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const DistSparseMatrix<Ring>& A,
  Base<Ring> beta,        DistSparseMatrix<Ring>& C );
template<typename Ring> 
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const DistSparseMatrix<Ring>& A,
                          DistSparseMatrix<Ring>& C );

// Her2k
// =====
template<typename Ring>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha,      const Matrix<Ring>& A, 
                   const Matrix<Ring>& B, 
  Base<Ring> beta,       Matrix<Ring>& C );

template<typename Ring>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B, 
                    Matrix<Ring>& C );

template<typename Ring>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha,      const AbstractDistMatrix<Ring>& A, 
                   const AbstractDistMatrix<Ring>& B,
  Base<Ring> beta,       AbstractDistMatrix<Ring>& C );

template<typename Ring>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
                    AbstractDistMatrix<Ring>& C );

// Multiply
// ========
// NOTE: The following routine multiplies a sparse matrix by a set of vectors
//       and is obviously not a BLAS routine. However, it is a basic linear
//       algebra routine making use of Elemental's core data structures, and
//       so this is the natural placement
template<typename Ring>
void Multiply
( Orientation orientation,
  Ring alpha, const SparseMatrix<Ring>& A, 
              const Matrix<Ring>& X,
  Ring beta,        Matrix<Ring>& Y );

template<typename Ring>
void Multiply
( Orientation orientation,
  Ring alpha, const DistSparseMatrix<Ring>& A, 
              const DistMultiVec<Ring>& X,
  Ring beta,        DistMultiVec<Ring>& Y );

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
template<typename Ring>
void Symm
( LeftOrRight side, UpperOrLower uplo,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B, 
  Ring beta,        Matrix<Ring>& C,
  bool conjugate=false );
template<typename Ring>
void Symm
( LeftOrRight side, UpperOrLower uplo,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,        AbstractDistMatrix<Ring>& C, 
  bool conjugate=false );

namespace symm {

template<typename Ring>
void LocalAccumulateLL
( Orientation orientation, Ring alpha,
  const DistMatrix<Ring>& A,
  const DistMatrix<Ring,MC,  STAR>& B_MC_STAR,
  const DistMatrix<Ring,STAR,MR  >& BTrans_STAR_MR,
        DistMatrix<Ring,MC,  STAR>& Z_MC_STAR,
        DistMatrix<Ring,MR,  STAR>& Z_MR_STAR );

template<typename Ring>
void LocalAccumulateLU
( Orientation orientation, Ring alpha,
  const DistMatrix<Ring>& A,
  const DistMatrix<Ring,MC,  STAR>& B_MC_STAR,
  const DistMatrix<Ring,STAR,MR  >& BTrans_STAR_MR,
        DistMatrix<Ring,MC,  STAR>& Z_MC_STAR,
        DistMatrix<Ring,MR,  STAR>& Z_MR_STAR );

template<typename Ring>
void LocalAccumulateRL
( Orientation orientation, Ring alpha,
  const DistMatrix<Ring>& A,
  const DistMatrix<Ring,STAR,MC  >& B_STAR_MC,
  const DistMatrix<Ring,MR,  STAR>& BTrans_MR_STAR,
        DistMatrix<Ring,MC,  STAR>& ZTrans_MC_STAR,
        DistMatrix<Ring,MR,  STAR>& ZTrans_MR_STAR );

template<typename Ring>
void LocalAccumulateRU
( Orientation orientation, Ring alpha,
  const DistMatrix<Ring,MC,  MR  >& A,
  const DistMatrix<Ring,STAR,MC  >& B_STAR_MC,
  const DistMatrix<Ring,MR,  STAR>& BTrans_MR_STAR,
        DistMatrix<Ring,MC,  STAR>& ZTrans_MC_STAR,
        DistMatrix<Ring,MR,  STAR>& ZTrans_MR_STAR );

} // namespace symm

// Syrk
// ====
template<typename Ring>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const Matrix<Ring>& A, 
  Ring beta,        Matrix<Ring>& C,
  bool conjugate=false );
template<typename Ring>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const Matrix<Ring>& A, 
                    Matrix<Ring>& C,
  bool conjugate=false );

template<typename Ring>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
  Ring beta,        AbstractDistMatrix<Ring>& C, 
  bool conjugate=false );
template<typename Ring>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
                    AbstractDistMatrix<Ring>& C,
  bool conjugate=false );

template<typename Ring>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const SparseMatrix<Ring>& A, 
  Ring beta,        SparseMatrix<Ring>& C, 
  bool conjugate=false );
template<typename Ring>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const SparseMatrix<Ring>& A, 
                    SparseMatrix<Ring>& C, 
  bool conjugate=false );

template<typename Ring>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const DistSparseMatrix<Ring>& A, 
  Ring beta,        DistSparseMatrix<Ring>& C, 
  bool conjugate=false );
template<typename Ring>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const DistSparseMatrix<Ring>& A, 
                    DistSparseMatrix<Ring>& C, 
  bool conjugate=false );

// Syr2k
// =====
template<typename Ring>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B, 
  Ring beta,        Matrix<Ring>& C,
  bool conjugate=false );

template<typename Ring>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B, 
                    Matrix<Ring>& C,
  bool conjugate=false );

template<typename Ring>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,        AbstractDistMatrix<Ring>& C,
  bool conjugate=false );

template<typename Ring>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
                    AbstractDistMatrix<Ring>& C,
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
template<typename Ring>
void Trmm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  Ring alpha, const Matrix<Ring>& A, 
                    Matrix<Ring>& B );
template<typename Ring>
void Trmm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
                    AbstractDistMatrix<Ring>& X );

template<typename Ring>
void LocalTrmm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  Ring alpha, const DistMatrix<Ring,STAR,STAR>& A, 
                    AbstractDistMatrix<Ring>& B );

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
  F alpha, const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B,
  bool checkIfSingular=false, TrsmAlgorithm alg=TRSM_DEFAULT );

template<typename F>
void LocalTrsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,STAR>& A, AbstractDistMatrix<F>& X,
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
template<typename Ring>
void Trtrmm( UpperOrLower uplo, Matrix<Ring>& A, bool conjugate=false );
template<typename Ring>
void Trtrmm
( UpperOrLower uplo, AbstractDistMatrix<Ring>& A, bool conjugate=false );
template<typename Ring>
void Trtrmm
( UpperOrLower uplo, DistMatrix<Ring,STAR,STAR>& A, bool conjugate=false );

// TwoSidedTrmm
// ============
template<typename Ring>
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  Matrix<Ring>& A, const Matrix<Ring>& B );
template<typename Ring>
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  AbstractDistMatrix<Ring>& A, const AbstractDistMatrix<Ring>& B );
template<typename Ring>
void LocalTwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag,
  DistMatrix<Ring,STAR,STAR>& A, const DistMatrix<Ring,STAR,STAR>& B );

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
  DistMatrix<F,STAR,STAR>& A, const DistMatrix<F,STAR,STAR>& B );

// Trrk
// ====
template<typename Ring>
void Trrk
( UpperOrLower uplo, 
  Orientation orientA, Orientation orientB,
  Ring alpha, const Matrix<Ring>& A, const Matrix<Ring>& B,
  Ring beta,        Matrix<Ring>& C );
template<typename Ring>
void Trrk
( UpperOrLower uplo, 
  Orientation orientA, Orientation orientB,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,        AbstractDistMatrix<Ring>& C );
template<typename Ring>
void LocalTrrk
( UpperOrLower uplo,
  Ring alpha, const DistMatrix<Ring,MC,  STAR>& A,
              const DistMatrix<Ring,STAR,MR  >& B,
  Ring beta,        DistMatrix<Ring,MC,  MR  >& C );
template<typename Ring>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientB,
  Ring alpha, const DistMatrix<Ring,MC,STAR>& A,
              const DistMatrix<Ring,MR,STAR>& B,
  Ring beta,        DistMatrix<Ring>& C );
template<typename Ring>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientA,
  Ring alpha, const DistMatrix<Ring,STAR,MC>& A,
              const DistMatrix<Ring,STAR,MR>& B,
  Ring beta,        DistMatrix<Ring,MC,  MR>& C );
template<typename Ring>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientA, Orientation orientB,
  Ring alpha, const DistMatrix<Ring,STAR,MC  >& A,
              const DistMatrix<Ring,MR,  STAR>& B,
  Ring beta,        DistMatrix<Ring,MC,  MR  >& C );

// Trr2k
// =====
/*
template<typename Ring>
void Trr2k
( UpperOrLower uplo, 
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B,
  Ring beta,  const Matrix<Ring>& C, 
              const Matrix<Ring>& D,
  Ringgamma,        Matrix<Ring>& E );
*/
template<typename Ring>
void Trr2k
( UpperOrLower uplo,
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,  const AbstractDistMatrix<Ring>& C, 
              const AbstractDistMatrix<Ring>& D,
  Ring gamma,       AbstractDistMatrix<Ring>& E );

// The distributions of the oriented matrices must match
template<typename Ring>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,  const AbstractDistMatrix<Ring>& C, 
              const AbstractDistMatrix<Ring>& D,
  Ring gamma,       AbstractDistMatrix<Ring>& E );

} // namespace El

#endif // ifndef EL_BLAS3_HPP
