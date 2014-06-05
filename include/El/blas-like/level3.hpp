/*
   Copyright (c) 2009-2014, Jack Poulson
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
enum GemmAlgorithm {
  GEMM_DEFAULT,
  GEMM_SUMMA_A,
  GEMM_SUMMA_B,
  GEMM_SUMMA_C,
  GEMM_SUMMA_DOT,
  GEMM_CANNON
};

template<typename T>
void Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

template<typename T>
void Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C );

template<typename T>
void Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C, GemmAlgorithm alg=GEMM_DEFAULT );

template<typename T>
void Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
                 DistMatrix<T>& C, GemmAlgorithm alg=GEMM_DEFAULT );

// Hemm
// ====
template<typename T>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

template<typename T>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C );

// Her2k
// =====
template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C );

template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C );

template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
                 DistMatrix<T>& C );

// Herk
// ====
template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C );

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, Matrix<T>& C );

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C );

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, DistMatrix<T>& C );

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
  F alpha, const DistMatrix<F>& A, const DistMatrix<F,VR,STAR>& shifts,
  DistMatrix<F>& B );

template<typename Real>
void MultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  Complex<Real> alpha,
  const DistMatrix<Real>& A,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
        DistMatrix<Real>& BReal, DistMatrix<Real>& BImag );

// MultiShiftTrsm
// ==============
template<typename F>
void MultiShiftTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, Matrix<F>& U, const Matrix<F>& shifts, Matrix<F>& X );

template<typename F>
void MultiShiftTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const DistMatrix<F>& U, const DistMatrix<F,VR,STAR>& shifts,
  DistMatrix<F>& X );

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
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& B,
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
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C,
  bool conjugate=false );

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
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C,
  bool conjugate=false );

template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
                 DistMatrix<T>& C,
  bool conjugate=false );

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
  T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C,
  bool conjugate=false );

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, DistMatrix<T>& C,
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
void Trdtrmm( UpperOrLower uplo, DistMatrix<F>& A, bool conjugate=false );

template<typename F>
void Trdtrmm
( UpperOrLower uplo,
  DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& dOff, bool conjugate=false );

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
  T alpha, const DistMatrix<T>& A, DistMatrix<T>& X );

// Trsm
// ====
template<typename F>
void Trsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& B,
  bool checkIfSingular=false );

// TODO: Greatly improve (and allow the user to modify) the mechanism for 
//       choosing between the different TRSM algorithms.
template<typename F>
void Trsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& B,
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
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& X,
  bool checkIfSingular=true );

// Trtrmm
// ======
template<typename T>
void Trtrmm( UpperOrLower uplo, Matrix<T>& A, bool conjugate=false );

template<typename T>
void Trtrmm( UpperOrLower uplo, DistMatrix<T>& A, bool conjugate=false );

// TwoSidedTrmm
// ============
template<typename T>
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<T>& A, const Matrix<T>& B );

template<typename T>
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  DistMatrix<T>& A, const DistMatrix<T>& B );

// TwoSidedTrsm
// ============
template<typename F>
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& B );

template<typename F>
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag,
  DistMatrix<F>& A, const DistMatrix<F>& B );

// Inline convenience functions
// ############################

// Gemm
// ====
template<typename T,Dist AColDist,Dist ARowDist,
                    Dist BColDist,Dist BRowDist,
                    Dist CColDist,Dist CRowDist>
inline void LocalGemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,AColDist,ARowDist>& A,
           const DistMatrix<T,BColDist,BRowDist>& B,
  T beta,        DistMatrix<T,CColDist,CRowDist>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalGemm");
        if( orientationOfA == NORMAL && orientationOfB == NORMAL )
        {
            if( AColDist != CColDist ||
                ARowDist != BColDist ||
                BRowDist != CRowDist )
                LogicError("C[X,Y] = A[X,Z] B[Z,Y]");
            if( A.ColAlign() != C.ColAlign() )
                LogicError("A's cols must align with C's rows");
            if( A.RowAlign() != B.ColAlign() )
                LogicError("A's rows must align with B's cols");
            if( B.RowAlign() != C.RowAlign() )
                LogicError("B's rows must align with C's rows");
            if( A.Height() != C.Height() ||
                A.Width() != B.Height() ||
                B.Width() != C.Width() )
                LogicError
                ("Nonconformal LocalGemmNN:\n",
                 "  A ~ ",A.Height()," x ",A.Width(),"\n",
                 "  B ~ ",B.Height()," x ",B.Width(),"\n",
                 "  C ~ ",C.Height()," x ",C.Width());
        }
        else if( orientationOfA == NORMAL )
        {
            if( AColDist != CColDist ||
                ARowDist != BRowDist ||
                BColDist != CRowDist )
                LogicError("C[X,Y] = A[X,Z] (B[Y,Z])^(T/H)");
            if( A.ColAlign() != C.ColAlign() )
                LogicError("A's cols must align with C's rows");
            if( A.RowAlign() != B.RowAlign() )
                LogicError("A's rows must align with B's rows");
            if( B.ColAlign() != C.RowAlign() )
                LogicError("B's cols must align with C's rows");
            if( A.Height() != C.Height() ||
                A.Width() != B.Width() ||
                B.Height() != C.Width() )
                LogicError
                ("Nonconformal LocalGemmNT:\n",
                 "  A ~ ",A.Height()," x ",A.Width(),"\n",
                 "  B ~ ",B.Height()," x ",B.Width(),"\n",
                 "  C ~ ",C.Height()," x ",C.Width());
        }
        else if( orientationOfB == NORMAL )
        {
            if( ARowDist != CColDist ||
                AColDist != BColDist ||
                BRowDist != CRowDist )
                LogicError("C[X,Y] = (A[Z,X])^(T/H) B[Z,Y]");
            if( A.RowAlign() != C.ColAlign() )
                LogicError("A's rows must align with C's cols");
            if( A.ColAlign() != B.ColAlign() )
                LogicError("A's cols must align with B's cols");
            if( B.RowAlign() != C.RowAlign() )
                LogicError("B's rows must align with C's rows");
            if( A.Width() != C.Height() ||
                A.Height() != B.Height() ||
                B.Width() != C.Width() )
                LogicError
                ("Nonconformal LocalGemmTN:\n",
                 "  A ~ ",A.Height()," x ",A.Width(),"\n",
                 "  B ~ ",B.Height()," x ",B.Width(),"\n",
                 "  C ~ ",C.Height()," x ",C.Width());
        }
        else
        {
            if( ARowDist != CColDist ||
                AColDist != BRowDist ||
                BColDist != CRowDist )
                LogicError("C[X,Y] = (A[Z,X])^(T/H) (B[Y,Z])^(T/H)");
            if( A.RowAlign() != C.ColAlign() )
                LogicError("A's rows must align with C's cols");
            if( A.ColAlign() != B.RowAlign() )
                LogicError("A's cols must align with B's rows");
            if( B.ColAlign() != C.RowAlign() )
                LogicError("B's cols must align with C's rows");
            if( A.Width() != C.Height() ||
                A.Height() != B.Width() ||
                B.Height() != C.Width() )
                LogicError
                ("Nonconformal LocalGemmTT:\n",
                 "  A ~ ",A.Height()," x ",A.Width(),"\n",
                 "  B ~ ",B.Height()," x ",B.Width(),"\n",
                 "  C ~ ",C.Height()," x ",C.Width());
        }
    )
    Gemm
    ( orientationOfA , orientationOfB,
      alpha, A.LockedMatrix(), B.LockedMatrix(), beta, C.Matrix() );
}

template<typename T,Dist AColDist,Dist ARowDist,
                    Dist BColDist,Dist BRowDist,
                    Dist CColDist,Dist CRowDist>
inline void LocalGemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,AColDist,ARowDist>& A,
           const DistMatrix<T,BColDist,BRowDist>& B,
                 DistMatrix<T,CColDist,CRowDist>& C )
{
    DEBUG_ONLY(CallStackEntry cse("LocalGemm"))
    const Int m = ( orientationOfA==NORMAL ? A.Height() : A.Width() );
    const Int n = ( orientationOfB==NORMAL ? B.Width() : B.Height() );
    Zeros( C, m, n );
    LocalGemm( orientationOfA, orientationOfB, alpha, A, B, T(0), C );
}

// MultiShiftQuasiTrsm
// ===================
template<typename F,Dist shiftColDist,
                    Dist     XColDist,Dist XRowDist>
inline void
LocalMultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const DistMatrix<F,STAR,STAR>& A,
           const DistMatrix<F,shiftColDist,STAR    >& shifts,
                 DistMatrix<F,    XColDist,XRowDist>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalMultiShiftQuasiTrsm");
        if( (side == LEFT &&  (     XColDist != STAR ||
                                shiftColDist != XRowDist) ) ||
            (side == RIGHT && (     XRowDist != STAR ||
                                shiftColDist != XColDist) ) )
            LogicError
            ("Dist of RHS and shifts must conform with that of triangle");
    )
    MultiShiftQuasiTrsm
    ( side, uplo, orientation,
      alpha, A.LockedMatrix(), shifts.LockedMatrix(), X.Matrix() );
}

template<typename Real,Dist shiftColDist,
                       Dist     XColDist,Dist XRowDist>
inline void
LocalMultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  Complex<Real> alpha,
  const DistMatrix<Real,         STAR,        STAR    >& A,
  const DistMatrix<Complex<Real>,shiftColDist,STAR    >& shifts,
        DistMatrix<Real,             XColDist,XRowDist>& XReal,
        DistMatrix<Real,             XColDist,XRowDist>& XImag )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalMultiShiftQuasiTrsm");
        if( (side == LEFT &&  (     XColDist != STAR ||
                                shiftColDist != XRowDist) ) ||
            (side == RIGHT && (     XRowDist != STAR ||
                                shiftColDist != XColDist) ) )
            LogicError
            ("Dist of RHS and shifts must conform with that of triangle");
    )
    MultiShiftQuasiTrsm
    ( side, uplo, orientation,
      alpha, A.LockedMatrix(), shifts.LockedMatrix(),
             XReal.Matrix(), XImag.Matrix() );
}

// QuasiTrsm
// =========
template<typename F,Dist XColDist,Dist XRowDist>
inline void
LocalQuasiTrsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation,
  F alpha, const DistMatrix<F,STAR,STAR>& A,
                 DistMatrix<F,XColDist,XRowDist>& X,
  bool checkIfSingular=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalQuasiTrsm");
        if( (side == LEFT && XColDist != STAR) ||
            (side == RIGHT && XRowDist != STAR) )
            LogicError
            ("Dist of RHS must conform with that of triangle");
    )
    QuasiTrsm
    ( side, uplo, orientation,
      alpha, A.LockedMatrix(), X.Matrix(), checkIfSingular );
}

// Trdtrmm
// =======
template<typename T>
inline void LocalTrdtrmm
( UpperOrLower uplo, DistMatrix<T,STAR,STAR>& A, bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("LocalTrdtrmm"))
    Trdtrmm( uplo, A.Matrix(), conjugate );
}

template<typename T>
inline void
LocalTrdtrmm
( UpperOrLower uplo,
  DistMatrix<T,STAR,STAR>& A, const DistMatrix<T,STAR,STAR>& dOff,
  bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("LocalTrdtrmm"))
    Trdtrmm( uplo, A.Matrix(), dOff.LockedMatrix(), conjugate );
}

// Trmm
// ====
template<typename T,Dist BColDist,Dist BRowDist>
inline void LocalTrmm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const DistMatrix<T,STAR,STAR>& A,
                 DistMatrix<T,BColDist,BRowDist>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrmm");
        if( (side == LEFT && BColDist != STAR) ||
            (side == RIGHT && BRowDist != STAR) )
            LogicError
            ("Dist of RHS must conform with that of triangle");
    )
    Trmm
    ( side, uplo, orientation, diag, alpha, A.LockedMatrix(), B.Matrix() );
}

// Trr2k
// =====
template<typename T>
void Trr2k
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
           const Matrix<T>& C, const Matrix<T>& D,
  T beta,        Matrix<T>& E );
template<typename T>
void Trr2k
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E );

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  );

// Trrk
// ====
template<typename T>
void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C );
template<typename T>
void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,  STAR>& A,
           const DistMatrix<T,STAR,MR  >& B,
  T beta,        DistMatrix<T,MC,  MR  >& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A,
           const DistMatrix<T,MR,STAR>& B,
  T beta,        DistMatrix<T>& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A,
           const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T,MC,  MR>& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A,
           const DistMatrix<T,MR,  STAR>& B,
  T beta,        DistMatrix<T,MC,  MR  >& C );

// Trsm
// ====
template<typename F,Dist XColDist,Dist XRowDist>
inline void
LocalTrsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,STAR>& A,
                 DistMatrix<F,XColDist,XRowDist>& X,
  bool checkIfSingular=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrsm");
        if( (side == LEFT && XColDist != STAR) ||
            (side == RIGHT && XRowDist != STAR) )
            LogicError
            ("Dist of RHS must conform with that of triangle");
    )
    // NOTE: Is this prototype available yet?!?
    Trsm
    ( side, uplo, orientation, diag,
      alpha, A.LockedMatrix(), X.Matrix(), checkIfSingular );
}

// TODO: Find a better name and/or home for this utility function
template<typename F,Dist colDist>
inline void AddInLocalData
( const DistMatrix<F,colDist,STAR>& X1, DistMatrix<F,STAR,STAR>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("AddInLocalData"))
    const Int width = X1.Width();
    const Int localHeight = X1.LocalHeight();
    const Int stride = X1.ColStride();
    const Int offset = X1.ColShift();
    for( Int j=0; j<width; ++j )
    {
        F* ZColBuffer = Z.Buffer(0,j);
        const F* X1ColBuffer = X1.LockedBuffer(0,j);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            ZColBuffer[offset+stride*iLoc] += X1ColBuffer[iLoc];
    }
}

// Trstrm
// ======
template<typename F>
inline void
LocalTrstrm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,STAR>& A,
                 DistMatrix<F,STAR,STAR>& X,
  bool checkIfSingular=true )
{
    DEBUG_ONLY(CallStackEntry cse("LocalTrstrm"))
    Trstrm
    ( side, uplo, orientation, diag,
      alpha, A.LockedMatrix(), X.Matrix(), checkIfSingular );
}

// Trtrmm
// ======
template<typename T>
inline void
LocalTrtrmm
( UpperOrLower uplo, DistMatrix<T,STAR,STAR>& A, bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("LocalTrtrmm"))
    Trtrmm( uplo, A.Matrix(), conjugate );
}

// TwoSidedTrmm
// ============
template<typename T>
inline void
LocalTwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag,
  DistMatrix<T,STAR,STAR>& A, const DistMatrix<T,STAR,STAR>& B )
{
    DEBUG_ONLY(CallStackEntry cse("LocalTwoSidedTrmm"))
    TwoSidedTrmm( uplo, diag, A.Matrix(), B.LockedMatrix() );
}

// TwoSidedTrsm
// ============
template<typename F>
inline void
LocalTwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag,
  DistMatrix<F,STAR,STAR>& A, const DistMatrix<F,STAR,STAR>& B )
{
    DEBUG_ONLY(CallStackEntry cse("LocalTwoSidedTrsm"))
    TwoSidedTrsm( uplo, diag, A.Matrix(), B.LockedMatrix() );
}

} // namespace El

#endif // ifndef EL_BLAS3_HPP
