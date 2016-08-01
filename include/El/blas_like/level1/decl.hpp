/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS1_DECL_HPP
#define EL_BLAS1_DECL_HPP

namespace El {

// TODO: More 'Contract' routines, e.g., {Contract,ContractedAxpy},
//       which sum results over the teams of processes that shared data in the
//       original distribution but do not in the final distribution. 
//       For example, a contraction of the form (U,Collect(V)) -> (U,V)
//       would perform the equivalent of an MPI_Reduce_scatter summation over 
//       the team of processes defining the 'V' row distribution.

// Adjoint
// =======
template<typename T>
void Adjoint( const Matrix<T>& A, Matrix<T>& B );
template<typename T>
void Adjoint( const ElementalMatrix<T>& A, ElementalMatrix<T>& B );
template<typename T>
void Adjoint( const BlockMatrix<T>& A, BlockMatrix<T>& B );
template<typename T>
void Adjoint( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );
template<typename T>
void Adjoint( const SparseMatrix<T>& A, SparseMatrix<T>& B );
template<typename T>
void Adjoint( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B );

// AdjointContract
// ===============
template<typename T>
void AdjointContract
( const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B );
template<typename T>
void AdjointContract
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B );

// AdjointAxpy
// ===========
template<typename T,typename S>
void AdjointAxpy
( S alpha, const Matrix<T>& X, Matrix<T>& Y );
template<typename T,typename S>
void AdjointAxpy
( S alpha, const SparseMatrix<T>& X, SparseMatrix<T>& Y );
template<typename T,typename S>
void AdjointAxpy
( S alpha, const ElementalMatrix<T>& X, ElementalMatrix<T>& Y );
template<typename T,typename S>
void AdjointAxpy
( S alpha, const BlockMatrix<T>& X, BlockMatrix<T>& Y );
template<typename T,typename S>
void AdjointAxpy
( S alpha, const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y );

// AdjointAxpyContract
// ===================
template<typename T>
void AdjointAxpyContract
( T alpha, const ElementalMatrix<T>& A,
                 ElementalMatrix<T>& B );
template<typename T>
void AdjointAxpyContract
( T alpha, const BlockMatrix<T>& A,
                 BlockMatrix<T>& B );

// AllReduce
// =========
template<typename T>
void AllReduce
( Matrix<T>& A, mpi::Comm comm, mpi::Op op=mpi::SUM );
template<typename T>
void AllReduce
( AbstractDistMatrix<T>& A, mpi::Comm comm, mpi::Op op=mpi::SUM );

// Axpy
// ====
template<typename T,typename S>
void Axpy( S alpha, const Matrix<T>& X, Matrix<T>& Y );
template<typename T,typename S>
void Axpy( S alpha, const ElementalMatrix<T>& X, ElementalMatrix<T>& Y );
template<typename T,typename S>
void Axpy( S alpha, const BlockMatrix<T>& X, BlockMatrix<T>& Y );
template<typename T,typename S>
void Axpy( S alpha, const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y );
template<typename T,typename S>
void Axpy( S alpha, const DistMultiVec<T>& X, DistMultiVec<T>& Y );
template<typename T,typename S>
void Axpy( S alpha, const SparseMatrix<T>& X, SparseMatrix<T>& Y );
template<typename T,typename S>
void Axpy( S alpha, const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y );

namespace axpy {
namespace util {

template<typename T>
void InterleaveMatrixUpdate
( T alpha, Int localHeight, Int localWidth,
  const T* A, Int colStrideA, Int rowStrideA,
        T* B, Int colStrideB, Int rowStrideB );
template<typename T>
void UpdateWithLocalData
( T alpha, const ElementalMatrix<T>& A, DistMatrix<T,STAR,STAR>& B );

} // namespace util
} // namespace axpy

// AxpyContract
// ============

template<typename T>
void AxpyContract
( T alpha, const ElementalMatrix<T>& A, ElementalMatrix<T>& B );
template<typename T>
void AxpyContract
( T alpha, const BlockMatrix<T>& A, BlockMatrix<T>& B );

// AxpyTrapezoid
// =============
template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alpha, const Matrix<T>& X, Matrix<T>& Y, Int offset=0 );
template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alpha,
  const ElementalMatrix<T>& X, ElementalMatrix<T>& Y, Int offset=0 );
template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alpha,
  const BlockMatrix<T>& X, BlockMatrix<T>& Y, Int offset=0 );
template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alpha, 
  const SparseMatrix<T>& X, SparseMatrix<T>& Y, Int offset=0 );
template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alpha, 
  const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y, Int offset=0 );

template<typename T>
void LocalAxpyTrapezoid
( UpperOrLower uplo, T alpha,
  const AbstractDistMatrix<T>& X,
        AbstractDistMatrix<T>& Y, Int offset=0 );

// Broadcast
// =========
template<typename T>
void Broadcast( Matrix<T>& A, mpi::Comm comm, int rank=0 );
template<typename T>
void Broadcast( AbstractDistMatrix<T>& A, mpi::Comm comm, int rank=0 );

// Column norms
// ============

// One norms
// ---------
// TODO

// Two norms
// ---------
template<typename F>
void ColumnTwoNorms
( const Matrix<F>& X, Matrix<Base<F>>& norms );
template<typename F,Dist U,Dist V>
void ColumnTwoNorms
( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,V,STAR>& norms );
template<typename F>
void ColumnTwoNorms
( const DistMultiVec<F>& X, Matrix<Base<F>>& norms );
template<typename F>
void ColumnTwoNorms
( const SparseMatrix<F>& X, Matrix<Base<F>>& norms );
template<typename F>
void ColumnTwoNorms
( const DistSparseMatrix<F>& X, DistMultiVec<Base<F>>& norms );

// Separated complex data
// ^^^^^^^^^^^^^^^^^^^^^^
template<typename Real,typename=EnableIf<IsReal<Real>>>
void ColumnTwoNorms
( const Matrix<Real>& XReal,
  const Matrix<Real>& XImag, 
        Matrix<Real>& norms );
template<typename Real,Dist U,Dist V,typename=EnableIf<IsReal<Real>>>
void ColumnTwoNorms
( const DistMatrix<Real,U,V>& XReal,
  const DistMatrix<Real,U,V>& XImag, 
        DistMatrix<Real,V,STAR>& norms );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void ColumnTwoNorms
( const DistMultiVec<Real>& XReal,
  const DistMultiVec<Real>& XImag, 
        Matrix<Real>& norms );

// Max norms
// ---------
template<typename F>
void ColumnMaxNorms
( const Matrix<F>& X, Matrix<Base<F>>& norms );
template<typename F,Dist U,Dist V>
void ColumnMaxNorms
( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,V,STAR>& norms );
template<typename F>
void ColumnMaxNorms
( const DistMultiVec<F>& X, Matrix<Base<F>>& norms );
template<typename F>
void ColumnMaxNorms
( const SparseMatrix<F>& X, Matrix<Base<F>>& norms );
template<typename F>
void ColumnMaxNorms
( const DistSparseMatrix<F>& X, DistMultiVec<Base<F>>& norms );

// Column minimum absolute values
// ==============================
// NOTE: While this is not a norm, it is often colloquially referred to as the
// "-infinity" norm
template<typename F>
void ColumnMinAbs
( const Matrix<F>& X, Matrix<Base<F>>& mins );
template<typename F,Dist U,Dist V>
void ColumnMinAbs
( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,V,STAR>& mins );
template<typename F>
void ColumnMinAbs
( const DistMultiVec<F>& X, Matrix<Base<F>>& mins );
template<typename F>
void ColumnMinAbs
( const SparseMatrix<F>& X, Matrix<Base<F>>& mins );
template<typename F>
void ColumnMinAbs
( const DistSparseMatrix<F>& X, DistMultiVec<Base<F>>& mins );

template<typename F>
void ColumnMinAbsNonzero
( const Matrix<F>& X, 
  const Matrix<Base<F>>& upperBounds,
        Matrix<Base<F>>& mins );
template<typename F,Dist U,Dist V>
void ColumnMinAbsNonzero
( const DistMatrix<F,U,V>& X, 
  const DistMatrix<Base<F>,V,STAR>& upperBounds,
        DistMatrix<Base<F>,V,STAR>& mins );
template<typename F>
void ColumnMinAbsNonzero
( const DistMultiVec<F>& X, 
  const Matrix<Base<F>>& upperBounds,
        Matrix<Base<F>>& mins );
template<typename F>
void ColumnMinAbsNonzero
( const SparseMatrix<F>& X, 
  const Matrix<Base<F>>& upperBounds,
        Matrix<Base<F>>& mins );
template<typename F>
void ColumnMinAbsNonzero
( const DistSparseMatrix<F>& X, 
  const DistMultiVec<Base<F>>& upperBounds,
        DistMultiVec<Base<F>>& mins );

// Row norms
// =========

// One norm
// --------
// TODO

// Two-norm
// --------
template<typename F>
void RowTwoNorms
( const Matrix<F>& X, Matrix<Base<F>>& norms );
template<typename F,Dist U,Dist V>
void RowTwoNorms
( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,U,STAR>& norms );
template<typename F>
void RowTwoNorms( const DistMultiVec<F>& X, DistMultiVec<Base<F>>& norms );
template<typename F>
void RowTwoNorms( const SparseMatrix<F>& X, Matrix<Base<F>>& norms );
template<typename F>
void RowTwoNorms( const DistSparseMatrix<F>& X, DistMultiVec<Base<F>>& norms );

// Max norm
// --------
template<typename F>
void RowMaxNorms
( const Matrix<F>& X, Matrix<Base<F>>& norms );
template<typename F,Dist U,Dist V>
void RowMaxNorms
( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,U,STAR>& norms );
template<typename F>
void RowMaxNorms( const DistMultiVec<F>& X, DistMultiVec<Base<F>>& norms );
template<typename F>
void RowMaxNorms( const SparseMatrix<F>& X, Matrix<Base<F>>& norms );
template<typename F>
void RowMaxNorms( const DistSparseMatrix<F>& X, DistMultiVec<Base<F>>& norms );

// Row minimum absolute values
// ===========================
// NOTE: While this is not a norm, it is often colloquially referred to as the
// "-infinity" norm
template<typename F>
void RowMinAbs
( const Matrix<F>& X, Matrix<Base<F>>& mins );
template<typename F,Dist U,Dist V>
void RowMinAbs
( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,U,STAR>& mins );
template<typename F>
void RowMinAbs
( const DistMultiVec<F>& X, DistMultiVec<Base<F>>& mins );
template<typename F>
void RowMinAbs
( const SparseMatrix<F>& X, Matrix<Base<F>>& mins );
template<typename F>
void RowMinAbs
( const DistSparseMatrix<F>& X, DistMultiVec<Base<F>>& mins );

template<typename F>
void RowMinAbsNonzero
( const Matrix<F>& X, 
  const Matrix<Base<F>>& upperBounds,
        Matrix<Base<F>>& mins );
template<typename F,Dist U,Dist V>
void RowMinAbsNonzero
( const DistMatrix<F,U,V>& X, 
  const DistMatrix<Base<F>,U,STAR>& upperBounds,
        DistMatrix<Base<F>,U,STAR>& mins );
template<typename F>
void RowMinAbsNonzero
( const DistMultiVec<F>& X, 
  const DistMultiVec<Base<F>>& upperBounds,
        DistMultiVec<Base<F>>& mins );
template<typename F>
void RowMinAbsNonzero
( const SparseMatrix<F>& X, 
  const Matrix<Base<F>>& upperBounds,
        Matrix<Base<F>>& mins );
template<typename F>
void RowMinAbsNonzero
( const DistSparseMatrix<F>& X, 
  const DistMultiVec<Base<F>>& upperBounds,
        DistMultiVec<Base<F>>& mins );

// Concatenation
// =============

// Horizontal concatenation: C := [A, B]
// -------------------------------------
template<typename T>
void HCat
( const Matrix<T>& A,
  const Matrix<T>& B, 
        Matrix<T>& C );
template<typename T>
void HCat
( const ElementalMatrix<T>& A,
  const ElementalMatrix<T>& B, 
        ElementalMatrix<T>& C );
template<typename T>
void HCat
( const SparseMatrix<T>& A,
  const SparseMatrix<T>& B, 
        SparseMatrix<T>& C );
template<typename T>
void HCat
( const DistSparseMatrix<T>& A,
  const DistSparseMatrix<T>& B, 
        DistSparseMatrix<T>& C );
template<typename T>
void HCat
( const DistMultiVec<T>& A,
  const DistMultiVec<T>& B, 
        DistMultiVec<T>& C );

// Vertical concatenation: C := [A; B]
// -----------------------------------
template<typename T>
void VCat
( const Matrix<T>& A,
  const Matrix<T>& B, 
        Matrix<T>& C );
template<typename T>
void VCat
( const ElementalMatrix<T>& A,
  const ElementalMatrix<T>& B, 
        ElementalMatrix<T>& C );
template<typename T>
void VCat
( const SparseMatrix<T>& A,
  const SparseMatrix<T>& B, 
        SparseMatrix<T>& C );
template<typename T>
void VCat
( const DistSparseMatrix<T>& A,
  const DistSparseMatrix<T>& B, 
        DistSparseMatrix<T>& C );
template<typename T>
void VCat
( const DistMultiVec<T>& A,
  const DistMultiVec<T>& B, 
        DistMultiVec<T>& C );

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
template<typename T>
void Conjugate( const ElementalMatrix<T>& A, ElementalMatrix<T>& B );

// ConjugateDiagonal
// =================
template<typename T>
void ConjugateDiagonal( Matrix<T>& A, Int offset=0 );
template<typename T>
void ConjugateDiagonal( AbstractDistMatrix<T>& A, Int offset=0 );

// ConjugateSubmatrix
// ==================
template<typename T>
void ConjugateSubmatrix
( Matrix<T>& A, const vector<Int>& I, const vector<Int>& J );
template<typename T>
void ConjugateSubmatrix
( AbstractDistMatrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J );

// Contract
// ========
template<typename T>
void Contract( const ElementalMatrix<T>& A, ElementalMatrix<T>& B );
template<typename T>
void Contract( const BlockMatrix<T>& A, BlockMatrix<T>& B );

// Copy
// ====

template<typename T>
void Copy( const Matrix<T>& A, Matrix<T>& B );
template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
void Copy( const Matrix<S>& A, Matrix<T>& B );

template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
void Copy( const ElementalMatrix<S>& A, ElementalMatrix<T>& B );

template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
void Copy( const BlockMatrix<S>& A, BlockMatrix<T>& B );

template<typename T>
void Copy( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );
template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
void Copy( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B );

template<typename T>
void CopyFromRoot
( const Matrix<T>& A, DistMatrix<T,CIRC,CIRC>& B,
  bool includingViewers=false );
template<typename T>
void CopyFromNonRoot( DistMatrix<T,CIRC,CIRC>& B,
  bool includingViewers=false );

template<typename T>
void CopyFromRoot
( const Matrix<T>& A, DistMatrix<T,CIRC,CIRC,BLOCK>& B,
  bool includingViewers=false );
template<typename T>
void CopyFromNonRoot( DistMatrix<T,CIRC,CIRC,BLOCK>& B,
  bool includingViewers=false );

void Copy( const Graph& A, Graph& B );
void Copy( const Graph& A, DistGraph& B );
void Copy( const DistGraph& A, Graph& B );
void Copy( const DistGraph& A, DistGraph& B );

void CopyFromRoot( const DistGraph& distGraph, Graph& graph );
void CopyFromNonRoot( const DistGraph& distGraph, int root=0 );

template<typename T>
void Copy( const SparseMatrix<T>& A, SparseMatrix<T>& B );

template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
void Copy( const SparseMatrix<S>& A, SparseMatrix<T>& B );

template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
void Copy( const SparseMatrix<S>& A, Matrix<T>& B );

template<typename T>
void Copy( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B );

template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
void Copy( const DistSparseMatrix<S>& A, DistSparseMatrix<T>& B );

template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
void Copy( const DistSparseMatrix<S>& A, AbstractDistMatrix<T>& B );

template<typename T>
void CopyFromRoot( const DistSparseMatrix<T>& ADist, SparseMatrix<T>& A );
template<typename T>
void CopyFromNonRoot( const DistSparseMatrix<T>& ADist, int root=0 );

template<typename T>
void Copy( const DistMultiVec<T>& A, DistMultiVec<T>& B );

template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
void Copy( const DistMultiVec<S>& A, DistMultiVec<T>& B );

template<typename T>
void Copy( const DistMultiVec<T>& A, AbstractDistMatrix<T>& B );
template<typename T>
void Copy( const AbstractDistMatrix<T>& A, DistMultiVec<T>& B );
template<typename T>
void CopyFromRoot( const DistMultiVec<T>& XDist, Matrix<T>& X );
template<typename T>
void CopyFromNonRoot( const DistMultiVec<T>& XDist, int root=0 );

namespace copy {
namespace util {

template<typename T>
void InterleaveMatrix
( Int height, Int width,
  const T* A, Int colStrideA, Int rowStrideA,
        T* B, Int colStrideB, Int rowStrideB );

template<typename T>
void ColStridedPack
( Int height, Int width,
  Int colAlign, Int colStride,
  const T* A,         Int ALDim,
        T* BPortions, Int portionSize );
template<typename T>
void ColStridedColumPack
( Int height, 
  Int colAlign, Int colStride,
  const T* A, 
        T* BPortions, Int portionSize );
template<typename T>
void ColStridedUnpack
( Int height, Int width,
  Int colAlign, Int colStride,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim );
template<typename T>
void PartialColStridedPack
( Int height, Int width,
  Int colAlign, Int colStride,
  Int colStrideUnion, Int colStridePart, Int colRankPart,
  Int colShiftA,
  const T* A,         Int ALDim,
        T* BPortions, Int portionSize );
template<typename T>
void PartialColStridedColumnPack
( Int height, 
  Int colAlign, Int colStride,
  Int colStrideUnion, Int colStridePart, Int colRankPart,
  Int colShiftA,
  const T* A,
        T* BPortions, Int portionSize );
template<typename T>
void PartialColStridedUnpack
( Int height, Int width,
  Int colAlign, Int colStride,
  Int colStrideUnion, Int colStridePart, Int colRankPart,
  Int colShiftB,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim );
template<typename T>
void PartialColStridedColumnUnpack
( Int height,
  Int colAlign, Int colStride,
  Int colStrideUnion, Int colStridePart, Int colRankPart,
  Int colShiftB,
  const T* APortions, Int portionSize,
        T* B );

template<typename T>
void RowStridedPack
( Int height, Int width,
  Int rowAlign, Int rowStride,
  const T* A,         Int ALDim,
        T* BPortions, Int portionSize );
template<typename T>
void RowStridedUnpack
( Int height, Int width,
  Int rowAlign, Int rowStride,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim );
template<typename T>
void PartialRowStridedPack
( Int height, Int width,
  Int rowAlign, Int rowStride,
  Int rowStrideUnion, Int rowStridePart, Int rowRankPart,
  Int rowShiftA,
  const T* A,         Int ALDim,
        T* BPortions, Int portionSize );
template<typename T>
void PartialRowStridedUnpack
( Int height, Int width,
  Int rowAlign, Int rowStride,
  Int rowStrideUnion, Int rowStridePart, Int rowRankPart,
  Int rowShiftB,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim );

template<typename T>
void StridedPack
( Int height, Int width,
  Int colAlign, Int colStride,
  Int rowAlign, Int rowStride,
  const T* A,         Int ALDim,
        T* BPortions, Int portionSize );
template<typename T>
void StridedUnpack
( Int height, Int width,
  Int colAlign, Int colStride,
  Int rowAlign, Int rowStride,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim );

} // namespace util
} // namespace copy

// DiagonalScale
// =============
template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<TDiag>& d, Matrix<T>& A );

template<typename TDiag,typename T,Dist U,Dist V,DistWrap wrapType=ELEMENT>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<TDiag>& d, DistMatrix<T,U,V,wrapType>& A );

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<TDiag>& d, AbstractDistMatrix<T>& A );

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<TDiag>& d, SparseMatrix<T>& A );

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<TDiag>& d, DistSparseMatrix<T>& A );

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<TDiag>& d, DistMultiVec<T>& X );

// DiagonalScaleTrapezoid
// ======================
template<typename TDiag,typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const Matrix<TDiag>& d, Matrix<T>& A, Int offset=0 );

template<typename TDiag,typename T,Dist U,Dist V>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const ElementalMatrix<TDiag>& d, DistMatrix<T,U,V>& A, Int offset=0 );

template<typename TDiag,typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const ElementalMatrix<TDiag>& d, ElementalMatrix<T>& A, Int offset=0 );

template<typename TDiag,typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const Matrix<TDiag>& d, SparseMatrix<T>& A, Int offset=0 );

template<typename TDiag,typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const DistMultiVec<TDiag>& d, DistSparseMatrix<T>& A, Int offset=0 );

// DiagonalSolve
// =============
template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<FDiag>& d, Matrix<F>& A, bool checkIfSingular=true );
template<typename F>
void SymmetricDiagonalSolve( const Matrix<Base<F>>& d, Matrix<F>& A );

template<typename FDiag,typename F,Dist U,Dist V>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const ElementalMatrix<FDiag>& d, DistMatrix<F,U,V>& A,
  bool checkIfSingular=true );

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const ElementalMatrix<FDiag>& d, ElementalMatrix<F>& A,
  bool checkIfSingular=true );

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<FDiag>& d, SparseMatrix<F>& A, 
  bool checkIfSingular=true );
template<typename F>
void SymmetricDiagonalSolve
( const Matrix<Base<F>>& d, SparseMatrix<F>& A );

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<FDiag>& d, DistSparseMatrix<F>& A, 
  bool checkIfSingular=true );
template<typename F>
void SymmetricDiagonalSolve
( const DistMultiVec<Base<F>>& d, DistSparseMatrix<F>& A );

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<FDiag>& d, DistMultiVec<F>& X, 
  bool checkIfSingular=true );

// Dot
// ===
template<typename T>
T Dot( const Matrix<T>& A, const Matrix<T>& B );
template<typename T>
T Dot( const ElementalMatrix<T>& A, const ElementalMatrix<T>& B );
template<typename T>
T Dot( const DistMultiVec<T>& A, const DistMultiVec<T>& B );

// Dotu
// ====
template<typename T>
T Dotu( const Matrix<T>& A, const Matrix<T>& B );
template<typename T>
T Dotu( const ElementalMatrix<T>& A, const ElementalMatrix<T>& B );
template<typename T>
T Dotu( const DistMultiVec<T>& A, const DistMultiVec<T>& B );

// EntrywiseFill
// =============
template<typename T>
void EntrywiseFill( Matrix<T>& A, function<T(void)> func );
template<typename T>
void EntrywiseFill( AbstractDistMatrix<T>& A, function<T(void)> func );
template<typename T>
void EntrywiseFill( DistMultiVec<T>& A, function<T(void)> func );

// EntrywiseMap
// ============
template<typename T>
void EntrywiseMap( Matrix<T>& A, function<T(T)> func );
template<typename T>
void EntrywiseMap( SparseMatrix<T>& A, function<T(T)> func );
template<typename T>
void EntrywiseMap( AbstractDistMatrix<T>& A, function<T(T)> func );
template<typename T>
void EntrywiseMap( DistSparseMatrix<T>& A, function<T(T)> func );
template<typename T>
void EntrywiseMap( DistMultiVec<T>& A, function<T(T)> func );

template<typename S,typename T>
void EntrywiseMap
( const Matrix<S>& A, Matrix<T>& B, function<T(S)> func );
template<typename S,typename T>
void EntrywiseMap
( const SparseMatrix<S>& A, SparseMatrix<T>& B, function<T(S)> func );
template<typename S,typename T>
void EntrywiseMap
( const ElementalMatrix<S>& A, ElementalMatrix<T>& B, 
  function<T(S)> func );
template<typename S,typename T>
void EntrywiseMap
( const BlockMatrix<S>& A, BlockMatrix<T>& B, 
  function<T(S)> func );
template<typename S,typename T>
void EntrywiseMap
( const DistSparseMatrix<S>& A, DistSparseMatrix<T>& B, 
  function<T(S)> func );
template<typename S,typename T>
void EntrywiseMap
( const DistMultiVec<S>& A, DistMultiVec<T>& B, 
  function<T(S)> func );

// Fill
// ====
template<typename T>
void Fill( Matrix<T>& A, T alpha );
template<typename T>
void Fill( AbstractDistMatrix<T>& A, T alpha );
template<typename T>
void Fill( DistMultiVec<T>& A, T alpha );
template<typename T>
void Fill( SparseMatrix<T>& A, T alpha );
template<typename T>
void Fill( DistSparseMatrix<T>& A, T alpha );

// FillDiagonal
// ============
template<typename T>
void FillDiagonal( Matrix<T>& A, T alpha, Int offset=0 );
template<typename T>
void FillDiagonal( AbstractDistMatrix<T>& A, T alpha, Int offset=0 );

// Full
// ====
template<typename T>
Matrix<T> Full( const SparseMatrix<T>& A );
// NOTE: A distributed version of the above does not exist because it is not
//       yet clear how Elemental should currently handle creating a grid within
//       a subroutine without causing a memory leak. DistMatrix may need to be
//       modified to allow for ownership of a grid.

// GetDiagonal
// ===========
template<typename T>
void GetDiagonal
( const Matrix<T>& A, Matrix<T>& d, Int offset=0 );
template<typename T>
void GetRealPartOfDiagonal
( const Matrix<T>& A, Matrix<Base<T>>& d, Int offset=0 );
template<typename T>
void GetImagPartOfDiagonal
( const Matrix<T>& A, Matrix<Base<T>>& d, Int offset=0 );

template<typename T>
Matrix<T> GetDiagonal( const Matrix<T>& A, Int offset=0 );
template<typename T>
Matrix<Base<T>> GetRealPartOfDiagonal( const Matrix<T>& A, Int offset=0 );
template<typename T>
Matrix<Base<T>> GetImagPartOfDiagonal( const Matrix<T>& A, Int offset=0 );

template<typename T,Dist U,Dist V>
void GetDiagonal
( const DistMatrix<T,U,V>& A, 
  ElementalMatrix<T>& d, Int offset=0 );
template<typename T,Dist U,Dist V>
void GetRealPartOfDiagonal
( const DistMatrix<T,U,V>& A, 
  ElementalMatrix<Base<T>>& d, Int offset=0 );
template<typename T,Dist U,Dist V>
void GetImagPartOfDiagonal
( const DistMatrix<T,U,V>& A, 
  ElementalMatrix<Base<T>>& d, Int offset=0 );
// Versions which will work for ElementalMatrix but which make use of a
// manual dynamic dispatch
template<typename T>
void GetDiagonal
( const ElementalMatrix<T>& A, 
  ElementalMatrix<T>& d, Int offset=0 );
template<typename T>
void GetRealPartOfDiagonal
( const ElementalMatrix<T>& A, 
  ElementalMatrix<Base<T>>& d, Int offset=0 );
template<typename T>
void GetImagPartOfDiagonal
( const ElementalMatrix<T>& A, 
  ElementalMatrix<Base<T>>& d, Int offset=0 );

template<typename T,Dist U,Dist V>
DistMatrix<T,DiagCol<U,V>(),DiagRow<U,V>()>
GetDiagonal( const DistMatrix<T,U,V>& A, Int offset=0 );
template<typename T,Dist U,Dist V>
DistMatrix<Base<T>,DiagCol<U,V>(),DiagRow<U,V>()>
GetRealPartOfDiagonal( const DistMatrix<T,U,V>& A, Int offset=0 );
template<typename T,Dist U,Dist V>
DistMatrix<Base<T>,DiagCol<U,V>(),DiagRow<U,V>()>
GetImagPartOfDiagonal( const DistMatrix<T,U,V>& A, Int offset=0 );

// GetMappedDiagonal
// =================
template<typename T,typename S>
void GetMappedDiagonal
( const Matrix<T>& A, Matrix<S>& d, function<S(T)> func, Int offset=0 );
template<typename T,typename S,Dist U,Dist V>
void GetMappedDiagonal
( const DistMatrix<T,U,V>& A, ElementalMatrix<S>& d, 
  function<S(T)> func, Int offset=0 );
template<typename T,typename S>
void GetMappedDiagonal
( const SparseMatrix<T>& A, Matrix<S>& d, 
  function<S(T)> func, Int offset=0 );
template<typename T,typename S>
void GetMappedDiagonal
( const DistSparseMatrix<T>& A, DistMultiVec<S>& d, 
  function<S(T)> func, Int offset=0 );

// GetSubgraph
// ===========
void GetSubgraph
( const Graph& graph,
        Range<Int> I,
        Range<Int> J,
        Graph& subgraph );
void GetSubgraph
( const Graph& graph,
  const vector<Int>& I,
        Range<Int> J,
        Graph& subgraph );
void GetSubgraph
( const Graph& graph,
        Range<Int> I,
  const vector<Int>& J,       
        Graph& subgraph );
void GetSubgraph
( const Graph& graph,
  const vector<Int>& I,
  const vector<Int>& J,       
        Graph& subgraph );

void GetSubgraph
( const DistGraph& graph,
        Range<Int> I,
        Range<Int> J,
        DistGraph& subgraph );
void GetSubgraph
( const DistGraph& graph,
  const vector<Int>& I,
        Range<Int> J,
        DistGraph& subgraph );
void GetSubgraph
( const DistGraph& graph,
        Range<Int> I,
  const vector<Int>& J,       
        DistGraph& subgraph );
void GetSubgraph
( const DistGraph& graph,
  const vector<Int>& I,
  const vector<Int>& J,       
        DistGraph& subgraph );

// GetSubmatrix
// ============

// Return a view
// ----------
template<typename T>
void GetSubmatrix
( const Matrix<T>& A,
        Range<Int> I,
        Range<Int> J,
        Matrix<T>& ASub );

template<typename T>
void GetSubmatrix
( const ElementalMatrix<T>& A,
        Range<Int> I,
        Range<Int> J, 
        ElementalMatrix<T>& ASub );


// Noncontiguous
// -------------
template<typename T>
void GetSubmatrix
( const Matrix<T>& A, 
        Range<Int> I,
  const vector<Int>& J, 
        Matrix<T>& ASub );
template<typename T>
void GetSubmatrix
( const Matrix<T>& A, 
  const vector<Int>& I,
  const Range<Int> J, 
        Matrix<T>& ASub );
template<typename T>
void GetSubmatrix
( const Matrix<T>& A, 
  const vector<Int>& I,
  const vector<Int>& J, 
        Matrix<T>& ASub );

template<typename T>
void GetSubmatrix
( const AbstractDistMatrix<T>& A, 
        Range<Int> I,
  const vector<Int>& J, 
        AbstractDistMatrix<T>& ASub );
template<typename T>
void GetSubmatrix
( const AbstractDistMatrix<T>& A, 
  const vector<Int>& I,
        Range<Int> J, 
        AbstractDistMatrix<T>& ASub );
template<typename T>
void GetSubmatrix
( const AbstractDistMatrix<T>& A, 
  const vector<Int>& I,
  const vector<Int>& J, 
        AbstractDistMatrix<T>& ASub );

template<typename T>
void GetSubmatrix
( const SparseMatrix<T>& A,
        Range<Int> I,
        Range<Int> J, 
        SparseMatrix<T>& ASub );
template<typename T>
void GetSubmatrix
( const SparseMatrix<T>& A,
        Range<Int> I,
  const vector<Int>& J, 
        SparseMatrix<T>& ASub );
template<typename T>
void GetSubmatrix
( const SparseMatrix<T>& A,
  const vector<Int>& I,
        Range<Int> J, 
        SparseMatrix<T>& ASub );
template<typename T>
void GetSubmatrix
( const SparseMatrix<T>& A,
  const vector<Int>& I,
  const vector<Int>& J, 
        SparseMatrix<T>& ASub );

template<typename T>
void GetSubmatrix
( const DistSparseMatrix<T>& A,
        Range<Int> I,
        Range<Int> J,
        DistSparseMatrix<T>& ASub );
template<typename T>
void GetSubmatrix
( const DistSparseMatrix<T>& A,
        Range<Int> I,
  const vector<Int>& J, 
        DistSparseMatrix<T>& ASub );
template<typename T>
void GetSubmatrix
( const DistSparseMatrix<T>& A,
  const vector<Int>& I,
        Range<Int> J,
        DistSparseMatrix<T>& ASub );
template<typename T>
void GetSubmatrix
( const DistSparseMatrix<T>& A,
  const vector<Int>& I,
  const vector<Int>& J, 
        DistSparseMatrix<T>& ASub );

template<typename T>
void GetSubmatrix
( const DistMultiVec<T>& A,
        Range<Int> I,
        Range<Int> J,
        DistMultiVec<T>& ASub );
template<typename T>
void GetSubmatrix
( const DistMultiVec<T>& A,
        Range<Int> I,
  const vector<Int>& J, 
        DistMultiVec<T>& ASub );
template<typename T>
void GetSubmatrix
( const DistMultiVec<T>& A,
  const vector<Int>& I,
        Range<Int> J,
        DistMultiVec<T>& ASub );
template<typename T>
void GetSubmatrix
( const DistMultiVec<T>& A,
  const vector<Int>& I,
  const vector<Int>& J, 
        DistMultiVec<T>& ASub );

// Hadamard
// ========
template<typename T>
void Hadamard( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C );
template<typename T>
void Hadamard
( const ElementalMatrix<T>& A,
  const ElementalMatrix<T>& B,
        ElementalMatrix<T>& C );
template<typename T>
void Hadamard
( const DistMultiVec<T>& A,
  const DistMultiVec<T>& B,
        DistMultiVec<T>& C );

// HilbertSchmidt
// ==============
template<typename T>
T HilbertSchmidt( const Matrix<T>& A, const Matrix<T>& B );
template<typename T>
T HilbertSchmidt
( const ElementalMatrix<T>& A, const ElementalMatrix<T>& C );
template<typename T>
T HilbertSchmidt( const DistMultiVec<T>& A, const DistMultiVec<T>& B );

// Imaginary part
// ==============
template<typename T>
void ImagPart
( const Matrix<T>& A, Matrix<Base<T>>& AImag );
template<typename T>
void ImagPart
( const ElementalMatrix<T>& A, ElementalMatrix<Base<T>>& AImag );
/* TODO: Sparse versions */

// IndexDependentFill
// ==================
template<typename T>
void IndexDependentFill( Matrix<T>& A, function<T(Int,Int)> func );
template<typename T>
void IndexDependentFill
( AbstractDistMatrix<T>& A, function<T(Int,Int)> func );

// IndexDependentMap
// =================
template<typename T>
void IndexDependentMap( Matrix<T>& A, function<T(Int,Int,T)> func );
template<typename T>
void IndexDependentMap
( AbstractDistMatrix<T>& A, function<T(Int,Int,T)> func );

template<typename S,typename T>
void IndexDependentMap
( const Matrix<S>& A, Matrix<T>& B, function<T(Int,Int,S)> func );
template<typename S,typename T>
void IndexDependentMap
( const ElementalMatrix<S>& A, ElementalMatrix<T>& B,
  function<T(Int,Int,S)> func );
template<typename S,typename T>
void IndexDependentMap
( const BlockMatrix<S>& A, BlockMatrix<T>& B,
  function<T(Int,Int,S)> func );

// Kronecker product
// =================
template<typename T>
void Kronecker( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C );
template<typename T>
void Kronecker
( const Matrix<T>& A, const Matrix<T>& B, ElementalMatrix<T>& C );
template<typename T>
void Kronecker
( const SparseMatrix<T>& A, const SparseMatrix<T>& B, SparseMatrix<T>& C );
template<typename T>
void Kronecker
( const SparseMatrix<T>& A, const Matrix<T>& B, SparseMatrix<T>& C );
template<typename T>
void Kronecker
( const Matrix<T>& A, const SparseMatrix<T>& B, SparseMatrix<T>& C );
template<typename T>
void Kronecker
( const SparseMatrix<T>& A, const SparseMatrix<T>& B, DistSparseMatrix<T>& C );
template<typename T>
void Kronecker
( const SparseMatrix<T>& A, const Matrix<T>& B, DistSparseMatrix<T>& C );
template<typename T>
void Kronecker
( const Matrix<T>& A, const SparseMatrix<T>& B, DistSparseMatrix<T>& C );

// MakeHermitian
// =============
template<typename T>
void MakeHermitian( UpperOrLower uplo, Matrix<T>& A );
template<typename T>
void MakeHermitian( UpperOrLower uplo, ElementalMatrix<T>& A );

template<typename T>
void MakeHermitian( UpperOrLower uplo, SparseMatrix<T>& A );
template<typename T>
void MakeHermitian( UpperOrLower uplo, DistSparseMatrix<T>& A );

// MakeDiagonalReal
// ================
template<typename T>
void MakeDiagonalReal( Matrix<T>& A, Int offset=0 );
template<typename T>
void MakeDiagonalReal( AbstractDistMatrix<T>& A, Int offset=0 );

// MakeReal
// ========
template<typename Real>
void MakeReal( Matrix<Real>& A );
template<typename Real>
void MakeReal( Matrix<Complex<Real>>& A );
template<typename T>
void MakeReal( AbstractDistMatrix<T>& A );

// MakeSubmatrixReal
// ================
template<typename T>
void MakeSubmatrixReal
( Matrix<T>& A, const vector<Int>& I, const vector<Int>& J );
template<typename T>
void MakeSubmatrixReal
( AbstractDistMatrix<T>& A, const vector<Int>& I, const vector<Int>& J );

// MakeSymmetric
// =============
template<typename T>
void MakeSymmetric( UpperOrLower uplo, Matrix<T>& A, bool conjugate=false );
template<typename T>
void MakeSymmetric
( UpperOrLower uplo, ElementalMatrix<T>& A, bool conjugate=false );

template<typename T>
void MakeSymmetric
( UpperOrLower uplo, SparseMatrix<T>& A, bool conjugate=false );
template<typename T>
void MakeSymmetric
( UpperOrLower uplo, DistSparseMatrix<T>& A, bool conjugate=false );

// MakeTrapezoidal
// ===============
template<typename T>
void MakeTrapezoidal( UpperOrLower uplo, Matrix<T>& A, Int offset=0 );
template<typename T>
void MakeTrapezoidal
( UpperOrLower uplo, AbstractDistMatrix<T>& A, Int offset=0 );

template<typename T>
void MakeTrapezoidal( UpperOrLower uplo, SparseMatrix<T>& A, Int offset=0 );
template<typename T>
void MakeTrapezoidal( UpperOrLower uplo, DistSparseMatrix<T>& A, Int offset=0 );

// Max
// ===
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Max( const Matrix<Real>& A );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Max( const AbstractDistMatrix<Real>& A );

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real SymmetricMax( UpperOrLower uplo, const Matrix<Real>& A );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real SymmetricMax( UpperOrLower uplo, const AbstractDistMatrix<Real>& A );

// MaxLoc
// ======
template<typename Real,typename=EnableIf<IsReal<Real>>>
Entry<Real> MaxLoc( const Matrix<Real>& A );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Entry<Real> MaxLoc( const AbstractDistMatrix<Real>& A );

template<typename Real,typename=EnableIf<IsReal<Real>>>
Entry<Real> SymmetricMaxLoc( UpperOrLower uplo, const Matrix<Real>& A );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Entry<Real>
SymmetricMaxLoc( UpperOrLower uplo, const AbstractDistMatrix<Real>& A );

template<typename Real,typename=EnableIf<IsReal<Real>>>
ValueInt<Real> VectorMaxLoc( const Matrix<Real>& x );
template<typename Real,typename=EnableIf<IsReal<Real>>>
ValueInt<Real> VectorMaxLoc( const AbstractDistMatrix<Real>& x );
template<typename Real,typename=EnableIf<IsReal<Real>>>
ValueInt<Real> VectorMaxLoc( const DistMultiVec<Real>& x );

// MaxAbs
// ======
template<typename T>
Base<T> MaxAbs( const Matrix<T>& A );
template<typename T>
Base<T> MaxAbs( const AbstractDistMatrix<T>& A );
template<typename T>
Base<T> MaxAbs( const SparseMatrix<T>& A );
template<typename T>
Base<T> MaxAbs( const DistSparseMatrix<T>& A );

template<typename T>
Base<T> SymmetricMaxAbs( UpperOrLower uplo, const Matrix<T>& A );
template<typename T>
Base<T> SymmetricMaxAbs( UpperOrLower uplo, const AbstractDistMatrix<T>& A );
template<typename T>
Base<T> SymmetricMaxAbs( UpperOrLower uplo, const SparseMatrix<T>& A );
template<typename T>
Base<T> SymmetricMaxAbs( UpperOrLower uplo, const DistSparseMatrix<T>& A );

// MaxAbsLoc
// =========
template<typename T>
ValueInt<Base<T>> VectorMaxAbsLoc( const Matrix<T>& x );
template<typename T>
ValueInt<Base<T>> VectorMaxAbsLoc( const AbstractDistMatrix<T>& x );

template<typename T>
Entry<Base<T>> MaxAbsLoc( const Matrix<T>& A );
template<typename T>
Entry<Base<T>> MaxAbsLoc( const AbstractDistMatrix<T>& A );
template<typename T>
Entry<Base<T>> MaxAbsLoc( const SparseMatrix<T>& A );
template<typename T>
Entry<Base<T>> MaxAbsLoc( const DistSparseMatrix<T>& A );

template<typename T>
Entry<Base<T>> 
SymmetricMaxAbsLoc( UpperOrLower uplo, const Matrix<T>& A );
template<typename T>
Entry<Base<T>> 
SymmetricMaxAbsLoc( UpperOrLower uplo, const AbstractDistMatrix<T>& A );
template<typename T>
Entry<Base<T>> 
SymmetricMaxAbsLoc( UpperOrLower uplo, const SparseMatrix<T>& A );
template<typename T>
Entry<Base<T>> 
SymmetricMaxAbsLoc( UpperOrLower uplo, const DistSparseMatrix<T>& A );

// Min
// ===
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Min( const Matrix<Real>& A );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Min( const AbstractDistMatrix<Real>& A );

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real SymmetricMin( UpperOrLower uplo, const Matrix<Real>& A );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real SymmetricMin( UpperOrLower uplo, const AbstractDistMatrix<Real>& A );

// MinLoc
// ======
template<typename Real,typename=EnableIf<IsReal<Real>>>
ValueInt<Real> VectorMinLoc( const Matrix<Real>& x );
template<typename Real,typename=EnableIf<IsReal<Real>>>
ValueInt<Real> VectorMinLoc( const AbstractDistMatrix<Real>& x );
template<typename Real,typename=EnableIf<IsReal<Real>>>
ValueInt<Real> VectorMinLoc( const DistMultiVec<Real>& x );

template<typename Real,typename=EnableIf<IsReal<Real>>>
Entry<Real> MinLoc( const Matrix<Real>& A );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Entry<Real> MinLoc( const AbstractDistMatrix<Real>& A );

template<typename Real,typename=EnableIf<IsReal<Real>>>
Entry<Real> SymmetricMinLoc( UpperOrLower uplo, const Matrix<Real>& A );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Entry<Real>
SymmetricMinLoc( UpperOrLower uplo, const AbstractDistMatrix<Real>& A );

// MinAbs
// ======
template<typename F>
Base<F> MinAbs( const Matrix<F>& A );
template<typename F>
Base<F> MinAbs( const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> SymmetricMinAbs( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricMinAbs( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// MinAbsLoc
// =========
template<typename F>
ValueInt<Base<F>> VectorMinAbsLoc( const Matrix<F>& x );
template<typename F>
ValueInt<Base<F>> VectorMinAbsLoc( const AbstractDistMatrix<F>& x );

template<typename F>
Entry<Base<F>> MinAbsLoc( const Matrix<F>& A );
template<typename F>
Entry<Base<F>> MinAbsLoc( const AbstractDistMatrix<F>& A );

template<typename F>
Entry<Base<F>> SymmetricMinAbsLoc( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Entry<Base<F>>
SymmetricMinAbsLoc( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Nrm2
// ====
template<typename F>
Base<F> Nrm2( const Matrix<F>& x );
template<typename F>
Base<F> Nrm2( const AbstractDistMatrix<F>& x );
template<typename F>
Base<F> Nrm2( const DistMultiVec<F>& x );

// QuasiDiagonalScale
// ==================
template<typename F,typename FMain>
void QuasiDiagonalScale
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<FMain>& d, const Matrix<F>& dSub,
  Matrix<F>& X, bool conjugated=false );
// TODO: Switch to full ElementalMatrix interface
template<typename F,typename FMain,Dist U,Dist V>
void QuasiDiagonalScale
( LeftOrRight side, UpperOrLower uplo,
  const ElementalMatrix<FMain>& d, const ElementalMatrix<F>& dSub,
  DistMatrix<F,U,V>& X, bool conjugated=false );

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
// TODO: Switch to full ElementalMatrix interface
template<typename F,typename FMain,Dist U,Dist V>
void
QuasiDiagonalSolve
( LeftOrRight side, UpperOrLower uplo,
  const ElementalMatrix<FMain>& d, const ElementalMatrix<F>& dSub,
  DistMatrix<F,U,V>& X, bool conjugated=false );

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

// Real part
// =========
template<typename T>
void RealPart( const Matrix<T>& A, Matrix<Base<T>>& AReal );
template<typename T>
void RealPart( const ElementalMatrix<T>& A, ElementalMatrix<Base<T>>& AReal );
/* TODO: Sparse versions */

// Reshape
// =======
template<typename T>
void Reshape( Int m, Int n, const Matrix<T>& A, Matrix<T>& B );
template<typename T>
Matrix<T> Reshape( Int m, Int n, const Matrix<T>& A );

template<typename T>
void Reshape
( Int m, Int n, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );
template<typename T>
DistMatrix<T> Reshape
( Int m, Int n, const AbstractDistMatrix<T>& A );

template<typename T>
void Reshape( Int m, Int n, const SparseMatrix<T>& A, SparseMatrix<T>& B );
template<typename T>
SparseMatrix<T> Reshape( Int m, Int n, const SparseMatrix<T>& A );

template<typename T>
void Reshape
( Int m, Int n, const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B );
template<typename T>
DistSparseMatrix<T> GetSubmatrix( Int m, Int n, const DistSparseMatrix<T>& A );

// Transform2x2 
// ============

// [a1,a2] := G [a1,a2], where G is 2x2
// ------------------------------------
template<typename T>
void Transform2x2
( const Matrix<T>& G,
        Matrix<T>& a1,
        Matrix<T>& a2 );
template<typename T>
void Transform2x2
( const AbstractDistMatrix<T>& G,
        AbstractDistMatrix<T>& a1,
        AbstractDistMatrix<T>& a2 );

// A([i1,i2],:) := G A([i1,i2],:), where G is 2x2
// ----------------------------------------------
template<typename T>
void Transform2x2Rows
( const Matrix<T>& G,
        Matrix<T>& A, Int i1, Int i2 );
template<typename T>
void Transform2x2Rows
( const AbstractDistMatrix<T>& G,
        AbstractDistMatrix<T>& A, Int i1, Int i2 );

// A(:,[j1,j2]) := A(:,[j1,j2]) G, where G is 2x2
// ----------------------------------------------
template<typename T>
void Transform2x2Cols
( const Matrix<T>& G,
        Matrix<T>& A, Int j1, Int j2 );
template<typename T>
void Transform2x2Cols
( const AbstractDistMatrix<T>& G,
        AbstractDistMatrix<T>& A, Int j1, Int j2 );

// TODO: SymmetricTransform2x2?

// Rotate (via Givens)
// ===================
// NOTE: BLAS calls this 'rot'

// [a1,a2] := [c s; -conj(s) c] [a1,a2]
// ------------------------------------
template<typename F>
void Rotate( Base<F> c, F s, Matrix<F>& a1, Matrix<F>& a2 );
template<typename F>
void Rotate
( Base<F> c, F s, AbstractDistMatrix<F>& a1, AbstractDistMatrix<F>& a2 );

// A([i1,i2],:) := [c s; -conj(s) c] A([i1,i2],:)
// ----------------------------------------------
template<typename F>
void RotateRows( Base<F> c, F s, Matrix<F>& A, Int i1, Int i2 );
template<typename F>
void RotateRows( Base<F> c, F s, AbstractDistMatrix<F>& A, Int i1, Int i2 );

// A(:,[j1,j2]) := A(:,[j1,j2]) [c s; -conj(s), c]
// -----------------------------------------------
template<typename F>
void RotateCols( Base<F> c, F s, Matrix<F>& A, Int j1, Int j2 );
template<typename F>
void RotateCols( Base<F> c, F s, AbstractDistMatrix<F>& A, Int j1, Int j2 );

// TODO: SymmetricRotation?

// Round
// =====
// Round each entry to the nearest integer
template<typename T>
void Round( Matrix<T>& A );
template<>
void Round( Matrix<Int>& A );
#ifdef EL_HAVE_MPC
template<>
void Round( Matrix<BigInt>& A );
#endif
template<typename T>
void Round( AbstractDistMatrix<T>& A );
template<typename T>
void Round( DistMultiVec<T>& A );
// TODO: Sparse matrix versions

// Scale
// =====
// TODO: Force S=T?
template<typename T,typename S>
void Scale( S alpha, Matrix<T>& A );
template<typename T,typename S>
void Scale( S alpha, AbstractDistMatrix<T>& A );
template<typename T,typename S>
void Scale( S alpha, SparseMatrix<T>& A );
template<typename T,typename S>
void Scale( S alpha, DistSparseMatrix<T>& A );
template<typename T,typename S>
void Scale( S alpha, DistMultiVec<T>& A );

template<typename Real,typename S,typename=EnableIf<IsReal<Real>>>
void Scale( S alpha, Matrix<Real>& AReal, Matrix<Real>& AImag );
template<typename Real,typename S,typename=EnableIf<IsReal<Real>>>
void Scale
( S alpha, AbstractDistMatrix<Real>& AReal, AbstractDistMatrix<Real>& AImag );

template<typename F>
void SafeScale
( Base<F> numerator, Base<F> denominator, Matrix<F>& A );
template<typename F>
void SafeScale
( Base<F> numerator, Base<F> denominator, AbstractDistMatrix<F>& A );
template<typename F>
void SafeScale
( Base<F> numerator, Base<F> denominator, SparseMatrix<F>& A );
template<typename F>
void SafeScale
( Base<F> numerator, Base<F> denominator, DistSparseMatrix<F>& A );
template<typename F>
void SafeScale
( Base<F> numerator, Base<F> denominator, DistMultiVec<F>& A );

template<typename F>
void SafeScaleTrapezoid
( Base<F> numerator, Base<F> denominator,
  UpperOrLower uplo, Matrix<F>& A, Int offset=0 );
template<typename F>
void SafeScaleTrapezoid
( Base<F> numerator, Base<F> denominator,
  UpperOrLower uplo, AbstractDistMatrix<F>& A, Int offset=0 );
template<typename F>
void SafeScaleTrapezoid
( Base<F> numerator, Base<F> denominator,
  UpperOrLower uplo, SparseMatrix<F>& A, Int offset=0 );
template<typename F>
void SafeScaleTrapezoid
( Base<F> numerator, Base<F> denominator,
  UpperOrLower uplo, DistSparseMatrix<F>& A, Int offset=0 );
template<typename F>
void SafeScaleTrapezoid
( Base<F> numerator, Base<F> denominator,
  UpperOrLower uplo, DistMultiVec<F>& A, Int offset=0 );

template<typename F>
void SafeScaleHermitianTridiag
( Base<F> numerator, Base<F> denominator, Matrix<Base<F>>& d, Matrix<F>& e );

// ScaleTrapezoid
// ==============
template<typename T,typename S>
void ScaleTrapezoid( S alpha, UpperOrLower uplo, Matrix<T>& A, Int offset=0 );
template<typename T,typename S>
void ScaleTrapezoid
( S alpha, UpperOrLower uplo, AbstractDistMatrix<T>& A, Int offset=0 );
template<typename T,typename S>
void ScaleTrapezoid
( S alpha, UpperOrLower uplo, SparseMatrix<T>& A, Int offset=0 );
template<typename T,typename S>
void ScaleTrapezoid
( S alpha, UpperOrLower uplo, DistSparseMatrix<T>& A, Int offset=0 );

// SetDiagonal
// ===========
template<typename T>
void SetDiagonal
( Matrix<T>& A, const Matrix<T>& d, Int offset=0 );
template<typename T>
void SetRealPartOfDiagonal
( Matrix<T>& A, const Matrix<Base<T>>& d, Int offset=0 );
template<typename T>
void SetImagPartOfDiagonal
( Matrix<T>& A, const Matrix<Base<T>>& d, Int offset=0 );

template<typename T,Dist U,Dist V>
void SetDiagonal
( DistMatrix<T,U,V>& A, const ElementalMatrix<T>& d, Int offset=0 );
template<typename T,Dist U,Dist V>
void SetRealPartOfDiagonal
( DistMatrix<T,U,V>& A, const ElementalMatrix<Base<T>>& d, Int offset=0 );
template<typename T,Dist U,Dist V>
void SetImagPartOfDiagonal
( DistMatrix<T,U,V>& A, const ElementalMatrix<Base<T>>& d, Int offset=0 );
// Versions which will work for ElementalMatrix but which make use of a
// manual dynamic dispatch
template<typename T>
void SetDiagonal
( ElementalMatrix<T>& A, const ElementalMatrix<T>& d, Int offset=0 );
template<typename T>
void SetRealPartOfDiagonal
( ElementalMatrix<T>& A, const ElementalMatrix<Base<T>>& d, Int offset=0 );
template<typename T>
void SetImagPartOfDiagonal
( ElementalMatrix<T>& A, const ElementalMatrix<Base<T>>& d, Int offset=0 );

// SetSubmatrix
// ============
template<typename T>
void SetSubmatrix
(       Matrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J, 
  const Matrix<T>& ASub );
template<typename T>
void SetSubmatrix
(       AbstractDistMatrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J, 
  const AbstractDistMatrix<T>& ASub );

// Swap
// ====
template<typename T>
void Swap( Orientation orientation, Matrix<T>& X, Matrix<T>& Y );
template<typename T>
void Swap
( Orientation orientation, AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y );

template<typename T>
void RowSwap( Matrix<T>& A, Int to, Int from );
template<typename T>
void RowSwap( AbstractDistMatrix<T>& A, Int to, Int from );

template<typename T>
void ColSwap( Matrix<T>& A, Int to, Int from );
template<typename T>
void ColSwap( AbstractDistMatrix<T>& A, Int to, Int from );

template<typename T>
void SymmetricSwap
( UpperOrLower uplo, Matrix<T>& A, Int to, Int from, bool conjugate=false );
template<typename T>
void SymmetricSwap
( UpperOrLower uplo, AbstractDistMatrix<T>& A, Int to, Int from,
  bool conjugate=false );

template<typename T>
void HermitianSwap( UpperOrLower uplo, Matrix<T>& A, Int to, Int from );
template<typename T>
void HermitianSwap
( UpperOrLower uplo, AbstractDistMatrix<T>& A, Int to, Int from );

// Symmetric2x2Inv
// ===============
template<typename F>
void Symmetric2x2Inv( UpperOrLower uplo, Matrix<F>& D, bool conjugate=false );

// Shift
// =====
template<typename T,typename S>
void Shift( Matrix<T>& A, S alpha );
template<typename T,typename S>
void Shift( AbstractDistMatrix<T>& A, S alpha );
template<typename T,typename S>
void Shift( DistMultiVec<T>& A, S alpha );

// ShiftDiagonal
// =============
template<typename T,typename S>
void ShiftDiagonal( Matrix<T>& A, S alpha, Int offset=0 );
template<typename T,typename S>
void ShiftDiagonal( AbstractDistMatrix<T>& A, S alpha, Int offset=0 );
template<typename T,typename S>
void ShiftDiagonal
( SparseMatrix<T>& A, S alpha, Int offset=0, bool existingDiag=false );
template<typename T,typename S>
void ShiftDiagonal
( DistSparseMatrix<T>& A, S alpha, Int offset=0, bool existingDiag=false );

// Transpose
// =========
template<typename T>
void Transpose
( const Matrix<T>& A,
        Matrix<T>& B,
  bool conjugate=false );
template<typename T>
void Transpose
( const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B,
  bool conjugate=false );
template<typename T>
void Transpose
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B,
  bool conjugate=false );
template<typename T>
void Transpose
( const AbstractDistMatrix<T>& A,
        AbstractDistMatrix<T>& B,
  bool conjugate=false );
template<typename T>
void Transpose
( const SparseMatrix<T>& A,
        SparseMatrix<T>& B,
  bool conjugate=false );
template<typename T>
void Transpose
( const DistSparseMatrix<T>& A,
        DistSparseMatrix<T>& B,
  bool conjugate=false );

// TransposeContract
// =================
template<typename T>
void TransposeContract
( const ElementalMatrix<T>& A, 
        ElementalMatrix<T>& B, bool conjugate=false );
template<typename T>
void TransposeContract
( const BlockMatrix<T>& A, 
        BlockMatrix<T>& B, bool conjugate=false );

// TransposeAxpy
// =============
template<typename T,typename S>
void TransposeAxpy
( S alpha, const Matrix<T>& X, Matrix<T>& Y, bool conjugate=false );
template<typename T,typename S>
void TransposeAxpy
( S alpha, const SparseMatrix<T>& X, SparseMatrix<T>& Y, bool conjugate=false );
template<typename T,typename S>
void TransposeAxpy
( S alpha, const ElementalMatrix<T>& X, ElementalMatrix<T>& Y, 
  bool conjugate=false );
template<typename T,typename S>
void TransposeAxpy
( S alpha, const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y, 
  bool conjugate=false );

// TransposeAxpyContract
// =====================
template<typename T>
void TransposeAxpyContract
( T alpha, const ElementalMatrix<T>& A, 
                 ElementalMatrix<T>& B, bool conjugate=false );
template<typename T>
void TransposeAxpyContract
( T alpha, const BlockMatrix<T>& A, 
                 BlockMatrix<T>& B, bool conjugate=false );

// UpdateDiagonal
// ==============
template<typename T>
void UpdateDiagonal
( Matrix<T>& A, T alpha, const Matrix<T>& d, Int offset=0 );
template<typename T>
void UpdateRealPartOfDiagonal
( Matrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset=0 );
template<typename T>
void UpdateImagPartOfDiagonal
( Matrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset=0 );

template<typename T,Dist U,Dist V>
void UpdateDiagonal
( DistMatrix<T,U,V>& A, T alpha, const ElementalMatrix<T>& d, 
  Int offset=0 );
template<typename T,Dist U,Dist V>
void UpdateRealPartOfDiagonal
( DistMatrix<T,U,V>& A, Base<T> alpha, const ElementalMatrix<Base<T>>& d, 
  Int offset=0 );
template<typename T,Dist U,Dist V>
void UpdateImagPartOfDiagonal
( DistMatrix<T,U,V>& A, Base<T> alpha, const ElementalMatrix<Base<T>>& d, 
  Int offset=0 );

template<typename T>
void UpdateDiagonal
( SparseMatrix<T>& A, T alpha, const Matrix<T>& d, Int offset=0, 
  bool diagExists=false );
template<typename T>
void UpdateRealPartOfDiagonal
( SparseMatrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset=0,
  bool diagExists=false );
template<typename T>
void UpdateImagPartOfDiagonal
( SparseMatrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset=0,
  bool diagExists=false );

template<typename T>
void UpdateDiagonal
( DistSparseMatrix<T>& A, T alpha, 
  const DistMultiVec<T>& d, Int offset=0, bool diagExists=false );
template<typename T>
void UpdateRealPartOfDiagonal
( DistSparseMatrix<T>& A, Base<T> alpha, 
  const DistMultiVec<Base<T>>& d, Int offset=0, bool diagExists=false );
template<typename T>
void UpdateImagPartOfDiagonal
( DistSparseMatrix<T>& A, Base<T> alpha, 
  const DistMultiVec<Base<T>>& d, Int offset=0, bool diagExists=false );

// UpdateMappedDiagonal
// ====================
template<typename T,typename S>
void UpdateMappedDiagonal
( Matrix<T>& A, const Matrix<S>& d, 
  function<void(T&,S)> func, Int offset=0 );
template<typename T,typename S,Dist U,Dist V>
void UpdateMappedDiagonal
( DistMatrix<T,U,V>& A, const ElementalMatrix<S>& d, 
  function<void(T&,S)> func, Int offset=0 );

template<typename T,typename S>
void UpdateMappedDiagonal
( SparseMatrix<T>& A, const Matrix<S>& d, 
  function<void(T&,S)> func, Int offset=0, bool diagExists=false );
template<typename T,typename S>
void UpdateMappedDiagonal
( DistSparseMatrix<T>& A, const DistMultiVec<S>& d, 
  function<void(T&,S)> func, Int offset=0, bool diagExists=false );

// UpdateSubmatrix
// ===============
template<typename T>
void UpdateSubmatrix
(       Matrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J, 
  T alpha, const Matrix<T>& ASub );
template<typename T>
void UpdateSubmatrix
(       AbstractDistMatrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J, 
  T alpha, const AbstractDistMatrix<T>& ASub );

// Zero
// ====
template<typename T>
void Zero( Matrix<T>& A );
template<typename T>
void Zero( AbstractDistMatrix<T>& A );
template<typename T>
void Zero( SparseMatrix<T>& A, bool clearMemory=true );
template<typename T>
void Zero( DistSparseMatrix<T>& A, bool clearMemory=true );
template<typename T>
void Zero( DistMultiVec<T>& A );

// Givens rotations
// ================
//
// Given phi and gamma, compute a Givens rotation such that
//
//  |       c   s | |   phi |  = | rho |, where c^2 + |s|^2 = 1
//  | -conj(s)  c | | gamma |    |  0  |
//
// This routine uses the stable approach suggested by Kahan and Demmel and
// returns the value rho.
//

template<typename Real>
Real Givens
( const Real& phi,
  const Real& gamma,
        Real& c,
        Real& s );
template<typename Real>
Complex<Real> Givens
( const Complex<Real>& phi,
  const Complex<Real>& gamma,
  Real& c,
  Complex<Real>& s );

} // namespace El

#endif // ifndef EL_BLAS1_DECL_HPP
