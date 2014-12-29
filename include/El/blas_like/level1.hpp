/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS1_HPP
#define EL_BLAS1_HPP

namespace El {

// Adjoint
// =======
template<typename T>
void Adjoint( const Matrix<T>& A, Matrix<T>& B );
template<typename T>
void Adjoint( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );
template<typename T>
void Adjoint
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );
template<typename T>
void Adjoint( const SparseMatrix<T>& A, SparseMatrix<T>& B );
template<typename T>
void Adjoint( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B );

namespace adjoint {

template<typename T,Dist U,Dist V>
void ColFilter
( const DistMatrix<T,V,Collect<U>()>& A, 
        DistMatrix<T,U,V           >& B );
template<typename T,Dist U,Dist V>
void ColFilter
( const BlockDistMatrix<T,V,Collect<U>()>& A, 
        BlockDistMatrix<T,U,V           >& B );

template<typename T,Dist U,Dist V>
void RowFilter
( const DistMatrix<T,Collect<V>(),U>& A,           
        DistMatrix<T,        U,   V>& B );
template<typename T,Dist U,Dist V>
void RowFilter
( const BlockDistMatrix<T,Collect<V>(),U>& A,           
        BlockDistMatrix<T,        U,   V>& B );

template<typename T,Dist U,Dist V>
void PartialColFilter
( const DistMatrix<T,V,Partial<U>()>& A, 
        DistMatrix<T,U,        V   >& B );
template<typename T,Dist U,Dist V>
void PartialColFilter
( const BlockDistMatrix<T,V,Partial<U>()>& A, 
        BlockDistMatrix<T,U,        V   >& B );

template<typename T,Dist U,Dist V>
void PartialRowFilter
( const DistMatrix<T,Partial<V>(),U>& A, 
        DistMatrix<T,U,           V>& B );
template<typename T,Dist U,Dist V>
void PartialRowFilter
( const BlockDistMatrix<T,Partial<V>(),U>& A, 
        BlockDistMatrix<T,U,           V>& B );

template<typename T,Dist U,Dist V>
void ColSumScatter
( const DistMatrix<T,V,Collect<U>()>& A, 
        DistMatrix<T,U,V           >& B );
template<typename T,Dist U,Dist V>
void ColSumScatter
( const BlockDistMatrix<T,V,Collect<U>()>& A, 
        BlockDistMatrix<T,U,V           >& B );

template<typename T,Dist U,Dist V>
void PartialColSumScatter
( const DistMatrix<T,V,Partial<U>()>& A,
        DistMatrix<T,U,V           >& B );
template<typename T,Dist U,Dist V>
void PartialColSumScatter
( const BlockDistMatrix<T,V,Partial<U>()>& A,
        BlockDistMatrix<T,U,V           >& B );

template<typename T,Dist U,Dist V>
void ColAllGather
( const DistMatrix<T,U,        V   >& A, 
        DistMatrix<T,V,Collect<U>()>& B );
template<typename T,Dist U,Dist V>
void ColAllGather
( const BlockDistMatrix<T,U,        V   >& A, 
        BlockDistMatrix<T,V,Collect<U>()>& B );

template<typename T,Dist U,Dist V>
void PartialColAllGather
( const DistMatrix<T,U,        V   >& A, 
        DistMatrix<T,V,Partial<U>()>& B );
template<typename T,Dist U,Dist V>
void PartialColAllGather
( const BlockDistMatrix<T,U,        V   >& A, 
        BlockDistMatrix<T,V,Partial<U>()>& B );

} // namespace adjoint

// AdjointAxpy
// =============
template<typename T,typename S>
void AdjointAxpy
( S alpha, const Matrix<T>& X, Matrix<T>& Y );
template<typename T,typename S>
void AdjointAxpy
( S alpha, const SparseMatrix<T>& X, SparseMatrix<T>& Y );
template<typename T,typename S>
void AdjointAxpy
( S alpha, const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y );
template<typename T,typename S>
void AdjointAxpy
( S alpha, const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y );

namespace adjoint_axpy {

template<typename T,Dist U,Dist V>
void ColSumScatter
( T alpha, const DistMatrix<T,V,Collect<U>()>& A,
                 DistMatrix<T,U,        V   >& B );
template<typename T,Dist U,Dist V>
void ColSumScatter
( T alpha, const BlockDistMatrix<T,V,Collect<U>()>& A,
                 BlockDistMatrix<T,U,        V   >& B );

template<typename T,Dist U,Dist V>
void PartialColSumScatter
( T alpha, const DistMatrix<T,V,Partial<U>()>& A,
                 DistMatrix<T,U,        V   >& B );
template<typename T,Dist U,Dist V>
void PartialColSumScatter
( T alpha, const BlockDistMatrix<T,V,Partial<U>()>& A,
                 BlockDistMatrix<T,U,        V   >& B );

} // namespace adjoint_axpy

// AllReduce
// =========
// TODO: Matrix<T> version?
template<typename T>
void AllReduce( AbstractDistMatrix<T>& A, mpi::Comm comm, mpi::Op op=mpi::SUM );
template<typename T>
void AllReduce( AbstractBlockDistMatrix<T>& A, mpi::Comm comm, mpi::Op op=mpi::SUM );

// Axpy
// ====
template<typename T,typename S>
void Axpy( S alpha, const Matrix<T>& X, Matrix<T>& Y );
template<typename T,typename S>
void Axpy( S alpha, const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y );
template<typename T,typename S>
void Axpy( S alpha, const DistMultiVec<T>& X, DistMultiVec<T>& Y );
template<typename T,typename S>
void Axpy( S alpha, const SparseMatrix<T>& X, SparseMatrix<T>& Y );
template<typename T,typename S>
void Axpy( S alpha, const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y );

namespace axpy {

template<typename T,Dist U,Dist V>
void SumScatter
( T alpha, const DistMatrix<T,Collect<U>(),Collect<V>()>& A, 
                 DistMatrix<T,        U,           V   >& B );
template<typename T,Dist U,Dist V>
void SumScatter
( T alpha, const BlockDistMatrix<T,Collect<U>(),Collect<V>()>& A, 
                 BlockDistMatrix<T,        U,           V   >& B );

template<typename T,Dist U,Dist V>
void ColSumScatter
( T alpha, const DistMatrix<T,Collect<U>(),V>& A, 
                 DistMatrix<T,        U,   V>& B );
template<typename T,Dist U,Dist V>
void ColSumScatter
( T alpha, const BlockDistMatrix<T,Collect<U>(),V>& A, 
                 BlockDistMatrix<T,        U,   V>& B );

template<typename T,Dist U,Dist V>
void RowSumScatter
( T alpha, const DistMatrix<T,U,Collect<V>()>& A, 
                 DistMatrix<T,U,        V   >& B );
template<typename T,Dist U,Dist V>
void RowSumScatter
( T alpha, const BlockDistMatrix<T,U,Collect<V>()>& A, 
                 BlockDistMatrix<T,U,        V   >& B );

template<typename T,Dist U,Dist V>
void PartialColSumScatter
( T alpha, const DistMatrix<T,Partial<U>(),V>& A, 
                 DistMatrix<T,        U,   V>& B );
template<typename T,Dist U,Dist V>
void PartialColSumScatter
( T alpha, const BlockDistMatrix<T,Partial<U>(),V>& A, 
                 BlockDistMatrix<T,        U,   V>& B );

template<typename T,Dist U,Dist V>
void PartialRowSumScatter
( T alpha, const DistMatrix<T,U,Partial<V>()>& A, 
                 DistMatrix<T,U,        V   >& B );
template<typename T,Dist U,Dist V>
void PartialRowSumScatter
( T alpha, const BlockDistMatrix<T,U,Partial<V>()>& A, 
                 BlockDistMatrix<T,U,        V   >& B );

namespace util {

template<typename T>
void InterleaveMatrixUpdate
( T alpha, Int localHeight, Int localWidth,
  const T* A, Int colStrideA, Int rowStrideA,
        T* B, Int colStrideB, Int rowStrideB );

} // namespace util

} // namespace axpy

// AxpyTrapezoid
// =============
template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alpha, const Matrix<T>& X, Matrix<T>& Y, Int offset=0 );
template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alpha,
  const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y, Int offset=0 );
template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alpha, 
  const SparseMatrix<T>& X, SparseMatrix<T>& Y, Int offset=0 );
template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alpha, 
  const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y, Int offset=0 );

// Broadcast
// =========
// TODO: Matrix<T> version?
template<typename T>
void Broadcast( AbstractDistMatrix<T>& A, mpi::Comm comm, Int rank=0 );
template<typename T>
void Broadcast( AbstractBlockDistMatrix<T>& A, mpi::Comm comm, Int rank=0 );

// Column norms
// ============
template<typename F>
void ColumnNorms( const Matrix<F>& X, Matrix<Base<F>>& norms );
template<typename Real>
void ColumnNorms
( const Matrix<Real>& XReal, const Matrix<Real>& XImag, 
  Matrix<Real>& norms );

template<typename F>
void ColumnNorms
( const AbstractDistMatrix<F>& X, Matrix<Base<F>>& norms );
template<typename F,Dist U,Dist V>
void ColumnNorms
( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,V,STAR>& norms );

template<typename Real>
void ColumnNorms
( const AbstractDistMatrix<Real>& XReal, const AbstractDistMatrix<Real>& XImag, 
  Matrix<Real>& norms );
template<typename Real,Dist U,Dist V>
void ColumnNorms
( const DistMatrix<Real,U,V>& XReal, const DistMatrix<Real,U,V>& XImag, 
  DistMatrix<Real,V,STAR>& norms );

template<typename F>
void ColumnNorms( const DistMultiVec<F>& X, Matrix<Base<F>>& norms );
template<typename Real>
void ColumnNorms
( const DistMultiVec<Real>& XReal, const DistMultiVec<Real>& XImag, 
  Matrix<Real>& norms );

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
void Conjugate( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );

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
( Matrix<T>& A, 
  const std::vector<Int>& rowInd, const std::vector<Int>& colInd );
template<typename T>
void ConjugateSubmatrix
( AbstractDistMatrix<T>& A, 
  const std::vector<Int>& rowInd, const std::vector<Int>& colInd );

// Copy
// ====
template<typename T>
void Copy( const Matrix<T>& A, Matrix<T>& B );
// TODO: A detailed description of which conversions are instantiated 
template<typename S,typename T>
void Copy( const Matrix<S>& A, Matrix<T>& B );
template<typename S,typename T>
void Copy( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B );
template<typename S,typename T>
void Copy( const AbstractBlockDistMatrix<S>& A, AbstractBlockDistMatrix<T>& B );

void Copy( const Graph& A, Graph& B );
void Copy( const Graph& A, DistGraph& B );
void Copy( const DistGraph& A, Graph& B );
void Copy( const DistGraph& A, DistGraph& B );

void CopyFromRoot( const DistGraph& distGraph, Graph& graph );
void CopyFromNonRoot( const DistGraph& distGraph, Int root=0 );

template<typename T>
void Copy( const SparseMatrix<T>& A, SparseMatrix<T>& B );
template<typename S,typename T>
void Copy( const SparseMatrix<S>& A, SparseMatrix<T>& B );
template<typename S,typename T>
void Copy( const SparseMatrix<S>& A, Matrix<T>& B );
template<typename T>
void Copy( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B );
template<typename S,typename T>
void Copy( const DistSparseMatrix<S>& A, AbstractDistMatrix<T>& B );
template<typename T>
void CopyFromRoot( const DistSparseMatrix<T>& ADist, SparseMatrix<T>& A );
template<typename T>
void CopyFromNonRoot( const DistSparseMatrix<T>& ADist, Int root=0 );

template<typename T>
void Copy( const DistMultiVec<T>& A, DistMultiVec<T>& B );
template<typename T>
void CopyFromRoot( const DistMultiVec<T>& XDist, Matrix<T>& X );
template<typename T>
void CopyFromNonRoot( const DistMultiVec<T>& XDist, Int root=0 );

namespace copy {

template<typename T,Dist U,Dist V>
void Translate( const DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B );
template<typename T,Dist U,Dist V>
void Translate( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,U,V>& B );

template<typename T>
void TranslateBetweenGrids
( const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B );
template<typename T>
void TranslateBetweenGrids
( const DistMatrix<T,STAR,STAR>& A, DistMatrix<T,STAR,STAR>& B );
// The fallback case that simply throws an exception
template<typename T,Dist U,Dist V>
void TranslateBetweenGrids
( const DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B );

// NOTE: Only instantiated for (U,V)=(MC,MR) and (U,V)=(MR,MC)
template<typename T,Dist U,Dist V>
void ColwiseVectorExchange
( const DistMatrix<T,ProductDist<U,V>(),STAR>& A,
        DistMatrix<T,ProductDist<V,U>(),STAR>& B );
template<typename T,Dist U,Dist V>
void RowwiseVectorExchange
( const DistMatrix<T,STAR,ProductDist<U,V>()>& A,
        DistMatrix<T,STAR,ProductDist<V,U>()>& B );

// NOTE: Only instantiated for (U,V)=(MC,MR) and (U,V)=(MR,MC)
template<typename T,Dist U,Dist V>
void TransposeDist( const DistMatrix<T,U,V>& A, DistMatrix<T,V,U>& B );

template<typename T,Dist U,Dist V>
void Filter
( const DistMatrix<T,Collect<U>(),Collect<V>()>& A,
        DistMatrix<T,        U,           V   >& B );
template<typename T,Dist U,Dist V>
void Filter
( const BlockDistMatrix<T,Collect<U>(),Collect<V>()>& A,
        BlockDistMatrix<T,        U,           V   >& B );

template<typename T,Dist U,Dist V>
void ColFilter
( const DistMatrix<T,Collect<U>(),V>& A,
        DistMatrix<T,        U,   V>& B );
template<typename T,Dist U,Dist V>
void ColFilter
( const BlockDistMatrix<T,Collect<U>(),V>& A,
        BlockDistMatrix<T,        U,   V>& B );

template<typename T,Dist U,Dist V>
void RowFilter
( const DistMatrix<T,U,Collect<V>()>& A,
        DistMatrix<T,U,        V   >& B );
template<typename T,Dist U,Dist V>
void RowFilter
( const BlockDistMatrix<T,U,Collect<V>()>& A,
        BlockDistMatrix<T,U,        V   >& B );

template<typename T,Dist U,Dist V>
void PartialColFilter
( const DistMatrix<T,Partial<U>(),V>& A,
        DistMatrix<T,        U,   V>& B );
template<typename T,Dist U,Dist V>
void PartialColFilter
( const BlockDistMatrix<T,Partial<U>(),V>& A,
        BlockDistMatrix<T,        U,   V>& B );

template<typename T,Dist U,Dist V>
void PartialRowFilter
( const DistMatrix<T,U,Partial<V>()>& A,
        DistMatrix<T,U,        V   >& B );
template<typename T,Dist U,Dist V>
void PartialRowFilter
( const BlockDistMatrix<T,U,Partial<V>()>& A,
        BlockDistMatrix<T,U,        V   >& B );

template<typename T,Dist U,Dist V>
void AllGather
( const DistMatrix<T,        U,           V   >& A, 
        DistMatrix<T,Collect<U>(),Collect<V>()>& B );
template<typename T,Dist U,Dist V>
void AllGather
( const BlockDistMatrix<T,        U,           V   >& A, 
        BlockDistMatrix<T,Collect<U>(),Collect<V>()>& B );

template<typename T,Dist U,Dist V>
void ColAllGather
( const DistMatrix<T,        U,   V>& A,
        DistMatrix<T,Collect<U>(),V>& B );
template<typename T,Dist U,Dist V>
void ColAllGather
( const BlockDistMatrix<T,        U,   V>& A,
        BlockDistMatrix<T,Collect<U>(),V>& B );

template<typename T,Dist U,Dist V>
void RowAllGather
( const DistMatrix<T,U,        V   >& A,
        DistMatrix<T,U,Collect<V>()>& B );
template<typename T,Dist U,Dist V>
void RowAllGather
( const BlockDistMatrix<T,U,        V   >& A,
        BlockDistMatrix<T,U,Collect<V>()>& B );

template<typename T,Dist U,Dist V>
void PartialColAllGather
( const DistMatrix<T,        U,   V>& A,
        DistMatrix<T,Partial<U>(),V>& B );
template<typename T,Dist U,Dist V>
void PartialColAllGather
( const BlockDistMatrix<T,        U,   V>& A,
        BlockDistMatrix<T,Partial<U>(),V>& B );

template<typename T,Dist U,Dist V>
void PartialRowAllGather
( const DistMatrix<T,U,        V   >& A,
        DistMatrix<T,U,Partial<V>()>& B );
template<typename T,Dist U,Dist V>
void PartialRowAllGather
( const BlockDistMatrix<T,U,        V   >& A,
        BlockDistMatrix<T,U,Partial<V>()>& B );

template<typename T,Dist U,Dist V>
void ColAllToAllDemote
( const DistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& A,
        DistMatrix<T,        U,                     V   >& B );
template<typename T,Dist U,Dist V>
void ColAllToAllDemote
( const BlockDistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& A,
        BlockDistMatrix<T,        U,                     V   >& B );

template<typename T,Dist U,Dist V>
void RowAllToAllDemote
( const DistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& A,
        DistMatrix<T,                U,             V   >& B );
template<typename T,Dist U,Dist V>
void RowAllToAllDemote
( const BlockDistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& A,
        BlockDistMatrix<T,                U,             V   >& B );

template<typename T,Dist U,Dist V>
void ColAllToAllPromote
( const DistMatrix<T,        U,                     V   >& A,
        DistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& B );
template<typename T,Dist U,Dist V>
void ColAllToAllPromote
( const BlockDistMatrix<T,        U,                     V   >& A,
        BlockDistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& B );

template<typename T,Dist U,Dist V>
void RowAllToAllPromote
( const DistMatrix<T,                U,             V   >& A,
        DistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& B );
template<typename T,Dist U,Dist V>
void RowAllToAllPromote
( const BlockDistMatrix<T,                U,             V   >& A,
        BlockDistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& B );

template<typename T,Dist U,Dist V>
void SumScatter
( const DistMatrix<T,Collect<U>(),Collect<V>()>& A, 
        DistMatrix<T,        U,           V   >& B );
template<typename T,Dist U,Dist V>
void SumScatter
( const BlockDistMatrix<T,Collect<U>(),Collect<V>()>& A, 
        BlockDistMatrix<T,        U,           V   >& B );

template<typename T,Dist U,Dist V>
void ColSumScatter
( const DistMatrix<T,Collect<U>(),V>& A, 
        DistMatrix<T,        U,   V>& B );
template<typename T,Dist U,Dist V>
void ColSumScatter
( const BlockDistMatrix<T,Collect<U>(),V>& A, 
        BlockDistMatrix<T,        U,   V>& B );

template<typename T,Dist U,Dist V>
void RowSumScatter
( const DistMatrix<T,U,Collect<V>()>& A, 
        DistMatrix<T,U,        V   >& B );
template<typename T,Dist U,Dist V>
void RowSumScatter
( const BlockDistMatrix<T,U,Collect<V>()>& A, 
        BlockDistMatrix<T,U,        V   >& B );

template<typename T,Dist U,Dist V>
void PartialColSumScatter
( const DistMatrix<T,Partial<U>(),V>& A, 
        DistMatrix<T,        U,   V>& B );
template<typename T,Dist U,Dist V>
void PartialColSumScatter
( const BlockDistMatrix<T,Partial<U>(),V>& A, 
        BlockDistMatrix<T,        U,   V>& B );

template<typename T,Dist U,Dist V>
void PartialRowSumScatter
( const DistMatrix<T,U,Partial<V>()>& A, 
        DistMatrix<T,U,        V   >& B );
template<typename T,Dist U,Dist V>
void PartialRowSumScatter
( const BlockDistMatrix<T,U,Partial<V>()>& A, 
        BlockDistMatrix<T,U,        V   >& B );

template<typename T>
void Gather
( const AbstractDistMatrix<T>& A,
        DistMatrix<T,CIRC,CIRC>& B );
template<typename T>
void Gather
( const AbstractBlockDistMatrix<T>& A,
        BlockDistMatrix<T,CIRC,CIRC>& B );

template<typename T>
void Scatter
( const DistMatrix<T,CIRC,CIRC>& A,
        AbstractDistMatrix<T>& B );
template<typename T>
void Scatter
( const BlockDistMatrix<T,CIRC,CIRC>& A,
        AbstractBlockDistMatrix<T>& B );

template<typename T>
void Scatter
( const DistMatrix<T,CIRC,CIRC>& A,
        DistMatrix<T,STAR,STAR>& B );
template<typename T>
void Scatter
( const BlockDistMatrix<T,CIRC,CIRC>& A,
        BlockDistMatrix<T,STAR,STAR>& B );

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

template<typename TDiag,typename T,Dist U,Dist V>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<TDiag>& d, DistMatrix<T,U,V>& A );

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
( Orientation orientation,
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
  const AbstractDistMatrix<TDiag>& d, DistMatrix<T,U,V>& A, Int offset=0 );

template<typename TDiag,typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<TDiag>& d, AbstractDistMatrix<T>& A, Int offset=0 );

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

template<typename FDiag,typename F,Dist U,Dist V>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<FDiag>& d, DistMatrix<F,U,V>& A,
  bool checkIfSingular=true );

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<FDiag>& d, AbstractDistMatrix<F>& A,
  bool checkIfSingular=true );

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<FDiag>& d, SparseMatrix<F>& A, 
  bool checkIfSingular=true );

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<FDiag>& d, DistSparseMatrix<F>& A, 
  bool checkIfSingular=true );

template<typename FDiag,typename F>
void DiagonalSolve
( Orientation orientation,
  const DistMultiVec<FDiag>& d, DistMultiVec<F>& X, 
  bool checkIfSingular=true );

// Dot
// ===
template<typename T>
T Dot( const Matrix<T>& A, const Matrix<T>& B );
template<typename T>
T Dot( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B );
template<typename T>
T Dot( const DistMultiVec<T>& A, const DistMultiVec<T>& B );

// Dotu
// ====
template<typename T>
T Dotu( const Matrix<T>& A, const Matrix<T>& B );
template<typename T>
T Dotu( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B );
template<typename T>
T Dotu( const DistMultiVec<T>& A, const DistMultiVec<T>& B );

// EntrywiseFill
// =============
template<typename T>
void EntrywiseFill( Matrix<T>& A, std::function<T(void)> func );
template<typename T>
void EntrywiseFill( AbstractDistMatrix<T>& A, std::function<T(void)> func );
template<typename T>
void EntrywiseFill
( AbstractBlockDistMatrix<T>& A, std::function<T(void)> func );
template<typename T>
void EntrywiseFill( DistMultiVec<T>& A, std::function<T(void)> func );

// EntrywiseMap
// ============
template<typename T>
void EntrywiseMap( Matrix<T>& A, std::function<T(T)> func );
template<typename T>
void EntrywiseMap( SparseMatrix<T>& A, std::function<T(T)> func );
template<typename T>
void EntrywiseMap( AbstractDistMatrix<T>& A, std::function<T(T)> func );
template<typename T>
void EntrywiseMap( AbstractBlockDistMatrix<T>& A, std::function<T(T)> func );
template<typename T>
void EntrywiseMap( DistSparseMatrix<T>& A, std::function<T(T)> func );
template<typename T>
void EntrywiseMap( DistMultiVec<T>& A, std::function<T(T)> func );

template<typename S,typename T>
void EntrywiseMap
( const Matrix<S>& A, Matrix<T>& B, std::function<T(S)> func );
template<typename S,typename T>
void EntrywiseMap
( const SparseMatrix<S>& A, SparseMatrix<T>& B, std::function<T(S)> func );
template<typename S,typename T>
void EntrywiseMap
( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B, 
  std::function<T(S)> func );
template<typename S,typename T>
void EntrywiseMap
( const AbstractBlockDistMatrix<S>& A, AbstractBlockDistMatrix<T>& B, 
  std::function<T(S)> func );
template<typename S,typename T>
void EntrywiseMap
( const DistSparseMatrix<S>& A, DistSparseMatrix<T>& B, 
  std::function<T(S)> func );
template<typename S,typename T>
void EntrywiseMap
( const DistMultiVec<S>& A, DistMultiVec<T>& B, 
  std::function<T(S)> func );

// Fill
// ====
template<typename T>
void Fill( Matrix<T>& A, T alpha );
template<typename T>
void Fill( AbstractDistMatrix<T>& A, T alpha );
template<typename T>
void Fill( AbstractBlockDistMatrix<T>& A, T alpha );
template<typename T>
void Fill( DistMultiVec<T>& A, T alpha );

// FillDiagonal
// ============
template<typename T,typename S>
void FillDiagonal( Matrix<T>& A, S alpha, Int offset=0 );
template<typename T,typename S>
void FillDiagonal( AbstractDistMatrix<T>& A, S alpha, Int offset=0 );
template<typename T,typename S>
void FillDiagonal( AbstractBlockDistMatrix<T>& A, S alpha, Int offset=0 );

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
  AbstractDistMatrix<T>& d, Int offset=0 );
template<typename T,Dist U,Dist V>
void GetRealPartOfDiagonal
( const DistMatrix<T,U,V>& A, 
  AbstractDistMatrix<Base<T>>& d, Int offset=0 );
template<typename T,Dist U,Dist V>
void GetImagPartOfDiagonal
( const DistMatrix<T,U,V>& A, 
  AbstractDistMatrix<Base<T>>& d, Int offset=0 );
// Versions which will work for AbstractDistMatrix but which make use of a
// manual dynamic dispatch
template<typename T>
void GetDiagonal
( const AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& d, Int offset=0 );
template<typename T>
void GetRealPartOfDiagonal
( const AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<Base<T>>& d, Int offset=0 );
template<typename T>
void GetImagPartOfDiagonal
( const AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<Base<T>>& d, Int offset=0 );

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
( const Matrix<T>& A, Matrix<S>& d, 
  std::function<S(T)> func, Int offset=0 );
template<typename T,typename S,Dist U,Dist V>
void GetMappedDiagonal
( const DistMatrix<T,U,V>& A, AbstractDistMatrix<S>& d, 
  std::function<S(T)> func, Int offset=0 );

// GetSubmatrix
// ============
template<typename T>
void GetSubmatrix
( const Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        Matrix<T>& ASub );
template<typename T>
void GetRealPartOfSubmatrix
( const Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        Matrix<Base<T>>& ASub );
template<typename T>
void GetImagPartOfSubmatrix
( const Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        Matrix<Base<T>>& ASub );

template<typename T>
Matrix<T> GetSubmatrix
( const Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J );
template<typename T>
Matrix<Base<T>> GetRealPartOfSubmatrix
( const Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J );
template<typename T>
Matrix<Base<T>> GetImagPartOfSubmatrix
( const Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J );

template<typename T>
void GetSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        AbstractDistMatrix<T>& ASub );
template<typename T>
void GetRealPartOfSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        AbstractDistMatrix<Base<T>>& ASub );
template<typename T>
void GetImagPartOfSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        AbstractDistMatrix<Base<T>>& ASub );

template<typename T>
DistMatrix<T,STAR,STAR> GetSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J );
template<typename T>
DistMatrix<Base<T>,STAR,STAR> GetRealPartOfSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J );
template<typename T>
DistMatrix<Base<T>,STAR,STAR> GetImagPartOfSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J );

// Hadamard
// ========
template<typename T>
void Hadamard( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C );
template<typename T>
void Hadamard
( const AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& C );

// HilbertSchmidt
// ==============
template<typename T>
T HilbertSchmidt( const Matrix<T>& A, const Matrix<T>& B );
template<typename T>
T HilbertSchmidt
( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& C );
template<typename T>
T HilbertSchmidt( const DistMultiVec<T>& A, const DistMultiVec<T>& B );

// Imaginary part
// ==============
template<typename T>
void ImagPart
( const Matrix<T>& A, Matrix<Base<T>>& AImag );
template<typename T>
void ImagPart
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<Base<T>>& AImag );

// IndexDependentFill
// ==================
template<typename T>
void IndexDependentFill( Matrix<T>& A, std::function<T(Int,Int)> func );
template<typename T>
void IndexDependentFill
( AbstractDistMatrix<T>& A, std::function<T(Int,Int)> func );
template<typename T>
void IndexDependentFill
( AbstractBlockDistMatrix<T>& A, std::function<T(Int,Int)> func );

// IndexDependentMap
// =================
template<typename T>
void IndexDependentMap( Matrix<T>& A, std::function<T(Int,Int,T)> func );
template<typename T>
void IndexDependentMap
( AbstractDistMatrix<T>& A, std::function<T(Int,Int,T)> func );
template<typename T>
void IndexDependentMap
( AbstractBlockDistMatrix<T>& A, std::function<T(Int,Int,T)> func );

template<typename S,typename T>
void IndexDependentMap
( const Matrix<S>& A, Matrix<T>& B, std::function<T(Int,Int,S)> func );
template<typename S,typename T>
void IndexDependentMap
( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B,
  std::function<T(Int,Int,S)> func );
template<typename S,typename T>
void IndexDependentMap
( const AbstractBlockDistMatrix<S>& A, AbstractBlockDistMatrix<T>& B,
  std::function<T(Int,Int,S)> func );

// MakeHermitian
// =============
template<typename T>
void MakeHermitian( UpperOrLower uplo, Matrix<T>& A );
template<typename T>
void MakeHermitian( UpperOrLower uplo, AbstractDistMatrix<T>& A );

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
template<typename T>
void MakeDiagonalReal( AbstractBlockDistMatrix<T>& A, Int offset=0 );

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
( Matrix<T>& A, 
  const std::vector<Int>& rowInd, const std::vector<Int>& colInd );
template<typename T>
void MakeSubmatrixReal
( AbstractDistMatrix<T>& A, 
  const std::vector<Int>& rowInd, const std::vector<Int>& colInd );

// MakeSymmetric
// =============
template<typename T>
void MakeSymmetric( UpperOrLower uplo, Matrix<T>& A, bool conjugate=false );
template<typename T>
void MakeSymmetric
( UpperOrLower uplo, AbstractDistMatrix<T>& A, bool conjugate=false );

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
void MakeTrapezoidal
( UpperOrLower uplo, AbstractBlockDistMatrix<T>& A, Int offset=0 );

template<typename T>
void MakeTrapezoidal( UpperOrLower uplo, SparseMatrix<T>& A, Int offset=0 );
template<typename T>
void MakeTrapezoidal( UpperOrLower uplo, DistSparseMatrix<T>& A, Int offset=0 );

// Max
// ===
template<typename Real>
ValueIntPair<Real> Max( const Matrix<Real>& A );
template<typename Real>
ValueIntPair<Real> Max( const AbstractDistMatrix<Real>& A );

template<typename Real>
ValueIntPair<Real> SymmetricMax( UpperOrLower uplo, const Matrix<Real>& A );
template<typename Real>
ValueIntPair<Real>
SymmetricMax( UpperOrLower uplo, const AbstractDistMatrix<Real>& A );

template<typename Real>
ValueInt<Real> VectorMax( const Matrix<Real>& x );
template<typename Real>
ValueInt<Real> VectorMax( const AbstractDistMatrix<Real>& x );
template<typename Real>
ValueInt<Real> VectorMax( const DistMultiVec<Real>& x );

// MaxAbs
// ======
template<typename F>
ValueInt<Base<F>> VectorMaxAbs( const Matrix<F>& x );
template<typename F>
ValueInt<Base<F>> VectorMaxAbs( const AbstractDistMatrix<F>& x );

template<typename F>
ValueIntPair<Base<F>> MaxAbs( const Matrix<F>& A );
template<typename F>
ValueIntPair<Base<F>> MaxAbs( const AbstractDistMatrix<F>& A );

template<typename F>
ValueIntPair<Base<F>> SymmetricMaxAbs( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
ValueIntPair<Base<F>> SymmetricMaxAbs
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Min
// ===
template<typename Real>
ValueInt<Real> VectorMin( const Matrix<Real>& x );
template<typename Real>
ValueInt<Real> VectorMin( const AbstractDistMatrix<Real>& x );
template<typename Real>
ValueInt<Real> VectorMin( const DistMultiVec<Real>& x );

template<typename Real>
ValueIntPair<Real> Min( const Matrix<Real>& A );
template<typename Real>
ValueIntPair<Real> Min( const AbstractDistMatrix<Real>& A );

template<typename Real>
ValueIntPair<Real> SymmetricMin( UpperOrLower uplo, const Matrix<Real>& A );
template<typename Real>
ValueIntPair<Real>
SymmetricMin( UpperOrLower uplo, const AbstractDistMatrix<Real>& A );

// MinAbs
// ======
template<typename F>
ValueInt<Base<F>> VectorMinAbs( const Matrix<F>& x );
template<typename F>
ValueInt<Base<F>> VectorMinAbs( const AbstractDistMatrix<F>& x );

template<typename F>
ValueIntPair<Base<F>> MinAbs( const Matrix<F>& A );
template<typename F>
ValueIntPair<Base<F>> MinAbs( const AbstractDistMatrix<F>& A );

template<typename F>
ValueIntPair<Base<F>> SymmetricMinAbs( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
ValueIntPair<Base<F>>
SymmetricMinAbs( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

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
// TODO: Switch to full AbstractDistMatrix interface
template<typename F,typename FMain,Dist U,Dist V>
void QuasiDiagonalScale
( LeftOrRight side, UpperOrLower uplo,
  const AbstractDistMatrix<FMain>& d, const AbstractDistMatrix<F>& dSub,
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
// TODO: Switch to full AbstractDistMatrix interface
template<typename F,typename FMain,Dist U,Dist V>
void
QuasiDiagonalSolve
( LeftOrRight side, UpperOrLower uplo,
  const AbstractDistMatrix<FMain>& d, const AbstractDistMatrix<F>& dSub,
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
void RealPart
( const Matrix<T>& A, Matrix<Base<T>>& AReal );
template<typename T>
void RealPart
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<Base<T>>& AReal );

// Scale
// =====
template<typename T,typename S>
void Scale( S alpha, Matrix<T>& A );
template<typename T,typename S>
void Scale( S alpha, AbstractDistMatrix<T>& A );
template<typename T,typename S>
void Scale( S alpha, AbstractBlockDistMatrix<T>& A );
template<typename T,typename S>
void Scale( S alpha, SparseMatrix<T>& A );
template<typename T,typename S>
void Scale( S alpha, DistSparseMatrix<T>& A );
template<typename T,typename S>
void Scale( S alpha, DistMultiVec<T>& A );

template<typename Real,typename S>
void Scale( S alpha, Matrix<Real>& AReal, Matrix<Real>& AImag );
template<typename Real,typename S>
void Scale
( S alpha, AbstractDistMatrix<Real>& AReal, AbstractDistMatrix<Real>& AImag );
template<typename Real,typename S>
void Scale
( S alpha, AbstractBlockDistMatrix<Real>& AReal,
           AbstractBlockDistMatrix<Real>& AImag );

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
( DistMatrix<T,U,V>& A, const AbstractDistMatrix<T>& d, Int offset=0 );
template<typename T,Dist U,Dist V>
void SetRealPartOfDiagonal
( DistMatrix<T,U,V>& A, const AbstractDistMatrix<Base<T>>& d, Int offset=0 );
template<typename T,Dist U,Dist V>
void SetImagPartOfDiagonal
( DistMatrix<T,U,V>& A, const AbstractDistMatrix<Base<T>>& d, Int offset=0 );
// Versions which will work for AbstractDistMatrix but which make use of a
// manual dynamic dispatch
template<typename T>
void SetDiagonal
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& d, Int offset=0 );
template<typename T>
void SetRealPartOfDiagonal
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Base<T>>& d, 
  Int offset=0 );
template<typename T>
void SetImagPartOfDiagonal
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Base<T>>& d, 
  Int offset=0 );

// SetSubmatrix
// ============
template<typename T>
void SetSubmatrix
(       Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const Matrix<T>& ASub );
template<typename T>
void SetRealPartOfSubmatrix
(       Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const Matrix<Base<T>>& ASub );
template<typename T>
void SetImagPartOfSubmatrix
(       Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const Matrix<Base<T>>& ASub );

template<typename T>
void SetSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const AbstractDistMatrix<T>& ASub );
template<typename T>
void SetRealPartOfSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const AbstractDistMatrix<Base<T>>& ASub );
template<typename T>
void SetImagPartOfSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const AbstractDistMatrix<Base<T>>& ASub );

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

// Symmetric2x2Scale
// =================
template<typename F>
void
Symmetric2x2Scale
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, Matrix<F>& A, bool conjugate=false );
template<typename F>
void
Symmetric2x2Scale
( LeftOrRight side, UpperOrLower uplo,
  const AbstractDistMatrix<F>& D, AbstractDistMatrix<F>& A,
  bool conjugate=false );

template<typename F>
void
FirstHalfOfSymmetric2x2Scale
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, Matrix<F>& a1, const Matrix<F>& a2,
  bool conjugate=false );

template<typename F>
void
SecondHalfOfSymmetric2x2Scale
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, const Matrix<F>& a1, Matrix<F>& a2,
  bool conjugate=false );

// Symmetric2x2Solve
// =================
template<typename F>
void
Symmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, Matrix<F>& A, bool conjugate=false );
template<typename F>
void
Symmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const AbstractDistMatrix<F>& D, AbstractDistMatrix<F>& A,
  bool conjugate=false );

template<typename F>
void
FirstHalfOfSymmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, Matrix<F>& a1, const Matrix<F>& a2,
  bool conjugate=false );

template<typename F>
void
SecondHalfOfSymmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, const Matrix<F>& a1, Matrix<F>& a2,
  bool conjugate=false );

// Shift
// =====
template<typename T,typename S>
void Shift( Matrix<T>& A, S alpha );
template<typename T,typename S>
void Shift( AbstractDistMatrix<T>& A, S alpha );
template<typename T,typename S>
void Shift( AbstractBlockDistMatrix<T>& A, S alpha );
template<typename T,typename S>
void Shift( DistMultiVec<T>& A, S alpha );

// ShiftDiagonal
// =============
template<typename T,typename S>
void ShiftDiagonal( Matrix<T>& A, S alpha, Int offset=0 );
template<typename T,typename S>
void ShiftDiagonal( AbstractDistMatrix<T>& A, S alpha, Int offset=0 );
template<typename T,typename S>
void ShiftDiagonal( AbstractBlockDistMatrix<T>& A, S alpha, Int offset=0 );
template<typename T,typename S>
void ShiftDiagonal( SparseMatrix<T>& A, S alpha, Int offset=0 );
template<typename T,typename S>
void ShiftDiagonal( DistSparseMatrix<T>& A, S alpha, Int offset=0 );

// Transpose
// =========
template<typename T>
void Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate=false );
template<typename T>
void Transpose
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B,
  bool conjugate=false );
template<typename T>
void Transpose
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B,
  bool conjugate=false );
template<typename T>
void Transpose
( const SparseMatrix<T>& A, SparseMatrix<T>& B, bool conjugate=false );
template<typename T>
void Transpose
( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B, bool conjugate=false );

namespace transpose {

template<typename T,Dist U,Dist V>
void ColFilter
( const DistMatrix<T,V,Collect<U>()>& A, 
        DistMatrix<T,U,V           >& B, bool conjugate=false );
template<typename T,Dist U,Dist V>
void ColFilter
( const BlockDistMatrix<T,V,Collect<U>()>& A, 
        BlockDistMatrix<T,U,V           >& B, bool conjugate=false );

template<typename T,Dist U,Dist V>
void RowFilter
( const DistMatrix<T,Collect<V>(),U>& A, 
        DistMatrix<T,        U,   V>& B, bool conjugate=false );
template<typename T,Dist U,Dist V>
void RowFilter
( const BlockDistMatrix<T,Collect<V>(),U>& A, 
        BlockDistMatrix<T,        U,   V>& B, bool conjugate=false );

template<typename T,Dist U,Dist V>
void PartialColFilter
( const DistMatrix<T,V,Partial<U>()>& A, 
        DistMatrix<T,U,        V   >& B, bool conjugate=false );
template<typename T,Dist U,Dist V>
void PartialColFilter
( const BlockDistMatrix<T,V,Partial<U>()>& A, 
        BlockDistMatrix<T,U,        V   >& B, bool conjugate=false );

template<typename T,Dist U,Dist V>
void PartialRowFilter
( const DistMatrix<T,Partial<V>(),U>& A, 
        DistMatrix<T,U,           V>& B, bool conjugate=false );
template<typename T,Dist U,Dist V>
void PartialRowFilter
( const BlockDistMatrix<T,Partial<V>(),U>& A, 
        BlockDistMatrix<T,U,           V>& B, bool conjugate=false );

template<typename T,Dist U,Dist V>
void ColSumScatter
( const DistMatrix<T,V,Collect<U>()>& A, 
        DistMatrix<T,U,V           >& B, bool conjugate=false );
template<typename T,Dist U,Dist V>
void ColSumScatter
( const BlockDistMatrix<T,V,Collect<U>()>& A, 
        BlockDistMatrix<T,U,V           >& B, bool conjugate=false );

template<typename T,Dist U,Dist V>
void PartialColSumScatter
( const DistMatrix<T,V,Partial<U>()>& A, 
        DistMatrix<T,U,V           >& B, bool conjugate=false );
template<typename T,Dist U,Dist V>
void PartialColSumScatter
( const BlockDistMatrix<T,V,Partial<U>()>& A, 
        BlockDistMatrix<T,U,V           >& B, bool conjugate=false );

template<typename T,Dist U,Dist V>
void ColAllGather
( const DistMatrix<T,U,        V   >& A, 
        DistMatrix<T,V,Collect<U>()>& B, bool conjugate=false );
template<typename T,Dist U,Dist V>
void ColAllGather
( const BlockDistMatrix<T,U,        V   >& A, 
        BlockDistMatrix<T,V,Collect<U>()>& B, bool conjugate=false );

template<typename T,Dist U,Dist V>
void PartialColAllGather
( const DistMatrix<T,U,        V   >& A, 
        DistMatrix<T,V,Partial<U>()>& B, bool conjugate=false );
template<typename T,Dist U,Dist V>
void PartialColAllGather
( const BlockDistMatrix<T,U,        V   >& A, 
        BlockDistMatrix<T,V,Partial<U>()>& B, bool conjugate=false );

} // namespace transpose

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
( S alpha, const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y, 
  bool conjugate=false );
template<typename T,typename S>
void TransposeAxpy
( S alpha, const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y, 
  bool conjugate=false );

namespace trans_axpy {

template<typename T,Dist U,Dist V>
void ColSumScatter
( T alpha, const DistMatrix<T,V,Collect<U>()>& A, 
                 DistMatrix<T,U,        V   >& B, bool conjugate=false );
template<typename T,Dist U,Dist V>
void ColSumScatter
( T alpha, const BlockDistMatrix<T,V,Collect<U>()>& A, 
                 BlockDistMatrix<T,U,        V   >& B, bool conjugate=false );

template<typename T,Dist U,Dist V>
void PartialColSumScatter
( T alpha, const DistMatrix<T,V,Partial<U>()>& A, 
                 DistMatrix<T,U,        V   >& B, bool conjugate=false );
template<typename T,Dist U,Dist V>
void PartialColSumScatter
( T alpha, const BlockDistMatrix<T,V,Partial<U>()>& A, 
                 BlockDistMatrix<T,U,        V   >& B, bool conjugate=false );

} // namespace trans_axpy

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
( DistMatrix<T,U,V>& A, T alpha, const AbstractDistMatrix<T>& d, 
  Int offset=0 );
template<typename T,Dist U,Dist V>
void UpdateRealPartOfDiagonal
( DistMatrix<T,U,V>& A, Base<T> alpha, const AbstractDistMatrix<Base<T>>& d, 
  Int offset=0 );
template<typename T,Dist U,Dist V>
void UpdateImagPartOfDiagonal
( DistMatrix<T,U,V>& A, Base<T> alpha, const AbstractDistMatrix<Base<T>>& d, 
  Int offset=0 );

// UpdateMappedDiagonal
// ====================
template<typename T,typename S>
void UpdateMappedDiagonal
( Matrix<T>& A, const Matrix<S>& d, 
  std::function<void(T&,S)> func, Int offset=0 );
template<typename T,typename S,Dist U,Dist V>
void UpdateMappedDiagonal
( DistMatrix<T,U,V>& A, const AbstractDistMatrix<S>& d, 
  std::function<void(T&,S)> func, Int offset=0 );

// UpdateSubmatrix
// ===============
template<typename T>
void UpdateSubmatrix
(       Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  T alpha, const Matrix<T>& ASub );
template<typename T>
void UpdateRealPartOfSubmatrix
(       Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  Base<T> alpha, const Matrix<Base<T>>& ASub );
template<typename T>
void UpdateImagPartOfSubmatrix
(       Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  Base<T> alpha, const Matrix<Base<T>>& ASub );

template<typename T>
void UpdateSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  T alpha, const AbstractDistMatrix<T>& ASub );
template<typename T>
void UpdateRealPartOfSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  Base<T> alpha, const AbstractDistMatrix<Base<T>>& ASub );
template<typename T>
void UpdateImagPartOfSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  Base<T> alpha, const AbstractDistMatrix<Base<T>>& ASub );

// Zero
// ====
template<typename T>
void Zero( Matrix<T>& A );
template<typename T>
void Zero( AbstractDistMatrix<T>& A );
template<typename T>
void Zero( AbstractBlockDistMatrix<T>& A );
template<typename T>
void Zero( SparseMatrix<T>& A );
template<typename T>
void Zero( DistSparseMatrix<T>& A );
template<typename T>
void Zero( DistMultiVec<T>& A );

} // namespace El

#endif // ifndef EL_BLAS1_HPP
