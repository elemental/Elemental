/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_HADAMARD_HPP
#define EL_BLAS_HADAMARD_HPP

// C(i,j) := A(i,j) B(i,j)

namespace El {

template<typename T>
void Hadamard( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C )
{
    EL_DEBUG_CSE
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Hadamard product requires equal dimensions");
    C.Resize( A.Height(), A.Width() );

    const Int height = A.Height();
    const Int width = A.Width();
    const T* ABuf = A.LockedBuffer();
    const T* BBuf = B.LockedBuffer();
    T* CBuf = C.Buffer();
    const Int ALDim = A.LDim();
    const Int BLDim = B.LDim();
    const Int CLDim = C.LDim();

    // Iterate over single loop if memory is contiguous. Otherwise
    // iterate over double loop.
    if( ALDim == height && BLDim == height && CLDim == height )
    {
        // Check if output matrix is equal to either input matrix
        if( CBuf == BBuf )
        {
            EL_PARALLEL_FOR
            for( Int i=0; i<height*width; ++i )
                CBuf[i] *= ABuf[i];
        }
        else if( CBuf == ABuf )
        {
            EL_PARALLEL_FOR
            for( Int i=0; i<height*width; ++i )
                CBuf[i] *= BBuf[i];
        }
        else
        {
            EL_PARALLEL_FOR
            for( Int i=0; i<height*width; ++i )
                CBuf[i] = ABuf[i] * BBuf[i];
        }
    }
    else
    {
        EL_PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            EL_SIMD
            for( Int i=0; i<height; ++i )
            {
                CBuf[i+j*CLDim] = ABuf[i+j*ALDim] * BBuf[i+j*BLDim];
            }
        }
    }
}

template<typename T>
void Hadamard
( const AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& C )
{
    EL_DEBUG_CSE
    const DistData& ADistData = A.DistData();
    const DistData& BDistData = B.DistData();
    DistData CDistData = C.DistData();
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Hadamard product requires equal dimensions");
    AssertSameGrids( A, B );
    if( ADistData.colDist != BDistData.colDist ||
        ADistData.rowDist != BDistData.rowDist ||
        BDistData.colDist != CDistData.colDist ||
        BDistData.rowDist != CDistData.rowDist )
        LogicError("A, B, and C must share the same distribution");
    if( A.ColAlign() != B.ColAlign() || A.RowAlign() != B.RowAlign() )
        LogicError("A and B must be aligned");
    if ( A.BlockHeight() != B.BlockHeight() ||
         A.BlockWidth() != B.BlockWidth())
      LogicError("A and B must have the same block size");
    C.AlignWith( A.DistData() );
    C.Resize( A.Height(), A.Width() );
    Hadamard( A.LockedMatrix(), B.LockedMatrix(), C.Matrix() );
}

template<typename T>
void Hadamard
( const DistMultiVec<T>& A, const DistMultiVec<T>& B, DistMultiVec<T>& C )
{
    EL_DEBUG_CSE
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Hadamard product requires equal dimensions");
    C.SetGrid( A.Grid() );
    C.Resize( A.Height(), A.Width() );
    Hadamard( A.LockedMatrix(), B.LockedMatrix(), C.Matrix() );
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void Hadamard \
  ( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C ); \
  EL_EXTERN template void Hadamard \
  ( const AbstractDistMatrix<T>& A, \
    const AbstractDistMatrix<T>& B, \
          AbstractDistMatrix<T>& C ); \
  EL_EXTERN template void Hadamard \
  ( const DistMultiVec<T>& A, \
    const DistMultiVec<T>& B, \
          DistMultiVec<T>& C );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_HADAMARD_HPP
