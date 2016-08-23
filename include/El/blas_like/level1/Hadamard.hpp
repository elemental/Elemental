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
    DEBUG_CSE
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Hadamard product requires equal dimensions");
    C.Resize( A.Height(), A.Width() );

    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            C(i,j) = A(i,j)*B(i,j);
}

template<typename T> 
void Hadamard
( const ElementalMatrix<T>& A,
  const ElementalMatrix<T>& B, 
        ElementalMatrix<T>& C )
{
    DEBUG_CSE
    const ElementalData ADistData = A.DistData();
    const ElementalData BDistData = B.DistData();
    ElementalData CDistData = C.DistData();
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
    C.AlignWith( A.DistData() );
    C.Resize( A.Height(), A.Width() );
    Hadamard( A.LockedMatrix(), B.LockedMatrix(), C.Matrix() );
}

template<typename T> 
void Hadamard
( const DistMultiVec<T>& A, const DistMultiVec<T>& B, DistMultiVec<T>& C )
{
    DEBUG_CSE
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Hadamard product requires equal dimensions");
    C.SetComm( A.Comm() );
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
  ( const ElementalMatrix<T>& A, \
    const ElementalMatrix<T>& B, \
          ElementalMatrix<T>& C ); \
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
