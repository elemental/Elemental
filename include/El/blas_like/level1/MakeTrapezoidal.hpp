/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_MAKETRAPEZOIDAL_HPP
#define EL_BLAS_MAKETRAPEZOIDAL_HPP

namespace El {

template<typename T>
void MakeTrapezoidal( UpperOrLower uplo, Matrix<T>& A, Int offset )
{
    EL_DEBUG_CSE
    const Int height = A.Height();
    const Int width = A.Width();
    const Int ldim = A.LDim();
    T* buffer = A.Buffer();

    if( uplo == LOWER )
    {
        EL_PARALLEL_FOR
        for( Int j=Max(0,offset+1); j<width; ++j )
        {
            const Int lastZeroRow = j-offset-1;
            const Int numZeroRows = Min( lastZeroRow+1, height );
            MemZero( &buffer[j*ldim], numZeroRows );
        }
    }
    else
    {
        EL_PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const Int firstZeroRow = Max(j-offset+1,0);
            if( firstZeroRow < height )
                MemZero( &buffer[firstZeroRow+j*ldim], height-firstZeroRow );
        }
    }
}

template<typename T>
void
MakeTrapezoidal( UpperOrLower uplo, AbstractDistMatrix<T>& A, Int offset )
{
    EL_DEBUG_CSE
    const Int height = A.Height();
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();

    T* buffer = A.Buffer();
    const Int ldim = A.LDim();

    if( uplo == LOWER )
    {
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int lastZeroRow = j-offset-1;
            if( lastZeroRow >= 0 )
            {
                const Int boundary = Min( lastZeroRow+1, height );
                const Int numZeroRows = A.LocalRowOffset(boundary);
                MemZero( &buffer[jLoc*ldim], numZeroRows );
            }
        }
    }
    else
    {
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int firstZeroRow = Max(j-offset+1,0);
            const Int numNonzeroRows = A.LocalRowOffset(firstZeroRow);
            if( numNonzeroRows < localHeight )
            {
                T* col = &buffer[numNonzeroRows+jLoc*ldim];
                MemZero( col, localHeight-numNonzeroRows );
            }
        }
    }
}

template<typename T>
void MakeTrapezoidal( UpperOrLower uplo, SparseMatrix<T>& A, Int offset )
{
    EL_DEBUG_CSE
    const Int numEntries = A.NumEntries();
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k);
        if( (uplo == LOWER && j-i > offset) ||
            (uplo == UPPER && j-i < offset) )
            A.QueueZero( i, j );
    }
    A.ProcessQueues();
}

template<typename T>
void MakeTrapezoidal( UpperOrLower uplo, DistSparseMatrix<T>& A, Int offset )
{
    EL_DEBUG_CSE
    const Int firstLocalRow = A.FirstLocalRow();
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k);
        if( (uplo == LOWER && j-i > offset) ||
            (uplo == UPPER && j-i < offset) )
            A.QueueLocalZero( i-firstLocalRow, j );
    }
    A.ProcessLocalQueues();
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void MakeTrapezoidal \
  ( UpperOrLower uplo, Matrix<T>& A, Int offset ); \
  EL_EXTERN template void MakeTrapezoidal \
  ( UpperOrLower uplo, AbstractDistMatrix<T>& A, Int offset ); \
  EL_EXTERN template void MakeTrapezoidal \
  ( UpperOrLower uplo, SparseMatrix<T>& A, Int offset ); \
  EL_EXTERN template void MakeTrapezoidal \
  ( UpperOrLower uplo, DistSparseMatrix<T>& A, Int offset );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_MAKETRAPEZOIDAL_HPP
