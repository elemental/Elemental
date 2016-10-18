/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_MAKESYMMETRIC_HPP
#define EL_BLAS_MAKESYMMETRIC_HPP

namespace El {

template<typename T>
void MakeSymmetric( UpperOrLower uplo, Matrix<T>& A, bool conjugate )
{
    DEBUG_CSE
    const Int n = A.Width();
    if( A.Height() != n )
        LogicError("Cannot make non-square matrix symmetric");

    if( conjugate )
        MakeDiagonalReal(A);

    T* ABuf = A.Buffer();
    const Int ldim = A.LDim();
    if( uplo == LOWER )
    {
        for( Int j=0; j<n; ++j )
        {
            for( Int i=j+1; i<n; ++i )
            {
                if( conjugate )
                    ABuf[j+i*ldim] = Conj(ABuf[i+j*ldim]); 
                else
                    ABuf[j+i*ldim] = ABuf[i+j*ldim];
            }    
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            for( Int i=0; i<j; ++i )
            {
                if( conjugate )
                    ABuf[j+i*ldim] = Conj(ABuf[i+j*ldim]);
                else
                    ABuf[j+i*ldim] = ABuf[i+j*ldim];
            }
        }
    }
}

template<typename T>
void MakeSymmetric
( UpperOrLower uplo, ElementalMatrix<T>& A, bool conjugate )
{
    DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    MakeTrapezoidal( uplo, A );
    if( conjugate )
        MakeDiagonalReal(A);

    unique_ptr<ElementalMatrix<T>> ATrans( A.Construct(A.Grid(),A.Root()) );
    Transpose( A, *ATrans, conjugate );
    if( uplo == LOWER )
        AxpyTrapezoid( UPPER, T(1), *ATrans, A, 1 );
    else
        AxpyTrapezoid( LOWER, T(1), *ATrans, A, -1 );
}

template<typename T>
void MakeSymmetric( UpperOrLower uplo, SparseMatrix<T>& A, bool conjugate )
{
    DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    MakeTrapezoidal( uplo, A );

    const Int m = A.Height();
    const Int numEntries = A.NumEntries();

    {
        const Int* sBuf = A.LockedSourceBuffer();
        const Int* tBuf = A.LockedTargetBuffer();
        T* vBuf = A.ValueBuffer();

        // Iterate over the diagonal entries
        Int numDiagonal = 0;
        for( Int i=0; i<m; ++i )
        {
            const Int e = A.Offset( i, i );
            if( e < numEntries && sBuf[e] == i && tBuf[e] == i )
            {
                ++numDiagonal;
                if( conjugate && IsComplex<T>::value )
                    vBuf[e] = RealPart(vBuf[e]);
            }
        }

        A.Reserve( numEntries-numDiagonal );
        // sBuf, tBuf, and vBuf are now invalidated due to reallocation
    }

    const Int* sBuf = A.LockedSourceBuffer();
    const Int* tBuf = A.LockedTargetBuffer();
    T* vBuf = A.ValueBuffer();

    for( Int k=0; k<numEntries; ++k ) 
    {
        if( sBuf[k] != tBuf[k] )
        {
            if( conjugate )
                A.QueueUpdate( tBuf[k], sBuf[k], Conj(vBuf[k]) );
            else
                A.QueueUpdate( tBuf[k], sBuf[k], vBuf[k] );
        }
    }
    A.ProcessQueues();
}

template<typename T>
void MakeSymmetric( UpperOrLower uplo, DistSparseMatrix<T>& A, bool conjugate )
{
    DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    MakeTrapezoidal( uplo, A );
    const Int numLocalEntries = A.NumLocalEntries();
    {
        T* vBuf = A.ValueBuffer();
        const Int* sBuf = A.LockedSourceBuffer();
        const Int* tBuf = A.LockedTargetBuffer();

        // Force the diagonal to be real
        // =============================
        if( conjugate && IsComplex<T>::value )
        {
            for( Int k=0; k<numLocalEntries; ++k )
                if( sBuf[k] == tBuf[k] )
                    vBuf[k] = RealPart(vBuf[k]);
        }

        // Compute the number of entries to send
        // =====================================
        Int numSend = 0;
        for( Int k=0; k<numLocalEntries; ++k )
        {
            const Int i = sBuf[k];
            const Int j = tBuf[k];
            if( i != j )
                ++numSend;
        }

        A.Reserve( numSend, numSend );
        // vBuf, sBuf, and tBuf are now invalidated due to reallocation
    }

    // Apply the updates
    // =================
    T* vBuf = A.ValueBuffer();
    const Int* sBuf = A.LockedSourceBuffer();
    const Int* tBuf = A.LockedTargetBuffer();
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = sBuf[k];
        const Int j = tBuf[k];
        if( i != j )
            A.QueueUpdate( j, i, ( conjugate ? Conj(vBuf[k]) : vBuf[k] ) );
    }
    A.ProcessQueues();
}

template<typename T>
void MakeHermitian( UpperOrLower uplo, Matrix<T>& A )
{
    DEBUG_CSE
    MakeSymmetric( uplo, A, true );
}

template<typename T>
void MakeHermitian( UpperOrLower uplo, ElementalMatrix<T>& A )
{
    DEBUG_CSE
    MakeSymmetric( uplo, A, true );
}

template<typename T>
void MakeHermitian( UpperOrLower uplo, SparseMatrix<T>& A )
{
    DEBUG_CSE
    MakeSymmetric( uplo, A, true );
}

template<typename T>
void MakeHermitian( UpperOrLower uplo, DistSparseMatrix<T>& A )
{
    DEBUG_CSE
    MakeSymmetric( uplo, A, true );
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void MakeSymmetric \
  ( UpperOrLower uplo, Matrix<T>& A, bool conjugate ); \
  EL_EXTERN template void MakeSymmetric \
  ( UpperOrLower uplo, ElementalMatrix<T>& A, bool conjugate ); \
  EL_EXTERN template void MakeSymmetric \
  ( UpperOrLower uplo, SparseMatrix<T>& A, bool conjugate ); \
  EL_EXTERN template void MakeSymmetric \
  ( UpperOrLower uplo, DistSparseMatrix<T>& A, bool conjugate ); \
  EL_EXTERN template void MakeHermitian \
  ( UpperOrLower uplo, Matrix<T>& A ); \
  EL_EXTERN template void MakeHermitian \
  ( UpperOrLower uplo, ElementalMatrix<T>& A ); \
  EL_EXTERN template void MakeHermitian \
  ( UpperOrLower uplo, SparseMatrix<T>& A ); \
  EL_EXTERN template void MakeHermitian \
  ( UpperOrLower uplo, DistSparseMatrix<T>& A );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_MAKESYMMETRIC_HPP
