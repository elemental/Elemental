/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_KRONECKER_HPP
#define EL_BLAS_KRONECKER_HPP

namespace El {

template<typename T> 
void Kronecker
( const Matrix<T>& A,
  const Matrix<T>& B,
        Matrix<T>& C )
{
    DEBUG_ONLY(CSE cse("Kronecker"))
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int mB = B.Height();
    const Int nB = B.Width();
    C.Resize( mA*mB, nA*nB );

    for( Int jA=0; jA<nA; ++jA )
    {
        for( Int iA=0; iA<mA; ++iA )
        {
            auto Cij = C( IR(iA*mB,(iA+1)*mB), IR(jA*nB,(jA+1)*nB) );
            Cij = B;
            Scale( A.Get(iA,jA), Cij );
        }
    }
}

template<typename T> 
void Kronecker
( const Matrix<T>& A,
  const Matrix<T>& B,
        ElementalMatrix<T>& CPre )
{
    DEBUG_ONLY(CSE cse("Kronecker"))

    DistMatrixWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& C = CProx.Get();

    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int mB = B.Height();
    const Int nB = B.Width();
    C.Resize( mA*mB, nA*nB );

    const Int localHeight = C.LocalHeight();
    const Int localWidth = C.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = C.GlobalCol(jLoc);
        const Int jA = j / nB;
        const Int jB = j % nB;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = C.GlobalRow(iLoc);
            const Int iA = i / mB;
            const Int iB = i % mB;
            C.SetLocal( iLoc, jLoc, A.Get(iA,jA)*B.Get(iB,jB) );
        }
    }
}

template<typename T> 
void Kronecker
( const SparseMatrix<T>& A,
  const SparseMatrix<T>& B,
        SparseMatrix<T>& C )
{
    DEBUG_ONLY(CSE cse("Kronecker"))
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int mB = B.Height();
    const Int nB = B.Width();
    const Int numEntriesA = A.NumEntries();
    const Int numEntriesB = B.NumEntries();
    C.Resize( mA*mB, nA*nB );
    C.Reserve( numEntriesA*numEntriesB );

    for( Int eA=0; eA<numEntriesA; ++eA )
    {
        const Int iA = A.Row(eA);
        const Int jA = A.Col(eA);
        const T valA = A.Value(eA);
        for( Int eB=0; eB<numEntriesB; ++eB )
        {
            const Int iB = B.Row(eB);
            const Int jB = B.Col(eB);
            const T valB = B.Value(eB);
            
            const Int i = iA*mB + iB;
            const Int j = jA*nB + jB;
            C.QueueUpdate( i, j, valA*valB );
        }
    }
    C.ProcessQueues();
}

template<typename T> 
void Kronecker
( const SparseMatrix<T>& A, const Matrix<T>& B, SparseMatrix<T>& C )
{
    DEBUG_ONLY(CSE cse("Kronecker"))
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int mB = B.Height();
    const Int nB = B.Width();
    const Int numEntriesA = A.NumEntries();
    C.Resize( mA*mB, nA*nB );
    C.Reserve( numEntriesA*mB*nB );

    for( Int eA=0; eA<numEntriesA; ++eA )
    {
        const Int iA = A.Row(eA);
        const Int jA = A.Col(eA);
        const T valA = A.Value(eA);

        for( Int jB=0; jB<nB; ++jB )
        {
            const Int j = jA*nB + jB;
            for( Int iB=0; iB<mB; ++iB )
            {
                const T valB = B.Get(iB,jB);
                const Int i = iA*mB + iB;
                C.QueueUpdate( i, j, valA*valB );
            }
        }
    }
    C.ProcessQueues();
}

template<typename T> 
void Kronecker
( const Matrix<T>& A, const SparseMatrix<T>& B, SparseMatrix<T>& C )
{
    DEBUG_ONLY(CSE cse("Kronecker"))
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int mB = B.Height();
    const Int nB = B.Width();
    const Int numEntriesB = B.NumEntries();
    C.Resize( mA*mB, nA*nB );
    C.Reserve( mA*nA*numEntriesB );

    for( Int jA=0; jA<nA; ++jA )
    {
        for( Int iA=0; iA<mA; ++iA )
        {
            const T valA = A.Get(iA,jA);
            for( Int eB=0; eB<numEntriesB; ++eB )
            {
                const Int iB = B.Row(eB);
                const Int jB = B.Col(eB);
                const T valB = B.Value(eB);
            
                const Int i = iA*mB + iB;
                const Int j = jA*nB + jB;
                C.QueueUpdate( i, j, valA*valB );
            }
        }
    }
    C.ProcessQueues();
}

template<typename T> 
void Kronecker
( const SparseMatrix<T>& A, const SparseMatrix<T>& B, DistSparseMatrix<T>& C )
{
    DEBUG_ONLY(CSE cse("Kronecker"))
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int mB = B.Height();
    const Int nB = B.Width();
    const Int numEntriesA = A.NumEntries();
    const Int numEntriesB = B.NumEntries();
    C.Resize( mA*mB, nA*nB );

    // Count the number of nonzeros in our local portion of C
    Int numEntriesC = 0;
    const Int localHeight = C.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = C.GlobalRow(iLoc);
        const Int iA = i / mB;
        const Int iB = i % mB;

        const Int numConnectA = A.NumConnections( iA );
        const Int numConnectB = B.NumConnections( iB );
        numEntriesC += numConnectA*numConnectB;
    }
    C.Reserve( numEntriesA*numEntriesB );

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = C.GlobalRow(iLoc);
        const Int iA = i / mB;
        const Int iB = i % mB;

        const Int offA = A.RowOffset( iA );
        const Int offB = A.RowOffset( iB );
        const Int numConnectA = A.NumConnections( iA );
        const Int numConnectB = B.NumConnections( iB );
        for( Int eA=offA; eA<offA+numConnectA; ++eA )
        {
            const Int jA = A.Col( eA );
            const T valA = A.Value( eA );
            for( Int eB=offB; eB<offB+numConnectB; ++eB )
            {
                const Int jB = B.Col( eB );
                const T valB = B.Value( eB );

                const Int j = jA*nB + jB;
                C.QueueLocalUpdate( iLoc, j, valA*valB );
            }
        }
    }
    C.ProcessLocalQueues();
}

template<typename T> 
void Kronecker
( const SparseMatrix<T>& A, const Matrix<T>& B, DistSparseMatrix<T>& C )
{
    DEBUG_ONLY(CSE cse("Kronecker"))
    LogicError("This routine is not yet written");
}

template<typename T> 
void Kronecker
( const Matrix<T>& A, const SparseMatrix<T>& B, DistSparseMatrix<T>& C )
{
    DEBUG_ONLY(CSE cse("Kronecker"))
    LogicError("This routine is not yet written");
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void Kronecker \
  ( const Matrix<T>& A, \
    const Matrix<T>& B, \
          Matrix<T>& C ); \
  EL_EXTERN template void Kronecker \
  ( const Matrix<T>& A, \
    const Matrix<T>& B, \
          ElementalMatrix<T>& C ); \
  EL_EXTERN template void Kronecker \
  ( const SparseMatrix<T>& A, \
    const SparseMatrix<T>& B, \
          SparseMatrix<T>& C ); \
  EL_EXTERN template void Kronecker \
  ( const Matrix<T>& A, \
    const SparseMatrix<T>& B, \
          SparseMatrix<T>& C ); \
  EL_EXTERN template void Kronecker \
  ( const SparseMatrix<T>& A, \
    const Matrix<T>& B, \
          SparseMatrix<T>& C ); \
  EL_EXTERN template void Kronecker \
  ( const SparseMatrix<T>& A, \
    const SparseMatrix<T>& B, \
          DistSparseMatrix<T>& C ); \
  EL_EXTERN template void Kronecker \
  ( const Matrix<T>& A, \
    const SparseMatrix<T>& B, \
          DistSparseMatrix<T>& C ); \
  EL_EXTERN template void Kronecker \
  ( const SparseMatrix<T>& A, \
    const Matrix<T>& B, \
          DistSparseMatrix<T>& C );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_KRONECKER_HPP
