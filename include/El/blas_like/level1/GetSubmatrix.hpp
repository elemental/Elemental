/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_GETSUBMATRIX_HPP
#define EL_BLAS_GETSUBMATRIX_HPP

namespace El {

// Contiguous
// ==========
template<typename T>
void GetSubmatrix
( const Matrix<T>& A,
        Range<Int> I,
        Range<Int> J,
        Matrix<T>& ASub )
{
    EL_DEBUG_CSE
    auto ASubView = A(I,J);
    ASub = ASubView;
}

template<typename T>
void GetSubmatrix
( const ElementalMatrix<T>& A,
        Range<Int> I,
        Range<Int> J,
        ElementalMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    unique_ptr<ElementalMatrix<T>>
      ASubView( A.Construct(A.Grid(),A.Root()) );
    LockedView( *ASubView, A, I, J );
    Copy( *ASubView, ASub );
}

template<typename T>
void GetSubmatrix
( const SparseMatrix<T>& A,
        Range<Int> I,
        Range<Int> J,
        SparseMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    if( I.end == END ) I.end = A.Height();
    if( J.end == END ) J.end = A.Width();
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    ASub.Resize( mSub, nSub );
    Zero( ASub );

    const Int* offsetBuf = A.LockedOffsetBuffer();
    const Int* colBuf = A.LockedTargetBuffer();
    const T* valBuf = A.LockedValueBuffer();

    // Reserve the number of nonzeros that live within the submatrix
    Int numNonzerosSub = 0;
    for( Int i=I.beg; i<I.end; ++i )
    {
        const Int rowOff = offsetBuf[i];
        const Int numConn = offsetBuf[i+1] - offsetBuf[i];
        for( Int e=rowOff; e<rowOff+numConn; ++e )
        {
            const Int j = colBuf[e];
            if( j >= J.beg && j < J.end )
                ++numNonzerosSub;
        }
    }
    ASub.Reserve( numNonzerosSub );

    // Insert the nonzeros
    for( Int i=I.beg; i<I.end; ++i )
    {
        const Int rowOff = offsetBuf[i];
        const Int numConn = offsetBuf[i+1] - offsetBuf[i];
        for( Int e=rowOff; e<rowOff+numConn; ++e )
        {
            const Int j = colBuf[e];
            if( j >= J.beg && j < J.end )
                ASub.QueueUpdate( i-I.beg, j-J.beg, valBuf[e] );
        }
    }
    ASub.ProcessQueues();
}

template<typename T>
void GetSubmatrix
( const SparseMatrix<T>& A,
        Range<Int> I,
  const vector<Int>& J,
        SparseMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    // TODO(poulson): Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

template<typename T>
void GetSubmatrix
( const SparseMatrix<T>& A,
  const vector<Int>& I,
        Range<Int> J,
        SparseMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    // TODO(poulson): Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

template<typename T>
void GetSubmatrix
( const SparseMatrix<T>& A,
  const vector<Int>& I,
  const vector<Int>& J,
        SparseMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    // TODO(poulson): Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

// TODO(poulson): Use lower-level access
template<typename T>
void GetSubmatrix
( const DistSparseMatrix<T>& A,
        Range<Int> I,
        Range<Int> J,
        DistSparseMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    if( I.end == END ) I.end = A.Height();
    if( J.end == END ) J.end = A.Width();

    ASub.SetGrid( A.Grid() );
    ASub.Resize( I.end-I.beg, J.end-J.beg );
    Zero( ASub );

    // Count the number of updates
    // ===========================
    Int numUpdates = 0;
    for( Int e=0; e<A.NumLocalEntries(); ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
            ++numUpdates;
    }

    // Queue and process the updates
    // =============================
    ASub.Reserve( numUpdates, numUpdates );
    for( Int e=0; e<A.NumLocalEntries(); ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
            ASub.QueueUpdate( i-I.beg, j-J.beg, A.Value(e) );
    }
    ASub.ProcessQueues();
}

template<typename T>
void GetSubmatrix
( const DistSparseMatrix<T>& A,
        Range<Int> I,
  const vector<Int>& J,
        DistSparseMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    // TODO(poulson): Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

template<typename T>
void GetSubmatrix
( const DistSparseMatrix<T>& A,
  const vector<Int>& I,
        Range<Int> J,
        DistSparseMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    // TODO(poulson): Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

template<typename T>
void GetSubmatrix
( const DistSparseMatrix<T>& A,
  const vector<Int>& I,
  const vector<Int>& J,
        DistSparseMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    // TODO(poulson): Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

template<typename T>
void GetSubmatrix
( const DistMultiVec<T>& A,
        Range<Int> I,
        Range<Int> J,
        DistMultiVec<T>& ASub )
{
    EL_DEBUG_CSE
    if( I.end == END ) I.end = A.Height();
    if( J.end == END ) J.end = A.Width();
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    const Int localHeight = A.LocalHeight();

    ASub.SetGrid( A.Grid() );
    ASub.Resize( mSub, nSub );
    Zero( ASub );

    const T* ABuf = A.LockedMatrix().LockedBuffer();
    const Int ALDim = A.LockedMatrix().LDim();

    // If no communication is necessary, take the easy and fast approach
    if( mSub == A.Height() )
    {
        Copy( A.LockedMatrix(), ASub.Matrix() );
        return;
    }

    // Compute the number of updates
    // =============================
    Int numUpdates = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        if( i >= I.end )
            break;
        else if( i >= I.beg )
            numUpdates += nSub;
    }

    // Queue and process the updates
    // =============================
    ASub.Reserve( numUpdates );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        if( i >= I.end )
            break;
        else if( i >= I.beg )
        {
            for( Int j=J.beg; j<J.end; ++j )
                ASub.QueueUpdate( i-I.beg, j-J.beg, ABuf[iLoc+j*ALDim] );
        }
    }
    ASub.ProcessQueues();
}

// Non-contiguous
// ==============
template<typename T>
void GetSubmatrix
( const Matrix<T>& A,
  const Range<Int> I,
  const vector<Int>& J,
        Matrix<T>& ASub )
{
    EL_DEBUG_CSE
    const Int mSub = I.end-I.beg;
    const Int nSub = J.size();
    ASub.Resize( mSub, nSub );

    T* ASubBuf = ASub.Buffer();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    const Int ASubLDim = ASub.LDim();

    for( Int jSub=0; jSub<nSub; ++jSub )
    {
        const Int j = J[jSub];
        MemCopy( &ASubBuf[jSub*ASubLDim], &ABuf[j*ALDim], mSub );
    }
}

template<typename T>
void GetSubmatrix
( const Matrix<T>& A,
  const vector<Int>& I,
  const Range<Int> J,
        Matrix<T>& ASub )
{
    EL_DEBUG_CSE
    const Int mSub = I.size();
    const Int nSub = J.end-J.beg;
    ASub.Resize( mSub, nSub );

    T* ASubBuf = ASub.Buffer();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    const Int ASubLDim = ASub.LDim();

    for( Int jSub=0; jSub<nSub; ++jSub )
    {
        const Int j = J.beg + jSub;
        for( Int iSub=0; iSub<mSub; ++iSub )
        {
            const Int i = I[iSub];
            ASubBuf[iSub+jSub*ASubLDim] = ABuf[i+j*ALDim];
        }
    }
}

template<typename T>
void GetSubmatrix
( const Matrix<T>& A,
  const vector<Int>& I,
  const vector<Int>& J,
        Matrix<T>& ASub )
{
    EL_DEBUG_CSE
    const Int mSub = I.size();
    const Int nSub = J.size();
    ASub.Resize( mSub, nSub );

    T* ASubBuf = ASub.Buffer();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    const Int ASubLDim = ASub.LDim();

    for( Int jSub=0; jSub<nSub; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<mSub; ++iSub )
        {
            const Int i = I[iSub];
            ASubBuf[iSub+jSub*ASubLDim] = ABuf[i+j*ALDim];
        }
    }
}

template<typename T>
void GetSubmatrix
( const AbstractDistMatrix<T>& A,
        Range<Int> I,
  const vector<Int>& J,
        AbstractDistMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    const Int mSub = I.end-I.beg;
    const Int nSub = J.size();
    const Grid& g = A.Grid();
    ASub.SetGrid( g );
    ASub.Resize( mSub, nSub );
    Zero( ASub );

    // TODO(poulson): Intelligently pick the redundant rank to pack from?

    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    // Count the number of updates
    // ===========================
    Int numUpdates = 0;
    if( A.RedundantRank() == 0 )
        for( Int i=I.beg; i<I.end; ++i )
            if( A.IsLocalRow(i) )
                for( auto& j : J )
                    if( A.IsLocalCol(j) )
                        ++numUpdates;

    // Queue and process the updates
    // =============================
    ASub.Reserve( numUpdates );
    if( A.RedundantRank() == 0 )
    {
        for( Int iSub=0; iSub<mSub; ++iSub )
        {
            const Int i = I.beg + iSub;
            if( A.IsLocalRow(i) )
            {
                const Int iLoc = A.LocalRow(i);
                for( Int jSub=0; jSub<nSub; ++jSub )
                {
                    const Int j = J[jSub];
                    if( A.IsLocalCol(j) )
                    {
                        const Int jLoc = A.LocalCol(j);
                        ASub.QueueUpdate( iSub, jSub, ABuf[iLoc+jLoc*ALDim] );
                    }
                }
            }
        }
    }
    ASub.ProcessQueues();
}

template<typename T>
void GetSubmatrix
( const AbstractDistMatrix<T>& A,
  const vector<Int>& I,
        Range<Int> J,
        AbstractDistMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    const Int mSub = I.size();
    const Int nSub = J.end-J.beg;
    const Grid& g = A.Grid();
    ASub.SetGrid( g );
    ASub.Resize( mSub, nSub );
    Zero( ASub );

    // TODO(poulson): Intelligently pick the redundant rank to pack from?

    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    // Count the number of updates
    // ===========================
    Int numUpdates = 0;
    if( A.RedundantRank() == 0 )
        for( auto& i : I )
            if( A.IsLocalRow(i) )
                for( Int j=J.beg; j<J.end; ++j )
                    if( A.IsLocalCol(j) )
                        ++numUpdates;

    // Queue and process the updates
    // =============================
    ASub.Reserve( numUpdates );
    if( A.RedundantRank() == 0 )
    {
        for( Int iSub=0; iSub<mSub; ++iSub )
        {
            const Int i = I[iSub];
            if( A.IsLocalRow(i) )
            {
                const Int iLoc = A.LocalRow(i);
                for( Int jSub=0; jSub<nSub; ++jSub )
                {
                    const Int j = J.beg + jSub;
                    if( A.IsLocalCol(j) )
                    {
                        const Int jLoc = A.LocalCol(j);
                        ASub.QueueUpdate( iSub, jSub, ABuf[iLoc+jLoc*ALDim] );
                    }
                }
            }
        }
    }
    ASub.ProcessQueues();
}

template<typename T>
void GetSubmatrix
( const AbstractDistMatrix<T>& A,
  const vector<Int>& I,
  const vector<Int>& J,
        AbstractDistMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    const Int mSub = I.size();
    const Int nSub = J.size();
    const Grid& g = A.Grid();
    ASub.SetGrid( g );
    ASub.Resize( mSub, nSub );
    Zero( ASub );

    // TODO(poulson): Intelligently pick the redundant rank to pack from?

    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    // Count the number of updates
    // ===========================
    Int numUpdates = 0;
    if( A.RedundantRank() == 0 )
        for( auto& i : I )
            if( A.IsLocalRow(i) )
                for( auto& j : J )
                    if( A.IsLocalCol(j) )
                        ++numUpdates;

    // Queue and process the updates
    // =============================
    ASub.Reserve( numUpdates );
    if( A.RedundantRank() == 0 )
    {
        for( Int iSub=0; iSub<mSub; ++iSub )
        {
            const Int i = I[iSub];
            if( A.IsLocalRow(i) )
            {
                const Int iLoc = A.LocalRow(i);
                for( Int jSub=0; jSub<nSub; ++jSub )
                {
                    const Int j = J[jSub];
                    if( A.IsLocalCol(j) )
                    {
                        const Int jLoc = A.LocalCol(j);
                        ASub.QueueUpdate( iSub, jSub, ABuf[iLoc+jLoc*ALDim] );
                    }
                }
            }
        }
    }
    ASub.ProcessQueues();
}

template<typename T>
void GetSubmatrix
( const DistMultiVec<T>& A,
        Range<Int> I,
  const vector<Int>& J,
        DistMultiVec<T>& ASub )
{
    EL_DEBUG_CSE
    // TODO(poulson): Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

template<typename T>
void GetSubmatrix
( const DistMultiVec<T>& A,
  const vector<Int>& I,
        Range<Int> J,
        DistMultiVec<T>& ASub )
{
    EL_DEBUG_CSE
    // TODO(poulson): Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

template<typename T>
void GetSubmatrix
( const DistMultiVec<T>& A,
  const vector<Int>& I,
  const vector<Int>& J,
        DistMultiVec<T>& ASub )
{
    EL_DEBUG_CSE
    // TODO(poulson): Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void GetSubmatrix \
  ( const Matrix<T>& A, \
          Range<Int> I, \
          Range<Int> J,  \
          Matrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const ElementalMatrix<T>& A, \
          Range<Int> I, \
          Range<Int> J, \
          ElementalMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const SparseMatrix<T>& A, \
          Range<Int> I, \
          Range<Int> J, \
          SparseMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const SparseMatrix<T>& A, \
          Range<Int> I, \
    const vector<Int>& J, \
          SparseMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const SparseMatrix<T>& A, \
    const vector<Int>& I, \
          Range<Int> J, \
          SparseMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const SparseMatrix<T>& A, \
    const vector<Int>& I, \
    const vector<Int>& J, \
          SparseMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const DistSparseMatrix<T>& A, \
          Range<Int> I, \
          Range<Int> J, \
          DistSparseMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const DistSparseMatrix<T>& A, \
          Range<Int> I, \
    const vector<Int>& J, \
          DistSparseMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const DistSparseMatrix<T>& A, \
    const vector<Int>& I, \
          Range<Int> J, \
          DistSparseMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const DistSparseMatrix<T>& A, \
    const vector<Int>& I, \
    const vector<Int>& J, \
          DistSparseMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const DistMultiVec<T>& A, \
          Range<Int> I, \
          Range<Int> J, \
          DistMultiVec<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const Matrix<T>& A, \
    const Range<Int> I, \
    const vector<Int>& J,  \
          Matrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const Matrix<T>& A, \
    const vector<Int>& I, \
    const Range<Int> J, \
          Matrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const Matrix<T>& A, \
    const vector<Int>& I, \
    const vector<Int>& J, \
          Matrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
          Range<Int> I, \
    const vector<Int>& J, \
          AbstractDistMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const vector<Int>& I, \
          Range<Int> J, \
          AbstractDistMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const vector<Int>& I, \
    const vector<Int>& J, \
          AbstractDistMatrix<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const DistMultiVec<T>& A, \
          Range<Int> I, \
    const vector<Int>& J, \
          DistMultiVec<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const DistMultiVec<T>& A, \
    const vector<Int>& I, \
          Range<Int> J, \
          DistMultiVec<T>& ASub ); \
  EL_EXTERN template void GetSubmatrix \
  ( const DistMultiVec<T>& A, \
    const vector<Int>& I, \
    const vector<Int>& J, \
          DistMultiVec<T>& ASub );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_GETSUBMATRIX_HPP
