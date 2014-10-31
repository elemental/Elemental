/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Multiply
( Orientation orientation, 
  T alpha, const SparseMatrix<T>& A, const Matrix<T>& X,
  T beta,                                  Matrix<T>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("Multiply");
        if( A.Height() != Y.Height() || A.Width() != X.Height() || 
            X.Width() != Y.Width() )
            LogicError("A, X, and Y did not conform");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int b = X.Width();

    // Y := beta Y
    Scale( beta, Y );

    // Accumulate
    if( orientation == NORMAL )
    {
        // Y := alpha A X + Y
        for( Int i=0; i<m; ++i )
        {
            const Int off = A.EntryOffset( i );
            const Int rowSize = A.NumConnections( i );
            for( Int k=0; k<rowSize; ++k )
            {
                const Int j = A.Col(k+off);
                const T AVal = A.Value(k+off);
                for( Int t=0; t<b; ++t )
                {
                    const T XVal = X.Get(j,t);
                    Y.Update( i, t, alpha*AVal*XVal );
                }
            }
        }
    }
    else
    {
        // Y := alpha A' X + Y
        const bool conjugate = ( orientation == ADJOINT );
        for( Int i=0; i<m; ++i ) 
        {
            const Int off = A.EntryOffset( i );
            const Int rowSize = A.NumConnections( i );
            for( Int k=0; k<rowSize; ++k )
            {
                const Int j = A.Col(k+off); 
                const T AVal = A.Value(k+off);
                for( Int t=0; t<b; ++t )
                {
                    const T XVal = X.Get(i,t);
                    if( conjugate )
                        Y.Update( j, t, alpha*Conj(AVal)*XVal );
                    else
                        Y.Update( j, t, alpha*AVal*XVal );
                }
            }
        }
    }
}

template<typename T>
void Multiply
( Orientation orientation, 
  T alpha, const DistSparseMatrix<T>& A, const DistMultiVec<T>& X,
  T beta,                                      DistMultiVec<T>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("Multiply");
        if( A.Height() != Y.Height() || A.Width() != X.Height() || 
            X.Width() != Y.Width() )
            LogicError("A, X, and Y did not conform");
        if( !mpi::Congruent( A.Comm(), X.Comm() ) || 
            !mpi::Congruent( X.Comm(), Y.Comm() ) )
            LogicError("Communicators did not match");
    )
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );

    // Y := beta Y
    Scale( beta, Y );

    if( orientation == NORMAL )
    {
        SparseNormalMultMeta<T>& meta = A.normalMultMeta;
        if( !meta.ready )
        {
            // Compute the set of row indices that we need from X
            const Int numLocalEntries = A.NumLocalEntries();
            std::set<Int> indexSet;
            for( Int e=0; e<numLocalEntries; ++e )
                indexSet.insert( A.Col(e) );
            const Int numRecvInds = indexSet.size();
            std::vector<Int> recvInds( numRecvInds );
            meta.recvSizes.clear();
            meta.recvSizes.resize( commSize, 0 );
            meta.recvOffs.resize( commSize );
            const Int blocksize = X.Blocksize();
            {
                Int off=0, lastOff=0, qPrev=0;
                std::set<Int>::const_iterator setIt;
                for( setIt=indexSet.begin(); setIt!=indexSet.end(); ++setIt )
                {
                    const Int j = *setIt;
                    const Int q = RowToProcess( j, blocksize, commSize );
                    while( qPrev != q )
                    {
                        meta.recvSizes[qPrev] = off - lastOff;
                        meta.recvOffs[qPrev+1] = off;
    
                        lastOff = off;
                        ++qPrev;
                    }
                    recvInds[off++] = j;
                }
                while( qPrev != commSize-1 )
                {
                    meta.recvSizes[qPrev] = off - lastOff;
                    meta.recvOffs[qPrev+1] = off;
                    lastOff = off;
                    ++qPrev;
                }
                meta.recvSizes[commSize-1] = off - lastOff;
            }

            // Coordinate
            meta.sendSizes.resize( commSize );
            mpi::AllToAll( &meta.recvSizes[0], 1, &meta.sendSizes[0], 1, comm );
            Int numSendInds=0;
            meta.sendOffs.resize( commSize );
            for( int q=0; q<commSize; ++q )
            {
                meta.sendOffs[q] = numSendInds;
                numSendInds += meta.sendSizes[q];
            }
            meta.sendInds.resize( numSendInds );
            mpi::AllToAll
            ( &recvInds[0],       &meta.recvSizes[0], &meta.recvOffs[0],
              &meta.sendInds[0], &meta.sendSizes[0], &meta.sendOffs[0], comm );

            meta.colOffs.resize( numLocalEntries );
            for( Int s=0; s<numLocalEntries; ++s )
                meta.colOffs[s] = Find( recvInds, A.Col(s) );
            meta.numRecvInds = numRecvInds;
            meta.ready = true;
        }

        // Convert the sizes and offsets to be compatible with the current width
        const Int b = X.Width();
        std::vector<int> recvSizes=meta.recvSizes,
                         recvOffs=meta.recvOffs,
                         sendSizes=meta.sendSizes,
                         sendOffs=meta.sendOffs;
        for( int q=0; q<commSize; ++q )
        {
            recvSizes[q] *= b;    
            recvOffs[q] *= b;
            sendSizes[q] *= b;
            sendOffs[q] *= b;
        }

        // Pack the send values
        const Int numSendInds = meta.sendInds.size();
        const Int firstLocalRow = X.FirstLocalRow();
        std::vector<T> sendVals( numSendInds*b );
        for( Int s=0; s<numSendInds; ++s )
        {
            const Int i = meta.sendInds[s];
            const Int iLocal = i - firstLocalRow;
            DEBUG_ONLY(
                if( iLocal < 0 || iLocal >= X.LocalHeight() )
                    LogicError("iLocal was out of bounds");
            )
            for( Int t=0; t<b; ++t )
                sendVals[s*b+t] = X.GetLocal( iLocal, t );
        }

        // Now send them
        std::vector<T> recvVals( meta.numRecvInds*b );
        mpi::AllToAll
        ( &sendVals[0], &sendSizes[0], &sendOffs[0],
          &recvVals[0], &recvSizes[0], &recvOffs[0], comm );
     
        // Perform the local multiply-accumulate, y := alpha A x + y
        const Int YLocalHeight = Y.LocalHeight();
        for( Int iLocal=0; iLocal<YLocalHeight; ++iLocal )
        {
            const Int off = A.EntryOffset( iLocal );
            const Int rowSize = A.NumConnections( iLocal );
            for( Int k=0; k<rowSize; ++k )
            {
                const Int colOff = meta.colOffs[k+off];
                const T AVal = A.Value(k+off);
                for( int t=0; t<b; ++t )
                {
                    const T XVal = recvVals[colOff*b+t];
                    Y.UpdateLocal( iLocal, t, alpha*AVal*XVal );
                }
            }
        }
    }
    else
    {
        LogicError("Distributed adjoint sparse multiply not yet supported");
    }
}

#define PROTO(T) \
    template void Multiply \
    ( Orientation orientation, \
      T alpha, const SparseMatrix<T>& A, const Matrix<T>& X, \
      T beta,                                  Matrix<T>& Y ); \
    template void Multiply \
    ( Orientation orientation, \
      T alpha, const DistSparseMatrix<T>& A, const DistMultiVec<T>& X, \
      T beta,                                      DistMultiVec<T>& Y );

#include "El/macros/Instantiate.h"

} // namespace El
