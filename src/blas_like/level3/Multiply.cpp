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
( T alpha, const DistSparseMatrix<T>& A, const DistMultiVec<T>& X,
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
    const int YLocalHeight = Y.LocalHeight();
    const int width = X.Width();
    const int numLocalEntries = A.NumLocalEntries();

    // Y := beta Y
    for( int j=0; j<width; ++j )
        for( int iLocal=0; iLocal<YLocalHeight; ++iLocal )
            Y.SetLocal( iLocal, j, beta*Y.GetLocal(iLocal,j) );

    SparseMultMeta<T>& meta = A.multMeta;
    if( !meta.ready )
    {
        // Compute the set of row indices that we need from X
        std::set<int> indexSet;
        for( int e=0; e<numLocalEntries; ++e )
            indexSet.insert( A.Col(e) );
        const int numRecvInds = indexSet.size();
        std::vector<int> recvInds( numRecvInds );
        meta.recvSizes.clear();
        meta.recvSizes.resize( commSize, 0 );
        meta.recvOffs.resize( commSize );
        const int blocksize = A.Blocksize();
        {
            int off=0, lastOff=0, qPrev=0;
            std::set<int>::const_iterator setIt;
            for( setIt=indexSet.begin(); setIt!=indexSet.end(); ++setIt )
            {
                const int j = *setIt;
                const int q = RowToProcess( j, blocksize, commSize );
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
        int numSendInds=0;
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
        for( int s=0; s<numLocalEntries; ++s )
            meta.colOffs[s] = Find( recvInds, A.Col(s) );
        meta.numRecvInds = numRecvInds;
        meta.ready = true;
    }

    // Convert the sizes and offsets to be compatible with the current width
    std::vector<int> recvSizes=meta.recvSizes,
                     recvOffs=meta.recvOffs,
                     sendSizes=meta.sendSizes,
                     sendOffs=meta.sendOffs;
    for( int q=0; q<commSize; ++q )
    {
        recvSizes[q] *= width;    
        recvOffs[q] *= width;
        sendSizes[q] *= width;
        sendOffs[q] *= width;
    }

    // Pack the send values
    const int numSendInds = meta.sendInds.size();
    const int firstLocalRow = A.FirstLocalRow();
    std::vector<T> sendVals( numSendInds*width );
    for( int s=0; s<numSendInds; ++s )
    {
        const int i = meta.sendInds[s];
        const int iLocal = i - firstLocalRow;
        DEBUG_ONLY(
            if( iLocal < 0 || iLocal >= X.LocalHeight() )
                LogicError("iLocal was out of bounds");
        )
        for( int j=0; j<width; ++j )
            sendVals[s*width+j] = X.GetLocal( iLocal, j );
    }

    // Now send them
    std::vector<T> recvVals( meta.numRecvInds*width );
    mpi::AllToAll
    ( &sendVals[0], &sendSizes[0], &sendOffs[0],
      &recvVals[0], &recvSizes[0], &recvOffs[0], comm );
     
    // Perform the local multiply-accumulate, y := alpha A x + y
    for( int iLocal=0; iLocal<YLocalHeight; ++iLocal )
    {
        const int off = A.LocalEntryOffset( iLocal );
        const int rowSize = A.NumConnections( iLocal );
        for( int k=0; k<rowSize; ++k )
        {
            const int colOff = meta.colOffs[k+off];
            const T AVal = A.Value(k+off);
            for( int j=0; j<width; ++j )
            {
                const T XVal = recvVals[colOff*width+j];
                const T update = alpha*AVal*XVal;
                Y.UpdateLocal( iLocal, j, update );
            }
        }
    }
}

#define PROTO(T) \
    template void Multiply \
    ( T alpha, const DistSparseMatrix<T>& A, const DistMultiVec<T>& X, \
      T beta,                                      DistMultiVec<T>& Y );

#include "El/macros/Instantiate.h"

} // namespace El
