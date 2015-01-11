/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
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
        if( X.Width() != Y.Width() )
            LogicError("X and Y must have the same width");
    )
    const Int m = A.Height();
    const Int b = X.Width();

    // Y := beta Y
    Scale( beta, Y );

    // Accumulate
    if( orientation == NORMAL )
    {
        // Y := alpha A X + Y
        if( A.Height() != Y.Height() )
            LogicError("A and Y must have the same height");
        if( A.Width() != X.Height() )
            LogicError("The width of A must match the height of X");
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
        if( A.Width() != Y.Height() )
            LogicError("The width of A must match the height of Y");
        if( A.Height() != X.Height() )
            LogicError("The height of A must match the height of X");
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
        if( X.Width() != Y.Width() )
            LogicError("X and Y must have the same width");
        if( !mpi::Congruent( A.Comm(), X.Comm() ) || 
            !mpi::Congruent( X.Comm(), Y.Comm() ) )
            LogicError("Communicators did not match");
    )
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );

    // Y := beta Y
    Scale( beta, Y );

    SparseMultMeta<T>& meta = A.multMeta;
    if( !meta.ready )
    {
        // Compute the set of row indices that we need from X in a normal
        // multiply or update of Y in the adjoint case
        const Int numLocalEntries = A.NumLocalEntries();
        std::set<Int> indexSet;
        for( Int e=0; e<numLocalEntries; ++e )
            indexSet.insert( A.Col(e) );
        const Int numRecvInds = indexSet.size();
        std::vector<Int> recvInds( numRecvInds );
        meta.recvSizes.clear();
        meta.recvSizes.resize( commSize, 0 );
        meta.recvOffs.resize( commSize );
        const Int blocksize = 
          ( orientation == NORMAL ? X.Blocksize() : Y.Blocksize() );;
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
        mpi::AllToAll
        ( meta.recvSizes.data(), 1, meta.sendSizes.data(), 1, comm );
        Int numSendInds=0;
        meta.sendOffs.resize( commSize );
        for( int q=0; q<commSize; ++q )
        {
            meta.sendOffs[q] = numSendInds;
            numSendInds += meta.sendSizes[q];
        }
        meta.sendInds.resize( numSendInds );
        mpi::AllToAll
        ( recvInds.data(),      meta.recvSizes.data(), meta.recvOffs.data(),
          meta.sendInds.data(), meta.sendSizes.data(), meta.sendOffs.data(), 
          comm );

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

    if( orientation == NORMAL )
    {
        if( A.Height() != Y.Height() )
            LogicError("A and Y must have the same height");
        if( A.Width() != X.Height() )
            LogicError("The width of A must match the height of X");

        // Pack the send values
        const Int numSendInds = meta.sendInds.size();
        const Int firstLocalRow = X.FirstLocalRow();
        std::vector<T> sendVals( numSendInds*b );
        for( Int s=0; s<numSendInds; ++s )
        {
            const Int i = meta.sendInds[s];
            const Int iLoc = i - firstLocalRow;
            for( Int t=0; t<b; ++t )
                sendVals[s*b+t] = X.GetLocal( iLoc, t );
        }

        // Now send them
        std::vector<T> recvVals( meta.numRecvInds*b );
        mpi::AllToAll
        ( sendVals.data(), sendSizes.data(), sendOffs.data(),
          recvVals.data(), recvSizes.data(), recvOffs.data(), comm );
     
        // Perform the local multiply-accumulate, y := alpha A x + y
        const Int ALocalHeight = A.LocalHeight();
        for( Int iLoc=0; iLoc<ALocalHeight; ++iLoc )
        {
            const Int off = A.EntryOffset( iLoc );
            const Int rowSize = A.NumConnections( iLoc );
            for( Int k=0; k<rowSize; ++k )
            {
                const Int colOff = meta.colOffs[k+off];
                const T AVal = A.Value(k+off);
                for( Int t=0; t<b; ++t )
                {
                    const T XVal = recvVals[colOff*b+t];
                    Y.UpdateLocal( iLoc, t, alpha*AVal*XVal );
                }
            }
        }
    }
    else
    {
        if( A.Width() != Y.Height() )
            LogicError("The width of A must match the height of Y");
        if( A.Height() != X.Height() )
            LogicError("The height of A must match the height of X");

        // Form and pack the updates to Y
        const bool conjugate = ( orientation == ADJOINT );
        std::vector<T> sendVals( meta.numRecvInds*b, 0 );
        const Int ALocalHeight = A.LocalHeight();
        for( Int iLoc=0; iLoc<ALocalHeight; ++iLoc )
        {
            const Int off = A.EntryOffset( iLoc );
            const Int rowSize = A.NumConnections( iLoc );
            for( Int k=0; k<rowSize; ++k )
            {
                const Int colOff = meta.colOffs[k+off];
                const T AVal = A.Value(k+off);
                for( Int t=0; t<b; ++t )
                {
                    const T XVal = X.GetLocal(iLoc,t);
                    if( conjugate )
                        sendVals[colOff*b+t] += alpha*Conj(AVal)*XVal;
                    else
                        sendVals[colOff*b+t] += alpha*AVal*XVal;
                }
            }
        }

        // Inject the updates to Y into the network
        const Int numRecvInds = meta.sendInds.size();
        std::vector<T> recvVals( numRecvInds*b );
        mpi::AllToAll
        ( sendVals.data(), recvSizes.data(), recvOffs.data(),
          recvVals.data(), sendSizes.data(), sendOffs.data(), comm );
     
        // Accumulate the received indices onto Y
        const Int firstLocalRow = Y.FirstLocalRow();
        for( Int s=0; s<numRecvInds; ++s )
        {
            const Int i = meta.sendInds[s];
            const Int iLoc = i - firstLocalRow;
            for( Int t=0; t<b; ++t )
                Y.UpdateLocal( iLoc, t, recvVals[s*b+t] );
        }
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
