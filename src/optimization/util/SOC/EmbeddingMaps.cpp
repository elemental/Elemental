/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace soc {

// TODO: Lower-level access

void EmbeddingMaps
( const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Int>& sparseOrders,
        Matrix<Int>& sparseFirstInds,
        Matrix<Int>& origToSparseOrders,
        Matrix<Int>& origToSparseFirstInds,
        Matrix<Int>& sparseToOrigOrders,
        Matrix<Int>& sparseToOrigFirstInds,
  Int cutoffSparse )
{
    DEBUG_ONLY(CSE cse("soc::EmbeddingMaps"))
    const Int k = orders.Height();

    // Form the metadata for the original index domain
    // -----------------------------------------------
    origToSparseOrders.Resize( k, 1 );
    origToSparseFirstInds.Resize( k, 1 );
    Int iSparse=0;
    for( Int i=0; i<k; )
    {
        const Int order = orders.Get(i,0);

        for( Int e=0; e<order; ++e )
            origToSparseFirstInds.Set( i+e, 0, iSparse );

        if( order > cutoffSparse )
        {
            for( Int e=0; e<order; ++e )
                origToSparseOrders.Set( i+e, 0, order+2 );
            iSparse += order+2;
        }
        else
        {
            for( Int e=0; e<order; ++e )
                origToSparseOrders.Set( i+e, 0, order );
            iSparse += order;
        }
        i += order;
    }
    const Int kSparse = iSparse;
    // Form the metadata for the sparsified index domain
    // -------------------------------------------------
    sparseToOrigOrders.Resize( kSparse, 1 );
    sparseToOrigFirstInds.Resize( kSparse, 1 );
    iSparse = 0;
    for( Int i=0; i<k; )
    {
        const Int order = orders.Get(i,0);
        if( order > cutoffSparse )
        {
            for( Int e=0; e<order+2; ++e )
                sparseToOrigOrders.Set( iSparse+e, 0, order );
            for( Int e=0; e<order+2; ++e )
                sparseToOrigFirstInds.Set( iSparse+e, 0, i );
            iSparse += order+2;
        }
        else
        {
            for( Int e=0; e<order; ++e )
                sparseToOrigOrders.Set( iSparse+e, 0, order );
            for( Int e=0; e<order; ++e )
                sparseToOrigFirstInds.Set( iSparse+e, 0, i );
            iSparse += order;
        }
        i += order;
    }
}

void EmbeddingMaps
( const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Int>& sparseOrders,
        DistMultiVec<Int>& sparseFirstInds,
        DistMultiVec<Int>& origToSparseOrders,
        DistMultiVec<Int>& origToSparseFirstInds,
        DistMultiVec<Int>& sparseToOrigOrders,
        DistMultiVec<Int>& sparseToOrigFirstInds,
  Int cutoffSparse )
{
    DEBUG_ONLY(CSE cse("soc::EmbeddingMaps"))
    const Int k = orders.Height();
    mpi::Comm comm = orders.Comm();
    const int commSize = mpi::Size( comm );

    const Int* orderBuf = orders.LockedMatrix().LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();
    const Int localHeight = orders.LocalHeight();

    // Allgather the list of cones with sufficiently large order
    // ---------------------------------------------------------
    // TODO: Send triplets instead?
    vector<Int> sendOrders, sendRoots;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = orders.GlobalRow(iLoc);
        const Int order = orderBuf[iLoc];
        const Int firstInd = firstIndBuf[iLoc];
        if( order > cutoffSparse && i == firstInd )
        {
            sendRoots.push_back( i );
            sendOrders.push_back( order );
        }
    }
    const int numSendRoots = sendRoots.size();
    vector<int> numRecvRoots(commSize);
    mpi::AllGather( &numSendRoots, 1, numRecvRoots.data(), 1, comm );
    vector<int> recvOffs;
    const int numRoots = Scan( numRecvRoots, recvOffs );
    // Receive the roots
    // ^^^^^^^^^^^^^^^^^
    vector<Int> recvRoots(numRoots);
    mpi::AllGather
    ( sendRoots.data(), numSendRoots,
      recvRoots.data(), numRecvRoots.data(), recvOffs.data(), comm );
    SwapClear( sendRoots );
    // Receive the orders
    // ^^^^^^^^^^^^^^^^^^
    vector<Int> recvOrders(numRoots);
    mpi::AllGather
    ( sendOrders.data(), numSendRoots,
      recvOrders.data(), numRecvRoots.data(), recvOffs.data(), comm );
    SwapClear( sendOrders );

    // TODO: Sort based upon the roots. The current distribution
    //       guarantees that they are already sorted.

    // Form the metadata for the original domain
    // -----------------------------------------

    origToSparseOrders.Resize( k, 1 );
    auto it = recvRoots.cbegin();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int order = orderBuf[iLoc];
        const Int root = firstIndBuf[iLoc];
        it = std::lower_bound( it, recvRoots.cend(), root );
        const bool embedded = ( it != recvRoots.cend() && *it == root );
        const Int sparseOrder = ( embedded ? order+2 : order );
        origToSparseOrders.SetLocal( iLoc, 0, sparseOrder );
    }

    origToSparseFirstInds.Resize( k, 1 );
    it = recvRoots.cbegin();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int root = firstIndBuf[iLoc];
        it = std::lower_bound( it, recvRoots.cend(), root );
        const Int numSparseBefore = Int(it-recvRoots.cbegin());
        const Int sparseOffset = 2*numSparseBefore;
        origToSparseFirstInds.SetLocal( iLoc, 0, root+sparseOffset );
    }

    // Form the metadata for the sparsified index domain
    // -------------------------------------------------
    Zeros( sparseOrders, k+2*numRoots, 1 );
    Zeros( sparseFirstInds, k+2*numRoots, 1 );
    Zeros( sparseToOrigOrders, k+2*numRoots, 1 );
    Zeros( sparseToOrigFirstInds, k+2*numRoots, 1 );

    // Count the number of updates for each
    // """"""""""""""""""""""""""""""""""""
    Int numQueues = 0;
    it = recvRoots.cbegin();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = firstInds.GlobalRow(iLoc);
        const Int root = firstIndBuf[iLoc];
        const Int order = orderBuf[iLoc];
        it = std::lower_bound( it, recvRoots.cend(), root );
        const bool embedded = ( it != recvRoots.cend() && root == *it );
        numQueues += 1;
        if( embedded && i - root == order-1 )
            numQueues += 2;
    }

    // Form sparseOrders
    // """""""""""""""""
    sparseOrders.Reserve( numQueues );
    it = recvRoots.cbegin();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = firstInds.GlobalRow(iLoc);
        const Int root = firstIndBuf[iLoc];
        const Int order = orderBuf[iLoc];
        it = std::lower_bound( it, recvRoots.cend(), root );
        const bool embedded = ( it != recvRoots.cend() && root == *it );
        const Int sparseOrder = ( embedded ? order+2 : order );
        const Int numSparseBefore = Int(it-recvRoots.cbegin());
        const Int sparseOffset = 2*numSparseBefore;
        const Int iSparse = i + sparseOffset;
        sparseOrders.QueueUpdate( iSparse, 0, sparseOrder );
        if( embedded && i - root == order-1 )
        {
            sparseOrders.QueueUpdate( iSparse+1, 0, sparseOrder );
           sparseOrders.QueueUpdate( iSparse+2, 0, sparseOrder );
        }
    }
    sparseOrders.ProcessQueues();

    // Form sparseFirstInds
    // """"""""""""""""""""
    sparseFirstInds.Reserve( numQueues );
    it = recvRoots.cbegin();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = firstInds.GlobalRow(iLoc);
        const Int root = firstIndBuf[iLoc];
        const Int order = orderBuf[iLoc];
        it = std::lower_bound( it, recvRoots.cend(), root );
        const bool embedded = ( it != recvRoots.cend() && root == *it );
        const Int numSparseBefore = Int(it-recvRoots.cbegin());
        const Int sparseOffset = 2*numSparseBefore;
        const Int iSparse = i + sparseOffset;
        const Int sparseRoot = root + sparseOffset;
        sparseFirstInds.QueueUpdate( iSparse, 0, sparseRoot );
        if( embedded && i - root == order-1 )
        {
            sparseFirstInds.QueueUpdate( iSparse+1, 0, sparseRoot );
            sparseFirstInds.QueueUpdate( iSparse+2, 0, sparseRoot );
        }
    }
    sparseFirstInds.ProcessQueues();

    // Form sparseToOrigOrders
    // """""""""""""""""""""""
    sparseToOrigOrders.Reserve( numQueues );
    it = recvRoots.cbegin();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = firstInds.GlobalRow(iLoc);
        const Int root = firstIndBuf[iLoc];
        const Int order = orderBuf[iLoc];
        it = std::lower_bound( it, recvRoots.cend(), root );
        const bool embedded = ( it != recvRoots.cend() && root == *it );
        const Int numSparseBefore = Int(it-recvRoots.cbegin());
        const Int sparseOffset = 2*numSparseBefore;
        const Int iSparse = i + sparseOffset;
        sparseToOrigOrders.QueueUpdate( iSparse, 0, order );
        if( embedded && i - root == order-1 )
        {
            sparseToOrigOrders.QueueUpdate( iSparse+1, 0, order );
            sparseToOrigOrders.QueueUpdate( iSparse+2, 0, order );
        }
    }
    sparseToOrigOrders.ProcessQueues();

    // Form sparseToOrigFirstInds
    // """"""""""""""""""""""""""
    sparseToOrigFirstInds.Reserve( numQueues );
    it = recvRoots.cbegin();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = firstInds.GlobalRow(iLoc);
        const Int root = firstIndBuf[iLoc];
        const Int order = orderBuf[iLoc];
        it = std::lower_bound( it, recvRoots.cend(), root );
        const bool embedded = ( it != recvRoots.cend() && root == *it );
        const Int numSparseBefore = Int(it-recvRoots.cbegin());
        const Int sparseOffset = 2*numSparseBefore;
        const Int iSparse = i + sparseOffset;
        sparseToOrigFirstInds.QueueUpdate( iSparse, 0, root );
        if( embedded && i - root == order-1 )
        {
            sparseToOrigFirstInds.QueueUpdate( iSparse+1, 0, root );
            sparseToOrigFirstInds.QueueUpdate( iSparse+2, 0, root );
        }
    }
    sparseToOrigFirstInds.ProcessQueues();
}

} // namespace soc
} // namespace El
