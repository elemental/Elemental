/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
void ConeMax
(       Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("ConeMax"))
    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    for( Int i=0; i<height; )
    {
        const Int order = orders.Get(i,0);
        const Int firstInd = firstInds.Get(i,0);
        if( i != firstInd )
            LogicError("Inconsistency in orders and firstInds");

        Real coneMax = x.Get(i,0);
        for( Int j=i+1; j<i+order; ++j )
            coneMax = Max(coneMax,x.Get(j,0));
        for( Int j=i; j<i+order; ++j )
            x.Set(j,0,coneMax);

        i += order;
    }
}

template<typename Real>
void ConeMax
(       AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("ConeMax"))
    AssertSameGrids( xPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    const Int localHeight = x.LocalHeight();
    mpi::Comm comm = x.DistComm();
    const int commSize = mpi::Size(comm);

    // Perform an mpi::AllToAll to scatter all of the cone roots of
    // order less than or equal to the cutoff 
    // TODO: Find a better strategy
    // A short-circuited ring algorithm would likely be significantly faster

    // Handle all cones with order <= cutoff
    // =====================================

    // Compute the maxima of the cones on the roots
    // --------------------------------------------
    {
        // For now, we will simply Gather each cone to each root rather than
        // performing a local reduction beforehand (TODO: do this)

        // Compute the metadata
        // ^^^^^^^^^^^^^^^^^^^^
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            if( order > cutoff )
                continue;

            const Int firstInd = firstInds.GetLocal(iLoc,0);
            if( i != firstInd )
            {
                const Int owner = firstInds.RowOwner(firstInd);
                // TODO: Don't count this if we also own the root
                ++sendCounts[owner];
            }
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll
        ( sendCounts.data(), 1, recvCounts.data(), 1, x.DistComm() );
        
        // Pack the data
        // ^^^^^^^^^^^^^
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        vector<Real> sendBuf(totalSend);
        auto offs = sendOffs;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            if( order > cutoff )
                continue;

            const Int firstInd = firstInds.GetLocal(iLoc,0);
            if( i != firstInd )
            {
                const Int owner = firstInds.RowOwner(firstInd);
                // TODO: Don't pack if we also own the root
                sendBuf[offs[owner]++] = x.GetLocal(iLoc,0);
            }
        }

        // Exchange the data
        // ^^^^^^^^^^^^^^^^^ 
        vector<int> recvOffs;
        const int totalRecv = Scan( recvCounts, recvOffs );
        vector<Real> recvBuf(totalRecv);
        mpi::AllToAll
        ( sendBuf.data(), sendCounts.data(), sendOffs.data(),
          recvBuf.data(), recvCounts.data(), recvOffs.data(), x.DistComm() );

        // Compute the maxima on the roots
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        offs = recvOffs;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            if( order > cutoff )
                continue;

            const Int firstInd = firstInds.GetLocal(iLoc,0);
            if( i == firstInd )
            {
                Real coneMax = x.GetLocal(iLoc,0);
                for( Int j=i+1; j<i+order; ++j )
                {
                    const Int owner = firstInds.RowOwner(j);
                    // TODO: Pull locally if locally owned
                    coneMax = Max(coneMax,recvBuf[offs[owner]++]);
                }
                x.SetLocal(iLoc,0,coneMax);
            }
        }
    }

    // Broadcast the maxima from the roots of the cones
    // ------------------------------------------------
    {
        // Count the number of remote updates (and set non-root entries to zero)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        Int numRemoteUpdates = 0;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            if( order > cutoff )
                continue;

            const Int firstInd = firstInds.GetLocal(iLoc,0);
            if( i == firstInd )
            {
                for( Int k=1; k<order; ++k )
                    if( !x.IsLocal(i+k,0) )
                        ++numRemoteUpdates;
            }
            else
                x.SetLocal( iLoc, 0, 0 );
        }
        // Queue and process the remote updates
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            if( order > cutoff )
                continue;

            const Int firstInd = firstInds.GetLocal(iLoc,0);
            if( i == firstInd )
                for( Int k=1; k<order; ++k )
                    x.QueueUpdate( i+k, 0, x.GetLocal(iLoc,0) );
        }
        x.ProcessQueues();
    }

    // Handle all of the cones with order > cutoff
    // ===========================================
    // Allgather the list of cones with sufficiently large order
    // ---------------------------------------------------------
    vector<Int> sendData;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int order = orders.GetLocal(iLoc,0);
        const Int firstInd = firstInds.GetLocal(iLoc,0);
        if( order > cutoff && i == firstInd )
        {
            sendData.push_back( i );
            sendData.push_back( order );
        }
    }
    const int numSendInts = sendData.size();
    vector<int> numRecvInts(commSize);
    mpi::AllGather( &numSendInts, 1, numRecvInts.data(), 1, comm );
    vector<int> recvOffs;
    const int totalRecv = Scan( numRecvInts, recvOffs );
    vector<Int> recvData(totalRecv);
    mpi::AllGather
    ( sendData.data(), numSendInts,
      recvData.data(), numRecvInts.data(), recvOffs.data(), comm );
    for( Int largeCone=0; largeCone<totalRecv/2; ++largeCone )
    {
        const Int i = recvData[2*largeCone+0];
        const Int order = recvData[2*largeCone+1];
        auto xCone = x( IR(i,i+order), ALL );
        const Real maxNorm = MaxNorm( xCone );
        Fill( xCone, maxNorm );
    }
}

template<typename Real>
void ConeMax
(       DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("ConeMax"))

    // TODO: Check that the communicators are congruent
    mpi::Comm comm = x.Comm();
    const int commSize = mpi::Size(comm);
    const int localHeight = x.LocalHeight();

    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    // TODO: Find a better strategy
    // A short-circuited ring algorithm would likely be significantly faster

    // Handle all cones with order <= cutoff
    // =====================================

    // Compute the maxima of the cones on the roots
    // --------------------------------------------
    {
        // For now, we will simply Gather each cone to each root rather than
        // performing a local reduction beforehand (TODO: do this)

        // Compute the metadata
        // ^^^^^^^^^^^^^^^^^^^^
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            if( order > cutoff )
                continue;

            const Int firstInd = firstInds.GetLocal(iLoc,0);
            if( i != firstInd )
            {
                const Int owner = firstInds.RowOwner(firstInd);
                // TODO: Don't count this if we also own the root
                ++sendCounts[owner];
            }
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, x.Comm() );
        
        // Pack the data
        // ^^^^^^^^^^^^^
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        vector<Real> sendBuf(totalSend);
        auto offs = sendOffs;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            if( order > cutoff )
                continue;

            const Int firstInd = firstInds.GetLocal(iLoc,0);
            if( i != firstInd )
            {
                const Int owner = firstInds.RowOwner(firstInd);
                // TODO: Don't pack if we also own the root
                sendBuf[offs[owner]++] = x.GetLocal(iLoc,0);
            }
        }

        // Exchange the data
        // ^^^^^^^^^^^^^^^^^ 
        vector<int> recvOffs;
        const int totalRecv = Scan( recvCounts, recvOffs );
        vector<Real> recvBuf(totalRecv);
        mpi::AllToAll
        ( sendBuf.data(), sendCounts.data(), sendOffs.data(),
          recvBuf.data(), recvCounts.data(), recvOffs.data(), x.Comm() );

        // Compute the maxima on the roots
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        offs = recvOffs;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            if( order > cutoff )
                continue;

            const Int firstInd = firstInds.GetLocal(iLoc,0);
            if( i == firstInd )
            {
                Real coneMax = x.GetLocal(iLoc,0);
                for( Int j=i+1; j<i+order; ++j )
                {
                    const Int owner = firstInds.RowOwner(j);
                    // TODO: Pull locally if locally owned
                    coneMax = Max(coneMax,recvBuf[offs[owner]++]);
                }
                x.SetLocal(iLoc,0,coneMax);
            }
        }
    }

    // Broadcast from the roots of the small cones
    // -------------------------------------------
    {
        // Count the number of remote updates (and set non-root entries to zero)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        Int numRemoteUpdates = 0;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            if( order > cutoff )
                continue;

            const Int firstInd = firstInds.GetLocal(iLoc,0);
            if( i == firstInd )
            {
                for( Int k=1; k<order; ++k )
                    if( !x.IsLocal(i+k,0) )
                        ++numRemoteUpdates;
            }
            else
                x.SetLocal( iLoc, 0, 0 );
        }
        // Queue and process the remote updates
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orders.GetLocal(iLoc,0);
            if( order > cutoff )
                continue;

            const Int firstInd = firstInds.GetLocal(iLoc,0);
            if( i == firstInd )
                for( Int k=1; k<order; ++k )
                    x.QueueUpdate( i+k, 0, x.GetLocal(iLoc,0) );
        }
        x.ProcessQueues();
    }

    // Handle all of the cones with order > cutoff
    // ===========================================
    // Allgather the list of cones with sufficiently large order
    // ---------------------------------------------------------
    vector<Int> sendData;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int order = orders.GetLocal(iLoc,0);
        const Int firstInd = firstInds.GetLocal(iLoc,0);
        if( order > cutoff && i == firstInd )
        {
            sendData.push_back( i );
            sendData.push_back( order );
        }
    }
    const int numSendInts = sendData.size();
    vector<int> numRecvInts(commSize);
    mpi::AllGather( &numSendInts, 1, numRecvInts.data(), 1, comm );
    vector<int> recvOffs;
    const int totalRecv = Scan( numRecvInts, recvOffs ); 
    vector<Int> recvData(totalRecv);
    mpi::AllGather
    ( sendData.data(), numSendInts,
      recvData.data(), numRecvInts.data(), recvOffs.data(), comm );
    for( Int largeCone=0; largeCone<totalRecv/2; ++largeCone )
    {
        const Int i = recvData[2*largeCone+0];
        const Int order = recvData[2*largeCone+1];

        // Compute our local maximum in this cone
        Real localConeMax = 0;
        const Int iFirst = x.FirstLocalRow();
        const Int iLast = iFirst + x.LocalHeight();
        for( Int j=Max(iFirst,i); j<Min(iLast,i+order); ++j )
            localConeMax = Max(localConeMax,x.GetLocal(j-iFirst,0));

        // Compute the maximum for this cone
        const Real coneMax = mpi::AllReduce( localConeMap, mpi::MAX, x.Comm() );
        
        // Overwrite the entire cone with this maximum
        for( Int j=Max(iFirst,i); j<Min(iLast,i+order); ++j )
            x.SetLocal(j-iFirst,0,coneMax);
    }
}

#define PROTO(Real) \
  template void ConeMax \
  (       Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void ConeMax \
  (       AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void ConeMax \
  (       DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
