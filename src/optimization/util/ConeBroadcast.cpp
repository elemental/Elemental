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
void ConeBroadcast
(       Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("ConeBroadcast"))
    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    Real* xBuf = x.Buffer();
    const Int* orderBuf = orders.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    for( Int i=0; i<height; )
    {
        const Int order = orderBuf[i];
        const Int firstInd = firstIndBuf[i];
        DEBUG_ONLY(
          if( i != firstInd )
              LogicError("Inconsistency in orders and firstInds");
        )

        const Real x0 = xBuf[i]; 
        for( Int j=i+1; j<i+order; ++j )
            xBuf[j] = x0;

        i += order;
    }
}

// TODO: Use lower-level access
template<typename Real>
void ConeBroadcast
(       ElementalMatrix<Real>& xPre, 
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("ConeBroadcast"))
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
    // Count the number of remote updates (and set non-root entries to zero)
    // ---------------------------------------------------------------------
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
    // ------------------------------------
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

    // Handle all of the cones with order > cutoff
    // ===========================================
    // Allgather the list of cones with sufficiently large order
    // ---------------------------------------------------------
    vector<Entry<Real>> sendData;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int order = orders.GetLocal(iLoc,0);
        const Int firstInd = firstInds.GetLocal(iLoc,0);
        if( order > cutoff && i == firstInd )
            sendData.push_back( Entry<Real>{i,order,x.GetLocal(iLoc,0)} );
    }
    const int numSendCones = sendData.size();
    vector<int> numRecvCones(commSize);
    mpi::AllGather( &numSendCones, 1, numRecvCones.data(), 1, comm );
    vector<int> recvOffs;
    const int totalRecv = Scan( numRecvCones, recvOffs );
    vector<Entry<Real>> recvData(totalRecv);
    mpi::AllGather
    ( sendData.data(), numSendCones,
      recvData.data(), numRecvCones.data(), recvOffs.data(), comm );
    for( Int largeCone=0; largeCone<totalRecv; ++largeCone )
    {
        const Int i = recvData[largeCone].i;
        const Int order = recvData[largeCone].j;
        const Real x0 = recvData[largeCone].value;
        auto xBot = x( IR(i+1,i+order), ALL );
        Fill( xBot, x0 );
    }
}

template<typename Real>
void ConeBroadcast
(       DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("ConeBroadcast"))

    // TODO: Check that the communicators are congruent
    mpi::Comm comm = x.Comm();
    const int commSize = mpi::Size(comm);
    const Int localHeight = x.LocalHeight();
    const Int firstLocalRow = x.FirstLocalRow();

    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    Real* xBuf = x.Matrix().Buffer();
    const Int* orderBuf = orders.LockedMatrix().LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();

    // TODO: Find a better strategy
    // A short-circuited ring algorithm would likely be significantly faster

    // Handle all cones with order <= cutoff
    // =====================================
    // Count the number of remote updates (and set non-root entries to zero)
    // ---------------------------------------------------------------------
    Int numRemoteUpdates = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = iLoc + firstLocalRow;
        const Int order = orderBuf[iLoc];
        if( order > cutoff )
            continue;

        const Int firstInd = firstIndBuf[iLoc];
        if( i == firstInd )
        {
            for( Int k=1; k<order; ++k )
                if( !x.IsLocal(i+k,0) )
                    ++numRemoteUpdates;
        }
        else
            xBuf[iLoc] = 0;
    }
    // Queue and process the remote updates
    // ------------------------------------
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = iLoc + firstLocalRow;
        const Int order = orderBuf[iLoc];
        if( order > cutoff )
            continue;

        const Int firstInd = firstIndBuf[iLoc];
        if( i == firstInd )
            for( Int k=1; k<order; ++k )
                x.QueueUpdate( i+k, 0, xBuf[iLoc] );
    }
    x.ProcessQueues();

    // Handle all of the cones with order > cutoff
    // ===========================================
    // Allgather the list of cones with sufficiently large order
    // ---------------------------------------------------------
    vector<Entry<Real>> sendData;
    // TODO: Count and reserve
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = iLoc + firstLocalRow;
        const Int order = orderBuf[iLoc];
        const Int firstInd = firstIndBuf[iLoc];
        if( order > cutoff && i == firstInd )
            sendData.push_back( Entry<Real>{i,order,xBuf[iLoc]} );
    }
    const int numSendCones = sendData.size();
    vector<int> numRecvCones(commSize);
    mpi::AllGather( &numSendCones, 1, numRecvCones.data(), 1, comm );
    vector<int> recvOffs;
    const int totalRecv = Scan( numRecvCones, recvOffs ); 
    vector<Entry<Real>> recvData(totalRecv);
    mpi::AllGather
    ( sendData.data(), numSendCones,
      recvData.data(), numRecvCones.data(), recvOffs.data(), comm );
    for( Int largeCone=0; largeCone<totalRecv; ++largeCone )
    {
        const Int i = recvData[largeCone].i;
        const Int order = recvData[largeCone].j;
        const Real x0 = recvData[largeCone].value;

        // Unpack the root entry
        const Int iFirst = firstLocalRow;
        const Int iLast = iFirst + localHeight;
        for( Int j=Max(iFirst,i+1); j<Min(iLast,i+order); ++j )
            xBuf[j-iFirst] = x0;
    }
}

#define PROTO(Real) \
  template void ConeBroadcast \
  (       Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void ConeBroadcast \
  (       ElementalMatrix<Real>& x, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, Int cutoff ); \
  template void ConeBroadcast \
  (       DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
