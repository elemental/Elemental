/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace {

template<typename Real>
EnableIf<IsReal<Real>,function<Real(Real,Real)>>
OpToReduce( mpi::Op op )
{
    function<Real(Real,Real)> reduce;
    if( op == mpi::SUM )
        reduce = []( Real alpha, Real beta ) { return alpha+beta; };
    else if( op == mpi::MAX )
        reduce = []( Real alpha, Real beta ) { return Max(alpha,beta); };
    else if( op == mpi::MIN )
        reduce = []( Real alpha, Real beta ) { return Min(alpha,beta); };
    else
        LogicError("Unsupported cone::AllReduce operation");
    return reduce;
}

template<typename F>
EnableIf<IsComplex<F>,function<F(F,F)>>
OpToReduce( mpi::Op op )
{
    function<F(F,F)> reduce;
    if( op == mpi::SUM )
        reduce = []( F alpha, F beta ) { return alpha+beta; };
    else
        LogicError("Unsupported cone::AllReduce operation");
    return reduce;
}

} // anonymous namespace

namespace cone {

template<typename F>
void AllReduce
(       Matrix<F>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
        mpi::Op op )
{
    DEBUG_CSE
    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    auto reduce = OpToReduce<F>( op );
    for( Int i=0; i<height; )
    {
        const Int order = orders(i);
        const Int firstInd = firstInds(i);
        if( i != firstInd )
            LogicError("Inconsistency in orders and firstInds");

        F coneRes = x(i);
        for( Int j=i+1; j<i+order; ++j )
            coneRes = reduce(coneRes,x(j));
        for( Int j=i; j<i+order; ++j )
            x(j) = coneRes;

        i += order;
    }
}

template<typename F>
void AllReduce
(       ElementalMatrix<F>& xPre, 
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  mpi::Op op,
  Int cutoff )
{
    DEBUG_CSE
    AssertSameGrids( xPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<F,F,VC,STAR>
      xProx( xPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    auto reduce = OpToReduce<F>( op );

    const Int localHeight = x.LocalHeight();
    mpi::Comm comm = x.DistComm();
    const int commSize = mpi::Size(comm);

          F* xBuf = x.Buffer();
    const Int* orderBuf = orders.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    // Perform an mpi::AllToAll to scatter all of the cone roots of
    // order less than or equal to the cutoff 
    // TODO: Find a better strategy
    // A short-circuited ring algorithm would likely be significantly faster

    // Handle all cones with order <= cutoff
    // =====================================

    // Compute the reduction of the cones on the roots
    // -----------------------------------------------
    {
        // For now, we will simply Gather each cone to each root rather than
        // performing a local reduction beforehand (TODO: do this)

        // Compute the metadata
        // ^^^^^^^^^^^^^^^^^^^^
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orderBuf[iLoc];
            if( order > cutoff )
                continue;

            const Int firstInd = firstIndBuf[iLoc];
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
        vector<F> sendBuf(totalSend);
        auto offs = sendOffs;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orderBuf[iLoc];
            if( order > cutoff )
                continue;

            const Int firstInd = firstIndBuf[iLoc];
            if( i != firstInd )
            {
                const Int owner = firstInds.RowOwner(firstInd);
                // TODO: Don't pack if we also own the root
                sendBuf[offs[owner]++] = xBuf[iLoc];
            }
        }

        // Exchange the data
        // ^^^^^^^^^^^^^^^^^ 
        vector<int> recvOffs;
        const int totalRecv = Scan( recvCounts, recvOffs );
        vector<F> recvBuf(totalRecv);
        mpi::AllToAll
        ( sendBuf.data(), sendCounts.data(), sendOffs.data(),
          recvBuf.data(), recvCounts.data(), recvOffs.data(), x.DistComm() );

        // Locally reduce on the roots
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^
        offs = recvOffs;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orderBuf[iLoc];
            if( order > cutoff )
                continue;

            const Int firstInd = firstIndBuf[iLoc];
            if( i == firstInd )
            {
                F coneRes = xBuf[iLoc];
                for( Int j=i+1; j<i+order; ++j )
                {
                    const Int owner = firstInds.RowOwner(j);
                    // TODO: Pull locally if locally owned
                    coneRes = reduce(coneRes,recvBuf[offs[owner]++]);
                }
                xBuf[iLoc] = coneRes;
            }
        }
    }

    // Broadcast the results from the roots of the cones
    // -------------------------------------------------
    {
        // Count the number of remote updates (and set non-root entries to zero)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        Int numRemoteUpdates = 0;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
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
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orderBuf[iLoc];
            if( order > cutoff )
                continue;

            const Int firstInd = firstIndBuf[iLoc];
            if( i == firstInd )
                for( Int k=1; k<order; ++k )
                    x.QueueUpdate( i+k, 0, xBuf[iLoc] );
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
        const Int order = orderBuf[iLoc];
        const Int firstInd = firstIndBuf[iLoc];
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
        F* xConeBuf = xCone.Buffer();

        // Compute our local result in this cone
        F localConeRes = 0;
        const Int xConeLocalHeight = xCone.LocalHeight();
        for( Int iLoc=0; iLoc<xConeLocalHeight; ++iLoc )
            localConeRes = reduce(localConeRes,xConeBuf[iLoc]);

        // Compute the maximum for this cone
        const F coneRes = mpi::AllReduce( localConeRes, op, x.DistComm() );
        for( Int iLoc=0; iLoc<xConeLocalHeight; ++iLoc )
            xConeBuf[iLoc] = coneRes;
    }
}

template<typename F>
void AllReduce
(       DistMultiVec<F>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, 
  mpi::Op op, Int cutoff )
{
    DEBUG_CSE

    // TODO: Check that the communicators are congruent
    mpi::Comm comm = x.Comm();
    const int commSize = mpi::Size(comm);
    const int localHeight = x.LocalHeight();

    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    auto reduce = OpToReduce<F>( op );

          F* xBuf = x.Matrix().Buffer();
    const Int* orderBuf = orders.LockedMatrix().LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();

    // TODO: Find a better strategy
    // A short-circuited ring algorithm would likely be significantly faster

    // Handle all cones with order <= cutoff
    // =====================================

    // Compute the reduction of the cones on the roots
    // -----------------------------------------------
    {
        // For now, we will simply Gather each cone to each root rather than
        // performing a local reduction beforehand (TODO: do this)

        // Compute the metadata
        // ^^^^^^^^^^^^^^^^^^^^
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orderBuf[iLoc];
            if( order > cutoff )
                continue;

            const Int firstInd = firstIndBuf[iLoc];
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
        vector<F> sendBuf(totalSend);
        auto offs = sendOffs;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orderBuf[iLoc];
            if( order > cutoff )
                continue;

            const Int firstInd = firstIndBuf[iLoc];
            if( i != firstInd )
            {
                const Int owner = firstInds.RowOwner(firstInd);
                // TODO: Don't pack if we also own the root
                sendBuf[offs[owner]++] = xBuf[iLoc];
            }
        }

        // Exchange the data
        // ^^^^^^^^^^^^^^^^^ 
        vector<int> recvOffs;
        const int totalRecv = Scan( recvCounts, recvOffs );
        vector<F> recvBuf(totalRecv);
        mpi::AllToAll
        ( sendBuf.data(), sendCounts.data(), sendOffs.data(),
          recvBuf.data(), recvCounts.data(), recvOffs.data(), x.Comm() );

        // Compute the maxima on the roots
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        offs = recvOffs;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orderBuf[iLoc];
            if( order > cutoff )
                continue;

            const Int firstInd = firstIndBuf[iLoc];
            if( i == firstInd )
            {
                F coneRes = xBuf[iLoc];
                for( Int j=i+1; j<i+order; ++j )
                {
                    const Int owner = firstInds.RowOwner(j);
                    // TODO: Pull locally if locally owned
                    coneRes = reduce(coneRes,recvBuf[offs[owner]++]);
                }
                xBuf[iLoc] = coneRes;
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
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = x.GlobalRow(iLoc);
            const Int order = orderBuf[iLoc];
            if( order > cutoff )
                continue;

            const Int firstInd = firstIndBuf[iLoc];
            if( i == firstInd )
                for( Int k=1; k<order; ++k )
                    x.QueueUpdate( i+k, 0, xBuf[iLoc] );
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
        const Int order = orderBuf[iLoc];
        const Int firstInd = firstIndBuf[iLoc];
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

        // Compute our local result in this cone
        F localConeRes = 0;
        const Int iFirst = x.FirstLocalRow();
        const Int iLast = iFirst + x.LocalHeight();
        for( Int j=Max(iFirst,i); j<Min(iLast,i+order); ++j )
            localConeRes = reduce(localConeRes,xBuf[j-iFirst]);

        // Compute the maximum for this cone
        const F coneRes = mpi::AllReduce( localConeRes, op, x.Comm() );
        for( Int j=Max(iFirst,i); j<Min(iLast,i+order); ++j )
            xBuf[j-iFirst] = coneRes;
    }
}

#define PROTO(F) \
  template void AllReduce \
  (       Matrix<F>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    mpi::Op op ); \
  template void AllReduce \
  (       ElementalMatrix<F>& x, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
    mpi::Op op, Int cutoff ); \
  template void AllReduce \
  (       DistMultiVec<F>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    mpi::Op op, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace cone
} // namespace El
