/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Members of second-order cones are stored contiguously within the column
// vector x, with the corresponding order of the cone each member belongs to
// stored in the same index of 'order', and the first index of the cone 
// being listed in the same index of 'firstInd'.
template<typename Real>
void SOCDots
( const Matrix<Real>& x, 
  const Matrix<Real>& y,
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCDots"))
    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");
    if( y.Height() != x.Height() || y.Width() != x.Width() )
        LogicError("x and y must be the same size");
    Zeros( z, x.Height(), x.Width() );

    for( Int i=0; i<height; )
    {
        const Int order = orders.Get(i,0);
        const Int firstInd = firstInds.Get(i,0);
        if( i != firstInd )
            LogicError("Inconsistency in orders and firstInds");

        // Compute the inner-product between two SOC members and broadcast
        // the result over an equivalently-sized z_i
        const Real dot = Dot( x(IR(i,i+order),ALL), y(IR(i,i+order),ALL) );
        z.Set( i, 0, dot );

        i += order;
    }
}

// TODO: An alternate, trivial implementation to benchmark against
template<typename Real>
void SOCDots
( const AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Real>& yPre,
        AbstractDistMatrix<Real>& zPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCDots"))
    AssertSameGrids( xPre, yPre, zPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto yPtr = ReadProxy<Real,VC,STAR>(&yPre,ctrl);
    auto zPtr = WriteProxy<Real,VC,STAR>(&zPre,ctrl);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& y = *yPtr;
    auto& z = *zPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");
    if( y.Height() != x.Height() || y.Width() != x.Width() )
        LogicError("x and y must be the same size");

    z.SetGrid( x.Grid() );
    Zeros( z, x.Height(), x.Width() );

    const Int localHeight = x.LocalHeight();
    mpi::Comm comm = x.DistComm();
    const int commSize = mpi::Size(comm);

    // TODO: Find a better strategy

    // Handle all second-order cones with order <= cutoff
    // ==================================================
    // Count the number of sends and handle the root portion of the work
    // -----------------------------------------------------------------
    Int numRemoteUpdates = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int order = orders.GetLocal(iLoc,0);
        if( order > cutoff )
            continue;

        const Int firstInd = firstInds.GetLocal(iLoc,0);
        if( i == firstInd )
            z.SetLocal( iLoc, 0, x.GetLocal(iLoc,0)*y.GetLocal(iLoc,0) );
        else if( !z.IsLocal(firstInd,0) )
            ++numRemoteUpdates;
    }
    // Queue and process the remote updates
    // ------------------------------------
    z.Reserve( numRemoteUpdates );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int order = orders.GetLocal(iLoc,0);
        if( order > cutoff )
            continue;

        const Int firstInd = firstInds.GetLocal(iLoc,0);
        if( i != firstInd )
            z.QueueUpdate( firstInd, 0, x.GetLocal(iLoc,0)*y.GetLocal(iLoc,0) );
    }
    z.ProcessQueues();

    // Handle all of the second-order cones with order > cutoff
    // ========================================================
    // Allgather the list of cones with sufficiently large order
    // ---------------------------------------------------------
    vector<Int> sendPairs;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int order = orders.GetLocal(iLoc,0);
        const Int firstInd = firstInds.GetLocal(iLoc,0);
        if( order > cutoff && i == firstInd )
        {
            sendPairs.push_back( i );
            sendPairs.push_back( order );
        }
    }
    const int numSendInts = sendPairs.size();
    vector<int> numRecvInts(commSize);
    mpi::AllGather( &numSendInts, 1, numRecvInts.data(), 1, comm );
    vector<int> recvOffs;
    const int totalRecv = Scan( numRecvInts, recvOffs );
    vector<Int> recvPairs(totalRecv);
    mpi::AllGather
    ( sendPairs.data(), numSendInts,
      recvPairs.data(), numRecvInts.data(), recvOffs.data(), comm );
    for( Int largeCone=0; largeCone<totalRecv/2; ++largeCone )
    {
        const Int i     = recvPairs[2*largeCone+0];
        const Int order = recvPairs[2*largeCone+1];
        auto xCone = x( IR(i,i+order), ALL );
        auto yCone = y( IR(i,i+order), ALL );
        z.Set( i, 0, Dot(xCone,yCone) );
    }
}

// TODO: An alternate, trivial implementation to benchmark against
template<typename Real>
void SOCDots
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCDots"))

    // TODO: Check that the communicators are congruent
    mpi::Comm comm = x.Comm();
    const int commSize = mpi::Size(comm);
    const int localHeight = x.LocalHeight();

    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");
    if( y.Height() != x.Height() || y.Width() != x.Width() )
        LogicError("x and y must be the same size");

    z.SetComm( x.Comm() );
    Zeros( z, x.Height(), x.Width() );

    // Perform an mpi::AllToAll to collect all of the second-order cones of
    // order less than or equal to the cutoff at the root locations and 
    // individually handle the remainder 
    // TODO: Find a better strategy

    // Handle all second-order cones with order <= cutoff
    // ==================================================
    // Count the number of sends and handle the local portion of the work
    // ------------------------------------------------------------------
    Int numRemoteUpdates = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int order = orders.GetLocal(iLoc,0);
        if( order > cutoff )
            continue;

        const Int firstInd = firstInds.GetLocal(iLoc,0);
        if( i == firstInd )
            z.SetLocal( iLoc, 0, x.GetLocal(iLoc,0)*y.GetLocal(iLoc,0) );
        else if( !z.IsLocal(firstInd,0) )
            ++numRemoteUpdates;
    }
    // Queue and process the remote updates
    // ------------------------------------
    z.Reserve( numRemoteUpdates );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int order = orders.GetLocal(iLoc,0);
        if( order > cutoff )
            continue;

        const Int firstInd = firstInds.GetLocal(iLoc,0);
        if( i != firstInd )
            z.QueueUpdate( firstInd, 0, x.GetLocal(iLoc,0)*y.GetLocal(iLoc,0) );
    }
    z.ProcessQueues();

    // Handle all of the second-order cones with order > cutoff
    // ========================================================
    // Allgather the list of cones with sufficiently large order
    // ---------------------------------------------------------
    vector<Int> sendPairs;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int order = orders.GetLocal(iLoc,0);
        const Int firstInd = firstInds.GetLocal(iLoc,0);
        if( order > cutoff && i == firstInd )
        {
            sendPairs.push_back( i );
            sendPairs.push_back( order );
        }
    }
    const int numSendInts = sendPairs.size();
    vector<int> numRecvInts(commSize);
    mpi::AllGather( &numSendInts, 1, numRecvInts.data(), 1, comm );
    vector<int> recvOffs;
    const int totalRecv = Scan( numRecvInts, recvOffs ); 
    vector<Int> recvPairs(totalRecv);
    mpi::AllGather
    ( sendPairs.data(), numSendInts,
      recvPairs.data(), numRecvInts.data(), recvOffs.data(), comm );
    for( Int largeCone=0; largeCone<totalRecv/2; ++largeCone )
    {
        const Int i     = recvPairs[2*largeCone+0];
        const Int order = recvPairs[2*largeCone+1];

        // Compute x(i:i+order)^T y(i:i+order)
        Real localDot = 0;
        const Int iFirst = x.FirstLocalRow();
        const Int iLast = iFirst + x.LocalHeight();
        for( Int j=Max(iFirst,i); j<Min(iLast,i+order); ++j )
            localDot += x.GetLocal(j-iFirst,0)*y.GetLocal(j-iFirst,0);
        const int owner = z.Owner(i,0);
        if( z.IsLocal(i,0) )
        {
            const Real dot = mpi::Reduce( localDot, owner, x.Comm() );
            z.Set( i, 0, dot );
        }
        else
            mpi::Reduce( localDot, owner, x.Comm() );
    }
}

#define PROTO(Real) \
  template void SOCDots \
  ( const Matrix<Real>& x, \
    const Matrix<Real>& y, \
          Matrix<Real>& z, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCDots \
  ( const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCDots \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
