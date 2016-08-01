/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace soc {

// Members of second-order cones are stored contiguously within the column
// vector x, with the corresponding order of the cone each member belongs to
// stored in the same index of 'order', and the first index of the cone 
// being listed in the same index of 'firstInd'.
template<typename Real,typename>
void Dots
( const Matrix<Real>& x, 
  const Matrix<Real>& y,
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_CSE
    const Int height = x.Height();
    DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
      if( y.Height() != x.Height() || y.Width() != x.Width() )
          LogicError("x and y must be the same size");
    )
    Zeros( z, x.Height(), x.Width() );

    for( Int i=0; i<height; )
    {
        const Int order = orders(i);
        const Int firstInd = firstInds(i);
        DEBUG_ONLY(
          if( i != firstInd )
              LogicError("Inconsistency in orders and firstInds");
        )

        // Compute the inner-product between two SOC members and broadcast
        // the result over an equivalently-sized z_i
        z(i) = blas::Dot( order, &x(i), 1, &y(i), 1 );
        i += order;
    }
}

// TODO: An alternate, trivial implementation to benchmark against
template<typename Real,typename>
void Dots
( const ElementalMatrix<Real>& xPre, 
  const ElementalMatrix<Real>& yPre,
        ElementalMatrix<Real>& zPre,
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_CSE
    AssertSameGrids( xPre, yPre, zPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Real,Real,VC,STAR>
      xProx( xPre, ctrl ),
      yProx( yPre, ctrl );
    DistMatrixWriteProxy<Real,Real,VC,STAR>
      zProx( zPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.GetLocked();
    auto& y = yProx.GetLocked();
    auto& z = zProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    const Int height = x.Height();
    DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
      if( y.Height() != x.Height() || y.Width() != x.Width() )
          LogicError("x and y must be the same size");
    )

    z.SetGrid( x.Grid() );
    Zeros( z, x.Height(), x.Width() );

    const Int localHeight = x.LocalHeight();
    mpi::Comm comm = x.DistComm();
    const int commSize = mpi::Size(comm);

    const Int* orderBuf = orders.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();
    const Real* xBuf = x.LockedBuffer();
    const Real* yBuf = y.LockedBuffer();
          Real* zBuf = z.Buffer();

    // TODO: Find a better strategy

    // Handle all second-order cones with order <= cutoff
    // ==================================================
    // Count the number of sends and handle the root portion of the work
    // -----------------------------------------------------------------
    Int numRemoteUpdates = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int order = orderBuf[iLoc];
        if( order > cutoff )
            continue;

        const Int firstInd = firstIndBuf[iLoc];
        if( i == firstInd )
            zBuf[iLoc] = xBuf[iLoc]*yBuf[iLoc];
        else if( !z.IsLocal(firstInd,0) )
            ++numRemoteUpdates;
    }
    // Queue and process the remote updates
    // ------------------------------------
    z.Reserve( numRemoteUpdates );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int order = orderBuf[iLoc];
        if( order > cutoff )
            continue;

        const Int firstInd = firstIndBuf[iLoc];
        if( i != firstInd )
            z.QueueUpdate( firstInd, 0, xBuf[iLoc]*yBuf[iLoc] );
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
        const Int order = orderBuf[iLoc];
        const Int firstInd = firstIndBuf[iLoc];
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
template<typename Real,typename>
void Dots
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
        Int cutoff )
{
    DEBUG_CSE

    // TODO: Check that the communicators are congruent
    mpi::Comm comm = x.Comm();
    const int commSize = mpi::Size(comm);
    const Int localHeight = x.LocalHeight();
    const Int firstLocalRow = x.FirstLocalRow();

    const Int height = x.Height();
    DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
      if( y.Height() != x.Height() || y.Width() != x.Width() )
          LogicError("x and y must be the same size");
    )

    z.SetComm( x.Comm() );
    Zeros( z, x.Height(), x.Width() );

    const Real* xBuf = x.LockedMatrix().LockedBuffer();
    const Real* yBuf = y.LockedMatrix().LockedBuffer();
          Real* zBuf = z.Matrix().Buffer();
    const Int* orderBuf = orders.LockedMatrix().LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();

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
        const Int i = iLoc + firstLocalRow;
        const Int order = orderBuf[iLoc];
        if( order > cutoff )
            continue;

        const Int firstInd = firstIndBuf[iLoc];
        if( i == firstInd )
            zBuf[iLoc] = xBuf[iLoc]*yBuf[iLoc];
        else if( !z.IsLocal(firstInd,0) )
            ++numRemoteUpdates;
    }
    // Queue and process the remote updates
    // ------------------------------------
    z.Reserve( numRemoteUpdates );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = iLoc + firstLocalRow;
        const Int order = orderBuf[iLoc];
        if( order > cutoff )
            continue;

        const Int firstInd = firstIndBuf[iLoc];
        if( i != firstInd )
            z.QueueUpdate( firstInd, 0, xBuf[iLoc]*yBuf[iLoc] );
    }
    z.ProcessQueues();

    // Handle all of the second-order cones with order > cutoff
    // ========================================================
    // Allgather the list of cones with sufficiently large order
    // ---------------------------------------------------------
    vector<Int> sendPairs;
    // TODO: Count and reserve so that the push_back's are all O(1)
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = iLoc + firstLocalRow;
        const Int order = orderBuf[iLoc];
        const Int firstInd = firstIndBuf[iLoc];
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
    const Int iFirst = firstLocalRow;
    const Int iLast = iFirst + localHeight;
    for( Int largeCone=0; largeCone<totalRecv/2; ++largeCone )
    {
        const Int i     = recvPairs[2*largeCone+0];
        const Int order = recvPairs[2*largeCone+1];

        // Compute x(i:i+order)^T y(i:i+order)
        Real localDot = 0;
        const Int jBeg = Max(iFirst,i);
        const Int jEnd = Min(iLast,i+order);
        for( Int j=jBeg; j<jEnd; ++j )
            localDot += xBuf[j-iFirst]*yBuf[j-iFirst];
        const int owner = z.Owner(i,0);
        if( z.IsLocal(i,0) )
        {
            const Real dot = mpi::Reduce( localDot, owner, x.Comm() );
            zBuf[i-firstLocalRow] = dot;
        }
        else
            mpi::Reduce( localDot, owner, x.Comm() );
    }
}

#define PROTO(Real) \
  template void Dots \
  ( const Matrix<Real>& x, \
    const Matrix<Real>& y, \
          Matrix<Real>& z, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void Dots \
  ( const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Real>& y, \
          ElementalMatrix<Real>& z, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
    Int cutoff ); \
  template void Dots \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace soc
} // namespace El
