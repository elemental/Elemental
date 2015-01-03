/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// Basis Pursuit is the search for a solution to A x = b which minimizes
// the one norm of x, i.e.,
//
//   min || x ||_1
//   s.t. A x = b.
//
// Real instances of the problem are expressable as a Linear Program [1] via 
// the transformation
//
// min 1^T [u;v]
// s.t. [A, -A] [u; v] = b, [u; v] >= 0.
//
// Complex instances of Basis Pursuit require Second-Order Cone Programming.
//

// TODO: Extend the existing LP control parameters
// TODO: Extend the (upcoming) SOCP control parameters

namespace El {

template<typename Real>
void BP
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("BP"))
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Real> c, xHat, AHat;

    // c := ones(2*n,1)
    // ================
    Ones( c, 2*n, 1 );

    // \hat A := [A, -A]
    // =================
    Zeros( AHat, m, 2*n );
    auto AHatL = AHat( IR(0,m), IR(0,n) );
    auto AHatR = AHat( IR(0,m), IR(n,2*n) );
    AHatL = A;
    AHatR = A;
    Scale( Real(-1), AHatR );

    // Solve the direct LP
    // ===================
    Matrix<Real> y, z;
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    x = xHat( IR(0,n), IR(0,1) );
    Axpy( Real(-1), xHat(IR(n,2*n),IR(0,1)), x );
}

template<typename Real>
void BP
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        AbstractDistMatrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("BP"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<Real> c(g), xHat(g), AHat(g);

    // c := ones(2*n,1)
    // ================
    Ones( c, 2*n, 1 );

    // \hat A := [A, -A]
    // =================
    Zeros( AHat, m, 2*n );
    auto AHatL = AHat( IR(0,m), IR(0,n) );
    auto AHatR = AHat( IR(0,m), IR(n,2*n) );
    AHatL = A;
    AHatR = A;
    Scale( Real(-1), AHatR );

    // Solve the direct LP
    // ===================
    DistMatrix<Real> y(g), z(g);
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    Copy( xHat( IR(0,n), IR(0,1) ), x );
    Axpy( Real(-1), xHat(IR(n,2*n),IR(0,1)), x );
}

template<typename Real>
void BP
( const SparseMatrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("BP"))
    const Int m = A.Height();
    const Int n = A.Width();
    SparseMatrix<Real> AHat;
    Matrix<Real> c, xHat;

    // c := ones(2*n,1)
    // ================
    Ones( c, 2*n, 1 );

    // \hat A := [A, -A]
    // =================
    const Int numEntriesA = A.NumEntries();
    Zeros( AHat, m, 2*n );
    AHat.Reserve( 2*numEntriesA );
    for( Int e=0; e<numEntriesA; ++e )
    {
        AHat.QueueUpdate( A.Row(e), A.Col(e),    A.Value(e) );
        AHat.QueueUpdate( A.Row(e), A.Col(e)+n, -A.Value(e) );
    }
    AHat.MakeConsistent();

    // Solve the direct LP
    // ===================
    Matrix<Real> y, z;
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    Zeros( x, n, 1 );
    for( Int i=0; i<n; ++i )
        x.Set( i, 0, xHat.Get(i,0)-xHat.Get(i+n,0) );
}

template<typename Real>
void BP
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        DistMultiVec<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("BP"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    DistSparseMatrix<Real> AHat(comm);
    DistMultiVec<Real> c(comm), xHat(comm);

    // c := ones(2*n,1)
    // ================
    Ones( c, 2*n, 1 );

    // \hat A := [A, -A]
    // ================
    // NOTE: Since A and \hat A are the same height and each distributed within
    //       columns, it is possible to form \hat A from A without communication
    const Int numLocalEntriesA = A.NumLocalEntries();
    Zeros( AHat, m, 2*n );
    AHat.Reserve( 2*numLocalEntriesA );
    for( Int e=0; e<numLocalEntriesA; ++e )
    {
        AHat.QueueLocalUpdate
        ( A.Row(e)-A.FirstLocalRow(), A.Col(e),    A.Value(e) );
        AHat.QueueLocalUpdate
        ( A.Row(e)-A.FirstLocalRow(), A.Col(e)+n, -A.Value(e) );
    }
    AHat.MakeConsistent();

    // Solve the direct LP
    // ===================
    DistMultiVec<Real> y(comm), z(comm);
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    Zeros( x, n, 1 );
    // Determine the send and recv counts/offsets
    // ------------------------------------------
    const Int commSize = mpi::Size(comm);
    std::vector<int> sendCounts(commSize,0);
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
    {
        const Int i = xHat.GlobalRow(iLoc);
        if( i < n )
            ++sendCounts[ x.RowOwner(i) ];
        else
            ++sendCounts[ x.RowOwner(i-n) ];
    }
    std::vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    std::vector<int> sendOffsets, recvOffsets;
    const int totalSend = Scan( sendCounts, sendOffsets );
    const int totalRecv = Scan( recvCounts, recvOffsets );
    // Pack the data 
    // -------------
    std::vector<Int> sSendBuf(totalSend);
    std::vector<Real> vSendBuf(totalSend);
    auto offsets = sendOffsets;
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
    {
        const Int i = xHat.GlobalRow(iLoc);
        if( i < n )
        {
            const int owner = x.RowOwner(i);
            sSendBuf[offsets[owner]] = i;
            vSendBuf[offsets[owner]] = xHat.GetLocal(iLoc,0);
            ++offsets[owner];
        }
        else
        {
            const int owner = x.RowOwner(i-n);
            sSendBuf[offsets[owner]] = i-n;
            vSendBuf[offsets[owner]] = -xHat.GetLocal(iLoc,0);
            ++offsets[owner];
        }
    }
    // Exchange the data
    // -----------------
    std::vector<Int> sRecvBuf(totalRecv);
    std::vector<Real> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    // Unpack the data
    // ---------------
    for( Int e=0; e<totalRecv; ++e )
        x.UpdateLocal( sRecvBuf[e]-x.FirstLocalRow(), 0, vRecvBuf[e] );
}

#define PROTO(Real) \
  template void BP \
  ( const Matrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void BP \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
          AbstractDistMatrix<Real>& x, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void BP \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void BP \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, \
          DistMultiVec<Real>& x, \
    const lp::direct::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namepace elem
