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
// decomposing x into its positive and negative parts, say (u,v), and posing
//
//   min 1^T [u;v]
//   s.t. [A, -A] [u; v] = b, [u; v] >= 0.
//
// After solving this LP, the solution is set to x := u - v.
//
// Complex instances of Basis Pursuit require Second-Order Cone Programming.
//
// [1] Scott S. Chen, David L. Donoho, and Michael A. Saunders,
//     "Atomic Decomposition by Basis Pursuit",
//     SIAM Review, Vol. 43, No. 1, pp. 129--159, 2001

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
    const Range<Int> uInd(0,n), vInd(n,2*n);
    Matrix<Real> c, AHat;

    // c := ones(2*n,1)
    // ================
    Ones( c, 2*n, 1 );

    // \hat A := [A, -A]
    // =================
    Zeros( AHat, m, 2*n );
    auto AHatu = AHat( IR(0,m), uInd );
    auto AHatv = AHat( IR(0,m), vInd );
    AHatu = A;
    Axpy( Real(-1), A, AHatv );

    // Solve the direct LP
    // ===================
    Matrix<Real> xHat, y, z;
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, IR(0,1) );
    Axpy( Real(-1), xHat(vInd,IR(0,1)), x );
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
    const Range<Int> uInd(0,n), vInd(n,2*n);
    DistMatrix<Real> c(g), AHat(g);

    // c := ones(2*n,1)
    // ================
    Ones( c, 2*n, 1 );

    // \hat A := [A, -A]
    // =================
    Zeros( AHat, m, 2*n );
    auto AHatu = AHat( IR(0,m), uInd );
    auto AHatv = AHat( IR(0,m), vInd );
    AHatu = A;
    Axpy( Real(-1), A, AHatv );

    // Solve the direct LP
    // ===================
    DistMatrix<Real> xHat(g), y(g), z(g);
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    Copy( xHat( uInd, IR(0,1) ), x );
    Axpy( Real(-1), xHat(vInd,IR(0,1)), x );
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
    const Range<Int> uInd(0,n), vInd(n,2*n);
    SparseMatrix<Real> AHat;
    Matrix<Real> c;

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
    Matrix<Real> xHat, y, z;
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, IR(0,1) );
    Axpy( Real(-1), xHat(vInd,IR(0,1)), x );
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
    DistMultiVec<Real> c(comm);

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
    DistMultiVec<Real> xHat(comm), y(comm), z(comm);
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    Zeros( x, n, 1 );
    // Determine the send and recv counts/offs
    // ------------------------------------------
    const Int commSize = mpi::Size(comm);
    vector<int> sendSizes(commSize,0);
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
    {
        const Int i = xHat.GlobalRow(iLoc);
        if( i < n )
            ++sendSizes[ x.RowOwner(i) ];
        else
            ++sendSizes[ x.RowOwner(i-n) ];
    }
    vector<int> recvSizes(commSize);
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, comm );
    vector<int> sendOffs, recvOffs;
    const int totalSend = Scan( sendSizes, sendOffs );
    const int totalRecv = Scan( recvSizes, recvOffs );
    // Pack the data 
    // -------------
    vector<Int> sSendBuf(totalSend);
    vector<Real> vSendBuf(totalSend);
    auto offs = sendOffs;
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
    {
        const Int i = xHat.GlobalRow(iLoc);
        if( i < n )
        {
            const int owner = x.RowOwner(i);
            sSendBuf[offs[owner]] = i;
            vSendBuf[offs[owner]] = xHat.GetLocal(iLoc,0);
            ++offs[owner];
        }
        else
        {
            const int owner = x.RowOwner(i-n);
            sSendBuf[offs[owner]] = i-n;
            vSendBuf[offs[owner]] = -xHat.GetLocal(iLoc,0);
            ++offs[owner];
        }
    }
    // Exchange the data
    // -----------------
    vector<Int> sRecvBuf(totalRecv);
    vector<Real> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendSizes.data(), sendOffs.data(),
      sRecvBuf.data(), recvSizes.data(), recvOffs.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendSizes.data(), sendOffs.data(),
      vRecvBuf.data(), recvSizes.data(), recvOffs.data(), comm );
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
