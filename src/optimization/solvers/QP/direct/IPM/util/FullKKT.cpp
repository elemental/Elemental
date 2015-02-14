/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "../../../affine/IPM/util.hpp"

namespace El {
namespace qp {
namespace direct {

// The full KKT system is of the form
//
//   |  Q A^T     -I    | | x |   |          -c            |
//   |  A 0        0    | | y |   |           b            |,
//   | -I 0   (-z <> x) | | z | = | - z <> (x o z + tau e) |
//
// and the particular system solved is of the form
//
//   |  Q A^T     -I    | | dx |   |   -rc    |
//   |  A 0        0    | | dy |   |   -rb    |,
//   | -I 0   (-z <> x) | | dz | = | z <> rmu |
//
// where 
//
//   rc = Q x + A^T y - z + c,
//   rb = A x - b,
//   rmu = x o z - tau e

template<typename Real>
void KKT
( const Matrix<Real>& Q, const Matrix<Real>& A, 
  const Matrix<Real>& x, const Matrix<Real>& z,
        Matrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    Zeros( J, 2*n+m, 2*n+m );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,2*n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd); 
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd); 
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd); 

    // Jxx := Q
    // ========
    Jxx = Q;

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := -I
    // =========
    Identity( Jzx, n, n );
    Scale( Real(-1), Jzx );

    // Jzz := - z <> x
    // ===============
    Matrix<Real> t;
    t = x;
    DiagonalSolve( LEFT, NORMAL, z, t );
    Scale( Real(-1), t );
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy ); 

        // Jxz := -I
        // =========
        Identity( Jxz, n, n );
        Scale( Real(-1), Jxz );
    }
}

template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& Q,    const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& xPre, const AbstractDistMatrix<Real>& zPre,
        AbstractDistMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    auto xPtr = ReadProxy<Real,STAR,STAR>(&xPre); auto& x = *xPtr;
    auto zPtr = ReadProxy<Real,STAR,STAR>(&zPre); auto& z = *zPtr;
    auto JPtr = WriteProxy<Real,MC,MR>(&JPre);    auto& J = *JPtr;

    Zeros( J, 2*n+m, 2*n+m );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,2*n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd);
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd);

    // Jxx := Q
    // ========
    Jxx = Q;

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := -I
    // =========
    Identity( Jzx, n, n );
    Scale( Real(-1), Jzx );

    // Jzz := - z <> x
    // ===============
    DistMatrix<Real,MC,STAR> t(x.Grid());
    t = x;
    DiagonalSolve( LEFT, NORMAL, z, t );
    Scale( Real(-1), t );
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy );

        // Jxz := -I
        // =========
        Identity( Jxz, n, n );
        Scale( Real(-1), Jxz );
    }
}

template<typename Real>
void KKT
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, 
  const Matrix<Real>& x,       const Matrix<Real>& z,
        SparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    Zeros( J, 2*n+m, 2*n+m );
    const Int numEntriesQ = Q.NumEntries();
    const Int numEntriesA = A.NumEntries();
    // Count the number of used entries of Q
    Int numUsedEntriesQ;
    if( onlyLower )
    {
        numUsedEntriesQ = 0;
        for( Int e=0; e<numEntriesQ; ++e )
            if( Q.Row(e) >= Q.Col(e) )
                ++numUsedEntriesQ;
    }
    else
        numUsedEntriesQ = numEntriesQ;

    if( onlyLower )
        J.Reserve( numUsedEntriesQ + numEntriesA + 2*n );
    else
        J.Reserve( numUsedEntriesQ + 2*numEntriesA + 3*n );

    // Jxx = Q
    // =======
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            J.QueueUpdate( Q.Row(e), Q.Col(e), Q.Value(e) );
    }

    // Jyx = A
    // =======
    for( Int e=0; e<numEntriesA; ++e )
        J.QueueUpdate( n+A.Row(e), A.Col(e), A.Value(e) );

    // Jzx = -I
    // ========
    for( Int e=0; e<n; ++e )
        J.QueueUpdate( n+m+e, e, Real(-1) );

    // Jzz = - z <> x
    // ==============
    for( Int e=0; e<n; ++e )
        J.QueueUpdate( n+m+e, n+m+e, -x.Get(e,0)/z.Get(e,0) );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        for( Int e=0; e<numEntriesA; ++e )
            J.QueueUpdate( A.Col(e), n+A.Row(e), A.Value(e) );

        // Jxz := -I
        // =========
        for( Int e=0; e<n; ++e )
            J.QueueUpdate( e, n+m+e, Real(-1) );
    }
    J.MakeConsistent();
}

template<typename Real>
void KKT
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    mpi::Comm comm = A.Comm();
    const Int commSize = mpi::Size( comm );

    J.SetComm( comm );
    Zeros( J, m+2*n, m+2*n );

    // Compute the number of entries to send to each process
    // =====================================================
    vector<int> sendCounts(commSize,0);
    // Jxx := Q
    // --------
    for( Int e=0; e<Q.NumLocalEntries(); ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            ++sendCounts[ J.RowOwner(Q.Row(e)) ];
    }
    // Jyx := A
    // --------
    for( Int e=0; e<A.NumLocalEntries(); ++e )
        ++sendCounts[ J.RowOwner(A.Row(e)+n) ];
    // Jxy := A^T
    // ----------
    DistSparseMatrix<Real> ATrans(comm);
    if( !onlyLower )
    {
        Transpose( A, ATrans );
        for( Int e=0; e<ATrans.NumLocalEntries(); ++e )
            ++sendCounts[ J.RowOwner(ATrans.Row(e)) ];
    }
    // Jzz := -z <> x
    // --------------
    for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
        ++sendCounts[ J.RowOwner( m+n + x.GlobalRow(iLoc) ) ];
    // Communicate to determine the number we receive from each process
    // ----------------------------------------------------------------
    vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    vector<int> sendOffsets, recvOffsets;
    const int totalSend = Scan( sendCounts, sendOffsets );
    const int totalRecv = Scan( recvCounts, recvOffsets );

    // Pack the triplets
    // =================
    vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
    vector<Real> vSendBuf(totalSend);
    auto offsets = sendOffsets;
    // Pack Q
    // ------
    for( Int e=0; e<Q.NumLocalEntries(); ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower ) 
        {
            const Real value = Q.Value(e);
            const int owner = J.RowOwner(i);
            sSendBuf[offsets[owner]] = i; 
            tSendBuf[offsets[owner]] = j;
            vSendBuf[offsets[owner]] = value;
            ++offsets[owner];
        }
    }
    // Pack A
    // ------
    for( Int e=0; e<A.NumLocalEntries(); ++e )
    {
        const Int i = A.Row(e) + n;
        const Int j = A.Col(e);
        const Real value = A.Value(e);
        const int owner = J.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        tSendBuf[offsets[owner]] = j;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }
    // Pack A^T
    // --------
    if( !onlyLower )
    {
        for( Int e=0; e<ATrans.NumLocalEntries(); ++e )
        {
            const Int i = ATrans.Row(e);
            const Int j = ATrans.Col(e) + n;
            const Real value = ATrans.Value(e);
            const int owner = J.RowOwner(i);
            sSendBuf[offsets[owner]] = i;
            tSendBuf[offsets[owner]] = j;
            vSendBuf[offsets[owner]] = value;
            ++offsets[owner];
        }
    }
    // Pack -z <> x
    // ------------
    for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
    {
        const Int i = m+n + x.GlobalRow(iLoc);
        const Int j = i;
        const Real value = -x.GetLocal(iLoc,0)/z.GetLocal(iLoc,0);
        const int owner = J.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        tSendBuf[offsets[owner]] = j;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }

    // Exchange the triplets
    // =====================
    vector<Int> sRecvBuf(totalRecv), tRecvBuf(totalRecv);
    vector<Real> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( tSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      tRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );

    // Unpack the triplets
    // ===================
    // Count the total number of entries for the negative identities
    // -------------------------------------------------------------
    Int negIdentUpdates = 0;
    for( Int iLoc=0; iLoc<J.LocalHeight(); ++iLoc )
    {
        const Int i = J.GlobalRow(iLoc);
        if( i < n && !onlyLower )
            ++negIdentUpdates;
        else if( i >= n+m )
            ++negIdentUpdates;
    }
    // Reserve the total number of local updates
    // -----------------------------------------
    J.Reserve( totalRecv+negIdentUpdates );
    // Append the local negative identity updates
    // ------------------------------------------
    for( Int iLoc=0; iLoc<J.LocalHeight(); ++iLoc )
    {
        const Int i = J.GlobalRow(iLoc);
        if( i < n && !onlyLower )
            J.QueueLocalUpdate( iLoc, i+(n+m), Real(-1) );
        else if( i >= n+m )
            J.QueueLocalUpdate( iLoc, i-(n+m), Real(-1) );
    }
    // Append the received local updates
    // ---------------------------------
    for( Int e=0; e<totalRecv; ++e )
        J.QueueLocalUpdate
        ( sRecvBuf[e]-J.FirstLocalRow(), tRecvBuf[e], vRecvBuf[e] );
    J.MakeConsistent();
}

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc,  const Matrix<Real>& rb, 
  const Matrix<Real>& rmu, const Matrix<Real>& z, 
        Matrix<Real>& d )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::KKTRHS"))
    const Int m = rb.Height();
    const Int n = rc.Height();
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,2*n+m);
    Zeros( d, 2*n+m, 1 );

    auto dx = d(xInd,IR(0,1));
    dx = rc;
    Scale( Real(-1), dx );

    auto dy = d(yInd,IR(0,1));
    dy = rb;
    Scale( Real(-1), dy );

    auto dz = d(zInd,IR(0,1));
    dz = rmu;
    DiagonalSolve( LEFT, NORMAL, z, dz );
}

template<typename Real>
void KKTRHS
( const AbstractDistMatrix<Real>& rc,  const AbstractDistMatrix<Real>& rb, 
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& z, 
        AbstractDistMatrix<Real>& dPre )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::KKTRHS"))

    auto dPtr = WriteProxy<Real,MC,MR>(&dPre); 
    auto& d = *dPtr;

    const Int m = rb.Height();
    const Int n = rc.Height();
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,2*n+m);
    Zeros( d, 2*n+m, 1 );

    auto dx = d(xInd,IR(0,1));
    Copy( rc, dx );
    Scale( Real(-1), dx );

    auto dy = d(yInd,IR(0,1));
    Copy( rb, dy );
    Scale( Real(-1), dy );

    auto dz = d(zInd,IR(0,1));
    Copy( rmu, dz );
    DiagonalSolve( LEFT, NORMAL, z, dz );
}

template<typename Real>
void KKTRHS
( const DistMultiVec<Real>& rc,  const DistMultiVec<Real>& rb, 
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& z, 
        DistMultiVec<Real>& d )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::KKTRHS"))
    const Int m = rb.Height();
    const Int n = rc.Height();
    Zeros( d, m+2*n, 1 );
    mpi::Comm comm = rmu.Comm();
    const Int commSize = mpi::Size( comm ); 

    // Compute the number of entries to send/recv from each process
    // ============================================================
    vector<int> sendCounts(commSize,0);
    for( Int iLoc=0; iLoc<rc.LocalHeight(); ++iLoc )
        ++sendCounts[ d.RowOwner( rc.GlobalRow(iLoc) ) ];
    for( Int iLoc=0; iLoc<rb.LocalHeight(); ++iLoc )
        ++sendCounts[ d.RowOwner( n + rb.GlobalRow(iLoc) ) ];
    for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
        ++sendCounts[ d.RowOwner( n+m + rmu.GlobalRow(iLoc) ) ];
    vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    vector<int> sendOffsets, recvOffsets;
    const int totalSend = Scan( sendCounts, sendOffsets );
    const int totalRecv = Scan( recvCounts, recvOffsets );

    // Pack the doublets
    // =================
    vector<Int> sSendBuf(totalSend);
    vector<Real> vSendBuf(totalSend);
    auto offsets = sendOffsets;
    for( Int iLoc=0; iLoc<rc.LocalHeight(); ++iLoc )
    {
        const Int i = rc.GlobalRow(iLoc);
        const Real value = -rc.GetLocal(iLoc,0);
        const int owner = d.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }
    for( Int iLoc=0; iLoc<rb.LocalHeight(); ++iLoc )
    {
        const Int i = n + rb.GlobalRow(iLoc);
        const Real value = -rb.GetLocal(iLoc,0);
        const int owner = d.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }
    for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
    {
        const Int i = n+m + rmu.GlobalRow(iLoc);
        const Real value = rmu.GetLocal(iLoc,0)/z.GetLocal(iLoc,0);
        const int owner = d.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }

    // Exchange and unpack the doublets
    // ================================
    vector<Int> sRecvBuf(totalRecv);
    vector<Real> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    for( Int e=0; e<totalRecv; ++e )
        d.UpdateLocal( sRecvBuf[e]-d.FirstLocalRow(), 0, vRecvBuf[e] );
}

template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const Matrix<Real>& d, 
        Matrix<Real>& dx, Matrix<Real>& dy, 
        Matrix<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::ExpandSolution"))
    qp::affine::ExpandCoreSolution( m, n, n, d, dx, dy, dz );
}

template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const AbstractDistMatrix<Real>& d, 
        AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, 
        AbstractDistMatrix<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::ExpandSolution"))
    qp::affine::ExpandCoreSolution( m, n, n, d, dx, dy, dz );
}

template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const DistMultiVec<Real>& d, 
        DistMultiVec<Real>& dx, DistMultiVec<Real>& dy, 
        DistMultiVec<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::ExpandSolution"))
    qp::affine::ExpandCoreSolution( m, n, n, d, dx, dy, dz );
}

#define PROTO(Real) \
  template void KKT \
  ( const Matrix<Real>& Q, const Matrix<Real>& A, \
    const Matrix<Real>& x, const Matrix<Real>& z, \
          Matrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, \
    const Matrix<Real>& x,       const Matrix<Real>& z, \
          SparseMatrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J, bool onlyLower ); \
  template void KKTRHS \
  ( const Matrix<Real>& rc,  const Matrix<Real>& rb, \
    const Matrix<Real>& rmu, const Matrix<Real>& z, \
          Matrix<Real>& d ); \
  template void KKTRHS \
  ( const AbstractDistMatrix<Real>& rc,  const AbstractDistMatrix<Real>& rb, \
    const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& d ); \
  template void KKTRHS \
  ( const DistMultiVec<Real>& rc,  const DistMultiVec<Real>& rb, \
    const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& z, \
          DistMultiVec<Real>& d ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const Matrix<Real>& d, \
          Matrix<Real>& dx, Matrix<Real>& dy, \
          Matrix<Real>& dz ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const AbstractDistMatrix<Real>& d, \
          AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, \
          AbstractDistMatrix<Real>& dz ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const DistMultiVec<Real>& d, \
          DistMultiVec<Real>& dx, DistMultiVec<Real>& dy, \
          DistMultiVec<Real>& dz );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace direct
} // namespace qp
} // namespace El
