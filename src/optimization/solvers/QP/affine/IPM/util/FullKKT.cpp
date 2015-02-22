/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "../../../direct/IPM/util.hpp"

namespace El {
namespace qp {
namespace affine {

// The full KKT system is of the form
//
//   | Q A^T      G     | | x |   |        -c             |
//   | A 0        0     | | y |   |         b             |,
//   | G 0    -(z <> s) | | z | = | -z <> (s o z + tau e) |
//
// and the particular system solved is of the form
//
//   | Q A^T      G     | | dx |   |     -rc        |
//   | A 0        0     | | dy |   |     -rb        |,
//   | G 0    -(z <> s) | | dz | = | -rh + z <> rmu |
//
// where 
//
//   rc  = Q x + A^T y + G^T z + c,
//   rb  = A x - b,
//   rh  = G x + s - h,
//   rmu = s o z - tau e

template<typename Real>
void KKT
( const Matrix<Real>& Q,
  const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& s, const Matrix<Real>& z,
        Matrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Zeros( J, n+m+k, n+m+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd); 
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd); 
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd); 

    // Jxx := Q
    // ========
    Jxx = Q;

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := G
    // ========
    Jzx = G;

    // Jzz := - z <> s
    // ===============
    Matrix<Real> t;
    t = s;
    DiagonalSolve( LEFT, NORMAL, z, t );
    Scale( Real(-1), t );
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy ); 

        // Jxz := G^T
        // ==========
        Transpose( G, Jxz );
    }
}

template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& Q,
  const AbstractDistMatrix<Real>& A,    const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& sPre, const AbstractDistMatrix<Real>& zPre,
        AbstractDistMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    auto sPtr = ReadProxy<Real,STAR,STAR>(&sPre); auto& s = *sPtr;
    auto zPtr = ReadProxy<Real,STAR,STAR>(&zPre); auto& z = *zPtr;
    auto JPtr = WriteProxy<Real,MC,MR>(&JPre);    auto& J = *JPtr;

    Zeros( J, n+m+k, n+m+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd);
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd);

    // Jxx := Q
    // ========
    Jxx = Q;

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := G
    // ========
    Jzx = G;

    // Jzz := - z <> s
    // ===============
    DistMatrix<Real,MC,STAR> t(s.Grid());
    t = s;
    DiagonalSolve( LEFT, NORMAL, z, t );
    Scale( Real(-1), t );
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy );

        // Jxz := G
        // ========
        Transpose( G, Jxz );
    }
}

template<typename Real>
void KKT
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& s,       const Matrix<Real>& z,
        SparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Zeros( J, n+m+k, n+m+k );
    const Int numEntriesQ = Q.NumEntries();
    const Int numEntriesA = A.NumEntries();
    const Int numEntriesG = G.NumEntries();
    // Count the number of entries of Q that we'll use 
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
        J.Reserve( numUsedEntriesQ + numEntriesA + numEntriesG + k );
    else
        J.Reserve( numUsedEntriesQ + 2*numEntriesA + 2*numEntriesG + k );

    // Jxx = Q
    // =======
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            J.QueueUpdate( i, j, Q.Value(e) );
    }

    // Jyx = A
    // =======
    for( Int e=0; e<numEntriesA; ++e )
        J.QueueUpdate( n+A.Row(e), A.Col(e), A.Value(e) );

    // Jzx = G
    // =======
    for( Int e=0; e<numEntriesG; ++e )
        J.QueueUpdate( n+m+G.Row(e), G.Col(e), G.Value(e) );

    // Jzz = -z <> s
    // =============
    for( Int e=0; e<k; ++e )
        J.QueueUpdate( n+m+e, n+m+e, -s.Get(e,0)/z.Get(e,0) );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        for( Int e=0; e<numEntriesA; ++e )
            J.QueueUpdate( A.Col(e), n+A.Row(e), A.Value(e) );

        // Jxz := G^T
        // ==========
        for( Int e=0; e<numEntriesG; ++e )
            J.QueueUpdate( G.Col(e), n+m+G.Row(e), G.Value(e) );
    }
    J.MakeConsistent();
}

template<typename Real>
void KKT
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& s,     const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    mpi::Comm comm = A.Comm();
    const Int commSize = mpi::Size( comm );

    J.SetComm( comm );
    Zeros( J, n+m+k, n+m+k );

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
            ++sendCounts[ J.RowOwner(i) ];
    }
    // Jyx := A
    // --------
    for( Int e=0; e<A.NumLocalEntries(); ++e )
        ++sendCounts[ J.RowOwner(A.Row(e)+n) ];
    // Jzx := G
    // --------
    for( Int e=0; e<G.NumLocalEntries(); ++e )
        ++sendCounts[ J.RowOwner(G.Row(e)+n+m) ];
    // Jxy := A^T
    // ----------
    DistSparseMatrix<Real> ATrans(comm);
    if( !onlyLower )
    {
        Transpose( A, ATrans );
        for( Int e=0; e<ATrans.NumLocalEntries(); ++e )
            ++sendCounts[ J.RowOwner(ATrans.Row(e)) ];
    }
    // Jxz := G^T
    // ----------
    DistSparseMatrix<Real> GTrans(comm);
    if( !onlyLower )
    {
        Transpose( G, GTrans );
        for( Int e=0; e<GTrans.NumLocalEntries(); ++e )
            ++sendCounts[ J.RowOwner(GTrans.Row(e)) ];
    }
    // Jzz := -z <> s
    // --------------
    for( Int iLoc=0; iLoc<s.LocalHeight(); ++iLoc )
        ++sendCounts[ J.RowOwner( m+n + s.GlobalRow(iLoc) ) ];
    // Communicate to determine the number we receive from each process
    // ----------------------------------------------------------------
    vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    vector<int> sendOffsets, recvOffsets;
    int totalSend = Scan( sendCounts, sendOffsets );
    int totalRecv = Scan( recvCounts, recvOffsets );

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
    // Pack G
    // ------
    for( Int e=0; e<G.NumLocalEntries(); ++e )
    {
        const Int i = G.Row(e) + n + m;
        const Int j = G.Col(e);
        const Real value = G.Value(e);
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
    // Pack G^T
    // --------
    if( !onlyLower )
    {
        for( Int e=0; e<GTrans.NumLocalEntries(); ++e )
        {
            const Int i = GTrans.Row(e);
            const Int j = GTrans.Col(e) + n + m;
            const Real value = GTrans.Value(e);
            const int owner = J.RowOwner(i);
            sSendBuf[offsets[owner]] = i;
            tSendBuf[offsets[owner]] = j;
            vSendBuf[offsets[owner]] = value;
            ++offsets[owner];
        }
    }
    // Pack -z <> s
    // ------------
    for( Int iLoc=0; iLoc<s.LocalHeight(); ++iLoc )
    {
        const Int i = m+n + s.GlobalRow(iLoc);
        const Int j = i;
        const Real value = -s.GetLocal(iLoc,0)/z.GetLocal(iLoc,0);
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
    J.Reserve( totalRecv );
    for( Int e=0; e<totalRecv; ++e )
        J.QueueLocalUpdate
        ( sRecvBuf[e]-J.FirstLocalRow(), tRecvBuf[e], vRecvBuf[e] );
    J.MakeConsistent();
}

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc, const Matrix<Real>& rb, 
  const Matrix<Real>& rh, const Matrix<Real>& rmu, 
  const Matrix<Real>& z, 
        Matrix<Real>& d )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::KKTRHS"))
    const Int n = rc.Height();
    const Int m = rb.Height();
    const Int k = rh.Height();
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    Zeros( d, n+m+k, 1 );

    auto dx = d(xInd,IR(0,1));
    dx = rc;
    Scale( Real(-1), dx );

    auto dy = d(yInd,IR(0,1));
    dy = rb;
    Scale( Real(-1), dy );

    auto dz = d(zInd,IR(0,1));
    dz = rmu;
    DiagonalSolve( LEFT, NORMAL, z, dz );
    Axpy( Real(-1), rh, dz );
}

template<typename Real>
void KKTRHS
( const AbstractDistMatrix<Real>& rc,  const AbstractDistMatrix<Real>& rb, 
  const AbstractDistMatrix<Real>& rh,  const AbstractDistMatrix<Real>& rmu, 
  const AbstractDistMatrix<Real>& z, 
        AbstractDistMatrix<Real>& dPre )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::KKTRHS"))

    auto dPtr = WriteProxy<Real,MC,MR>(&dPre);
    auto& d = *dPtr;

    const Int n = rc.Height();
    const Int m = rb.Height();
    const Int k = rh.Height();
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    Zeros( d, n+m+k, 1 );

    auto dx = d(xInd,IR(0,1));
    Copy( rc, dx );
    Scale( Real(-1), dx );

    auto dy = d(yInd,IR(0,1));
    Copy( rb, dy );
    Scale( Real(-1), dy );

    auto dz = d(zInd,IR(0,1));
    Copy( rmu, dz );
    DiagonalSolve( LEFT, NORMAL, z, dz );
    Axpy( Real(-1), rh, dz );
}

template<typename Real>
void KKTRHS
( const DistMultiVec<Real>& rc, const DistMultiVec<Real>& rb, 
  const DistMultiVec<Real>& rh, const DistMultiVec<Real>& rmu, 
  const DistMultiVec<Real>& z, 
        DistMultiVec<Real>& d )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::KKTRHS"))
    const Int n = rc.Height();
    const Int m = rb.Height();
    const int k = rh.Height();
    Zeros( d, n+m+k, 1 );
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
        const Real value = rmu.GetLocal(iLoc,0)/z.GetLocal(iLoc,0) - 
                           rh.GetLocal(iLoc,0);
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
void ExpandCoreSolution
( Int m, Int n, Int k,
  const Matrix<Real>& d,
        Matrix<Real>& dx, Matrix<Real>& dy,
        Matrix<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::ExpandCoreSolution"))
    if( d.Height() != n+m+k || d.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    dx = d(xInd,IR(0,1));
    dy = d(yInd,IR(0,1));
    dz = d(zInd,IR(0,1));
}

template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const AbstractDistMatrix<Real>& dPre,
        AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy,
        AbstractDistMatrix<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::ExpandCoreSolution"))

    auto dPtr = ReadProxy<Real,MC,MR>(&dPre);
    auto& d = *dPtr;

    if( d.Height() != n+m+k || d.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    Copy( d(xInd,IR(0,1)), dx );
    Copy( d(yInd,IR(0,1)), dy );
    Copy( d(zInd,IR(0,1)), dz );
}

template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const DistMultiVec<Real>& d,
        DistMultiVec<Real>& dx, DistMultiVec<Real>& dy,
        DistMultiVec<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::ExpandCoreSolution"))
    if( d.Height() != n+m+k || d.Width() != 1 )
        LogicError("Right-hand side was the wrong size");
    mpi::Comm comm = d.Comm();
    const Int commSize = mpi::Size(comm);

    dx.SetComm( comm );
    dy.SetComm( comm );
    dz.SetComm( comm );
    dx.Resize( n, 1 );
    dy.Resize( m, 1 );
    dz.Resize( k, 1 );

    // Compute the metadata for the AllToAll
    // =====================================
    vector<int> sendCounts(commSize,0);
    for( Int iLoc=0; iLoc<d.LocalHeight(); ++iLoc )
    {
        const Int i = d.GlobalRow(iLoc);
        if( i < n )
            ++sendCounts[ dx.RowOwner(i) ];
        else if( i < n+m )
            ++sendCounts[ dy.RowOwner(i-n) ];
        else
            ++sendCounts[ dz.RowOwner(i-(n+m)) ];
    }
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
    for( Int iLoc=0; iLoc<d.LocalHeight(); ++iLoc )
    {
        const Int i = d.GlobalRow(iLoc);
        int owner;
        if( i < n )
            owner = dx.RowOwner(i);
        else if( i < n+m )
            owner = dy.RowOwner(i-n);
        else
            owner = dz.RowOwner(i-(n+m));
        sSendBuf[offsets[owner]] = i;
        vSendBuf[offsets[owner]] = d.GetLocal(iLoc,0);
        ++offsets[owner];
    }

    // Exchange the doublets
    // =====================
    vector<Int> sRecvBuf(totalRecv);
    vector<Real> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );

    // Unpack the doublets
    // ===================
    for( Int e=0; e<totalRecv; ++e )
    {
        const Int i = sRecvBuf[e];
        if( i < n )
            dx.SetLocal( i-dx.FirstLocalRow(), 0, vRecvBuf[e] );
        else if( i < n+m )
            dy.SetLocal( i-n-dy.FirstLocalRow(), 0, vRecvBuf[e] );
        else
            dz.SetLocal( i-(n+m)-dz.FirstLocalRow(), 0, vRecvBuf[e] );
    }
}

template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const Matrix<Real>& d,  const Matrix<Real>& rmu,
  const Matrix<Real>& s,  const Matrix<Real>& z, 
        Matrix<Real>& dx,       Matrix<Real>& dy, 
        Matrix<Real>& dz,       Matrix<Real>& ds )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::ExpandSolution"))
    const Int k = s.Height();
    ExpandCoreSolution( m, n, k, d, dx, dy, dz );
    // ds := - z <> ( rmu + s o dz )
    // =============================
    ds = dz;
    DiagonalScale( LEFT, NORMAL, s, ds );
    Axpy( Real(1), rmu, ds );
    DiagonalSolve( LEFT, NORMAL, z, ds );
    Scale( Real(-1), ds );
}

template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const AbstractDistMatrix<Real>& d,  const AbstractDistMatrix<Real>& rmu,
  const AbstractDistMatrix<Real>& s,  const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& dx,       AbstractDistMatrix<Real>& dy, 
        AbstractDistMatrix<Real>& dz,       AbstractDistMatrix<Real>& ds )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::ExpandSolution"))
    const int k = s.Height();
    ExpandCoreSolution( m, n, k, d, dx, dy, dz );
    // ds := - z <> ( rmu + s o dz )
    // =============================
    Copy( dz, ds );
    DiagonalScale( LEFT, NORMAL, s, ds );
    Axpy( Real(1), rmu, ds );
    DiagonalSolve( LEFT, NORMAL, z, ds );
    Scale( Real(-1), ds );
}

template<typename Real>
void ExpandSolution
( Int m, Int n,
  const DistMultiVec<Real>& d, const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& s, const DistMultiVec<Real>& z,
        DistMultiVec<Real>& dx,      DistMultiVec<Real>& dy, 
        DistMultiVec<Real>& dz,      DistMultiVec<Real>& ds )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::ExpandSolution"))
    const Int k = s.Height();
    ExpandCoreSolution( m, n, k, d, dx, dy, dz );
    // ds := - z <> ( rmu + s o dz )
    // =============================
    ds = dz;
    DiagonalScale( LEFT, NORMAL, s, ds );
    Axpy( Real(1), rmu, ds );
    DiagonalSolve( LEFT, NORMAL, z, ds );
    Scale( Real(-1), ds );
}

#define PROTO(Real) \
  template void KKT \
  ( const Matrix<Real>& Q, \
    const Matrix<Real>& A, const Matrix<Real>& G, \
    const Matrix<Real>& s, const Matrix<Real>& z, \
          Matrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const AbstractDistMatrix<Real>& Q, \
    const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const SparseMatrix<Real>& Q, \
    const SparseMatrix<Real>& A, const SparseMatrix<Real>& G, \
    const Matrix<Real>& s,       const Matrix<Real>& z, \
          SparseMatrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const DistSparseMatrix<Real>& Q, \
    const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& s,     const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J, bool onlyLower ); \
  template void KKTRHS \
  ( const Matrix<Real>& rc, const Matrix<Real>& rb, \
    const Matrix<Real>& rh, const Matrix<Real>& rmu, \
    const Matrix<Real>& z, \
          Matrix<Real>& d ); \
  template void KKTRHS \
  ( const AbstractDistMatrix<Real>& rc, const AbstractDistMatrix<Real>& rb, \
    const AbstractDistMatrix<Real>& rh, const AbstractDistMatrix<Real>& rmu, \
    const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& d ); \
  template void KKTRHS \
  ( const DistMultiVec<Real>& rc, const DistMultiVec<Real>& rb, \
    const DistMultiVec<Real>& rh, const DistMultiVec<Real>& rmu, \
    const DistMultiVec<Real>& z, \
          DistMultiVec<Real>& d ); \
  template void ExpandCoreSolution \
  ( Int m, Int n, Int k, \
    const Matrix<Real>& d, \
          Matrix<Real>& dx, Matrix<Real>& dy, \
          Matrix<Real>& dz ); \
  template void ExpandCoreSolution \
  ( Int m, Int n, Int k, \
    const AbstractDistMatrix<Real>& d, \
          AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, \
          AbstractDistMatrix<Real>& dz ); \
  template void ExpandCoreSolution \
  ( Int m, Int n, Int k, \
    const DistMultiVec<Real>& d, \
          DistMultiVec<Real>& dx, DistMultiVec<Real>& dy, \
          DistMultiVec<Real>& dz ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const Matrix<Real>& d,  const Matrix<Real>& rmu, \
    const Matrix<Real>& s,  const Matrix<Real>& z, \
          Matrix<Real>& dx,       Matrix<Real>& dy, \
          Matrix<Real>& dz,       Matrix<Real>& ds ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const AbstractDistMatrix<Real>& d,  const AbstractDistMatrix<Real>& rmu, \
    const AbstractDistMatrix<Real>& s,  const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& dx,       AbstractDistMatrix<Real>& dy, \
          AbstractDistMatrix<Real>& dz,       AbstractDistMatrix<Real>& ds ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const DistMultiVec<Real>& d,  const DistMultiVec<Real>& rmu, \
    const DistMultiVec<Real>& s,  const DistMultiVec<Real>& z, \
          DistMultiVec<Real>& dx,       DistMultiVec<Real>& dy, \
          DistMultiVec<Real>& dz,       DistMultiVec<Real>& ds );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace affine
} // namespace qp
} // namespace El
