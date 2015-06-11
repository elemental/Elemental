/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

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
    DEBUG_ONLY(CSE cse("qp::affine::KKT"))
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
  const AbstractDistMatrix<Real>& s,    const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("qp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

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
    DistMatrix<Real,MC,STAR> t(s);
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
    DEBUG_ONLY(CSE cse("qp::affine::KKT"))
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
    J.ProcessQueues();
}

template<typename Real>
void KKT
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& s,     const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("qp::affine::KKT"))
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

    // Pack the triplets
    // =================
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    vector<Entry<Real>> sendBuf(totalSend);
    auto offs = sendOffs;
    // Pack Q
    // ------
    for( Int e=0; e<Q.NumLocalEntries(); ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            sendBuf[offs[J.RowOwner(i)]++] = Entry<Real>{ i, j, Q.Value(e) };
    }
    // Pack A
    // ------
    for( Int e=0; e<A.NumLocalEntries(); ++e )
    {
        const Int i = A.Row(e) + n;
        const Int j = A.Col(e);
        sendBuf[offs[J.RowOwner(i)]++] = Entry<Real>{ i, j, A.Value(e) };
    }
    // Pack G
    // ------
    for( Int e=0; e<G.NumLocalEntries(); ++e )
    {
        const Int i = G.Row(e) + n + m;
        const Int j = G.Col(e);
        sendBuf[offs[J.RowOwner(i)]++] = Entry<Real>{ i, j, G.Value(e) };
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
            sendBuf[offs[J.RowOwner(i)]++] = Entry<Real>{ i, j, value };
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
            sendBuf[offs[J.RowOwner(i)]++] = Entry<Real>{ i, j, value };
        }
    }
    // Pack -z <> s
    // ------------
    for( Int iLoc=0; iLoc<s.LocalHeight(); ++iLoc )
    {
        const Int i = m+n + s.GlobalRow(iLoc);
        const Int j = i;
        const Real value = -s.GetLocal(iLoc,0)/z.GetLocal(iLoc,0);
        sendBuf[offs[J.RowOwner(i)]++] = Entry<Real>{ i, j, value };
    }

    // Exchange and unpack the triplets
    // ================================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    J.Reserve( recvBuf.size() );
    for( auto& entry : recvBuf )
        J.QueueUpdate( entry );
    J.ProcessQueues();
}

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc, const Matrix<Real>& rb, 
  const Matrix<Real>& rh, const Matrix<Real>& rmu, 
  const Matrix<Real>& z, 
        Matrix<Real>& d )
{
    DEBUG_ONLY(CSE cse("qp::affine::KKTRHS"))
    const Int n = rc.Height();
    const Int m = rb.Height();
    const Int k = rh.Height();
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    Zeros( d, n+m+k, 1 );

    auto dx = d(xInd,ALL);
    dx = rc;
    Scale( Real(-1), dx );

    auto dy = d(yInd,ALL);
    dy = rb;
    Scale( Real(-1), dy );

    auto dz = d(zInd,ALL);
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
    DEBUG_ONLY(CSE cse("qp::affine::KKTRHS"))

    auto dPtr = WriteProxy<Real,MC,MR>(&dPre);
    auto& d = *dPtr;

    const Int n = rc.Height();
    const Int m = rb.Height();
    const Int k = rh.Height();
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    Zeros( d, n+m+k, 1 );

    auto dx = d(xInd,ALL);
    Copy( rc, dx );
    Scale( Real(-1), dx );

    auto dy = d(yInd,ALL);
    Copy( rb, dy );
    Scale( Real(-1), dy );

    auto dz = d(zInd,ALL);
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
    DEBUG_ONLY(CSE cse("qp::affine::KKTRHS"))
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

    // Pack the doublets
    // =================
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    vector<ValueInt<Real>> sendBuf(totalSend);
    auto offs = sendOffs;
    for( Int iLoc=0; iLoc<rc.LocalHeight(); ++iLoc )
    {
        const Int i = rc.GlobalRow(iLoc);
        const Real value = -rc.GetLocal(iLoc,0);
        sendBuf[offs[d.RowOwner(i)]++] = ValueInt<Real>{ value, i };
    }
    for( Int iLoc=0; iLoc<rb.LocalHeight(); ++iLoc )
    {
        const Int i = n + rb.GlobalRow(iLoc);
        const Real value = -rb.GetLocal(iLoc,0);
        sendBuf[offs[d.RowOwner(i)]++] = ValueInt<Real>{ value, i };
    }
    for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
    {
        const Int i = n+m + rmu.GlobalRow(iLoc);
        const Real value = rmu.GetLocal(iLoc,0)/z.GetLocal(iLoc,0) - 
                           rh.GetLocal(iLoc,0);
        sendBuf[offs[d.RowOwner(i)]++] = ValueInt<Real>{ value, i };
    }

    // Exchange and unpack the doublets
    // ================================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    for( auto& entry : recvBuf )
        d.Update( entry.index, 0, entry.value );
}

template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const Matrix<Real>& d,
        Matrix<Real>& dx, Matrix<Real>& dy,
        Matrix<Real>& dz )
{
    DEBUG_ONLY(CSE cse("qp::affine::ExpandCoreSolution"))
    if( d.Height() != n+m+k || d.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    dx = d(xInd,ALL);
    dy = d(yInd,ALL);
    dz = d(zInd,ALL);
}

template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const AbstractDistMatrix<Real>& dPre,
        AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy,
        AbstractDistMatrix<Real>& dz )
{
    DEBUG_ONLY(CSE cse("qp::affine::ExpandCoreSolution"))

    auto dPtr = ReadProxy<Real,MC,MR>(&dPre);
    auto& d = *dPtr;

    if( d.Height() != n+m+k || d.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    Copy( d(xInd,ALL), dx );
    Copy( d(yInd,ALL), dy );
    Copy( d(zInd,ALL), dz );
}

template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const DistMultiVec<Real>& d,
        DistMultiVec<Real>& dx, DistMultiVec<Real>& dy,
        DistMultiVec<Real>& dz )
{
    DEBUG_ONLY(CSE cse("qp::affine::ExpandCoreSolution"))
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

    // Pack the doublets
    // =================
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    vector<ValueInt<Real>> sendBuf(totalSend);
    auto offs = sendOffs;
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
        sendBuf[offs[owner]++] = ValueInt<Real>{ d.GetLocal(iLoc,0), i };
    }

    // Exchange and unpack the doublets
    // ================================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    for( auto& entry : recvBuf )
    {
        const Int i = entry.index;
        if( i < n )
            dx.Set( i, 0, entry.value );
        else if( i < n+m )
            dy.Set( i-n, 0, entry.value );
        else
            dz.Set( i-(n+m), 0, entry.value );
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
    DEBUG_ONLY(CSE cse("qp::affine::ExpandSolution"))
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
    DEBUG_ONLY(CSE cse("qp::affine::ExpandSolution"))
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
    DEBUG_ONLY(CSE cse("qp::affine::ExpandSolution"))
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
