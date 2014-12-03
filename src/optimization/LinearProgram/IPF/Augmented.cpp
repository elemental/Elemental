/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace lin_prog {

// Form 
//    J = | -S*inv(X) A^T | and y = | -r_c + r_mu / X |,
//        |  A        0   |         | -r_b            |
// where 
//    S   = diag(s),
//    X   = diag(x),
//    e   = ones(n,1),
//    r_b = A x - b,
//    r_c = A^T l + s - c, and
//    r_mu = X S e - tau e.
//
// The implied system is of the form
//   J | dx | = y,     ds = -(r_mu + S dx) / X.
//     | dl |
//

template<typename Real>
void AugmentedKKT
( const Matrix<Real>& A, 
  const Matrix<Real>& s, const Matrix<Real>& x,
  Matrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::AugmentedKKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    Zeros( J, m+n, m+n );
    const IR xInd(0,n), lInd(n,n+m);
    auto Jxx = J(xInd,xInd); auto Jxl = J(xInd,lInd);
    auto Jlx = J(lInd,xInd); auto Jll = J(lInd,lInd);
    Matrix<Real> d( s );
    Scale( Real(-1), d );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( Jxx, d );
    Jlx = A;
    if( !onlyLower )
        Transpose( A, Jxl );
}

template<typename Real>
void AugmentedKKT
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& x,
  AbstractDistMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::AugmentedKKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    auto JPtr = WriteProxy<Real,MC,MR>(&JPre);
    auto& J = *JPtr;

    Zeros( J, m+n, m+n );
    const IR xInd(0,n), lInd(n,n+m);
    auto Jxx = J(xInd,xInd); auto Jxl = J(xInd,lInd);
    auto Jlx = J(lInd,xInd); auto Jll = J(lInd,lInd);
    DistMatrix<Real,STAR,STAR> d( s );
    Scale( Real(-1), d );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( Jxx, d.Matrix() );
    Jlx = A;
    if( !onlyLower )
        Transpose( A, Jxl );
}

template<typename Real>
void AugmentedKKT
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& s, const Matrix<Real>& x,
  SparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::AugmentedKKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();

    Zeros( J, m+n, m+n );
    J.Reserve( 2*numEntries + n ); 
    // -S*inv(X) updates
    for( Int j=0; j<n; ++j )
        J.Update( j, j, -s.Get(j,0)/x.Get(j,0) );
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k);
        const Real value = A.Value(k);
        // A update
        J.Update( i+n, j, value );
        // A^T update
        if( !onlyLower )
            J.Update( j, i+n, value );
    }
    J.MakeConsistent();
}

template<typename Real>
void AugmentedKKT
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& s, const DistMultiVec<Real>& x,
  DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::AugmentedKKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    mpi::Comm comm = A.Comm();
    const Int commSize = mpi::Size( comm );

    DistSparseMatrix<Real> ATrans(comm);
    Transpose( A, ATrans );

    J.SetComm( comm );
    Zeros( J, m+n, m+n );

    // Compute the number of entries to send to each process
    // =====================================================
    std::vector<int> sendCounts(commSize,0);
    // For placing A into the bottom-left corner
    // -----------------------------------------
    for( Int k=0; k<A.NumLocalEntries(); ++k )
        ++sendCounts[ J.RowOwner(A.Row(k)+n) ];
    // For placing A^T into the top-right corner
    // -----------------------------------------
    if( !onlyLower )
        for( Int k=0; k<ATrans.NumLocalEntries(); ++k )
            ++sendCounts[ J.RowOwner(ATrans.Row(k)) ];
    // For placing -S*inv(X) into the top-left corner
    // ----------------------------------------------
    for( Int k=0; k<x.LocalHeight(); ++k )
        ++sendCounts[ J.RowOwner( k+x.FirstLocalRow() ) ];
    // Communicate to determine the number we receive from each process
    // ----------------------------------------------------------------
    std::vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );

    // Convert the send/recv counts into offsets and total sizes
    // =========================================================
    std::vector<int> sendOffsets, recvOffsets;
    int totalSend = Scan( sendCounts, sendOffsets );
    int totalRecv = Scan( recvCounts, recvOffsets );

    // Pack the triplets
    // =================
    std::vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
    std::vector<Real> vSendBuf(totalSend);
    auto offsets = sendOffsets;
    // Pack A
    // ------
    for( Int k=0; k<A.NumLocalEntries(); ++k )
    {
        const Int i = A.Row(k) + n;
        const Int j = A.Col(k);
        const Real value = A.Value(k);
        const Int owner = J.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        tSendBuf[offsets[owner]] = j;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }
    // Pack A^T
    // --------
    if( !onlyLower )
    {
        for( Int k=0; k<ATrans.NumLocalEntries(); ++k )
        {
            const Int i = ATrans.Row(k);
            const Int j = ATrans.Col(k) + n;
            const Real value = ATrans.Value(k);
            const Int owner = J.RowOwner(i);
            sSendBuf[offsets[owner]] = i;
            tSendBuf[offsets[owner]] = j;
            vSendBuf[offsets[owner]] = value;
            ++offsets[owner];
        }
    }
    // Pack -S inv(X)
    // --------------
    for( Int k=0; k<x.LocalHeight(); ++k )
    {
        const Int i = k + x.FirstLocalRow();
        const Int j = i;
        const Real value = -s.GetLocal(k,0)/x.GetLocal(k,0);
        const Int owner = J.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        tSendBuf[offsets[owner]] = j;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }

    // Exchange and unpack the triplets
    // ================================
    std::vector<Int> sRecvBuf(totalRecv), tRecvBuf(totalRecv);
    std::vector<Real> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( tSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      tRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    J.Reserve( totalRecv );
    for( Int k=0; k<totalRecv; ++k )
        J.QueueLocalUpdate
        ( sRecvBuf[k]-J.FirstLocalRow(), tRecvBuf[k], vRecvBuf[k] );
    J.MakeConsistent();
}

template<typename Real>
void AugmentedKKTRHS
( const Matrix<Real>& x, 
  const Matrix<Real>& rmu, const Matrix<Real>& rc, const Matrix<Real>& rb,
  Matrix<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::AugmentedKKTRHS"))
    const Int m = rb.Height();
    const Int n = rmu.Height();
    const IR xInd(0,n), lInd(n,n+m);
    Zeros( y, m+n, 1 );

    auto yx = y(xInd,IR(0,1));
    yx = rmu;
    for( Int i=0; i<n; ++i )
        yx.Set( i, 0, yx.Get(i,0)/x.Get(i,0) );
    Axpy( Real(-1), rc, yx );

    auto yl = y(lInd,IR(0,1));
    yl = rb;
    Scale( Real(-1), yl );
}

template<typename Real>
void AugmentedKKTRHS
( const AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc, 
  const AbstractDistMatrix<Real>& rb,
  AbstractDistMatrix<Real>& yPre )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::AugmentedKKTRHS"))

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.rowConstrain = true;
    ctrl.colAlign = 0;
    ctrl.rowAlign = 0;
    
    auto xPtr = ReadProxy<Real,MC,MR>(&xPre,ctrl);  auto& x = *xPtr;
    auto yPtr = WriteProxy<Real,MC,MR>(&yPre,ctrl); auto& y = *yPtr;

    const Int m = rb.Height();
    const Int n = rmu.Height();
    const IR xInd(0,n), lInd(n,n+m);
    Zeros( y, m+n, 1 );

    auto yx = y(xInd,IR(0,1));
    yx = rmu;
    for( Int iLoc=0; iLoc<yx.LocalHeight(); ++iLoc )
        yx.SetLocal( iLoc, 0, yx.GetLocal(iLoc,0)/x.GetLocal(iLoc,0) );
    Axpy( Real(-1), rc, yx );

    auto yl = y(lInd,IR(0,1));
    yl = rb;
    Scale( Real(-1), yl );
}

template<typename Real>
void AugmentedKKTRHS
( const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rc, 
  const DistMultiVec<Real>& rb, DistMultiVec<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::FormAugmentedSystem"))
    const Int m = rb.Height();
    const Int n = x.Height();

    mpi::Comm comm = x.Comm();
    const Int commSize = mpi::Size( comm );

    Zeros( y, m+n, 1 );
    DistMultiVec<Real> yx(comm), yl(comm);

    yx = rc;
    Scale( Real(-1), yx );
    for( Int k=0; k<x.LocalHeight(); ++k )
        yx.UpdateLocal( k, 0, rmu.GetLocal(k,0)/x.GetLocal(k,0) );    

    yl = rb;
    Scale( Real(-1), yl );

    // Compute the number of entries to send to each process
    // =====================================================
    std::vector<int> sendCounts(commSize), recvCounts(commSize);
    for( Int q=0; q<commSize; ++q )
        sendCounts[q] = 0;
    for( Int k=0; k<yx.LocalHeight(); ++k )
        ++sendCounts[ y.RowOwner( k+yx.FirstLocalRow() ) ];
    for( Int k=0; k<yl.LocalHeight(); ++k )
        ++sendCounts[ y.RowOwner( k+yl.FirstLocalRow()+n ) ];
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    
    // Convert the send/recv counts into offsets and total sizes
    // =========================================================
    std::vector<int> sendOffsets, recvOffsets;
    const int totalSend = Scan( sendCounts, sendOffsets );
    const int totalRecv = Scan( recvCounts, recvOffsets );

    // Pack the triplets
    // =================
    std::vector<int> sSendBuf(totalSend);
    std::vector<Real> vSendBuf(totalSend);
    auto offsets = sendOffsets;
    for( Int k=0; k<yx.LocalHeight(); ++k )
    {
        const Int i = k + yx.FirstLocalRow();
        const Real value = yx.GetLocal(k,0);
        const Int owner = y.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }
    for( Int k=0; k<yl.LocalHeight(); ++k )
    {
        const Int i = k + yl.FirstLocalRow() + n;
        const Real value = yl.GetLocal(k,0);
        const Int owner = y.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }

    // Exchange and unpack the triplets
    // ================================
    std::vector<int> sRecvBuf(totalRecv);
    std::vector<Real> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    Zeros( y, m+n, 1 );
    for( Int k=0; k<totalRecv; ++k )
        y.UpdateLocal( sRecvBuf[k]-y.FirstLocalRow(), 0, vRecvBuf[k] );
}

template<typename Real>
void ExpandAugmentedSolution
( const Matrix<Real>& s, const Matrix<Real>& x,
  const Matrix<Real>& rmu, const Matrix<Real>& y,
  Matrix<Real>& ds, Matrix<Real>& dx, Matrix<Real>& dl )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::ExpandAugmentedSolution"))
    const Int n = rmu.Height();
    const Int m = y.Height() - n;

    // Extract dx and dl from [dx; dl]
    // ===============================
    dx.Resize( n, 1 );
    dl.Resize( m, 1 );
    const IR xInd(0,n), lInd(n,n+m);
    auto yx = y(xInd,IR(0,1));
    auto yl = y(lInd,IR(0,1));
    dx = yx;
    dl = yl;

    // ds := -(r_mu + S dx) / X
    // ========================
    ds.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
    {
        const Real x_i = x.Get(i,0);
        const Real s_i = s.Get(i,0);
        const Real dx_i = dx.Get(i,0);
        const Real rmu_i = rmu.Get(i,0);
        ds.Set( i, 0, -(rmu_i + s_i*dx_i)/x_i );
    }
}

template<typename Real>
void ExpandAugmentedSolution
( const AbstractDistMatrix<Real>& sPre, const AbstractDistMatrix<Real>& xPre,
  const AbstractDistMatrix<Real>& rmuPre, const AbstractDistMatrix<Real>& yPre,
  AbstractDistMatrix<Real>& dsPre, AbstractDistMatrix<Real>& dxPre, 
  AbstractDistMatrix<Real>& dl )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::ExpandAugmentedSolution"))

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.rowConstrain = true;
    ctrl.colAlign = 0;
    ctrl.rowAlign = 0;

    auto sPtr = ReadProxy<Real,MC,MR>(&sPre,ctrl); auto& s = *sPtr;
    auto xPtr = ReadProxy<Real,MC,MR>(&xPre,ctrl); auto& x = *xPtr;

    auto rmuPtr = ReadProxy<Real,MC,MR>(&rmuPre); auto& rmu = *rmuPtr;
    auto yPtr   = ReadProxy<Real,MC,MR>(&yPre);   auto& y   = *yPtr;

    auto dsPtr = WriteProxy<Real,MC,MR>(&dsPre,ctrl); auto& ds = *dsPtr;
    auto dxPtr = WriteProxy<Real,MC,MR>(&dxPre,ctrl); auto& dx = *dxPtr;

    const Int n = rmu.Height();
    const Int m = y.Height() - n;

    // Extract dx and dl from [dx; dl]
    // ===============================
    dx.Resize( n, 1 );
    dl.Resize( m, 1 );
    const IR xInd(0,n), lInd(n,n+m);
    auto yx = y(xInd,IR(0,1));
    auto yl = y(lInd,IR(0,1));
    dx = yx;
    Copy( yl, dl );

    // ds := -(r_mu + S dx) / X
    // ========================
    ds.Resize( n, 1 );
    for( Int iLoc=0; iLoc<ds.LocalHeight(); ++iLoc )
    {
        const Real x_i = x.GetLocal(iLoc,0);
        const Real s_i = s.GetLocal(iLoc,0);
        const Real dx_i = dx.GetLocal(iLoc,0);
        const Real rmu_i = rmu.GetLocal(iLoc,0);
        ds.SetLocal( iLoc, 0, -(rmu_i + s_i*dx_i)/x_i );
    }
}

template<typename Real>
void ExpandAugmentedSolution
( const DistMultiVec<Real>& s, const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& y,
  DistMultiVec<Real>& ds, DistMultiVec<Real>& dx, DistMultiVec<Real>& dl )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::ExpandAugmentedSolution"))
    const Int n = rmu.Height();
    const Int m = y.Height() - n;
    mpi::Comm comm = s.Comm();
    const Int commSize = mpi::Size(comm);

    // Extract dx and dl from [dx; dl]
    // ===============================
    dx.Resize( n, 1 );
    dl.Resize( m, 1 );
    // Compute the number of entries to send to each process
    // -----------------------------------------------------
    std::vector<int> sendCounts(commSize,0);
    for( Int iLoc=0; iLoc<y.LocalHeight(); ++iLoc )
    {
        const Int i = y.FirstLocalRow() + iLoc;
        if( i < n )
            ++sendCounts[ dx.RowOwner(i) ];
        else
            ++sendCounts[ dl.RowOwner(i-n) ];
    }
    // Communicate to determine the number we receive from each process
    // ----------------------------------------------------------------
    std::vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    std::vector<int> sendOffsets, recvOffsets;
    int totalSend = Scan( sendCounts, sendOffsets );
    int totalRecv = Scan( recvCounts, recvOffsets );
    // Pack the entries and row indices of dx and dl
    // ---------------------------------------------
    std::vector<Int> sSendBuf(totalSend);
    std::vector<Real> vSendBuf(totalSend);
    auto offsets = sendOffsets;
    for( Int iLoc=0; iLoc<y.LocalHeight(); ++iLoc )
    {
        const Int i = y.FirstLocalRow() + iLoc;
        if( i < n )
        {
            const Int owner = dx.RowOwner(i); 
            sSendBuf[offsets[owner]] = i; 
            vSendBuf[offsets[owner]] = y.GetLocal(iLoc,0);
            ++offsets[owner]; 
        }
        else
        {
            const Int owner = dl.RowOwner(i-n);
            sSendBuf[offsets[owner]] = i;
            vSendBuf[offsets[owner]] = y.GetLocal(iLoc,0);
            ++offsets[owner];
        }
    }
    // Exchange and unpack the entries and indices
    // -------------------------------------------
    std::vector<Int> sRecvBuf(totalRecv);
    std::vector<Real> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    for( Int k=0; k<totalRecv; ++k )
    {
        const Int i = sRecvBuf[k];
        if( i < n )
            dx.SetLocal( i-dx.FirstLocalRow(), 0, vRecvBuf[k] );
        else
            dl.SetLocal( i-n-dl.FirstLocalRow(), 0, vRecvBuf[k] );
    }

    // ds := -(r_mu + S dx) / X
    // ========================
    ds.Resize( n, 1 );
    for( Int iLoc=0; iLoc<ds.LocalHeight(); ++iLoc )
    {
        const Real x_i = x.GetLocal(iLoc,0);
        const Real s_i = s.GetLocal(iLoc,0);
        const Real dx_i = dx.GetLocal(iLoc,0);
        const Real rmu_i = rmu.GetLocal(iLoc,0);
        ds.SetLocal( iLoc, 0, -(rmu_i + s_i*dx_i)/x_i );
    }
}

#define PROTO(Real) \
  template void AugmentedKKT \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& s, const Matrix<Real>& x, \
    Matrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& x, \
    AbstractDistMatrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& s, const Matrix<Real>& x, \
    SparseMatrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, \
    DistSparseMatrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKTRHS \
  ( const Matrix<Real>& x, \
    const Matrix<Real>& rmu, const Matrix<Real>& rc, const Matrix<Real>& rb, \
    Matrix<Real>& y ); \
  template void AugmentedKKTRHS \
  ( const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc, \
    const AbstractDistMatrix<Real>& rb, AbstractDistMatrix<Real>& y ); \
  template void AugmentedKKTRHS \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rc, \
    const DistMultiVec<Real>& rb, DistMultiVec<Real>& y ); \
  template void ExpandAugmentedSolution \
  ( const Matrix<Real>& s, const Matrix<Real>& x, \
    const Matrix<Real>& rmu, const Matrix<Real>& y, \
    Matrix<Real>& ds, Matrix<Real>& dx, Matrix<Real>& dl ); \
  template void ExpandAugmentedSolution \
  ( const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& y, \
    AbstractDistMatrix<Real>& ds, AbstractDistMatrix<Real>& dx, \
    AbstractDistMatrix<Real>& dl ); \
  template void ExpandAugmentedSolution \
  ( const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& y, \
    DistMultiVec<Real>& ds, DistMultiVec<Real>& dx, DistMultiVec<Real>& dl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace lin_prog
} // namespace El
