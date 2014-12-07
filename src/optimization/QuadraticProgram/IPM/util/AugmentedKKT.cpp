/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace quad_prog {

// Form 
//    J = | -Q-Z*inv(X) A^T | and y = | -r_c + r_mu / X |,
//        |  A          0   |         | -r_b            |
// where 
//    Z   = diag(z),
//    X   = diag(x),
//    e   = ones(n,1),
//    r_b = A x - b,
//    r_c = A^T y + z - Q x - c, and
//    r_mu = X Z e - tau e.
//
// The implied system is of the form
//   J | dx | = rhs,   dz = -(r_mu + Z dx) / X.
//     | dy |
//

template<typename Real>
void AugmentedKKT
( const Matrix<Real>& Q, const Matrix<Real>& A, 
  const Matrix<Real>& x, const Matrix<Real>& z,
        Matrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::AugmentedKKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    Zeros( J, m+n, m+n );
    const IR xInd(0,n), yInd(n,n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd);
    Matrix<Real> d( z );
    Scale( Real(-1), d );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( Jxx, d );
    Axpy( Real(-1), Q, Jxx );
    Jyx = A;
    if( !onlyLower )
        Transpose( A, Jxy );
}

template<typename Real>
void AugmentedKKT
( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z,
  AbstractDistMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::AugmentedKKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    auto JPtr = WriteProxy<Real,MC,MR>(&JPre);
    auto& J = *JPtr;

    Zeros( J, m+n, m+n );
    const IR xInd(0,n), yInd(n,n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd);
    DistMatrix<Real,STAR,STAR> d( z );
    Scale( Real(-1), d );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( Jxx, d.Matrix() );
    Axpy( Real(-1), Q, Jxx );
    Jyx = A;
    if( !onlyLower )
        Transpose( A, Jxy );
}

template<typename Real>
void AugmentedKKT
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, 
  const Matrix<Real>& x, const Matrix<Real>& z,
  SparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::AugmentedKKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    Zeros( J, m+n, m+n );
    if( onlyLower )
        J.Reserve( A.NumEntries() + Q.NumEntries() + n );
    else
        J.Reserve( 2*A.NumEntries() + Q.NumEntries() + n ); 

    // -Z*inv(X) updates
    for( Int j=0; j<n; ++j )
        J.Update( j, j, -z.Get(j,0)/x.Get(j,0) );
    // -Q update
    for( Int k=0; k<Q.NumEntries(); ++k )
        J.Update( Q.Row(k), Q.Col(k), -Q.Value(k) );
    // A and A^T updates
    for( Int k=0; k<A.NumEntries(); ++k )
    {
        J.Update( A.Row(k)+n, A.Col(k), A.Value(k) );
        if( !onlyLower )
            J.Update( A.Col(k), A.Row(k)+n, A.Value(k) );
    }
    J.MakeConsistent();
}

template<typename Real>
void AugmentedKKT
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& x, const DistMultiVec<Real>& z,
  DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::AugmentedKKT"))
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
    // For placing -Q into the top-left corner
    // ---------------------------------------
    for( Int k=0; k<Q.NumLocalEntries(); ++k )
        ++sendCounts[ J.RowOwner(Q.Row(k)) ];
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
    // Pack -Z inv(X)
    // --------------
    for( Int k=0; k<x.LocalHeight(); ++k )
    {
        const Int i = k + x.FirstLocalRow();
        const Int j = i;
        const Real value = -z.GetLocal(k,0)/x.GetLocal(k,0);
        const Int owner = J.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        tSendBuf[offsets[owner]] = j;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }
    // Pack -Q
    // -------
    for( Int k=0; k<Q.NumLocalEntries(); ++k )
    {
        const Int i = Q.Row(k);
        const Int j = Q.Col(k);
        const Real value = -Q.Value(k);
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
  Matrix<Real>& rhs )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::AugmentedKKTRHS"))
    const Int m = rb.Height();
    const Int n = rmu.Height();
    const IR xInd(0,n), yInd(n,n+m);
    Zeros( rhs, m+n, 1 );

    auto rhs_x = rhs(xInd,IR(0,1));
    rhs_x = rmu;
    for( Int i=0; i<n; ++i )
        rhs_x.Set( i, 0, rhs_x.Get(i,0)/x.Get(i,0) );
    Axpy( Real(-1), rc, rhs_x );

    auto rhs_y = rhs(yInd,IR(0,1));
    rhs_y = rb;
    Scale( Real(-1), rhs_y );
}

template<typename Real>
void AugmentedKKTRHS
( const AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc, 
  const AbstractDistMatrix<Real>& rb,
  AbstractDistMatrix<Real>& rhsPre )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::AugmentedKKTRHS"))

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.rowConstrain = true;
    ctrl.colAlign = 0;
    ctrl.rowAlign = 0;
    
    auto xPtr   = ReadProxy<Real,MC,MR>(&xPre,ctrl);    auto& x   = *xPtr;
    auto rhsPtr = WriteProxy<Real,MC,MR>(&rhsPre,ctrl); auto& rhs = *rhsPtr;

    const Int m = rb.Height();
    const Int n = rmu.Height();
    const IR xInd(0,n), yInd(n,n+m);
    Zeros( rhs, m+n, 1 );

    auto rhs_x = rhs(xInd,IR(0,1));
    rhs_x = rmu;
    for( Int iLoc=0; iLoc<rhs_x.LocalHeight(); ++iLoc )
        rhs_x.SetLocal( iLoc, 0, rhs_x.GetLocal(iLoc,0)/x.GetLocal(iLoc,0) );
    Axpy( Real(-1), rc, rhs_x );

    auto rhs_y = rhs(yInd,IR(0,1));
    rhs_y = rb;
    Scale( Real(-1), rhs_y );
}

template<typename Real>
void AugmentedKKTRHS
( const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rc, 
  const DistMultiVec<Real>& rb, DistMultiVec<Real>& rhs )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::FormAugmentedSystem"))
    const Int m = rb.Height();
    const Int n = x.Height();

    mpi::Comm comm = x.Comm();
    const Int commSize = mpi::Size( comm );

    Zeros( rhs, m+n, 1 );
    DistMultiVec<Real> rhs_x(comm), rhs_y(comm);

    rhs_x = rc;
    Scale( Real(-1), rhs_x );
    for( Int k=0; k<x.LocalHeight(); ++k )
        rhs_x.UpdateLocal( k, 0, rmu.GetLocal(k,0)/x.GetLocal(k,0) );    

    rhs_y = rb;
    Scale( Real(-1), rhs_y );

    // Compute the number of entries to send to each process
    // =====================================================
    std::vector<int> sendCounts(commSize), recvCounts(commSize);
    for( Int q=0; q<commSize; ++q )
        sendCounts[q] = 0;
    for( Int k=0; k<rhs_x.LocalHeight(); ++k )
        ++sendCounts[ rhs.RowOwner( k+rhs_x.FirstLocalRow() ) ];
    for( Int k=0; k<rhs_y.LocalHeight(); ++k )
        ++sendCounts[ rhs.RowOwner( k+rhs_y.FirstLocalRow()+n ) ];
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
    for( Int k=0; k<rhs_x.LocalHeight(); ++k )
    {
        const Int i = k + rhs_x.FirstLocalRow();
        const Real value = rhs_x.GetLocal(k,0);
        const Int owner = rhs.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }
    for( Int k=0; k<rhs_y.LocalHeight(); ++k )
    {
        const Int i = k + rhs_y.FirstLocalRow() + n;
        const Real value = rhs_y.GetLocal(k,0);
        const Int owner = rhs.RowOwner(i);
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
    Zeros( rhs, m+n, 1 );
    for( Int k=0; k<totalRecv; ++k )
        rhs.UpdateLocal( sRecvBuf[k]-rhs.FirstLocalRow(), 0, vRecvBuf[k] );
}

template<typename Real>
void ExpandAugmentedSolution
( const Matrix<Real>& x, const Matrix<Real>& z,
  const Matrix<Real>& rmu, const Matrix<Real>& rhs,
  Matrix<Real>& dx, Matrix<Real>& dy, Matrix<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::ExpandAugmentedSolution"))
    const Int n = rmu.Height();
    const Int m = rhs.Height() - n;

    // Extract dx and dy from [dx; dy]
    // ===============================
    dx.Resize( n, 1 );
    dy.Resize( m, 1 );
    const IR xInd(0,n), yInd(n,n+m);
    auto rhs_x = rhs(xInd,IR(0,1));
    auto rhs_y = rhs(yInd,IR(0,1));
    dx = rhs_x;
    dy = rhs_y;

    // dz := -(r_mu + Z dx) / X
    // ========================
    dz.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
    {
        const Real x_i = x.Get(i,0);
        const Real z_i = z.Get(i,0);
        const Real dx_i = dx.Get(i,0);
        const Real rmu_i = rmu.Get(i,0);
        dz.Set( i, 0, -(rmu_i + z_i*dx_i)/x_i );
    }
}

template<typename Real>
void ExpandAugmentedSolution
( const AbstractDistMatrix<Real>& xPre, const AbstractDistMatrix<Real>& zPre,
  const AbstractDistMatrix<Real>& rmuPre, 
  const AbstractDistMatrix<Real>& rhsPre,
        AbstractDistMatrix<Real>& dxPre, 
        AbstractDistMatrix<Real>& dy, 
        AbstractDistMatrix<Real>& dzPre )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::ExpandAugmentedSolution"))

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.rowConstrain = true;
    ctrl.colAlign = 0;
    ctrl.rowAlign = 0;

    auto xPtr = ReadProxy<Real,MC,MR>(&xPre,ctrl); auto& x = *xPtr;
    auto zPtr = ReadProxy<Real,MC,MR>(&zPre,ctrl); auto& z = *zPtr;

    auto rmuPtr = ReadProxy<Real,MC,MR>(&rmuPre); auto& rmu = *rmuPtr;
    auto rhsPtr = ReadProxy<Real,MC,MR>(&rhsPre); auto& rhs = *rhsPtr;

    auto dxPtr = WriteProxy<Real,MC,MR>(&dxPre,ctrl); auto& dx = *dxPtr;
    auto dzPtr = WriteProxy<Real,MC,MR>(&dzPre,ctrl); auto& dz = *dzPtr;

    const Int n = rmu.Height();
    const Int m = rhs.Height() - n;

    // Extract dx and dy from [dx; dy]
    // ===============================
    dx.Resize( n, 1 );
    dy.Resize( m, 1 );
    const IR xInd(0,n), yInd(n,n+m);
    auto rhs_x = rhs(xInd,IR(0,1));
    auto rhs_y = rhs(yInd,IR(0,1));
    dx = rhs_x;
    Copy( rhs_y, dy );

    // dz := -(r_mu + Z dx) / X
    // ========================
    dz.Resize( n, 1 );
    for( Int iLoc=0; iLoc<dz.LocalHeight(); ++iLoc )
    {
        const Real x_i = x.GetLocal(iLoc,0);
        const Real z_i = z.GetLocal(iLoc,0);
        const Real dx_i = dx.GetLocal(iLoc,0);
        const Real rmu_i = rmu.GetLocal(iLoc,0);
        dz.SetLocal( iLoc, 0, -(rmu_i + z_i*dx_i)/x_i );
    }
}

template<typename Real>
void ExpandAugmentedSolution
( const DistMultiVec<Real>& x, const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rhs,
  DistMultiVec<Real>& dx, DistMultiVec<Real>& dy, DistMultiVec<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::ExpandAugmentedSolution"))
    const Int n = rmu.Height();
    const Int m = rhs.Height() - n;
    mpi::Comm comm = z.Comm();
    const Int commSize = mpi::Size(comm);

    // Extract dx and dy from [dx; dy]
    // ===============================
    dx.Resize( n, 1 );
    dy.Resize( m, 1 );
    // Compute the number of entries to send to each process
    // -----------------------------------------------------
    std::vector<int> sendCounts(commSize,0);
    for( Int iLoc=0; iLoc<rhs.LocalHeight(); ++iLoc )
    {
        const Int i = rhs.FirstLocalRow() + iLoc;
        if( i < n )
            ++sendCounts[ dx.RowOwner(i) ];
        else
            ++sendCounts[ dy.RowOwner(i-n) ];
    }
    // Communicate to determine the number we receive from each process
    // ----------------------------------------------------------------
    std::vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    std::vector<int> sendOffsets, recvOffsets;
    int totalSend = Scan( sendCounts, sendOffsets );
    int totalRecv = Scan( recvCounts, recvOffsets );
    // Pack the entries and row indices of dx and dy
    // ---------------------------------------------
    std::vector<Int> sSendBuf(totalSend);
    std::vector<Real> vSendBuf(totalSend);
    auto offsets = sendOffsets;
    for( Int iLoc=0; iLoc<rhs.LocalHeight(); ++iLoc )
    {
        const Int i = rhs.FirstLocalRow() + iLoc;
        if( i < n )
        {
            const Int owner = dx.RowOwner(i); 
            sSendBuf[offsets[owner]] = i; 
            vSendBuf[offsets[owner]] = rhs.GetLocal(iLoc,0);
            ++offsets[owner]; 
        }
        else
        {
            const Int owner = dy.RowOwner(i-n);
            sSendBuf[offsets[owner]] = i;
            vSendBuf[offsets[owner]] = rhs.GetLocal(iLoc,0);
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
            dy.SetLocal( i-n-dy.FirstLocalRow(), 0, vRecvBuf[k] );
    }

    // dz := -(r_mu + Z dx) / X
    // ========================
    dz.Resize( n, 1 );
    for( Int iLoc=0; iLoc<dz.LocalHeight(); ++iLoc )
    {
        const Real x_i = x.GetLocal(iLoc,0);
        const Real z_i = z.GetLocal(iLoc,0);
        const Real dx_i = dx.GetLocal(iLoc,0);
        const Real rmu_i = rmu.GetLocal(iLoc,0);
        dz.SetLocal( iLoc, 0, -(rmu_i + z_i*dx_i)/x_i );
    }
}

#define PROTO(Real) \
  template void AugmentedKKT \
  ( const Matrix<Real>& Q, const Matrix<Real>& A, \
    const Matrix<Real>& x, const Matrix<Real>& z, \
    Matrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z, \
    AbstractDistMatrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, \
    const Matrix<Real>& x, const Matrix<Real>& z, \
    SparseMatrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& x, const DistMultiVec<Real>& z, \
    DistSparseMatrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKTRHS \
  ( const Matrix<Real>& x, \
    const Matrix<Real>& rmu, const Matrix<Real>& rc, const Matrix<Real>& rb, \
    Matrix<Real>& rhs ); \
  template void AugmentedKKTRHS \
  ( const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc, \
    const AbstractDistMatrix<Real>& rb, AbstractDistMatrix<Real>& rhs ); \
  template void AugmentedKKTRHS \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rc, \
    const DistMultiVec<Real>& rb, DistMultiVec<Real>& rhs ); \
  template void ExpandAugmentedSolution \
  ( const Matrix<Real>& x, const Matrix<Real>& z, \
    const Matrix<Real>& rmu, const Matrix<Real>& rhs, \
    Matrix<Real>& dx, Matrix<Real>& dy, Matrix<Real>& dz ); \
  template void ExpandAugmentedSolution \
  ( const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rhs, \
    AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, \
    AbstractDistMatrix<Real>& dz ); \
  template void ExpandAugmentedSolution \
  ( const DistMultiVec<Real>& x, const DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rhs, \
    DistMultiVec<Real>& dx, DistMultiVec<Real>& dy, DistMultiVec<Real>& dz );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace quad_prog
} // namespace El
