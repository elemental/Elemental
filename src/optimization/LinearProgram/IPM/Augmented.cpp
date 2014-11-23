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
//    J = | -S*inv(X) A^T | and y = | -r_c + s - tau inv(X) e |,
//        |  A        0   |         | -r_b                    |
// where 
//    S   = diag(s),
//    X   = diag(x),
//    e   = ones(n,1),
//    r_b = A x - b, and
//    r_c = A^T l + s - c.
//
// The implied system is of the form
//   J | \Delta x | = y,
//     | \Delta l |
// and \Delta s = -s + tau inv(X) e - inv(X) S \Delta x.
//

template<typename Real>
void FormAugmentedSystem
( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& s, const Matrix<Real>& x, const Matrix<Real>& l,
  Real tau, Matrix<Real>& J, Matrix<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::FormAugmentedSystem"))
    const Int m = A.Height();
    const Int n = A.Width();

    // Form the Jacobian, J
    // ====================
    Zeros( J, m+n, m+n );
    IR xInd(0,n);
    IR lInd(n,n+m);
    auto Jxx = J(xInd,xInd); auto Jxl = J(xInd,lInd);
    auto Jlx = J(lInd,xInd); auto Jll = J(lInd,lInd);
    Matrix<Real> d( s );
    Scale( Real(-1), d );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( Jxx, d );
    Transpose( A, Jxl );
    Jlx = A;

    // Form the right-hand side, y
    // ===========================
    Zeros( y, m+n, 1 );
    auto yx = y(xInd,IR(0,1));
    auto yl = y(lInd,IR(0,1));
    yx = c;
    Gemv( TRANSPOSE, Real(-1), A, l, Real(1), yx );
    for( Int j=0; j<n; ++j )
        yx.Update( j, 0, -tau/x.Get(j,0) );
    yl = b;
    Gemv( NORMAL, Real(-1), A, x, Real(1), yl );
}

template<typename Real>
void FormAugmentedSystem
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Real>& l,
  Real tau, AbstractDistMatrix<Real>& JPre, AbstractDistMatrix<Real>& yPre )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::FormAugmentedSystem"))
    const Int m = A.Height();
    const Int n = A.Width();

    auto JPtr = WriteProxy<Real,MC,MR>(&JPre); auto& J = *JPtr;
    auto yPtr = WriteProxy<Real,MC,MR>(&yPre); auto& y = *yPtr;

    // Form the Jacobian, J
    // ====================
    Zeros( J, m+n, m+n );
    IR xInd(0,n), lInd(n,n+m);
    auto Jxx = J(xInd,xInd); auto Jxl = J(xInd,lInd);
    auto Jlx = J(lInd,xInd); auto Jll = J(lInd,lInd);
    DistMatrix<Real,STAR,STAR> d( s );
    Scale( Real(-1), d );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( Jxx, d.Matrix() );
    Transpose( A, Jxl );
    Jlx = A;

    // Form the right-hand side, y
    // ===========================
    Zeros( y, m+n, 1 );

    auto yx = y(xInd,IR(0,1));
    yx = c;
    Gemv( TRANSPOSE, Real(-1), A, l, Real(1), yx );
    DistMatrix<Real> g(x);
    auto map = [&]( Real alpha ) { return -tau/alpha; };
    EntrywiseMap( g, std::function<Real(Real)>(map) );
    Axpy( Real(1), g, yx );

    auto yl = y(lInd,IR(0,1));
    yl = b;
    Gemv( NORMAL, Real(-1), A, x, Real(1), yl );
}

template<typename Real>
void FormAugmentedSystem
( const SparseMatrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& s, const Matrix<Real>& x, const Matrix<Real>& l,
  Real tau, SparseMatrix<Real>& J, Matrix<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::FormAugmentedSystem"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();

    // Form the Jacobian, J
    // ====================
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
        J.Update( j, i+n, value );
    }
    J.MakeConsistent();

    // Form the right-hand side, y
    // ===========================
    Zeros( y, m+n, 1 );
    IR xInd(0,n), lInd(n,n+m);

    auto yx = y(xInd,IR(0,1));
    yx = c;
    Multiply( TRANSPOSE, Real(-1), A, l, Real(1), yx );
    for( Int j=0; j<n; ++j )
        yx.Update( j, 0, -tau/x.Get(j,0) );

    auto yl = y(lInd,IR(0,1));
    yl = b;
    Multiply( NORMAL, Real(-1), A, x, Real(1), yl );
}

template<typename Real>
void FormAugmentedSystem
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& l,
  Real tau, DistSparseMatrix<Real>& J, DistMultiVec<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::FormAugmentedSystem"))
    const Int m = A.Height();
    const Int n = A.Width();

    mpi::Comm comm = A.Comm();
    const Int commSize = mpi::Size( comm );

    DistSparseMatrix<Real> ATrans(comm);
    Transpose( A, ATrans );

    // Form the Jacobian
    // =================
    J.SetComm( comm );
    Zeros( J, m+n, m+n );
    
    // Compute the number of entries to send to each process
    // -----------------------------------------------------
    std::vector<int> sendCounts(commSize,0);
    // For placing A into the bottom-left corner
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for( Int k=0; k<A.NumLocalEntries(); ++k )
        ++sendCounts[ J.RowOwner(A.Row(k)+n) ];
    // For placing A^T into the top-right corner
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for( Int k=0; k<ATrans.NumLocalEntries(); ++k )
        ++sendCounts[ J.RowOwner(ATrans.Row(k)) ];
    // For placing -S*inv(X) into the top-left corner
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for( Int k=0; k<x.LocalHeight(); ++k )
        ++sendCounts[ J.RowOwner( k+x.FirstLocalRow() ) ];
    // Communicate to determine the number we receive from each process
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    
    // Convert the send/recv counts into offsets and total sizes
    // ---------------------------------------------------------
    std::vector<int> sendOffsets, recvOffsets;
    int totalSend = Scan( sendCounts, sendOffsets );
    int totalRecv = Scan( recvCounts, recvOffsets );

    // Pack the triplets
    // -----------------
    std::vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
    std::vector<Real> vSendBuf(totalSend);
    auto offsets = sendOffsets;
    // Pack A
    // ^^^^^^
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
    // ^^^^^^^^
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
    // Pack -S*inv(X)
    // ^^^^^^^^^^^^^^
    for( Int k=0; k<x.LocalHeight(); ++k )
    {
        const Int i = k + x.FirstLocalRow();
        const Int j = i;
        const Int value = -s.GetLocal(k,0)/x.GetLocal(k,0);
        const Int owner = J.RowOwner(i);
        sSendBuf[offsets[owner]] = i; 
        tSendBuf[offsets[owner]] = j;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }

    // Exchange and unpack the triplets
    // --------------------------------
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

    // Form the two halves of the right-hand side
    // ==========================================
    Zeros( y, m+n, 1 );
    DistMultiVec<Real> yx(comm), yl(comm);

    yx = c;
    Multiply( NORMAL, Real(-1), ATrans, l, Real(1), yx );
    for( Int k=0; k<x.LocalHeight(); ++k )
        yx.UpdateLocal( k, 0, -tau/x.GetLocal(k,0) );

    yl = b;
    Multiply( NORMAL, Real(-1), A, x, Real(1), yl );

    // Compute the number of entries to send to each process
    // -----------------------------------------------------
    for( Int q=0; q<commSize; ++q )
        sendCounts[commSize] = 0;
    for( Int k=0; k<yx.LocalHeight(); ++k )
        ++sendCounts[ y.RowOwner( k+yx.FirstLocalRow() ) ];
    for( Int k=0; k<yl.LocalHeight(); ++k )
        ++sendCounts[ y.RowOwner( k+yl.FirstLocalRow()+n ) ];
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    
    // Convert the send/recv counts into offsets and total sizes
    // ---------------------------------------------------------
    totalSend = Scan( sendCounts, sendOffsets );
    totalRecv = Scan( recvCounts, recvOffsets );

    // Pack the triplets
    // -----------------
    sSendBuf.resize(totalSend);
    vSendBuf.resize(totalSend);
    offsets = sendOffsets;
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
    // --------------------------------
    sRecvBuf.resize(totalRecv);
    vRecvBuf.resize(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    Zeros( y, m+n, 1 );
    for( Int k=0; k<totalRecv; ++k )
        y.UpdateLocal( sRecvBuf[k]-J.FirstLocalRow(), 0, vRecvBuf[k] );
}

#define PROTO(Real) \
  template void FormAugmentedSystem \
  ( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& s, const Matrix<Real>& x, const Matrix<Real>& l, \
    Real tau, Matrix<Real>& J, Matrix<Real>& y ); \
  template void FormAugmentedSystem \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Real>& l, \
    Real tau, AbstractDistMatrix<Real>& J, AbstractDistMatrix<Real>& y ); \
  template void FormAugmentedSystem \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& s, const Matrix<Real>& x, const Matrix<Real>& l, \
    Real tau, SparseMatrix<Real>& J, Matrix<Real>& y ); \
  template void FormAugmentedSystem \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& l, \
    Real tau, DistSparseMatrix<Real>& J, DistMultiVec<Real>& y );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace lin_prog
} // namespace El
