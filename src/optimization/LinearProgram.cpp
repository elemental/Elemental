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
//
//    J = | 0    A      | and y = |         -r_b            |,
//        | A^T -D^{-2} |         | -r_c + s - tau X^{-1} e |
//
// where 
//
//    D   = diag(x)^{1/2} / diag(s)^{1/2},
//    e   = ones(n,1),
//    r_b = A x - b, and
//    r_c = A^T l + s - c.
//
// The implied system is of the form
//
//   J | \Delta l | = y,
//     | \Delta x |
//
// and \Delta s = -s + tau diag(x)^{-1} e - diag(x)^{-1} diag(s) \Delta x.
//

template<typename Real>
void FormAugmentedSystem
( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& x, const Matrix<Real>& l, const Matrix<Real>& s,
  Real tau, Matrix<Real>& J, Matrix<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::FormAugmentedSystem"))
    const Int m = A.Height();
    const Int n = A.Width();

    // Form the Jacobian, J
    // ====================
    Zeros( J, m+n, m+n );
    auto JTR = J(IR(0,m  ),IR(m,m+n));
    auto JBL = J(IR(m,m+n),IR(0,m  ));
    auto JBR = J(IR(m,m+n),IR(m,m+n));
    JTR = A;
    Transpose( A, JBL );
    Matrix<Real> d( s );
    Scale( Real(-1), d );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( JBR, d );

    // Form the right-hand side, y
    // ===========================
    Zeros( y, m+n, 1 );
    auto yT = y(IR(0,m),IR(0,1));
    yT = b;
    Gemv( NORMAL, Real(-1), A, x, Real(1), yT );
    auto yB = y(IR(m,m+n),IR(0,1));
    yB = c;
    Gemv( TRANSPOSE, Real(-1), A, l, Real(1), yB );
    for( Int j=0; j<n; ++j )
        yB.Update( j, 0, -tau/x.Get(j,0) );
}

template<typename Real>
void FormAugmentedSystem
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& l, 
  const AbstractDistMatrix<Real>& s,
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
    auto JTR = J(IR(0,m  ),IR(m,m+n));
    auto JBL = J(IR(m,m+n),IR(0,m  ));
    auto JBR = J(IR(m,m+n),IR(m,m+n));
    JTR = A;
    Transpose( A, JBL );
    DistMatrix<Real,STAR,STAR> d( s );
    Scale( Real(-1), d );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( JBR, d.Matrix() );

    // Form the right-hand side, y
    // ===========================
    Zeros( y, m+n, 1 );
    auto yT = y(IR(0,m),IR(0,1));
    yT = b;
    Gemv( NORMAL, Real(-1), A, x, Real(1), yT );
    auto yB = y(IR(m,m+n),IR(0,1));
    yB = c;
    Gemv( TRANSPOSE, Real(-1), A, l, Real(1), yB );
    DistMatrix<Real> g( x );
    auto lambda = [&]( Real alpha ) { return -tau/alpha; };
    EntrywiseMap( g, std::function<Real(Real)>(lambda) );
    Axpy( Real(1), g, yB );
}

template<typename Real>
void FormAugmentedSystem
( const SparseMatrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& x, const Matrix<Real>& l, const Matrix<Real>& s,
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
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k);
        const Real value = A.Value(k);
        // A update
        J.Update( i, j+m, value );
        // A^T update
        J.Update( j+m, i, value );
    }
    // -D^{-2} updates
    for( Int j=0; j<n; ++j )
        J.Update( m+j, m+j, -s.Get(j,0)/x.Get(j,0) );
    J.MakeConsistent();

    // Form the right-hand side, y
    // ===========================
    Zeros( y, m+n, 1 );
    auto yT = y(IR(0,m),IR(0,1));
    yT = b;
    Multiply( NORMAL, Real(-1), A, x, Real(1), yT );
    auto yB = y(IR(m,m+n),IR(0,1));
    yB = c;
    Multiply( TRANSPOSE, Real(-1), A, l, Real(1), yB );
    for( Int j=0; j<n; ++j )
        yB.Update( j, 0, -tau/x.Get(j,0) );
}

template<typename Real>
void FormAugmentedSystem
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& x, const DistMultiVec<Real>& l, 
  const DistMultiVec<Real>& s,
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
    // For placing A into the top-right corner
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for( Int k=0; k<A.NumLocalEntries(); ++k )
        ++sendCounts[ J.RowOwner(A.Row(k)) ];
    // For placing A^T into the bottom-left corner
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for( Int k=0; k<ATrans.NumLocalEntries(); ++k )
        ++sendCounts[ J.RowOwner( ATrans.Row(k)+m ) ];
    // For placing -diag(s)/diag(x) into the bottom-right corner
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for( Int k=0; k<x.LocalHeight(); ++k )
        ++sendCounts[ J.RowOwner( k+x.FirstLocalRow()+m ) ];
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
        const Int i = A.Row(k);
        const Int j = A.Col(k) + m;
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
        const Int i = ATrans.Row(k) + m;
        const Int j = ATrans.Col(k);
        const Real value = ATrans.Value(k);
        const Int owner = J.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        tSendBuf[offsets[owner]] = j;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }
    // Pack -D^{-2}
    // ^^^^^^^^^^^^
    for( Int k=0; k<x.LocalHeight(); ++k )
    {
        const Int i = k + x.FirstLocalRow() + m;
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
    DistMultiVec<Real> yT(comm), yB(comm);
    yT = b;
    yB = c;
    Multiply( NORMAL, Real(-1), A, x, Real(1), yT );
    Multiply( NORMAL, Real(-1), ATrans, l, Real(1), yB );
    for( Int k=0; k<x.LocalHeight(); ++k )
        yB.UpdateLocal( k, 0, -tau/x.GetLocal(k,0) );

    // Compute the number of entries to send to each process
    // -----------------------------------------------------
    for( Int q=0; q<commSize; ++q )
        sendCounts[commSize] = 0;
    for( Int k=0; k<yT.LocalHeight(); ++k )
        ++sendCounts[ y.RowOwner( k+yT.FirstLocalRow() ) ];
    for( Int k=0; k<yB.LocalHeight(); ++k )
        ++sendCounts[ y.RowOwner( k+yB.FirstLocalRow()+m ) ];
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
    for( Int k=0; k<yT.LocalHeight(); ++k )
    {
        const Int i = k + yT.FirstLocalRow();
        const Real value = yT.GetLocal(k,0);
        const Int owner = y.RowOwner(i);
        sSendBuf[offsets[owner]] = i;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }
    for( Int k=0; k<yB.LocalHeight(); ++k )
    {
        const Int i = k + yB.FirstLocalRow() + m;
        const Real value = yB.GetLocal(k,0);
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

// Form 
//
//    J = | A D^2 A^T |, and 
//    y = [ b - A diag(s)^{-1} ( diag(x) r_c + tau e ) ]
//
// where 
//
//    D   = diag(x)^{1/2} / diag(s)^{1/2},
//    e   = ones(n,1),
//    r_b = A x - b, and
//    r_c = A^T l + s - c.
//
// The implied system is of the form
//
//   J | \Delta l | = y,
// 
//  \Delta s = -r_c - A^T \Delta l, and
//  \Delta x = -x + diag(s)^{-1} (tau e - diag(x) \Delta s).
//

template<typename Real>
void FormNormalSystem
( const SparseMatrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& x, const Matrix<Real>& l, const Matrix<Real>& s,
  Real tau, SparseMatrix<Real>& J, Matrix<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::FormNormalSystem"))
    const Int n = A.Width();

    // Form D = diag(d) = diag(x)^{1/2} / diag(s)^{1/2}
    // ================================================
    Matrix<Real> d;
    d.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        d.Set( i, 0, Sqrt(x.Get(i,0))/Sqrt(s.Get(i,0)) );

    // Form the Jacobian, J
    // ====================
    SparseMatrix<Real> G;
    Transpose( A, G );
    DiagonalScale( LEFT, NORMAL, d, G );
    Syrk( LOWER, TRANSPOSE, Real(1), G, J );
    MakeSymmetric( LOWER, J );

    // Form the right-hand side, y
    // ===========================
    // Form the 'c' residual, A^T l + s - c
    // ------------------------------------
    Matrix<Real> cResid;
    cResid = s;
    Axpy( Real(-1), c, cResid );
    Multiply( TRANSPOSE, Real(1), A, l, Real(1), cResid );
    // Form the portion of the right-hand side to be multiplied by A
    // -------------------------------------------------------------
    Matrix<Real> g;
    g = cResid;
    for( Int i=0; i<n; ++i )
    {
        const Real gamma = g.Get(i,0); 
        const Real si = s.Get(i,0);
        const Real xi = x.Get(i,0);
        g.Set( i, 0, (gamma*xi + tau)/si );
    }
    // Form the right-hand side, y
    // ---------------------------
    y = b;
    Multiply( NORMAL, Real(-1), A, g, Real(1), y );
}

template<typename Real>
void SolveNormalSystem
( const SparseMatrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& x, const Matrix<Real>& l, const Matrix<Real>& s,
  Real tau, const SparseMatrix<Real>& J, const Matrix<Real>& y,
  Matrix<Real>& dx, Matrix<Real>& dl, Matrix<Real>& ds )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::SolveNormalSystem"))
    const Int n = A.Width();

    // NOTE: SymmetricSolve not yet supported for sequential matrices
    /*
    // Compute the proposed change in the Lagrange multiplier
    // ======================================================
    dl = y;
    SymmetricSolve( J, dl );

    // Compute the proposed change in the dual variable
    // ================================================
    // ds := c - s
    // -----------
    ds = c; 
    Axpy( Real(-1), s, ds );
    // g := l + dl
    // -----------
    Matrix<Real> g;
    g = l;
    Axpy( Real(1), dl, g );
    // ds := ds - A^T g = c - s - A^T (l + dl)
    // ---------------------------------------
    Multiply( TRANSPOSE, Real(-1), A, g, Real(1), ds );

    // Compute the proposed change in the primal variable
    // ==================================================
    Zeros( dx, n, 1 );
    for( Int i=0; i<n; ++i )
    {
        const Real xi = x.Get(i,0);
        const Real si = s.Get(i,0);
        const Real dsi = ds.Get(i,0);
        dx.Set( i, 0, -xi + (tau - dsi*xi)/si );
    }
    */
}

template<typename Real>
void FormNormalSystem
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& x, const DistMultiVec<Real>& l, 
  const DistMultiVec<Real>& s,
  Real tau, DistSparseMatrix<Real>& J, DistMultiVec<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::FormNormalSystem"))
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    if( !mpi::Congruent( comm, b.Comm() ) )
        LogicError("Communicators of A and b must match");
    if( !mpi::Congruent( comm, c.Comm() ) )
        LogicError("Communicators of A and c must match");
    if( !mpi::Congruent( comm, x.Comm() ) )
        LogicError("Communicators of A and x must match");
    if( !mpi::Congruent( comm, l.Comm() ) )
        LogicError("Communicators of A and l must match");
    if( !mpi::Congruent( comm, s.Comm() ) )
        LogicError("Communicators of A and s must match");

    // Form D = diag(d) = diag(x)^{1/2} / diag(s)^{1/2}
    // ================================================
    DistMultiVec<Real> d(comm);
    d.Resize( n, 1 );
    for( Int iLoc=0; iLoc<d.LocalHeight(); ++iLoc )
        d.SetLocal
        ( iLoc, 0, Sqrt(x.GetLocal(iLoc,0))/Sqrt(s.GetLocal(iLoc,0)) );

    // Form the Jacobian, J
    // ====================
    DistSparseMatrix<Real> G(comm);
    Transpose( A, G );
    DiagonalScale( LEFT, NORMAL, d, G );
    Syrk( LOWER, TRANSPOSE, Real(1), G, J );
    MakeSymmetric( LOWER, J );

    // Form the right-hand side, y
    // ===========================
    // Form the 'c' residual, A^T l + s - c
    // ------------------------------------
    DistMultiVec<Real> cResid(comm);
    cResid = s;
    Axpy( Real(-1), c, cResid );
    Multiply( TRANSPOSE, Real(1), A, l, Real(1), cResid );
    // Form the portion of the right-hand side to be multiplied by A
    // -------------------------------------------------------------
    DistMultiVec<Real> g(comm);
    g = cResid;
    for( Int iLoc=0; iLoc<g.LocalHeight(); ++iLoc )
    {
        const Real gamma = g.GetLocal(iLoc,0); 
        const Real si = s.GetLocal(iLoc,0);
        const Real xi = x.GetLocal(iLoc,0);
        g.SetLocal( iLoc, 0, (gamma*xi + tau)/si );
    }
    // Form the right-hand side, y
    // ---------------------------
    y = b;
    Multiply( NORMAL, Real(-1), A, g, Real(1), y );
}

template<typename Real>
void SolveNormalSystem
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& x, const DistMultiVec<Real>& l, 
  const DistMultiVec<Real>& s,
  Real tau, const DistSparseMatrix<Real>& J, const DistMultiVec<Real>& y,
  DistMultiVec<Real>& dx, DistMultiVec<Real>& dl, DistMultiVec<Real>& ds )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::SolveNormalSystem"))
    // TODO: Check that the communicators are congruent

    // Compute the proposed change in the Lagrange multiplier
    // ======================================================
    dl = y;
    SymmetricSolve( J, dl );

    // Compute the proposed change in the dual variable
    // ================================================
    // ds := c - s
    // -----------
    ds = c; 
    Axpy( Real(-1), s, ds );
    // g := l + dl
    // -----------
    DistMultiVec<Real> g(A.Comm());
    g = l;
    Axpy( Real(1), dl, g );
    // ds := ds - A^T g = c - s - A^T (l + dl)
    // ---------------------------------------
    Multiply( TRANSPOSE, Real(-1), A, g, Real(1), ds );

    // Compute the proposed change in the primal variable
    // ==================================================
    const Int n = A.Width();
    Zeros( dx, n, 1 );
    const Int nLoc = dx.LocalHeight();
    for( Int iLoc=0; iLoc<nLoc; ++iLoc )
    {
        const Real xi = x.GetLocal(iLoc,0);
        const Real si = s.GetLocal(iLoc,0);
        const Real dsi = ds.GetLocal(iLoc,0);
        dx.SetLocal( iLoc, 0, -xi + (tau - dsi*xi)/si );
    }
}

template<typename Real>
Real IPFLineSearch
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& x,  const Matrix<Real>& l,  const Matrix<Real>& s,
  const Matrix<Real>& dx, const Matrix<Real>& dl, const Matrix<Real>& ds,
  Real gamma, Real beta, Real psi, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPFLineSearch"))
    if( gamma <= Real(0) || gamma >= Real(1) )
        LogicError("gamma must be in (0,1)");
    if( beta < Real(1) )
        LogicError("beta must be at least one");
    // TODO: Ensure the dimensions match
    if( b.Width() != 1 || c.Width() != 1 ||
        x.Width() != 1 || l.Width() != 1 || s.Width() != 1 ||
        dx.Width() != 1 || dl.Width() != 1 || ds.Width() != 1 )
        LogicError("{b,c,x,l,s,dx,dl,ds} must be column vectors");
    const Int m = A.Height();
    const Int n = A.Width();

    // Continually halve the step size until it is small enough so that the new
    // state lies within the "-\infty" neighborhood of the central path, i.e.,
    //  (a) || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu,
    //  (b) || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu,
    //  (c) x(\alpha), s(\alpha) > 0, and, for all i,
    //  (d) x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha),
    // where 
    //    x(\alpha) = x + \alpha dx,
    //    l(\alpha) = l + \alpha dl,
    //    s(\alpha) = s + \alpha ds,
    //    r_b(\alpha) = A x(\alpha) - b, and
    //    r_c(\alpha) = A^T l(\alpha) + s(\alpha) - c,
    // and the Armijo condition,
    //   \mu(\alpha) \le (1 - \alpha/\psi) \mu,
    // is also satisfied.
    // ===============================================
    // Setup
    // -----
    Matrix<Real> A_x, A_dx, AT_l, AT_dl, rb, rc;
    Zeros( A_x,   m, 1 );
    Zeros( A_dx,  m, 1 );
    Zeros( AT_l,  n, 1 );
    Zeros( AT_dl, n, 1 );
    Multiply( NORMAL,    Real(1), A, x,  Real(0), A_x   );
    Multiply( NORMAL,    Real(1), A, dx, Real(0), A_dx  );
    Multiply( TRANSPOSE, Real(1), A, l,  Real(0), AT_l  );
    Multiply( TRANSPOSE, Real(1), A, dl, Real(0), AT_dl );
    rb = A_x; 
    Axpy( Real(-1), b, rb );
    rc = AT_l;
    Axpy( Real(1), s, rc );
    Axpy( Real(-1), c, rc );
    const Real rbNrm2 = Nrm2( rb );
    const Real rcNrm2 = Nrm2( rc );
    const Real mu = Dot(x,s) / n;
    // Perform the line search using the cached data
    // ---------------------------------------------
    Real alpha = 1;
    Matrix<Real> x_alpha, s_alpha, rb_alpha, rc_alpha;
    for( Int k=0; k<100; ++k, alpha=alpha/2 )
    {
        // x(\alpha) = x + \alpha dx
        // ^^^^^^^^^^^^^^^^^^^^^^^^^
        x_alpha = x;
        Axpy( alpha, dx, x_alpha );

        // s(\alpha) = s + \alpha ds
        // ^^^^^^^^^^^^^^^^^^^^^^^^^ 
        s_alpha = s;
        Axpy( alpha, ds, s_alpha );

        // \mu(\alpha) = x(\alpha)^T s / n
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        const Real mu_alpha = Dot(x_alpha,s_alpha) / n;
        if( mu_alpha > (1-alpha/psi)*mu )
        {
            if( print )
                std::cout << "Armijo condition not satisfied" << std::endl;
            continue;
        }

        // Check 
        //    x(\alpha), s(\alpha) > 0, and 
        //    x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        bool balanced = true;
        for( Int i=0; i<n; ++i )
        {
            const Real xi_alpha = x_alpha.Get(i,0);
            const Real si_alpha = s_alpha.Get(i,0);
            if( xi_alpha <= Real(0) || si_alpha <= Real(0) )
                balanced = false;
            if( xi_alpha*si_alpha < gamma*mu_alpha )
                balanced = false;
        }
        if( !balanced )
        {
            if( print )
                std::cout << "  unbalanced entries" << std::endl;
            continue;
        }
        // Check || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rb_alpha = rb;
        Axpy( alpha, A_dx, rb_alpha );
        const Real rb_alphaNrm2 = Nrm2( rb_alpha );
        if( rb_alphaNrm2 > rbNrm2*beta*mu_alpha/mu )
        {
            if( print )
                std::cout << "  r_b failure: " << rb_alphaNrm2 << " > "
                          << rbNrm2*beta*mu_alpha/mu << std::endl;
            continue;
        }
        // Check || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rc_alpha = rc;
        Axpy( alpha, AT_dl, rc_alpha );
        Axpy( alpha, ds, rc_alpha );
        const Real rc_alphaNrm2 = Nrm2( rc_alpha );
        if( rc_alphaNrm2 > rcNrm2*beta*mu_alpha/mu )
        {
            if( print )
                std::cout << "  r_c failure: " << rc_alphaNrm2 << " > "
                          << rcNrm2*beta*mu_alpha/mu << std::endl;
        }
        else
            break;
    }
    return alpha;
}

template<typename Real>
Real IPFLineSearch
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, 
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& l, 
  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& dx, 
  const DistMultiVec<Real>& dl, 
  const DistMultiVec<Real>& ds,
  Real gamma, Real beta, Real psi, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPFLineSearch"))
    if( gamma <= Real(0) || gamma >= Real(1) )
        LogicError("gamma must be in (0,1)");
    if( beta < Real(1) )
        LogicError("beta must be at least one");
    // TODO: Ensure communicators are congruent
    // TODO: Ensure the dimensions match
    if( b.Width() != 1 || c.Width() != 1 ||
        x.Width() != 1 || l.Width() != 1 || s.Width() != 1 ||
        dx.Width() != 1 || dl.Width() != 1 || ds.Width() != 1 )
        LogicError("{b,c,x,l,s,dx,dl,ds} must be column vectors");
    const Int m = A.Height();
    const Int n = A.Width();
    const Int nLocal = x.LocalHeight();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);

    // Continually halve the step size until it is small enough so that the new
    // state lies within the "-\infty" neighborhood of the central path, i.e.,
    //  (a) || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu,
    //  (b) || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu,
    //  (c) x(\alpha), s(\alpha) > 0, and, for all i,
    //  (d) x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha),
    // where 
    //    x(\alpha) = x + \alpha dx,
    //    l(\alpha) = l + \alpha dl,
    //    s(\alpha) = s + \alpha ds,
    //    r_b(\alpha) = A x(\alpha) - b, and
    //    r_c(\alpha) = A^T l(\alpha) + s(\alpha) - c,
    // and the Armijo condition,
    //   \mu(\alpha) \le (1 - \alpha/\psi) \mu,
    // is also satisfied.
    // ===============================================
    // Setup
    // -----
    DistMultiVec<Real> A_x(comm), A_dx(comm), AT_l(comm), AT_dl(comm), 
                       rb(comm), rc(comm);
    Zeros( A_x,   m, 1 );
    Zeros( A_dx,  m, 1 );
    Zeros( AT_l,  n, 1 );
    Zeros( AT_dl, n, 1 );
    Multiply( NORMAL,    Real(1), A, x,  Real(0), A_x   );
    Multiply( NORMAL,    Real(1), A, dx, Real(0), A_dx  );
    Multiply( TRANSPOSE, Real(1), A, l,  Real(0), AT_l  );
    Multiply( TRANSPOSE, Real(1), A, dl, Real(0), AT_dl );
    rb = A_x; 
    Axpy( Real(-1), b, rb );
    rc = AT_l;
    Axpy( Real(1), s, rc );
    Axpy( Real(-1), c, rc );
    const Real rbNrm2 = Nrm2( rb );
    const Real rcNrm2 = Nrm2( rc );
    const Real mu = Dot(x,s) / n;
    // Perform the line search using the cached data
    // ---------------------------------------------
    Real alpha = 1;
    DistMultiVec<Real> x_alpha(comm), s_alpha(comm), 
                       rb_alpha(comm), rc_alpha(comm);
    for( Int k=0; k<100; ++k, alpha=alpha/2 )
    {
        // x(\alpha) = x + \alpha dx
        // ^^^^^^^^^^^^^^^^^^^^^^^^^
        x_alpha = x;
        Axpy( alpha, dx, x_alpha );

        // s(\alpha) = s + \alpha ds
        // ^^^^^^^^^^^^^^^^^^^^^^^^^ 
        s_alpha = s;
        Axpy( alpha, ds, s_alpha );

        // \mu(\alpha) = x(\alpha)^T s / n
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        const Real mu_alpha = Dot(x_alpha,s_alpha) / n;
        if( mu_alpha > (1-alpha/psi)*mu )
        {
            if( print && commRank == 0 )
                std::cout << "Armijo condition not satisfied" << std::endl;
            continue;
        }

        // Check 
        //    x(\alpha), s(\alpha) > 0, and 
        //    x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        byte locallyBalanced = true;
        for( Int iLoc=0; iLoc<nLocal; ++iLoc )
        {
            const Real xi_alpha = x_alpha.GetLocal(iLoc,0);
            const Real si_alpha = s_alpha.GetLocal(iLoc,0);
            if( xi_alpha <= Real(0) || si_alpha <= Real(0) )
                locallyBalanced = false;
            if( xi_alpha*si_alpha < gamma*mu_alpha )
                locallyBalanced = false;
        }
        const byte balanced = 
            mpi::AllReduce( locallyBalanced, mpi::BINARY_AND, comm ); 
        if( !balanced )
        {
            if( print && commRank == 0 )
                std::cout << "  unbalanced entries" << std::endl;
            continue;
        }
        // Check || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rb_alpha = rb;
        Axpy( alpha, A_dx, rb_alpha );
        const Real rb_alphaNrm2 = Nrm2( rb_alpha );
        if( rb_alphaNrm2 > rbNrm2*beta*mu_alpha/mu )
        {
            if( print && commRank == 0 )
                std::cout << "  r_b failure: " << rb_alphaNrm2 << " > "
                          << rbNrm2*beta*mu_alpha/mu << std::endl;
            continue;
        }
        // Check || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rc_alpha = rc;
        Axpy( alpha, AT_dl, rc_alpha );
        Axpy( alpha, ds, rc_alpha );
        const Real rc_alphaNrm2 = Nrm2( rc_alpha );
        if( rc_alphaNrm2 > rcNrm2*beta*mu_alpha/mu )
        {
            if( print && commRank == 0 )
                std::cout << "  r_c failure: " << rc_alphaNrm2 << " > "
                          << rcNrm2*beta*mu_alpha/mu << std::endl;
        }
        else
            break;
    }
    return alpha;
}

template<typename Real>
void IPF
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,  const Matrix<Real>& c,
  Matrix<Real>& x, Matrix<Real>& l, Matrix<Real>& s,
  Real muTol, Real rbTol, Real rcTol, Int maxIts,
  Real sigma, Real gamma, Real beta, Real psi, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPF"))    
    const Int n = A.Width();
    SparseMatrix<Real> J;
    Matrix<Real> y, rb, rc, dx, dl, ds;
    for( Int numIts=0; numIts<maxIts; ++numIts )
    {
        // Check for convergence
        // =====================
        // mu = x^T s / n
        // ---------------
        const Real mu = Dot(x,s) / n;
        // || r_b ||_2 = || A x - b ||_2
        // -----------------------------
        rb = b;
        Multiply( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        // || r_c ||_2 = || A^T l + s - c ||_2
        // -----------------------------------
        rc = c;
        Multiply( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
        const Real rcNrm2 = Nrm2( rc );
        // Now check the pieces
        // --------------------
        if( mu <= muTol && rbNrm2 <= rbTol && rcNrm2 <= rcTol )
            break;
        else if( print )
            std::cout << " iter " << numIts << ":\n"
                      << "  || r_b ||_2  = " << rbNrm2
                      << "  || r_c ||_2  = " << rcNrm2
                      << "  mu = " << mu << std::endl;

        // Construct the reduced KKT system, J dl = y
        // ==========================================
        FormNormalSystem( A, b, c, x, l, s, sigma*mu, J, y );

        // Compute the proposed step from the KKT system
        // =============================================
        SolveNormalSystem( A, b, c, x, l, s, sigma*mu, J, y, dx, dl, ds );

#ifndef RELEASE
          // Sanity checks
          // =============
          Matrix<Real> dsError, dxError, dlError;

          dsError.Resize( n, 1 );
          for( Int i=0; i<n; ++i )
          {
              const Real xi = x.Get(i,0);
              const Real si = s.Get(i,0);
              dsError.Set( i, 0, xi*si - sigma*mu );
          }
          const Real rmuNrm2 = Nrm2( dsError );
          for( Int i=0; i<n; ++i )
          {
              const Real xi = x.Get(i,0);
              const Real si = s.Get(i,0);
              const Real dxi = dx.Get(i,0);
              const Real dsi = ds.Get(i,0);
              dsError.Update( i, 0, xi*dsi + si*dxi );
          }
          const Real dsErrorNrm2 = Nrm2( dsError );

          dlError = ds;
          Multiply( TRANSPOSE, Real(1), A, dl, Real(1), dlError );
          Axpy( Real(1), rc, dlError );
          const Real dlErrorNrm2 = Nrm2( dlError );

          dxError = rb;
          Multiply( NORMAL, Real(1), A, dx, Real(1), dxError );
          const Real dxErrorNrm2 = Nrm2( dxError );
  
          if( print )
              std::cout << "  || dsError ||_2 / || r_mu ||_2 = " 
                        << dsErrorNrm2/rmuNrm2 << "\n"
                        << "  || dxError ||_2 / || r_c ||_2 = " 
                        << dxErrorNrm2/rcNrm2 << "\n"
                        << "  || dlError ||_2 / || r_b ||_2 = " 
                        << dlErrorNrm2/rbNrm2 << std::endl;
#endif

        // Decide on the step length
        // =========================
        const Real alpha =
          IPFLineSearch
          ( A, b, c, x, l, s, dx, dl, ds, gamma, beta, psi, print );
        if( print )
            std::cout << "  alpha = " << alpha << std::endl;

        // Update the state by stepping a distance of alpha
        // ================================================
        Axpy( alpha, dx, x );
        Axpy( alpha, dl, l );
        Axpy( alpha, ds, s );
    }
}

template<typename Real>
void IPF
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b,  const DistMultiVec<Real>& c,
  DistMultiVec<Real>& x, DistMultiVec<Real>& l, DistMultiVec<Real>& s,
  Real muTol, Real rbTol, Real rcTol, Int maxIts,
  Real sigma, Real gamma, Real beta, Real psi, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPF"))    
    // TODO: Input checks
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);

    DistSparseMatrix<Real> J(comm);
    DistMultiVec<Real> y(comm), rb(comm), rc(comm), 
                       dx(comm), dl(comm), ds(comm);
    for( Int numIts=0; numIts<maxIts; ++numIts )
    {
        // Check for convergence
        // =====================
        // mu = x^T s / n
        // ---------------
        const Real mu = Dot(x,s) / n;
        // || r_b ||_2 = || A x - b ||_2
        // -----------------------------
        rb = b; 
        Multiply( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        // || r_c ||_2 = || A^T l + s - c ||_2
        // -----------------------------------
        rc = c;
        Multiply( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
        const Real rcNrm2 = Nrm2( rc );
        // Now check the pieces
        // --------------------
        if( mu <= muTol && rbNrm2 <= rbTol && rcNrm2 <= rcTol )
            break;
        else if( print && commRank == 0 )
            std::cout << "  iter " << numIts << ": \n"
                      << "  || r_b ||_2 = " << rbNrm2 
                      << "  || r_c ||_2 = " << rcNrm2 
                      << "  mu = " << mu << std::endl;

        // Construct the reduced KKT system, J dl = y
        // ==========================================
        FormNormalSystem( A, b, c, x, l, s, sigma*mu, J, y );
       
        // Compute the proposed step from the KKT system
        // =============================================
        SolveNormalSystem( A, b, c, x, l, s, sigma*mu, J, y, dx, dl, ds );

#ifndef RELEASE
          // Sanity checks
          // =============
          DistMultiVec<Real> dsError(comm), dxError(comm), dlError(comm);

          dsError.Resize( n, 1 );
          for( Int iLoc=0; iLoc<dsError.LocalHeight(); ++iLoc )
          {
              const Real xi = x.GetLocal(iLoc,0);
              const Real si = s.GetLocal(iLoc,0);
              dsError.SetLocal( iLoc, 0, xi*si - sigma*mu );
          }
          const Real rmuNrm2 = Nrm2( dsError );
          for( Int iLoc=0; iLoc<dsError.LocalHeight(); ++iLoc )
          {
              const Real xi = x.GetLocal(iLoc,0);
              const Real si = s.GetLocal(iLoc,0);
              const Real dxi = dx.GetLocal(iLoc,0);
              const Real dsi = ds.GetLocal(iLoc,0);
              dsError.UpdateLocal( iLoc, 0, xi*dsi + si*dxi );
          }
          const Real dsErrorNrm2 = Nrm2( dsError );

          dlError = ds;
          Multiply( TRANSPOSE, Real(1), A, dl, Real(1), dlError );
          Axpy( Real(1), rc, dlError );
          const Real dlErrorNrm2 = Nrm2( dlError );

          dxError = rb;
          Multiply( NORMAL, Real(1), A, dx, Real(1), dxError );
          const Real dxErrorNrm2 = Nrm2( dxError );

          if( print && commRank == 0 )
              std::cout << "  || dsError ||_2 / || r_mu ||_2 = " 
                        << dsErrorNrm2/rmuNrm2 << "\n"
                        << "  || dxError ||_2 / || r_c ||_2 = " 
                        << dxErrorNrm2/rcNrm2 << "\n"
                        << "  || dlError ||_2 / || r_b ||_2 = " 
                        << dlErrorNrm2/rbNrm2 << std::endl;
#endif

        // Decide on the step length
        // =========================
        const Real alpha = 
          IPFLineSearch 
          ( A, b, c, x, l, s, dx, dl, ds, gamma, beta, psi, print );
        if( print && commRank == 0 )
            std::cout << "  alpha = " << alpha << std::endl;

        // Update the state by stepping a distance of alpha
        // ================================================
        Axpy( alpha, dx, x );
        Axpy( alpha, dl, l );
        Axpy( alpha, ds, s );
    }
}

} // namespace lin_prog

// These implementations are adaptations of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/linprog/linprog.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This ADMM attempts to solve the following linear program:
//     minimize    c' x
//     subject to  A x = b, x >= 0
//

template<typename Real>
Int LinearProgram
( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, 
  Matrix<Real>& z,
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("LinearProgram"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");
    
    // Cache a custom partially-pivoted LU factorization of 
    //    |  rho*I   A^H | = | B11  B12 |
    //    |  A       0   |   | B21  B22 |
    // by (justifiably) avoiding pivoting in the first n steps of
    // the factorization, so that
    //    [I,rho*I] = lu(rho*I).
    // The factorization would then proceed with 
    //    B21 := B21 U11^{-1} = A (rho*I)^{-1} = A/rho
    //    B12 := L11^{-1} B12 = I A^H = A^H.
    // The Schur complement would then be
    //    B22 := B22 - B21 B12 = 0 - (A*A^H)/rho.
    // We then factor said matrix with LU with partial pivoting and
    // swap the necessary rows of B21 in order to implicitly commute
    // the row pivots with the Gauss transforms in the manner standard
    // for GEPP. Unless A A' is singular, pivoting should not be needed,
    // as Cholesky factorization of the negative matrix should be valid.
    //
    // The result is the factorization
    //   | I 0   | | rho*I A^H | = | I   0   | | rho*I U12 |,
    //   | 0 P22 | | A     0   |   | L21 L22 | | 0     U22 |
    // where [L22,U22] are stored within B22.
    Matrix<Real> U12, L21, B22, bPiv;
    Adjoint( A, U12 );
    L21 = A; Scale( 1/rho, L21 );
    Herk( LOWER, NORMAL, -1/rho, A, B22 );
    MakeHermitian( LOWER, B22 );
    // TODO: Replace with sparse-direct Cholesky version?
    Matrix<Int> perm2;
    LU( B22, perm2 );
    PermuteRows( L21, perm2 );
    bPiv = b;
    PermuteRows( bPiv, perm2 );

    // Possibly form the inverse of L22 U22
    Matrix<Real> X22;
    if( inv )
    {
        X22 = B22;
        MakeTrapezoidal( LOWER, X22 );
        SetDiagonal( X22, Real(1) );
        TriangularInverse( LOWER, UNIT, X22 );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, Real(1), B22, X22 );
    }

    Int numIter=0;
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Real> g, xTmp, y, t;
    Zeros( g, m+n, 1 );
    PartitionDown( g, xTmp, y, n );
    Matrix<Real> x, u, zOld, xHat;
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    Zeros( t, n, 1 );
    while( numIter < maxIter )
    {
        zOld = z;

        // Find x from
        //  | rho*I  A^H | | x | = | rho*(z-u)-c | 
        //  | A      0   | | y |   | b           |
        // via our cached custom factorization:
        // 
        // |x| = inv(U) inv(L) P' |rho*(z-u)-c|
        // |y|                    |b          |
        //     = |rho*I U12|^{-1} |I   0  | |I 0   | |rho*(z-u)-c|
        //     = |0     U22|      |L21 L22| |0 P22'| |b          |
        //     = "                        " |rho*(z-u)-c|
        //                                  | P22' b    |
        xTmp = z;
        Axpy( Real(-1), u, xTmp );
        Scale( rho, xTmp );
        Axpy( Real(-1), c, xTmp );
        y = bPiv;
        Gemv( NORMAL, Real(-1), L21, xTmp, Real(1), y );
        if( inv )
        {
            Gemv( NORMAL, Real(1), X22, y, t );
            y = t;
        }
        else
        {
            Trsv( LOWER, NORMAL, UNIT, B22, y );
            Trsv( UPPER, NORMAL, NON_UNIT, B22, y );
        }
        Gemv( NORMAL, Real(-1), U12, y, Real(1), xTmp );
        Scale( 1/rho, xTmp );

        // xHat := alpha*x + (1-alpha)*zOld
        xHat = xTmp;
        Scale( alpha, xHat );
        Axpy( 1-alpha, zOld, xHat );

        // z := pos(xHat+u)
        z = xHat;
        Axpy( Real(1), u, z );
        LowerClip( z, Real(0) );

        // u := u + (xHat-z)
        Axpy( Real(1),  xHat, u );
        Axpy( Real(-1), z,    u );

        const Real objective = Dot( c, xTmp );

        // rNorm := || x - z ||_2
        t = xTmp;
        Axpy( Real(-1), z, t );
        const Real rNorm = FrobeniusNorm( t );
        // sNorm := |rho| || z - zOld ||_2
        t = z;
        Axpy( Real(-1), zOld, t );
        const Real sNorm = Abs(rho)*FrobeniusNorm( t );

        const Real epsPri = Sqrt(Real(n))*absTol +
            relTol*Max(FrobeniusNorm(xTmp),FrobeniusNorm(z));
        const Real epsDual = Sqrt(Real(n))*absTol +
            relTol*Abs(rho)*FrobeniusNorm(u);

        if( progress )
        {
            t = xTmp;
            LowerClip( t, Real(0) );
            Axpy( Real(-1), xTmp, t );
            const Real clipDist = FrobeniusNorm( t );
            std::cout << numIter << ": "
              << "||x-z||_2=" << rNorm << ", "
              << "epsPri=" << epsPri << ", "
              << "|rho| ||z-zOld||_2=" << sNorm << ", "
              << "epsDual=" << epsDual << ", "
              << "||x-Pos(x)||_2=" << clipDist << ", "
              << "c'x=" << objective << std::endl;
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( maxIter == numIter )
        std::cout << "ADMM failed to converge" << std::endl;
    x = xTmp;
    return numIter;
}

template<typename Real>
Int LinearProgram
( const AbstractDistMatrix<Real>& APre, const AbstractDistMatrix<Real>& bPre,
  const AbstractDistMatrix<Real>& cPre,       AbstractDistMatrix<Real>& zPre, 
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("LinearProgram"))

    auto APtr = ReadProxy<Real,MC,MR>( &APre );  auto& A = *APtr;
    auto bPtr = ReadProxy<Real,MC,MR>( &bPre );  auto& b = *bPtr;
    auto cPtr = ReadProxy<Real,MC,MR>( &cPre );  auto& c = *cPtr;
    auto zPtr = WriteProxy<Real,MC,MR>( &zPre ); auto& z = *zPtr;

    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");
    
    // Cache a custom partially-pivoted LU factorization of 
    //    |  rho*I   A^H | = | B11  B12 |
    //    |  A       0   |   | B21  B22 |
    // by (justifiably) avoiding pivoting in the first n steps of
    // the factorization, so that
    //    [I,rho*I] = lu(rho*I).
    // The factorization would then proceed with 
    //    B21 := B21 U11^{-1} = A (rho*I)^{-1} = A/rho
    //    B12 := L11^{-1} B12 = I A^H = A^H.
    // The Schur complement would then be
    //    B22 := B22 - B21 B12 = 0 - (A*A^H)/rho.
    // We then factor said matrix with LU with partial pivoting and
    // swap the necessary rows of B21 in order to implicitly commute
    // the row pivots with the Gauss transforms in the manner standard
    // for GEPP. Unless A A' is singular, pivoting should not be needed,
    // as Cholesky factorization of the negative matrix should be valid.
    //
    // The result is the factorization
    //   | I 0   | | rho*I A^H | = | I   0   | | rho*I U12 |,
    //   | 0 P22 | | A     0   |   | L21 L22 | | 0     U22 |
    // where [L22,U22] are stored within B22.
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    DistMatrix<Real> U12(grid), L21(grid), B22(grid), bPiv(grid);
    U12.Align( 0,                 n%U12.RowStride() );
    L21.Align( n%L21.ColStride(), 0                 );
    B22.Align( n%B22.ColStride(), n%B22.RowStride() );
    Adjoint( A, U12 );
    L21 = A; Scale( 1/rho, L21 );
    Herk( LOWER, NORMAL, -1/rho, A, B22 );
    MakeHermitian( LOWER, B22 );
    DistMatrix<Int,VC,STAR> perm2(grid);
    LU( B22, perm2 );
    PermuteRows( L21, perm2 );
    bPiv = b;
    PermuteRows( bPiv, perm2 );

    // Possibly form the inverse of L22 U22
    DistMatrix<Real> X22(grid);
    if( inv )
    {
        X22 = B22;
        MakeTrapezoidal( LOWER, X22 );
        SetDiagonal( X22, Real(1) );
        TriangularInverse( LOWER, UNIT, X22 );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, Real(1), B22, X22 );
    }

    Int numIter=0;
    DistMatrix<Real> g(grid), xTmp(grid), y(grid), t(grid);
    Zeros( g, m+n, 1 );
    PartitionDown( g, xTmp, y, n );
    DistMatrix<Real> x(grid), u(grid), zOld(grid), xHat(grid);
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    Zeros( t, n, 1 );
    while( numIter < maxIter )
    {
        zOld = z;

        // Find x from
        //  | rho*I  A^H | | x | = | rho*(z-u)-c | 
        //  | A      0   | | y |   | b           |
        // via our cached custom factorization:
        // 
        // |x| = inv(U) inv(L) P' |rho*(z-u)-c|
        // |y|                    |b          |
        //     = |rho*I U12|^{-1} |I   0  | |I 0   | |rho*(z-u)-c|
        //     = |0     U22|      |L21 L22| |0 P22'| |b          |
        //     = "                        " |rho*(z-u)-c|
        //                                  | P22' b    |
        xTmp = z;
        Axpy( Real(-1), u, xTmp );
        Scale( rho, xTmp );
        Axpy( Real(-1), c, xTmp );
        y = bPiv;
        Gemv( NORMAL, Real(-1), L21, xTmp, Real(1), y );
        if( inv )
        {
            Gemv( NORMAL, Real(1), X22, y, t );
            y = t;
        }
        else
        {
            Trsv( LOWER, NORMAL, UNIT, B22, y );
            Trsv( UPPER, NORMAL, NON_UNIT, B22, y );
        }
        Gemv( NORMAL, Real(-1), U12, y, Real(1), xTmp );
        Scale( 1/rho, xTmp );

        // xHat := alpha*x + (1-alpha)*zOld
        xHat = xTmp;
        Scale( alpha, xHat );
        Axpy( 1-alpha, zOld, xHat );

        // z := pos(xHat+u)
        z = xHat;
        Axpy( Real(1), u, z );
        LowerClip( z, Real(0) );

        // u := u + (xHat-z)
        Axpy( Real(1),  xHat, u );
        Axpy( Real(-1), z,    u );

        const Real objective = Dot( c, xTmp );

        // rNorm := || x - z ||_2
        t = xTmp;
        Axpy( Real(-1), z, t );
        const Real rNorm = FrobeniusNorm( t );
        // sNorm := |rho| || z - zOld ||_2
        t = z;
        Axpy( Real(-1), zOld, t );
        const Real sNorm = Abs(rho)*FrobeniusNorm( t );

        const Real epsPri = Sqrt(Real(n))*absTol +
            relTol*Max(FrobeniusNorm(xTmp),FrobeniusNorm(z));
        const Real epsDual = Sqrt(Real(n))*absTol +
            relTol*Abs(rho)*FrobeniusNorm(u);

        if( progress )
        {
            t = xTmp;
            LowerClip( t, Real(0) );
            Axpy( Real(-1), xTmp, t );
            const Real clipDist = FrobeniusNorm( t );
            if( grid.Rank() == 0 )
                std::cout << numIter << ": "
                  << "||x-z||_2=" << rNorm << ", "
                  << "epsPri=" << epsPri << ", "
                  << "|rho| ||z-zOld||_2=" << sNorm << ", "
                  << "epsDual=" << epsDual << ", "
                  << "||x-Pos(x)||_2=" << clipDist << ", "
                  << "c'x=" << objective << std::endl;
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( maxIter == numIter && grid.Rank() == 0 )
        std::cout << "ADMM failed to converge" << std::endl;
    x = xTmp;
    return numIter;
}

#define PROTO(Real) \
  template void lin_prog::FormAugmentedSystem \
  ( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& x, const Matrix<Real>& l, const Matrix<Real>& s, \
    Real tau, Matrix<Real>& J, Matrix<Real>& y ); \
  template void lin_prog::FormAugmentedSystem \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& l, \
    const AbstractDistMatrix<Real>& s, \
    Real tau, AbstractDistMatrix<Real>& J, AbstractDistMatrix<Real>& y ); \
  template void lin_prog::FormAugmentedSystem \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& x, const Matrix<Real>& l, const Matrix<Real>& s, \
    Real tau, SparseMatrix<Real>& J, Matrix<Real>& y ); \
  template void lin_prog::FormAugmentedSystem \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& x, const DistMultiVec<Real>& l, \
    const DistMultiVec<Real>& s, \
    Real tau, DistSparseMatrix<Real>& J, DistMultiVec<Real>& y ); \
  template void lin_prog::FormNormalSystem \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& x, const Matrix<Real>& l, const Matrix<Real>& s, \
    Real tau, SparseMatrix<Real>& J, Matrix<Real>& y ); \
  template void lin_prog::SolveNormalSystem \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& x, const Matrix<Real>& l, const Matrix<Real>& s, \
    Real tau, const SparseMatrix<Real>& J, const Matrix<Real>& y, \
    Matrix<Real>& dx, Matrix<Real>& dl, Matrix<Real>& ds ); \
  template void lin_prog::FormNormalSystem \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& x, const DistMultiVec<Real>& l, \
    const DistMultiVec<Real>& s, \
    Real tau, DistSparseMatrix<Real>& J, DistMultiVec<Real>& y ); \
  template void lin_prog::SolveNormalSystem \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& x, const DistMultiVec<Real>& l, \
    const DistMultiVec<Real>& s, \
    Real tau, const DistSparseMatrix<Real>& J, const DistMultiVec<Real>& y, \
    DistMultiVec<Real>& dx, DistMultiVec<Real>& dl, DistMultiVec<Real>& ds ); \
  template Real lin_prog::IPFLineSearch \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b,  const Matrix<Real>& c, \
    const Matrix<Real>& x,  const Matrix<Real>& l,  const Matrix<Real>& s, \
    const Matrix<Real>& dx, const Matrix<Real>& dl, const Matrix<Real>& ds, \
    Real gamma, Real beta, Real psi, bool print ); \
  template Real lin_prog::IPFLineSearch \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b,  const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& x,  const DistMultiVec<Real>& l, \
    const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& dx, const DistMultiVec<Real>& dl, \
    const DistMultiVec<Real>& ds, \
    Real gamma, Real beta, Real psi, bool print ); \
  template void lin_prog::IPF \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b,  const Matrix<Real>& c, \
    Matrix<Real>& x, Matrix<Real>& l, Matrix<Real>& s, \
    Real muTol, Real rbTol, Real rcTol, Int maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print ); \
  template void lin_prog::IPF \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b,  const DistMultiVec<Real>& c, \
    DistMultiVec<Real>& x, DistMultiVec<Real>& l, DistMultiVec<Real>& s, \
    Real muTol, Real rbTol, Real rcTol, Int maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print ); \
  template Int LinearProgram \
  ( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, \
    Matrix<Real>& z, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress ); \
  template Int LinearProgram \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
    const AbstractDistMatrix<Real>& c,       AbstractDistMatrix<Real>& z, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
