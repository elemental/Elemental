/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F> 
void LeastSquares
( Orientation orientation, Matrix<F>& A, const Matrix<F>& B, 
  Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LeastSquares"))

    Matrix<F> t;
    Matrix<Base<F>> d;

    const Int m = A.Height();
    const Int n = A.Width();
    if( m >= n )
    {
        QR( A, t, d );
        qr::SolveAfter( orientation, A, t, d, B, X );
    }
    else
    {
        LQ( A, t, d );
        lq::SolveAfter( orientation, A, t, d, B, X );
    }
}

template<typename F> 
void LeastSquares
( Orientation orientation, AbstractDistMatrix<F>& APre,
  const AbstractDistMatrix<F>& B, AbstractDistMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LeastSquares"))

    auto APtr = ReadProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());

    const Int m = A.Height();
    const Int n = A.Width();
    if( m >= n )
    {
        QR( A, t, d );
        qr::SolveAfter( orientation, A, t, d, B, X );
    }
    else
    {
        LQ( A, t, d );
        lq::SolveAfter( orientation, A, t, d, B, X );
    }
}

// The following routines solve either
//
//   Minimum length: min || X ||_F s.t. op(A) X = B, or
//   Least squares:  min || op(A) X - B ||_F,
//
// where op(A) is either A, A^T, or A^H, via forming a Hermitian 
// quasi-semidefinite system J D = \hat{B}, where J is 
//
//    | I   A | when height(A) >= width(A), or
//    | A^H 0 |
//     
//    | I A^H | when height(A) < width(A).
//    | A  0  |
//
// When height(op(A)) < width(op(A)), the system
//
//     | I     op(A)^H | | X | = | 0 |
//     | op(A)    0    | | Y |   | B |
//
// guarantees that op(A) X = B and X is in range(op(A)^H), which shows that 
// X solves the minimum length problem. Otherwise, the system
//
//     | I        op(A) | | R | = | B |
//     | op(A)^H    0   | | X |   | 0 |
// 
// guarantees that R = B - op(A) X and R in null(op(A)^H), which is equivalent
// to solving min || op(A) X - B ||_F.
//
// The Hermitian quasi-semidefinite systems are solved by converting them into
// Hermitian quasi-definite form via a priori regularization, applying an 
// LDL^H factorization with static pivoting to the regularized system, and
// using the iteratively-refined solution of with the regularized factorization
// as a preconditioner for the original problem (defaulting to Flexible GMRES
// for now).
//
// This approach originated within 
//
//    Michael Saunders, 
//   "Chapter 8, Cholesky-based Methods for Sparse Least Squares:
//    The Benefits of Regularization",
//    in L. Adams and J.L. Nazareth (eds.), Linear and Nonlinear Conjugate
//    Gradient-Related Methods, SIAM, Philadelphia, 92--100 (1996).
//
// But note that SymmLQ and LSQR were used rather than flexible GMRES, and 
// iteratively refining *within* the preconditioner was not discussed.
//

template<typename F>
void LeastSquares
( Orientation orientation,
  const SparseMatrix<F>& A, const Matrix<F>& B, Matrix<F>& X,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(
      CallStackEntry cse("LeastSquares");
      if( orientation == NORMAL && A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
      if( orientation != NORMAL && A.Width() != B.Height() )
          LogicError("Width of A and height of B must match");
    )
    typedef Base<F> Real;

    // ABar := inv(D_r) op(A) inv(D_c)
    // ===============================
    SparseMatrix<F> ABar;
    if( orientation == NORMAL )
        ABar = A;
    else if( orientation == TRANSPOSE )
        Transpose( A, ABar );
    else
        Adjoint( A, ABar );
    const Int m = ABar.Height();
    const Int n = ABar.Width();
    const Int numEntriesA = ABar.NumEntries();
    const bool progress = false;
    Matrix<Real> dRow, dCol;
    GeomEquil( ABar, dRow, dCol, progress );

    // BBar := inv(D_r) B
    // ==================
    auto BBar = B;
    DiagonalSolve( LEFT, NORMAL, dRow, BBar );
    const Int k = BBar.Width();

    SparseMatrix<F> J;
    Zeros( J, m+n, m+n );
    J.Reserve( 2*numEntriesA + Max(m,n) );
    if( m >= n )
    {
        // Form J = [D_r^{-2}, ABar; ABar^H, 0]
        // =============================
        for( Int e=0; e<numEntriesA; ++e )
        {
            J.QueueUpdate( ABar.Row(e),   ABar.Col(e)+m,      ABar.Value(e)  );
            J.QueueUpdate( ABar.Col(e)+m, ABar.Row(e),   Conj(ABar.Value(e)) );
        }
        for( Int e=0; e<m; ++e )
            J.QueueUpdate( e, e, Pow(dRow.Get(e,0),Real(-2)) );
    }
    else
    {
        // Form J = [D_c^{-2}, ABar^H; ABar, 0]
        // ====================================
        for( Int e=0; e<numEntriesA; ++e )
        {
            J.QueueUpdate( ABar.Col(e),   ABar.Row(e)+n, Conj(ABar.Value(e)) );
            J.QueueUpdate( ABar.Row(e)+n, ABar.Col(e),        ABar.Value(e)  );
        }
        for( Int e=0; e<n; ++e )
            J.QueueUpdate( e, e, Pow(dCol.Get(e,0),Real(-2)) );
    }
    J.MakeConsistent();

    Matrix<F> D;
    Zeros( D, m+n, k );
    if( m >= n )
    {
        // Form D = [BBar; 0]
        // ==================
        auto DT = D( IR(0,m), IR(0,k) );
        DT = BBar;
    }
    else
    {
        // Form D = [0; BBar]
        // ==================
        auto DB = D( IR(n,m+n), IR(0,k) );
        DB = BBar;
    }

    // Compute the regularized quasi-semidefinite fact of J
    // ====================================================
    const Real minReductionFactor = 2;
    const Int maxRefineIts = 50;
    bool print = false;
    const Real epsilon = lapack::MachineEpsilon<Real>();
    bool aPriori = true;
    const Real regMagPrimal = Pow(epsilon,Real(0.5));
    const Real regMagDual = Pow(epsilon,Real(0.5));
    const Real pivTol = MaxNorm(J)*epsilon;

    Matrix<Real> regCand, reg;
    regCand.Resize( m+n, 1 );
    for( Int i=0; i<Max(m,n); ++i )
        regCand.Set( i, 0, regMagPrimal );
    for( Int i=Max(m,n); i<m+n; ++i )
        regCand.Set( i, 0, -regMagDual );
    Zeros( reg, m+n, 1 );

    vector<Int> map, invMap;
    SymmNodeInfo info;
    Separator rootSep;
    NestedDissection( J.LockedGraph(), map, rootSep, info );
    InvertMap( map, invMap );
    SymmFront<F> JFront( J, map, info );

    MatrixNode<Real> regCandNodal(invMap,info,regCand),
                     regNodal(invMap,info,reg);
    RegularizedQSDLDL
    ( info, JFront, pivTol, regCandNodal, regNodal, aPriori );
    regNodal.Push( invMap, info, reg );

    // Successively solve each of the k linear systems
    // ===============================================
    // TODO: Extend the iterative refinement to handle multiple RHS
    Matrix<F> u;
    Zeros( u, m+n, 1 );
    for( Int j=0; j<k; ++j )
    {
        auto d = D( IR(0,m+n), IR(j,j+1) );
        u = d;
        reg_qsd_ldl::SolveAfter
        ( J, reg, invMap, info, JFront, u,
          REG_REFINE_FGMRES,
          minReductionFactor, maxRefineIts, print );
        d = u;
    }

    Zeros( X, n, k );
    if( m >= n )
    {
        // Extract XBar from [R; XBar]
        // ===========================
        auto DB = D( IR(m,m+n), IR(0,k) );
        X = DB;
    }
    else
    {
        // Extract XBar from [XBar; Y]
        // ===========================
        auto DT = D( IR(0,n), IR(0,k) );
        X = DT;
    }
    DiagonalSolve( LEFT, NORMAL, dCol, X );
}

template<typename F>
void LeastSquares
( Orientation orientation,
  const DistSparseMatrix<F>& A, const DistMultiVec<F>& B, DistMultiVec<F>& X,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(
      CallStackEntry cse("LeastSquares");
      if( orientation == NORMAL && A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
      if( orientation != NORMAL && A.Width() != B.Height() )
          LogicError("Width of A and height of B must match");
    )
    typedef Base<F> Real;
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size(comm);

    // ABar := inv(D_r) op(A) inv(D_c) 
    // ===============================
    DistSparseMatrix<F> ABar(comm);
    if( orientation == NORMAL )
        ABar = A;
    else if( orientation == TRANSPOSE )
        Transpose( A, ABar );
    else
        Adjoint( A, ABar );
    const Int m = ABar.Height();
    const Int n = ABar.Width();
    const bool progress = false;
    DistMultiVec<Real> dRow(comm), dCol(comm);
    GeomEquil( ABar, dRow, dCol, progress );

    // BBar := inv(D_r) B
    // ==================
    const Int k = B.Width();
    auto BBar = B;
    DiagonalSolve( LEFT, NORMAL, dRow, BBar );

    // Set J to [D_r^{-2}, ABar; ABar^H, 0] or [D_c^{-2}, ABar^H; ABar, 0] 
    // ===================================================================
    DistSparseMatrix<F> J(comm);
    Zeros( J, m+n, m+n );
    const Int numLocalEntriesA = ABar.NumLocalEntries();
    {
        // Compute metadata
        // ----------------
        vector<int> sendCounts(commSize,0);
        for( Int e=0; e<numLocalEntriesA; ++e )
        {
            const Int i = ABar.Row(e);
            const Int j = ABar.Col(e);
            if( m >= n )
            {
                // Sending ABar
                ++sendCounts[ J.RowOwner(i) ];
                // Sending ABar^H
                ++sendCounts[ J.RowOwner(j+m) ];
            }
            else
            {
                // Sending ABar
                ++sendCounts[ J.RowOwner(i+n) ];
                // Sending ABar^H
                ++sendCounts[ J.RowOwner(j) ];
            }
        }
        if( m >= n )
        {
            for( Int iLoc=0; iLoc<dRow.LocalHeight(); ++iLoc )
                ++sendCounts[ J.RowOwner(dRow.GlobalRow(iLoc)) ];
        }
        else
        {
            for( Int iLoc=0; iLoc<dCol.LocalHeight(); ++iLoc )
                ++sendCounts[ J.RowOwner(dCol.GlobalRow(iLoc)) ];
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        vector<int> sendOffsets(commSize), recvOffsets(commSize);
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack
        // ----
        vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
        vector<F> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int e=0; e<numLocalEntriesA; ++e )
        {
            const Int i = ABar.Row(e);
            const Int j = ABar.Col(e);
            const F value = ABar.Value(e);

            if( m >= n )
            {
                // Sending ABar
                int owner = J.RowOwner(i);
                sSendBuf[offsets[owner]] = i;
                tSendBuf[offsets[owner]] = j+m;
                vSendBuf[offsets[owner]] = value;
                ++offsets[owner];

                // Sending ABar^H
                owner = J.RowOwner(j+m);
                sSendBuf[offsets[owner]] = j+m;
                tSendBuf[offsets[owner]] = i;
                vSendBuf[offsets[owner]] = Conj(value);
                ++offsets[owner];
            }
            else
            {
                // Sending ABar
                int owner = J.RowOwner(i+n);
                sSendBuf[offsets[owner]] = i+n;
                tSendBuf[offsets[owner]] = j;
                vSendBuf[offsets[owner]] = value;
                ++offsets[owner];

                // Sending ABar^H
                owner = J.RowOwner(j);
                sSendBuf[offsets[owner]] = j;
                tSendBuf[offsets[owner]] = i+n;
                vSendBuf[offsets[owner]] = Conj(value);
                ++offsets[owner];
            }
        }
        if( m >= n )
        {
            for( Int iLoc=0; iLoc<dRow.LocalHeight(); ++iLoc )
            {
                const Int i = dRow.GlobalRow(iLoc);
                const int owner = J.RowOwner(i);
                sSendBuf[offsets[owner]] = i; 
                tSendBuf[offsets[owner]] = i;
                vSendBuf[offsets[owner]] = Pow(dRow.GetLocal(iLoc,0),Real(-2));
                ++offsets[owner];
            }
        }
        else
        {
            for( Int iLoc=0; iLoc<dCol.LocalHeight(); ++iLoc )
            {
                const Int i = dCol.GlobalRow(iLoc);
                const int owner = J.RowOwner(i);
                sSendBuf[offsets[owner]] = i; 
                tSendBuf[offsets[owner]] = i;
                vSendBuf[offsets[owner]] = Pow(dCol.GetLocal(iLoc,0),Real(-2));
                ++offsets[owner];
            }
        }

        // Exchange
        // --------
        vector<Int> sRecvBuf(totalRecv), tRecvBuf(totalRecv);
        vector<F> vRecvBuf(totalRecv);
        mpi::AllToAll
        ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( tSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          tRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        // Unpack
        // ------
        J.Reserve( totalRecv );
        for( Int e=0; e<totalRecv; ++e )
            J.QueueLocalUpdate
            ( sRecvBuf[e]-J.FirstLocalRow(), tRecvBuf[e], vRecvBuf[e] );
        J.MakeConsistent();
    }

    // Set D to [BBar; 0] or [0; BBar]
    // ===============================
    DistMultiVec<F> D(comm);
    Zeros( D, m+n, k );
    {
        // Compute metadata
        // ----------------
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<BBar.LocalHeight(); ++iLoc )
        {
            const Int i = BBar.GlobalRow(iLoc);
            if( m >= n )
                sendCounts[ D.RowOwner(i) ] += k;
            else
                sendCounts[ D.RowOwner(i+n) ] += k;
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        vector<int> sendOffsets(commSize), recvOffsets(commSize);
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack
        // ----
        vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
        vector<F> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int iLoc=0; iLoc<BBar.LocalHeight(); ++iLoc )
        {
            const Int i = BBar.GlobalRow(iLoc);

            if( m >= n )
            {
                int owner = D.RowOwner(i);
                for( Int j=0; j<k; ++j )
                {
                    const F value = BBar.GetLocal(iLoc,j);
                    sSendBuf[offsets[owner]] = i;
                    tSendBuf[offsets[owner]] = j;
                    vSendBuf[offsets[owner]] = value;
                    ++offsets[owner];
                }
            }
            else
            {
                int owner = D.RowOwner(i+n);
                for( Int j=0; j<k; ++j )
                {
                    const F value = BBar.GetLocal(iLoc,j);
                    sSendBuf[offsets[owner]] = i+n;
                    tSendBuf[offsets[owner]] = j;
                    vSendBuf[offsets[owner]] = value;
                    ++offsets[owner];
                }
            }
        }
        // Exchange
        // --------
        vector<Int> sRecvBuf(totalRecv), tRecvBuf(totalRecv);
        vector<F> vRecvBuf(totalRecv);
        mpi::AllToAll
        ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( tSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          tRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        // Unpack
        // ------
        for( Int e=0; e<totalRecv; ++e )
            D.UpdateLocal
            ( sRecvBuf[e]-D.FirstLocalRow(), tRecvBuf[e], vRecvBuf[e] );
    }

    // Compute the dynamically-regularized quasi-semidefinite fact of J
    // ================================================================
    const Real minReductionFactor = 2;
    const Int maxRefineIts = 50;
    bool print = false;
    const Real epsilon = lapack::MachineEpsilon<Real>();
    bool aPriori = true;
    const Real regMagPrimal = Pow(epsilon,Real(0.5));
    const Real regMagDual = Pow(epsilon,Real(0.5));
    const Real pivTol = MaxNorm(J)*epsilon;

    DistMultiVec<Real> regCand(comm), reg(comm);
    regCand.Resize( m+n, 1 );
    for( Int iLoc=0; iLoc<regCand.LocalHeight(); ++iLoc )
    {
        const Int i = regCand.GlobalRow(iLoc);
        if( i < Max(m,n) )
            regCand.SetLocal( iLoc, 0, regMagPrimal );
        else
            regCand.SetLocal( iLoc, 0, -regMagDual );
    }
    Zeros( reg, m+n, 1 );

    DistMap map, invMap;
    DistSymmNodeInfo info;
    DistSeparator rootSep;
    NestedDissection( J.LockedDistGraph(), map, rootSep, info );
    InvertMap( map, invMap );
    DistSymmFront<F> JFront( J, map, rootSep, info );

    DistMultiVecNode<Real> regCandNodal(invMap,info,regCand),
                           regNodal(invMap,info,reg);
    RegularizedQSDLDL
    ( info, JFront, pivTol, regCandNodal, regNodal, aPriori, LDL_1D );
    regNodal.Push( invMap, info, reg );

    // Successively solve each of the k linear systems
    // ===============================================
    // TODO: Extend the iterative refinement to handle multiple right-hand sides
    DistMultiVec<F> u(comm);
    Zeros( u, m+n, 1 );
    auto& DLoc = D.Matrix();
    auto& uLoc = u.Matrix();
    const Int DLocHeight = DLoc.Height();
    for( Int j=0; j<k; ++j )
    {
        auto dLoc = DLoc( IR(0,DLocHeight), IR(j,j+1) );
        Copy( dLoc, uLoc );
        reg_qsd_ldl::SolveAfter
        ( J, reg, invMap, info, JFront, u,
          REG_REFINE_FGMRES,
          minReductionFactor, maxRefineIts, print );
        Copy( uLoc, dLoc );
    }


    // Extract XBar from [R; XBar] or [XBar; Y] and then rescale
    // =========================================================
    Zeros( X, n, k );
    {
        // Compute metadata
        // ----------------
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<DLocHeight; ++iLoc )
        {
            const Int i = D.GlobalRow(iLoc);
            if( m > n )
            {
                if( i >= m )
                    sendCounts[ X.RowOwner(i-m) ] += k;
            }
            else
            {
                if( i < n )
                    sendCounts[ X.RowOwner(i) ] += k;
                else
                    break;
            }
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        vector<int> sendOffsets(commSize), recvOffsets(commSize);
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack
        // ----
        vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
        vector<F> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int iLoc=0; iLoc<DLocHeight; ++iLoc )
        {
            const Int i = D.GlobalRow(iLoc);
            if( m >= n )
            {
                if( i >= m )
                {
                    int owner = X.RowOwner(i-m);
                    for( Int j=0; j<k; ++j )
                    {
                        const F value = D.GetLocal(iLoc,j);
                        sSendBuf[offsets[owner]] = i-m;
                        tSendBuf[offsets[owner]] = j;
                        vSendBuf[offsets[owner]] = value;
                        ++offsets[owner];
                    }
                }
            }
            else
            {
                if( i < n )
                {
                    int owner = X.RowOwner(i);
                    for( Int j=0; j<k; ++j )
                    {
                        const F value = D.GetLocal(iLoc,j);
                        sSendBuf[offsets[owner]] = i;
                        tSendBuf[offsets[owner]] = j;
                        vSendBuf[offsets[owner]] = value;
                        ++offsets[owner];
                    }
                }
                else
                    break;
            }
        }
        // Exchange
        // --------
        vector<Int> sRecvBuf(totalRecv), tRecvBuf(totalRecv);
        vector<F> vRecvBuf(totalRecv);
        mpi::AllToAll
        ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( tSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          tRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        // Unpack
        // ------
        for( Int e=0; e<totalRecv; ++e )
            X.SetLocal
            ( sRecvBuf[e]-X.FirstLocalRow(), tRecvBuf[e], vRecvBuf[e] );
    }
    DiagonalSolve( LEFT, NORMAL, dCol, X );
}

#define PROTO(F) \
  template void LeastSquares \
  ( Orientation orientation, Matrix<F>& A, const Matrix<F>& B, \
    Matrix<F>& X ); \
  template void LeastSquares \
  ( Orientation orientation, AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& B, AbstractDistMatrix<F>& X ); \
  template void LeastSquares \
  ( Orientation orientation, \
    const SparseMatrix<F>& A, const Matrix<F>& B, \
    Matrix<F>& X, const BisectCtrl& ctrl ); \
  template void LeastSquares \
  ( Orientation orientation, \
    const DistSparseMatrix<F>& A, const DistMultiVec<F>& B, \
    DistMultiVec<F>& X, const BisectCtrl& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
