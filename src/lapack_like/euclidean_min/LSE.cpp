/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// This file implements both dense and sparse-direct solutions of 
// Equality-constrained Least Squares (LSE):
//
//     min_x || A x - c ||_2 subject to B x = d.
//
// For dense instances of the problem, a Generalized RQ factorization can be
// employed as long as A is m x n, B is p x n, and p <= n <= m+p. It is 
// assumed that B has full row rank, p, and [A;B] has full column rank, n.
//
// A Generalized RQ factorization of (B,A),
//    B = T Q = | 0 T12 | Q,  A = Z | R11 R12 | Q,
//                                  |   0 R22 |
// where Q and Z are unitary and R and T are upper-trapezoidal, allows us to
// re-express the constraint
//     T Q x = d,
// as
//     | 0 T12 | | y1 | = d,
//               | y2 |   
// where y = Q x, which only requires the solution of the upper-triangular 
// system
//     T12 y2 = d.
//
// The objective can be rewritten as
//     || A x - c ||_2 = || Z^H A x - Z^H c ||_2
//                     = ||   R Q x - Z^H c ||_2
// which, defining g = Z^H c, can be partitioned as
//     | R11 R12 | | y1 | - | g1 | = | R11 y1 + R12 y2 - g1 |.
//     |   0 R22 | | y2 |   | g2 |   |          R22 y2 - g2 |
// Since y2 is fixed by the constraint, the norm is minimized by setting the
// top term to zero, which involves solving the upper-triangular system
//     R11 y1 = g1 - R12 y2.
//       
// On exit of the internal dense implementation, A and B are overwritten with 
// their implicit Generalized RQ factorization of (B,A), and, optionally, C is 
// overwritten with the rotated residual matrix
//     Z^H (A X - C) = (R Q X - Z^H C) = |           0 |,
//                                       | R22 Y2 - G1 |
// where R22 is an upper-trapezoidal (not necessarily triangular) matrix.
// D is overwritten with arbitrary values.
//
// Note that essentially the same scheme is used in LAPACK's {S,D,C,Z}GGLSE.
//
// For sparse instances of the LSE problem, the symmetric quasi-semidefinite
// augmented system
//
//     | 0    A^H    B^H | |     x    |   | 0 |
//     | A -alpha*I   0  | | -r/alpha | = | c |
//     | B     0      0  | |  y/alpha |   | d |
//
// is formed, equilibrated, and then a priori regularization is added in order
// to make the system sufficiently quasi-definite. A Cholesky-like factorization
// of this regularized system is then used as a preconditioner for FGMRES(k).
//

namespace El {

namespace lse {

template<typename F> 
void Overwrite
( Matrix<F>& A, Matrix<F>& B, 
  Matrix<F>& C, Matrix<F>& D, 
  Matrix<F>& X, bool computeResidual )
{
    DEBUG_ONLY(CSE cse("lse::Overwrite"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = B.Height();
    const Int numRhs = D.Width();
    if( m != C.Height() )
        LogicError("A and C must be the same height");
    if( p != D.Height() )
        LogicError("B and D must be the same height");
    if( numRhs != C.Width() )
        LogicError("C and D must be the same width");
    if( n < p )
        LogicError("LSE requires width(A) >= height(B)");
    if( m+p < n )
        LogicError("LSE requires height(A)+height(B) >= width(A)");
    const bool checkIfSingular = true;

    // Compute the implicit Generalized RQ decomposition of (B,A)
    Matrix<F> tA, tB;
    Matrix<Base<F>> dA, dB;
    GRQ( B, tB, dB, A, tA, dA );

    // G := Z^H C
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, C );

    // Partition the relevant matrices
    Zeros( X, n, numRhs );
    Matrix<F> Y1, Y2;
    PartitionUp( X, Y1, Y2, p );
    Matrix<F> T11, T12;
    PartitionLeft( B, T11, T12, p );
    Matrix<F> R11, R12, R21, R22;
    PartitionDownDiagonal( A, R11, R12, R21, R22, n-p );
    Matrix<F> G1, G2;
    PartitionDown( C, G1, G2, n-p );

    // Solve T12 Y2 = D
    Y2 = D; 
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), T12, Y2, checkIfSingular );

    // G1 := G1 - R12 Y2
    Gemm( NORMAL, NORMAL, F(-1), R12, Y2, F(1), G1 );

    // Solve R11 Y1 = G1
    Y1 = G1;
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R11, Y1, checkIfSingular );

    if( computeResidual )
    {
        // R22 is upper-trapezoidal, and so it is best to decompose it in terms
        // of its upper-left triangular block and either its bottom zero 
        // block or right non-zero block. Putting k=Min(p,m-(n-p)), then
        // the k x k upper-left block is upper-triangular. If m >= n, the
        // bottom m-(n-p) - k = m-n rows are zero, otherwise the right 
        // p - k = n-m.columns are nonzero.
        if( m < n )
        {
            Matrix<F> R22L, R22R;
            PartitionLeft( R22, R22L, R22R, n-m );
            Matrix<F> DT, DB;
            PartitionUp( D, DT, DB, n-m );
            Gemm( NORMAL, NORMAL, F(-1), R22R, DB, F(1), G2 );
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R22L, DT );
            G2 -= DT;
        }
        else
        {
            Matrix<F> R22T, R22B;
            PartitionUp( R22, R22T, R22B, m-n );
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R22T, D );
            Matrix<F> G2T, G2B;
            PartitionUp( G2, G2T, G2B, m-n );
            G2T -= D;
        }
        Zero( G1 );
    }

    // X := Q^H Y
    rq::ApplyQ( LEFT, ADJOINT, B, tB, dB, X );
}

template<typename F> 
void Overwrite
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& BPre, 
  AbstractDistMatrix<F>& CPre, AbstractDistMatrix<F>& DPre, 
  AbstractDistMatrix<F>& XPre, bool computeResidual )
{
    DEBUG_ONLY(CSE cse("lse::Overwrite"))

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre ); auto& A = *APtr;
    auto BPtr = ReadWriteProxy<F,MC,MR>( &BPre ); auto& B = *BPtr;
    auto CPtr = ReadWriteProxy<F,MC,MR>( &CPre ); auto& C = *CPtr;
    auto DPtr = ReadWriteProxy<F,MC,MR>( &DPre ); auto& D = *DPtr;
    auto XPtr = WriteProxy<F,MC,MR>( &XPre );     auto& X = *XPtr;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = B.Height();
    const Int numRhs = D.Width();
    if( m != C.Height() )
        LogicError("A and C must be the same height");
    if( p != D.Height() )
        LogicError("B and D must be the same height");
    if( numRhs != C.Width() )
        LogicError("C and D must be the same width");
    if( n < p )
        LogicError("LSE requires width(A) >= height(B)");
    if( m+p < n )
        LogicError("LSE requires height(A)+height(B) >= width(A)");
    const Grid& g = A.Grid();
    if( g != B.Grid() || g != C.Grid() || g != D.Grid() )
        LogicError("All matrices must be distributed over the same grid");
    X.SetGrid( g );
    const bool checkIfSingular = true;

    // Compute the implicit Generalized RQ decomposition of (B,A)
    DistMatrix<F,MD,STAR> tA(g), tB(g);
    DistMatrix<Base<F>,MD,STAR> dA(g), dB(g);
    GRQ( B, tB, dB, A, tA, dA );

    // G := Z^H C
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, C );

    // Partition the relevant matrices
    Zeros( X, n, numRhs );
    DistMatrix<F> Y1(g), Y2(g);
    PartitionUp( X, Y1, Y2, p );
    DistMatrix<F> T11(g), T12(g);
    PartitionLeft( B, T11, T12, p );
    DistMatrix<F> R11(g), R12(g), R21(g), R22(g);
    PartitionDownDiagonal( A, R11, R12, R21, R22, n-p );
    DistMatrix<F> G1(g), G2(g);
    PartitionDown( C, G1, G2, n-p );

    // Solve T12 Y2 = D
    Y2 = D; 
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), T12, Y2, checkIfSingular );

    // G1 := G1 - R12 Y2
    Gemm( NORMAL, NORMAL, F(-1), R12, Y2, F(1), G1 );

    // Solve R11 Y1 = G1
    Y1 = G1;
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R11, Y1, checkIfSingular );

    if( computeResidual )
    {
        // R22 is upper-trapezoidal, and so it is best to decompose it in terms
        // of its upper-left triangular block and either its bottom zero 
        // block or right non-zero block. Putting k=Min(p,m-(n-p)), then
        // the k x k upper-left block is upper-triangular. If m >= n, the
        // bottom m-(n-p) - k = m-n rows are zero, otherwise the right 
        // p - k = n-m.columns are nonzero.
        if( m < n )
        {
            DistMatrix<F> R22L(g), R22R(g);
            PartitionLeft( R22, R22L, R22R, n-m );
            DistMatrix<F> DT(g), DB(g);
            PartitionUp( D, DT, DB, n-m );
            Gemm( NORMAL, NORMAL, F(-1), R22R, DB, F(1), G2 );
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R22L, DT );
            G2 -= DT;
        }
        else
        {
            DistMatrix<F> R22T(g), R22B(g);
            PartitionUp( R22, R22T, R22B, m-n );
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R22T, D );
            DistMatrix<F> G2T(g), G2B(g);
            PartitionUp( G2, G2T, G2B, m-n );
            G2T -= D;
        }
        Zero( G1 );
    }

    // X := Q^H Y
    rq::ApplyQ( LEFT, ADJOINT, B, tB, dB, X );
}

} // namespace lse

template<typename F> 
void LSE
( const Matrix<F>& A, const Matrix<F>& B, 
  const Matrix<F>& C, const Matrix<F>& D, 
        Matrix<F>& X )
{
    DEBUG_ONLY(CSE cse("LSE"))
    Matrix<F> ACopy( A ), BCopy( B ), CCopy( C ), DCopy( D );
    lse::Overwrite( ACopy, BCopy, CCopy, DCopy, X );
}

template<typename F> 
void LSE
( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B, 
  const AbstractDistMatrix<F>& C, const AbstractDistMatrix<F>& D, 
        AbstractDistMatrix<F>& X )
{
    DEBUG_ONLY(CSE cse("LSE"))
    DistMatrix<F> ACopy( A ), BCopy( B ), CCopy( C ), DCopy( D );
    lse::Overwrite( ACopy, BCopy, CCopy, DCopy, X );
}

template<typename F> 
void LSE
( const SparseMatrix<F>& A, const SparseMatrix<F>& B, 
  const Matrix<F>& C,       const Matrix<F>& D, 
        Matrix<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LSE"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Height();
    const Int numRHS = C.Width();

    // Rescale W = [ A; B ]
    // ====================
    Matrix<Real> dC;
    SparseMatrix<F> W;
    VCat( A, B, W );
    if( ctrl.equilibrate )
    {
        ColumnTwoNorms( W, dC );
        auto normMap = []( Real beta )
          { return ( beta < Sqrt(Epsilon<Real>()) ? Real(1) : beta ); };
        EntrywiseMap( dC, function<Real(Real)>(normMap) );
        DiagonalSolve( RIGHT, NORMAL, dC, W );
    }
    else
    {
        Ones( dC, n, 1 );
    }
    Real normScale = 1;
    if( ctrl.scaleTwoNorm )
    {
        normScale = TwoNormEstimate( W, ctrl.basisSize );
        if( ctrl.progress )
            cout << "Estimated || [ A; B ] ||_2 ~= " << normScale << endl;
        W *= F(1)/normScale;
        dC *= normScale;
    }

    // Form the augmented RHS
    // ======================
    //   G = [ 0; C; D ]
    Matrix<F> G;
    Zeros( G, n+m+k, numRHS );
    {
        auto Gr = G( IR(n,n+m),     ALL );
        auto Gy = G( IR(n+m,n+m+k), ALL );
        Gr = C;
        Gy = D;
    }

    // Form the augmented matrix
    // =========================
    //
    //         | 0     A^H    B^H |
    //     J = | A  -alpha*I   0  |
    //         | B      0      0  |
    //
    const Int numEntriesW = W.NumEntries();
    SparseMatrix<F> J; 
    Zeros( J, n+m+k, n+m+k );
    J.Reserve( 2*numEntriesW+m ); 
    for( Int e=0; e<numEntriesW; ++e )
    {
        J.QueueUpdate( W.Row(e)+n, W.Col(e),        W.Value(e)  );
        J.QueueUpdate( W.Col(e),   W.Row(e)+n, Conj(W.Value(e)) );
    }
    for( Int e=0; e<m; ++e )
        J.QueueUpdate( e+n, e+n, -ctrl.alpha );
    J.ProcessQueues();
    
    // Add the a priori regularization
    // ===============================
    Matrix<Real> reg;
    Zeros( reg, n+m+k, 1 );
    for( Int i=0; i<n; ++i )
        reg.Set( i, 0, ctrl.qsdCtrl.regPrimal );
    for( Int i=n; i<n+m+k; ++i )
        reg.Set( i, 0, -ctrl.qsdCtrl.regDual );
    SparseMatrix<F> JOrig;
    JOrig = J;
    UpdateRealPartOfDiagonal( J, Real(1), reg );

    // Factor the regularized system
    // =============================
    vector<Int> map, invMap; 
    ldl::NodeInfo info;
    ldl::Separator rootSep;
    ldl::NestedDissection( J.LockedGraph(), map, rootSep, info );
    InvertMap( map, invMap );
    ldl::Front<F> JFront( J, map, info );
    LDL( info, JFront );    

    // Successively solve each of the numRHS linear systems
    // ====================================================
    Matrix<F> u;
    Zeros( u, n+m+k, 1 );
    for( Int j=0; j<numRHS; ++j )
    {
        auto g = G( ALL, IR(j) );
        u = g;
        reg_qsd_ldl::SolveAfter
        ( JOrig, reg, invMap, info, JFront, u, ctrl.qsdCtrl );
        g = u;
    }

    // Extract X from G = [ Dc*X; -R/alpha; Y/alpha ]
    // ==============================================
    X = G( IR(0,n), ALL );
    DiagonalSolve( LEFT, NORMAL, dC, X );
}

template<typename F> 
void LSE
( const DistSparseMatrix<F>& A, const DistSparseMatrix<F>& B, 
  const DistMultiVec<F>& C,     const DistMultiVec<F>& D, 
        DistMultiVec<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LSE"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Height();
    const Int numRHS = C.Width();
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

    // Rescale W = [ A; B ]
    // ====================
    DistMultiVec<Real> dC(comm);
    DistSparseMatrix<F> W(comm);
    VCat( A, B, W );
    if( ctrl.equilibrate )
    {
        ColumnTwoNorms( W, dC );
        auto normMap = []( Real beta )
          { return ( beta < Sqrt(Epsilon<Real>()) ? Real(1) : beta ); };
        EntrywiseMap( dC, function<Real(Real)>(normMap) );
        DiagonalSolve( RIGHT, NORMAL, dC, W );
    }
    else
    {
        Ones( dC, n, 1 );
    }
    Real normScale = 1;
    if( ctrl.scaleTwoNorm )
    {
        normScale = TwoNormEstimate( W, ctrl.basisSize );
        if( ctrl.progress && commRank == 0 )
            cout << "Estimated || [ A; B ] ||_2 ~= " << normScale << endl;
        W *= F(1)/normScale;
        dC *= normScale;
    }

    // Form the augmented RHS
    // ======================
    //   G = [ 0; C; D ]
    DistMultiVec<F> G(comm);
    Zeros( G, n+m+k, numRHS );
    {
        // Compute the metadata
        // --------------------
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<C.LocalHeight(); ++iLoc )
        {
            const Int i = C.GlobalRow(iLoc);
            sendCounts[ G.RowOwner(i+n) ] += numRHS;
        }
        for( Int iLoc=0; iLoc<D.LocalHeight(); ++iLoc )
        {
            const Int i = D.GlobalRow(iLoc);
            sendCounts[ G.RowOwner(i+n+m) ] += numRHS;
        }
        // Pack the data
        // -------------
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        auto offs = sendOffs;
        vector<Entry<F>> sendBuf(totalSend);
        for( Int iLoc=0; iLoc<C.LocalHeight(); ++iLoc )
        {
            const Int i = C.GlobalRow(iLoc);
            const int owner = G.RowOwner(i+n);
            for( Int j=0; j<numRHS; ++j )
                sendBuf[offs[owner]++] = Entry<F>{i+n,j,C.GetLocal(iLoc,j)};
        }
        for( Int iLoc=0; iLoc<D.LocalHeight(); ++iLoc )
        {
            const Int i = D.GlobalRow(iLoc);
            const int owner = G.RowOwner(i+n+m);
            for( Int j=0; j<numRHS; ++j )
                sendBuf[offs[owner]++] = Entry<F>{i+n+m,j,D.GetLocal(iLoc,j)};
        }
        // Exchange and unpack the data
        // ----------------------------
        auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
        for( auto& entry : recvBuf )
            G.Set( entry );
    }

    // Form the augmented matrix
    // =========================
    //
    //         | 0     A^H    B^H |
    //     J = | A  -alpha*I   0  |
    //         | B      0      0  |
    //
    const Int numEntriesW = W.NumLocalEntries();
    DistSparseMatrix<F> J(comm); 
    Zeros( J, n+m+k, n+m+k );
    {
        // Compute the metadata
        // --------------------
        vector<int> sendCounts(commSize,0);
        for( Int e=0; e<numEntriesW; ++e )
        {
            ++sendCounts[ J.RowOwner(W.Row(e)+n) ];
            ++sendCounts[ J.RowOwner(W.Col(e))   ];
        }
        // Pack W
        // ------
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        auto offs = sendOffs;
        vector<Entry<F>> sendBuf(totalSend);
        for( Int e=0; e<numEntriesW; ++e )
        {
            const Int i = W.Row(e);
            const Int j = W.Col(e);
            const F value = W.Value(e);
            // Send this entry of W into its normal position
            int owner = J.RowOwner(i+n);
            sendBuf[offs[owner]++] = Entry<F>{i+n,j,value};
            // Send this entry of W into its adjoint position
            owner = J.RowOwner(j);
            sendBuf[offs[owner]++] = Entry<F>{j,i+n,Conj(value)};
        }
        // Exchange and unpack the data
        // ----------------------------
        auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
        // Count the total number of negative alpha updates
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        Int numNegAlphaUpdates = 0;
        for( Int iLoc=0; iLoc<J.LocalHeight(); ++iLoc )
        {
            const Int i = J.GlobalRow(iLoc);
            if( i >= n+m )
                break;
            else if( i >= n )
                ++numNegAlphaUpdates;
        }
        // Unpack
        // ^^^^^^
        J.Reserve( recvBuf.size() + numNegAlphaUpdates );
        for( auto& entry : recvBuf )
            J.QueueUpdate( entry );
        for( Int iLoc=0; iLoc<J.LocalHeight(); ++iLoc )
        {
            const Int i = J.GlobalRow(iLoc);
            if( i >= n+m )
                break;
            else if( i >= n )
                J.QueueUpdate( i, i, -ctrl.alpha );
        }
        J.ProcessQueues();
    }

    // Add the a priori regularization
    // ===============================
    DistMultiVec<Real> reg(comm);
    Zeros( reg, n+m+k, 1 );
    for( Int iLoc=0; iLoc<reg.LocalHeight(); ++iLoc )
    {
        const Int i = reg.GlobalRow(iLoc);
        if( i < n )
            reg.SetLocal( iLoc, 0, ctrl.qsdCtrl.regPrimal );
        else
            reg.SetLocal( iLoc, 0, -ctrl.qsdCtrl.regDual );
    }
    DistSparseMatrix<F> JOrig(comm);
    JOrig = J;
    UpdateRealPartOfDiagonal( J, Real(1), reg );

    // Factor the regularized system
    // =============================
    DistMap map, invMap; 
    ldl::DistNodeInfo info;
    ldl::DistSeparator rootSep;
    ldl::NestedDissection( J.LockedDistGraph(), map, rootSep, info );
    InvertMap( map, invMap );
    ldl::DistFront<F> JFront( J, map, rootSep, info );
    LDL( info, JFront );    

    // Successively solve each of the numRHS linear systems
    // ====================================================
    DistMultiVec<F> u(comm);
    Zeros( u, n+m+k, 1 );
    auto& GLoc = G.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<numRHS; ++j )
    {
        auto gLoc = GLoc( ALL, IR(j) );
        uLoc = gLoc;
        reg_qsd_ldl::SolveAfter
        ( JOrig, reg, invMap, info, JFront, u, ctrl.qsdCtrl );
        gLoc = uLoc;
    }

    // Extract X from G = [ Dc*X; -R/alpha; Y/alpha ]
    // ==============================================
    X = G( IR(0,n), ALL );
    DiagonalSolve( LEFT, NORMAL, dC, X );
}

#define PROTO(F) \
  template void lse::Overwrite \
  ( Matrix<F>& A, Matrix<F>& B, \
    Matrix<F>& C, Matrix<F>& D, \
    Matrix<F>& X, bool computeResidual ); \
  template void lse::Overwrite \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, \
    AbstractDistMatrix<F>& C, AbstractDistMatrix<F>& D, \
    AbstractDistMatrix<F>& X, bool computeResidual ); \
  template void LSE \
  ( const Matrix<F>& A, const Matrix<F>& B, \
    const Matrix<F>& C, const Matrix<F>& D, \
          Matrix<F>& X ); \
  template void LSE \
  ( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B, \
    const AbstractDistMatrix<F>& C, const AbstractDistMatrix<F>& D, \
          AbstractDistMatrix<F>& X ); \
  template void LSE \
  ( const SparseMatrix<F>& A, const SparseMatrix<F>& B, \
    const Matrix<F>& C,       const Matrix<F>& D, \
          Matrix<F>& X, \
    const LeastSquaresCtrl<Base<F>>& ctrl ); \
  template void LSE \
  ( const DistSparseMatrix<F>& A, const DistSparseMatrix<F>& B, \
    const DistMultiVec<F>& C,     const DistMultiVec<F>& D, \
          DistMultiVec<F>& X, \
    const LeastSquaresCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
