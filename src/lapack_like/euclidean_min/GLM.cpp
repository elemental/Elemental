/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// This file implements both dense and sparse-direct solutions of 
// General (Gauss-Markov) Linear Model (GLM):
//
//     min_{x,y} || y ||_2 such that A x + B y = d. 
//
// For dense instances of the problem, where A is m x n and B is m x p, 
// we assume that n <= m <= n+p as well as that A has full column rank, n, and 
// [A B] has full row rank, m.
//
// A Generalized QR factorization of (A,B),
//     A = Q R = Q | R11 |, B = Q T Z = Q | T11 T12 | Z,
//                 | 0   |                |   0 T22 |
// where Q and Z are unitary and R and T are upper-trapezoidal, allows us to
// re-express the constraint as 
//     (Q^H d) = | R11 | x + | T11 T12 | (Z y).
//               |   0 |     |   0 T22 |
// which is re-written as
//      | g1 | = | R11 x + T11 c1 + T12 c2 |
//      | g2 |   |                  T22 c2 |.
// Since || c ||_2 == || Z y ||_2 = || y ||_2 is to be minimized, and c2 is 
// fixed, our only freedom is in the choice of c1, which we set to zero.
// Then all that is left is to solve
//      R11 x = g1 - T12 c2
// for x.
//
// Note that essentially the same scheme is used in LAPACK's {S,D,C,Z}GGGLM.
//
// For sparse instances of the GLM problem, the symmetric quasi-semidefinite
// augmented system
//
//     |  0   A      B    | |    z    |   | d/alpha |
//     | A^H  0      0    | | x/alpha | = |    0    |
//     | B^H  0  -alpha*I | | y/alpha |   |    0    |
//
// is formed, equilibrated, and then a priori regularization is added in order
// to make the system sufficiently quasi-definite. A Cholesky-like factorization
// of this regularized system is then used as a preconditioner for FGMRES(k).
//

namespace El {

namespace glm {

// For the following two routines, on exit, A and B are overwritten with their 
// implicit Generalized QR factorization and D is overwritten with X
template<typename F> 
void Overwrite( Matrix<F>& A, Matrix<F>& B, Matrix<F>& D, Matrix<F>& Y )
{
    DEBUG_ONLY(CSE cse("glm::Overwrite"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = B.Width();
    const Int numRhs = D.Width();
    if( m != B.Height() || m != D.Height() )
        LogicError("A, B, and D must be the same height");
    if( m < n )
        LogicError("GLM requires height(A) >= width(A)");
    if( n+p < m )
        LogicError("GLM requires width(A)+width(B) >= height(A)");
    const bool checkIfSingular = true;

    // Compute the implicit Generalized QR decomposition of (A,B)
    Matrix<F> tA, tB;
    Matrix<Base<F>> dA, dB;
    GQR( A, tA, dA, B, tB, dB );

    // G := Q^H D
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, D );

    // Partition the relevant matrices
    Matrix<F> G1, G2;
    PartitionDown( D, G1, G2, n );
    Matrix<F> R11, R21;
    PartitionDown( A, R11, R21, n );
    Matrix<F> T11, T12, T21, T22;
    PartitionUpOffsetDiagonal
    ( p-m,
      B, T11, T12,
         T21, T22, m-n );
    Zeros( Y, p, numRhs );
    Matrix<F> C1, C2;
    PartitionDown( Y, C1, C2, n+p-m );

    // Solve T22 C2 = G2
    C2 = G2;
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), T22, C2, checkIfSingular );

    // G1 := G1 - T12 C2
    Gemm( NORMAL, NORMAL, F(-1), T12, C2, F(1), G1 );
    
    // Solve R11 X = G1 
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R11, G1, checkIfSingular );
    D.Resize( n, numRhs );

    // Y := Z^H C
    rq::ApplyQ( LEFT, ADJOINT, B, tB, dB, Y );
}

template<typename F> 
void Overwrite
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& BPre, 
  AbstractDistMatrix<F>& DPre, AbstractDistMatrix<F>& YPre )
{
    DEBUG_ONLY(CSE cse("glm::Overwrite"))

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre ); auto& A = *APtr;
    auto BPtr = ReadWriteProxy<F,MC,MR>( &BPre ); auto& B = *BPtr;
    auto DPtr = ReadWriteProxy<F,MC,MR>( &DPre ); auto& D = *DPtr;
    auto YPtr = WriteProxy<F,MC,MR>( &YPre );     auto& Y = *YPtr;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = B.Width();
    const Int numRhs = D.Width();
    if( m != B.Height() || m != D.Height() )
        LogicError("A, B, and D must be the same height");
    if( m < n )
        LogicError("GLM requires height(A) >= width(A)");
    if( n+p < m )
        LogicError("GLM requires width(A)+width(B) >= height(A)");
    const Grid& g = A.Grid();
    if( g != B.Grid() || g != D.Grid() )
        LogicError("All matrices must have the same grid");
    Y.SetGrid( g );
    const bool checkIfSingular = true;

    // Compute the implicit Generalized QR decomposition of (A,B)
    DistMatrix<F,MD,STAR> tA(g), tB(g);
    DistMatrix<Base<F>,MD,STAR> dA(g), dB(g);
    GQR( A, tA, dA, B, tB, dB );

    // G := Q^H D
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, D );

    // Partition the relevant matrices
    DistMatrix<F> G1(g), G2(g);
    PartitionDown( D, G1, G2, n );
    DistMatrix<F> R11(g), R21(g);
    PartitionDown( A, R11, R21, n );
    DistMatrix<F> T11(g), T12(g), T21(g), T22(g);
    PartitionUpOffsetDiagonal
    ( p-m,
      B, T11, T12,
         T21, T22, m-n );
    Zeros( Y, p, numRhs );
    DistMatrix<F> C1(g), C2(g);
    PartitionDown( Y, C1, C2, n+p-m );

    // Solve T22 C2 = G2
    C2 = G2;
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), T22, C2, checkIfSingular );

    // G1 := G1 - T12 C2
    Gemm( NORMAL, NORMAL, F(-1), T12, C2, F(1), G1 );
    
    // Solve R11 X = G1
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R11, G1, checkIfSingular );
    D.Resize( n, numRhs );

    // Y := Z^H C
    rq::ApplyQ( LEFT, ADJOINT, B, tB, dB, Y );
}

} // namespace glm

template<typename F> 
void GLM
( const Matrix<F>& A, const Matrix<F>& B, 
  const Matrix<F>& D, 
        Matrix<F>& X, Matrix<F>& Y )
{
    DEBUG_ONLY(CSE cse("GLM"))
    Matrix<F> ACopy( A ), BCopy( B );
    X = D;
    glm::Overwrite( ACopy, BCopy, X, Y );
}

template<typename F> 
void GLM
( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B, 
  const AbstractDistMatrix<F>& D, 
        AbstractDistMatrix<F>& X,       AbstractDistMatrix<F>& Y )
{
    DEBUG_ONLY(CSE cse("GLM"))
    DistMatrix<F> ACopy( A ), BCopy( B );
    Copy( D, X );
    glm::Overwrite( ACopy, BCopy, X, Y );
}

template<typename F>
void GLM
( const SparseMatrix<F>& A, const SparseMatrix<F>& B,
  const Matrix<F>& D,
        Matrix<F>& X,             Matrix<F>& Y, 
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("GLM"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    const Int numRHS = D.Width();

    // Rescale the rows of W := [ A, B ]
    // ===============================
    Matrix<Real> dR;
    SparseMatrix<F> W;
    HCat( A, B, W );
    if( ctrl.equilibrate )
    {
        RowTwoNorms( W, dR ); 
        auto normMap = []( Real beta )
          { return ( beta < Sqrt(Epsilon<Real>()) ? Real(1) : beta ); };
        EntrywiseMap( dR, function<Real(Real)>(normMap) );
        DiagonalSolve( LEFT, NORMAL, dR, W );
    }
    else
    {
        Ones( dR, m, 1 );
    }
    Real normScale = 1;
    if( ctrl.scaleTwoNorm )
    {
        normScale = TwoNormEstimate( W, ctrl.basisSize );
        if( ctrl.progress )
            cout << "Estimated || [ A, B ] ||_2 ~= " << normScale << endl;
        W *= F(1)/normScale;
        dR *= normScale;
    }

    // Form the augmented RHS (and rescale it)
    // =======================================
    //   G = [ D/alpha; 0; 0 ]
    Matrix<F> G;
    Zeros( G, m+n+k, numRHS );
    {
        auto Gz = G( IR(0,m), ALL );
        // Use X as a temporary
        X = D;
        X *= F(1)/ctrl.alpha;
        DiagonalSolve( LEFT, NORMAL, dR, X );
        Gz = X;
        X.Empty();
    }

    // Form the augmented matrix
    // =========================
    //
    //         | 0    A      B    |
    //     J = | A^H  0      0    |
    //         | B^H  0  -alpha*I |
    //
    const Int numEntriesW = W.NumEntries();
    SparseMatrix<F> J;
    Zeros( J, m+n+k, m+n+k );
    J.Reserve( 2*numEntriesW+m );
    for( Int e=0; e<numEntriesW; ++e )
    {
        J.QueueUpdate( W.Row(e),   W.Col(e)+m, W.Value(e)       );
        J.QueueUpdate( W.Col(e)+m, W.Row(e),   Conj(W.Value(e)) );
    }
    for( Int e=0; e<k; ++e )
        J.QueueUpdate( e+m+n, e+m+n, -ctrl.alpha );
    J.ProcessQueues();

    // Add the a priori regularization
    // ===============================
    Matrix<Real> reg;
    Zeros( reg, m+n+k, 1 );
    for( Int i=0; i<m; ++i )
        reg.Set( i, 0, ctrl.qsdCtrl.regPrimal );
    for( Int i=m; i<m+n+k; ++i )
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
    Zeros( u, m+n+k, 1 );
    for( Int j=0; j<numRHS; ++j )
    {
        auto g = G( ALL, IR(j) );
        u = g;
        reg_qsd_ldl::SolveAfter
        ( JOrig, reg, invMap, info, JFront, u, ctrl.qsdCtrl );
        g = u;
    }

    // Extract X and Y from G = [ Z; X/alpha; Y/alpha ]
    // ================================================
    // Well, actually, the solution has been equilibrated, but the division
    // by alpha commutes.
    X = G( IR(m,m+n),     ALL );
    Y = G( IR(m+n,m+n+k), ALL );
    X *= ctrl.alpha;
    Y *= ctrl.alpha;
}

template<typename F>
void GLM
( const DistSparseMatrix<F>& A, const DistSparseMatrix<F>& B,
  const DistMultiVec<F>& D,
        DistMultiVec<F>& X,           DistMultiVec<F>& Y,
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("GLM"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    const Int numRHS = D.Width();
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

    // Rescale the rows of W := [ A, B ]
    // ===============================
    DistMultiVec<Real> dR(comm);
    DistSparseMatrix<F> W(comm);
    HCat( A, B, W );
    if( ctrl.equilibrate )
    {
        RowTwoNorms( W, dR );
        auto normMap = []( Real beta )
          { return ( beta < Sqrt(Epsilon<Real>()) ? Real(1) : beta ); };
        EntrywiseMap( dR, function<Real(Real)>(normMap) );
        DiagonalSolve( LEFT, NORMAL, dR, W );
    }
    else
    {
        Ones( dR, m, 1 );
    }
    Real normScale = 1;
    if( ctrl.scaleTwoNorm )
    {
        normScale = TwoNormEstimate( W, ctrl.basisSize );
        if( ctrl.progress && commRank == 0 )
            cout << "Estimated || [ A, B ] ||_2 ~= " << normScale << endl;
        W *= F(1)/normScale;
        dR *= normScale;
    }

    // Form the augmented RHS
    // ======================
    //   G = [ D/alpha; 0; 0 ]
    DistMultiVec<F> G(comm);
    Zeros( G, m+n+k, numRHS );
    X = D;
    X *= F(1)/ctrl.alpha;
    DiagonalSolve( LEFT, NORMAL, dR, X );
    {
        // Compute the metadata
        // --------------------
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<X.LocalHeight(); ++iLoc )
            sendCounts[ G.RowOwner(X.GlobalRow(iLoc)) ] += numRHS;
        // Pack the data
        // -------------
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        auto offs = sendOffs;
        vector<Entry<F>> sendBuf(totalSend);
        for( Int iLoc=0; iLoc<X.LocalHeight(); ++iLoc )
        {
            const Int i = X.GlobalRow(iLoc);
            const int owner = G.RowOwner(i);
            for( Int j=0; j<numRHS; ++j )
                sendBuf[offs[owner]++] = Entry<F>{ i, j, X.GetLocal(iLoc,j) };
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
    //         | 0    A      B    |
    //     J = | A^H  0      0    |
    //         | B^H  0  -alpha*I |
    //
    const Int numEntriesW = W.NumLocalEntries();
    DistSparseMatrix<F> J(comm);
    Zeros( J, m+n+k, m+n+k );
    {
        // Compute the metadata
        // --------------------
        vector<int> sendCounts(commSize,0);
        for( Int e=0; e<numEntriesW; ++e )
        {
            ++sendCounts[ J.RowOwner(W.Row(e))   ];
            ++sendCounts[ J.RowOwner(W.Col(e)+m) ];
        }
        // Pack the A and B data
        // ---------------------
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
            int owner = J.RowOwner(i);
            sendBuf[offs[owner]++] = Entry<F>{i,j+m,value};
            // Send this entry of W into its adjoint position
            owner = J.RowOwner(j+m);
            sendBuf[offs[owner]++] = Entry<F>{j+m,i,Conj(value)};
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
            if( i >= m+n )
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
            if( i >= m+n )
                J.QueueUpdate( i, i, -ctrl.alpha );
        }
        J.ProcessQueues();
    }

    // Add the a priori regularization
    // ===============================
    DistMultiVec<Real> reg(comm);
    Zeros( reg, m+n+k, 1 );
    for( Int iLoc=0; iLoc<reg.LocalHeight(); ++iLoc )
    {
        const Int i = reg.GlobalRow(iLoc);
        if( i < m )
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
    Zeros( u, m+n+k, 1 );
    auto& GLoc = G.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<numRHS; ++j )
    {
        auto gLoc = GLoc( ALL, IR(j) );
        Copy( gLoc, uLoc );
        reg_qsd_ldl::SolveAfter
        ( JOrig, reg, invMap, info, JFront, u, ctrl.qsdCtrl );
        Copy( uLoc, gLoc );
    }

    // Extract X and Y from G = [ Z; X/alpha; Y/alpha ]
    // ================================================
    // Well, actually, the solution has been equilibrated, but the division
    // by alpha commutes.
    X = G( IR(m,m+n),     ALL );
    Y = G( IR(m+n,m+n+k), ALL );
    X *= ctrl.alpha;
    Y *= ctrl.alpha;
}

#define PROTO(F) \
  template void glm::Overwrite \
  ( Matrix<F>& A, Matrix<F>& B, Matrix<F>& D, Matrix<F>& Y ); \
  template void glm::Overwrite \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, \
    AbstractDistMatrix<F>& D, AbstractDistMatrix<F>& Y ); \
  template void GLM \
  ( const Matrix<F>& A, const Matrix<F>& B, \
    const Matrix<F>& D, \
          Matrix<F>& X,       Matrix<F>& Y ); \
  template void GLM \
  ( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B, \
    const AbstractDistMatrix<F>& D, \
          AbstractDistMatrix<F>& X, AbstractDistMatrix<F>& Y ); \
  template void GLM \
  ( const SparseMatrix<F>& A, const SparseMatrix<F>& B, \
    const Matrix<F>& D, \
          Matrix<F>& X,             Matrix<F>& Y, \
    const LeastSquaresCtrl<Base<F>>& ctrl ); \
  template void GLM \
  ( const DistSparseMatrix<F>& A, const DistSparseMatrix<F>& B, \
    const DistMultiVec<F>& D, \
          DistMultiVec<F>& X,           DistMultiVec<F>& Y, \
    const LeastSquaresCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
