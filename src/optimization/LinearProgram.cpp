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
    const Int blocksizeJ = J.Blocksize();
    // For placing A into the top-right corner
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for( Int k=0; k<A.NumLocalEntries(); ++k )
    {
        const Int i = A.Row(k);
        const Int owner = RowToProcess( i, blocksizeJ, commSize );
        ++sendCounts[owner];
    }
    // For placing A^T into the bottom-left corner
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for( Int k=0; k<ATrans.NumLocalEntries(); ++k )
    {
        const Int i = ATrans.Row(k) + m;
        const Int owner = RowToProcess( i, blocksizeJ, commSize );
        ++sendCounts[owner];
    }
    // For placing -diag(s)/diag(x) into the bottom-right corner
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for( Int k=0; k<x.LocalHeight(); ++k )
    {
        const Int i = k + x.FirstLocalRow() + m;
        const Int owner = RowToProcess( i, blocksizeJ, commSize );
        ++sendCounts[owner];
    }

    // Communicate to determine the number we receive from each process
    // ----------------------------------------------------------------
    std::vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    
    // Convert the send/recv counts into offsets and total sizes
    // ---------------------------------------------------------
    Int totalSend=0, totalRecv=0;
    std::vector<int> sendOffsets(commSize), recvOffsets(commSize);
    for( Int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = totalSend;
        recvOffsets[q] = totalRecv;
        totalSend += sendCounts[q];
        totalRecv += recvCounts[q];
    }

    // Pack the triplets
    // -----------------
    std::vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
    std::vector<Real> vSendBuf(totalSend);
    std::vector<int> offsets = sendOffsets;
    // Pack A
    // ^^^^^^
    for( Int k=0; k<A.NumLocalEntries(); ++k ) 
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k) + m;
        const Real value = A.Value(k);
        const Int owner = RowToProcess( i, blocksizeJ, commSize );
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
        const Real value = A.Value(k);
        const Int owner = RowToProcess( i, blocksizeJ, commSize );
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
        const Int owner = RowToProcess( i, blocksizeJ, commSize );
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
    {
        const Int i = k + yT.FirstLocalRow();
        const Int owner = RowToProcess( i, blocksizeJ, commSize );
        ++sendCounts[owner];
    }
    for( Int k=0; k<yB.LocalHeight(); ++k )
    {
        const Int i = k + yB.FirstLocalRow() + m;
        const Int owner = RowToProcess( i, blocksizeJ, commSize );
        ++sendCounts[owner];
    }

    // Communicate to determine the number we receive from each process
    // ----------------------------------------------------------------
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    
    // Convert the send/recv counts into offsets and total sizes
    // ---------------------------------------------------------
    totalSend=0, totalRecv=0;
    for( Int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = totalSend;
        recvOffsets[q] = totalRecv;
        totalSend += sendCounts[q];
        totalRecv += recvCounts[q];
    }

    // Pack the triplets
    // -----------------
    sSendBuf.resize(totalSend);
    vSendBuf.resize(totalSend);
    offsets = sendOffsets;
    for( Int k=0; k<yT.LocalHeight(); ++k )
    {
        const Int i = k + yT.FirstLocalRow();
        const Real value = yT.GetLocal(k,0);
        const Int owner = RowToProcess( i, blocksizeJ, commSize );
        sSendBuf[offsets[owner]] = i;
        vSendBuf[offsets[owner]] = value;
        ++offsets[owner];
    }
    for( Int k=0; k<yB.LocalHeight(); ++k )
    {
        const Int i = k + yB.FirstLocalRow() + m;
        const Real value = yB.GetLocal(k,0);
        const Int owner = RowToProcess( i, blocksizeJ, commSize );
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
