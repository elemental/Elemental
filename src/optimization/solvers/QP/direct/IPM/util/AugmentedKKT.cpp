/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace qp {
namespace direct {

// Form 
//
//    | Q + (x <> z)  A^T | | dx | = | -r_c - x <> r_mu |,
//    |      A         0  | | dy |   | -r_b             |
//
// where 
//
//    r_b  = A x - b,
//    r_c  = Q x + A^T y - z + c,
//    r_mu = x o z - tau e,
//
// and dz can be computed using
//
//   dz = - x <> (r_mu + z o dx)
//

template<typename Real>
void AugmentedKKT
( const Matrix<Real>& Q,
  const Matrix<Real>& A, 
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        Matrix<Real>& J, bool onlyLower )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();

    Zeros( J, m+n, m+n );
    const IR xInd(0,n), yInd(n,n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd);
    Matrix<Real> d( z );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( Jxx, d );
    Jxx += Q;
    Jyx = A;
    if( !onlyLower )
        Transpose( A, Jxy );
}

template<typename Real>
void AugmentedKKT
( const ElementalMatrix<Real>& Q,
  const ElementalMatrix<Real>& A, 
  const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();

    DistMatrixWriteProxy<Real,Real,MC,MR> JProx( JPre );
    auto& J = JProx.Get();

    Zeros( J, m+n, m+n );
    const IR xInd(0,n), yInd(n,n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd);
    DistMatrix<Real,MC,STAR> d( z );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( Jxx, d );
    Jxx += Q;
    Jyx = A;
    if( !onlyLower )
        Transpose( A, Jxy );
}

template<typename Real>
void AugmentedKKT
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A, 
        Real gamma,
        Real delta,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntriesQ = Q.NumEntries();
    const Int numEntriesA = A.NumEntries();

    // Count the number of used entries of Q
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

    Zeros( J, m+n, m+n );
    if( onlyLower )
        J.Reserve( numEntriesA + numUsedEntriesQ + n+m );
    else
        J.Reserve( 2*numEntriesA + numUsedEntriesQ + n+m ); 

    // x o inv(z) + gamma^2*I updates
    for( Int j=0; j<n; ++j )
        J.QueueUpdate( j, j, z(j)/x(j)+gamma*gamma );

    // Q update
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j )
            J.QueueUpdate( i, j, Q.Value(e) );
        else if( !onlyLower )
            J.QueueUpdate( i, j, Q.Value(e) );
    }

    // A and A^T updates
    for( Int e=0; e<numEntriesA; ++e )
    {
        J.QueueUpdate( A.Row(e)+n, A.Col(e), A.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( A.Col(e), A.Row(e)+n, A.Value(e) );
    }

    // -delta^2*I 
    for( Int i=0; i<m; ++i )
        J.QueueUpdate( i+n, i+n, -delta*delta );

    J.ProcessQueues();
    J.FreezeSparsity();
}

template<typename Real>
void AugmentedKKT
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
        Real gamma,
        Real delta,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntriesQ = Q.NumLocalEntries();
    const Int numEntriesA = A.NumLocalEntries();
    auto& xLoc = x.LockedMatrix();
    auto& zLoc = z.LockedMatrix();

    J.SetComm( A.Comm() );
    Zeros( J, m+n, m+n );
    const Int JLocalHeight = J.LocalHeight();

    // Compute the number of entries to send
    // =====================================
    Int numEntries = 0;
    numEntries += numEntriesA;
    if( !onlyLower ) 
        numEntries += numEntriesA;
    numEntries += x.LocalHeight(); 
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            ++numEntries;
    }
    for( Int iLoc=0; iLoc<JLocalHeight; ++iLoc )
    {
        const Int i = J.GlobalRow(iLoc);
        if( i > n )
            ++numEntries;
    }

    // Queue the entries
    // =================
    J.Reserve( numEntries, numEntries );
    // Pack A
    // ------
    for( Int e=0; e<numEntriesA; ++e )
    {
        const Int i = A.Row(e) + n;
        const Int j = A.Col(e);
        J.QueueUpdate( i, j, A.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( j, i, A.Value(e) );
    }
    // Pack x o inv(z) + gamma^2*I
    // ---------------------------
    for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Real value = zLoc(iLoc)/xLoc(iLoc)+gamma*gamma;
        J.QueueUpdate( i, i, value );
    }
    // Pack Q
    // ------
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            J.QueueUpdate( i, j, Q.Value(e) );
    }
    // Pack -delta^2*I
    // ---------------
    for( Int iLoc=0; iLoc<JLocalHeight; ++iLoc )
    {
        const Int i = J.GlobalRow(iLoc);
        if( i > n )
            J.QueueUpdate( i, i, -delta*delta );
    }

    J.ProcessQueues();
    J.FreezeSparsity();
}

template<typename Real>
void AugmentedKKTRHS
( const Matrix<Real>& x, 
  const Matrix<Real>& rc,
  const Matrix<Real>& rb, 
  const Matrix<Real>& rmu,
        Matrix<Real>& d )
{
    DEBUG_CSE
    const Int m = rb.Height();
    const Int n = rmu.Height();
    const IR xInd(0,n), yInd(n,n+m);
    Zeros( d, m+n, 1 );

    // dx := - (r_c + x <> r_mu)
    // =========================
    auto dx = d(xInd,ALL);
    dx = rmu;
    DiagonalSolve( LEFT, NORMAL, x, dx );
    dx += rc;
    dx *= -1;

    // dy := -r_b
    // ==========
    auto dy = d(yInd,ALL);
    dy = rb;
    dy *= -1;
}

template<typename Real>
void AugmentedKKTRHS
( const ElementalMatrix<Real>& x, 
  const ElementalMatrix<Real>& rc,
  const ElementalMatrix<Real>& rb, 
  const ElementalMatrix<Real>& rmu,
        ElementalMatrix<Real>& dPre )
{
    DEBUG_CSE

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.rowConstrain = true;
    ctrl.colAlign = 0;
    ctrl.rowAlign = 0;

    DistMatrixWriteProxy<Real,Real,MC,MR> dProx( dPre, ctrl );
    auto& d = dProx.Get();

    const Int m = rb.Height();
    const Int n = rmu.Height();
    const IR xInd(0,n), yInd(n,n+m);
    Zeros( d, m+n, 1 );

    // dx := - (r_c + x <> r_mu)
    // =========================
    auto dx = d(xInd,ALL);
    dx = rmu;
    DiagonalSolve( LEFT, NORMAL, x, dx );
    dx += rc;
    dx *= -1;

    // dy := -r_b
    // ==========
    auto dy = d(yInd,ALL);
    dy = rb;
    dy *= -1;
}

template<typename Real>
void AugmentedKKTRHS
( const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rb, 
  const DistMultiVec<Real>& rmu, 
        DistMultiVec<Real>& d )
{
    DEBUG_CSE
    const Int m = rb.Height();
    const Int n = x.Height();
    auto& xLoc = x.LockedMatrix();
    auto& rcLoc = rc.LockedMatrix();
    auto& rbLoc = rb.LockedMatrix();
    auto& rmuLoc = rmu.LockedMatrix();

    d.SetComm( x.Comm() );
    Zeros( d, m+n, 1 );

    d.Reserve( rc.LocalHeight()+rb.LocalHeight() );
    for( Int iLoc=0; iLoc<rc.LocalHeight(); ++iLoc )
    {
        const Int i = rc.GlobalRow(iLoc);
        const Real value = -rcLoc(iLoc) - rmuLoc(iLoc)/xLoc(iLoc);
        d.QueueUpdate( i, 0, value );
    }
    for( Int iLoc=0; iLoc<rb.LocalHeight(); ++iLoc )
    {
        const Int i = rb.GlobalRow(iLoc) + n;
        const Real value = -rbLoc(iLoc);
        d.QueueUpdate( i, 0, value );
    }
    d.ProcessQueues();
}

template<typename Real>
void ExpandAugmentedSolution
( const Matrix<Real>& x,
  const Matrix<Real>& z,
  const Matrix<Real>& rmu,
  const Matrix<Real>& d,
        Matrix<Real>& dx,
        Matrix<Real>& dy, 
        Matrix<Real>& dz )
{
    DEBUG_CSE
    const Int n = rmu.Height();
    const Int m = d.Height() - n;

    // Extract dx and dy from [dx; dy]
    // ===============================
    dx = d(IR(0,n  ),ALL);
    dy = d(IR(n,n+m),ALL);

    // dz := - x <> (r_mu + z o dx)
    // ============================
    dz = dx;
    DiagonalScale( LEFT, NORMAL, z, dz );
    dz += rmu;
    DiagonalSolve( LEFT, NORMAL, x, dz );
    dz *= -1;
}

template<typename Real>
void ExpandAugmentedSolution
( const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& z,
  const ElementalMatrix<Real>& rmu,
  const ElementalMatrix<Real>& dPre,
        ElementalMatrix<Real>& dx,
        ElementalMatrix<Real>& dy, 
        ElementalMatrix<Real>& dz )
{
    DEBUG_CSE

    DistMatrixReadProxy<Real,Real,MC,MR> dProx( dPre );
    auto& d = dProx.GetLocked();

    const Int n = rmu.Height();
    const Int m = d.Height() - n;

    // Extract dx and dy from [dx; dy]
    // ===============================
    dx = d(IR(0,n  ),ALL);
    dy = d(IR(n,n+m),ALL);

    // dz := - x <> (r_mu + z o dx)
    // ============================
    dz = dx;
    DiagonalScale( LEFT, NORMAL, z, dz );
    dz += rmu;
    DiagonalSolve( LEFT, NORMAL, x, dz );
    dz *= -1;
}

template<typename Real>
void ExpandAugmentedSolution
( const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& d,
        DistMultiVec<Real>& dx,
        DistMultiVec<Real>& dy, 
        DistMultiVec<Real>& dz )
{
    DEBUG_CSE
    const Int n = rmu.Height();
    const Int m = d.Height() - n;
    mpi::Comm comm = z.Comm();
    dx.SetComm( comm );
    dy.SetComm( comm );
    dz.SetComm( comm );

    auto& dLoc = d.LockedMatrix();

    // Extract dx and dy from [dx; dy]
    // ===============================
    Zeros( dx, n, 1 );
    Zeros( dy, m, 1 );
    Int dxNumEntries=0, dyNumEntries=0;
    for( Int iLoc=0; iLoc<d.LocalHeight(); ++iLoc )
    {
        const Int i = d.GlobalRow(iLoc);
        if( i < n )
            ++dxNumEntries;
        else
            ++dyNumEntries;
    }
    dx.Reserve( dxNumEntries );
    dy.Reserve( dyNumEntries );
    for( Int iLoc=0; iLoc<d.LocalHeight(); ++iLoc )
    {
        const Int i = d.GlobalRow(iLoc);
        const Real value = dLoc(iLoc);
        if( i < n )
            dx.QueueUpdate( i, 0, value );
        else
            dy.QueueUpdate( i-n, 0, value );
    }
    dx.ProcessQueues();
    dy.ProcessQueues();

    // dz := - x <> (r_mu + z o dx)
    // ============================
    dz = dx;
    DiagonalScale( LEFT, NORMAL, z, dz );
    dz += rmu;
    DiagonalSolve( LEFT, NORMAL, x, dz );
    dz *= -1;
}

#define PROTO(Real) \
  template void AugmentedKKT \
  ( const Matrix<Real>& Q, \
    const Matrix<Real>& A, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
          Matrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const ElementalMatrix<Real>& Q, \
    const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Real>& z, \
          ElementalMatrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const SparseMatrix<Real>& Q, \
    const SparseMatrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
          SparseMatrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const DistSparseMatrix<Real>& Q, \
    const DistSparseMatrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void AugmentedKKTRHS \
  ( const Matrix<Real>& x, \
    const Matrix<Real>& rc, \
    const Matrix<Real>& rb, \
    const Matrix<Real>& rmu, \
          Matrix<Real>& d ); \
  template void AugmentedKKTRHS \
  ( const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Real>& rc, \
    const ElementalMatrix<Real>& rb, \
    const ElementalMatrix<Real>& rmu, \
          ElementalMatrix<Real>& d ); \
  template void AugmentedKKTRHS \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& rc, \
    const DistMultiVec<Real>& rb, \
    const DistMultiVec<Real>& rmu, \
          DistMultiVec<Real>& d ); \
  template void ExpandAugmentedSolution \
  ( const Matrix<Real>& x, \
    const Matrix<Real>& z, \
    const Matrix<Real>& rmu, \
    const Matrix<Real>& d, \
          Matrix<Real>& dx, \
          Matrix<Real>& dy, \
          Matrix<Real>& dz ); \
  template void ExpandAugmentedSolution \
  ( const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Real>& z, \
    const ElementalMatrix<Real>& rmu, \
    const ElementalMatrix<Real>& d, \
          ElementalMatrix<Real>& dx, \
          ElementalMatrix<Real>& dy, \
          ElementalMatrix<Real>& dz ); \
  template void ExpandAugmentedSolution \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& rmu, \
    const DistMultiVec<Real>& d, \
          DistMultiVec<Real>& dx, \
          DistMultiVec<Real>& dy, \
          DistMultiVec<Real>& dz );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace direct
} // namespace qp
} // namespace El
