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
( const Matrix<Real>& Q, const Matrix<Real>& A, 
  const Matrix<Real>& x, const Matrix<Real>& z,
        Matrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("qp::direct::AugmentedKKT"))
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
( const AbstractDistMatrix<Real>& Q,    const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& x,    const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("qp::direct::AugmentedKKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    auto JPtr = WriteProxy<Real,MC,MR>(&JPre);
    auto& J = *JPtr;

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
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, 
  const Matrix<Real>& x,       const Matrix<Real>& z,
        SparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("qp::direct::AugmentedKKT"))
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
        J.Reserve( numEntriesA + numUsedEntriesQ + n );
    else
        J.Reserve( 2*numEntriesA + numUsedEntriesQ + n ); 

    // x <> z updates
    for( Int j=0; j<n; ++j )
        J.QueueUpdate( j, j, z.Get(j,0)/x.Get(j,0) );
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
    J.ProcessQueues();
}

template<typename Real>
void AugmentedKKT
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("qp::direct::AugmentedKKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    J.SetComm( A.Comm() );
    Zeros( J, m+n, m+n );

    // Compute the number of entries to send
    // =====================================
    Int numEntries = 0;
    numEntries += A.NumLocalEntries();
    if( !onlyLower ) 
        numEntries += A.NumLocalEntries();
    numEntries += x.LocalHeight(); 
    for( Int e=0; e<Q.NumLocalEntries(); ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            ++numEntries;
    }

    // Queue the entries
    // =================
    J.Reserve( numEntries, numEntries );
    // Pack A
    // ------
    for( Int e=0; e<A.NumLocalEntries(); ++e )
    {
        const Int i = A.Row(e) + n;
        const Int j = A.Col(e);
        J.QueueUpdate( i, j, A.Value(e), false );
        if( !onlyLower )
            J.QueueUpdate( j, i, A.Value(e), false );
    }
    // Pack x <> z
    // -----------
    for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Real value = z.GetLocal(iLoc,0)/x.GetLocal(iLoc,0);
        J.QueueUpdate( i, i, value, false );
    }
    // Pack Q
    // ------
    for( Int e=0; e<Q.NumLocalEntries(); ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            J.QueueUpdate( i, j, Q.Value(e), false );
    }
    J.ProcessQueues();
}

template<typename Real>
void AugmentedKKTRHS
( const Matrix<Real>& x, 
  const Matrix<Real>& rc,  const Matrix<Real>& rb, 
  const Matrix<Real>& rmu,
        Matrix<Real>& d )
{
    DEBUG_ONLY(CSE cse("qp::direct::AugmentedKKTRHS"))
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
( const AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Real>& rc,   const AbstractDistMatrix<Real>& rb, 
  const AbstractDistMatrix<Real>& rmu,
        AbstractDistMatrix<Real>& dPre )
{
    DEBUG_ONLY(CSE cse("qp::direct::AugmentedKKTRHS"))

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.rowConstrain = true;
    ctrl.colAlign = 0;
    ctrl.rowAlign = 0;
    auto dPtr = WriteProxy<Real,MC,MR>(&dPre,ctrl); 
    auto& d = *dPtr;

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
  const DistMultiVec<Real>& rc,  const DistMultiVec<Real>& rb, 
  const DistMultiVec<Real>& rmu, 
        DistMultiVec<Real>& d )
{
    DEBUG_ONLY(CSE cse("qp::direct::FormAugmentedSystem"))
    const Int m = rb.Height();
    const Int n = x.Height();
    d.SetComm( x.Comm() );
    Zeros( d, m+n, 1 );

    d.Reserve( rc.LocalHeight()+rb.LocalHeight() );
    for( Int iLoc=0; iLoc<rc.LocalHeight(); ++iLoc )
    {
        const Int i = rc.GlobalRow(iLoc);
        const Real value = -rc.GetLocal(iLoc,0) -
                            rmu.GetLocal(iLoc,0)/x.GetLocal(iLoc,0);
        d.QueueUpdate( i, 0, value );
    }
    for( Int iLoc=0; iLoc<rb.LocalHeight(); ++iLoc )
    {
        const Int i = rb.GlobalRow(iLoc) + n;
        const Real value = -rb.GetLocal(iLoc,0);
        d.QueueUpdate( i, 0, value );
    }
    d.ProcessQueues();
}

template<typename Real>
void ExpandAugmentedSolution
( const Matrix<Real>& x,   const Matrix<Real>& z,
  const Matrix<Real>& rmu, const Matrix<Real>& d,
        Matrix<Real>& dx,        Matrix<Real>& dy, 
        Matrix<Real>& dz )
{
    DEBUG_ONLY(CSE cse("qp::direct::ExpandAugmentedSolution"))
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
( const AbstractDistMatrix<Real>& x,   const AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& dPre,
        AbstractDistMatrix<Real>& dx,        AbstractDistMatrix<Real>& dy, 
        AbstractDistMatrix<Real>& dz )
{
    DEBUG_ONLY(CSE cse("qp::direct::ExpandAugmentedSolution"))

    auto dPtr = ReadProxy<Real,MC,MR>(&dPre);
    auto& d = *dPtr;

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
( const DistMultiVec<Real>& x,   const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& d,
        DistMultiVec<Real>& dx,        DistMultiVec<Real>& dy, 
        DistMultiVec<Real>& dz )
{
    DEBUG_ONLY(CSE cse("qp::direct::ExpandAugmentedSolution"))
    const Int n = rmu.Height();
    const Int m = d.Height() - n;
    mpi::Comm comm = z.Comm();
    dx.SetComm( comm );
    dy.SetComm( comm );
    dz.SetComm( comm );

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
        const Real value = d.GetLocal(iLoc,0);
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
  ( const Matrix<Real>& Q, const Matrix<Real>& A, \
    const Matrix<Real>& x, const Matrix<Real>& z, \
    Matrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z, \
    AbstractDistMatrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, \
    const Matrix<Real>& x,       const Matrix<Real>& z, \
    SparseMatrix<Real>& J,       bool onlyLower ); \
  template void AugmentedKKT \
  ( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z, \
    DistSparseMatrix<Real>& J,       bool onlyLower ); \
  template void AugmentedKKTRHS \
  ( const Matrix<Real>& x, \
    const Matrix<Real>& rc, const Matrix<Real>& rb, \
    const Matrix<Real>& rmu, \
          Matrix<Real>& d ); \
  template void AugmentedKKTRHS \
  ( const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Real>& rc,  const AbstractDistMatrix<Real>& rb, \
    const AbstractDistMatrix<Real>& rmu, \
          AbstractDistMatrix<Real>& d ); \
  template void AugmentedKKTRHS \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& rc,  const DistMultiVec<Real>& rb, \
    const DistMultiVec<Real>& rmu, \
          DistMultiVec<Real>& d ); \
  template void ExpandAugmentedSolution \
  ( const Matrix<Real>& x,   const Matrix<Real>& z, \
    const Matrix<Real>& rmu, const Matrix<Real>& d, \
          Matrix<Real>& dx,        Matrix<Real>& dy, \
          Matrix<Real>& dz ); \
  template void ExpandAugmentedSolution \
  ( const AbstractDistMatrix<Real>& x,   const AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& d, \
          AbstractDistMatrix<Real>& dx,        AbstractDistMatrix<Real>& dy, \
          AbstractDistMatrix<Real>& dz ); \
  template void ExpandAugmentedSolution \
  ( const DistMultiVec<Real>& x,   const DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& d, \
          DistMultiVec<Real>& dx,        DistMultiVec<Real>& dy, \
          DistMultiVec<Real>& dz );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace direct
} // namespace qp
} // namespace El
