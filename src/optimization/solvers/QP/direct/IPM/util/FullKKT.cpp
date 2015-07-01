/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "../../../affine/IPM/util.hpp"

namespace El {
namespace qp {
namespace direct {

// The full KKT system is of the form
//
//   |  Q A^T     -I    | | x |   |          -c            |
//   |  A 0        0    | | y |   |           b            |,
//   | -I 0   (-z <> x) | | z | = | - z <> (x o z + tau e) |
//
// and the particular system solved is of the form
//
//   |  Q A^T     -I    | | dx |   |   -rc    |
//   |  A 0        0    | | dy |   |   -rb    |,
//   | -I 0   (-z <> x) | | dz | = | z <> rmu |
//
// where 
//
//   rc = Q x + A^T y - z + c,
//   rb = A x - b,
//   rmu = x o z - tau e

template<typename Real>
void KKT
( const Matrix<Real>& Q, const Matrix<Real>& A, 
  const Matrix<Real>& x, const Matrix<Real>& z,
        Matrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("qp::direct::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    Zeros( J, 2*n+m, 2*n+m );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,2*n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd); 
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd); 
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd); 

    // Jxx := Q
    // ========
    Jxx = Q;

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := -I
    // =========
    Identity( Jzx, n, n );
    Jzx *= -1;

    // Jzz := - z <> x
    // ===============
    Matrix<Real> t;
    t = x;
    DiagonalSolve( LEFT, NORMAL, z, t );
    t *= -1;
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy ); 

        // Jxz := -I
        // =========
        Identity( Jxz, n, n );
        Jxz *= -1;
    }
}

template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& Q,    const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& x,    const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("qp::direct::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    auto JPtr = WriteProxy<Real,MC,MR>(&JPre); 
    auto& J = *JPtr;

    Zeros( J, 2*n+m, 2*n+m );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,2*n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd);
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd);

    // Jxx := Q
    // ========
    Jxx = Q;

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := -I
    // =========
    Identity( Jzx, n, n );
    Jzx *= -1;

    // Jzz := - z <> x
    // ===============
    DistMatrix<Real,MC,STAR> t(x);
    DiagonalSolve( LEFT, NORMAL, z, t );
    t *= -1;
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy );

        // Jxz := -I
        // =========
        Identity( Jxz, n, n );
        Jxz *= -1;
    }
}

template<typename Real>
void KKT
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, 
  const Matrix<Real>& x,       const Matrix<Real>& z,
        SparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("qp::direct::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    Zeros( J, 2*n+m, 2*n+m );
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

    if( onlyLower )
        J.Reserve( numUsedEntriesQ + numEntriesA + 2*n );
    else
        J.Reserve( numUsedEntriesQ + 2*numEntriesA + 3*n );

    // Jxx = Q
    // =======
    for( Int e=0; e<numEntriesQ; ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            J.QueueUpdate( Q.Row(e), Q.Col(e), Q.Value(e) );
    }

    // Jyx = A
    // =======
    for( Int e=0; e<numEntriesA; ++e )
        J.QueueUpdate( n+A.Row(e), A.Col(e), A.Value(e) );

    // Jzx = -I
    // ========
    for( Int e=0; e<n; ++e )
        J.QueueUpdate( n+m+e, e, Real(-1) );

    // Jzz = - z <> x
    // ==============
    for( Int e=0; e<n; ++e )
        J.QueueUpdate( n+m+e, n+m+e, -x.Get(e,0)/z.Get(e,0) );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        for( Int e=0; e<numEntriesA; ++e )
            J.QueueUpdate( A.Col(e), n+A.Row(e), A.Value(e) );

        // Jxz := -I
        // =========
        for( Int e=0; e<n; ++e )
            J.QueueUpdate( e, n+m+e, Real(-1) );
    }
    J.ProcessQueues();
}

template<typename Real>
void KKT
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("qp::direct::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    J.SetComm( A.Comm() );
    Zeros( J, m+2*n, m+2*n );

    // Count the number of entries to send
    // ===================================
    Int numEntries = 0;
    for( Int e=0; e<Q.NumLocalEntries(); ++e )
    {
        const Int i = Q.Row(e);
        const Int j = Q.Col(e);
        if( i >= j || !onlyLower )
            ++numEntries;
    }
    numEntries += A.NumLocalEntries();
    if( !onlyLower )
        numEntries += A.NumLocalEntries();
    numEntries += x.LocalHeight();
    // Count the total number of entries for the negative identities
    // -------------------------------------------------------------
    Int negIdentUpdates = 0;
    for( Int iLoc=0; iLoc<J.LocalHeight(); ++iLoc )
    {
        const Int i = J.GlobalRow(iLoc);
        if( i < n && !onlyLower )
            ++negIdentUpdates;
        else if( i >= n+m )
            ++negIdentUpdates;
    }

    // Pack and process the updates
    // ============================
    J.Reserve( numEntries+negIdentUpdates, numEntries );
    // Append the local negative identity updates
    // ------------------------------------------
    for( Int iLoc=0; iLoc<J.LocalHeight(); ++iLoc )
    {
        const Int i = J.GlobalRow(iLoc);
        if( (i < n && !onlyLower) || i >= n+m )
            J.QueueUpdate( i, i+(n+m), Real(-1) );
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
    // Pack -z <> x
    // ------------
    for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
    {
        const Int i = m+n + x.GlobalRow(iLoc);
        const Real value = -x.GetLocal(iLoc,0)/z.GetLocal(iLoc,0);
        J.QueueUpdate( i, i, value );
    }
    J.ProcessQueues();
}

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc,  const Matrix<Real>& rb, 
  const Matrix<Real>& rmu, const Matrix<Real>& z, 
        Matrix<Real>& d )
{
    DEBUG_ONLY(CSE cse("qp::direct::KKTRHS"))
    const Int m = rb.Height();
    const Int n = rc.Height();
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,2*n+m);
    Zeros( d, 2*n+m, 1 );

    auto dx = d(xInd,ALL);
    dx = rc;
    dx *= -1;

    auto dy = d(yInd,ALL);
    dy = rb;
    dy *= -1;

    auto dz = d(zInd,ALL);
    dz = rmu;
    DiagonalSolve( LEFT, NORMAL, z, dz );
}

template<typename Real>
void KKTRHS
( const AbstractDistMatrix<Real>& rc,  const AbstractDistMatrix<Real>& rb, 
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& z, 
        AbstractDistMatrix<Real>& dPre )
{
    DEBUG_ONLY(CSE cse("qp::direct::KKTRHS"))

    auto dPtr = WriteProxy<Real,MC,MR>(&dPre); 
    auto& d = *dPtr;

    const Int m = rb.Height();
    const Int n = rc.Height();
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,2*n+m);
    Zeros( d, 2*n+m, 1 );

    auto dx = d(xInd,ALL);
    dx = rc;
    dx *= -1;

    auto dy = d(yInd,ALL);
    dy = rb;
    dy *= -1;

    auto dz = d(zInd,ALL);
    dz = rmu;
    DiagonalSolve( LEFT, NORMAL, z, dz );
}

template<typename Real>
void KKTRHS
( const DistMultiVec<Real>& rc,  const DistMultiVec<Real>& rb, 
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& z, 
        DistMultiVec<Real>& d )
{
    DEBUG_ONLY(CSE cse("qp::direct::KKTRHS"))
    const Int m = rb.Height();
    const Int n = rc.Height();
    d.SetComm( rmu.Comm() );
    Zeros( d, m+2*n, 1 );

    d.Reserve( rc.LocalHeight()+rb.LocalHeight()+rmu.LocalHeight() );
    for( Int iLoc=0; iLoc<rc.LocalHeight(); ++iLoc )
    {
        const Int i = rc.GlobalRow(iLoc);
        const Real value = -rc.GetLocal(iLoc,0);
        d.QueueUpdate( i, 0, value );
    }
    for( Int iLoc=0; iLoc<rb.LocalHeight(); ++iLoc )
    {
        const Int i = n + rb.GlobalRow(iLoc);
        const Real value = -rb.GetLocal(iLoc,0);
        d.QueueUpdate( i, 0, value );
    }
    for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
    {
        const Int i = n+m + rmu.GlobalRow(iLoc);
        const Real value = rmu.GetLocal(iLoc,0)/z.GetLocal(iLoc,0);
        d.QueueUpdate( i, 0, value );
    }
    d.ProcessQueues();
}

template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const Matrix<Real>& d, 
        Matrix<Real>& dx, Matrix<Real>& dy, 
        Matrix<Real>& dz )
{
    DEBUG_ONLY(CSE cse("qp::direct::ExpandSolution"))
    qp::affine::ExpandCoreSolution( m, n, n, d, dx, dy, dz );
}

template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const AbstractDistMatrix<Real>& d, 
        AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, 
        AbstractDistMatrix<Real>& dz )
{
    DEBUG_ONLY(CSE cse("qp::direct::ExpandSolution"))
    qp::affine::ExpandCoreSolution( m, n, n, d, dx, dy, dz );
}

template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const DistMultiVec<Real>& d, 
        DistMultiVec<Real>& dx, DistMultiVec<Real>& dy, 
        DistMultiVec<Real>& dz )
{
    DEBUG_ONLY(CSE cse("qp::direct::ExpandSolution"))
    qp::affine::ExpandCoreSolution( m, n, n, d, dx, dy, dz );
}

#define PROTO(Real) \
  template void KKT \
  ( const Matrix<Real>& Q, const Matrix<Real>& A, \
    const Matrix<Real>& x, const Matrix<Real>& z, \
          Matrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, \
    const Matrix<Real>& x,       const Matrix<Real>& z, \
          SparseMatrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J, bool onlyLower ); \
  template void KKTRHS \
  ( const Matrix<Real>& rc,  const Matrix<Real>& rb, \
    const Matrix<Real>& rmu, const Matrix<Real>& z, \
          Matrix<Real>& d ); \
  template void KKTRHS \
  ( const AbstractDistMatrix<Real>& rc,  const AbstractDistMatrix<Real>& rb, \
    const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& d ); \
  template void KKTRHS \
  ( const DistMultiVec<Real>& rc,  const DistMultiVec<Real>& rb, \
    const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& z, \
          DistMultiVec<Real>& d ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const Matrix<Real>& d, \
          Matrix<Real>& dx, Matrix<Real>& dy, \
          Matrix<Real>& dz ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const AbstractDistMatrix<Real>& d, \
          AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, \
          AbstractDistMatrix<Real>& dz ); \
  template void ExpandSolution \
  ( Int m, Int n, \
    const DistMultiVec<Real>& d, \
          DistMultiVec<Real>& dx, DistMultiVec<Real>& dy, \
          DistMultiVec<Real>& dz );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace direct
} // namespace qp
} // namespace El
