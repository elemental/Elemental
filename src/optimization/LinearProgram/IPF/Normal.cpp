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
//    J = | A X inv(S) A^T |, and 
//    y = [ -r_b + A inv(S) (r_mu - X r_c) ]
// where 
//    X   = diag(x),
//    S   = diag(s),
//    e   = ones(n,1),
//    r_b = A x - b, 
//    r_c = A^T l + s - c, and
//    r_mu = X S e - tau e.
//
// The implied system is of the form
//
//   J | dl | = y,
// 
//  ds = -r_c - A^T dl, and
//  dx = -(r_mu - X ds) / S
//

template<typename Real>
void NormalKKT
( const SparseMatrix<Real>& A,
  const Matrix<Real>& s, const Matrix<Real>& x,
  SparseMatrix<Real>& J )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::NormalKKT"))

    // Form X^{1/2} / S^{1/2}
    // ======================
    const Int n = A.Width();
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
}

template<typename Real>
void NormalKKTRHS
( const SparseMatrix<Real>& A,
  const Matrix<Real>& s, const Matrix<Real>& x, 
  const Matrix<Real>& rmu, const Matrix<Real>& rc, const Matrix<Real>& rb,
        Matrix<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::NormalKKTRHS"))

    // Form the portion of the right-hand side to be multiplied by A
    // =============================================================
    const Int n = A.Width();
    Matrix<Real> g;
    g.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        g.Set( i, 0, (rmu.Get(i,0)-x.Get(i,0)*rc.Get(i,0))/s.Get(i,0) );

    // Form the right-hand side, y
    // ===========================
    y = rb;
    Multiply( NORMAL, Real(1), A, g, Real(-1), y );
}

template<typename Real>
void NormalKKT
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, 
  DistSparseMatrix<Real>& J )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::NormalKKT"))
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    if( !mpi::Congruent( comm, x.Comm() ) )
        LogicError("Communicators of A and x must match");
    if( !mpi::Congruent( comm, s.Comm() ) )
        LogicError("Communicators of A and s must match");

    // Form X^{1/2} / S^{1/2}
    // ======================
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
}

template<typename Real>
void NormalKKTRHS
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rc, 
  const DistMultiVec<Real>& rb,
        DistMultiVec<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::NormalKKTRHS"))
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    if( !mpi::Congruent( comm, rmu.Comm() ) )
        LogicError("Communicators of A and r_mu must match");
    if( !mpi::Congruent( comm, rc.Comm() ) )
        LogicError("Communicators of A and r_c must match");
    if( !mpi::Congruent( comm, rb.Comm() ) )
        LogicError("Communicators of A and r_b must match");
    if( !mpi::Congruent( comm, s.Comm() ) )
        LogicError("Communicators of A and s must match");
    if( !mpi::Congruent( comm, x.Comm() ) )
        LogicError("Communicators of A and x must match");

    // Form the portion of the right-hand side to be multiplied by A
    // =============================================================
    DistMultiVec<Real> g(comm);
    g.Resize( n, 1 );
    for( Int iLoc=0; iLoc<g.LocalHeight(); ++iLoc )
    {
        const Real rmu_i = rmu.GetLocal(iLoc,0);
        const Real rc_i = rc.GetLocal(iLoc,0);
        const Real x_i = x.GetLocal(iLoc,0);
        const Real s_i = s.GetLocal(iLoc,0);
        g.SetLocal( iLoc, 0, (rmu_i-x_i*rc_i)/s_i );
    }

    // Form the right-hand side, y
    // ===========================
    y = rb;
    Multiply( NORMAL, Real(1), A, g, Real(-1), y );
}

template<typename Real>
void ExpandNormalSolution
( const SparseMatrix<Real>& A, const Matrix<Real>& c,
  const Matrix<Real>& s,       const Matrix<Real>& x,
  const Matrix<Real>& rmu,     const Matrix<Real>& rc,
  const Matrix<Real>& dl, Matrix<Real>& ds, Matrix<Real>& dx )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::ExpandNormalSolution"))

    // ds := -r_c - A^T dl
    // ===================
    ds = rc;
    Multiply( TRANSPOSE, Real(-1), A, dl, Real(-1), ds );

    // dx := -(r_mu + X ds) / S
    // ========================
    const Int n = A.Width();
    Zeros( dx, n, 1 );
    for( Int i=0; i<n; ++i )
    {
        const Real s_i = s.Get(i,0);
        const Real x_i = x.Get(i,0);
        const Real rmu_i = rmu.Get(i,0);
        const Real ds_i = ds.Get(i,0);
        dx.Set( i, 0, -(rmu_i+x_i*ds_i)/s_i );
    }
}

template<typename Real>
void ExpandNormalSolution
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& s,     const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& rmu,   const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& dl,
  DistMultiVec<Real>& ds, DistMultiVec<Real>& dx )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::ExpandNormalSolution"))

    // ds := -r_c - A^T dl
    // ===================
    ds = rc;
    Multiply( TRANSPOSE, Real(-1), A, dl, Real(-1), ds );

    // dx := -(r_mu + X ds) / S
    // ========================
    const Int n = A.Width();
    Zeros( dx, n, 1 );
    for( Int iLoc=0; iLoc<dx.LocalHeight(); ++iLoc )
    {
        const Real s_i = s.GetLocal(iLoc,0);
        const Real x_i = x.GetLocal(iLoc,0);
        const Real rmu_i = rmu.GetLocal(iLoc,0);
        const Real ds_i = ds.GetLocal(iLoc,0);
        dx.SetLocal( iLoc, 0, -(rmu_i+x_i*ds_i)/s_i );
    }
}

#define PROTO(Real) \
  template void NormalKKT \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& s, const Matrix<Real>& x, \
    SparseMatrix<Real>& J ); \
  template void NormalKKT \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, \
    DistSparseMatrix<Real>& J ); \
  template void NormalKKTRHS \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& s, const Matrix<Real>& x, \
    const Matrix<Real>& rmu, const Matrix<Real>& rc, const Matrix<Real>& rb, \
          Matrix<Real>& y ); \
  template void NormalKKTRHS \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rc, \
    const DistMultiVec<Real>& rb, DistMultiVec<Real>& y ); \
  template void ExpandNormalSolution \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& c, \
    const Matrix<Real>& s,       const Matrix<Real>& x, \
    const Matrix<Real>& rmu,     const Matrix<Real>& rc, \
    const Matrix<Real>& dl, Matrix<Real>& ds, Matrix<Real>& dx ); \
  template void ExpandNormalSolution \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& s,     const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& rmu,  const DistMultiVec<Real>& rc, \
    const DistMultiVec<Real>& dl, \
    DistMultiVec<Real>& ds, DistMultiVec<Real>& dx );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace lin_prog
} // namespace El
