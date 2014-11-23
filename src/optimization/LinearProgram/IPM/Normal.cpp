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
//    J = | A X*inv(S) A^T |, and 
//    y = [ b - A inv(S) ( X r_c + tau e ) ]
// where 
//    X   = diag(x),
//    S   = diag(s),
//    e   = ones(n,1),
//    r_b = A x - b, and
//    r_c = A^T l + s - c.
//
// The implied system is of the form
//
//   J | \Delta l | = y,
// 
//  \Delta s = -r_c - A^T \Delta l, and
//  \Delta x = -x + inv(S)(tau e - X \Delta s).
//

template<typename Real>
void FormNormalSystem
( const SparseMatrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& s, const Matrix<Real>& x, const Matrix<Real>& l,
  Real tau, SparseMatrix<Real>& J, Matrix<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::FormNormalSystem"))
    const Int n = A.Width();

    // Form X^{1/2} / S^{1/2}
    // ======================
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
void FormNormalSystem
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& l,
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
( const SparseMatrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& s, const Matrix<Real>& x, const Matrix<Real>& l,
  Real tau, const SparseMatrix<Real>& J, const Matrix<Real>& y,
  Matrix<Real>& ds, Matrix<Real>& dx, Matrix<Real>& dl )
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
void SolveNormalSystem
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& l,
  Real tau, const DistSparseMatrix<Real>& J, const DistMultiVec<Real>& y,
  DistMultiVec<Real>& ds, DistMultiVec<Real>& dx, DistMultiVec<Real>& dl )
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

#define PROTO(Real) \
  template void FormNormalSystem \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& s, const Matrix<Real>& x, const Matrix<Real>& l, \
    Real tau, SparseMatrix<Real>& J, Matrix<Real>& y ); \
  template void FormNormalSystem \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& l, \
    Real tau, DistSparseMatrix<Real>& J, DistMultiVec<Real>& y ); \
  template void SolveNormalSystem \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& s, const Matrix<Real>& x, const Matrix<Real>& l, \
    Real tau, const SparseMatrix<Real>& J, const Matrix<Real>& y, \
    Matrix<Real>& dx, Matrix<Real>& dl, Matrix<Real>& ds ); \
  template void SolveNormalSystem \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& s, const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& l, \
    Real tau, const DistSparseMatrix<Real>& J, const DistMultiVec<Real>& y, \
    DistMultiVec<Real>& ds, DistMultiVec<Real>& dx, DistMultiVec<Real>& dl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace lin_prog
} // namespace El
