/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace lp {
namespace direct {

// Form 
//    J = | A (z <> x) A^T |, and 
//    y = [ r_b - A (z <> (x o r_c + r_mu)) ]
// where 
//    r_b  = A x - b, 
//    r_c  = A^T y - z + c, and
//    r_mu = x o z - tau e.
//
// The implied system is of the form
//
//   J | dy | = d,
// 
//  dz = r_c + A^T dy, and
//  dx = -z <> (r_mu + x o dz)
//

template<typename Real>
void NormalKKT
( const Matrix<Real>& A,
  const Matrix<Real>& x, const Matrix<Real>& z,
        Matrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::NormalKKT"))

    // d := sqrt(z) <> sqrt(x)
    // =======================
    const Int n = A.Width();
    Matrix<Real> d;
    d.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        d.Set( i, 0, Sqrt(x.Get(i,0))/Sqrt(z.Get(i,0)) );

    // Form the Jacobian, J := A D^2 A^T
    // =================================
    auto AD( A );
    DiagonalScale( RIGHT, NORMAL, d, AD );
    Syrk( LOWER, NORMAL, Real(1), AD, J );
    if( !onlyLower )
        MakeSymmetric( LOWER, J );
}

template<typename Real>
void NormalKKT
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::NormalKKT"))

    // From sqrt(x) and sqrt(z)
    // ========================
    DistMatrix<Real,MR,STAR> xSqrt(A.Grid()), zSqrt(A.Grid());
    xSqrt.Align(0,0);
    zSqrt.Align(0,0);
    xSqrt = x;
    zSqrt = z;
    EntrywiseMap( xSqrt, function<Real(Real)>(sqrt) );
    EntrywiseMap( zSqrt, function<Real(Real)>(sqrt) );

    // Form the Jacobian, J := A D^2 A^T
    // =================================
    DistMatrix<Real,MC,MR> AD( A );
    DiagonalScale( RIGHT, NORMAL, xSqrt, AD );
    DiagonalSolve( RIGHT, NORMAL, zSqrt, AD );
    Syrk( LOWER, NORMAL, Real(1), AD, J );
    if( !onlyLower )
        MakeSymmetric( LOWER, J );
}

template<typename Real>
void NormalKKT
( const SparseMatrix<Real>& A,
  const Matrix<Real>& x,       const Matrix<Real>& z,
        SparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::NormalKKT"))

    // d := sqrt(z) <> sqrt(x)
    // =======================
    const Int n = A.Width();
    Matrix<Real> d;
    d.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        d.Set( i, 0, Sqrt(x.Get(i,0))/Sqrt(z.Get(i,0)) );

    // Form the Jacobian, J := A D^2 A^T
    // =================================
    SparseMatrix<Real> G;
    Transpose( A, G );
    DiagonalScale( LEFT, NORMAL, d, G );
    Syrk( LOWER, TRANSPOSE, Real(1), G, J );
    if( !onlyLower )
        MakeSymmetric( LOWER, J );
}

template<typename Real>
void NormalKKT
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z, 
        DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::NormalKKT"))
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    if( !mpi::Congruent( comm, x.Comm() ) )
        LogicError("Communicators of A and x must match");
    if( !mpi::Congruent( comm, z.Comm() ) )
        LogicError("Communicators of A and z must match");

    // d := sqrt(z) <> sqrt(x)
    // =======================
    DistMultiVec<Real> d(comm);
    d.Resize( n, 1 );
    for( Int iLoc=0; iLoc<d.LocalHeight(); ++iLoc )
        d.SetLocal
        ( iLoc, 0, Sqrt(x.GetLocal(iLoc,0))/Sqrt(z.GetLocal(iLoc,0)) );

    // Form the Jacobian, J := A D^2 A^T
    // =================================
    DistSparseMatrix<Real> G(comm);
    Transpose( A, G );
    DiagonalScale( LEFT, NORMAL, d, G );
    Syrk( LOWER, TRANSPOSE, Real(1), G, J );
    if( !onlyLower )
        MakeSymmetric( LOWER, J );
}

template<typename Real>
void NormalKKTRHS
( const Matrix<Real>& A,
  const Matrix<Real>& x,  const Matrix<Real>& z, 
  const Matrix<Real>& rc, const Matrix<Real>& rb, 
  const Matrix<Real>& rmu,
        Matrix<Real>& d )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::NormalKKTRHS"))

    // Form the portion of the right-hand side to be multiplied by A
    // =============================================================
    Matrix<Real> g( rc );
    DiagonalScale( LEFT, NORMAL, x, g );
    Axpy( Real(1), rmu, g );
    DiagonalSolve( LEFT, NORMAL, z, g );

    // Form the right-hand side, d
    // ===========================
    d = rb;
    Gemv( NORMAL, Real(-1), A, g, Real(1), d );
}

template<typename Real>
void NormalKKTRHS
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& x,  const AbstractDistMatrix<Real>& z, 
  const AbstractDistMatrix<Real>& rc, const AbstractDistMatrix<Real>& rb, 
  const AbstractDistMatrix<Real>& rmu,        
        AbstractDistMatrix<Real>& d )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::NormalKKTRHS"))

    // Form the portion of the right-hand side to be multiplied by A
    // =============================================================
    DistMatrix<Real,MC,MR> g( rc );
    DiagonalScale( LEFT, NORMAL, x, g );
    Axpy( Real(1), rmu, g );
    DiagonalSolve( LEFT, NORMAL, z, g );

    // Form the right-hand side, d 
    // ===========================
    Copy( rb, d );
    Gemv( NORMAL, Real(-1), A, g, Real(1), d );
}

template<typename Real>
void NormalKKTRHS
( const SparseMatrix<Real>& A,
  const Matrix<Real>& x,       const Matrix<Real>& z, 
  const Matrix<Real>& rc,      const Matrix<Real>& rb, 
  const Matrix<Real>& rmu,
        Matrix<Real>& d )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::NormalKKTRHS"))

    // Form the portion of the right-hand side to be multiplied by A
    // =============================================================
    Matrix<Real> g( rc );
    DiagonalScale( LEFT, NORMAL, x, g );
    Axpy( Real(1), rmu, g );
    DiagonalSolve( LEFT, NORMAL, z, g );

    // Form the right-hand side, d 
    // ===========================
    d = rb;
    Multiply( NORMAL, Real(-1), A, g, Real(1), d );
}

template<typename Real>
void NormalKKTRHS
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z, 
  const DistMultiVec<Real>& rc,    const DistMultiVec<Real>& rb, 
  const DistMultiVec<Real>& rmu,
        DistMultiVec<Real>& d )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::NormalKKTRHS"))
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    if( !mpi::Congruent( comm, rmu.Comm() ) )
        LogicError("Communicators of A and r_mu must match");
    if( !mpi::Congruent( comm, rc.Comm() ) )
        LogicError("Communicators of A and r_c must match");
    if( !mpi::Congruent( comm, rb.Comm() ) )
        LogicError("Communicators of A and r_b must match");
    if( !mpi::Congruent( comm, x.Comm() ) )
        LogicError("Communicators of A and x must match");
    if( !mpi::Congruent( comm, z.Comm() ) )
        LogicError("Communicators of A and z must match");

    // Form the portion of the right-hand side to be multiplied by A
    // =============================================================
    DistMultiVec<Real> g(comm);
    g.Resize( n, 1 );
    for( Int iLoc=0; iLoc<g.LocalHeight(); ++iLoc )
    {
        const Real rmu_i = rmu.GetLocal(iLoc,0);
        const Real rc_i = rc.GetLocal(iLoc,0);
        const Real x_i = x.GetLocal(iLoc,0);
        const Real z_i = z.GetLocal(iLoc,0);
        g.SetLocal( iLoc, 0, (x_i*rc_i+rmu_i)/z_i );
    }

    // Form the right-hand side, d
    // ===========================
    d = rb;
    Multiply( NORMAL, Real(-1), A, g, Real(1), d );
}

template<typename Real>
void ExpandNormalSolution
( const Matrix<Real>& A,  const Matrix<Real>& c,
  const Matrix<Real>& x,  const Matrix<Real>& z,
  const Matrix<Real>& rc, const Matrix<Real>& rmu,
        Matrix<Real>& dx, const Matrix<Real>& dy, 
        Matrix<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::ExpandNormalSolution"))

    // dz := r_c + A^T dy
    // ==================
    dz = rc;
    Gemv( TRANSPOSE, Real(1), A, dy, Real(1), dz );

    // dx := -z <> (r_mu + x o dz)
    // ===========================
    dx = dz;
    DiagonalScale( LEFT, NORMAL, x, dx );
    Axpy( Real(1), rmu, dx );
    DiagonalSolve( LEFT, NORMAL, z, dx );
    Scale( Real(-1), dx );
}

template<typename Real>
void ExpandNormalSolution
( const AbstractDistMatrix<Real>& A,  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& x,  const AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Real>& rc, const AbstractDistMatrix<Real>& rmu,
        AbstractDistMatrix<Real>& dx, const AbstractDistMatrix<Real>& dy, 
        AbstractDistMatrix<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::ExpandNormalSolution"))

    // dz := r_c + A^T dy
    // ==================
    Copy( rc, dz );
    Gemv( TRANSPOSE, Real(1), A, dy, Real(1), dz );

    // dx := -z <> (r_mu + x o dz)
    // ===========================
    Copy( dz, dx );
    DiagonalScale( LEFT, NORMAL, x, dx );
    Axpy( Real(1), rmu, dx );
    DiagonalSolve( LEFT, NORMAL, z, dx );
    Scale( Real(-1), dx );
}

template<typename Real>
void ExpandNormalSolution
( const SparseMatrix<Real>& A, const Matrix<Real>& c,
  const Matrix<Real>& x,       const Matrix<Real>& z,
  const Matrix<Real>& rc,      const Matrix<Real>& rmu,
        Matrix<Real>& dx,      const Matrix<Real>& dy, 
        Matrix<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::ExpandNormalSolution"))

    // dz := r_c + A^T dy
    // ==================
    dz = rc;
    Multiply( TRANSPOSE, Real(1), A, dy, Real(1), dz );

    // dx := -z <> (r_mu + x o dz)
    // ===========================
    dx = dz;
    DiagonalScale( LEFT, NORMAL, x, dx );
    Axpy( Real(1), rmu, dx );
    DiagonalSolve( LEFT, NORMAL, z, dx );
    Scale( Real(-1), dx );
}

template<typename Real>
void ExpandNormalSolution
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rc,    const DistMultiVec<Real>& rmu,
        DistMultiVec<Real>& dx,    const DistMultiVec<Real>& dy, 
        DistMultiVec<Real>& dz )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::ExpandNormalSolution"))

    // dz := r_c + A^T dy
    // ==================
    dz = rc;
    Multiply( TRANSPOSE, Real(1), A, dy, Real(1), dz );

    // dx := -z <> (r_mu + x o dz)
    // ===========================
    dx = dz;
    DiagonalScale( LEFT, NORMAL, x, dx );
    Axpy( Real(1), rmu, dx );
    DiagonalSolve( LEFT, NORMAL, z, dx );
    Scale( Real(-1), dx );
}

#define PROTO(Real) \
  template void NormalKKT \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& x, const Matrix<Real>& z, \
          Matrix<Real>& J, bool onlyLower ); \
  template void NormalKKT \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& J, bool onlyLower ); \
  template void NormalKKT \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& x,       const Matrix<Real>& z, \
          SparseMatrix<Real>& J, bool onlyLower ); \
  template void NormalKKT \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J, bool onlyLower ); \
  template void NormalKKTRHS \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& x,  const Matrix<Real>& z, \
    const Matrix<Real>& rc, const Matrix<Real>& rb, \
    const Matrix<Real>& rmu, \
          Matrix<Real>& d ); \
  template void NormalKKTRHS \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& x,  const AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Real>& rc, const AbstractDistMatrix<Real>& rb, \
    const AbstractDistMatrix<Real>& rmu, \
          AbstractDistMatrix<Real>& d ); \
  template void NormalKKTRHS \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& x,       const Matrix<Real>& z, \
    const Matrix<Real>& rc,      const Matrix<Real>& rb, \
    const Matrix<Real>& rmu, \
          Matrix<Real>& d ); \
  template void NormalKKTRHS \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& x,   const DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& rc,  const DistMultiVec<Real>& rb, \
    const DistMultiVec<Real>& rmu, \
          DistMultiVec<Real>& d ); \
  template void ExpandNormalSolution \
  ( const Matrix<Real>& A,  const Matrix<Real>& c, \
    const Matrix<Real>& x,  const Matrix<Real>& z, \
    const Matrix<Real>& rc, const Matrix<Real>& rmu, \
          Matrix<Real>& dx, const Matrix<Real>& dy, \
          Matrix<Real>& dz ); \
  template void ExpandNormalSolution \
  ( const AbstractDistMatrix<Real>& A,  const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& x,  const AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Real>& rc, const AbstractDistMatrix<Real>& rmu, \
          AbstractDistMatrix<Real>& dx, const AbstractDistMatrix<Real>& dy, \
          AbstractDistMatrix<Real>& dz ); \
  template void ExpandNormalSolution \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& c, \
    const Matrix<Real>& x,       const Matrix<Real>& z, \
    const Matrix<Real>& rc,      const Matrix<Real>& rmu, \
          Matrix<Real>& dx,      const Matrix<Real>& dy, \
          Matrix<Real>& dz ); \
  template void ExpandNormalSolution \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& rc,    const DistMultiVec<Real>& rmu, \
          DistMultiVec<Real>& dx,    const DistMultiVec<Real>& dy, \
          DistMultiVec<Real>& dz );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace direct
} // namespace lp
} // namespace El
