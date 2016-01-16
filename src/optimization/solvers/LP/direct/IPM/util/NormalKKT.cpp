/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace lp {
namespace direct {

// Starting with
//
//   | gamma^2 I      A^T          -I    | | dx |   |     r_1     |
//   |     A      -delta^2 I        0    | | dy | = |     r_2     |,
//   |    -I           0       -inv(Z) X | | dz |   | -inv(Z) r_3 |
//
// using -inv(Z) X as a block pivot yields the system
//
//   | inv(X) Z + gamma^2 I      A^T    | | dx | = | r_1 + inv(X) r_3 | ,
//   |          A            -delta^2 I | | dy |   |        r_2       |
//
// where we can recover 
//
//   dz = gamma^2 dx + A^T dy - r_1. 
//
// Subsequently using inv(X) Z + gamma^2 I, which is convenient to define as 
// inv(D)^2, as a block pivot yields
//
//   (A D^2 A^T + delta^2 I) dy = A D^2 (r_1 + inv(X) r_3) - r_2,
//
// where we can recover 
//
//   dx = D^2 (-A^T dy + r_1 + inv(X) r_3).
//
// Please see [1] for detailed notes on this derivation. 
//
// NOTE: At this time, r_1 = -r_c, r_2 = -r_b, r_3 = -r_mu.
//
// Also, in order to compensate for the floating-point error in the formation
// of A D^2 A^T, the diagonal entries of A D^2 A^T are set to 1+1e-14 times 
// their absolute value. This approach is due to a personal communication with
// Tor Mykle
//
// [1] Michael Saunders, "Notes 7: PDCO - Primal-Dual Interior Methods",
//     from the course notes of Large-Scale Numerical Optimization, 
//     Stanford University, Management Science and Engineering (and ICME), 2015.
//     Last accessed from: 
//     http://web.stanford.edu/class/msande318/notes/notes07-PDinterior.pdf
//

// TODO: Avoid reforming D

template<typename Real>
void NormalKKT
( const Matrix<Real>& A,
        Real gamma, 
        Real delta,
  const Matrix<Real>& x, 
  const Matrix<Real>& z,
        Matrix<Real>& J, 
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::direct::NormalKKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    // dInv := sqrt( (z ./ x) .+ gamma^2 )
    // ===================================
    Matrix<Real> dInv;
    dInv.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        dInv.Set( i, 0, Sqrt(z.Get(i,0)/x.Get(i,0) + gamma*gamma) );

    // Form A D^2 A^T + delta^2 I
    // ==========================
    auto AD( A );
    DiagonalSolve( RIGHT, NORMAL, dInv, AD );
    Zeros( J, m, m );
    ShiftDiagonal( J, delta*delta );
    Syrk( LOWER, NORMAL, Real(1), AD, Real(1), J );

    // TODO: Inflate diagonal?

    if( !onlyLower )
        MakeSymmetric( LOWER, J );
}

template<typename Real>
void NormalKKT
( const ElementalMatrix<Real>& A,
        Real gamma, 
        Real delta,
  const ElementalMatrix<Real>& xPre,
  const ElementalMatrix<Real>& zPre,
        ElementalMatrix<Real>& J,
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::direct::NormalKKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    ElementalProxyCtrl ctrl;
    ctrl.colAlign = 0;
    ctrl.colConstrain = true;

    DistMatrixReadProxy<Real,Real,MR,STAR>
      xProx( xPre, ctrl ),
      zProx( zPre, ctrl );
    auto& x = xProx.GetLocked();
    auto& z = zProx.GetLocked();

    DistMatrix<Real,MR,STAR> dInv(A.Grid());
    dInv.Resize( n, 1 );
    {
        const Int nLocal = dInv.LocalHeight();
        for( Int iLoc=0; iLoc<nLocal; ++iLoc )
            dInv.SetLocal
            ( iLoc, 0, 
              Sqrt(z.GetLocal(iLoc,0)/x.GetLocal(iLoc,0)+gamma*gamma) );
    }

    // Form A D^2 A^T + delta^2 I
    // ==========================
    DistMatrix<Real,MC,MR> AD( A );
    DiagonalSolve( RIGHT, NORMAL, dInv, AD );
    Zeros( J, m, m );
    ShiftDiagonal( J, delta*delta );
    Syrk( LOWER, NORMAL, Real(1), AD, Real(1), J );

    // TODO: Inflate diagonal?

    if( !onlyLower )
        MakeSymmetric( LOWER, J );
}

template<typename Real>
void NormalKKT
( const SparseMatrix<Real>& A,
        Real gamma, 
        Real delta,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J, 
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::direct::NormalKKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    // TODO: Expose this value as a parameter
    const Real inflateRatio = Pow(limits::Epsilon<Real>(),Real(0.83));

    // dInv := sqrt( (z ./ x) .+ gamma^2 )
    // ===================================
    Matrix<Real> dInv;
    dInv.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        dInv.Set( i, 0, Sqrt(z.Get(i,0)/x.Get(i,0) + gamma*gamma) );

    // Form A D^2 A^T + delta^2 I
    // ==========================
    SparseMatrix<Real> G;
    // TODO: Avoid forming this within an inner loop...
    Transpose( A, G );
    DiagonalSolve( LEFT, NORMAL, dInv, G );
    Zeros( J, m, m );
    ShiftDiagonal( J, delta*delta );
    Syrk( LOWER, TRANSPOSE, Real(1), G, Real(1), J );

    // Inflate the diagonal in a small relative sense
    // ==============================================
    // TODO: Create EntrywiseMapDiagonal and replace this with it
    Real* valBuf = J.ValueBuffer();
    for( Int i=0; i<m; ++i )
    {
        const Int e = J.Offset( i, i );
        const Real diagAbs = Abs(valBuf[e]);
        valBuf[e] = (1+inflateRatio)*diagAbs;
    }

    if( !onlyLower )
        MakeSymmetric( LOWER, J );
}

template<typename Real>
void NormalKKT
( const DistSparseMatrix<Real>& A, 
        Real gamma,
        Real delta,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z, 
        DistSparseMatrix<Real>& J, 
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::direct::NormalKKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    if( !mpi::Congruent( comm, x.Comm() ) )
        LogicError("Communicators of A and x must match");
    if( !mpi::Congruent( comm, z.Comm() ) )
        LogicError("Communicators of A and z must match");

    // TODO: Expose this value as a parameter
    const Real inflateRatio = Pow(limits::Epsilon<Real>(),Real(0.83));

    // dInv := sqrt( (z ./ x) .+ gamma^2 )
    // ===================================
    DistMultiVec<Real> dInv(comm);
    dInv.Resize( n, 1 );
    const Int dInvLocalHeight = dInv.LocalHeight();
    for( Int iLoc=0; iLoc<dInvLocalHeight; ++iLoc )
        dInv.SetLocal
        ( iLoc, 0, Sqrt(z.GetLocal(iLoc,0)/x.GetLocal(iLoc,0) + gamma*gamma) );

    // Form A D^2 A^T + delta^2 I
    // ==========================
    DistSparseMatrix<Real> G(comm);
    // TODO: Avoid forming this within an inner loop...
    Transpose( A, G );
    DiagonalSolve( LEFT, NORMAL, dInv, G );
    Zeros( J, m, m );
    ShiftDiagonal( J, delta*delta );
    Syrk( LOWER, TRANSPOSE, Real(1), G, Real(1), J );

    // Inflate the diagonal in a small relative sense
    // ==============================================
    // TODO: Create EntrywiseMapDiagonal and replace this with it
    Real* valBuf = J.ValueBuffer();
    const Int JLocalHeight = J.LocalHeight();
    for( Int iLoc=0; iLoc<JLocalHeight; ++iLoc )
    {
        const Int i = J.GlobalRow(iLoc);
        const Int e = J.Offset( iLoc, i );
        const Real diagAbs = Abs(valBuf[e]);
        valBuf[e] = (1+inflateRatio)*diagAbs;
    }

    if( !onlyLower )
        MakeSymmetric( LOWER, J );
}

template<typename Real>
void NormalKKTRHS
( const Matrix<Real>& A,
        Real gamma,
  const Matrix<Real>& x,
  const Matrix<Real>& z, 
  const Matrix<Real>& rc,
  const Matrix<Real>& rb, 
  const Matrix<Real>& rmu,
        Matrix<Real>& d )
{
    DEBUG_ONLY(CSE cse("lp::direct::NormalKKTRHS"))
    const Int n = A.Width();

    // dInv := sqrt( (z ./ x) .+ gamma^2 )
    // ===================================
    Matrix<Real> dInv;
    dInv.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        dInv.Set( i, 0, Sqrt(z.Get(i,0)/x.Get(i,0) + gamma*gamma) );

    // g := D^2( r_1 + inv(X) r_3 ) = D^2 ( -r_c - inv(X) r_mu )
    // =========================================================
    Matrix<Real> g( rmu );
    DiagonalSolve( LEFT, NORMAL, x, g );
    g += rc;
    g *= -1;
    DiagonalSolve( LEFT, NORMAL, dInv, g );
    DiagonalSolve( LEFT, NORMAL, dInv, g );

    // d := A D^2 ( r_1 + inv(X) r_3 ) - r_2 = A g + r_b
    // =================================================
    d = rb;
    Gemv( NORMAL, Real(1), A, g, Real(1), d );
}

template<typename Real>
void NormalKKTRHS
( const ElementalMatrix<Real>& A,
        Real gamma,
  const ElementalMatrix<Real>& xPre,
  const ElementalMatrix<Real>& zPre, 
  const ElementalMatrix<Real>& rc,
  const ElementalMatrix<Real>& rb, 
  const ElementalMatrix<Real>& rmu,        
        ElementalMatrix<Real>& d )
{
    DEBUG_ONLY(CSE cse("lp::direct::NormalKKTRHS"))
    const Int n = A.Width();

    ElementalProxyCtrl ctrl;
    ctrl.colAlign = 0;
    ctrl.colConstrain = true;

    DistMatrixReadProxy<Real,Real,MR,STAR>
      xProx( xPre, ctrl ),
      zProx( zPre, ctrl );
    auto& x = xProx.GetLocked();
    auto& z = zProx.GetLocked();

    DistMatrix<Real,MR,STAR> dInv(A.Grid());
    dInv.Resize( n, 1 );
    {
        const Int nLocal = dInv.LocalHeight();
        for( Int iLoc=0; iLoc<nLocal; ++iLoc )
            dInv.SetLocal
            ( iLoc, 0, 
              Sqrt(z.GetLocal(iLoc,0)/x.GetLocal(iLoc,0)+gamma*gamma) );
    }

    // g := D^2( r_1 + inv(X) r_3 ) = D^2 ( -r_c - inv(X) r_mu )
    // =========================================================
    DistMatrix<Real,MC,MR> g( rmu );
    DiagonalSolve( LEFT, NORMAL, x, g );
    g += rc;
    g *= -1;
    DiagonalSolve( LEFT, NORMAL, dInv, g );
    DiagonalSolve( LEFT, NORMAL, dInv, g );

    // d := A D^2 ( r_1 + inv(X) r_3 ) - r_2 = A g + r_b
    // =================================================
    d = rb;
    Gemv( NORMAL, Real(1), A, g, Real(1), d );
}

template<typename Real>
void NormalKKTRHS
( const SparseMatrix<Real>& A,
        Real gamma,
  const Matrix<Real>& x,
  const Matrix<Real>& z, 
  const Matrix<Real>& rc,
  const Matrix<Real>& rb, 
  const Matrix<Real>& rmu,
        Matrix<Real>& d )
{
    DEBUG_ONLY(CSE cse("lp::direct::NormalKKTRHS"))
    const Int n = A.Width();

    // dInv := sqrt( (z ./ x) .+ gamma^2 )
    // ===================================
    Matrix<Real> dInv;
    dInv.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        dInv.Set( i, 0, Sqrt(z.Get(i,0)/x.Get(i,0) + gamma*gamma) );

    // g := D^2( r_1 + inv(X) r_3 ) = D^2 ( -r_c - inv(X) r_mu )
    // =========================================================
    Matrix<Real> g( rmu );
    DiagonalSolve( LEFT, NORMAL, x, g );
    g += rc;
    g *= -1;
    DiagonalSolve( LEFT, NORMAL, dInv, g );
    DiagonalSolve( LEFT, NORMAL, dInv, g );

    // d := A D^2 ( r_1 + inv(X) r_3 ) - r_2 = A g + r_b
    // =================================================
    d = rb;
    Multiply( NORMAL, Real(1), A, g, Real(1), d );
}

template<typename Real>
void NormalKKTRHS
( const DistSparseMatrix<Real>& A,
        Real gamma,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z, 
  const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rb, 
  const DistMultiVec<Real>& rmu,
        DistMultiVec<Real>& d )
{
    DEBUG_ONLY(CSE cse("lp::direct::NormalKKTRHS"))
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

    DistMultiVec<Real> dInv(comm);
    dInv.Resize( n, 1 );
    {
        const Int nLocal = dInv.LocalHeight();
        for( Int iLoc=0; iLoc<nLocal; ++iLoc )
            dInv.SetLocal
            ( iLoc, 0, 
              Sqrt(z.GetLocal(iLoc,0)/x.GetLocal(iLoc,0)+gamma*gamma) );
    }

    // g := D^2( r_1 + inv(X) r_3 ) = D^2 ( -r_c - inv(X) r_mu )
    // =========================================================
    DistMultiVec<Real> g( rmu );
    DiagonalSolve( LEFT, NORMAL, x, g );
    g += rc;
    g *= -1;
    DiagonalSolve( LEFT, NORMAL, dInv, g );
    DiagonalSolve( LEFT, NORMAL, dInv, g );

    // d := A D^2 ( r_1 + inv(X) r_3 ) - r_2 = A g + r_b
    // =================================================
    d = rb;
    Multiply( NORMAL, Real(1), A, g, Real(1), d );
}

template<typename Real>
void ExpandNormalSolution
( const Matrix<Real>& A,
        Real gamma,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
  const Matrix<Real>& rc,
  const Matrix<Real>& rmu,
        Matrix<Real>& dx,
  const Matrix<Real>& dy, 
        Matrix<Real>& dz )
{
    DEBUG_ONLY(CSE cse("lp::direct::ExpandNormalSolution"))
    const Int n = A.Width();

    // dInv := sqrt( (z ./ x) .+ gamma^2 )
    // ===================================
    Matrix<Real> dInv;
    dInv.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        dInv.Set( i, 0, Sqrt(z.Get(i,0)/x.Get(i,0) + gamma*gamma) );

    // Use dz as a temporary for storing A^T dy
    // ========================================
    Zeros( dz, n, 1 );
    Gemv( TRANSPOSE, Real(1), A, dy, Real(0), dz );

    // dx = D^2 (-A^T dy + r_1 + inv(X) r_3) 
    //    = D^2 (-A^T dy - r_c - inv(X) r_mu)
    // ======================================
    dx = rmu;
    DiagonalSolve( LEFT, NORMAL, x, dx );
    dx += rc;
    dx += dz;
    dx *= -1;
    DiagonalSolve( LEFT, NORMAL, dInv, dx );
    DiagonalSolve( LEFT, NORMAL, dInv, dx );

    // dz := r_c + gamma^2 dx + A^T dy
    // ===============================
    Axpy( gamma*gamma, dx, dz );
    dz += rc;
}

template<typename Real>
void ExpandNormalSolution
( const ElementalMatrix<Real>& A,
        Real gamma,
  const ElementalMatrix<Real>& xPre,
  const ElementalMatrix<Real>& zPre,
  const ElementalMatrix<Real>& rc,
  const ElementalMatrix<Real>& rmu,
        ElementalMatrix<Real>& dx,
  const ElementalMatrix<Real>& dy, 
        ElementalMatrix<Real>& dz )
{
    DEBUG_ONLY(CSE cse("lp::direct::ExpandNormalSolution"))
    const Int n = A.Width();

    ElementalProxyCtrl ctrl;
    ctrl.colAlign = 0;
    ctrl.colConstrain = true;

    DistMatrixReadProxy<Real,Real,MR,STAR>
      xProx( xPre, ctrl ),
      zProx( zPre, ctrl );
    auto& x = xProx.GetLocked();
    auto& z = zProx.GetLocked();

    DistMatrix<Real,MR,STAR> dInv(A.Grid());
    dInv.Resize( n, 1 );
    {
        const Int nLocal = dInv.LocalHeight();
        for( Int iLoc=0; iLoc<nLocal; ++iLoc )
            dInv.SetLocal
            ( iLoc, 0, 
              Sqrt(z.GetLocal(iLoc,0)/x.GetLocal(iLoc,0)+gamma*gamma) );
    }

    // Use dz as a temporary for storing A^T dy
    // ========================================
    Zeros( dz, n, 1 );
    Gemv( TRANSPOSE, Real(1), A, dy, Real(0), dz );

    // dx = D^2 (-A^T dy + r_1 + inv(X) r_3) 
    //    = D^2 (-A^T dy - r_c - inv(X) r_mu)
    // ======================================
    dx = rmu;
    DiagonalSolve( LEFT, NORMAL, x, dx );
    dx += rc;
    dx += dz;
    dx *= -1;
    DiagonalSolve( LEFT, NORMAL, dInv, dx );
    DiagonalSolve( LEFT, NORMAL, dInv, dx );

    // dz := r_c + gamma^2 dx + A^T dy
    // ===============================
    Axpy( gamma*gamma, dx, dz );
    dz += rc;
}

template<typename Real>
void ExpandNormalSolution
( const SparseMatrix<Real>& A,
        Real gamma,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
  const Matrix<Real>& rc,
  const Matrix<Real>& rmu,
        Matrix<Real>& dx,
  const Matrix<Real>& dy, 
        Matrix<Real>& dz )
{
    DEBUG_ONLY(CSE cse("lp::direct::ExpandNormalSolution"))
    const Int n = A.Width();

    // dInv := sqrt( (z ./ x) .+ gamma^2 )
    // ===================================
    Matrix<Real> dInv;
    dInv.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        dInv.Set( i, 0, Sqrt(z.Get(i,0)/x.Get(i,0) + gamma*gamma) );

    // Use dz as a temporary for storing A^T dy
    // ========================================
    Zeros( dz, n, 1 );
    Multiply( TRANSPOSE, Real(1), A, dy, Real(0), dz );

    // dx = D^2 (-A^T dy + r_1 + inv(X) r_3) 
    //    = D^2 (-A^T dy - r_c - inv(X) r_mu)
    // ======================================
    dx = rmu;
    DiagonalSolve( LEFT, NORMAL, x, dx );
    dx += rc;
    dx += dz;
    dx *= -1;
    DiagonalSolve( LEFT, NORMAL, dInv, dx );
    DiagonalSolve( LEFT, NORMAL, dInv, dx );

    // dz := r_c + gamma^2 dx + A^T dy
    // ===============================
    Axpy( gamma*gamma, dx, dz );
    dz += rc;
}

template<typename Real>
void ExpandNormalSolution
( const DistSparseMatrix<Real>& A,
        Real gamma,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rmu,
        DistMultiVec<Real>& dx,
  const DistMultiVec<Real>& dy, 
        DistMultiVec<Real>& dz )
{
    DEBUG_ONLY(CSE cse("lp::direct::ExpandNormalSolution"))
    const Int n = A.Width();

    DistMultiVec<Real> dInv(A.Comm());
    dInv.Resize( n, 1 );
    {
        const Int nLocal = dInv.LocalHeight();
        for( Int iLoc=0; iLoc<nLocal; ++iLoc )
            dInv.SetLocal
            ( iLoc, 0, 
              Sqrt(z.GetLocal(iLoc,0)/x.GetLocal(iLoc,0)+gamma*gamma) );
    }

    // Use dz as a temporary for storing A^T dy
    // ========================================
    Zeros( dz, n, 1 );
    Multiply( TRANSPOSE, Real(1), A, dy, Real(0), dz );

    // dx = D^2 (-A^T dy + r_1 + inv(X) r_3) 
    //    = D^2 (-A^T dy - r_c - inv(X) r_mu)
    // ======================================
    dx = rmu;
    DiagonalSolve( LEFT, NORMAL, x, dx );
    dx += rc;
    dx += dz;
    dx *= -1;
    DiagonalSolve( LEFT, NORMAL, dInv, dx );
    DiagonalSolve( LEFT, NORMAL, dInv, dx );

    // dz := r_c + gamma^2 dx + A^T dy
    // ===============================
    Axpy( gamma*gamma, dx, dz );
    dz += rc;
}

#define PROTO(Real) \
  template void NormalKKT \
  ( const Matrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
          Matrix<Real>& J, bool onlyLower ); \
  template void NormalKKT \
  ( const ElementalMatrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Real>& z, \
          ElementalMatrix<Real>& J, bool onlyLower ); \
  template void NormalKKT \
  ( const SparseMatrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
          SparseMatrix<Real>& J, bool onlyLower ); \
  template void NormalKKT \
  ( const DistSparseMatrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J, bool onlyLower ); \
  template void NormalKKTRHS \
  ( const Matrix<Real>& A, \
          Real gamma, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
    const Matrix<Real>& rc, \
    const Matrix<Real>& rb, \
    const Matrix<Real>& rmu, \
          Matrix<Real>& d ); \
  template void NormalKKTRHS \
  ( const ElementalMatrix<Real>& A, \
          Real gamma, \
    const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Real>& z, \
    const ElementalMatrix<Real>& rc, \
    const ElementalMatrix<Real>& rb, \
    const ElementalMatrix<Real>& rmu, \
          ElementalMatrix<Real>& d ); \
  template void NormalKKTRHS \
  ( const SparseMatrix<Real>& A, \
          Real gamma, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
    const Matrix<Real>& rc, \
    const Matrix<Real>& rb, \
    const Matrix<Real>& rmu, \
          Matrix<Real>& d ); \
  template void NormalKKTRHS \
  ( const DistSparseMatrix<Real>& A, \
          Real gamma, \
    const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& rc, \
    const DistMultiVec<Real>& rb, \
    const DistMultiVec<Real>& rmu, \
          DistMultiVec<Real>& d ); \
  template void ExpandNormalSolution \
  ( const Matrix<Real>& A, \
          Real gamma, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
    const Matrix<Real>& rc, \
    const Matrix<Real>& rmu, \
          Matrix<Real>& dx, \
    const Matrix<Real>& dy, \
          Matrix<Real>& dz ); \
  template void ExpandNormalSolution \
  ( const ElementalMatrix<Real>& A, \
          Real gamma, \
    const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Real>& z, \
    const ElementalMatrix<Real>& rc, \
    const ElementalMatrix<Real>& rmu, \
          ElementalMatrix<Real>& dx, \
    const ElementalMatrix<Real>& dy, \
          ElementalMatrix<Real>& dz ); \
  template void ExpandNormalSolution \
  ( const SparseMatrix<Real>& A, \
          Real gamma, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
    const Matrix<Real>& rc, \
    const Matrix<Real>& rmu, \
          Matrix<Real>& dx, \
    const Matrix<Real>& dy, \
          Matrix<Real>& dz ); \
  template void ExpandNormalSolution \
  ( const DistSparseMatrix<Real>& A, \
          Real gamma, \
    const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& rc, \
    const DistMultiVec<Real>& rmu, \
          DistMultiVec<Real>& dx, \
    const DistMultiVec<Real>& dy, \
          DistMultiVec<Real>& dz );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace direct
} // namespace lp
} // namespace El
