/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "../util.hpp"

namespace El {
namespace lp {
namespace primal {

// Solve
//
//    | I  A^T | | x | = | -c |,
//    | A   0  | | y |   |  b |
//
// then compute 
//
//    alpha_p = -min x,
//    alpha_d = -max x,
//
// and set 
//
//    x := ( alpha_p < 0 ?  x :  x + (1+alpha_p)e ),
//    z := ( alpha_d < 0 ? -x : -x + (1+alpha_d)e ).
//
// Note that this is a simplification of the "dual" conic form initialization
// procedure described in 
//
//   L. Vandenberghe
//   "The CVXOPT linear and quadratic cone program solvers"
//   <http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf>
//

template<typename Real>
void Initialize
( const Matrix<Real>& A, 
  const Matrix<Real>& b, const Matrix<Real>& c,
        Matrix<Real>& x,       Matrix<Real>& y,
        Matrix<Real>& z )
{
    DEBUG_ONLY(CallStackEntry cse("lp::primal::Initialize"))
    const Int n = A.Width();

    Matrix<Real> J, ones;
    Ones( ones, n, 1 );
    AugmentedKKT( A, ones, ones, J );

    Matrix<Real> rc, rb, rmu, d;
    rc = c;
    rb = b;
    Scale( Real(-1), rb );
    Zeros( rmu, n, 1 );
    AugmentedKKTRHS( ones, rc, rb, rmu, d );

    SymmetricSolve( LOWER, NORMAL, J, d );

    ExpandAugmentedSolution( ones, ones, rmu, d, x, y, z );
    z = x;
    Scale( Real(-1), z );

    const auto xMinPair = VectorMin( x );
    const auto zMinPair = VectorMin( z );
    const Real alphaPrimal = -xMinPair.value;
    const Real alphaDual   = -zMinPair.value;
    if( alphaPrimal >= Real(0) )
        Shift( x, alphaPrimal+1 );
    if( alphaDual >= Real(0) )
        Shift( z, alphaDual+1 );
}

template<typename Real>
void Initialize
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z )
{
    DEBUG_ONLY(CallStackEntry cse("lp::primal::Initialize"))
    const Int n = A.Width();
    const Grid& g = A.Grid();

    DistMatrix<Real> J(g), ones(g);
    Ones( ones, n, 1 );
    AugmentedKKT( A, ones, ones, J );

    DistMatrix<Real> rc(g), rb(g), rmu(g), d(g);
    rc = c;
    rb = b;
    Scale( Real(-1), rb );
    Zeros( rmu, n, 1 );
    AugmentedKKTRHS( ones, rc, rb, rmu, d );

    SymmetricSolve( LOWER, NORMAL, J, d );

    ExpandAugmentedSolution( ones, ones, rmu, d, x, y, z );
    Copy( x, z );
    Scale( Real(-1), z );

    const auto xMinPair = VectorMin( x );
    const auto zMinPair = VectorMin( z );
    const Real alphaPrimal = -xMinPair.value;
    const Real alphaDual   = -zMinPair.value;
    if( alphaPrimal >= Real(0) )
        Shift( x, alphaPrimal+1 );
    if( alphaDual >= Real(0) )
        Shift( z, alphaDual+1 );
}

template<typename Real>
void Initialize
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,       const Matrix<Real>& c,
        Matrix<Real>& x,             Matrix<Real>& y,
        Matrix<Real>& z )
{
    DEBUG_ONLY(CallStackEntry cse("lp::primal::Initialize"))
    const Int n = A.Width();

    SparseMatrix<Real> J;
    Matrix<Real> ones;
    Ones( ones, n, 1 );
    AugmentedKKT( A, ones, ones, J, false );

    Matrix<Real> rc, rb, rmu, d;
    rc = c;
    rb = b;
    Scale( Real(-1), rb );
    Zeros( rmu, n, 1 );
    AugmentedKKTRHS( ones, rc, rb, rmu, d );

    // TODO: Solve with regularization
    //SymmetricSolve( LOWER, NORMAL, J, d );
    LogicError("Sequential sparse-direct solves not yet supported");

    ExpandAugmentedSolution( ones, ones, rmu, d, x, y, z );
    z = x;
    Scale( Real(-1), z );

    const auto xMinPair = VectorMin( x );
    const auto zMinPair = VectorMin( z );
    const Real alphaPrimal = -xMinPair.value;
    const Real alphaDual   = -zMinPair.value;
    if( alphaPrimal >= Real(0) )
        Shift( x, alphaPrimal+1 );
    if( alphaDual >= Real(0) )
        Shift( z, alphaDual+1 );
}

template<typename Real>
void Initialize
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b,      const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,            DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMap& map,                     DistMap& invMap, 
        DistSeparatorTree& sepTree,       DistSymmInfo& info )
{
    DEBUG_ONLY(CallStackEntry cse("lp::primal::Initialize"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();

    DistSparseMatrix<Real> J(comm);
    DistMultiVec<Real> ones(comm);
    Ones( ones, n, 1 );
    AugmentedKKT( A, ones, ones, J, false );

    DistMultiVec<Real> rc(comm), rb(comm), rmu(comm), d(comm);
    rc = c;
    rb = b;
    Scale( Real(-1), rb );
    Zeros( rmu, n, 1 );
    AugmentedKKTRHS( ones, rc, rb, rmu, d );

    // Approximately solve the system using dynamic regularization
    // ===========================================================
    DistMultiVec<Real> regCand(comm), reg(comm);
    DistNodalMultiVec<Real> regCandNodal, regNodal;
    const Real epsilon = lapack::MachineEpsilon<Real>();
    const Real pivTol = MaxNorm(J)*epsilon;
    const Real regMagPrimal = Pow(epsilon,Real(0.75));
    const Real regMagLagrange = Pow(epsilon,Real(0.5));
    regCand.Resize( n+m, 1 );
    for( Int iLoc=0; iLoc<regCand.LocalHeight(); ++iLoc )
    {
        const Int i = regCand.FirstLocalRow() + iLoc;
        if( i < n )
            regCand.SetLocal( iLoc, 0, regMagPrimal );
        else
            regCand.SetLocal( iLoc, 0, -regMagLagrange );
    }
    // Do not use any a priori regularization
    Zeros( reg, m+n, 1 );
    // Compute the proposed step from the KKT system
    // ---------------------------------------------
    NestedDissection( J.LockedDistGraph(), map, sepTree, info );
    map.FormInverse( invMap );

    DistSymmFrontTree<Real> JFrontTree;
    JFrontTree.Initialize( J, map, sepTree, info );
    regCandNodal.Pull( invMap, info, regCand );
    regNodal.Pull( invMap, info, reg );
    RegularizedLDL
    ( info, JFrontTree, pivTol, regCandNodal, regNodal, LDL_1D );
    regNodal.Push( invMap, info, reg );
    // TODO: Iterative refinement
    /*
    SolveWithIterativeRefinement
    ( J, invMap, info, JFrontTree, d, 
      minReductionFactor, maxRefineIts );
    */
    DistNodalMultiVec<Real> dNodal;
    dNodal.Pull( invMap, info, d );
    Solve( info, JFrontTree, dNodal );
    dNodal.Push( invMap, info, d );

    ExpandAugmentedSolution( ones, ones, rmu, d, x, y, z );
    z = x;
    Scale( Real(-1), z );

    const auto xMinPair = VectorMin( x );
    const auto zMinPair = VectorMin( z );
    const Real alphaPrimal = -xMinPair.value;
    const Real alphaDual   = -zMinPair.value;
    if( alphaPrimal >= Real(0) )
        Shift( x, alphaPrimal+1 );
    if( alphaDual >= Real(0) )
        Shift( z, alphaDual+1 );
}

#define PROTO(Real) \
  template void Initialize \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
          Matrix<Real>& x,       Matrix<Real>& y, \
          Matrix<Real>& z ); \
  template void Initialize \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
          AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z ); \
  template void Initialize \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b,       const Matrix<Real>& c, \
          Matrix<Real>& x,             Matrix<Real>& y, \
          Matrix<Real>& z ); \
  template void Initialize \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b,      const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x,            DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
          DistMap& map,                     DistMap& invMap, \
          DistSeparatorTree& sepTree,       DistSymmInfo& info );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace primal
} // namespace lp
} // namespace El
