/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "./LP/direct/IPM.hpp"
#include "./LP/affine/IPM.hpp"
#include "./LP/MPS.hpp"

namespace El {

template<typename Real>
void LP
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == LP_ADMM )
        lp::direct::ADMM
        ( problem.A, problem.b, problem.c, solution.x, ctrl.admmCtrl );
    else if( ctrl.approach == LP_MEHROTRA )
        lp::direct::Mehrotra( problem, solution, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

// This interface is now deprecated.
template<typename Real>
void LP
( const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    DirectLPProblem<Matrix<Real>,Matrix<Real>> problem;
    DirectLPSolution<Matrix<Real>> solution;
    LockedView( problem.c, c );
    LockedView( problem.A, A );
    LockedView( problem.b, b );
    solution.x = x;
    solution.y = y;
    solution.z = z;
    LP( problem, solution, ctrl );
    x = solution.x;
    y = solution.y;
    z = solution.z;
}

template<typename Real>
void LP
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == LP_MEHROTRA )
        lp::affine::Mehrotra( problem, solution, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

// This interface is now deprecated.
template<typename Real>
void LP
( const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    AffineLPProblem<Matrix<Real>,Matrix<Real>> problem;
    LockedView( problem.c, c );
    LockedView( problem.A, A );
    LockedView( problem.b, b );
    LockedView( problem.G, G );
    LockedView( problem.h, h );
    AffineLPSolution<Matrix<Real>> solution;
    solution.x = x;
    solution.s = s;
    solution.y = y;
    solution.z = z;
    LP( problem, solution, ctrl );
    x = solution.x;
    s = solution.s;
    y = solution.y;
    z = solution.z;
}

template<typename Real>
void LP
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        DirectLPSolution<DistMatrix<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == LP_ADMM )
        lp::direct::ADMM
        ( problem.A, problem.b, problem.c, solution.x, ctrl.admmCtrl );
    else if( ctrl.approach == LP_MEHROTRA )
        lp::direct::Mehrotra( problem, solution, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

// This interface is now deprecated.
template<typename Real>
void LP
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>> problem;
    DirectLPSolution<DistMatrix<Real>> solution;
    ForceSimpleAlignments( problem, A.Grid() );
    ForceSimpleAlignments( solution, A.Grid() );
    Copy( c, problem.c );
    Copy( A, problem.A );
    Copy( b, problem.b );
    Copy( x, solution.x );
    Copy( y, solution.y );
    Copy( z, solution.z );
    LP( problem, solution, ctrl );
    Copy( solution.x, x );
    Copy( solution.y, y );
    Copy( solution.z, z );
}

template<typename Real>
void LP
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        AffineLPSolution<DistMatrix<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == LP_MEHROTRA )
        lp::affine::Mehrotra( problem, solution, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

// This interface is now deprecated.
template<typename Real>
void LP
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>> problem;
    AffineLPSolution<DistMatrix<Real>> solution;
    ForceSimpleAlignments( problem, A.Grid() );
    ForceSimpleAlignments( solution, A.Grid() );
    Copy( c, problem.c );
    Copy( A, problem.A );
    Copy( b, problem.b );
    Copy( G, problem.G );
    Copy( h, problem.h );
    Copy( x, solution.x );
    Copy( y, solution.y );
    Copy( z, solution.z );
    Copy( s, solution.s );
    LP( problem, solution, ctrl );
    Copy( solution.x, x );
    Copy( solution.y, y );
    Copy( solution.z, z );
    Copy( solution.s, s );
}

template<typename Real>
void LP
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == LP_MEHROTRA )
        lp::direct::Mehrotra( problem, solution, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

// This interface is now deprecated.
template<typename Real>
void LP
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    DirectLPProblem<SparseMatrix<Real>,Matrix<Real>> problem;
    DirectLPSolution<Matrix<Real>> solution;
    problem.c = c;
    problem.A = A;
    problem.b = b;
    solution.x = x;
    solution.y = y;
    solution.z = z;
    LP( problem, solution, ctrl );
    x = solution.x;
    y = solution.y;
    z = solution.z;
}

template<typename Real>
void LP
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == LP_MEHROTRA )
        lp::affine::Mehrotra( problem, solution, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

// This interface is now deprecated.
template<typename Real>
void LP
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    AffineLPProblem<SparseMatrix<Real>,Matrix<Real>> problem;
    AffineLPSolution<Matrix<Real>> solution;
    problem.c = c;
    problem.A = A;
    problem.b = b;
    problem.G = G;
    problem.h = h;
    solution.x = x;
    solution.y = y;
    solution.z = z;
    solution.s = s;
    LP( problem, solution, ctrl );
    x = solution.x;
    y = solution.y;
    z = solution.z;
    s = solution.s;
}

template<typename Real>
void LP
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        DirectLPSolution<DistMultiVec<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == LP_MEHROTRA )
        lp::direct::Mehrotra( problem, solution, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

// This interface is now deprecated.
template<typename Real>
void LP
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = A.Grid();

    DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>> problem;
    DirectLPSolution<DistMultiVec<Real>> solution;
    ForceSimpleAlignments( problem, grid );
    ForceSimpleAlignments( solution, grid );

    problem.c = c;
    problem.A = A;
    problem.b = b;
    solution.x = x;
    solution.y = y;
    solution.z = z;
    LP( problem, solution, ctrl );
    x = solution.x;
    y = solution.y;
    z = solution.z;
}

template<typename Real>
void LP
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        AffineLPSolution<DistMultiVec<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == LP_MEHROTRA )
        lp::affine::Mehrotra( problem, solution, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

// This interface is now deprecated.
template<typename Real>
void LP
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = A.Grid();

    AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>> problem;
    AffineLPSolution<DistMultiVec<Real>> solution;
    ForceSimpleAlignments( problem, grid );
    ForceSimpleAlignments( solution, grid );

    problem.c = c;
    problem.A = A;
    problem.b = b;
    problem.G = G;
    problem.h = h;
    solution.x = x;
    solution.y = y;
    solution.z = z;
    solution.s = s;
    LP( problem, solution, ctrl );
    x = solution.x;
    y = solution.y;
    z = solution.z;
    s = solution.s;
}

#define PROTO(Real) \
  template void LP \
  ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem, \
          DirectLPSolution<Matrix<Real>>& solution, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem, \
          AffineLPSolution<Matrix<Real>>& solution, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem, \
          DirectLPSolution<DistMatrix<Real>>& solution, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, \
    const AbstractDistMatrix<Real>& c, \
          AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem, \
          AffineLPSolution<DistMatrix<Real>>& solution, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& b, \
    const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& h, \
          AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& s, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem, \
          DirectLPSolution<Matrix<Real>>& solution, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem, \
          AffineLPSolution<Matrix<Real>>& solution, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem, \
          DirectLPSolution<DistMultiVec<Real>>& solution, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem, \
          AffineLPSolution<DistMultiVec<Real>>& solution, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
          DistMultiVec<Real>& s, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void ReadMPS \
  ( AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem, \
    const string& filename, \
    bool compressed, \
    bool minimize, \
    bool keepNonnegativeWithZeroUpperBounds, \
    bool metadataSummary ); \
  template void ReadMPS \
  ( AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem, \
    const string& filename, \
    bool compressed, \
    bool minimize, \
    bool keepNonnegativeWithZeroUpperBounds, \
    bool metadataSummary ); \
  template void ReadMPS \
  ( AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem, \
    const string& filename, \
    bool compressed, \
    bool minimize, \
    bool keepNonnegativeWithZeroUpperBounds, \
    bool metadataSummary ); \
  template void ReadMPS \
  ( AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem, \
    const string& filename, \
    bool compressed, \
    bool minimize, \
    bool keepNonnegativeWithZeroUpperBounds, \
    bool metadataSummary ); \
  template void WriteMPS \
  ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem, \
    const string& filename, \
    bool compressed ); \
  template void WriteMPS \
  ( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem, \
    const string& filename, \
    bool compressed ); \
  template void WriteMPS \
  ( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem, \
    const string& filename, \
    bool compressed ); \
  template void WriteMPS \
  ( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem, \
    const string& filename, \
    bool compressed ); \
  template void WriteMPS \
  ( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem, \
    const string& filename, \
    bool compressed ); \
  template void WriteMPS \
  ( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem, \
    const string& filename, \
    bool compressed ); \
  template void WriteMPS \
  ( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem, \
    const string& filename, \
    bool compressed ); \
  template void WriteMPS \
  ( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem, \
    const string& filename, \
    bool compressed );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
