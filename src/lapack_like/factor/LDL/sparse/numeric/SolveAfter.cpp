/*
   Copyright (c) 2013-2015, Jack Poulson.
   All rights reserved.

   Copyright (c) 2011-2013, Jack Poulson and Lexing Ying.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

namespace El {
namespace ldl {

template<typename F>
void SolveAfter
( const vector<Int>& invMap,
  const NodeInfo& info, 
  const Front<F>& front,
        Matrix<F>& X )
{
    DEBUG_ONLY(CSE cse("ldl::SolveAfter"))

    MatrixNode<F> XNodal( invMap, info, X );
    SolveAfter( info, front, XNodal );
    XNodal.Push( invMap, info, X );
}

template<typename F>
void SolveAfter
( const NodeInfo& info,
  const Front<F>& front,
        MatrixNode<F>& X )
{
    DEBUG_ONLY(CSE cse("ldl::SolveAfter"))

    const Orientation orientation = ( front.isHermitian ? ADJOINT : TRANSPOSE );
    if( BlockFactorization(front.type) )
    {
        // Solve against block diagonal factor, L D
        LowerSolve( NORMAL, info, front, X );
        // Solve against the (conjugate-)transpose of the block unit diagonal L
        LowerSolve( orientation, info, front, X );
    }
    else
    {
        // Solve against unit diagonal L
        LowerSolve( NORMAL, info, front, X );
        // Solve against diagonal
        DiagonalSolve( info, front, X );
        // Solve against the (conjugate-)transpose of the unit diagonal L
        LowerSolve( orientation, info, front, X );
    }
}

template<typename F>
void SolveAfter
( const DistMap& invMap,
  const DistNodeInfo& info, 
  const DistFront<F>& front,
        DistMultiVec<F>& X )
{
    DEBUG_ONLY(CSE cse("ldl::SolveAfter"))

    if( FrontIs1D(front.type) )
    {
        DistMultiVecNode<F> XNodal( invMap, info, X );
        SolveAfter( info, front, XNodal );
        XNodal.Push( invMap, info, X );
    }
    else
    {
        DistMatrixNode<F> XNodal( invMap, info, X );
        SolveAfter( info, front, XNodal );
        XNodal.Push( invMap, info, X );
    }
}

template<typename F>
void SolveAfter
( const DistNodeInfo& info, 
  const DistFront<F>& front,
        DistMultiVecNode<F>& X )
{
    DEBUG_ONLY(CSE cse("ldl::SolveAfter"))

    // TODO: Only perform the switch if there are a sufficient 
    //       number of right-hand sides?
    /*
    if( !FrontIs1D(front.type) )
    {
        // TODO: Add warning?
        DistMatrixNode<F> XMat( X );
        SolveAfter( info, front, XMat );
        X = XMat;
        return;
    }
    */

    const Orientation orientation = ( front.isHermitian ? ADJOINT : TRANSPOSE );
    if( BlockFactorization(front.type) )
    {
        // Solve against block diagonal factor, L D
        LowerSolve( NORMAL, info, front, X );
        // Solve against the (conjugate-)transpose of the block unit diagonal L
        LowerSolve( orientation, info, front, X );
    }
    else
    {
        // Solve against unit diagonal L
        LowerSolve( NORMAL, info, front, X );
        // Solve against diagonal
        DiagonalSolve( info, front, X );
        // Solve against the (conjugate-)transpose of the unit diagonal L
        LowerSolve( orientation, info, front, X );
    }
}

template<typename F>
void SolveAfter
( const DistNodeInfo& info, 
  const DistFront<F>& front,
        DistMatrixNode<F>& X )
{
    DEBUG_ONLY(CSE cse("ldl::SolveAfter"))

    if( FrontIs1D(front.type) )
    {
        // TODO: Add warning?
        DistMultiVecNode<F> XMV( X );
        SolveAfter( info, front, XMV );
        X = XMV;
        return;
    }

    const Orientation orientation = ( front.isHermitian ? ADJOINT : TRANSPOSE );
    if( BlockFactorization(front.type) )
    {
        // Solve against block diagonal factor, L D
        LowerSolve( NORMAL, info, front, X );
        // Solve against the (conjugate-)transpose of the block unit diagonal L
        LowerSolve( orientation, info, front, X );
    }
    else
    {
        // Solve against unit diagonal L
        LowerSolve( NORMAL, info, front, X );
        // Solve against diagonal
        DiagonalSolve( info, front, X );
        // Solve against the (conjugate-)transpose of the unit diagonal L
        LowerSolve( orientation, info, front, X );
    }
} 

// TODO: Improve these implementations 
//       (e.g., limit maxRefineIts to 3 by default)
template<typename F>
Int SolveWithIterativeRefinement
( const SparseMatrix<F>& A,
  const vector<Int>& invMap,
  const NodeInfo& info,
  const Front<F>& front,
        Matrix<F>& y,
  Base<F> minReductionFactor, Int maxRefineIts )
{
    DEBUG_ONLY(CSE cse("ldl::SolveWithIterativeRefinement"))
    auto yOrig = y;

    // Compute the initial guess
    // =========================
    Matrix<F> x;
    MatrixNode<F> xNodal( invMap, info, y );
    SolveAfter( info, front, xNodal );
    xNodal.Push( invMap, info, x );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        Matrix<F> dx, xCand; 
        Multiply( NORMAL, F(-1), A, x, F(1), y );
        Base<F> errorNorm = Nrm2( y );
        for( ; refineIt<maxRefineIts; ++refineIt )
        {
            // Compute the proposed update to the solution
            // -------------------------------------------
            xNodal.Pull( invMap, info, y );
            SolveAfter( info, front, xNodal );
            xNodal.Push( invMap, info, dx );
            xCand = x;
            xCand += dx;

            // If the proposed update lowers the residual, accept it
            // -----------------------------------------------------
            y = yOrig;
            Multiply( NORMAL, F(-1), A, xCand, F(1), y );
            Base<F> newErrorNorm = Nrm2( y );
            if( minReductionFactor*newErrorNorm < errorNorm )
            {
                x = xCand;
                errorNorm = newErrorNorm;
            }
            else if( newErrorNorm < errorNorm )
            {
                x = xCand;
                errorNorm = newErrorNorm;
                break;
            }
            else
                break;
        }
    }
    // Store the final result
    // ======================
    y = x;
    return refineIt;
}

template<typename F>
Int SolveWithIterativeRefinement
( const DistSparseMatrix<F>& A,
  const DistMap& invMap,
  const DistNodeInfo& info,
  const DistFront<F>& front,
        DistMultiVec<F>& y,
  Base<F> minReductionFactor, Int maxRefineIts )
{
    DEBUG_ONLY(CSE cse("ldl::SolveWithIterativeRefinement"))
    mpi::Comm comm = y.Comm();

    DistMultiVec<F> yOrig(comm);
    yOrig = y;

    ldl::DistMultiVecNodeMeta meta;

    // Compute the initial guess
    // =========================
    DistMultiVec<F> x(comm);
    DistMultiVecNode<F> xNodal;
    xNodal.Pull( invMap, info, y, meta );
    SolveAfter( info, front, xNodal );
    xNodal.Push( invMap, info, x, meta );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<F> dx(comm), xCand(comm); 
        Multiply( NORMAL, F(-1), A, x, F(1), y );
        Base<F> errorNorm = Nrm2( y );
        for( ; refineIt<maxRefineIts; ++refineIt )
        {
            // Compute the proposed update to the solution
            // -------------------------------------------
            xNodal.Pull( invMap, info, y, meta );
            SolveAfter( info, front, xNodal );
            xNodal.Push( invMap, info, dx, meta );
            xCand = x;
            xCand += dx;

            // If the proposed update lowers the residual, accept it
            // -----------------------------------------------------
            y = yOrig;
            Multiply( NORMAL, F(-1), A, xCand, F(1), y );
            Base<F> newErrorNorm = Nrm2( y );
            if( minReductionFactor*newErrorNorm < errorNorm )
            {
                x = xCand;
                errorNorm = newErrorNorm;
            }
            else if( newErrorNorm < errorNorm )
            {
                x = xCand;
                errorNorm = newErrorNorm;
                break;
            }
            else
                break;
        }
    }
    // Store the final result
    // ======================
    y = x;
    return refineIt;
}

#define PROTO(F) \
  template void SolveAfter \
  ( const vector<Int>& invMap, \
    const NodeInfo& info, \
    const Front<F>& front, \
          Matrix<F>& X ); \
  template void SolveAfter \
  ( const DistMap& invMap, \
    const DistNodeInfo& info, \
    const DistFront<F>& front, \
          DistMultiVec<F>& X ); \
  template void SolveAfter \
  ( const NodeInfo& info, \
    const Front<F>& front, \
          MatrixNode<F>& X ); \
  template void SolveAfter \
  ( const DistNodeInfo& info, \
    const DistFront<F>& front, \
          DistMultiVecNode<F>& X ); \
  template void SolveAfter \
  ( const DistNodeInfo& info, \
    const DistFront<F>& front, \
          DistMatrixNode<F>& X ); \
  template Int SolveWithIterativeRefinement \
  ( const SparseMatrix<F>& A, \
    const vector<Int>& invMap, \
    const NodeInfo& info, \
    const Front<F>& front, \
          Matrix<F>& y, \
    Base<F> minReductionFactor, Int maxRefineIts ); \
  template Int SolveWithIterativeRefinement \
  ( const DistSparseMatrix<F>& A, \
    const DistMap& invMap, \
    const DistNodeInfo& info, \
    const DistFront<F>& front, \
          DistMultiVec<F>& y, \
    Base<F> minReductionFactor, Int maxRefineIts );
 
#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace ldl
} // namespace El
