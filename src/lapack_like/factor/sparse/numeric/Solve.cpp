/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

namespace El {

template<typename F>
void Solve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("Solve"))
    if( !FrontsAre1d(L.frontType) )
        LogicError("Invalid front type for 1D solve");
    const Orientation orientation = ( L.isHermitian ? ADJOINT : TRANSPOSE );
    if( BlockFactorization(L.frontType) )
    {
        // Solve against block diagonal factor, L D
        LowerSolve( NORMAL, info, L, X );
        // Solve against the (conjugate-)transpose of the block unit diagonal L
        LowerSolve( orientation, info, L, X );
    }
    else
    {
        // Solve against unit diagonal L
        LowerSolve( NORMAL, info, L, X );
        // Solve against diagonal
        DiagonalSolve( info, L, X );
        // Solve against the (conjugate-)transpose of the unit diagonal L
        LowerSolve( orientation, info, L, X );
    }
}

template<typename F>
void Solve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("Solve"))
    if( FrontsAre1d(L.frontType) )
        LogicError("Invalid front type for 2D solve");
    const Orientation orientation = ( L.isHermitian ? ADJOINT : TRANSPOSE );
    if( BlockFactorization(L.frontType) )
    {
        // Solve against block diagonal factor, L D
        LowerSolve( NORMAL, info, L, X );
        // Solve against the (conjugate-)transpose of the block unit diagonal L
        LowerSolve( orientation, info, L, X );
    }
    else
    {
        // Solve against unit diagonal L
        LowerSolve( NORMAL, info, L, X );
        // Solve against diagonal
        DiagonalSolve( info, L, X );
        // Solve against the (conjugate-)transpose of the unit diagonal L
        LowerSolve( orientation, info, L, X );
    }
} 

template<typename F>
void SolveWithIterativeRefinement
( const DistSparseMatrix<F>& A,
  const DistMap& invMap, const DistSymmInfo& info,
  const DistSymmFrontTree<F>& AFact, DistMultiVec<F>& y,
  Base<F> minReductionFactor, Int maxRefineIts )
{
    DEBUG_ONLY(CallStackEntry cse("IterativeRefinement"))
    mpi::Comm comm = y.Comm();

    DistMultiVec<F> yOrig(comm);
    yOrig = y;

    // Compute the initial guess
    // =========================
    DistMultiVec<F> x(comm);
    DistNodalMultiVec<F> xNodal;
    xNodal.Pull( invMap, info, y );
    Solve( info, AFact, xNodal );
    xNodal.Push( invMap, info, x );

    if( maxRefineIts > 0 )
    {
        DistMultiVec<F> dx(comm), xCand(comm); 
        Multiply( NORMAL, F(-1), A, x, F(1), y );
        Base<F> errorNorm = Nrm2( y );
        for( Int refineIt=0; refineIt<maxRefineIts; ++refineIt )
        {
            // Compute the proposed update to the solution
            // -------------------------------------------
            xNodal.Pull( invMap, info, y );
            Solve( info, AFact, xNodal );
            xNodal.Push( invMap, info, dx );
            xCand = x;
            Axpy( F(1), dx, xCand );

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
            if( newErrorNorm < errorNorm )
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
}

// TODO: Add iterative refinement parameter
template<typename F>
void SymmetricSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, 
  bool conjugate, const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricSolve"))
    DistSymmInfo info;
    DistSeparatorTree sepTree;
    DistMap map, inverseMap;
    NestedDissection( A.LockedDistGraph(), map, sepTree, info, ctrl );
    map.FormInverse( inverseMap );

    DistSymmFrontTree<F> frontTree( A, map, sepTree, info, conjugate );
    LDL( info, frontTree, LDL_INTRAPIV_1D );

    DistNodalMultiVec<F> XNodal;
    XNodal.Pull( inverseMap, info, X );
    Solve( info, frontTree, XNodal );
    XNodal.Push( inverseMap, info, X );
}

// TODO: Add iterative refinement parameter
template<typename F>
void HermitianSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSolve"))
    SymmetricSolve( A, X, true, ctrl );
}

#define PROTO(F) \
  template void Solve \
  ( const DistSymmInfo& info, \
    const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X ); \
  template void Solve \
  ( const DistSymmInfo& info, \
    const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X ); \
  template void SolveWithIterativeRefinement \
  ( const DistSparseMatrix<F>& A, \
    const DistMap& invMap, const DistSymmInfo& info, \
    const DistSymmFrontTree<F>& AFact, DistMultiVec<F>& y, \
    Base<F> minReductionFactor, Int maxRefineIts ); \
  template void SymmetricSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, bool conjugate, \
    const BisectCtrl& ctrl ); \
  template void HermitianSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, \
    const BisectCtrl& ctrl );
 
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
