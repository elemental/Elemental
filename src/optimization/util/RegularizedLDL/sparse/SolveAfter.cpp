/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace reg_ldl {

template<typename F>
Int SolveAfter
( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,             const DistSymmInfo& info,
  const DistSymmFrontTree<F>& AFact,       DistMultiVec<F>& y,
  Base<F> minReductionFactor,              Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_ldl::SolveAfter"))
    mpi::Comm comm = y.Comm();
    const Int commRank = mpi::Rank(comm);

    DistMultiVec<F> yOrig(comm);
    yOrig = y;

    // Compute the initial guess
    // =========================
    DistMultiVec<F> x(comm);
    DistNodalMultiVec<F> xNodal;
    xNodal.Pull( invMap, info, y );
    Solve( info, AFact, xNodal );
    xNodal.Push( invMap, info, x );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<F> dx(comm), xCand(comm), xRegScaled(comm);
        xRegScaled = x;
        DiagonalScale( NORMAL, reg, xRegScaled );
        Axpy( F(1), xRegScaled, y );
        Multiply( NORMAL, F(-1), A, x, F(1), y );
        Base<F> errorNorm = Nrm2( y );
        if( progress && commRank == 0 )
            std::cout << "    original error norm: " << errorNorm << std::endl;
        for( ; refineIt<maxRefineIts; ++refineIt )
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
            xRegScaled = xCand;
            DiagonalScale( NORMAL, reg, xRegScaled );
            Axpy( F(-1), xRegScaled, y );
            Multiply( NORMAL, F(-1), A, xCand, F(1), y );
            Base<F> newErrorNorm = Nrm2( y );
            if( progress && commRank == 0 )
                std::cout << "    reduced by factor " 
                          << errorNorm/newErrorNorm << std::endl;
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
  template Int SolveAfter \
  ( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg, \
    const DistMap& invMap,             const DistSymmInfo& info, \
    const DistSymmFrontTree<F>& AFact,       DistMultiVec<F>& y, \
    Base<F> minReductionFactor,              Int maxRefineIts, \
    bool progress );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace reg_ldl
} // namespace El
