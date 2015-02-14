/*
   Copyright (c) 2009-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2013-2014, Jack Poulson and 
   The Georgia Institute of Technology.
   All rights reserved.

   Copyright (c) 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

namespace El {

// TODO: Add iterative refinement parameter
template<typename F>
void SymmetricSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, 
  bool conjugate, const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricSolve"))
    DistSymmNodeInfo info;
    DistSeparator rootSep;
    DistMap map, invMap;
    NestedDissection( A.LockedDistGraph(), map, rootSep, info, ctrl );
    map.FormInverse( invMap );

    DistSymmFront<F> front( A, map, rootSep, info, conjugate );
    LDL( info, front, LDL_INTRAPIV_1D );

    // TODO: Extend ldl::SolveWithIterativeRefinement to support multiple
    //       right-hand sides
    /*
    ldl::SolveWithIterativeRefinement
    ( A, invMap, info, front, X, minReductionFactor, maxRefineIts );
    */
    DistMultiVecNode<F> XNodal( invMap, info, X );
    ldl::SolveAfter( info, front, XNodal );
    XNodal.Push( invMap, info, X );
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
  template void SymmetricSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, bool conjugate, \
    const BisectCtrl& ctrl ); \
  template void HermitianSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, \
    const BisectCtrl& ctrl );
 
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
