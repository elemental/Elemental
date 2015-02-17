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
( const SparseMatrix<F>& A, Matrix<F>& X, 
  bool conjugate, const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricSolve"))
    SymmNodeInfo info;
    Separator rootSep;
    vector<Int> map, invMap;
    NestedDissection( A.LockedGraph(), map, rootSep, info, ctrl );
    InvertMap( map, invMap );

    SymmFront<F> front( A, map, info, conjugate );
    LDL( info, front );

    // TODO: Extend ldl::SolveWithIterativeRefinement to support multiple
    //       right-hand sides
    /*
    ldl::SolveWithIterativeRefinement
    ( A, invMap, info, front, X, minReductionFactor, maxRefineIts );
    */
    ldl::SolveAfter( invMap, info, front, X );
}

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
    InvertMap( map, invMap );

    DistSymmFront<F> front( A, map, rootSep, info, conjugate );
    LDL( info, front, LDL_INTRAPIV_1D );

    // TODO: Extend ldl::SolveWithIterativeRefinement to support multiple
    //       right-hand sides
    /*
    ldl::SolveWithIterativeRefinement
    ( A, invMap, info, front, X, minReductionFactor, maxRefineIts );
    */
    ldl::SolveAfter( invMap, info, front, X );
}

// TODO: Add iterative refinement parameter
template<typename F>
void HermitianSolve
( const SparseMatrix<F>& A, Matrix<F>& X, const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSolve"))
    SymmetricSolve( A, X, true, ctrl );
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
  ( const SparseMatrix<F>& A, Matrix<F>& X, bool conjugate, \
    const BisectCtrl& ctrl ); \
  template void SymmetricSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, bool conjugate, \
    const BisectCtrl& ctrl ); \
  template void HermitianSolve \
  ( const SparseMatrix<F>& A, Matrix<F>& X, \
    const BisectCtrl& ctrl ); \
  template void HermitianSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, \
    const BisectCtrl& ctrl );
 
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
