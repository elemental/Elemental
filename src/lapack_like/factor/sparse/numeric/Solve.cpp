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
  template void SymmetricSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, bool conjugate, \
    const BisectCtrl& ctrl ); \
  template void HermitianSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, \
    const BisectCtrl& ctrl );
 
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
