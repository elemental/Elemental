/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./LowerSolve/Forward.hpp"
#include "./LowerSolve/Backward.hpp"

namespace El {

template<typename F>
void LowerSolve
( Orientation orientation, const SymmNodeInfo& info, 
  const SymmFront<F>& front, MatrixNode<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LowerSolve"))
    if( orientation == NORMAL )
        LowerForwardSolve( info, front, X );
    else
        LowerBackwardSolve( info, front, X, orientation==ADJOINT );
}

template<typename F>
void LowerSolve
( Orientation orientation, const DistSymmNodeInfo& info, 
  const DistSymmFront<F>& front, DistMultiVecNode<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LowerSolve"))
    if( orientation == NORMAL )
        LowerForwardSolve( info, front, X );
    else
        LowerBackwardSolve( info, front, X, orientation==ADJOINT );
}

template<typename F>
void LowerSolve
( Orientation orientation, const DistSymmNodeInfo& info, 
  const DistSymmFront<F>& front, DistMatrixNode<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LowerSolve"))
    if( orientation == NORMAL )
        LowerForwardSolve( info, front, X );
    else
        LowerBackwardSolve( info, front, X, orientation==ADJOINT );
}

#define PROTO(F) \
  template void LowerSolve \
  ( Orientation orientation, const SymmNodeInfo& info, \
    const SymmFront<F>& front, MatrixNode<F>& X ); \
  template void LowerSolve \
  ( Orientation orientation, const DistSymmNodeInfo& info, \
    const DistSymmFront<F>& front, DistMultiVecNode<F>& X ); \
  template void LowerSolve \
  ( Orientation orientation, const DistSymmNodeInfo& info, \
    const DistSymmFront<F>& front, DistMatrixNode<F>& X );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
