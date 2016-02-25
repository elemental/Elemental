/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./LowerSolve/Forward.hpp"
#include "./LowerSolve/Backward.hpp"

namespace El {
namespace ldl {

template<typename F>
void LowerSolve
( Orientation orientation,
  const NodeInfo& info, 
  const Front<F>& front,
        MatrixNode<F>& X )
{
    DEBUG_ONLY(CSE cse("LowerSolve"))
    if( orientation == NORMAL )
        LowerForwardSolve( info, front, X );
    else
        LowerBackwardSolve( info, front, X, orientation==ADJOINT );
}

template<typename F>
void LowerSolve
( Orientation orientation,
  const DistNodeInfo& info, 
  const DistFront<F>& front,
        DistMultiVecNode<F>& X )
{
    DEBUG_ONLY(CSE cse("LowerSolve"))
    if( orientation == NORMAL )
        LowerForwardSolve( info, front, X );
    else
        LowerBackwardSolve( info, front, X, orientation==ADJOINT );
}

template<typename F>
void LowerSolve
( Orientation orientation,
  const DistNodeInfo& info, 
  const DistFront<F>& front,
        DistMatrixNode<F>& X )
{
    DEBUG_ONLY(CSE cse("LowerSolve"))
    if( orientation == NORMAL )
        LowerForwardSolve( info, front, X );
    else
        LowerBackwardSolve( info, front, X, orientation==ADJOINT );
}

#define PROTO(F) \
  template void LowerSolve \
  ( Orientation orientation, \
    const NodeInfo& info, \
    const Front<F>& front, \
          MatrixNode<F>& X ); \
  template void LowerSolve \
  ( Orientation orientation, \
    const DistNodeInfo& info, \
    const DistFront<F>& front, \
          DistMultiVecNode<F>& X ); \
  template void LowerSolve \
  ( Orientation orientation, \
    const DistNodeInfo& info, \
    const DistFront<F>& front, \
          DistMatrixNode<F>& X );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace ldl
} // namespace El
