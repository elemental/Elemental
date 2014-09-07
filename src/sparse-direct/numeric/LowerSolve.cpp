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

#include "./LowerSolve/Front.hpp"
#include "./LowerSolve/FrontBlock.hpp"

#include "./LowerSolve/Local.hpp"
#include "./LowerSolve/Dist.hpp"

namespace El {

template<typename F>
void LowerSolve
( Orientation orientation, const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LowerSolve"))
    if( orientation == NORMAL )
    {
        LocalLowerForwardSolve( info, L, X );
        DistLowerForwardSolve( info, L, X );
    }
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        DistLowerBackwardSolve( info, L, X, conjugate );
        LocalLowerBackwardSolve( info, L, X, conjugate );
    }
}

template<typename F>
void LowerSolve
( Orientation orientation, const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LowerSolve"))
    if( orientation == NORMAL )
    {
        LocalLowerForwardSolve( info, L, X );
        DistLowerForwardSolve( info, L, X );
    }
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        DistLowerBackwardSolve( info, L, X, conjugate );
        LocalLowerBackwardSolve( info, L, X, conjugate );
    }
}

#define PROTO(F) \
  template void LowerSolve \
  ( Orientation orientation, const DistSymmInfo& info, \
    const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X ); \
  template void LowerSolve \
  ( Orientation orientation, const DistSymmInfo& info, \
    const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
