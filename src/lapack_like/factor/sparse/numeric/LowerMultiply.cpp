/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./LowerMultiply/Front.hpp"

#include "./LowerMultiply/Local.hpp"
#include "./LowerMultiply/Dist.hpp"

namespace El {

template<typename T>
void LowerMultiply
( Orientation orientation, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, 
  DistNodalMultiVec<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LowerMultiply"))
    if( orientation == NORMAL )
    {
        LocalLowerMultiplyNormal( diagOffset, info, L, X );
        DistLowerMultiplyNormal( diagOffset, info, L, X );
    }
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        DistLowerMultiplyTranspose( diagOffset, info, L, X, conjugate );
        LocalLowerMultiplyTranspose( diagOffset, info, L, X, conjugate );
    }
}

#define PROTO(T) \
  template void LowerMultiply \
  ( Orientation orientation, int diagOffset, \
    const DistSymmInfo& info, const DistSymmFrontTree<T>& L, \
    DistNodalMultiVec<T>& X );

// NOTE: Trmm is not currently implemented for Int
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
