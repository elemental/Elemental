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

#include "./LDL/ProcessFront.hpp"
#include "./LDL/ProcessFrontBlock.hpp"

#include "./LDL/ProcessLocalTree.hpp"
#include "./LDL/ProcessDistTree.hpp"

namespace El {

template<typename T>
void InitializeDistLeaf
( const DistSymmInfo& info, DistSymmFrontTree<T>& L )
{
    DEBUG_ONLY(CallStackEntry cse("InitializeDistLeaf"))
    const DistSymmNodeInfo& node = info.distNodes[0];
    Matrix<T>& topLocalFrontL = L.localFronts.back().frontL;
    DistMatrix<T>& front2dL = L.distFronts[0].front2dL;

    front2dL.LockedAttach
    ( topLocalFrontL.Height(), topLocalFrontL.Width(), *node.grid, 0, 0,
      topLocalFrontL.LockedBuffer(), topLocalFrontL.LDim() );
}

template<typename F>
void 
LDL( DistSymmInfo& info, DistSymmFrontTree<F>& L, SymmFrontType newFrontType )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    if( !Unfactored(L.frontType) )
        LogicError("Matrix is already factored");

    // Convert from 1D to 2D if necessary
    ChangeFrontType( L, SYMM_2D );

    // Perform the initial factorization
    L.frontType = InitialFactorType(newFrontType);
    ldl::ProcessLocalTree( info, L );
    ldl::ProcessDistTree( info, L );

    // Convert the fronts from the initial factorization to the requested form
    ChangeFrontType( L, newFrontType );
}

#define PROTO_INT(T) \
  template void InitializeDistLeaf \
  ( const DistSymmInfo& info, DistSymmFrontTree<T>& L );

#define PROTO(F) \
  PROTO_INT(F) \
  template void LDL \
  ( DistSymmInfo& info, DistSymmFrontTree<F>& L, SymmFrontType newFrontType );

#include "El/macros/Instantiate.h"

} // namespace El
