/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./RegularizedQSDLDL/ProcessFront.hpp"
#include "./RegularizedQSDLDL/ProcessLocalTree.hpp"
#include "./RegularizedQSDLDL/ProcessDistTree.hpp"

namespace El {

template<typename F>
void RegularizedQSDLDL
( DistSymmInfo& info, DistSymmFrontTree<F>& L, Base<F> pivTol, 
  const DistNodalMultiVec<Base<F>>& regCand, 
        DistNodalMultiVec<Base<F>>& reg,
  bool aPriori,
  SymmFrontType newFrontType )
{
    DEBUG_ONLY(CallStackEntry cse("RegularizedQSDLDL"))
    if( !Unfactored(L.frontType) )
        LogicError("Matrix is already factored");

    // Convert from 1D to 2D if necessary
    ChangeFrontType( L, SYMM_2D );

    // Perform the initial factorization
    L.frontType = InitialFactorType(newFrontType);
    reg_qsd_ldl::ProcessLocalTree( info, L, pivTol, regCand, reg, aPriori );
    reg_qsd_ldl::ProcessDistTree( info, L, pivTol, regCand, reg, aPriori );

    // Convert the fronts from the initial factorization to the requested form
    ChangeFrontType( L, newFrontType );
}

#define PROTO(F) \
  template void RegularizedQSDLDL \
  ( DistSymmInfo& info, DistSymmFrontTree<F>& L, \
    Base<F> pivTol, \
    const DistNodalMultiVec<Base<F>>& regCand, \
          DistNodalMultiVec<Base<F>>& reg, \
    bool aPriori, \
    SymmFrontType newFrontType );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
