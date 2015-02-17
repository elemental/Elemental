/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./RegularizedQSDLDL/Process.hpp"

namespace El {

template<typename F>
void RegularizedQSDLDL
( const SymmNodeInfo& info, SymmFront<F>& front, Base<F> pivTol, 
  const MatrixNode<Base<F>>& regCand, 
        MatrixNode<Base<F>>& reg,
  bool aPriori,
  SymmFrontType newType )
{
    DEBUG_ONLY(CallStackEntry cse("RegularizedQSDLDL"))
    if( !Unfactored(front.type) )
        LogicError("Matrix is already factored");

    // Convert from 1D to 2D if necessary
    ChangeFrontType( front, SYMM_2D );

    // Perform the initial factorization
    reg_qsd_ldl::Process
    ( info, front, pivTol, regCand, reg, aPriori, InitialFactorType(newType) );

    // Convert the fronts from the initial factorization to the requested form
    ChangeFrontType( front, newType );
}

template<typename F>
void RegularizedQSDLDL
( const DistSymmNodeInfo& info, DistSymmFront<F>& front, Base<F> pivTol, 
  const DistMultiVecNode<Base<F>>& regCand, 
        DistMultiVecNode<Base<F>>& reg,
  bool aPriori,
  SymmFrontType newType )
{
    DEBUG_ONLY(CallStackEntry cse("RegularizedQSDLDL"))
    if( !Unfactored(front.type) )
        LogicError("Matrix is already factored");

    // Convert from 1D to 2D if necessary
    ChangeFrontType( front, SYMM_2D );

    // Perform the initial factorization
    reg_qsd_ldl::Process
    ( info, front, pivTol, regCand, reg, aPriori, InitialFactorType(newType) );

    // Convert the fronts from the initial factorization to the requested form
    ChangeFrontType( front, newType );
}

#define PROTO(F) \
  template void RegularizedQSDLDL \
  ( const SymmNodeInfo& info, SymmFront<F>& front, \
    Base<F> pivTol, \
    const MatrixNode<Base<F>>& regCand, \
          MatrixNode<Base<F>>& reg, \
    bool aPriori, \
    SymmFrontType newType ); \
  template void RegularizedQSDLDL \
  ( const DistSymmNodeInfo& info, DistSymmFront<F>& front, \
    Base<F> pivTol, \
    const DistMultiVecNode<Base<F>>& regCand, \
          DistMultiVecNode<Base<F>>& reg, \
    bool aPriori, \
    SymmFrontType newType );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
