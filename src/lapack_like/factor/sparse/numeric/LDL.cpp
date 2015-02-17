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

#include "./LDL/Process.hpp"

namespace El {

template<typename F>
void LDL
( const SymmNodeInfo& info, SymmFront<F>& front, SymmFrontType newType )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    if( !Unfactored(front.type) )
        LogicError("Matrix is already factored");

    // Convert from 1D to 2D if necessary
    ChangeFrontType( front, SYMM_2D );

    // Perform the initial factorization
    ldl::Process( info, front, InitialFactorType(newType) );

    // Convert the fronts from the initial factorization to the requested form
    ChangeFrontType( front, newType );
}

template<typename F>
void LDL
( const DistSymmNodeInfo& info, DistSymmFront<F>& front, SymmFrontType newType )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    if( !Unfactored(front.type) )
        LogicError("Matrix is already factored");

    // Convert from 1D to 2D if necessary
    ChangeFrontType( front, SYMM_2D );

    // Perform the initial factorization
    ldl::Process( info, front, InitialFactorType(newType) );

    // Convert the fronts from the initial factorization to the requested form
    ChangeFrontType( front, newType );
}

#define PROTO(F) \
  template void LDL \
  ( const SymmNodeInfo& info, SymmFront<F>& front, \
    SymmFrontType newType ); \
  template void LDL \
  ( const DistSymmNodeInfo& info, DistSymmFront<F>& front, \
    SymmFrontType newType );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
