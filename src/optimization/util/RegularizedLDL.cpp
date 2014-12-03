/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./RegularizedLDL/dense/Var3.hpp"

#include "./RegularizedLDL/sparse/ProcessFront.hpp"
#include "./RegularizedLDL/sparse/ProcessLocalTree.hpp"
#include "./RegularizedLDL/sparse/ProcessDistTree.hpp"

namespace El {

template<typename F>
void RegularizedLDL
( Matrix<F>& A, Base<F> pivTol, Base<F> regMag,
  const Matrix<Int>& pivSign, Matrix<Base<F>>& reg )
{
    DEBUG_ONLY(CallStackEntry cse("RegularizedLDL"))
    reg_ldl::Var3( A, pivTol, regMag, pivSign, reg );
}

template<typename F>
void RegularizedLDL
( AbstractDistMatrix<F>& A, Base<F> pivTol, Base<F> regMag,
  const AbstractDistMatrix<Int>& pivSign, AbstractDistMatrix<Base<F>>& reg )
{
    DEBUG_ONLY(CallStackEntry cse("RegularizedLDL"))
    reg_ldl::Var3( A, pivTol, regMag, pivSign, reg );
}

// TODO: Modify this routine to incorporate regularization
template<typename F>
void RegularizedLDL
( DistSymmInfo& info, DistSymmFrontTree<F>& L, SymmFrontType newFrontType )
{
    DEBUG_ONLY(CallStackEntry cse("RegularizedLDL"))
    LogicError("No regularization has yet been incorporated");
    if( !Unfactored(L.frontType) )
        LogicError("Matrix is already factored");

    // Convert from 1D to 2D if necessary
    ChangeFrontType( L, SYMM_2D );

    // Perform the initial factorization
    L.frontType = InitialFactorType(newFrontType);
    reg_ldl::ProcessLocalTree( info, L );
    reg_ldl::ProcessDistTree( info, L );

    // Convert the fronts from the initial factorization to the requested form
    ChangeFrontType( L, newFrontType );
}

#define PROTO(F) \
  template void RegularizedLDL \
  ( Matrix<F>& A, Base<F> pivTol, Base<F> regMag, \
    const Matrix<Int>& pivSign, Matrix<Base<F>>& reg ); \
  template void RegularizedLDL \
  ( AbstractDistMatrix<F>& A, Base<F> pivTol, Base<F> regMag, \
    const AbstractDistMatrix<Int>& pivSign, \
          AbstractDistMatrix<Base<F>>& reg ); \
  template void RegularizedLDL \
  ( DistSymmInfo& info, DistSymmFrontTree<F>& L, SymmFrontType newFrontType );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
