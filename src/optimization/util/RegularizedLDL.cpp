/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./RegularizedLDL/Var3.hpp"

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

#define PROTO(F) \
  template void RegularizedLDL \
  ( Matrix<F>& A, Base<F> pivTol, Base<F> regMag, \
    const Matrix<Int>& pivSign, Matrix<Base<F>>& reg ); \
  template void RegularizedLDL \
  ( AbstractDistMatrix<F>& A, Base<F> pivTol, Base<F> regMag, \
    const AbstractDistMatrix<Int>& pivSign, \
          AbstractDistMatrix<Base<F>>& reg );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
