/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
void HermitianInverse
( UpperOrLower uplo, Matrix<F>& A, const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianInverse"))
    SymmetricInverse( uplo, A, true, ctrl );
}

template<typename F>
void HermitianInverse
( UpperOrLower uplo, AbstractDistMatrix<F>& A, 
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianInverse"))
    SymmetricInverse( uplo, A, true, ctrl );
}

template<typename F>
void LocalHermitianInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LocalHermitianInverse"))
    SymmetricInverse( uplo, A.Matrix(), true, ctrl );
}

#define PROTO(F) \
  template void HermitianInverse \
  ( UpperOrLower uplo, Matrix<F>& A, const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void HermitianInverse \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void LocalHermitianInverse \
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, \
    const LDLPivotCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
