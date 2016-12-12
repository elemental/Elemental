/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
void HermitianInverse
( UpperOrLower uplo, Matrix<Field>& A, const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    SymmetricInverse( uplo, A, true, ctrl );
}

template<typename Field>
void HermitianInverse
( UpperOrLower uplo, AbstractDistMatrix<Field>& A,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    SymmetricInverse( uplo, A, true, ctrl );
}

template<typename Field>
void LocalHermitianInverse
( UpperOrLower uplo, DistMatrix<Field,STAR,STAR>& A,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    SymmetricInverse( uplo, A.Matrix(), true, ctrl );
}

#define PROTO(Field) \
  template void HermitianInverse \
  ( UpperOrLower uplo, Matrix<Field>& A, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template void HermitianInverse \
  ( UpperOrLower uplo, AbstractDistMatrix<Field>& A, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template void LocalHermitianInverse \
  ( UpperOrLower uplo, DistMatrix<Field,STAR,STAR>& A, \
    const LDLPivotCtrl<Base<Field>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
