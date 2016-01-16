/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
InertiaType Inertia
( UpperOrLower uplo,
  Matrix<F>& A,
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("Inertia"))
    if( uplo == UPPER )
        LogicError("This option not yet supported");
    Permutation p;
    Matrix<F> dSub;
    LDL( A, dSub, p, true, ctrl );
    return ldl::Inertia( GetRealPartOfDiagonal(A), dSub );
}

template<typename F>
InertiaType Inertia
( UpperOrLower uplo,
  ElementalMatrix<F>& APre, 
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("Inertia"))
    if( uplo == UPPER )
        LogicError("This option not yet supported");

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    DistPermutation p( A.Grid() );
    DistMatrix<F,MD,STAR> dSub( A.Grid() );
    LDL( A, dSub, p, true, ctrl );
    return ldl::Inertia( GetRealPartOfDiagonal(A), dSub );
}

#define PROTO(F) \
  template InertiaType Inertia \
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template InertiaType Inertia \
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    const LDLPivotCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
