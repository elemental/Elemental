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
InertiaType Inertia
( UpperOrLower uplo, Matrix<F>& A, const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("Inertia"))
    if( uplo == UPPER )
        LogicError("This option not yet supported");
    Matrix<Int> p;
    Matrix<F> dSub;
    LDL( A, dSub, p, true, ctrl );
    return ldl::Inertia( GetRealPartOfDiagonal(A), dSub );
}

template<typename F>
InertiaType Inertia
( UpperOrLower uplo, AbstractDistMatrix<F>& APre, 
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("Inertia"))
    if( uplo == UPPER )
        LogicError("This option not yet supported");

    auto APtr = ReadProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    DistMatrix<Int,VC,STAR> p( A.Grid() );
    DistMatrix<F,MD,STAR> dSub( A.Grid() );
    LDL( A, dSub, p, true, ctrl );
    return ldl::Inertia( GetRealPartOfDiagonal(A), dSub );
}

#define PROTO(F) \
  template InertiaType Inertia \
  ( UpperOrLower uplo, Matrix<F>& A, const LDLPivotCtrl<Base<F>>& ctrl ); \
  template InertiaType Inertia \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, \
    const LDLPivotCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
