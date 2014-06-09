/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename F>
InertiaType Inertia
( UpperOrLower uplo, Matrix<F>& A, LDLPivotType pivotType )
{
    DEBUG_ONLY(CallStackEntry cse("Inertia"))
    if( uplo == UPPER )
        LogicError("This option not yet supported");
    Matrix<Int> pPerm;
    Matrix<F> dSub;
    LDL( A, dSub, pPerm, true, pivotType );
    return ldl::Inertia( A.GetRealPartOfDiagonal(), dSub );
}

template<typename F>
InertiaType Inertia
( UpperOrLower uplo, DistMatrix<F>& A, LDLPivotType pivotType )
{
    DEBUG_ONLY(CallStackEntry cse("Inertia"))
    if( uplo == UPPER )
        LogicError("This option not yet supported");
    DistMatrix<Int,VC,STAR> pPerm( A.Grid() );
    DistMatrix<F,MD,STAR> dSub( A.Grid() );
    LDL( A, dSub, pPerm, true, pivotType );
    return ldl::Inertia( A.GetRealPartOfDiagonal(), dSub );
}

#define PROTO(F) \
  template InertiaType Inertia \
  ( UpperOrLower uplo, Matrix<F>& A, LDLPivotType pivotType ); \
  template InertiaType Inertia \
  ( UpperOrLower uplo, DistMatrix<F>& A, LDLPivotType pivotType );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
