/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// NOTE: This overwrites both triangles of the inverse.
template<typename F>
void SymmetricInverse
( UpperOrLower uplo, Matrix<F>& A, bool conjugate, 
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricInverse"))
    if( uplo == LOWER )
    {
        Matrix<Int> p;
        Matrix<F> dSub;
        LDL( A, dSub, p, conjugate, ctrl );
        TriangularInverse( LOWER, UNIT, A ); 
        Trdtrmm( LOWER, A, dSub, conjugate );

        // NOTE: Fill in both triangles of the inverse
        Matrix<Int> pInv;
        InvertPermutation( p, pInv );
        MakeSymmetric( LOWER, A, conjugate );
        PermuteRows( A, pInv, p );
        PermuteCols( A, pInv, p ); 
    }
    else
        LogicError("This option is not yet supported");
}

template<typename F>
void SymmetricInverse
( UpperOrLower uplo, AbstractDistMatrix<F>& APre, bool conjugate, 
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricInverse"))

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    if( uplo == LOWER )
    {
        DistMatrix<Int,VC,STAR> p( A.Grid() );
        DistMatrix<F,MD,STAR> dSub( A.Grid() );

        LDL( A, dSub, p, conjugate, ctrl );
        TriangularInverse( LOWER, UNIT, A ); 
        Trdtrmm( LOWER, A, dSub, conjugate );

        // NOTE: Fill in both triangles of the inverse
        DistMatrix<Int,VC,STAR> pInv(p.Grid());
        InvertPermutation( p, pInv );
        MakeSymmetric( LOWER, A, conjugate );
        PermuteRows( A, pInv, p );
        PermuteCols( A, pInv, p );
    }
    else
        LogicError("This option is not yet supported");
}

template<typename F>
void LocalSymmetricInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, bool conjugate, 
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LocalSymmetricInverse"))
    SymmetricInverse( uplo, A.Matrix(), conjugate, ctrl );
}

#define PROTO(F) \
  template void SymmetricInverse \
  ( UpperOrLower uplo, Matrix<F>& A, bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void SymmetricInverse \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void LocalSymmetricInverse \
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
