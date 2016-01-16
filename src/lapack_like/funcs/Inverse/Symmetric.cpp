/*
   Copyright (c) 2009-2016, Jack Poulson
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
( UpperOrLower uplo,
  Matrix<F>& A,
  bool conjugate, 
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SymmetricInverse"))
    if( uplo == LOWER )
    {
        Permutation P;
        Matrix<F> dSub;
        LDL( A, dSub, P, conjugate, ctrl );
        TriangularInverse( LOWER, UNIT, A ); 
        Trdtrmm( LOWER, A, dSub, conjugate );

        // NOTE: Fill in both triangles of the inverse
        MakeSymmetric( LOWER, A, conjugate );
        P.InversePermuteRows( A );
        P.InversePermuteCols( A );
    }
    else
        LogicError("This option is not yet supported");
}

template<typename F>
void SymmetricInverse
( UpperOrLower uplo,
  ElementalMatrix<F>& APre,
  bool conjugate, 
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SymmetricInverse"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    if( uplo == LOWER )
    {
        DistPermutation P( A.Grid() );
        DistMatrix<F,MD,STAR> dSub( A.Grid() );

        LDL( A, dSub, P, conjugate, ctrl );
        TriangularInverse( LOWER, UNIT, A ); 
        Trdtrmm( LOWER, A, dSub, conjugate );

        // NOTE: Fill in both triangles of the inverse
        MakeSymmetric( LOWER, A, conjugate );
        P.InversePermuteRows( A );
        P.InversePermuteCols( A );
    }
    else
        LogicError("This option is not yet supported");
}

template<typename F>
void LocalSymmetricInverse
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A,
  bool conjugate, 
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LocalSymmetricInverse"))
    SymmetricInverse( uplo, A.Matrix(), conjugate, ctrl );
}

#define PROTO(F) \
  template void SymmetricInverse \
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void SymmetricInverse \
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void LocalSymmetricInverse \
  ( UpperOrLower uplo, \
    DistMatrix<F,STAR,STAR>& A, \
    bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
