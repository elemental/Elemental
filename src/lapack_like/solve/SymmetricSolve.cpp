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
void SymmetricSolve
( UpperOrLower uplo, Orientation orientation, Matrix<F>& A, 
  Matrix<F>& B, bool conjugate, const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricSolve"))
    if( uplo == UPPER )
        LogicError("Upper Bunch-Kaufman is not yet supported");
    Matrix<Int> p; 
    Matrix<F> dSub;
    LDL( A, dSub, p, conjugate, ctrl );
    const bool conjFlip = ( (orientation == ADJOINT && conjugate == false) ||
                            (orientation == TRANSPOSE && conjugate == true) );
    if( conjFlip )
        Conjugate( B );
    ldl::SolveAfter( A, dSub, p, B, conjugate );
    if( conjFlip )
        Conjugate( B );
}

template<typename F>
void SymmetricSolve
( UpperOrLower uplo, Orientation orientation, AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& BPre, bool conjugate, 
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricSolve"))
    if( uplo == UPPER )
        LogicError("Upper Bunch-Kaufman is not yet supported");

    auto APtr = ReadProxy<F,MC,MR>( &APre );      auto& A = *APtr;
    auto BPtr = ReadWriteProxy<F,MC,MR>( &BPre ); auto& B = *BPtr;

    DistMatrix<Int,VC,STAR> p(A.Grid()); 
    DistMatrix<F,MD,STAR> dSub(A.Grid());
    LDL( A, dSub, p, conjugate, ctrl );
    const bool conjFlip = ( (orientation == ADJOINT && conjugate == false) ||
                            (orientation == TRANSPOSE && conjugate == true) );
    if( conjFlip )
        Conjugate( B );
    ldl::SolveAfter( A, dSub, p, B, conjugate );
    if( conjFlip )
        Conjugate( B );
}

#define PROTO(F) \
  template void SymmetricSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    Matrix<F>& A, Matrix<F>& B, bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void SymmetricSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
