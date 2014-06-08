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
void SymmetricSolve
( UpperOrLower uplo, Orientation orientation, Matrix<F>& A, 
  Matrix<F>& B, bool conjugate, LDLPivotType pivotType )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricSolve"))
    if( uplo == UPPER )
        LogicError("Upper Bunch-Kaufman is not yet supported");
    Matrix<Int> pPerm; 
    Matrix<F> dSub;
    LDL( A, dSub, pPerm, conjugate, pivotType );
    const bool conjFlip = ( (orientation == ADJOINT && conjugate == false) ||
                            (orientation == TRANSPOSE && conjugate == true) );
    if( conjFlip )
        Conjugate( B );
    ldl::SolveAfter( A, dSub, pPerm, B, conjugate );
    if( conjFlip )
        Conjugate( B );
}

template<typename F>
void SymmetricSolve
( UpperOrLower uplo, Orientation orientation, DistMatrix<F>& A, 
  DistMatrix<F>& B, bool conjugate, LDLPivotType pivotType )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricSolve"))
    if( uplo == UPPER )
        LogicError("Upper Bunch-Kaufman is not yet supported");
    DistMatrix<Int,VC,STAR> pPerm(A.Grid()); 
    DistMatrix<F,MD,STAR> dSub(A.Grid());
    LDL( A, dSub, pPerm, conjugate, pivotType );
    const bool conjFlip = ( (orientation == ADJOINT && conjugate == false) ||
                            (orientation == TRANSPOSE && conjugate == true) );
    if( conjFlip )
        Conjugate( B );
    ldl::SolveAfter( A, dSub, pPerm, B, conjugate );
    if( conjFlip )
        Conjugate( B );
}

#define PROTO(F) \
  template void SymmetricSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    Matrix<F>& A, Matrix<F>& B, bool conjugate, \
    LDLPivotType pivotType ); \
  template void SymmetricSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    DistMatrix<F>& A, DistMatrix<F>& B, bool conjugate, \
    LDLPivotType pivotType );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
