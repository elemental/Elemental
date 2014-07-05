/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void MakeSymmetric( UpperOrLower uplo, Matrix<T>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("MakeSymmetric"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    Matrix<T> d;
    A.GetDiagonal( d );
    if( conjugate )
        MakeReal( d );

    if( uplo == LOWER )
        MakeTriangular( LOWER, A );
    else
        MakeTriangular( UPPER, A );
    SetDiagonal( A, T(0) );
    Matrix<T> ATrans;
    Transpose( A, ATrans, conjugate );
    Axpy( T(1), ATrans, A );

    A.SetDiagonal( d );
}

template<typename T>
void MakeSymmetric
( UpperOrLower uplo, AbstractDistMatrix<T>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("MakeSymmetric"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    std::unique_ptr<AbstractDistMatrix<T>> ATrans( A.Construct(A.Grid()) );
    Transpose( A, *ATrans, conjugate );
    if( uplo == LOWER )
    {
        MakeTriangular( LOWER, A );
        MakeTrapezoidal( UPPER, *ATrans, 1 );
    }
    else
    {
        MakeTriangular( UPPER, A );
        MakeTrapezoidal( LOWER, *ATrans, -1 );
    }
    Axpy( T(1), *ATrans, A );
    A.MakeDiagonalReal();
}

#define PROTO(F) \
  template void MakeSymmetric \
  ( UpperOrLower uplo, Matrix<F>& A, bool conjugate ); \
  template void MakeSymmetric \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool conjugate ); 

#include "El/macros/Instantiate.h"

} // namespace El
