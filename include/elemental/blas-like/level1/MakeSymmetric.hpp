/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MAKESYMMETRIC_HPP
#define ELEM_MAKESYMMETRIC_HPP

#include "./Axpy.hpp"
#include "./MakeReal.hpp"
#include "./MakeTriangular.hpp"
#include "./SetDiagonal.hpp"
#include "./Transpose.hpp"

namespace elem {

template<typename T>
inline void
MakeSymmetric( UpperOrLower uplo, Matrix<T>& A, bool conjugate=false )
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
inline void
MakeSymmetric( UpperOrLower uplo, DistMatrix<T>& A, bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("MakeSymmetric"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    const Grid& g = A.Grid();
    DistMatrix<T,MD,STAR> d(g);
    A.GetDiagonal( d );
    if( conjugate )
        MakeReal( d );

    if( uplo == LOWER )
        MakeTriangular( LOWER, A );
    else
        MakeTriangular( UPPER, A );
    SetDiagonal( A, T(0) );
    DistMatrix<T> ATrans(g);
    Transpose( A, ATrans, conjugate );
    Axpy( T(1), ATrans, A );

    A.SetDiagonal( d );
}

} // namespace elem

#endif // ifndef ELEM_MAKESYMMETRIC_HPP
