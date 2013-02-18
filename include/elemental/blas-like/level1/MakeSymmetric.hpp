/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_MAKESYMMETRIC_HPP
#define BLAS_MAKESYMMETRIC_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/MakeReal.hpp"
#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level1/SetDiagonal.hpp"
#include "elemental/blas-like/level1/Transpose.hpp"

namespace elem {

template<typename T>
inline void
MakeSymmetric( UpperOrLower uplo, Matrix<T>& A, bool conjugate=false )
{
#ifndef RELEASE
    PushCallStack("MakeSymmetric");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make non-square matrix symmetric");

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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
MakeSymmetric( UpperOrLower uplo, DistMatrix<T>& A, bool conjugate=false )
{
#ifndef RELEASE
    PushCallStack("MakeSymmetric");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make non-square matrix symmetric");

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
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef BLAS_MAKESYMMETRIC_HPP
