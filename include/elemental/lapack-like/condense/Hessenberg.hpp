/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HESSENBERG_HPP
#define ELEM_HESSENBERG_HPP

#include ELEM_MAKETRAPEZOIDAL_INC
#include ELEM_MAKETRIANGULAR_INC

#include "./Hessenberg/L.hpp"
#include "./Hessenberg/U.hpp"
#include "./Hessenberg/ApplyQ.hpp"

namespace elem {

template<typename F>
inline void Hessenberg( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("Hessenberg"))
    if( uplo == UPPER )
        hessenberg::U( A, t );
    else
        hessenberg::L( A, t );
}

template<typename F> 
inline void Hessenberg
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
    DEBUG_ONLY(CallStackEntry cse("Hessenberg"))
    if( uplo == UPPER )
        hessenberg::U( A, t );
    else
        hessenberg::L( A, t );
}

template<typename F>
inline void Hessenberg( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Hessenberg"))
    Matrix<F> t;
    Hessenberg( uplo, A, t );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

template<typename F> 
inline void Hessenberg( UpperOrLower uplo, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Hessenberg"))
    DistMatrix<F,STAR,STAR> t(A.Grid());
    Hessenberg( uplo, A, t );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

} // namespace elem

#endif // ifndef ELEM_HESSENBERG_HPP
