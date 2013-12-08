/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_CHOLESKY_HPP
#define ELEM_LAPACK_CHOLESKY_HPP

// TODO: Reorganize Cholesky implementation?
namespace elem {
template<typename F>
void LocalCholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A );
template<typename F>
void LocalReverseCholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A );
} // namespace elem

#include "./Cholesky/LVar3.hpp"
#include "./Cholesky/LVar3Square.hpp"
#include "./Cholesky/LVar3Pivoted.hpp"
#include "./Cholesky/UVar3.hpp"
#include "./Cholesky/UVar3Square.hpp"
#include "./Cholesky/UVar3Pivoted.hpp"
#include "./Cholesky/SolveAfter.hpp"

namespace elem {

template<typename F>
inline void
LocalCholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LocalCholesky");
#endif
    Cholesky( uplo, A.Matrix() );
}

template<typename F>
inline void
LocalReverseCholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LocalReverseCholesky");
#endif
    ReverseCholesky( uplo, A.Matrix() );
}

// TODO: Pivoted Reverse Cholesky?

template<typename F>
inline void
Cholesky( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("Cholesky");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    if( uplo == LOWER )
        cholesky::LVar3( A );
    else
        cholesky::UVar3( A );
}

template<typename F>
inline void
Cholesky( UpperOrLower uplo, Matrix<F>& A, Matrix<Int>& p )
{
#ifndef RELEASE
    CallStackEntry cse("Cholesky");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    if( uplo == LOWER )
        cholesky::LVar3( A, p );
    else
        cholesky::UVar3( A, p );
}

template<typename F>
inline void
ReverseCholesky( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("ReverseCholesky");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    if( uplo == LOWER )
        cholesky::ReverseLVar3( A );
    else
        cholesky::ReverseUVar3( A );
}

template<typename F> 
inline void
Cholesky( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("Cholesky");
#endif
    const Grid& g = A.Grid();
    if( g.Height() == g.Width() )
    {
        if( uplo == LOWER )
            cholesky::LVar3Square( A );
        else
            cholesky::UVar3Square( A );
    }
    else
    {
        if( uplo == LOWER )
            cholesky::LVar3( A );
        else
            cholesky::UVar3( A );
    }
}

template<typename F> 
inline void
Cholesky( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Int,VC,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry cse("Cholesky");
#endif
    if( uplo == LOWER )
        cholesky::LVar3( A, p );
    else
        cholesky::UVar3( A, p );
}

template<typename F> 
inline void
ReverseCholesky( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("ReverseCholesky");
#endif
    if( uplo == LOWER )
        cholesky::ReverseLVar3( A );
    else
        cholesky::ReverseUVar3( A );
    /*
    const Grid& g = A.Grid();
    if( g.Height() == g.Width() )
    {
        if( uplo == LOWER )
            cholesky::ReverseLVar3Square( A );
        else
            cholesky::ReverseUVar3Square( A );
    }
    else
    {
        if( uplo == LOWER )
            cholesky::ReverseLVar3( A );
        else
            cholesky::ReverseUVar3( A );
    }
    */
}

} // namespace elem

#include "elemental/blas-like/level1/MakeHermitian.hpp"
#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/lapack-like/SquareRoot.hpp"
#include "elemental/lapack-like/LQ.hpp"
#include "elemental/lapack-like/QR.hpp"

namespace elem {

template<typename F>
inline void
HPSDCholesky( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("HPSDCholesky");
#endif
    HPSDSquareRoot( uplo, A );
    MakeHermitian( uplo, A );

    if( uplo == LOWER )
    {
        LQ( A );
        MakeTriangular( LOWER, A );
    }
    else
    {
        QR( A );
        MakeTriangular( UPPER, A );
    }
}

template<typename F>
inline void
HPSDCholesky( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("HPSDCholesky");
#endif
    EnsurePMRRR();

    HPSDSquareRoot( uplo, A );
    MakeHermitian( uplo, A );

    if( uplo == LOWER )
    {
        LQ( A );
        MakeTriangular( LOWER, A );
    }
    else
    {
        QR( A );
        MakeTriangular( UPPER, A );
    }
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_CHOLESKY_HPP
