/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CHOLESKY_HPP
#define ELEM_CHOLESKY_HPP

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

#include "./Cholesky/LMod.hpp"
#include "./Cholesky/UMod.hpp"

namespace elem {

template<typename F>
inline void
LocalCholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalCholesky"))
    Cholesky( uplo, A.Matrix() );
}

template<typename F>
inline void
LocalReverseCholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalReverseCholesky"))
    ReverseCholesky( uplo, A.Matrix() );
}

// TODO: Pivoted Reverse Cholesky?

template<typename F>
inline void
Cholesky( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("Cholesky");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        cholesky::LVar3( A );
    else
        cholesky::UVar3( A );
}

template<typename F>
inline void
Cholesky( UpperOrLower uplo, Matrix<F>& A, Matrix<Int>& pPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("Cholesky");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        cholesky::LVar3( A, pPerm );
    else
        cholesky::UVar3( A, pPerm );
}

template<typename F>
inline void
ReverseCholesky( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("ReverseCholesky");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        cholesky::ReverseLVar3( A );
    else
        cholesky::ReverseUVar3( A );
}

template<typename F> 
inline void
Cholesky( UpperOrLower uplo, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Cholesky"))
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

template<typename F,Dist UPerm> 
inline void
Cholesky
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm )
{
    DEBUG_ONLY(CallStackEntry cse("Cholesky"))
    if( uplo == LOWER )
        cholesky::LVar3( A, pPerm );
    else
        cholesky::UVar3( A, pPerm );
}

template<typename F> 
inline void
ReverseCholesky( UpperOrLower uplo, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("ReverseCholesky"))
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

// Either 
//         L' L'^H := L L^H + alpha V V^H
// or
//         U'^H U' := U^H U + alpha V V^H
template<typename F>
inline void
CholeskyMod( UpperOrLower uplo, Matrix<F>& T, Base<F> alpha, Matrix<F>& V )
{
    DEBUG_ONLY(CallStackEntry cse("CholeskyMod"))
    if( alpha == Base<F>(0) )
        return;
    if( uplo == LOWER )
        cholesky::LMod( T, alpha, V );
    else
        cholesky::UMod( T, alpha, V );
}

template<typename F>
inline void
CholeskyMod
( UpperOrLower uplo, DistMatrix<F>& T, Base<F> alpha, DistMatrix<F>& V )
{
    DEBUG_ONLY(CallStackEntry cse("CholeskyMod"))
    if( alpha == Base<F>(0) )
        return;
    if( uplo == LOWER )
        cholesky::LMod( T, alpha, V );
    else
        cholesky::UMod( T, alpha, V );
}

} // namespace elem

#include ELEM_MAKEHERMITIAN_INC
#include ELEM_MAKETRIANGULAR_INC

#include ELEM_LQ_INC
#include ELEM_QR_INC

#include ELEM_SQUAREROOT_INC

namespace elem {

template<typename F>
inline void
HPSDCholesky( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HPSDCholesky"))
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
    DEBUG_ONLY(CallStackEntry cse("HPSDCholesky"))

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

#endif // ifndef ELEM_CHOLESKY_HPP
