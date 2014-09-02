/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Cholesky/LVar3.hpp"
#include "./Cholesky/LVar3Square.hpp"
#include "./Cholesky/LVar3Pivoted.hpp"
#include "./Cholesky/UVar3.hpp"
#include "./Cholesky/UVar3Square.hpp"
#include "./Cholesky/UVar3Pivoted.hpp"
#include "./Cholesky/SolveAfter.hpp"

#include "./Cholesky/LMod.hpp"
#include "./Cholesky/UMod.hpp"

namespace El {

// TODO: Pivoted Reverse Cholesky?

template<typename F>
void Cholesky( UpperOrLower uplo, Matrix<F>& A )
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
void Cholesky( UpperOrLower uplo, Matrix<F>& A, Matrix<Int>& pPerm )
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
void ReverseCholesky( UpperOrLower uplo, Matrix<F>& A )
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
void Cholesky( UpperOrLower uplo, AbstractDistMatrix<F>& A )
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

template<typename F> 
void Cholesky
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& pPerm )
{
    DEBUG_ONLY(CallStackEntry cse("Cholesky"))
    if( uplo == LOWER )
        cholesky::LVar3( A, pPerm );
    else
        cholesky::UVar3( A, pPerm );
}

template<typename F> 
void ReverseCholesky( UpperOrLower uplo, AbstractDistMatrix<F>& A )
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
void CholeskyMod( UpperOrLower uplo, Matrix<F>& T, Base<F> alpha, Matrix<F>& V )
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
void CholeskyMod
( UpperOrLower uplo, AbstractDistMatrix<F>& T, 
  Base<F> alpha, AbstractDistMatrix<F>& V )
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
void HPSDCholesky( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HPSDCholesky"))
    HPSDSquareRoot( uplo, A );
    MakeHermitian( uplo, A );

    if( uplo == LOWER )
        lq::ExplicitTriang( A );
    else
        qr::ExplicitTriang( A );
}

template<typename F>
void HPSDCholesky( UpperOrLower uplo, AbstractDistMatrix<F>& APre )
{
    DEBUG_ONLY(CallStackEntry cse("HPSDCholesky"))

    // NOTE: This should be removed once HPSD, LQ, and QR have been generalized
    auto APtr = ReadWriteProxy( &APre );
    auto& A = *APtr;

    HPSDSquareRoot( uplo, A );
    MakeHermitian( uplo, A );

    if( uplo == LOWER )
        lq::ExplicitTriang( A );
    else
        qr::ExplicitTriang( A );
}

#define PROTO(F) \
  template void Cholesky( UpperOrLower uplo, Matrix<F>& A ); \
  template void Cholesky( UpperOrLower uplo, AbstractDistMatrix<F>& A ); \
  template void ReverseCholesky( UpperOrLower uplo, Matrix<F>& A ); \
  template void ReverseCholesky \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A ); \
  template void Cholesky( UpperOrLower uplo, Matrix<F>& A, Matrix<Int>& p ); \
  template void Cholesky \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p ); \
  template void CholeskyMod \
  ( UpperOrLower uplo, Matrix<F>& T, Base<F> alpha, Matrix<F>& V ); \
  template void CholeskyMod \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& T, \
    Base<F> alpha, AbstractDistMatrix<F>& V ); \
  template void HPSDCholesky( UpperOrLower uplo, Matrix<F>& A ); \
  template void HPSDCholesky( UpperOrLower uplo, AbstractDistMatrix<F>& A ); \
  template void cholesky::SolveAfter \
  ( UpperOrLower uplo, Orientation orientation, \
    const Matrix<F>& A, Matrix<F>& B ); \
  template void cholesky::SolveAfter \
  ( UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B ); \
  template void cholesky::SolveAfter \
  ( UpperOrLower uplo, Orientation orientation, \
    const Matrix<F>& A, const Matrix<Int>& pPerm, Matrix<F>& B ); \
  template void cholesky::SolveAfter \
  ( UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<F>& A, const AbstractDistMatrix<Int>& pPerm, \
          AbstractDistMatrix<F>& B ); 

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
