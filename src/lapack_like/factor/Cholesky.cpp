/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Cholesky/LVar3.hpp"
#include "./Cholesky/LVar3Pivoted.hpp"
#include "./Cholesky/UVar3.hpp"
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
      CSE cse("Cholesky");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    if( uplo == LOWER )
        cholesky::LVar3( A );
    else
        cholesky::UVar3( A );
}

template<typename F>
void Cholesky( UpperOrLower uplo, Matrix<F>& A, Permutation& p )
{
    DEBUG_ONLY(
      CSE cse("Cholesky");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    if( uplo == LOWER )
        cholesky::LVar3( A, p );
    else
        cholesky::UVar3( A, p );
}

template<typename F>
void ReverseCholesky( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(
      CSE cse("ReverseCholesky");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    if( uplo == LOWER )
        cholesky::ReverseLVar3( A );
    else
        cholesky::ReverseUVar3( A );
}

namespace cholesky {

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
void ScaLAPACKHelper( UpperOrLower uplo, AbstractDistMatrix<F>& A )
{
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    // TODO: Add support for optionally timing the proxy redistribution
    DistMatrixReadWriteProxy<F,F,MC,MR,BLOCK> ABlockProx( A );
    auto& ABlock = ABlockProx.Get();

    const Int n = ABlock.Height();
    const int bHandle = blacs::Handle( ABlock );
    const int context = blacs::GridInit( bHandle, ABlock );
    auto descA = FillDesc( ABlock, context );

    const char uploChar = UpperOrLowerToChar( uplo );
    scalapack::Cholesky( uploChar, n, ABlock.Buffer(), descA.data() );

    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
void ScaLAPACKHelper( UpperOrLower uplo, AbstractDistMatrix<F>& A )
{
    RuntimeError("There is no ScaLAPACK support for this datatype");
}

} // anonymous namespace

template<typename F> 
void Cholesky( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool scalapack )
{
    DEBUG_ONLY(CSE cse("Cholesky"))
    if( scalapack )
    {
        cholesky::ScaLAPACKHelper( uplo, A );
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
( UpperOrLower uplo, AbstractDistMatrix<F>& A, DistPermutation& p )
{
    DEBUG_ONLY(CSE cse("Cholesky"))
    if( uplo == LOWER )
        cholesky::LVar3( A, p );
    else
        cholesky::UVar3( A, p );
}

template<typename F>
void Cholesky
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{ Cholesky( uplo, A.Matrix() ); }

template<typename F> 
void ReverseCholesky( UpperOrLower uplo, AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("ReverseCholesky"))
    if( uplo == LOWER )
        cholesky::ReverseLVar3( A );
    else
        cholesky::ReverseUVar3( A );
}

template<typename F>
void ReverseCholesky
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{ ReverseCholesky( uplo, A.Matrix() ); }

// Either 
//         L' L'^H := L L^H + alpha V V^H
// or
//         U'^H U' := U^H U + alpha V V^H
template<typename F>
void CholeskyMod( UpperOrLower uplo, Matrix<F>& T, Base<F> alpha, Matrix<F>& V )
{
    DEBUG_ONLY(CSE cse("CholeskyMod"))
    if( alpha == Base<F>(0) )
        return;
    if( uplo == LOWER )
        cholesky::LMod( T, alpha, V );
    else
        cholesky::UMod( T, alpha, V );
}

template<typename F>
void CholeskyMod
( UpperOrLower uplo,
  AbstractDistMatrix<F>& T, 
  Base<F> alpha,
  AbstractDistMatrix<F>& V )
{
    DEBUG_ONLY(CSE cse("CholeskyMod"))
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
    DEBUG_ONLY(CSE cse("HPSDCholesky"))
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
    DEBUG_ONLY(CSE cse("HPSDCholesky"))

    // NOTE: This should be removed once HPSD, LQ, and QR have been generalized
    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    HPSDSquareRoot( uplo, A );
    MakeHermitian( uplo, A );

    if( uplo == LOWER )
        lq::ExplicitTriang( A );
    else
        qr::ExplicitTriang( A );
}

#define PROTO_BASE(F) \
  template void Cholesky( UpperOrLower uplo, Matrix<F>& A ); \
  template void Cholesky \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool scalapack ); \
  template void Cholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A ); \
  template void ReverseCholesky( UpperOrLower uplo, Matrix<F>& A ); \
  template void ReverseCholesky \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A ); \
  template void ReverseCholesky \
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A ); \
  template void Cholesky( UpperOrLower uplo, Matrix<F>& A, Permutation& p ); \
  template void Cholesky \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<F>& A, \
    DistPermutation& p ); \
  template void CholeskyMod \
  ( UpperOrLower uplo, Matrix<F>& T, Base<F> alpha, Matrix<F>& V ); \
  template void CholeskyMod \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<F>& T, \
    Base<F> alpha, \
    AbstractDistMatrix<F>& V ); \
  template void cholesky::SolveAfter \
  ( UpperOrLower uplo, Orientation orientation, \
    const Matrix<F>& A, \
          Matrix<F>& B ); \
  template void cholesky::SolveAfter \
  ( UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<F>& A, \
          AbstractDistMatrix<F>& B ); \
  template void cholesky::SolveAfter \
  ( UpperOrLower uplo, Orientation orientation, \
    const Matrix<F>& A, \
    const Permutation& p, \
          Matrix<F>& B ); \
  template void cholesky::SolveAfter \
  ( UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<F>& A, \
    const DistPermutation& p, \
          AbstractDistMatrix<F>& B ); 

#define PROTO(F) \
  PROTO_BASE(F) \
  template void HPSDCholesky( UpperOrLower uplo, Matrix<F>& A ); \
  template void HPSDCholesky( UpperOrLower uplo, AbstractDistMatrix<F>& A );

#define PROTO_DOUBLEDOUBLE PROTO_BASE(DoubleDouble)
#define PROTO_QUADDOUBLE PROTO_BASE(QuadDouble)
#define PROTO_QUAD PROTO_BASE(Quad)
#define PROTO_COMPLEX_QUAD PROTO_BASE(Complex<Quad>)
#define PROTO_BIGFLOAT PROTO_BASE(BigFloat)

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
