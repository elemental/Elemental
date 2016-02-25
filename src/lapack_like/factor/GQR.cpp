/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F> 
void GQR
( Matrix<F>& A,
  Matrix<F>& tA,
  Matrix<Base<F>>& dA, 
  Matrix<F>& B,
  Matrix<F>& tB,
  Matrix<Base<F>>& dB )
{
    DEBUG_ONLY(CSE cse("GQR"))
    QR( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    RQ( B, tB, dB );
}

template<typename F> 
void GQR
( ElementalMatrix<F>& APre, 
  ElementalMatrix<F>& tA,
  ElementalMatrix<Base<F>>& dA,
  ElementalMatrix<F>& BPre, 
  ElementalMatrix<F>& tB,
  ElementalMatrix<Base<F>>& dB )
{
    DEBUG_ONLY(CSE cse("GQR"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre ), BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    QR( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    RQ( B, tB, dB );
}

namespace gqr {

template<typename F> 
void ExplicitTriang( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CSE cse("gqr::ExplicitTriang"))
    Matrix<F> tA;
    Matrix<Base<F>> dA;
    QR( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    MakeTrapezoidal( UPPER, A );
    rq::ExplicitTriang( B );
}

template<typename F> 
void ExplicitTriang( ElementalMatrix<F>& APre, ElementalMatrix<F>& BPre )
{
    DEBUG_ONLY(CSE cse("gqr::ExplicitTriang"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre ), BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> tA(g);
    DistMatrix<Base<F>,MD,STAR> dA(g);
    QR( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    MakeTrapezoidal( UPPER, A );
    rq::ExplicitTriang( B );
}

} // namespace gqr

#define PROTO(F) \
  template void GQR \
  ( Matrix<F>& A, \
    Matrix<F>& tA, \
    Matrix<Base<F>>& dA, \
    Matrix<F>& B, \
    Matrix<F>& tB, \
    Matrix<Base<F>>& dB ); \
  template void GQR \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& tA, \
    ElementalMatrix<Base<F>>& dA, \
    ElementalMatrix<F>& B, \
    ElementalMatrix<F>& tB, \
    ElementalMatrix<Base<F>>& dB ); \
  template void gqr::ExplicitTriang \
  ( Matrix<F>& A, Matrix<F>& B ); \
  template void gqr::ExplicitTriang \
  ( ElementalMatrix<F>& A, ElementalMatrix<F>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
