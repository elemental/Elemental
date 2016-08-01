/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename F> 
void GRQ
( Matrix<F>& A,
  Matrix<F>& phaseA,
  Matrix<Base<F>>& signatureA, 
  Matrix<F>& B,
  Matrix<F>& phaseB,
  Matrix<Base<F>>& signatureB )
{
    DEBUG_CSE
    RQ( A, phaseA, signatureA );
    rq::ApplyQ( RIGHT, ADJOINT, A, phaseA, signatureA, B );
    QR( B, phaseB, signatureB );
}

template<typename F> 
void GRQ
( ElementalMatrix<F>& APre, 
  ElementalMatrix<F>& phaseA,
  ElementalMatrix<Base<F>>& signatureA,
  ElementalMatrix<F>& BPre, 
  ElementalMatrix<F>& phaseB,
  ElementalMatrix<Base<F>>& signatureB )
{
    DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre ), BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    RQ( A, phaseA, signatureA );
    rq::ApplyQ( RIGHT, ADJOINT, A, phaseA, signatureA, B );
    QR( B, phaseB, signatureB );
}

namespace grq {

template<typename F> 
void ExplicitTriang( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_CSE
    Matrix<F> phaseA;
    Matrix<Base<F>> signatureA;
    RQ( A, phaseA, signatureA );
    rq::ApplyQ( RIGHT, ADJOINT, A, phaseA, signatureA, B );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
    qr::ExplicitTriang( B );
}

template<typename F> 
void ExplicitTriang( ElementalMatrix<F>& APre, ElementalMatrix<F>& BPre )
{
    DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre ), BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> phaseA(g);
    DistMatrix<Base<F>,MD,STAR> signatureA(g);
    RQ( A, phaseA, signatureA );
    rq::ApplyQ( RIGHT, ADJOINT, A, phaseA, signatureA, B );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
    qr::ExplicitTriang( B );
}

} // namespace grq


#define PROTO(F) \
  template void GRQ \
  ( Matrix<F>& A, \
    Matrix<F>& phaseA, \
    Matrix<Base<F>>& signatureA, \
    Matrix<F>& B, \
    Matrix<F>& phaseB, \
    Matrix<Base<F>>& signatureB ); \
  template void GRQ \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& phaseA, \
    ElementalMatrix<Base<F>>& signatureA, \
    ElementalMatrix<F>& B, \
    ElementalMatrix<F>& phaseB, \
    ElementalMatrix<Base<F>>& signatureB ); \
  template void grq::ExplicitTriang \
  ( Matrix<F>& A, Matrix<F>& B ); \
  template void grq::ExplicitTriang \
  ( ElementalMatrix<F>& A, ElementalMatrix<F>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
