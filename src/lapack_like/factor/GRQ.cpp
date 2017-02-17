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
  Matrix<F>& householderScalarsA,
  Matrix<Base<F>>& signatureA,
  Matrix<F>& B,
  Matrix<F>& householderScalarsB,
  Matrix<Base<F>>& signatureB )
{
    EL_DEBUG_CSE
    RQ( A, householderScalarsA, signatureA );
    rq::ApplyQ( RIGHT, ADJOINT, A, householderScalarsA, signatureA, B );
    QR( B, householderScalarsB, signatureB );
}

template<typename F>
void GRQ
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& householderScalarsA,
  AbstractDistMatrix<Base<F>>& signatureA,
  AbstractDistMatrix<F>& BPre,
  AbstractDistMatrix<F>& householderScalarsB,
  AbstractDistMatrix<Base<F>>& signatureB )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre ), BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    RQ( A, householderScalarsA, signatureA );
    rq::ApplyQ( RIGHT, ADJOINT, A, householderScalarsA, signatureA, B );
    QR( B, householderScalarsB, signatureB );
}

namespace grq {

template<typename F>
void ExplicitTriang( Matrix<F>& A, Matrix<F>& B )
{
    EL_DEBUG_CSE
    Matrix<F> householderScalarsA;
    Matrix<Base<F>> signatureA;
    RQ( A, householderScalarsA, signatureA );
    rq::ApplyQ( RIGHT, ADJOINT, A, householderScalarsA, signatureA, B );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
    qr::ExplicitTriang( B );
}

template<typename F>
void ExplicitTriang( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& BPre )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre ), BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> householderScalarsA(g);
    DistMatrix<Base<F>,MD,STAR> signatureA(g);
    RQ( A, householderScalarsA, signatureA );
    rq::ApplyQ( RIGHT, ADJOINT, A, householderScalarsA, signatureA, B );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
    qr::ExplicitTriang( B );
}

} // namespace grq


#define PROTO(F) \
  template void GRQ \
  ( Matrix<F>& A, \
    Matrix<F>& householderScalarsA, \
    Matrix<Base<F>>& signatureA, \
    Matrix<F>& B, \
    Matrix<F>& householderScalarsB, \
    Matrix<Base<F>>& signatureB ); \
  template void GRQ \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& householderScalarsA, \
    AbstractDistMatrix<Base<F>>& signatureA, \
    AbstractDistMatrix<F>& B, \
    AbstractDistMatrix<F>& householderScalarsB, \
    AbstractDistMatrix<Base<F>>& signatureB ); \
  template void grq::ExplicitTriang \
  ( Matrix<F>& A, Matrix<F>& B ); \
  template void grq::ExplicitTriang \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
