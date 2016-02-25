/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_QR_EXPLICIT_HPP
#define EL_QR_EXPLICIT_HPP

namespace El {
namespace qr {

template<typename F>
void ExplicitTriang( Matrix<F>& A, const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("qr::ExplicitTriang"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    if( ctrl.colPiv )
    {
        Permutation Omega;
        BusingerGolub( A, t, d, Omega, ctrl );
    }
    else
        Householder( A, t, d );

    A.Resize( t.Height(), A.Width() );
    MakeTrapezoidal( UPPER, A );
}

template<typename F>
void ExplicitTriang( ElementalMatrix<F>& A, const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("qr::ExplicitTriang"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());
    if( ctrl.colPiv )
    {
        DistPermutation Omega(A.Grid());
        BusingerGolub( A, t, d, Omega, ctrl );
    }
    else
        Householder( A, t, d );

    A.Resize( t.Height(), A.Width() );
    MakeTrapezoidal( UPPER, A );
}

template<typename F>
void ExplicitUnitary
( Matrix<F>& A, bool thinQR, const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("qr::ExplicitUnitary"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    if( ctrl.colPiv )
    {
        Permutation Omega;
        QR( A, t, d, Omega, ctrl );
    }
    else
        QR( A, t, d );

    if( thinQR ) 
    {
        A.Resize( A.Height(), t.Height() );
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
        DiagonalScale( RIGHT, NORMAL, d, A );
    }
    else
    {
        auto ACopy = A;
        // TODO: Use an extension of ExpandPackedReflectors to make this faster
        Identity( A, A.Height(), A.Height() );
        qr::ApplyQ( LEFT, NORMAL, ACopy, t, d, A );
    }
}

template<typename F>
void ExplicitUnitary
( ElementalMatrix<F>& APre, bool thinQR, const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("qr::ExplicitUnitary"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    if( ctrl.colPiv )
    {
        DistPermutation Omega(g);
        QR( A, t, d, Omega, ctrl );
    }
    else
        QR( A, t, d );

    if( thinQR )
    {
        A.Resize( A.Height(), t.Height() );
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
        DiagonalScale( RIGHT, NORMAL, d, A );
    }
    else
    {
        auto ACopy = A;
        // TODO: Use an extension of ExpandPackedReflectors to make this faster
        Identity( A, A.Height(), A.Height() );
        qr::ApplyQ( LEFT, NORMAL, ACopy, t, d, A );
    }
}

template<typename F>
void Explicit
( Matrix<F>& A,
  Matrix<F>& R,
  bool thinQR,
  const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("qr::Explicit"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    if( ctrl.colPiv )
    {
        Permutation Omega;
        QR( A, t, d, Omega, ctrl );
    }
    else
        QR( A, t, d );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numIts = t.Height();

    auto AT = A( IR(0,numIts), IR(0,n) );
    R = AT;
    MakeTrapezoidal( UPPER, R );

    if( thinQR )
    {
        A.Resize( m, numIts );
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
        DiagonalScale( RIGHT, NORMAL, d, A );
    }
    else
    {
        auto ACopy = A;
        // TODO: Use an extension of ExpandPackedReflectors to make this faster
        Identity( A, A.Height(), A.Height() );
        qr::ApplyQ( LEFT, NORMAL, ACopy, t, d, A );
    }
}

template<typename F>
void Explicit
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& R, 
  bool thinQR,
  const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("qr::Explicit"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    if( ctrl.colPiv )
    {
        DistPermutation Omega(g);
        QR( A, t, d, Omega, ctrl );
    }
    else
        QR( A, t, d );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numIts = t.Height();

    auto AT = A( IR(0,numIts), IR(0,n) );
    Copy( AT, R );
    MakeTrapezoidal( UPPER, R );

    if( thinQR )
    {
        A.Resize( m, numIts );
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
        DiagonalScale( RIGHT, NORMAL, d, A );
    }
    else
    {
        auto ACopy = A;
        // TODO: Use an extension of ExpandPackedReflectors to make this faster
        Identity( A, A.Height(), A.Height() );
        qr::ApplyQ( LEFT, NORMAL, ACopy, t, d, A );
    }
}

template<typename F>
void Explicit
( Matrix<F>& A,
  Matrix<F>& R,
  Matrix<Int>& OmegaFull,
  bool thinQR,
  const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("qr::Explicit"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    Permutation Omega;
    QR( A, t, d, Omega, ctrl );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numIts = t.Height();

    auto AT = A( IR(0,numIts), IR(0,n) );
    R = AT;
    MakeTrapezoidal( UPPER, R );

    if( thinQR )
    {
        A.Resize( m, numIts );
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
        DiagonalScale( RIGHT, NORMAL, d, A );
    }
    else
    {
        auto ACopy = A;
        // TODO: Use an extension of ExpandPackedReflectors to make this faster
        Identity( A, A.Height(), A.Height() );
        qr::ApplyQ( LEFT, NORMAL, ACopy, t, d, A );
    }

    Omega.ExplicitMatrix( OmegaFull );
} 

template<typename F>
void Explicit
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& R, 
  ElementalMatrix<Int>& OmegaFull,
  bool thinQR,
  const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("qr::Explicit"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    DistPermutation Omega(g);
    QR( A, t, d, Omega, ctrl );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numIts = t.Height();

    auto AT = A( IR(0,numIts), IR(0,n) );
    Copy( AT, R );
    MakeTrapezoidal( UPPER, R );

    if( thinQR )
    {
        A.Resize( m, numIts );
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
        DiagonalScale( RIGHT, NORMAL, d, A );
    }
    else
    {
        auto ACopy = A;
        // TODO: Use an extension of ExpandPackedReflectors to make this faster
        Identity( A, A.Height(), A.Height() );
        qr::ApplyQ( LEFT, NORMAL, ACopy, t, d, A );
    }

    Omega.ExplicitMatrix( OmegaFull );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_EXPLICIT_HPP
