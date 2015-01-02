/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_QR_EXPLICIT_HPP
#define EL_QR_EXPLICIT_HPP

namespace El {
namespace qr {

template<typename F>
void ExplicitTriang( Matrix<F>& A, const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ExplicitTriang"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    if( ctrl.colPiv )
    {
        Matrix<Int> p;
        BusingerGolub( A, t, d, p, ctrl );
    }
    else
        Householder( A, t, d );

    A.Resize( t.Height(), A.Width() );
    MakeTrapezoidal( UPPER, A );
}

template<typename F>
void ExplicitTriang( AbstractDistMatrix<F>& A, const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ExplicitTriang"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());
    if( ctrl.colPiv )
    {
        DistMatrix<Int,VC,STAR> p(A.Grid());
        BusingerGolub( A, t, d, p, ctrl );
    }
    else
        Householder( A, t, d );

    A.Resize( t.Height(), A.Width() );
    MakeTrapezoidal( UPPER, A );
}

template<typename F>
void ExplicitUnitary( Matrix<F>& A, const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ExplicitUnitary"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    if( ctrl.colPiv )
    {
        Matrix<Int> p;
        QR( A, t, d, p, ctrl );
    }
    else
        QR( A, t, d );

    A.Resize( A.Height(), t.Height() );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
}

template<typename F>
void ExplicitUnitary( AbstractDistMatrix<F>& APre, const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ExplicitUnitary"))

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    if( ctrl.colPiv )
    {
        DistMatrix<Int,VR,STAR> p(g);
        QR( A, t, d, p, ctrl );
    }
    else
        QR( A, t, d );

    A.Resize( A.Height(), t.Height() );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
}

template<typename F>
void Explicit( Matrix<F>& A, Matrix<F>& R, const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    if( ctrl.colPiv )
    {
        Matrix<Int> p;
        QR( A, t, d, p, ctrl );
    }
    else
        QR( A, t, d );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numIts = t.Height();

    auto AT = A( IR(0,numIts), IR(0,n) );
    R = AT;
    MakeTrapezoidal( UPPER, R );

    A.Resize( m, numIts );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
}

template<typename F>
void Explicit
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& R, 
  const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    if( ctrl.colPiv )
    {
        DistMatrix<Int,VR,STAR> p(g);
        QR( A, t, d, p, ctrl );
    }
    else
        QR( A, t, d );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numIts = t.Height();

    auto AT = A( IR(0,numIts), IR(0,n) );
    Copy( AT, R );
    MakeTrapezoidal( UPPER, R );

    A.Resize( m, numIts );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
}

template<typename F>
void Explicit
( Matrix<F>& A, Matrix<F>& R, Matrix<Int>& P, const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    Matrix<Int> p;
    QR( A, t, d, p, ctrl );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numIts = t.Height();

    auto AT = A( IR(0,numIts), IR(0,n) );
    R = AT;
    MakeTrapezoidal( UPPER, R );

    A.Resize( m, numIts );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );

    ExplicitPermutation( p, P );
} 

template<typename F>
void Explicit
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& R, 
  AbstractDistMatrix<Int>& P, const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    DistMatrix<Int,VC,STAR> p(g);
    QR( A, t, d, p, ctrl );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numIts = t.Height();

    auto AT = A( IR(0,numIts), IR(0,n) );
    Copy( AT, R );
    MakeTrapezoidal( UPPER, R );

    A.Resize( m, numIts );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );

    ExplicitPermutation( p, P );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_EXPLICIT_HPP
