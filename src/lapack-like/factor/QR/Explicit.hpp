/*
   Copyright (c) 2009-2014, Jack Poulson
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
void Explicit( Matrix<F>& A, bool colPiv )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    Matrix<Int> p;
    Matrix<F> t;
    Matrix<Base<F>> d;
    if( colPiv )
        QR( A, t, d, p );
    else
        QR( A, t, d );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
}

template<typename F>
void Explicit( AbstractDistMatrix<F>& APre, bool colPiv )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    const Grid& g = APre.Grid();

    DistMatrix<F> A(g);
    Copy( APre, A, READ_WRITE_PROXY );

    DistMatrix<Int,VR,STAR> p(g);
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    if( colPiv )
        QR( A, t, d, p );
    else
        QR( A, t, d );

    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );

    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
}

template<typename F>
void Explicit( Matrix<F>& A, Matrix<F>& R, bool colPiv )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    Matrix<Int> p;
    if( colPiv )
        QR( A, t, d, p );
    else
        QR( A, t, d );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    auto AT = View( A, IndexRange(0,minDim), IndexRange(0,n) );
    R = AT;
    MakeTriangular( UPPER, R );

    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
}

template<typename F>
void Explicit
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& R, bool colPiv )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    const Grid& g = APre.Grid();

    DistMatrix<F> A(g);
    Copy( APre, A, READ_WRITE_PROXY );

    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    DistMatrix<Int,VR,STAR> p(g);
    if( colPiv )
        QR( A, t, d, p );
    else
        QR( A, t, d );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    auto AT = View( A, IndexRange(0,minDim), IndexRange(0,n) );
    Copy( AT, R );
    MakeTriangular( UPPER, R );

    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );

    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
}

template<typename F>
void Explicit( Matrix<F>& A, Matrix<F>& R, Matrix<Int>& p )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    QR( A, t, d, p );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    auto AT = View( A, IndexRange(0,minDim), IndexRange(0,n) );
    R = AT;
    MakeTriangular( UPPER, R );

    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
} 

template<typename F>
void Explicit
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& R, 
  AbstractDistMatrix<Int>& p )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    const Grid& g = APre.Grid();

    DistMatrix<F> A(g);
    Copy( APre, A, READ_WRITE_PROXY );

    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    QR( A, t, d, p );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    auto AT = View( A, IndexRange(0,minDim), IndexRange(0,n) );
    Copy( AT, R );
    MakeTriangular( UPPER, R );

    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );

    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_EXPLICIT_HPP
