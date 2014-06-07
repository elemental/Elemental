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
    Matrix<Int> pPerm;
    Matrix<F> t;
    Matrix<Base<F>> d;
    if( colPiv )
        QR( A, t, d, pPerm );
    else
        QR( A, t, d );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
}

template<typename F>
void Explicit( DistMatrix<F>& A, bool colPiv )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    const Grid& g = A.Grid();
    DistMatrix<Int,VR,STAR> pPerm(g);
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    if( colPiv )
        QR( A, t, d, pPerm );
    else
        QR( A, t, d );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
}

template<typename F>
void Explicit( Matrix<F>& A, Matrix<F>& R, bool colPiv )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    Matrix<Int> pPerm;
    if( colPiv )
        QR( A, t, d, pPerm );
    else
        QR( A, t, d );
    Matrix<F> AT, AB;
    PartitionDown( A, AT, AB, Min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
}

template<typename F>
void Explicit( DistMatrix<F>& A, DistMatrix<F>& R, bool colPiv )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    DistMatrix<Int,VR,STAR> pPerm(g);
    if( colPiv )
        QR( A, t, d, pPerm );
    else
        QR( A, t, d );
    DistMatrix<F> AT(g), AB(g);
    PartitionDown( A, AT, AB, Min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
}

template<typename F>
void Explicit( Matrix<F>& A, Matrix<F>& R, Matrix<Int>& pPerm )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    QR( A, t, d, pPerm );
    Matrix<F> AT, AB;
    PartitionDown( A, AT, AB, Min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
} 

template<typename F,Dist UPerm>
void Explicit
( DistMatrix<F>& A, DistMatrix<F>& R, DistMatrix<Int,UPerm,STAR>& pPerm )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    QR( A, t, d, pPerm );
    DistMatrix<F> AT(g), AB(g);
    PartitionDown( A, AT, AB, Min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
    DiagonalScale( RIGHT, NORMAL, d, A );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_EXPLICIT_HPP
