/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QR_EXPLICIT_HPP
#define ELEM_QR_EXPLICIT_HPP

#include ELEM_MAKETRIANGULAR_INC
#include ELEM_EXPANDPACKEDREFLECTORS_INC
#include ELEM_QR_INC

namespace elem {
namespace qr {

template<typename F>
inline void
Explicit( Matrix<F>& A, bool colPiv=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    Matrix<F> t;
    if( colPiv )
    {
        Matrix<Int> p;
        QR( A, t, p );
    }
    else
    {
        QR( A, t );
    }
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
}

template<typename F>
inline void
Explicit( DistMatrix<F>& A, bool colPiv=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    if( colPiv )
    {
        DistMatrix<Int,VR,STAR> p(g);
        QR( A, t, p );
    }
    else
    {
        QR( A, t );
    }
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
}

template<typename F>
inline void
Explicit( Matrix<F>& A, Matrix<F>& R, bool colPiv=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    Matrix<F> t;
    if( colPiv )
    {
        Matrix<Int> p;
        QR( A, t, p );
    }
    else
    {
        QR( A, t );
    }
    Matrix<F> AT,
              AB;
    PartitionDown
    ( A, AT,
         AB, Min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
}

template<typename F>
inline void
Explicit( DistMatrix<F>& A, DistMatrix<F>& R, bool colPiv=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Explicit"))
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    if( colPiv )
    {
        DistMatrix<Int,VR,STAR> p(g);
        QR( A, t, p );
    }
    else
    {
        QR( A, t );
    }
    DistMatrix<F> AT(g),
                  AB(g);
    PartitionDown
    ( A, AT,
         AB, Min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, A, t );
}

} // namespace qr
} // namespace elem

#endif // ifndef ELEM_QR_EXPLICIT_HPP
