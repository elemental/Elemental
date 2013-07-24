/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_QR_EXPLICIT_HPP
#define ELEM_LAPACK_QR_EXPLICIT_HPP

#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/lapack-like/ExpandPackedReflectors.hpp"
#include "elemental/lapack-like/QR.hpp"

namespace elem {
namespace qr {

template<typename F>
inline void
Explicit( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("qr::Explicit");
#endif
    Matrix<F> t;
    QR( A, t );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

template<typename F>
inline void
Explicit( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("qr::Explicit");
#endif
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t( g );
    QR( A, t );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

template<typename F>
inline void
Explicit( Matrix<F>& A, Matrix<F>& R )
{
#ifndef RELEASE
    CallStackEntry cse("qr::Explicit");
#endif
    Matrix<F> t;
    QR( A, t );
    Matrix<F> AT,
              AB;
    PartitionDown
    ( A, AT,
         AB, std::min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

template<typename F>
inline void
Explicit( DistMatrix<F>& A, DistMatrix<F>& R )
{
#ifndef RELEASE
    CallStackEntry cse("qr::Explicit");
#endif
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t( g );
    QR( A, t );
    DistMatrix<F> AT(g),
                  AB(g);
    PartitionDown
    ( A, AT,
         AB, std::min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

} // namespace qr
} // namespace elem

#endif // ifndef ELEM_LAPACK_QR_EXPLICIT_HPP
