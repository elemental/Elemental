/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_QR_EXPLICIT_HPP
#define LAPACK_QR_EXPLICIT_HPP

#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/lapack-like/ExpandPackedReflectors.hpp"
#include "elemental/lapack-like/QR.hpp"

namespace elem {
namespace qr {

template<typename Real>
inline void
ExplicitHelper( Matrix<Real>& A )
{
    QR( A );
    ExpandPackedReflectors( LOWER, VERTICAL, 0, A );
}

template<typename Real>
inline void
ExplicitHelper( DistMatrix<Real>& A )
{
    QR( A );
    ExpandPackedReflectors( LOWER, VERTICAL, 0, A );
}

template<typename Real>
inline void
ExplicitHelper( Matrix<Complex<Real> >& A )
{
    Matrix<Complex<Real> > t;
    QR( A, t );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

template<typename Real>
inline void
ExplicitHelper( DistMatrix<Complex<Real> >& A )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>,MD,STAR> t( g );
    QR( A, t );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

template<typename Real>
inline void
ExplicitHelper( Matrix<Real>& A, Matrix<Real>& R )
{
    QR( A );
    Matrix<Real> AT,
                 AB;
    PartitionDown
    ( A, AT,
         AB, std::min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, 0, A );
}

template<typename Real>
inline void
ExplicitHelper( DistMatrix<Real>& A, DistMatrix<Real>& R )
{
    const Grid& g = A.Grid();
    QR( A );
    DistMatrix<Real> AT(g),
                     AB(g);
    PartitionDown
    ( A, AT,
         AB, std::min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, 0, A );
}

template<typename Real>
inline void
ExplicitHelper( Matrix<Complex<Real> >& A, Matrix<Complex<Real> >& R )
{
    Matrix<Complex<Real> > t;
    QR( A, t );
    Matrix<Complex<Real> > AT,
                           AB;
    PartitionDown
    ( A, AT,
         AB, std::min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

template<typename Real>
inline void
ExplicitHelper
( DistMatrix<Complex<Real> >& A, DistMatrix<Complex<Real> >& R )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>,MD,STAR> t( g );
    QR( A, t );
    DistMatrix<Complex<Real> > AT(g),
                               AB(g);
    PartitionDown
    ( A, AT,
         AB, std::min(A.Height(),A.Width()) );
    R = AT;
    MakeTriangular( UPPER, R );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

template<typename F> 
inline void
Explicit( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("qr::Explicit");
#endif
    ExplicitHelper( A );
}

template<typename F> 
inline void
Explicit( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("qr::Explicit");
#endif
    ExplicitHelper( A );
}

template<typename F> 
inline void
Explicit( Matrix<F>& A, Matrix<F>& R )
{
#ifndef RELEASE
    CallStackEntry entry("qr::Explicit");
#endif
    ExplicitHelper( A, R );
}

template<typename F> 
inline void
Explicit( DistMatrix<F>& A, DistMatrix<F>& R )
{
#ifndef RELEASE
    CallStackEntry entry("qr::Explicit");
#endif
    ExplicitHelper( A, R );
}

} // namespace qr
} // namespace elem

#endif // ifndef LAPACK_QR_EXPLICIT_HPP
