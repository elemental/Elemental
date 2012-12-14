/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

namespace internal {

template<typename Real>
inline void
ExplicitQRHelper( Matrix<Real>& A )
{
    QR( A );
    ExpandPackedReflectors( LOWER, VERTICAL, 0, A );
}

template<typename Real>
inline void
ExplicitQRHelper( DistMatrix<Real>& A )
{
    QR( A );
    ExpandPackedReflectors( LOWER, VERTICAL, 0, A );
}

template<typename Real>
inline void
ExplicitQRHelper( Matrix<Complex<Real> >& A )
{
    Matrix<Complex<Real> > t;
    QR( A, t );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

template<typename Real>
inline void
ExplicitQRHelper( DistMatrix<Complex<Real> >& A )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>,MD,STAR> t( g );
    QR( A, t );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

template<typename Real>
inline void
ExplicitQRHelper( Matrix<Real>& A, Matrix<Real>& R )
{
    QR( A );
    Matrix<Real> AT,
                 AB;
    PartitionDown
    ( A, AT,
         AB, std::min(A.Height(),A.Width()) );
    R = AT;
    MakeTrapezoidal( LEFT, UPPER, 0, R );
    ExpandPackedReflectors( LOWER, VERTICAL, 0, A );
}

template<typename Real>
inline void
ExplicitQRHelper( DistMatrix<Real>& A, DistMatrix<Real>& R )
{
    const Grid& g = A.Grid();
    QR( A );
    DistMatrix<Real> AT(g),
                     AB(g);
    PartitionDown
    ( A, AT,
         AB, std::min(A.Height(),A.Width()) );
    R = AT;
    MakeTrapezoidal( LEFT, UPPER, 0, R );
    ExpandPackedReflectors( LOWER, VERTICAL, 0, A );
}

template<typename Real>
inline void
ExplicitQRHelper( Matrix<Complex<Real> >& A, Matrix<Complex<Real> >& R )
{
    Matrix<Complex<Real> > t;
    QR( A, t );
    Matrix<Complex<Real> > AT,
                           AB;
    PartitionDown
    ( A, AT,
         AB, std::min(A.Height(),A.Width()) );
    R = AT;
    MakeTrapezoidal( LEFT, UPPER, 0, R );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

template<typename Real>
inline void
ExplicitQRHelper
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
    MakeTrapezoidal( LEFT, UPPER, 0, R );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, A, t );
}

} // namespace internal

template<typename F> 
inline void
ExplicitQR( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("ExplicitQR");
#endif
    internal::ExplicitQRHelper( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
ExplicitQR( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("ExplicitQR");
#endif
    internal::ExplicitQRHelper( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
ExplicitQR( Matrix<F>& A, Matrix<F>& R )
{
#ifndef RELEASE
    PushCallStack("ExplicitQR");
#endif
    internal::ExplicitQRHelper( A, R );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
ExplicitQR( DistMatrix<F>& A, DistMatrix<F>& R )
{
#ifndef RELEASE
    PushCallStack("ExplicitQR");
#endif
    internal::ExplicitQRHelper( A, R );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
