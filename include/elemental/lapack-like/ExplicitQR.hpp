/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {

namespace internal {

template<typename Real>
inline void
ExplicitQRHelper( Matrix<Real>& A )
{
    QR( A );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Real> Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitQRHelper( DistMatrix<Real,MC,MR>& A )
{
    QR( A );

    // TODO: Replace this with an in-place expansion of Q
    const Grid& g = A.Grid();
    DistMatrix<Real,MC,MR> Q( g );
    Q.AlignWith( A );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitQRHelper( Matrix<Complex<Real> >& A )
{
    Matrix<Complex<Real> > t;
    QR( A, t );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Complex<Real> > Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitQRHelper( DistMatrix<Complex<Real>,MC,MR>& A )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>,MD,STAR> t( g );
    QR( A, t );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<Complex<Real>,MC,MR> Q( g );
    Q.AlignWith( A );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitQRHelper( Matrix<Real>& A, Matrix<Real>& R )
{
    QR( A );
    R = A;
    MakeTrapezoidal( LEFT, UPPER, 0, R );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Real> Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitQRHelper( DistMatrix<Real,MC,MR>& A, DistMatrix<Real,MC,MR>& R )
{
    QR( A );
    R = A;
    MakeTrapezoidal( LEFT, UPPER, 0, R );

    // TODO: Replace this with an in-place expansion of Q
    const Grid& g = A.Grid();
    DistMatrix<Real,MC,MR> Q( g );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitQRHelper( Matrix<Complex<Real> >& A, Matrix<Complex<Real> >& R )
{
    Matrix<Complex<Real> > t;
    QR( A, t );
    R = A;
    MakeTrapezoidal( LEFT, UPPER, 0, R );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Complex<Real> > Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitQRHelper
( DistMatrix<Complex<Real>,MC,MR>& A, DistMatrix<Complex<Real>,MC,MR>& R )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>,MD,STAR> t( g );
    QR( A, t );
    R = A;
    MakeTrapezoidal( LEFT, UPPER, 0, R );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<Complex<Real>,MC,MR> Q( g );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
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
ExplicitQR( DistMatrix<F,MC,MR>& A )
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
ExplicitQR( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& R )
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
