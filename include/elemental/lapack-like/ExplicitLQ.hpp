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
ExplicitLQHelper( Matrix<Real>& A )
{
    LQ( A );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Real> Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( DistMatrix<Real,MC,MR>& A )
{
    LQ( A );

    // TODO: Replace this with an in-place expansion of Q
    const Grid& g = A.Grid();
    DistMatrix<Real,MC,MR> Q( g );
    Q.AlignWith( A );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( Matrix<Complex<Real> >& A )
{
    Matrix<Complex<Real> > t;
    LQ( A, t );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Complex<Real> > Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( DistMatrix<Complex<Real>,MC,MR>& A )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>,MD,STAR> t( g );
    LQ( A, t );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<Complex<Real>,MC,MR> Q( g );
    Q.AlignWith( A );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( Matrix<Real>& L, Matrix<Real>& A )
{
    LQ( A );
    L = A;
    MakeTrapezoidal( LEFT, LOWER, 0, L );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Real> Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( DistMatrix<Real,MC,MR>& L, DistMatrix<Real,MC,MR>& A )
{
    LQ( A );
    L = A;
    MakeTrapezoidal( LEFT, LOWER, 0, L );

    // TODO: Replace this with an in-place expansion of Q
    const Grid& g = A.Grid();
    DistMatrix<Real,MC,MR> Q( g );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( Matrix<Complex<Real> >& L, Matrix<Complex<Real> >& A )
{
    Matrix<Complex<Real> > t;
    LQ( A, t );
    L = A;
    MakeTrapezoidal( LEFT, LOWER, 0, L );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Complex<Real> > Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper
( DistMatrix<Complex<Real>,MC,MR>& L, DistMatrix<Complex<Real>,MC,MR>& A )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>,MD,STAR> t( g );
    LQ( A, t );
    L = A;
    MakeTrapezoidal( LEFT, LOWER, 0, L );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<Complex<Real>,MC,MR> Q( g );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

} // namespace internal

template<typename F> 
inline void
ExplicitLQ( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("ExplicitLQ");
#endif
    internal::ExplicitLQHelper( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
ExplicitLQ( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("ExplicitLQ");
#endif
    internal::ExplicitLQHelper( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
ExplicitLQ( Matrix<F>& L, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("ExplicitLQ");
#endif
    internal::ExplicitLQHelper( L, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
ExplicitLQ( DistMatrix<F,MC,MR>& L, DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("ExplicitLQ");
#endif
    internal::ExplicitLQHelper( L, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
