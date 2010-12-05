/*
   Copyright (C) 1992-2008 The University of Tennessee
   All rights reserved.

   Copyright (c) 2009-2010, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is partially based upon the LAPACK 
   routines dlarfg.f and zlarfg.f.

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
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

template<typename R>
R
elemental::lapack::Reflector
( Matrix<R>& chi, Matrix<R>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::Reflector");
#endif
    if( x.Height() == 0 )
    {
        chi(0,0) *= (R)-1;
#ifndef RELEASE
        PopCallStack();
#endif
        return (R)2;
    }

    R norm = blas::Nrm2( x );
    R alpha = chi(0,0);

    R beta;
    if( alpha <= 0 )
        beta = wrappers::lapack::SafeNorm( alpha, norm );
    else
        beta = -wrappers::lapack::SafeNorm( alpha, norm );

    R safeMin = numeric_limits<R>::min() / numeric_limits<R>::epsilon();
    int count = 0;
    if( Abs( beta ) < safeMin )
    {
        R invOfSafeMin = static_cast<R>(1) / safeMin;
        do
        {
            ++count;
            blas::Scal( invOfSafeMin, x );
            alpha *= invOfSafeMin;
            beta *= invOfSafeMin;
        } while( Abs( beta ) < safeMin );

        norm = blas::Nrm2( x );
        if( alpha <= 0 )
            beta = wrappers::lapack::SafeNorm( alpha, norm );
        else
            beta = -wrappers::lapack::SafeNorm( alpha, norm );
    }

    R tau = ( beta - alpha ) / beta;
    blas::Scal( static_cast<R>(1)/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    chi(0,0) = beta;
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
complex<R>
elemental::lapack::Reflector
( Matrix< complex<R> >& chi, Matrix< complex<R> >& x )
{
#ifndef RELEASE
    PushCallStack("lapack::Reflector");
#endif
    typedef complex<R> C;

    R norm = blas::Nrm2( x );
    C alpha = chi(0,0);

    if( norm == 0 && imag(alpha) == (R)0 )
    {
        chi(0,0) *= (R)-1;
#ifndef RELEASE
        PopCallStack();
#endif
        return (C)2;
    }

    R beta;
    if( real(alpha) <= 0 )
        beta = wrappers::lapack::SafeNorm( real(alpha), imag(alpha), norm );
    else
        beta = -wrappers::lapack::SafeNorm( real(alpha), imag(alpha), norm );

    R safeMin = numeric_limits<R>::min() / numeric_limits<R>::epsilon();
    int count = 0;
    if( Abs( beta ) < safeMin )
    {
        R invOfSafeMin = static_cast<R>(1) / safeMin;
        do
        {
            ++count;
            blas::Scal( (C)invOfSafeMin, x );
            alpha *= invOfSafeMin;
            beta *= invOfSafeMin;
        } while( Abs( beta ) < safeMin );

        norm = blas::Nrm2( x );
        if( real(alpha) <= 0 )
            beta = wrappers::lapack::SafeNorm
                   ( real(alpha), imag(alpha), norm );
        else
            beta = -wrappers::lapack::SafeNorm
                    ( real(alpha), imag(alpha), norm );
    }

    C tau = C( (beta-real(alpha))/beta, -imag(alpha)/beta );
    blas::Scal( static_cast<C>(1)/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    chi(0,0) = beta;
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}
#endif // WITHOUT_COMPLEX

template<typename T>
T
elemental::lapack::Reflector
( DistMatrix<T,MC,MR>& chi, DistMatrix<T,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::Reflector");
    if( chi.Grid() != x.Grid() )
        throw logic_error( "chi and x must be distributed over the same grid" );
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw logic_error( "chi must be a scalar." );
    if( x.Height() != 1 && x.Width() != 1 )
        throw logic_error( "x must be a vector." );
#endif
    if( max( x.Height(), x.Width() ) == 0 )
        return (T)0;

    const Grid& g = x.Grid();
    T tau;
    if( x.Width() == 1 )
    {
        const bool thisIsMyColumn = ( g.MRRank() == x.RowAlignment() );
        if( thisIsMyColumn )
            tau = lapack::internal::ColReflector( chi, x );
        wrappers::mpi::Broadcast( &tau, 1, x.RowAlignment(), g.MRComm() );
    }
    else
    {
        const bool thisIsMyRow = ( g.MCRank() == x.ColAlignment() );
        if( thisIsMyRow )
            tau = lapack::internal::RowReflector( chi, x );
        wrappers::mpi::Broadcast( &tau, 1, x.ColAlignment(), g.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

template float elemental::lapack::Reflector
( Matrix<float>& chi, Matrix<float>& x );

template float elemental::lapack::Reflector
( DistMatrix<float,MC,MR>& chi, DistMatrix<float,MC,MR>& x );

template double elemental::lapack::Reflector
( Matrix<double>& chi, Matrix<double>& x );

template double elemental::lapack::Reflector
( DistMatrix<double,MC,MR>& chi, DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template scomplex elemental::lapack::Reflector
( Matrix<scomplex>& chi, Matrix<scomplex>& x );

template scomplex elemental::lapack::Reflector
( DistMatrix<scomplex,MC,MR>& chi, DistMatrix<scomplex,MC,MR>& x );

template dcomplex elemental::lapack::Reflector
( Matrix<dcomplex>& chi, Matrix<dcomplex>& x );

template dcomplex elemental::lapack::Reflector
( DistMatrix<dcomplex,MC,MR>& chi, DistMatrix<dcomplex,MC,MR>& x );
#endif

