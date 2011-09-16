/*
   Copyright (C) 1992-2008 The University of Tennessee
   All rights reserved.

   Copyright (c) 2009-2011, Jack Poulson
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

#include "./Reflector/ColReflector.hpp"
#include "./Reflector/RowReflector.hpp"

template<typename R> // representation of a real number
inline R
elemental::advanced::Reflector
( Matrix<R>& chi, Matrix<R>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::Reflector");
#endif
    if( x.Height() == 0 )
    {
        chi.Set(0,0,-chi.Get(0,0));
#ifndef RELEASE
        PopCallStack();
#endif
        return (R)2;
    }

    R norm = basic::Nrm2( x );
    R alpha = chi.Get(0,0);

    R beta;
    if( alpha <= 0 )
        beta = lapack::SafeNorm( alpha, norm );
    else
        beta = -lapack::SafeNorm( alpha, norm );

    const R safeMin = lapack::MachineSafeMin<R>() /
                      lapack::MachineEpsilon<R>();
    int count = 0;
    if( Abs( beta ) < safeMin )
    {
        R invOfSafeMin = static_cast<R>(1) / safeMin;
        do
        {
            ++count;
            basic::Scal( invOfSafeMin, x );
            alpha *= invOfSafeMin;
            beta *= invOfSafeMin;
        } while( Abs( beta ) < safeMin );

        norm = basic::Nrm2( x );
        if( alpha <= 0 )
            beta = lapack::SafeNorm( alpha, norm );
        else
            beta = -lapack::SafeNorm( alpha, norm );
    }

    R tau = ( beta - alpha ) / beta;
    basic::Scal( static_cast<R>(1)/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    chi.Set(0,0,beta);
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
inline std::complex<R>
elemental::advanced::Reflector
( Matrix<std::complex<R> >& chi, Matrix<std::complex<R> >& x )
{
#ifndef RELEASE
    PushCallStack("lapack::Reflector");
#endif
    typedef std::complex<R> C;

    R norm = basic::Nrm2( x );
    C alpha = chi.Get(0,0);

    if( norm == 0 && imag(alpha) == (R)0 )
    {
        chi.Set(0,0,-chi.Get(0,0));
#ifndef RELEASE
        PopCallStack();
#endif
        return (C)2;
    }

    R beta;
    if( real(alpha) <= 0 )
        beta = lapack::SafeNorm( real(alpha), imag(alpha), norm );
    else
        beta = -lapack::SafeNorm( real(alpha), imag(alpha), norm );

    const R safeMin = lapack::MachineSafeMin<R>() /
                      lapack::MachineEpsilon<R>();
    int count = 0;
    if( Abs( beta ) < safeMin )
    {
        R invOfSafeMin = static_cast<R>(1) / safeMin;
        do
        {
            ++count;
            basic::Scal( (C)invOfSafeMin, x );
            alpha *= invOfSafeMin;
            beta *= invOfSafeMin;
        } while( Abs( beta ) < safeMin );

        norm = basic::Nrm2( x );
        if( real(alpha) <= 0 )
            beta = lapack::SafeNorm
                   ( real(alpha), imag(alpha), norm );
        else
            beta = -lapack::SafeNorm
                    ( real(alpha), imag(alpha), norm );
    }

    C tau = C( (beta-real(alpha))/beta, -imag(alpha)/beta );
    basic::Scal( static_cast<C>(1)/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    chi.Set(0,0,beta);
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}
#endif // WITHOUT_COMPLEX

template<typename F> // represents a real or complex number
inline F
elemental::advanced::Reflector
( DistMatrix<F,MC,MR>& chi, DistMatrix<F,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::Reflector");
    if( chi.Grid() != x.Grid() )
        throw std::logic_error
        ("chi and x must be distributed over the same grid");
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw std::logic_error("chi must be a scalar");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("x must be a vector");
#endif
    if( std::max( x.Height(), x.Width() ) == 0 )
        return (F)0;

    const Grid& g = x.Grid();
    F tau;
    if( x.Width() == 1 && x.RowAlignment() == chi.RowAlignment() )
    {
        const bool thisIsMyColumn = ( g.MRRank() == x.RowAlignment() );
        if( thisIsMyColumn )
            tau = advanced::internal::ColReflector( chi, x );
        mpi::Broadcast( &tau, 1, x.RowAlignment(), g.MRComm() );
    }
    else
    {
        const bool thisIsMyRow = ( g.MCRank() == x.ColAlignment() );
        if( thisIsMyRow )
            tau = advanced::internal::RowReflector( chi, x );
        mpi::Broadcast( &tau, 1, x.ColAlignment(), g.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}
