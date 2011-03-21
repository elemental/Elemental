/*
   Copyright (c) 1992-2008 The University of Tennessee. 
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
#include "elemental/advanced_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::imports::mpi;

template<typename R> // representation of a real number
R
elemental::advanced::internal::RowReflector
( DistMatrix<R,MC,MR>& chi, DistMatrix<R,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::RowReflector");
    if( chi.Grid() != x.Grid() )
        throw logic_error( "chi and x must be distributed over the same grid" );
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw logic_error( "chi must be a scalar." );
    if( x.Height() != 1 )
        throw logic_error( "x must be a row vector." );
    if( chi.Grid().MCRank() != chi.ColAlignment() )
        throw logic_error( "Reflecting with incorrect row of processes." );
    if( x.Grid().MCRank() != x.ColAlignment() )
        throw logic_error( "Reflecting with incorrect row of processes." );
#endif
    
    const Grid& g = x.Grid();
    const int c = g.Width();
    const int myCol = g.MRRank();

    if( x.Width() == 0 )
    {
        if( myCol == chi.RowAlignment() )
            chi.SetLocalEntry(0,0,-chi.GetLocalEntry(0,0));
#ifndef RELEASE
        PopCallStack();
#endif
        return (R)2;
    }

    vector<R> localNorms(c);
    R localNorm = basic::Nrm2( x.LockedLocalMatrix() ); 
    AllGather( &localNorm, 1, &localNorms[0], 1, g.MRComm() );
    R norm = imports::blas::Nrm2( c, &localNorms[0], 1 );

    R alpha;
    if( myCol == chi.RowAlignment() )
        alpha = chi.GetLocalEntry(0,0);
    Broadcast( &alpha, 1, chi.RowAlignment(), g.MRComm() );

    R beta;
    if( alpha <= 0 )
        beta = imports::lapack::SafeNorm( alpha, norm );
    else
        beta = -imports::lapack::SafeNorm( alpha, norm );

    const R safeMin = imports::lapack::MachineSafeMin<R>() /
                      imports::lapack::MachineEpsilon<R>();
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

        localNorm = basic::Nrm2( x.LockedLocalMatrix() );
        AllGather( &localNorm, 1, &localNorms[0], 1, g.MRComm() );
        norm = imports::blas::Nrm2( c, &localNorms[0], 1 );
        if( alpha <= 0 )
            beta = imports::lapack::SafeNorm( alpha, norm );
        else
            beta = -imports::lapack::SafeNorm( alpha, norm );
    }

    R tau = ( beta-alpha ) / beta;
    basic::Scal( static_cast<R>(1)/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    if( myCol == chi.RowAlignment() )
        chi.SetLocalEntry(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
complex<R>
elemental::advanced::internal::RowReflector
( DistMatrix<complex<R>,MC,MR>& chi, DistMatrix<complex<R>,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::RowReflector");    
    if( chi.Grid() != x.Grid() )
        throw logic_error( "chi and x must be distributed over the same grid" );
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw logic_error( "chi must be a scalar." );
    if( x.Height() != 1 )
        throw logic_error( "x must be a row vector." );
    if( chi.Grid().MCRank() != chi.ColAlignment() )
        throw logic_error( "Reflecting with incorrect row of processes." );
    if( x.Grid().MCRank() != x.ColAlignment() )
        throw logic_error( "Reflecting with incorrect row of processes." );
#endif
    typedef complex<R> C;

    const Grid& g = x.Grid();
    const int c = g.Width();
    const int myCol = g.MRRank();

    vector<R> localNorms(c);
    R localNorm = basic::Nrm2( x.LockedLocalMatrix() ); 
    AllGather( &localNorm, 1, &localNorms[0], 1, g.MRComm() );
    R norm = imports::blas::Nrm2( c, &localNorms[0], 1 );

    C alpha;
    if( myCol == chi.RowAlignment() )
        alpha = chi.GetLocalEntry(0,0);
    Broadcast( &alpha, 1, chi.RowAlignment(), g.MRComm() );

    if( norm == (R)0 && imag(alpha) == (R)0 )
    {
        if( myCol == chi.RowAlignment() )
            chi.SetLocalEntry(0,0,-chi.GetLocalEntry(0,0));
#ifndef RELEASE
        PopCallStack();
#endif
        return (C)2;
    }

    R beta;
    if( real(alpha) <= 0 )
        beta = imports::lapack::SafeNorm( real(alpha), imag(alpha), norm );
    else
        beta = -imports::lapack::SafeNorm( real(alpha), imag(alpha), norm );

    const R safeMin = imports::lapack::MachineSafeMin<R>() /
                      imports::lapack::MachineEpsilon<R>();
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

        localNorm = basic::Nrm2( x.LockedLocalMatrix() );
        AllGather( &localNorm, 1, &localNorms[0], 1, g.MRComm() );
        norm = imports::blas::Nrm2( c, &localNorms[0], 1 );
        if( real(alpha) <= 0 )
        {
            beta = imports::lapack::SafeNorm
                   ( real(alpha), imag(alpha), norm );
        }
        else
        {
            beta = -imports::lapack::SafeNorm
                    ( real(alpha), imag(alpha), norm );
        }
    }

    C tau = C( (beta-real(alpha))/beta, -imag(alpha)/beta );
    basic::Scal( static_cast<C>(1)/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    if( myCol == chi.RowAlignment() )
        chi.SetLocalEntry(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}
#endif // WITHOUT_COMPLEX

template float elemental::advanced::internal::RowReflector
( DistMatrix<float,MC,MR>& chi, DistMatrix<float,MC,MR>& x );

template double elemental::advanced::internal::RowReflector
( DistMatrix<double,MC,MR>& chi, DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template scomplex elemental::advanced::internal::RowReflector
( DistMatrix<scomplex,MC,MR>& chi, DistMatrix<scomplex,MC,MR>& x );

template dcomplex elemental::advanced::internal::RowReflector
( DistMatrix<dcomplex,MC,MR>& chi, DistMatrix<dcomplex,MC,MR>& x );
#endif

