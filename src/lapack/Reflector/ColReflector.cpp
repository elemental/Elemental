/*
   Copyright (c) 1992-2008 The University of Tennessee.
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
using namespace elemental::wrappers::mpi;

template<typename R>
R
elemental::lapack::internal::ColReflector
( DistMatrix<R,MC,MR>& chi, DistMatrix<R,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::ColReflector");
    if( chi.GetGrid() != x.GetGrid() )
        throw logic_error( "chi and x must be distributed over the same grid" );
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw logic_error( "chi must be a scalar." );
    if( x.Width() != 1 )
        throw logic_error( "x must be a column vector." );
    if( chi.GetGrid().MRRank() != chi.RowAlignment() )
        throw logic_error( "Reflecting with incorrect column of processes." );
    if( x.GetGrid().MRRank() != x.RowAlignment() )
        throw logic_error( "Reflecting with incorrect column of processes." );
#endif

    const Grid& g = x.GetGrid();
    const int r = g.Height();
    const int myRow = g.MCRank();

    if( x.Height() == 0 )
    {
        if( myRow == chi.ColAlignment() )
            chi.LocalEntry(0,0) *= (R)-1;
#ifndef RELEASE
        PopCallStack();
#endif
        return (R)2;
    }

    vector<R> localNorms(r);
    R localNorm = blas::Nrm2( x.LockedLocalMatrix() ); 
    AllGather( &localNorm, 1, &localNorms[0], 1, g.MCComm() );
    R norm = wrappers::blas::Nrm2( r, &localNorms[0], 1 );

    R alpha;
    if( myRow == chi.ColAlignment() )
        alpha = chi.LocalEntry(0,0);
    Broadcast( &alpha, 1, chi.ColAlignment(), g.MCComm() );

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

        localNorm = blas::Nrm2( x.LockedLocalMatrix() );
        AllGather( &localNorm, 1, &localNorms[0], 1, g.MCComm() );
        norm = wrappers::blas::Nrm2( r, &localNorms[0], 1 );
        if( alpha <= 0 )
            beta = wrappers::lapack::SafeNorm( alpha, norm );
        else
            beta = -wrappers::lapack::SafeNorm( alpha, norm );
    }

    R tau = ( beta-alpha ) / beta;
    blas::Scal( static_cast<R>(1)/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    if( myRow == chi.ColAlignment() )
        chi.LocalEntry(0,0) = beta;
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
complex<R>
elemental::lapack::internal::ColReflector
( DistMatrix<complex<R>,MC,MR>& chi, DistMatrix<complex<R>,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::ColReflector");
    if( chi.GetGrid() != x.GetGrid() )
        throw logic_error( "chi and x must be distributed over the same grid" );
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw logic_error( "chi must be a scalar." );
    if( x.Width() != 1 )
        throw logic_error( "x must be a column vector." );
    if( chi.GetGrid().MRRank() != chi.RowAlignment() )
        throw logic_error( "Reflecting with incorrect column of processes." );
    if( x.GetGrid().MRRank() != x.RowAlignment() )
        throw logic_error( "Reflecting with incorrect column of processes." );
#endif
    typedef complex<R> C;

    const Grid& g = x.GetGrid();
    const int r = g.Height();
    const int myRow = g.MCRank();

    vector<R> localNorms(r);
    R localNorm = blas::Nrm2( x.LockedLocalMatrix() ); 
    AllGather( &localNorm, 1, &localNorms[0], 1, g.MCComm() );
    R norm = wrappers::blas::Nrm2( r, &localNorms[0], 1 );

    C alpha;
    if( myRow == chi.ColAlignment() )
        alpha = chi.LocalEntry(0,0);
    Broadcast( &alpha, 1, chi.ColAlignment(), g.MCComm() );

    if( norm == (R)0 && imag(alpha) == (R)0 )
    {
        if( myRow == chi.ColAlignment() )
            chi.LocalEntry(0,0) *= (R)-1;
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

        localNorm = blas::Nrm2( x.LockedLocalMatrix() );
        AllGather( &localNorm, 1, &localNorms[0], 1, g.MCComm() );
        norm = wrappers::blas::Nrm2( r, &localNorms[0], 1 );
        if( real(alpha) <= 0 )
        {
            beta = wrappers::lapack::SafeNorm
                   ( real(alpha), imag(alpha), norm );
        }
        else
        {
            beta = -wrappers::lapack::SafeNorm
                   ( real(alpha), imag(alpha), norm );
        }
    }

    C tau = C( (beta-real(alpha))/beta, -imag(alpha)/beta );
    blas::Scal( static_cast<C>(1)/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    if( myRow == chi.ColAlignment() )
        chi.LocalEntry(0,0) = beta;
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}
#endif // WITHOUT_COMPLEX

template float elemental::lapack::internal::ColReflector
( DistMatrix<float,MC,MR>& chi, DistMatrix<float,MC,MR>& x );

template double elemental::lapack::internal::ColReflector
( DistMatrix<double,MC,MR>& chi, DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template scomplex elemental::lapack::internal::ColReflector
( DistMatrix<scomplex,MC,MR>& chi, DistMatrix<scomplex,MC,MR>& x );

template dcomplex elemental::lapack::internal::ColReflector
( DistMatrix<dcomplex,MC,MR>& chi, DistMatrix<dcomplex,MC,MR>& x );
#endif 

