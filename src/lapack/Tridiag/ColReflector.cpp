/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
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
    if( x.Height() == 0 )
        return (R)0;

    const Grid& g = x.GetGrid();
    const int r = g.Height();
    const int myRow = g.MCRank();

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
        return (C)0;

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

