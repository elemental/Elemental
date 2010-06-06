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
elemental::lapack::internal::LocalColReflector
( DistMatrix<R,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::LocalColReflector");
    if( x.Width() != 1 )
        throw "x must be a column vector.";
    if( x.GetGrid().MRRank() != x.RowAlignment() )
        throw "x is not aligned correctly.";
#endif
    if( x.Height() <= 1 )
        return (R)0;

    const Grid& grid = x.GetGrid();
    const int r = grid.Height();
    const int myRow = grid.MCRank();

    // For partitioning x
    DistMatrix<R,MC,MR> chi1(grid);
    DistMatrix<R,MC,MR> x2(grid);

    PartitionDown
    ( x, chi1,
         x2,  1 );

    vector<R> localNorms(r);
    R localNorm = blas::Nrm2( x2.LockedLocalMatrix() ); 
    AllGather( &localNorm, 1, &localNorms[0], 1, grid.MCComm() );
    R norm = wrappers::blas::Nrm2( r, &localNorms[0], 1 );

    R alpha;
    if( myRow == chi1.ColAlignment() )
        alpha = chi1.LocalEntry(0,0);
    Broadcast( &alpha, 1, chi1.ColAlignment(), grid.MCComm() );

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
            blas::Scal( invOfSafeMin, x2 );
            alpha *= invOfSafeMin;
            beta *= invOfSafeMin;
        } while( Abs( beta ) < safeMin );

        localNorm = blas::Nrm2( x2.LockedLocalMatrix() );
        AllGather( &localNorm, 1, &localNorms[0], 1, grid.MCComm() );
        norm = wrappers::blas::Nrm2( r, &localNorms[0], 1 );
        if( alpha <= 0 )
            beta = wrappers::lapack::SafeNorm( alpha, norm );
        else
            beta = -wrappers::lapack::SafeNorm( alpha, norm );
    }

    R tau = ( beta-alpha ) / beta;
    blas::Scal( static_cast<R>(1)/(alpha-beta), x2 );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    if( myRow == chi1.ColAlignment() )
        chi1.LocalEntry(0,0) = beta;
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
complex<R>
elemental::lapack::internal::LocalColReflector
( DistMatrix<complex<R>,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::LocalColReflector");
    if( x.Width() != 1 )
        throw "x must be a column vector.";
    if( x.GetGrid().MRRank() != x.RowAlignment() )
        throw "x is not aligned correctly.";
#endif
    typedef complex<R> C;

    if( x.Height() <= 1 )
        return (C)0;

    const Grid& grid = x.GetGrid();
    const int r = grid.Height();
    const int myRow = grid.MCRank();

    // For partitioning x
    DistMatrix<C,MC,MR> chi1(grid);
    DistMatrix<C,MC,MR> x2(grid);

    PartitionDown
    ( x, chi1,
         x2,   1 );

    vector<R> localNorms(r);
    R localNorm = blas::Nrm2( x2.LockedLocalMatrix() ); 
    AllGather( &localNorm, 1, &localNorms[0], 1, grid.MCComm() );
    R norm = wrappers::blas::Nrm2( r, &localNorms[0], 1 );

    C alpha;
    if( myRow == chi1.ColAlignment() )
        alpha = chi1.LocalEntry(0,0);
    Broadcast( &alpha, 1, chi1.ColAlignment(), grid.MCComm() );

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
            blas::Scal( (C)invOfSafeMin, x2 );
            alpha *= invOfSafeMin;
            beta *= invOfSafeMin;
        } while( Abs( beta ) < safeMin );

        localNorm = blas::Nrm2( x2.LockedLocalMatrix() );
        AllGather( &localNorm, 1, &localNorms[0], 1, grid.MCComm() );
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
    blas::Scal( static_cast<C>(1)/(alpha-beta), x2 );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    if( myRow == chi1.ColAlignment() )
        chi1.LocalEntry(0,0) = beta;
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}
#endif // WITHOUT_COMPLEX

template float elemental::lapack::internal::LocalColReflector
( DistMatrix<float,MC,MR>& x );

template double elemental::lapack::internal::LocalColReflector
( DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template scomplex elemental::lapack::internal::LocalColReflector
( DistMatrix<scomplex,MC,MR>& x );

template dcomplex elemental::lapack::internal::LocalColReflector
( DistMatrix<dcomplex,MC,MR>& x );
#endif 

