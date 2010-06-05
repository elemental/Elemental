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

template<typename R>
R
elemental::lapack::internal::Reflector
( DistMatrix<R,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::Reflector");
    if( x.Height() != 1 && x.Width() != 1 )
        throw "x must be a vector.";
#endif
    if( max( x.Height(), x.Width() ) <= 1 )
        return (R)0;

    const Grid& grid = x.GetGrid();

    // For partitioning x
    DistMatrix<R,MC,MR> chi1(grid);
    DistMatrix<R,MC,MR> x2(grid);

    if( x.Width() == 1 )
    {
        PartitionDown
        ( x, chi1,
             x2,   1 );
    }
    else
    {
        PartitionRight( x, chi1, x2, 1 );
    }

    R x2_Norm = blas::Nrm2( x2 );
    R alpha = chi1.Get(0,0);
    R beta;
    if( alpha <= 0 )
        beta = wrappers::lapack::SafeNorm( alpha, x2_Norm );
    else
        beta = -wrappers::lapack::SafeNorm( alpha, x2_Norm );

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

        x2_Norm = blas::Nrm2( x2 );
        if( alpha <= 0 )
            beta = wrappers::lapack::SafeNorm( alpha, x2_Norm );
        else
            beta = -wrappers::lapack::SafeNorm( alpha, x2_Norm );
    }
    R tau = ( beta-alpha ) / beta;
    blas::Scal( static_cast<R>(1)/(alpha-beta), x2 );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    chi1.Set(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
complex<R>
elemental::lapack::internal::Reflector
( DistMatrix<complex<R>,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::Reflector");
    if( x.Height() != 1 && x.Width() != 1 )
        throw "x must be a vector.";
#endif
    typedef complex<R> C;

    if( max( x.Height(), x.Width() ) <= 1 )
        return (C)0;

    const Grid& grid = x.GetGrid();

    // For partitioning x
    DistMatrix<C,MC,MR> chi1(grid);
    DistMatrix<C,MC,MR> x2(grid);

    if( x.Width() == 1 )
    {
        PartitionDown
        ( x, chi1,
             x2,   1 );
    }
    else
    {
        PartitionRight( x, chi1, x2, 1 );
    }

    R x2_Norm = blas::Nrm2( x2 );
    C alpha = chi1.Get(0,0);
    R beta;
    if( real(alpha) <= 0 )
        beta = wrappers::lapack::SafeNorm( real(alpha), imag(alpha), x2_Norm );
    else
        beta = -wrappers::lapack::SafeNorm( real(alpha), imag(alpha), x2_Norm );

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

        x2_Norm = blas::Nrm2( x2 );
        if( real(alpha) <= 0 )
        {
            beta = wrappers::lapack::SafeNorm
                   ( real(alpha), imag(alpha), x2_Norm );
        }
        else
        {
            beta = -wrappers::lapack::SafeNorm
                    ( real(alpha), imag(alpha), x2_Norm );
        }
    }
    C tau = C( (beta-real(alpha))/beta, -imag(alpha)/beta );
    blas::Scal( static_cast<C>(1)/(alpha-beta), x2 );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    chi1.Set(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}
#endif // WITHOUT_COMPLEX

template float elemental::lapack::internal::Reflector
( DistMatrix<float,MC,MR>& x );

template double elemental::lapack::internal::Reflector
( DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template scomplex elemental::lapack::internal::Reflector
( DistMatrix<scomplex,MC,MR>& x );

template dcomplex elemental::lapack::internal::Reflector
( DistMatrix<dcomplex,MC,MR>& x );
#endif

