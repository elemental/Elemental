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
( DistMatrix<R,MC,MR>& chi, DistMatrix<R,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::Reflector");
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw logic_error( "chi must be a scalar." );
    if( x.Height() != 1 && x.Width() != 1 )
        throw logic_error( "x must be a vector." );
#endif
    if( max( x.Height(), x.Width() ) == 0 )
        return (R)0;

    R norm = blas::Nrm2( x );
    R alpha = chi.Get(0,0);
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
    R tau = ( beta-alpha ) / beta;
    blas::Scal( static_cast<R>(1)/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    chi.Set(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
complex<R>
elemental::lapack::internal::Reflector
( DistMatrix<complex<R>,MC,MR>& chi, DistMatrix<complex<R>,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::Reflector");
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw logic_error( "chi must be a scalar." );
    if( x.Height() != 1 && x.Width() != 1 )
        throw logic_error( "x must be a vector." );
#endif
    typedef complex<R> C;

    if( max( x.Height(), x.Width() ) == 0 )
        return (C)0;

    R norm = blas::Nrm2( x );
    C alpha = chi.Get(0,0);

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

        norm = blas::Nrm2( x );
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
    chi.Set(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}
#endif // WITHOUT_COMPLEX

template float elemental::lapack::internal::Reflector
( DistMatrix<float,MC,MR>& chi, DistMatrix<float,MC,MR>& x );

template double elemental::lapack::internal::Reflector
( DistMatrix<double,MC,MR>& chi, DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template scomplex elemental::lapack::internal::Reflector
( DistMatrix<scomplex,MC,MR>& chi, DistMatrix<scomplex,MC,MR>& x );

template dcomplex elemental::lapack::internal::Reflector
( DistMatrix<dcomplex,MC,MR>& chi, DistMatrix<dcomplex,MC,MR>& x );
#endif

