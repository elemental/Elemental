/*
   Copyright (C) 1992-2008 The University of Tennessee
   All rights reserved.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.

   This file is partially based upon the LAPACK routines dlarfg.f and zlarfg.f.
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
    if( chi.GetGrid() != x.GetGrid() )
        throw logic_error( "chi and x must be distributed over the same grid" );
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw logic_error( "chi must be a scalar." );
    if( x.Height() != 1 && x.Width() != 1 )
        throw logic_error( "x must be a vector." );
#endif
    if( max( x.Height(), x.Width() ) == 0 )
        return (T)0;

    const Grid& g = x.GetGrid();
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

