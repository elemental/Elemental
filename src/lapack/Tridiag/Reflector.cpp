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

template<typename T>
T
elemental::lapack::internal::Reflector
( DistMatrix<T,MC,MR>& chi, DistMatrix<T,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::Reflector");
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

