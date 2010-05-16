/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

template<typename R>
R
elemental::lapack::internal::LocalRowReflector
( DistMatrix<R,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::LocalRowReflector");
    if( x.Height() != 1 )
        throw "x must be a row vector.";
    if( x.GetGrid().MCRank() != x.ColAlignment() )
        throw "x is not aligned correctly.";
#endif
    if( x.Width() <= 1 )
        return (R)0;

    const Grid& grid = x.GetGrid();
    const int c = grid.Width();
    const int myCol = grid.MRRank();

    // For partitioning x
    DistMatrix<R,MC,MR> chi1(grid);
    DistMatrix<R,MC,MR> x2(grid);

    PartitionRight( x,  chi1, x2, 1 );

    R localNorm = blas::Nrm2( x2.LockedLocalMatrix() ); 
    R* localNorms = new R[c];
    AllGather( &localNorm, 1, localNorms, 1, grid.MRComm() );
    R norm = wrappers::blas::Nrm2( c, localNorms, 1 );

    R alpha;
    if( myCol == chi1.RowAlignment() )
        alpha = chi1.LocalEntry(0,0);
    Broadcast( &alpha, 1, chi1.RowAlignment(), grid.MRComm() );

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
        AllGather( &localNorm, 1, localNorms, 1, grid.MRComm() );
        norm = wrappers::blas::Nrm2( c, localNorms, 1 );
        if( alpha <= 0 )
            beta = wrappers::lapack::SafeNorm( alpha, norm );
        else
            beta = -wrappers::lapack::SafeNorm( alpha, norm );
    }
    delete localNorms;

    R tau = ( beta-alpha ) / beta;
    blas::Scal( static_cast<R>(1)/(alpha-beta), x2 );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    if( myCol == chi1.RowAlignment() )
        chi1.LocalEntry(0,0) = beta;
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

template float elemental::lapack::internal::LocalRowReflector
( DistMatrix<float,MC,MR>& x );

template double elemental::lapack::internal::LocalRowReflector
( DistMatrix<double,MC,MR>& x );

