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
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

template<typename R>
R
elemental::blas::Nrm2
( const DistMatrix<R,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("blas::Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw "x must be a vector.";
#endif
    R norm;
    const Grid& grid = x.GetGrid();

    if( x.Width() == 1 )
    {
        const int ownerCol = x.RowAlignment();
        if( grid.MRRank() == ownerCol )
        {
            R localNorm = Nrm2( x.LockedLocalMatrix() ); 
            
            const int r = grid.Height();
            vector<R> localNorms(r);
            R* localNormsPtr = &localNorms[0];
            AllGather( &localNorm, 1, localNormsPtr, 1, grid.MCComm() );
            norm = wrappers::blas::Nrm2( r, localNormsPtr, 1 );
        }
        Broadcast( &norm, 1, ownerCol, grid.MRComm() );
    }
    else
    {
        const int ownerRow = x.ColAlignment();
        if( grid.MCRank() == ownerRow )
        {
            R localNorm = Nrm2( x.LockedLocalMatrix() );

            const int c = grid.Width();
            vector<R> localNorms(c);
            R* localNormsPtr = &localNorms[0];
            AllGather( &localNorm, 1, localNormsPtr, 1, grid.MRComm() );
            norm = wrappers::blas::Nrm2( c, localNormsPtr, 1 );
        }
        Broadcast( &norm, 1, ownerRow, grid.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
R
elemental::blas::Nrm2
( const DistMatrix< complex<R>, MC, MR >& x )
{
#ifndef RELEASE
    PushCallStack("blas::Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw "x must be a vector.";
#endif
    R norm;
    const Grid& grid = x.GetGrid();

    if( x.Width() == 1 )
    {
        const int ownerCol = x.RowAlignment();
        if( grid.MRRank() == ownerCol )
        {
            R localNorm = Nrm2( x.LockedLocalMatrix() ); 
            
            const int r = grid.Height();
            vector<R> localNorms(r);
            R* localNormsPtr = &localNorms[0];
            AllGather( &localNorm, 1, localNormsPtr, 1, grid.MCComm() );
            norm = wrappers::blas::Nrm2( r, localNormsPtr, 1 );
        }
        Broadcast( &norm, 1, ownerCol, grid.MRComm() );
    }
    else
    {
        const int ownerRow = x.ColAlignment();
        if( grid.MCRank() == ownerRow )
        {
            R localNorm = Nrm2( x.LockedLocalMatrix() );

            const int c = grid.Width();
            vector<R> localNorms(c);
            R* localNormsPtr = &localNorms[0];
            AllGather( &localNorm, 1, localNormsPtr, 1, grid.MRComm() );
            norm = wrappers::blas::Nrm2( c, localNormsPtr, 1 );
        }
        Broadcast( &norm, 1, ownerRow, grid.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}
#endif

template float elemental::blas::Nrm2
( const DistMatrix<float,MC,MR>& x );

template double elemental::blas::Nrm2
( const DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template float elemental::blas::Nrm2
( const DistMatrix< complex<float>, MC, MR >& x );

template double elemental::blas::Nrm2
( const DistMatrix< complex<double>, MC, MR >& x );
#endif

