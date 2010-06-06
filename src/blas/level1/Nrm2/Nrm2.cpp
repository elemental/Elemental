/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
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
        throw logic_error( "x must be a vector." );
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
        throw logic_error( "x must be a vector." );
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

