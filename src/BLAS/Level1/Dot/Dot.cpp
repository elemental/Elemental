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
#include "ElementalBLAS_Internal.h"
using namespace std;
using namespace Elemental;
using namespace Elemental::wrappers::MPI;

/* 
   C++ does not (currently) allow for partial function template specialization,
   so implementing Dot will be unnecessarily obfuscated. Sorry.
  
   The compromise is to have the user-level routine, 
  
     template<typename T, Distribution U, Distribution V,
                          Distribution W, Distribution Z >
     T
     BLAS::Dot( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

   simply route to overloaded pseudo-partial-specializations within the 
   BLAS::Internal namespace. For example,

     template<typename T, Distribution U, Distribution V>
     T
     BLAS::Internal::Dot
     ( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,MR>& y );
*/
template<typename T, Distribution U, Distribution V,
                     Distribution W, Distribution Z >
inline T
Elemental::BLAS::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Dot");
#endif
    T dotProduct = BLAS::Internal::Dot( x, y );
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename T, Distribution U, Distribution V>
inline T
Elemental::BLAS::Internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Dot");
#endif
    const Grid& grid = x.GetGrid();
#ifndef RELEASE
    if( x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y must be distributed over the same grid." << endl;
    }
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be the same length." << endl;
        DumpCallStack();
        throw exception();    
    }
#endif
    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,MR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( grid.MRRank() == ownerCol )
        { 
            T localDot = BLAS::Dot
                         ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            AllReduce
            ( &localDot, &globalDot, 1, MPI_SUM, grid.MCComm() );
        }
        Broadcast( &globalDot, 1, ownerCol, grid.MRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MR,MC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( grid.MCRank() == ownerRow )
        {
            T localDot = BLAS::Dot
                         ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            AllReduce
            ( &localDot, &globalDot, 1, MPI_SUM, grid.MRComm() );
        }
        Broadcast( &globalDot, 1, ownerRow, grid.MCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MR,MC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( grid.MRRank() == ownerCol )
        {
            T localDot = BLAS::Dot
                         ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            AllReduce
            ( &localDot, &globalDot, 1, MPI_SUM, grid.MCComm() );
        }
        Broadcast( &globalDot, 1, ownerCol, grid.MRComm() );
    }
    else
    {
        DistMatrix<T,MC,MR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( grid.MCRank() == ownerRow )
        {
            T localDot = BLAS::Dot
                         ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            AllReduce
            ( &localDot, &globalDot, 1, MPI_SUM, grid.MRComm() );
        }
        Broadcast( &globalDot, 1, ownerRow, grid.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
inline T
Elemental::BLAS::Internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,Star>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Dot");
#endif
    const Grid& grid = x.GetGrid();
#ifndef RELEASE
    if( x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y must be distributed over the same grid." << endl;
    }
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be the same length." << endl;
        DumpCallStack();
        throw exception();    
    }
#endif
    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.MCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,Star,MC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( grid.MCRank() == ownerRow )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, ownerRow, grid.MCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,Star,MC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.MCComm() );
    }
    else
    {
        DistMatrix<T,MC,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( grid.MCRank() == ownerRow )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, ownerRow, grid.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
inline T
Elemental::BLAS::Internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,MR>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Dot");
#endif
    const Grid& grid = x.GetGrid();
#ifndef RELEASE
    if( x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y must be distributed over the same grid." << endl;
    }
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be the same length." << endl;
        DumpCallStack();
        throw exception();    
    }
#endif
    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,Star,MR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( grid.MRRank() == ownerCol )
        { 
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, ownerCol, grid.MRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MR,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.MRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MR,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( grid.MRRank() == ownerCol )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, ownerCol, grid.MRComm() );
    }
    else
    {
        DistMatrix<T,Star,MR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.MRComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
inline T
Elemental::BLAS::Internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MR,MC>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Dot");
#endif
    const Grid& grid = x.GetGrid();
#ifndef RELEASE
    if( x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y must be distributed over the same grid." << endl;
    }
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be the same length." << endl;
        DumpCallStack();
        throw exception();    
    }
#endif
    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MR,MC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( grid.MCRank() == ownerRow )
        { 
            T localDot = BLAS::Dot
                         ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            AllReduce
            ( &localDot, &globalDot, 1, MPI_SUM, grid.MRComm() );
        }
        Broadcast( &globalDot, 1, ownerRow, grid.MCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,MR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( grid.MRRank() == ownerCol )
        {
            T localDot = BLAS::Dot
                         ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            AllReduce
            ( &localDot, &globalDot, 1, MPI_SUM, grid.MCComm() );
        }
        Broadcast( &globalDot, 1, ownerCol, grid.MRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MC,MR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( grid.MCRank() == ownerRow )
        {
            T localDot = BLAS::Dot
                         ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            AllReduce
            ( &localDot, &globalDot, 1, MPI_SUM, grid.MRComm() );
        }
        Broadcast( &globalDot, 1, ownerRow, grid.MCComm() );
    }
    else
    {
        DistMatrix<T,MR,MC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( grid.MRRank() == ownerCol )
        {
            T localDot = BLAS::Dot
                         ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            AllReduce
            ( &localDot, &globalDot, 1, MPI_SUM, grid.MCComm() );
        }
        Broadcast( &globalDot, 1, ownerCol, grid.MRComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
inline T
Elemental::BLAS::Internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MR,Star>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Dot");
#endif
    const Grid& grid = x.GetGrid();
#ifndef RELEASE
    if( x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y must be distributed over the same grid." << endl;
    }
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be the same length." << endl;
        DumpCallStack();
        throw exception();    
    }
#endif
    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MR,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.MRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,Star,MR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( grid.MRRank() == ownerCol )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, ownerCol, grid.MRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,Star,MR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.MRComm() );
    }
    else
    {
        DistMatrix<T,MR,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( grid.MRRank() == ownerCol )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, ownerCol, grid.MRComm() );
    }
#ifndef RELEASE 
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
inline T
Elemental::BLAS::Internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,MC>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Dot");
#endif
    const Grid& grid = x.GetGrid();
#ifndef RELEASE
    if( x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y must be distributed over the same grid." << endl;
    }
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be the same length." << endl;
        DumpCallStack();
        throw exception();    
    }
#endif
    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,Star,MC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( grid.MCRank() == ownerRow )
        { 
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, ownerRow, grid.MCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.MCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MC,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( grid.MCRank() == ownerRow )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, ownerRow, grid.MCComm() );
    }
    else
    {
        DistMatrix<T,Star,MC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
inline T
Elemental::BLAS::Internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,VC,Star>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Dot");
#endif
    const Grid& grid = x.GetGrid();
#ifndef RELEASE
    if( x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y must be distributed over the same grid." << endl;
    }
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be the same length." << endl;
        DumpCallStack();
        throw exception();    
    }
#endif
    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,VC,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.VCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,Star,VC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( grid.VCRank() == owner )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, owner, grid.VCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,Star,VC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.VCComm() );
    }
    else
    {
        DistMatrix<T,VC,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( grid.VCRank() == owner )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, owner, grid.VCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
inline T
Elemental::BLAS::Internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,VC>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Dot");
#endif
    const Grid& grid = x.GetGrid();
#ifndef RELEASE
    if( x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y must be distributed over the same grid." << endl;
    }
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        if( grid.VCRank() == 0 )
            cerr << "Dot requires x and y to be the same length." << endl;
        DumpCallStack();
        throw exception();    
    }
#endif
    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,Star,VC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( grid.VCRank() == owner )
        { 
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, owner, grid.VCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,VC,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.VCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,VC,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( grid.VCRank() == owner )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, owner, grid.VCComm() );
    }
    else
    {
        DistMatrix<T,Star,VC> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.VCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
inline T
Elemental::BLAS::Internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,VR,Star>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Dot");
#endif
    const Grid& grid = x.GetGrid();
#ifndef RELEASE
    if( x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y must be distributed over the same grid." << endl;
    }
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        if( grid.VRRank() == 0 )
            cerr << "Dot requires x and y to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        if( grid.VRRank() == 0 )
            cerr << "Dot requires x and y to be the same length." << endl;
        DumpCallStack();
        throw exception();    
    }
#endif
    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,VR,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.VRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,Star,VR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( grid.VRRank() == owner )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, owner, grid.VRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,Star,VR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.VRComm() );
    }
    else
    {
        DistMatrix<T,VR,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( grid.VRRank() == owner )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, owner, grid.VRComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
inline T
Elemental::BLAS::Internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,VR>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Dot");
#endif
    const Grid& grid = x.GetGrid();
#ifndef RELEASE
    if( x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y must be distributed over the same grid." << endl;
    }
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        if( grid.VRRank() == 0 )
            cerr << "Dot requires x and y to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        if( grid.VRRank() == 0 )
            cerr << "Dot requires x and y to be the same length." << endl;
        DumpCallStack();
        throw exception();    
    }
#endif
    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,Star,VR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( grid.VRRank() == owner )
        { 
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, owner, grid.VRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,VR,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.VRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,VR,Star> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( grid.VRRank() == owner )
        {
            globalDot = BLAS::Dot
                        ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        Broadcast( &globalDot, 1, owner, grid.VRComm() );
    }
    else
    {
        DistMatrix<T,Star,VR> xRedist(grid);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = BLAS::Dot
                     ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        AllReduce( &localDot, &globalDot, 1, MPI_SUM, grid.VRComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
inline T
Elemental::BLAS::Internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,Star>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Dot");
#endif
    const Grid& grid = x.GetGrid();
#ifndef RELEASE
    if( x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y must be distributed over the same grid." << endl;
    }
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        if( grid.VRRank() == 0 )
            cerr << "Dot requires x and y to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        if( grid.VRRank() == 0 )
            cerr << "Dot requires x and y to be the same length." << endl;
        DumpCallStack();
        throw exception();    
    }
#endif
    DistMatrix<T,Star,Star> xRedist(grid);
    xRedist = x;

    T globalDot = BLAS::Dot
                  ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

// We only need to explicitly instantiate the wrapper, since the underlying 
// routines in BLAS::Internal will be immediately called 
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,MC,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,MC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,Star,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,MR,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,MR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,Star,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,VC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,Star,VC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,VR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,Star,VR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,Star,Star>& y );

template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,Star>& x,
  const DistMatrix<float,MC,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,Star>& x,
  const DistMatrix<float,MC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,Star>& x,
  const DistMatrix<float,Star,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,Star>& x,
  const DistMatrix<float,MR,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,Star>& x,
  const DistMatrix<float,MR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,Star>& x,
  const DistMatrix<float,Star,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,Star>& x,
  const DistMatrix<float,VC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,Star>& x,
  const DistMatrix<float,Star,VC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,Star>& x,
  const DistMatrix<float,VR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,Star>& x,
  const DistMatrix<float,Star,VR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MC,Star>& x,
  const DistMatrix<float,Star,Star>& y );

template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MR>& x,
  const DistMatrix<float,MC,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MR>& x,
  const DistMatrix<float,MC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MR>& x,
  const DistMatrix<float,Star,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MR>& x,
  const DistMatrix<float,MR,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MR>& x,
  const DistMatrix<float,MR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MR>& x,
  const DistMatrix<float,Star,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MR>& x,
  const DistMatrix<float,VC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MR>& x,
  const DistMatrix<float,Star,VC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MR>& x,
  const DistMatrix<float,VR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MR>& x,
  const DistMatrix<float,Star,VR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MR>& x,
  const DistMatrix<float,Star,Star>& y );

template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,MC,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,MC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,Star,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,MR,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,MR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,Star,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,VC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,Star,VC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,VR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,Star,VR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,Star,Star>& y );

template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,Star>& x,
  const DistMatrix<float,MC,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,Star>& x,
  const DistMatrix<float,MC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,Star>& x,
  const DistMatrix<float,Star,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,Star>& x,
  const DistMatrix<float,MR,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,Star>& x,
  const DistMatrix<float,MR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,Star>& x,
  const DistMatrix<float,Star,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,Star>& x,
  const DistMatrix<float,VC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,Star>& x,
  const DistMatrix<float,Star,VC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,Star>& x,
  const DistMatrix<float,VR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,Star>& x,
  const DistMatrix<float,Star,VR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,MR,Star>& x,
  const DistMatrix<float,Star,Star>& y );

template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MC>& x,
  const DistMatrix<float,MC,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MC>& x,
  const DistMatrix<float,MC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MC>& x,
  const DistMatrix<float,Star,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MC>& x,
  const DistMatrix<float,MR,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MC>& x,
  const DistMatrix<float,MR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MC>& x,
  const DistMatrix<float,Star,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MC>& x,
  const DistMatrix<float,VC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MC>& x,
  const DistMatrix<float,Star,VC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MC>& x,
  const DistMatrix<float,VR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MC>& x,
  const DistMatrix<float,Star,VR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,MC>& x,
  const DistMatrix<float,Star,Star>& y );

template float Elemental::BLAS::Dot
( const DistMatrix<float,VC,Star>& x,
  const DistMatrix<float,MC,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VC,Star>& x,
  const DistMatrix<float,MC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VC,Star>& x,
  const DistMatrix<float,Star,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VC,Star>& x,
  const DistMatrix<float,MR,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VC,Star>& x,
  const DistMatrix<float,MR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VC,Star>& x,
  const DistMatrix<float,Star,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VC,Star>& x,
  const DistMatrix<float,VC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VC,Star>& x,
  const DistMatrix<float,Star,VC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VC,Star>& x,
  const DistMatrix<float,VR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VC,Star>& x,
  const DistMatrix<float,Star,VR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VC,Star>& x,
  const DistMatrix<float,Star,Star>& y );

template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VC>& x,
  const DistMatrix<float,MC,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VC>& x,
  const DistMatrix<float,MC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VC>& x,
  const DistMatrix<float,Star,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VC>& x,
  const DistMatrix<float,MR,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VC>& x,
  const DistMatrix<float,MR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VC>& x,
  const DistMatrix<float,Star,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VC>& x,
  const DistMatrix<float,VC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VC>& x,
  const DistMatrix<float,Star,VC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VC>& x,
  const DistMatrix<float,VR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VC>& x,
  const DistMatrix<float,Star,VR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VC>& x,
  const DistMatrix<float,Star,Star>& y );

template float Elemental::BLAS::Dot
( const DistMatrix<float,VR,Star>& x,
  const DistMatrix<float,MC,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VR,Star>& x,
  const DistMatrix<float,MC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VR,Star>& x,
  const DistMatrix<float,Star,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VR,Star>& x,
  const DistMatrix<float,MR,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VR,Star>& x,
  const DistMatrix<float,MR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VR,Star>& x,
  const DistMatrix<float,Star,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VR,Star>& x,
  const DistMatrix<float,VC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VR,Star>& x,
  const DistMatrix<float,Star,VC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VR,Star>& x,
  const DistMatrix<float,VR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VR,Star>& x,
  const DistMatrix<float,Star,VR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,VR,Star>& x,
  const DistMatrix<float,Star,Star>& y );

template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VR>& x,
  const DistMatrix<float,MC,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VR>& x,
  const DistMatrix<float,MC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VR>& x,
  const DistMatrix<float,Star,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VR>& x,
  const DistMatrix<float,MR,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VR>& x,
  const DistMatrix<float,MR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VR>& x,
  const DistMatrix<float,Star,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VR>& x,
  const DistMatrix<float,VC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VR>& x,
  const DistMatrix<float,Star,VC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VR>& x,
  const DistMatrix<float,VR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VR>& x,
  const DistMatrix<float,Star,VR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,VR>& x,
  const DistMatrix<float,Star,Star>& y );

template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,Star>& x,
  const DistMatrix<float,MC,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,Star>& x,
  const DistMatrix<float,MC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,Star>& x,
  const DistMatrix<float,Star,MR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,Star>& x,
  const DistMatrix<float,MR,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,Star>& x,
  const DistMatrix<float,MR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,Star>& x,
  const DistMatrix<float,Star,MC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,Star>& x,
  const DistMatrix<float,VC,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,Star>& x,
  const DistMatrix<float,Star,VC>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,Star>& x,
  const DistMatrix<float,VR,Star>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,Star>& x,
  const DistMatrix<float,Star,VR>& y );
template float Elemental::BLAS::Dot
( const DistMatrix<float,Star,Star>& x,
  const DistMatrix<float,Star,Star>& y );


template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,MC,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,MC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,Star,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,MR,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,MR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,Star,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,VC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,Star,VC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,VR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,Star,VR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,Star,Star>& y );

template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,Star>& x,
  const DistMatrix<double,MC,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,Star>& x,
  const DistMatrix<double,MC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,Star>& x,
  const DistMatrix<double,Star,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,Star>& x,
  const DistMatrix<double,MR,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,Star>& x,
  const DistMatrix<double,MR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,Star>& x,
  const DistMatrix<double,Star,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,Star>& x,
  const DistMatrix<double,VC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,Star>& x,
  const DistMatrix<double,Star,VC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,Star>& x,
  const DistMatrix<double,VR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,Star>& x,
  const DistMatrix<double,Star,VR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MC,Star>& x,
  const DistMatrix<double,Star,Star>& y );

template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MR>& x,
  const DistMatrix<double,MC,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MR>& x,
  const DistMatrix<double,MC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MR>& x,
  const DistMatrix<double,Star,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MR>& x,
  const DistMatrix<double,MR,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MR>& x,
  const DistMatrix<double,MR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MR>& x,
  const DistMatrix<double,Star,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MR>& x,
  const DistMatrix<double,VC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MR>& x,
  const DistMatrix<double,Star,VC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MR>& x,
  const DistMatrix<double,VR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MR>& x,
  const DistMatrix<double,Star,VR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MR>& x,
  const DistMatrix<double,Star,Star>& y );

template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,MC,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,MC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,Star,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,MR,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,MR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,Star,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,VC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,Star,VC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,VR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,Star,VR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,Star,Star>& y );

template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,Star>& x,
  const DistMatrix<double,MC,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,Star>& x,
  const DistMatrix<double,MC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,Star>& x,
  const DistMatrix<double,Star,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,Star>& x,
  const DistMatrix<double,MR,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,Star>& x,
  const DistMatrix<double,MR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,Star>& x,
  const DistMatrix<double,Star,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,Star>& x,
  const DistMatrix<double,VC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,Star>& x,
  const DistMatrix<double,Star,VC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,Star>& x,
  const DistMatrix<double,VR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,Star>& x,
  const DistMatrix<double,Star,VR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,MR,Star>& x,
  const DistMatrix<double,Star,Star>& y );

template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MC>& x,
  const DistMatrix<double,MC,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MC>& x,
  const DistMatrix<double,MC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MC>& x,
  const DistMatrix<double,Star,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MC>& x,
  const DistMatrix<double,MR,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MC>& x,
  const DistMatrix<double,MR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MC>& x,
  const DistMatrix<double,Star,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MC>& x,
  const DistMatrix<double,VC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MC>& x,
  const DistMatrix<double,Star,VC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MC>& x,
  const DistMatrix<double,VR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MC>& x,
  const DistMatrix<double,Star,VR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,MC>& x,
  const DistMatrix<double,Star,Star>& y );

template double Elemental::BLAS::Dot
( const DistMatrix<double,VC,Star>& x,
  const DistMatrix<double,MC,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VC,Star>& x,
  const DistMatrix<double,MC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VC,Star>& x,
  const DistMatrix<double,Star,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VC,Star>& x,
  const DistMatrix<double,MR,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VC,Star>& x,
  const DistMatrix<double,MR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VC,Star>& x,
  const DistMatrix<double,Star,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VC,Star>& x,
  const DistMatrix<double,VC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VC,Star>& x,
  const DistMatrix<double,Star,VC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VC,Star>& x,
  const DistMatrix<double,VR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VC,Star>& x,
  const DistMatrix<double,Star,VR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VC,Star>& x,
  const DistMatrix<double,Star,Star>& y );

template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VC>& x,
  const DistMatrix<double,MC,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VC>& x,
  const DistMatrix<double,MC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VC>& x,
  const DistMatrix<double,Star,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VC>& x,
  const DistMatrix<double,MR,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VC>& x,
  const DistMatrix<double,MR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VC>& x,
  const DistMatrix<double,Star,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VC>& x,
  const DistMatrix<double,VC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VC>& x,
  const DistMatrix<double,Star,VC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VC>& x,
  const DistMatrix<double,VR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VC>& x,
  const DistMatrix<double,Star,VR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VC>& x,
  const DistMatrix<double,Star,Star>& y );

template double Elemental::BLAS::Dot
( const DistMatrix<double,VR,Star>& x,
  const DistMatrix<double,MC,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VR,Star>& x,
  const DistMatrix<double,MC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VR,Star>& x,
  const DistMatrix<double,Star,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VR,Star>& x,
  const DistMatrix<double,MR,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VR,Star>& x,
  const DistMatrix<double,MR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VR,Star>& x,
  const DistMatrix<double,Star,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VR,Star>& x,
  const DistMatrix<double,VC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VR,Star>& x,
  const DistMatrix<double,Star,VC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VR,Star>& x,
  const DistMatrix<double,VR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VR,Star>& x,
  const DistMatrix<double,Star,VR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,VR,Star>& x,
  const DistMatrix<double,Star,Star>& y );

template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VR>& x,
  const DistMatrix<double,MC,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VR>& x,
  const DistMatrix<double,MC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VR>& x,
  const DistMatrix<double,Star,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VR>& x,
  const DistMatrix<double,MR,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VR>& x,
  const DistMatrix<double,MR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VR>& x,
  const DistMatrix<double,Star,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VR>& x,
  const DistMatrix<double,VC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VR>& x,
  const DistMatrix<double,Star,VC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VR>& x,
  const DistMatrix<double,VR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VR>& x,
  const DistMatrix<double,Star,VR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,VR>& x,
  const DistMatrix<double,Star,Star>& y );

template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,Star>& x,
  const DistMatrix<double,MC,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,Star>& x,
  const DistMatrix<double,MC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,Star>& x,
  const DistMatrix<double,Star,MR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,Star>& x,
  const DistMatrix<double,MR,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,Star>& x,
  const DistMatrix<double,MR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,Star>& x,
  const DistMatrix<double,Star,MC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,Star>& x,
  const DistMatrix<double,VC,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,Star>& x,
  const DistMatrix<double,Star,VC>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,Star>& x,
  const DistMatrix<double,VR,Star>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,Star>& x,
  const DistMatrix<double,Star,VR>& y );
template double Elemental::BLAS::Dot
( const DistMatrix<double,Star,Star>& x,
  const DistMatrix<double,Star,Star>& y );

#ifndef WITHOUT_COMPLEX
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,MC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,Star,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,MR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,Star,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,VC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,Star,VC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,VR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,Star,VR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,Star,Star>& y );

template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,Star>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,Star>& x,
  const DistMatrix<scomplex,MC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,Star>& x,
  const DistMatrix<scomplex,Star,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,Star>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,Star>& x,
  const DistMatrix<scomplex,MR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,Star>& x,
  const DistMatrix<scomplex,Star,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,Star>& x,
  const DistMatrix<scomplex,VC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,Star>& x,
  const DistMatrix<scomplex,Star,VC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,Star>& x,
  const DistMatrix<scomplex,VR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,Star>& x,
  const DistMatrix<scomplex,Star,VR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MC,Star>& x,
  const DistMatrix<scomplex,Star,Star>& y );

template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MR>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MR>& x,
  const DistMatrix<scomplex,MC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MR>& x,
  const DistMatrix<scomplex,Star,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MR>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MR>& x,
  const DistMatrix<scomplex,MR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MR>& x,
  const DistMatrix<scomplex,Star,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MR>& x,
  const DistMatrix<scomplex,VC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MR>& x,
  const DistMatrix<scomplex,Star,VC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MR>& x,
  const DistMatrix<scomplex,VR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MR>& x,
  const DistMatrix<scomplex,Star,VR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MR>& x,
  const DistMatrix<scomplex,Star,Star>& y );

template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,MC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,Star,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,MR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,Star,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,VC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,Star,VC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,VR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,Star,VR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,Star,Star>& y );

template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,Star>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,Star>& x,
  const DistMatrix<scomplex,MC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,Star>& x,
  const DistMatrix<scomplex,Star,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,Star>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,Star>& x,
  const DistMatrix<scomplex,MR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,Star>& x,
  const DistMatrix<scomplex,Star,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,Star>& x,
  const DistMatrix<scomplex,VC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,Star>& x,
  const DistMatrix<scomplex,Star,VC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,Star>& x,
  const DistMatrix<scomplex,VR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,Star>& x,
  const DistMatrix<scomplex,Star,VR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,MR,Star>& x,
  const DistMatrix<scomplex,Star,Star>& y );

template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MC>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MC>& x,
  const DistMatrix<scomplex,MC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MC>& x,
  const DistMatrix<scomplex,Star,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MC>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MC>& x,
  const DistMatrix<scomplex,MR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MC>& x,
  const DistMatrix<scomplex,Star,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MC>& x,
  const DistMatrix<scomplex,VC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MC>& x,
  const DistMatrix<scomplex,Star,VC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MC>& x,
  const DistMatrix<scomplex,VR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MC>& x,
  const DistMatrix<scomplex,Star,VR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,MC>& x,
  const DistMatrix<scomplex,Star,Star>& y );

template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VC,Star>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VC,Star>& x,
  const DistMatrix<scomplex,MC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VC,Star>& x,
  const DistMatrix<scomplex,Star,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VC,Star>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VC,Star>& x,
  const DistMatrix<scomplex,MR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VC,Star>& x,
  const DistMatrix<scomplex,Star,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VC,Star>& x,
  const DistMatrix<scomplex,VC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VC,Star>& x,
  const DistMatrix<scomplex,Star,VC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VC,Star>& x,
  const DistMatrix<scomplex,VR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VC,Star>& x,
  const DistMatrix<scomplex,Star,VR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VC,Star>& x,
  const DistMatrix<scomplex,Star,Star>& y );

template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VC>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VC>& x,
  const DistMatrix<scomplex,MC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VC>& x,
  const DistMatrix<scomplex,Star,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VC>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VC>& x,
  const DistMatrix<scomplex,MR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VC>& x,
  const DistMatrix<scomplex,Star,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VC>& x,
  const DistMatrix<scomplex,VC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VC>& x,
  const DistMatrix<scomplex,Star,VC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VC>& x,
  const DistMatrix<scomplex,VR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VC>& x,
  const DistMatrix<scomplex,Star,VR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VC>& x,
  const DistMatrix<scomplex,Star,Star>& y );

template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VR,Star>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VR,Star>& x,
  const DistMatrix<scomplex,MC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VR,Star>& x,
  const DistMatrix<scomplex,Star,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VR,Star>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VR,Star>& x,
  const DistMatrix<scomplex,MR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VR,Star>& x,
  const DistMatrix<scomplex,Star,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VR,Star>& x,
  const DistMatrix<scomplex,VC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VR,Star>& x,
  const DistMatrix<scomplex,Star,VC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VR,Star>& x,
  const DistMatrix<scomplex,VR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VR,Star>& x,
  const DistMatrix<scomplex,Star,VR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,VR,Star>& x,
  const DistMatrix<scomplex,Star,Star>& y );

template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VR>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VR>& x,
  const DistMatrix<scomplex,MC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VR>& x,
  const DistMatrix<scomplex,Star,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VR>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VR>& x,
  const DistMatrix<scomplex,MR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VR>& x,
  const DistMatrix<scomplex,Star,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VR>& x,
  const DistMatrix<scomplex,VC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VR>& x,
  const DistMatrix<scomplex,Star,VC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VR>& x,
  const DistMatrix<scomplex,VR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VR>& x,
  const DistMatrix<scomplex,Star,VR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,VR>& x,
  const DistMatrix<scomplex,Star,Star>& y );

template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,Star>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,Star>& x,
  const DistMatrix<scomplex,MC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,Star>& x,
  const DistMatrix<scomplex,Star,MR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,Star>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,Star>& x,
  const DistMatrix<scomplex,MR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,Star>& x,
  const DistMatrix<scomplex,Star,MC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,Star>& x,
  const DistMatrix<scomplex,VC,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,Star>& x,
  const DistMatrix<scomplex,Star,VC>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,Star>& x,
  const DistMatrix<scomplex,VR,Star>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,Star>& x,
  const DistMatrix<scomplex,Star,VR>& y );
template scomplex Elemental::BLAS::Dot
( const DistMatrix<scomplex,Star,Star>& x,
  const DistMatrix<scomplex,Star,Star>& y );

template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,MC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,Star,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,MR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,Star,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,VC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,Star,VC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,VR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,Star,VR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,Star,Star>& y );

template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,Star>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,Star>& x,
  const DistMatrix<dcomplex,MC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,Star>& x,
  const DistMatrix<dcomplex,Star,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,Star>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,Star>& x,
  const DistMatrix<dcomplex,MR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,Star>& x,
  const DistMatrix<dcomplex,Star,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,Star>& x,
  const DistMatrix<dcomplex,VC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,Star>& x,
  const DistMatrix<dcomplex,Star,VC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,Star>& x,
  const DistMatrix<dcomplex,VR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,Star>& x,
  const DistMatrix<dcomplex,Star,VR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MC,Star>& x,
  const DistMatrix<dcomplex,Star,Star>& y );

template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MR>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MR>& x,
  const DistMatrix<dcomplex,MC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MR>& x,
  const DistMatrix<dcomplex,Star,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MR>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MR>& x,
  const DistMatrix<dcomplex,MR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MR>& x,
  const DistMatrix<dcomplex,Star,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MR>& x,
  const DistMatrix<dcomplex,VC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MR>& x,
  const DistMatrix<dcomplex,Star,VC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MR>& x,
  const DistMatrix<dcomplex,VR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MR>& x,
  const DistMatrix<dcomplex,Star,VR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MR>& x,
  const DistMatrix<dcomplex,Star,Star>& y );

template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,MC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,Star,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,MR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,Star,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,VC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,Star,VC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,VR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,Star,VR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,Star,Star>& y );

template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,Star>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,Star>& x,
  const DistMatrix<dcomplex,MC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,Star>& x,
  const DistMatrix<dcomplex,Star,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,Star>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,Star>& x,
  const DistMatrix<dcomplex,MR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,Star>& x,
  const DistMatrix<dcomplex,Star,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,Star>& x,
  const DistMatrix<dcomplex,VC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,Star>& x,
  const DistMatrix<dcomplex,Star,VC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,Star>& x,
  const DistMatrix<dcomplex,VR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,Star>& x,
  const DistMatrix<dcomplex,Star,VR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,MR,Star>& x,
  const DistMatrix<dcomplex,Star,Star>& y );

template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MC>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MC>& x,
  const DistMatrix<dcomplex,MC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MC>& x,
  const DistMatrix<dcomplex,Star,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MC>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MC>& x,
  const DistMatrix<dcomplex,MR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MC>& x,
  const DistMatrix<dcomplex,Star,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MC>& x,
  const DistMatrix<dcomplex,VC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MC>& x,
  const DistMatrix<dcomplex,Star,VC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MC>& x,
  const DistMatrix<dcomplex,VR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MC>& x,
  const DistMatrix<dcomplex,Star,VR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,MC>& x,
  const DistMatrix<dcomplex,Star,Star>& y );

template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VC,Star>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VC,Star>& x,
  const DistMatrix<dcomplex,MC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VC,Star>& x,
  const DistMatrix<dcomplex,Star,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VC,Star>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VC,Star>& x,
  const DistMatrix<dcomplex,MR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VC,Star>& x,
  const DistMatrix<dcomplex,Star,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VC,Star>& x,
  const DistMatrix<dcomplex,VC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VC,Star>& x,
  const DistMatrix<dcomplex,Star,VC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VC,Star>& x,
  const DistMatrix<dcomplex,VR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VC,Star>& x,
  const DistMatrix<dcomplex,Star,VR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VC,Star>& x,
  const DistMatrix<dcomplex,Star,Star>& y );

template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VC>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VC>& x,
  const DistMatrix<dcomplex,MC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VC>& x,
  const DistMatrix<dcomplex,Star,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VC>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VC>& x,
  const DistMatrix<dcomplex,MR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VC>& x,
  const DistMatrix<dcomplex,Star,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VC>& x,
  const DistMatrix<dcomplex,VC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VC>& x,
  const DistMatrix<dcomplex,Star,VC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VC>& x,
  const DistMatrix<dcomplex,VR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VC>& x,
  const DistMatrix<dcomplex,Star,VR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VC>& x,
  const DistMatrix<dcomplex,Star,Star>& y );

template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VR,Star>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VR,Star>& x,
  const DistMatrix<dcomplex,MC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VR,Star>& x,
  const DistMatrix<dcomplex,Star,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VR,Star>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VR,Star>& x,
  const DistMatrix<dcomplex,MR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VR,Star>& x,
  const DistMatrix<dcomplex,Star,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VR,Star>& x,
  const DistMatrix<dcomplex,VC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VR,Star>& x,
  const DistMatrix<dcomplex,Star,VC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VR,Star>& x,
  const DistMatrix<dcomplex,VR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VR,Star>& x,
  const DistMatrix<dcomplex,Star,VR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,VR,Star>& x,
  const DistMatrix<dcomplex,Star,Star>& y );

template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VR>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VR>& x,
  const DistMatrix<dcomplex,MC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VR>& x,
  const DistMatrix<dcomplex,Star,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VR>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VR>& x,
  const DistMatrix<dcomplex,MR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VR>& x,
  const DistMatrix<dcomplex,Star,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VR>& x,
  const DistMatrix<dcomplex,VC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VR>& x,
  const DistMatrix<dcomplex,Star,VC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VR>& x,
  const DistMatrix<dcomplex,VR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VR>& x,
  const DistMatrix<dcomplex,Star,VR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,VR>& x,
  const DistMatrix<dcomplex,Star,Star>& y );

template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,Star>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,Star>& x,
  const DistMatrix<dcomplex,MC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,Star>& x,
  const DistMatrix<dcomplex,Star,MR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,Star>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,Star>& x,
  const DistMatrix<dcomplex,MR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,Star>& x,
  const DistMatrix<dcomplex,Star,MC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,Star>& x,
  const DistMatrix<dcomplex,VC,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,Star>& x,
  const DistMatrix<dcomplex,Star,VC>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,Star>& x,
  const DistMatrix<dcomplex,VR,Star>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,Star>& x,
  const DistMatrix<dcomplex,Star,VR>& y );
template dcomplex Elemental::BLAS::Dot
( const DistMatrix<dcomplex,Star,Star>& x,
  const DistMatrix<dcomplex,Star,Star>& y );
#endif
