/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/basic_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::imports;

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

/* 
   C++ does not (currently) allow for partial function template specialization,
   so implementing Dot will be unnecessarily obfuscated. Sorry.
  
   The compromise is to have the user-level routine, 
  
     template<typename T, Distribution U, Distribution V,
                          Distribution W, Distribution Z >
     T
     basic::Dot( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

   simply route to overloaded pseudo-partial-specializations within the 
   basic::internal namespace. For example,

     template<typename T, Distribution U, Distribution V>
     T
     basic::internal::Dot
     ( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,MR>& y );
*/
template<typename T, Distribution U, Distribution V,
                     Distribution W, Distribution Z >
T
elemental::basic::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )
{
#ifndef RELEASE
    PushCallStack("basic::Dot");
#endif
    T dotProduct = basic::internal::Dot( x, y );
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename T, Distribution U, Distribution V>
T
elemental::basic::internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::Dot");
    if( x.Grid() != y.Grid() )
        throw logic_error( "x and y must be distributed over the same grid." );
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw logic_error( "Dot requires x and y to be vectors." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw logic_error( "Dot requires x and y to be the same length." );
#endif
    const Grid& g = x.Grid();

    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( g.MRRank() == ownerCol )
        { 
            T localDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce
            ( &localDot, &globalDot, 1, mpi::SUM, g.MCComm() );
        }
        mpi::Broadcast( &globalDot, 1, ownerCol, g.MRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( g.MCRank() == ownerRow )
        {
            T localDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce
            ( &localDot, &globalDot, 1, mpi::SUM, g.MRComm() );
        }
        mpi::Broadcast( &globalDot, 1, ownerRow, g.MCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( g.MRRank() == ownerCol )
        {
            T localDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce
            ( &localDot, &globalDot, 1, mpi::SUM, g.MCComm() );
        }
        mpi::Broadcast( &globalDot, 1, ownerCol, g.MRComm() );
    }
    else
    {
        DistMatrix<T,MC,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( g.MCRank() == ownerRow )
        {
            T localDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce
            ( &localDot, &globalDot, 1, mpi::SUM, g.MRComm() );
        }
        mpi::Broadcast( &globalDot, 1, ownerRow, g.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

#ifdef ENABLE_ALL_DISTRIBUTED_DOT
template<typename T, Distribution U, Distribution V>
T
elemental::basic::internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::Dot");
    if( x.Grid() != y.Grid() )
        throw logic_error( "x and y must be distributed over the same grid." );
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw logic_error( "Dot requires x and y to be vectors." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw logic_error( "Dot requires x and y to be the same length." );
#endif
    const Grid& g = x.Grid();

    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.MCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,STAR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( g.MCRank() == ownerRow )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, ownerRow, g.MCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.MCComm() );
    }
    else
    {
        DistMatrix<T,MC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( g.MCRank() == ownerRow )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, ownerRow, g.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
T
elemental::basic::internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,STAR,MR>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::Dot");
    if( x.Grid() != y.Grid() )
        throw logic_error( "x and y must be distributed over the same grid." );
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw logic_error( "Dot requires x and y to be vectors." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw logic_error( "Dot requires x and y to be the same length." );
#endif
    const Grid& g = x.Grid();

    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,STAR,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( g.MRRank() == ownerCol )
        { 
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, ownerCol, g.MRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.MRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( g.MRRank() == ownerCol )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, ownerCol, g.MRComm() );
    }
    else
    {
        DistMatrix<T,STAR,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.MRComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
T
elemental::basic::internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MR,MC>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::Dot");
    if( x.Grid() != y.Grid() )
        throw logic_error( "x and y must be distributed over the same grid." );
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw logic_error( "Dot requires x and y to be vectors." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw logic_error( "Dot requires x and y to be the same length." );
#endif
    const Grid& g = x.Grid();

    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( g.MCRank() == ownerRow )
        { 
            T localDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce
            ( &localDot, &globalDot, 1, mpi::SUM, g.MRComm() );
        }
        mpi::Broadcast( &globalDot, 1, ownerRow, g.MCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( g.MRRank() == ownerCol )
        {
            T localDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce
            ( &localDot, &globalDot, 1, mpi::SUM, g.MCComm() );
        }
        mpi::Broadcast( &globalDot, 1, ownerCol, g.MRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MC,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( g.MCRank() == ownerRow )
        {
            T localDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce
            ( &localDot, &globalDot, 1, mpi::SUM, g.MRComm() );
        }
        mpi::Broadcast( &globalDot, 1, ownerRow, g.MCComm() );
    }
    else
    {
        DistMatrix<T,MR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( g.MRRank() == ownerCol )
        {
            T localDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce
            ( &localDot, &globalDot, 1, mpi::SUM, g.MCComm() );
        }
        mpi::Broadcast( &globalDot, 1, ownerCol, g.MRComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
T
elemental::basic::internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MR,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::Dot");
    if( x.Grid() != y.Grid() )
        throw logic_error( "x and y must be distributed over the same grid." );
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw logic_error( "Dot requires x and y to be vectors." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw logic_error( "Dot requires x and y to be the same length." );
#endif
    const Grid& g = x.Grid();

    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.MRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,STAR,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( g.MRRank() == ownerCol )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, ownerCol, g.MRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.MRComm() );
    }
    else
    {
        DistMatrix<T,MR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( g.MRRank() == ownerCol )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, ownerCol, g.MRComm() );
    }
#ifndef RELEASE 
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
T
elemental::basic::internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,STAR,MC>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::Dot");
    if( x.Grid() != y.Grid() )
        throw logic_error( "x and y must be distributed over the same grid." );
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw logic_error( "Dot requires x and y to be vectors." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw logic_error( "Dot requires x and y to be the same length." );
#endif
    const Grid& g = x.Grid();

    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,STAR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( g.MCRank() == ownerRow )
        { 
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, ownerRow, g.MCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.MCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( g.MCRank() == ownerRow )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, ownerRow, g.MCComm() );
    }
    else
    {
        DistMatrix<T,STAR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
T
elemental::basic::internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,VC,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::Dot");
    if( x.Grid() != y.Grid() )
        throw logic_error( "x and y must be distributed over the same grid." );
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw logic_error( "Dot requires x and y to be vectors." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw logic_error( "Dot requires x and y to be the same length." );
#endif
    const Grid& g = x.Grid();

    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,VC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.VCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,STAR,VC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( g.VCRank() == owner )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, owner, g.VCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,VC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.VCComm() );
    }
    else
    {
        DistMatrix<T,VC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( g.VCRank() == owner )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, owner, g.VCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
T
elemental::basic::internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,STAR,VC>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::Dot");
    if( x.Grid() != y.Grid() )
        throw logic_error( "x and y must be distributed over the same grid." );
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw logic_error( "Dot requires x and y to be vectors." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw logic_error( "Dot requires x and y to be the same length." );
#endif
    const Grid& g = x.Grid();

    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,STAR,VC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( g.VCRank() == owner )
        { 
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, owner, g.VCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,VC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.VCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,VC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( g.VCRank() == owner )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, owner, g.VCComm() );
    }
    else
    {
        DistMatrix<T,STAR,VC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.VCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
T
elemental::basic::internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,VR,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::Dot");
    if( x.Grid() != y.Grid() )
        throw logic_error( "x and y must be distributed over the same grid." );
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw logic_error( "Dot requires x and y to be vectors." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw logic_error( "Dot requires x and y to be the same length." );
#endif
    const Grid& g = x.Grid();

    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,VR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.VRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,STAR,VR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( g.VRRank() == owner )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, owner, g.VRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,VR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.VRComm() );
    }
    else
    {
        DistMatrix<T,VR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( g.VRRank() == owner )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, owner, g.VRComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
T
elemental::basic::internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,STAR,VR>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::Dot");
    if( x.Grid() != y.Grid() )
        throw logic_error( "x and y must be distributed over the same grid." );
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw logic_error( "Dot requires x and y to be vectors." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw logic_error( "Dot requires x and y to be the same length." );
#endif
    const Grid& g = x.Grid();

    T globalDot;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,STAR,VR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( g.VRRank() == owner )
        { 
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, owner, g.VRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,VR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.VRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,VR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( g.VRRank() == owner )
        {
            globalDot = 
                basic::Dot
                ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDot, 1, owner, g.VRComm() );
    }
    else
    {
        DistMatrix<T,STAR,VR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDot = 
            basic::Dot
            ( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDot, &globalDot, 1, mpi::SUM, g.VRComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Distribution U, Distribution V>
T
elemental::basic::internal::Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,STAR,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::Dot");
    if( x.Grid() != y.Grid() )
        throw logic_error( "x and y must be distributed over the same grid." );
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw logic_error( "Dot requires x and y to be vectors." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw logic_error( "Dot requires x and y to be the same length." );
#endif
    const Grid& g = x.Grid();

    DistMatrix<T,STAR,STAR> xRedist(g);
    xRedist = x;

    T globalDot = 
        basic::Dot( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}
#endif // ENABLE_ALL_DISTRIBUTED_DOT

// We only need to explicitly instantiate the wrapper, since the underlying 
// routines in basic::internal will be immediately called 
template float elemental::basic::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,MC,MR>& y );
#ifdef ENABLE_ALL_DISTRIBUTED_DOT
template float elemental::basic::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,MC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,STAR,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,MR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,MR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,STAR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,VC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,STAR,VC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,VR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,STAR,VR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,MR>& x,
  const DistMatrix<float,STAR,STAR>& y );

template float elemental::basic::Dot
( const DistMatrix<float,MC,STAR>& x,
  const DistMatrix<float,MC,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,STAR>& x,
  const DistMatrix<float,MC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,STAR>& x,
  const DistMatrix<float,STAR,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,STAR>& x,
  const DistMatrix<float,MR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,STAR>& x,
  const DistMatrix<float,MR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,STAR>& x,
  const DistMatrix<float,STAR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,STAR>& x,
  const DistMatrix<float,VC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,STAR>& x,
  const DistMatrix<float,STAR,VC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,STAR>& x,
  const DistMatrix<float,VR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,STAR>& x,
  const DistMatrix<float,STAR,VR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MC,STAR>& x,
  const DistMatrix<float,STAR,STAR>& y );

template float elemental::basic::Dot
( const DistMatrix<float,STAR,MR>& x,
  const DistMatrix<float,MC,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MR>& x,
  const DistMatrix<float,MC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MR>& x,
  const DistMatrix<float,STAR,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MR>& x,
  const DistMatrix<float,MR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MR>& x,
  const DistMatrix<float,MR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MR>& x,
  const DistMatrix<float,STAR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MR>& x,
  const DistMatrix<float,VC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MR>& x,
  const DistMatrix<float,STAR,VC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MR>& x,
  const DistMatrix<float,VR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MR>& x,
  const DistMatrix<float,STAR,VR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MR>& x,
  const DistMatrix<float,STAR,STAR>& y );

template float elemental::basic::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,MC,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,MC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,STAR,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,MR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,MR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,STAR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,VC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,STAR,VC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,VR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,STAR,VR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,MC>& x,
  const DistMatrix<float,STAR,STAR>& y );

template float elemental::basic::Dot
( const DistMatrix<float,MR,STAR>& x,
  const DistMatrix<float,MC,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,STAR>& x,
  const DistMatrix<float,MC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,STAR>& x,
  const DistMatrix<float,STAR,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,STAR>& x,
  const DistMatrix<float,MR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,STAR>& x,
  const DistMatrix<float,MR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,STAR>& x,
  const DistMatrix<float,STAR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,STAR>& x,
  const DistMatrix<float,VC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,STAR>& x,
  const DistMatrix<float,STAR,VC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,STAR>& x,
  const DistMatrix<float,VR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,STAR>& x,
  const DistMatrix<float,STAR,VR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,MR,STAR>& x,
  const DistMatrix<float,STAR,STAR>& y );

template float elemental::basic::Dot
( const DistMatrix<float,STAR,MC>& x,
  const DistMatrix<float,MC,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MC>& x,
  const DistMatrix<float,MC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MC>& x,
  const DistMatrix<float,STAR,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MC>& x,
  const DistMatrix<float,MR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MC>& x,
  const DistMatrix<float,MR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MC>& x,
  const DistMatrix<float,STAR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MC>& x,
  const DistMatrix<float,VC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MC>& x,
  const DistMatrix<float,STAR,VC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MC>& x,
  const DistMatrix<float,VR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MC>& x,
  const DistMatrix<float,STAR,VR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,MC>& x,
  const DistMatrix<float,STAR,STAR>& y );

template float elemental::basic::Dot
( const DistMatrix<float,VC,STAR>& x,
  const DistMatrix<float,MC,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VC,STAR>& x,
  const DistMatrix<float,MC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VC,STAR>& x,
  const DistMatrix<float,STAR,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VC,STAR>& x,
  const DistMatrix<float,MR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VC,STAR>& x,
  const DistMatrix<float,MR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VC,STAR>& x,
  const DistMatrix<float,STAR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VC,STAR>& x,
  const DistMatrix<float,VC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VC,STAR>& x,
  const DistMatrix<float,STAR,VC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VC,STAR>& x,
  const DistMatrix<float,VR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VC,STAR>& x,
  const DistMatrix<float,STAR,VR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VC,STAR>& x,
  const DistMatrix<float,STAR,STAR>& y );

template float elemental::basic::Dot
( const DistMatrix<float,STAR,VC>& x,
  const DistMatrix<float,MC,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VC>& x,
  const DistMatrix<float,MC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VC>& x,
  const DistMatrix<float,STAR,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VC>& x,
  const DistMatrix<float,MR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VC>& x,
  const DistMatrix<float,MR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VC>& x,
  const DistMatrix<float,STAR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VC>& x,
  const DistMatrix<float,VC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VC>& x,
  const DistMatrix<float,STAR,VC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VC>& x,
  const DistMatrix<float,VR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VC>& x,
  const DistMatrix<float,STAR,VR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VC>& x,
  const DistMatrix<float,STAR,STAR>& y );

template float elemental::basic::Dot
( const DistMatrix<float,VR,STAR>& x,
  const DistMatrix<float,MC,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VR,STAR>& x,
  const DistMatrix<float,MC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VR,STAR>& x,
  const DistMatrix<float,STAR,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VR,STAR>& x,
  const DistMatrix<float,MR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VR,STAR>& x,
  const DistMatrix<float,MR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VR,STAR>& x,
  const DistMatrix<float,STAR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VR,STAR>& x,
  const DistMatrix<float,VC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VR,STAR>& x,
  const DistMatrix<float,STAR,VC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VR,STAR>& x,
  const DistMatrix<float,VR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VR,STAR>& x,
  const DistMatrix<float,STAR,VR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,VR,STAR>& x,
  const DistMatrix<float,STAR,STAR>& y );

template float elemental::basic::Dot
( const DistMatrix<float,STAR,VR>& x,
  const DistMatrix<float,MC,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VR>& x,
  const DistMatrix<float,MC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VR>& x,
  const DistMatrix<float,STAR,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VR>& x,
  const DistMatrix<float,MR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VR>& x,
  const DistMatrix<float,MR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VR>& x,
  const DistMatrix<float,STAR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VR>& x,
  const DistMatrix<float,VC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VR>& x,
  const DistMatrix<float,STAR,VC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VR>& x,
  const DistMatrix<float,VR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VR>& x,
  const DistMatrix<float,STAR,VR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,VR>& x,
  const DistMatrix<float,STAR,STAR>& y );

template float elemental::basic::Dot
( const DistMatrix<float,STAR,STAR>& x,
  const DistMatrix<float,MC,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,STAR>& x,
  const DistMatrix<float,MC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,STAR>& x,
  const DistMatrix<float,STAR,MR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,STAR>& x,
  const DistMatrix<float,MR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,STAR>& x,
  const DistMatrix<float,MR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,STAR>& x,
  const DistMatrix<float,STAR,MC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,STAR>& x,
  const DistMatrix<float,VC,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,STAR>& x,
  const DistMatrix<float,STAR,VC>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,STAR>& x,
  const DistMatrix<float,VR,STAR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,STAR>& x,
  const DistMatrix<float,STAR,VR>& y );
template float elemental::basic::Dot
( const DistMatrix<float,STAR,STAR>& x,
  const DistMatrix<float,STAR,STAR>& y );
#endif // ENABLE_ALL_DISTRIBUTED_DOT

template double elemental::basic::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,MC,MR>& y );
#ifdef ENABLE_ALL_DISTRIBUTED_DOT
template double elemental::basic::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,MC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,STAR,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,MR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,MR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,STAR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,VC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,STAR,VC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,VR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,STAR,VR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,MR>& x,
  const DistMatrix<double,STAR,STAR>& y );

template double elemental::basic::Dot
( const DistMatrix<double,MC,STAR>& x,
  const DistMatrix<double,MC,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,STAR>& x,
  const DistMatrix<double,MC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,STAR>& x,
  const DistMatrix<double,STAR,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,STAR>& x,
  const DistMatrix<double,MR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,STAR>& x,
  const DistMatrix<double,MR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,STAR>& x,
  const DistMatrix<double,STAR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,STAR>& x,
  const DistMatrix<double,VC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,STAR>& x,
  const DistMatrix<double,STAR,VC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,STAR>& x,
  const DistMatrix<double,VR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,STAR>& x,
  const DistMatrix<double,STAR,VR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MC,STAR>& x,
  const DistMatrix<double,STAR,STAR>& y );

template double elemental::basic::Dot
( const DistMatrix<double,STAR,MR>& x,
  const DistMatrix<double,MC,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MR>& x,
  const DistMatrix<double,MC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MR>& x,
  const DistMatrix<double,STAR,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MR>& x,
  const DistMatrix<double,MR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MR>& x,
  const DistMatrix<double,MR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MR>& x,
  const DistMatrix<double,STAR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MR>& x,
  const DistMatrix<double,VC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MR>& x,
  const DistMatrix<double,STAR,VC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MR>& x,
  const DistMatrix<double,VR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MR>& x,
  const DistMatrix<double,STAR,VR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MR>& x,
  const DistMatrix<double,STAR,STAR>& y );

template double elemental::basic::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,MC,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,MC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,STAR,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,MR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,MR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,STAR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,VC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,STAR,VC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,VR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,STAR,VR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,MC>& x,
  const DistMatrix<double,STAR,STAR>& y );

template double elemental::basic::Dot
( const DistMatrix<double,MR,STAR>& x,
  const DistMatrix<double,MC,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,STAR>& x,
  const DistMatrix<double,MC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,STAR>& x,
  const DistMatrix<double,STAR,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,STAR>& x,
  const DistMatrix<double,MR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,STAR>& x,
  const DistMatrix<double,MR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,STAR>& x,
  const DistMatrix<double,STAR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,STAR>& x,
  const DistMatrix<double,VC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,STAR>& x,
  const DistMatrix<double,STAR,VC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,STAR>& x,
  const DistMatrix<double,VR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,STAR>& x,
  const DistMatrix<double,STAR,VR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,MR,STAR>& x,
  const DistMatrix<double,STAR,STAR>& y );

template double elemental::basic::Dot
( const DistMatrix<double,STAR,MC>& x,
  const DistMatrix<double,MC,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MC>& x,
  const DistMatrix<double,MC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MC>& x,
  const DistMatrix<double,STAR,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MC>& x,
  const DistMatrix<double,MR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MC>& x,
  const DistMatrix<double,MR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MC>& x,
  const DistMatrix<double,STAR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MC>& x,
  const DistMatrix<double,VC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MC>& x,
  const DistMatrix<double,STAR,VC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MC>& x,
  const DistMatrix<double,VR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MC>& x,
  const DistMatrix<double,STAR,VR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,MC>& x,
  const DistMatrix<double,STAR,STAR>& y );

template double elemental::basic::Dot
( const DistMatrix<double,VC,STAR>& x,
  const DistMatrix<double,MC,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VC,STAR>& x,
  const DistMatrix<double,MC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VC,STAR>& x,
  const DistMatrix<double,STAR,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VC,STAR>& x,
  const DistMatrix<double,MR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VC,STAR>& x,
  const DistMatrix<double,MR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VC,STAR>& x,
  const DistMatrix<double,STAR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VC,STAR>& x,
  const DistMatrix<double,VC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VC,STAR>& x,
  const DistMatrix<double,STAR,VC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VC,STAR>& x,
  const DistMatrix<double,VR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VC,STAR>& x,
  const DistMatrix<double,STAR,VR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VC,STAR>& x,
  const DistMatrix<double,STAR,STAR>& y );

template double elemental::basic::Dot
( const DistMatrix<double,STAR,VC>& x,
  const DistMatrix<double,MC,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VC>& x,
  const DistMatrix<double,MC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VC>& x,
  const DistMatrix<double,STAR,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VC>& x,
  const DistMatrix<double,MR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VC>& x,
  const DistMatrix<double,MR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VC>& x,
  const DistMatrix<double,STAR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VC>& x,
  const DistMatrix<double,VC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VC>& x,
  const DistMatrix<double,STAR,VC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VC>& x,
  const DistMatrix<double,VR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VC>& x,
  const DistMatrix<double,STAR,VR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VC>& x,
  const DistMatrix<double,STAR,STAR>& y );

template double elemental::basic::Dot
( const DistMatrix<double,VR,STAR>& x,
  const DistMatrix<double,MC,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VR,STAR>& x,
  const DistMatrix<double,MC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VR,STAR>& x,
  const DistMatrix<double,STAR,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VR,STAR>& x,
  const DistMatrix<double,MR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VR,STAR>& x,
  const DistMatrix<double,MR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VR,STAR>& x,
  const DistMatrix<double,STAR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VR,STAR>& x,
  const DistMatrix<double,VC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VR,STAR>& x,
  const DistMatrix<double,STAR,VC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VR,STAR>& x,
  const DistMatrix<double,VR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VR,STAR>& x,
  const DistMatrix<double,STAR,VR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,VR,STAR>& x,
  const DistMatrix<double,STAR,STAR>& y );

template double elemental::basic::Dot
( const DistMatrix<double,STAR,VR>& x,
  const DistMatrix<double,MC,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VR>& x,
  const DistMatrix<double,MC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VR>& x,
  const DistMatrix<double,STAR,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VR>& x,
  const DistMatrix<double,MR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VR>& x,
  const DistMatrix<double,MR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VR>& x,
  const DistMatrix<double,STAR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VR>& x,
  const DistMatrix<double,VC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VR>& x,
  const DistMatrix<double,STAR,VC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VR>& x,
  const DistMatrix<double,VR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VR>& x,
  const DistMatrix<double,STAR,VR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,VR>& x,
  const DistMatrix<double,STAR,STAR>& y );

template double elemental::basic::Dot
( const DistMatrix<double,STAR,STAR>& x,
  const DistMatrix<double,MC,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,STAR>& x,
  const DistMatrix<double,MC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,STAR>& x,
  const DistMatrix<double,STAR,MR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,STAR>& x,
  const DistMatrix<double,MR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,STAR>& x,
  const DistMatrix<double,MR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,STAR>& x,
  const DistMatrix<double,STAR,MC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,STAR>& x,
  const DistMatrix<double,VC,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,STAR>& x,
  const DistMatrix<double,STAR,VC>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,STAR>& x,
  const DistMatrix<double,VR,STAR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,STAR>& x,
  const DistMatrix<double,STAR,VR>& y );
template double elemental::basic::Dot
( const DistMatrix<double,STAR,STAR>& x,
  const DistMatrix<double,STAR,STAR>& y );
#endif // ENABLE_ALL_DISTRIBUTED_DOT

#ifndef WITHOUT_COMPLEX
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,MC,MR>& y );
#ifdef ENABLE_ALL_DISTRIBUTED_DOT
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,MC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,STAR,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,MR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,STAR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,VC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,STAR,VC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,VR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,STAR,VR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,MR>& x,
  const DistMatrix<scomplex,STAR,STAR>& y );

template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,STAR>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,STAR>& x,
  const DistMatrix<scomplex,MC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,STAR>& x,
  const DistMatrix<scomplex,STAR,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,STAR>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,STAR>& x,
  const DistMatrix<scomplex,MR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,STAR>& x,
  const DistMatrix<scomplex,STAR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,STAR>& x,
  const DistMatrix<scomplex,VC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,STAR>& x,
  const DistMatrix<scomplex,STAR,VC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,STAR>& x,
  const DistMatrix<scomplex,VR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,STAR>& x,
  const DistMatrix<scomplex,STAR,VR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MC,STAR>& x,
  const DistMatrix<scomplex,STAR,STAR>& y );

template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MR>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MR>& x,
  const DistMatrix<scomplex,MC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MR>& x,
  const DistMatrix<scomplex,STAR,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MR>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MR>& x,
  const DistMatrix<scomplex,MR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MR>& x,
  const DistMatrix<scomplex,STAR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MR>& x,
  const DistMatrix<scomplex,VC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MR>& x,
  const DistMatrix<scomplex,STAR,VC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MR>& x,
  const DistMatrix<scomplex,VR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MR>& x,
  const DistMatrix<scomplex,STAR,VR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MR>& x,
  const DistMatrix<scomplex,STAR,STAR>& y );

template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,MC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,STAR,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,MR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,STAR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,VC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,STAR,VC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,VR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,STAR,VR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,MC>& x,
  const DistMatrix<scomplex,STAR,STAR>& y );

template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,STAR>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,STAR>& x,
  const DistMatrix<scomplex,MC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,STAR>& x,
  const DistMatrix<scomplex,STAR,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,STAR>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,STAR>& x,
  const DistMatrix<scomplex,MR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,STAR>& x,
  const DistMatrix<scomplex,STAR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,STAR>& x,
  const DistMatrix<scomplex,VC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,STAR>& x,
  const DistMatrix<scomplex,STAR,VC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,STAR>& x,
  const DistMatrix<scomplex,VR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,STAR>& x,
  const DistMatrix<scomplex,STAR,VR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,MR,STAR>& x,
  const DistMatrix<scomplex,STAR,STAR>& y );

template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MC>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MC>& x,
  const DistMatrix<scomplex,MC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MC>& x,
  const DistMatrix<scomplex,STAR,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MC>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MC>& x,
  const DistMatrix<scomplex,MR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MC>& x,
  const DistMatrix<scomplex,STAR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MC>& x,
  const DistMatrix<scomplex,VC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MC>& x,
  const DistMatrix<scomplex,STAR,VC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MC>& x,
  const DistMatrix<scomplex,VR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MC>& x,
  const DistMatrix<scomplex,STAR,VR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,MC>& x,
  const DistMatrix<scomplex,STAR,STAR>& y );

template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VC,STAR>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VC,STAR>& x,
  const DistMatrix<scomplex,MC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VC,STAR>& x,
  const DistMatrix<scomplex,STAR,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VC,STAR>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VC,STAR>& x,
  const DistMatrix<scomplex,MR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VC,STAR>& x,
  const DistMatrix<scomplex,STAR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VC,STAR>& x,
  const DistMatrix<scomplex,VC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VC,STAR>& x,
  const DistMatrix<scomplex,STAR,VC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VC,STAR>& x,
  const DistMatrix<scomplex,VR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VC,STAR>& x,
  const DistMatrix<scomplex,STAR,VR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VC,STAR>& x,
  const DistMatrix<scomplex,STAR,STAR>& y );

template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VC>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VC>& x,
  const DistMatrix<scomplex,MC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VC>& x,
  const DistMatrix<scomplex,STAR,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VC>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VC>& x,
  const DistMatrix<scomplex,MR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VC>& x,
  const DistMatrix<scomplex,STAR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VC>& x,
  const DistMatrix<scomplex,VC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VC>& x,
  const DistMatrix<scomplex,STAR,VC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VC>& x,
  const DistMatrix<scomplex,VR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VC>& x,
  const DistMatrix<scomplex,STAR,VR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VC>& x,
  const DistMatrix<scomplex,STAR,STAR>& y );

template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VR,STAR>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VR,STAR>& x,
  const DistMatrix<scomplex,MC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VR,STAR>& x,
  const DistMatrix<scomplex,STAR,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VR,STAR>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VR,STAR>& x,
  const DistMatrix<scomplex,MR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VR,STAR>& x,
  const DistMatrix<scomplex,STAR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VR,STAR>& x,
  const DistMatrix<scomplex,VC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VR,STAR>& x,
  const DistMatrix<scomplex,STAR,VC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VR,STAR>& x,
  const DistMatrix<scomplex,VR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VR,STAR>& x,
  const DistMatrix<scomplex,STAR,VR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,VR,STAR>& x,
  const DistMatrix<scomplex,STAR,STAR>& y );

template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VR>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VR>& x,
  const DistMatrix<scomplex,MC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VR>& x,
  const DistMatrix<scomplex,STAR,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VR>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VR>& x,
  const DistMatrix<scomplex,MR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VR>& x,
  const DistMatrix<scomplex,STAR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VR>& x,
  const DistMatrix<scomplex,VC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VR>& x,
  const DistMatrix<scomplex,STAR,VC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VR>& x,
  const DistMatrix<scomplex,VR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VR>& x,
  const DistMatrix<scomplex,STAR,VR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,VR>& x,
  const DistMatrix<scomplex,STAR,STAR>& y );

template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,STAR>& x,
  const DistMatrix<scomplex,MC,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,STAR>& x,
  const DistMatrix<scomplex,MC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,STAR>& x,
  const DistMatrix<scomplex,STAR,MR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,STAR>& x,
  const DistMatrix<scomplex,MR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,STAR>& x,
  const DistMatrix<scomplex,MR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,STAR>& x,
  const DistMatrix<scomplex,STAR,MC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,STAR>& x,
  const DistMatrix<scomplex,VC,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,STAR>& x,
  const DistMatrix<scomplex,STAR,VC>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,STAR>& x,
  const DistMatrix<scomplex,VR,STAR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,STAR>& x,
  const DistMatrix<scomplex,STAR,VR>& y );
template scomplex elemental::basic::Dot
( const DistMatrix<scomplex,STAR,STAR>& x,
  const DistMatrix<scomplex,STAR,STAR>& y );
#endif // ENABLE_ALL_DISTRIBUTED_DOT

template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
#ifdef ENABLE_ALL_DISTRIBUTED_DOT
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,MC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,STAR,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,MR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,STAR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,VC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,STAR,VC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,VR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,STAR,VR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,MR>& x,
  const DistMatrix<dcomplex,STAR,STAR>& y );

template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,STAR>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,STAR>& x,
  const DistMatrix<dcomplex,MC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,STAR>& x,
  const DistMatrix<dcomplex,STAR,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,STAR>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,STAR>& x,
  const DistMatrix<dcomplex,MR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,STAR>& x,
  const DistMatrix<dcomplex,STAR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,STAR>& x,
  const DistMatrix<dcomplex,VC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,STAR>& x,
  const DistMatrix<dcomplex,STAR,VC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,STAR>& x,
  const DistMatrix<dcomplex,VR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,STAR>& x,
  const DistMatrix<dcomplex,STAR,VR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MC,STAR>& x,
  const DistMatrix<dcomplex,STAR,STAR>& y );

template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MR>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MR>& x,
  const DistMatrix<dcomplex,MC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MR>& x,
  const DistMatrix<dcomplex,STAR,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MR>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MR>& x,
  const DistMatrix<dcomplex,MR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MR>& x,
  const DistMatrix<dcomplex,STAR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MR>& x,
  const DistMatrix<dcomplex,VC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MR>& x,
  const DistMatrix<dcomplex,STAR,VC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MR>& x,
  const DistMatrix<dcomplex,VR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MR>& x,
  const DistMatrix<dcomplex,STAR,VR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MR>& x,
  const DistMatrix<dcomplex,STAR,STAR>& y );

template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,MC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,STAR,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,MR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,STAR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,VC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,STAR,VC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,VR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,STAR,VR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,MC>& x,
  const DistMatrix<dcomplex,STAR,STAR>& y );

template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,STAR>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,STAR>& x,
  const DistMatrix<dcomplex,MC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,STAR>& x,
  const DistMatrix<dcomplex,STAR,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,STAR>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,STAR>& x,
  const DistMatrix<dcomplex,MR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,STAR>& x,
  const DistMatrix<dcomplex,STAR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,STAR>& x,
  const DistMatrix<dcomplex,VC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,STAR>& x,
  const DistMatrix<dcomplex,STAR,VC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,STAR>& x,
  const DistMatrix<dcomplex,VR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,STAR>& x,
  const DistMatrix<dcomplex,STAR,VR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,MR,STAR>& x,
  const DistMatrix<dcomplex,STAR,STAR>& y );

template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MC>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MC>& x,
  const DistMatrix<dcomplex,MC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MC>& x,
  const DistMatrix<dcomplex,STAR,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MC>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MC>& x,
  const DistMatrix<dcomplex,MR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MC>& x,
  const DistMatrix<dcomplex,STAR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MC>& x,
  const DistMatrix<dcomplex,VC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MC>& x,
  const DistMatrix<dcomplex,STAR,VC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MC>& x,
  const DistMatrix<dcomplex,VR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MC>& x,
  const DistMatrix<dcomplex,STAR,VR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,MC>& x,
  const DistMatrix<dcomplex,STAR,STAR>& y );

template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VC,STAR>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VC,STAR>& x,
  const DistMatrix<dcomplex,MC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VC,STAR>& x,
  const DistMatrix<dcomplex,STAR,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VC,STAR>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VC,STAR>& x,
  const DistMatrix<dcomplex,MR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VC,STAR>& x,
  const DistMatrix<dcomplex,STAR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VC,STAR>& x,
  const DistMatrix<dcomplex,VC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VC,STAR>& x,
  const DistMatrix<dcomplex,STAR,VC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VC,STAR>& x,
  const DistMatrix<dcomplex,VR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VC,STAR>& x,
  const DistMatrix<dcomplex,STAR,VR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VC,STAR>& x,
  const DistMatrix<dcomplex,STAR,STAR>& y );

template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VC>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VC>& x,
  const DistMatrix<dcomplex,MC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VC>& x,
  const DistMatrix<dcomplex,STAR,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VC>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VC>& x,
  const DistMatrix<dcomplex,MR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VC>& x,
  const DistMatrix<dcomplex,STAR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VC>& x,
  const DistMatrix<dcomplex,VC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VC>& x,
  const DistMatrix<dcomplex,STAR,VC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VC>& x,
  const DistMatrix<dcomplex,VR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VC>& x,
  const DistMatrix<dcomplex,STAR,VR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VC>& x,
  const DistMatrix<dcomplex,STAR,STAR>& y );

template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VR,STAR>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VR,STAR>& x,
  const DistMatrix<dcomplex,MC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VR,STAR>& x,
  const DistMatrix<dcomplex,STAR,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VR,STAR>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VR,STAR>& x,
  const DistMatrix<dcomplex,MR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VR,STAR>& x,
  const DistMatrix<dcomplex,STAR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VR,STAR>& x,
  const DistMatrix<dcomplex,VC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VR,STAR>& x,
  const DistMatrix<dcomplex,STAR,VC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VR,STAR>& x,
  const DistMatrix<dcomplex,VR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VR,STAR>& x,
  const DistMatrix<dcomplex,STAR,VR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,VR,STAR>& x,
  const DistMatrix<dcomplex,STAR,STAR>& y );

template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VR>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VR>& x,
  const DistMatrix<dcomplex,MC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VR>& x,
  const DistMatrix<dcomplex,STAR,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VR>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VR>& x,
  const DistMatrix<dcomplex,MR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VR>& x,
  const DistMatrix<dcomplex,STAR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VR>& x,
  const DistMatrix<dcomplex,VC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VR>& x,
  const DistMatrix<dcomplex,STAR,VC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VR>& x,
  const DistMatrix<dcomplex,VR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VR>& x,
  const DistMatrix<dcomplex,STAR,VR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,VR>& x,
  const DistMatrix<dcomplex,STAR,STAR>& y );

template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,STAR>& x,
  const DistMatrix<dcomplex,MC,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,STAR>& x,
  const DistMatrix<dcomplex,MC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,STAR>& x,
  const DistMatrix<dcomplex,STAR,MR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,STAR>& x,
  const DistMatrix<dcomplex,MR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,STAR>& x,
  const DistMatrix<dcomplex,MR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,STAR>& x,
  const DistMatrix<dcomplex,STAR,MC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,STAR>& x,
  const DistMatrix<dcomplex,VC,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,STAR>& x,
  const DistMatrix<dcomplex,STAR,VC>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,STAR>& x,
  const DistMatrix<dcomplex,VR,STAR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,STAR>& x,
  const DistMatrix<dcomplex,STAR,VR>& y );
template dcomplex elemental::basic::Dot
( const DistMatrix<dcomplex,STAR,STAR>& x,
  const DistMatrix<dcomplex,STAR,STAR>& y );
#endif // ENABLE_ALL_DISTRIBUTED_DOT
#endif // WITHOUT_COMPLEX
