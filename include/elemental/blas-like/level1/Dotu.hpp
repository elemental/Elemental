/*
   Copyright (c) 2009-2012, Jack Poulson
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

/* 
   C++ does not (currently) allow for partial function template specialization,
   so implementing Dotu will be unnecessarily obfuscated. Sorry.
  
   The compromise is to have the user-level routine, 
  
     template<typename T,Distribution U,Distribution V,
                         Distribution W,Distribution Z>
     T Dotu( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

   simply route to overloaded pseudo-partial-specializations within the 
   internal namespace. For example,

     template<typename T,Distribution U,Distribution V>
     T internal::DotuHelper
     ( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,MR>& y );
*/
namespace elem {

namespace internal {

template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("internal::DotuHelper");
    if( x.Grid() != y.Grid() )
        throw std::logic_error("{x,y} must be distributed over the same grid");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Dotu requires x and y to be vectors");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("Dotu requires x and y to be the same length");
#endif
    const Grid& g = x.Grid();

    T globalDotu;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( g.Col() == ownerCol )
        { 
            T localDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.ColComm() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerCol, g.RowComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( g.Row() == ownerRow )
        {
            T localDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.RowComm() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerRow, g.ColComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( g.Col() == ownerCol )
        {
            T localDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.ColComm() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerCol, g.RowComm() );
    }
    else
    {
        DistMatrix<T,MC,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( g.Row() == ownerRow )
        {
            T localDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.RowComm() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerRow, g.ColComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDotu;
}

template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("internal::DotuHelper");
    if( x.Grid() != y.Grid() )
        throw std::logic_error("{x,y} must be distributed over the same grid");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Dotu requires x and y to be vectors");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("Dotu requires x and y to be the same length");
#endif
    const Grid& g = x.Grid();

    T globalDotu;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.ColComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,STAR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( g.Row() == ownerRow )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerRow, g.ColComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.ColComm() );
    }
    else
    {
        DistMatrix<T,MC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.ColAlignment();
        if( g.Row() == ownerRow )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerRow, g.ColComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDotu;
}

template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T,STAR,MR>& y )
{
#ifndef RELEASE
    PushCallStack("internal::DotuHelper");
    if( x.Grid() != y.Grid() )
        throw std::logic_error("{x,y} must be distributed over the same grid");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Dotu requires x and y to be vectors");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("Dotu requires x and y to be the same length");
#endif
    const Grid& g = x.Grid();

    T globalDotu;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,STAR,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( g.Col() == ownerCol )
        { 
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerCol, g.RowComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.RowComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.RowAlignment();
        if( g.Col() == ownerCol )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerCol, g.RowComm() );
    }
    else
    {
        DistMatrix<T,STAR,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.RowComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDotu;
}

template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T,MR,MC>& y )
{
#ifndef RELEASE
    PushCallStack("internal::DotuHelper");
    if( x.Grid() != y.Grid() )
        throw std::logic_error("{x,y} must be distributed over the same grid");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Dotu requires x and y to be vectors");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("Dotu requires x and y to be the same length");
#endif
    const Grid& g = x.Grid();

    T globalDotu;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( g.Row() == ownerRow )
        { 
            T localDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.RowComm() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerRow, g.ColComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( g.Col() == ownerCol )
        {
            T localDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.ColComm() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerCol, g.RowComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MC,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( g.Row() == ownerRow )
        {
            T localDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.RowComm() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerRow, g.ColComm() );
    }
    else
    {
        DistMatrix<T,MR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( g.Col() == ownerCol )
        {
            T localDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
            mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.ColComm() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerCol, g.RowComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDotu;
}

template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T,MR,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("internal::DotuHelper");
    if( x.Grid() != y.Grid() )
        throw std::logic_error("{x,y} must be distributed over the same grid");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Dotu requires x and y to be vectors");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("Dotu requires x and y to be the same length");
#endif
    const Grid& g = x.Grid();

    T globalDotu;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.RowComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,STAR,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( g.Col() == ownerCol )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerCol, g.RowComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.RowComm() );
    }
    else
    {
        DistMatrix<T,MR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerCol = y.ColAlignment();
        if( g.Col() == ownerCol )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerCol, g.RowComm() );
    }
#ifndef RELEASE 
    PopCallStack();
#endif
    return globalDotu;
}

template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T,STAR,MC>& y )
{
#ifndef RELEASE
    PushCallStack("internal::DotuHelper");
    if( x.Grid() != y.Grid() )
        throw std::logic_error("{x,y} must be distributed over the same grid");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Dotu requires x and y to be vectors");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("Dotu requires x and y to be the same length");
#endif
    const Grid& g = x.Grid();

    T globalDotu;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,STAR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( g.Row() == ownerRow )
        { 
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerRow, g.ColComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.ColComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int ownerRow = y.RowAlignment();
        if( g.Row() == ownerRow )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, ownerRow, g.ColComm() );
    }
    else
    {
        DistMatrix<T,STAR,MC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.ColComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDotu;
}

template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T,VC,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("internal::DotuHelper");
    if( x.Grid() != y.Grid() )
        throw std::logic_error("{x,y} must be distributed over the same grid");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Dotu requires x and y to be vectors");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("Dotu requires x and y to be the same length");
#endif
    const Grid& g = x.Grid();

    T globalDotu;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,VC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.VCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,STAR,VC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( g.VCRank() == owner )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, owner, g.VCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,VC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.VCComm() );
    }
    else
    {
        DistMatrix<T,VC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( g.VCRank() == owner )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, owner, g.VCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDotu;
}

template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T,STAR,VC>& y )
{
#ifndef RELEASE
    PushCallStack("internal::DotuHelper");
    if( x.Grid() != y.Grid() )
        throw std::logic_error("{x,y} must be distributed over the same grid");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Dotu requires x and y to be vectors");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("Dotu requires x and y to be the same length");
#endif
    const Grid& g = x.Grid();

    T globalDotu;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,STAR,VC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( g.VCRank() == owner )
        { 
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, owner, g.VCComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,VC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.VCComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,VC,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( g.VCRank() == owner )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, owner, g.VCComm() );
    }
    else
    {
        DistMatrix<T,STAR,VC> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.VCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDotu;
}

template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T,VR,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("internal::DotuHelper");
    if( x.Grid() != y.Grid() )
        throw std::logic_error("{x,y} must be distributed over the same grid");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Dotu requires x and y to be vectors");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("Dotu requires x and y to be the same length");
#endif
    const Grid& g = x.Grid();

    T globalDotu;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,VR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.VRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,STAR,VR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( g.VRRank() == owner )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, owner, g.VRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,VR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.VRComm() );
    }
    else
    {
        DistMatrix<T,VR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.ColAlignment();
        if( g.VRRank() == owner )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, owner, g.VRComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDotu;
}

template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T,STAR,VR>& y )
{
#ifndef RELEASE
    PushCallStack("internal::DotuHelper");
    if( x.Grid() != y.Grid() )
        throw std::logic_error("{x,y} must be distributed over the same grid");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Dotu requires x and y to be vectors");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("Dotu requires x and y to be the same length");
#endif
    const Grid& g = x.Grid();

    T globalDotu;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,STAR,VR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( g.VRRank() == owner )
        { 
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, owner, g.VRComm() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,VR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.VRComm() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,VR,STAR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        int owner = y.RowAlignment();
        if( g.VRRank() == owner )
        {
            globalDotu = 
                Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        }
        mpi::Broadcast( &globalDotu, 1, owner, g.VRComm() );
    }
    else
    {
        DistMatrix<T,STAR,VR> xRedist(g);
        xRedist.AlignWith( y );
        xRedist = x;

        T localDotu = 
            Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
        mpi::AllReduce( &localDotu, &globalDotu, 1, mpi::SUM, g.VRComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDotu;
}

template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T,STAR,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("internal::DotuHelper");
    if( x.Grid() != y.Grid() )
        throw std::logic_error("{x,y} must be distributed over the same grid");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Dotu requires x and y to be vectors");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("Dotu requires x and y to be the same length");
#endif
    const Grid& g = x.Grid();

    DistMatrix<T,STAR,STAR> xRedist(g);
    xRedist = x;

    T globalDotu = 
        Dotu( xRedist.LockedLocalMatrix(), y.LockedLocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDotu;
}

} // namespace internal

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline T
Dotu( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )
{
#ifndef RELEASE
    PushCallStack("Dotu");
#endif
    T dotProduct = internal::DotuHelper( x, y );
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

} // namespace elem
