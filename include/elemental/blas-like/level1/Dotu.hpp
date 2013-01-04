/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename T>
inline T
Dotu( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Dotu");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Expected vector inputs");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("x and y must be the same length");
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = blas::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = blas::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = blas::Dotu
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = blas::Dotu
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

namespace internal {

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
     ( const DistMatrix<T,U,V>& x, const DistMatrix<T>& y );
*/
template<typename T,Distribution U,Distribution V>
inline T
DotuHelper( const DistMatrix<T,U,V>& x, const DistMatrix<T>& y )
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
        DistMatrix<T> xRedist(g);
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
        DistMatrix<T> xRedist(g);
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
        DistMatrix<T> xRedist(g);
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
        DistMatrix<T> xRedist(g);
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
