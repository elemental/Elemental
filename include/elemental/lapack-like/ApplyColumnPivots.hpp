/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_APPLYCOLUMNPIVOTS_HPP
#define LAPACK_APPLYCOLUMNPIVOTS_HPP

#include "elemental/lapack-like/ComposePivots.hpp"

namespace elem {

template<typename F>
inline void
ApplyColumnPivots( Matrix<F>& A, const Matrix<int>& p )
{
#ifndef RELEASE
    PushCallStack("ApplyColumnPivots");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
    if( p.Height() != A.Width() )
        throw std::logic_error("p must be the same length as the width of A");
#endif
    const int height = A.Height();
    const int width = A.Width();
    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    for( int j=0; j<width; ++j )
    {
        const int k = p.Get(j,0);
        F* Aj = A.Buffer(0,j);
        F* Ak = A.Buffer(0,k);
        for( int i=0; i<height; ++i )
        {
            F temp = Aj[i];
            Aj[i] = Ak[i];
            Ak[i] = temp;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
ApplyInverseColumnPivots( Matrix<F>& A, const Matrix<int>& p )
{
#ifndef RELEASE
    PushCallStack("ApplyInverseColumnPivots");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
    if( p.Height() != A.Width() )
        throw std::logic_error("p must be the same length as the width of A");
#endif
    const int height = A.Height();
    const int width = A.Width();
    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    for( int j=width-1; j>=0; --j )
    {
        const int k = p.Get(j,0);
        F* Aj = A.Buffer(0,j);
        F* Ak = A.Buffer(0,k);
        for( int i=0; i<height; ++i )
        {
            F temp = Aj[i];
            Aj[i] = Ak[i];
            Ak[i] = temp;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
ApplyColumnPivots( DistMatrix<F>& A, const DistMatrix<int,VC,STAR>& p )
{
#ifndef RELEASE
    PushCallStack("ApplyColumnPivots");
#endif
    DistMatrix<int,STAR,STAR> p_STAR_STAR( p );
    ApplyColumnPivots( A, p_STAR_STAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
ApplyInverseColumnPivots
( DistMatrix<F>& A, const DistMatrix<int,VC,STAR>& p )
{
#ifndef RELEASE
    PushCallStack("ApplyInverseColumnPivots");
#endif
    DistMatrix<int,STAR,STAR> p_STAR_STAR( p );
    ApplyInverseColumnPivots( A, p_STAR_STAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
ApplyColumnPivots( DistMatrix<F>& A, const DistMatrix<int,STAR,STAR>& p )
{
#ifndef RELEASE
    PushCallStack("ApplyColumnPivots");
#endif
    std::vector<int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyColumnPivots( A, image, preimage );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
ApplyInverseColumnPivots
( DistMatrix<F>& A, const DistMatrix<int,STAR,STAR>& p )
{
#ifndef RELEASE
    PushCallStack("ApplyInverseColumnPivots");
#endif
    std::vector<int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyColumnPivots( A, preimage, image );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
ApplyColumnPivots
( Matrix<F>& A, 
  const std::vector<int>& image,
  const std::vector<int>& preimage )
{
    const int b = image.size();
#ifndef RELEASE
    PushCallStack("ApplyColumnPivots");
    if( A.Width() < b || b != preimage.size() )
        throw std::logic_error
        ("image and preimage must be vectors of equal length that are not "
         "wider than A.");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( m == 0 || n == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // TODO: Optimize this routine

    // Make a copy of the first b columns
    Matrix<F> AColPanView;
    LockedView( AColPanView, A, 0, 0, m, b );
    Matrix<F> AColPanCopy( AColPanView );

    // Make a copy of the preimage columns
    Matrix<F> APreimageCopy( m, b );
    for( int j=0; j<b; ++j )
    {
        const int jPre = preimage[j];
        if( jPre >= b )
            MemCopy( APreimageCopy.Buffer(0,j), A.LockedBuffer(0,jPre), m );
    }

    // Apply the permutations
    for( int j=0; j<b; ++j )
    {
        const int jPre = preimage[j];
        const int jPost = image[j];
        // Move row[i] into row[image[i]]
        MemCopy( A.Buffer(0,jPost), AColPanCopy.LockedBuffer(0,j), m );
        // Move row[preimage[i]] into row[i]
        if( jPre >= b )
            MemCopy( A.Buffer(0,j), APreimageCopy.LockedBuffer(0,j), m );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
ApplyColumnPivots
( DistMatrix<F>& A, 
  const std::vector<int>& image,
  const std::vector<int>& preimage )
{
    const int b = image.size();
#ifndef RELEASE
    PushCallStack("ApplyColumnPivots");
    if( A.Width() < b || b != preimage.size() )
        throw std::logic_error
        ("image and preimage must be vectors of equal length that are not "
         "wider than A.");
#endif
    const int localHeight = A.LocalHeight();
    if( A.Height() == 0 || A.Width() == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Extract the relevant process grid information
    const Grid& g = A.Grid();
    const int c = g.Width();
    const int rowAlignment = A.RowAlignment();
    const int rowShift = A.RowShift();
    const int myCol = g.Col();

    // Extract the send and recv counts from the image and preimage.
    // This process's sends may be logically partitioned into two sets:
    //   (a) sends from rows [0,...,b-1]
    //   (b) sends from rows [b,...]
    // The latter is analyzed with preimage, the former deduced with image.
    std::vector<int> sendCounts(c,0), recvCounts(c,0);
    for( int j=rowShift; j<b; j+=c )
    {
        const int sendCol = image[j];         
        const int sendTo = (rowAlignment+sendCol) % c; 
        sendCounts[sendTo] += localHeight;

        const int recvCol = preimage[j];
        const int recvFrom = (rowAlignment+recvCol) % c;
        recvCounts[recvFrom] += localHeight;
    }
    for( int j=0; j<b; ++j )
    {
        const int sendCol = image[j];
        if( sendCol >= b )
        {
            const int sendTo = (rowAlignment+sendCol) % c;
            if( sendTo == myCol )
            {
                const int sendFrom = (rowAlignment+j) % c;
                recvCounts[sendFrom] += localHeight;
            }
        }

        const int recvCol = preimage[j];
        if( recvCol >= b )
        {
            const int recvFrom = (rowAlignment+recvCol) % c;
            if( recvFrom == myCol )
            {
                const int recvTo = (rowAlignment+j) % c;
                sendCounts[recvTo] += localHeight;
            }
        }
    }

    // Construct the send and recv displacements from the counts
    std::vector<int> sendDispls(c), recvDispls(c);
    int totalSend=0, totalRecv=0;
    for( int i=0; i<c; ++i )
    {
        sendDispls[i] = totalSend;
        recvDispls[i] = totalRecv;
        totalSend += sendCounts[i];
        totalRecv += recvCounts[i];
    }
#ifndef RELEASE
    if( totalSend != totalRecv )
    {
        std::ostringstream msg;
        msg << "Send and recv counts do not match: (send,recv)=" 
             << totalSend << "," << totalRecv;
        throw std::logic_error( msg.str().c_str() );
    }
#endif

    // Fill vectors with the send data
    std::vector<F> sendData(std::max(1,totalSend));
    std::vector<int> offsets(c,0);
    const int localWidth = LocalLength( b, rowShift, c );
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int sendCol = image[rowShift+jLocal*c];
        const int sendTo = (rowAlignment+sendCol) % c;
        const int offset = sendDispls[sendTo]+offsets[sendTo];
        MemCopy( &sendData[offset], A.LocalBuffer(0,jLocal), localHeight );
        offsets[sendTo] += localHeight;
    }
    for( int j=0; j<b; ++j )
    {
        const int recvCol = preimage[j];
        if( recvCol >= b )
        {
            const int recvFrom = (rowAlignment+recvCol) % c; 
            if( recvFrom == myCol )
            {
                const int recvTo = (rowAlignment+j) % c;
                const int jLocal = (recvCol-rowShift) / c;
                const int offset = sendDispls[recvTo]+offsets[recvTo];
                MemCopy
                ( &sendData[offset], A.LocalBuffer(0,jLocal), localHeight );
                offsets[recvTo] += localHeight;
            }
        }
    }

    // Communicate all pivot rows
    std::vector<F> recvData(std::max(1,totalRecv));
    mpi::AllToAll
    ( &sendData[0], &sendCounts[0], &sendDispls[0],
      &recvData[0], &recvCounts[0], &recvDispls[0], g.RowComm() );

    // Unpack the recv data
    for( int k=0; k<c; ++k )
    {
        offsets[k] = 0;
        int thisRowShift = Shift( k, rowAlignment, c );
        for( int j=thisRowShift; j<b; j+=c )
        {
            const int sendCol = image[j];
            const int sendTo = (rowAlignment+sendCol) % c;
            if( sendTo == myCol )
            {
                const int offset = recvDispls[k]+offsets[k];
                const int jLocal = (sendCol-rowShift) / c;
                MemCopy
                ( A.LocalBuffer(0,jLocal), &recvData[offset], localHeight );
                offsets[k] += localHeight;
            }
        }
    }
    for( int j=0; j<b; ++j )
    {
        const int recvCol = preimage[j];
        if( recvCol >= b )
        {
            const int recvTo = (rowAlignment+j) % c;
            if( recvTo == myCol )
            {
                const int recvFrom = (rowAlignment+recvCol) % c; 
                const int jLocal = (j-rowShift) / c;
                const int offset = recvDispls[recvFrom]+offsets[recvFrom];
                MemCopy
                ( A.LocalBuffer(0,jLocal), &recvData[offset], localHeight );
                offsets[recvFrom] += localHeight;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_APPLYCOLUMNPIVOTS_HPP
