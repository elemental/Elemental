/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_APPLYROWPIVOTS_HPP
#define ELEM_LAPACK_APPLYROWPIVOTS_HPP

#include "elemental/lapack-like/ComposePivots.hpp"

namespace elem {

template<typename F>
inline void
ApplyRowPivots( Matrix<F>& A, const Matrix<int>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyRowPivots");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
    if( p.Height() > A.Height() )
        throw std::logic_error("p cannot be larger than height of A");
#endif
    const int height = A.Height();
    const int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const int numPivots = p.Height();
    const int ldim = A.LDim();
    for( int i=0; i<numPivots; ++i )
    {
        const int k = p.Get(i,0);
        F* Ai = A.Buffer(i,0);
        F* Ak = A.Buffer(k,0);
        for( int j=0; j<width; ++j )
        {
            F temp = Ai[j*ldim];
            Ai[j*ldim] = Ak[j*ldim];
            Ak[j*ldim] = temp;
        }
    }
}

template<typename F>
inline void
ApplyInverseRowPivots( Matrix<F>& A, const Matrix<int>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyInverseRowPivots");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
    if( p.Height() > A.Height() )
        throw std::logic_error("p cannot be larger than height of A");
#endif
    const int height = A.Height();
    const int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const int numPivots = p.Height();
    const int ldim = A.LDim();
    for( int i=numPivots-1; i>=0; --i )
    {
        const int k = p.Get(i,0);
        F* Ai = A.Buffer(i,0);
        F* Ak = A.Buffer(k,0);
        for( int j=0; j<width; ++j )
        {
            F temp = Ai[j*ldim];
            Ai[j*ldim] = Ak[j*ldim];
            Ak[j*ldim] = temp;
        }
    }
}

template<typename F,Distribution U,Distribution V>
inline void
ApplyRowPivots( DistMatrix<F>& A, const DistMatrix<int,U,V>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyRowPivots");
#endif
    DistMatrix<int,STAR,STAR> p_STAR_STAR( p );
    ApplyRowPivots( A, p_STAR_STAR );
}

template<typename F,Distribution U,Distribution V>
inline void
ApplyInverseRowPivots
( DistMatrix<F>& A, const DistMatrix<int,U,V>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyInverseRowPivots");
#endif
    DistMatrix<int,STAR,STAR> p_STAR_STAR( p );
    ApplyInverseRowPivots( A, p_STAR_STAR );
}

template<typename F>
inline void
ApplyRowPivots( DistMatrix<F>& A, const DistMatrix<int,STAR,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyRowPivots");
#endif
    std::vector<int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyRowPivots( A, image, preimage );
}

template<typename F>
inline void
ApplyInverseRowPivots
( DistMatrix<F>& A, const DistMatrix<int,STAR,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyInverseRowPivots");
#endif
    std::vector<int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyRowPivots( A, preimage, image );
}

template<typename F> 
inline void
ApplyRowPivots
( Matrix<F>& A, 
  const std::vector<int>& image,
  const std::vector<int>& preimage )
{
    const int b = image.size();
#ifndef RELEASE
    CallStackEntry entry("ApplyRowPivots");
    if( A.Height() < b || b != (int)preimage.size() )
        throw std::logic_error
        ("image and preimage must be vectors of equal length that are not "
         "taller than A.");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( m == 0 || n == 0 )
        return;

    // TODO: Optimize this routine

    // Make a copy of the first b rows
    Matrix<F> ARowPanView;
    LockedView( ARowPanView, A, 0, 0, b, n );
    Matrix<F> ARowPanCopy( ARowPanView );

    // Make a copy of the preimage rows
    Matrix<F> APreimageCopy( b, n );
    for( int i=0; i<b; ++i ) 
    {
        const int iPre = preimage[i];
        if( iPre >= b )
            for( int j=0; j<n; ++j )
                APreimageCopy.Set(i,j,A.Get(iPre,j));
    }

    // Apply the permutations
    for( int i=0; i<b; ++i )
    {
        const int iPre = preimage[i];
        const int iPost = image[i];
        // Move row[i] into row[image[i]]
        for( int j=0; j<n; ++j )
            A.Set(iPost,j,ARowPanCopy.Get(i,j));
        if( iPre >= b )
        {
            // Move row[preimage[i]] into row[i]
            for( int j=0; j<n; ++j )
                A.Set(i,j,APreimageCopy.Get(i,j));
        }
    }
}

template<typename F> 
inline void
ApplyRowPivots
( DistMatrix<F>& A, 
  const std::vector<int>& image,
  const std::vector<int>& preimage )
{
    const int b = image.size();
#ifndef RELEASE
    CallStackEntry entry("ApplyRowPivots");
    if( A.Height() < b || b != (int)preimage.size() )
        throw std::logic_error
        ("image and preimage must be vectors of equal length that are not "
         "taller than A.");
#endif
    const int localWidth = A.LocalWidth();
    if( A.Height() == 0 || A.Width() == 0 )
        return;
    
    // Extract the relevant process grid information
    const Grid& g = A.Grid();
    const int r = g.Height();
    const int colAlignment = A.ColAlignment();
    const int colShift = A.ColShift();
    const int myRow = g.Row();

    // Extract the send and recv counts from the image and preimage.
    // This process's sends may be logically partitioned into two sets:
    //   (a) sends from rows [0,...,b-1]
    //   (b) sends from rows [b,...]
    // The latter is analyzed with preimage, the former deduced with image.
    std::vector<int> sendCounts(r,0), recvCounts(r,0);
    for( int i=colShift; i<b; i+=r )
    {
        const int sendRow = image[i];         
        const int sendTo = (colAlignment+sendRow) % r; 
        sendCounts[sendTo] += localWidth;

        const int recvRow = preimage[i];
        const int recvFrom = (colAlignment+recvRow) % r;
        recvCounts[recvFrom] += localWidth;
    }
    for( int i=0; i<b; ++i )
    {
        const int sendRow = image[i];
        if( sendRow >= b )
        {
            const int sendTo = (colAlignment+sendRow) % r;
            if( sendTo == myRow )
            {
                const int sendFrom = (colAlignment+i) % r;
                recvCounts[sendFrom] += localWidth;
            }
        }

        const int recvRow = preimage[i];
        if( recvRow >= b )
        {
            const int recvFrom = (colAlignment+recvRow) % r;
            if( recvFrom == myRow )
            {
                const int recvTo = (colAlignment+i) % r;
                sendCounts[recvTo] += localWidth;
            }
        }
    }

    // Construct the send and recv displacements from the counts
    std::vector<int> sendDispls(r), recvDispls(r);
    int totalSend=0, totalRecv=0;
    for( int i=0; i<r; ++i )
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
    const int ALDim = A.LDim();
    std::vector<F> sendData(std::max(1,totalSend));
    std::vector<int> offsets(r,0);
    const int localHeight = Length( b, colShift, r );
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int sendRow = image[colShift+iLoc*r];
        const int sendTo = (colAlignment+sendRow) % r;
        const int offset = sendDispls[sendTo]+offsets[sendTo];
        const F* ABuffer = A.Buffer(iLoc,0);
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
            sendData[offset+jLoc] = ABuffer[jLoc*ALDim];
        offsets[sendTo] += localWidth;
    }
    for( int i=0; i<b; ++i )
    {
        const int recvRow = preimage[i];
        if( recvRow >= b )
        {
            const int recvFrom = (colAlignment+recvRow) % r; 
            if( recvFrom == myRow )
            {
                const int recvTo = (colAlignment+i) % r;
                const int iLoc = (recvRow-colShift) / r;
                const int offset = sendDispls[recvTo]+offsets[recvTo];
                const F* ABuffer = A.Buffer(iLoc,0);
                for( int jLoc=0; jLoc<localWidth; ++jLoc )
                    sendData[offset+jLoc] = ABuffer[jLoc*ALDim];
                offsets[recvTo] += localWidth;
            }
        }
    }

    // Communicate all pivot rows
    std::vector<F> recvData(std::max(1,totalRecv));
    mpi::AllToAll
    ( &sendData[0], &sendCounts[0], &sendDispls[0],
      &recvData[0], &recvCounts[0], &recvDispls[0], g.ColComm() );

    // Unpack the recv data
    for( int k=0; k<r; ++k )
    {
        offsets[k] = 0;
        int thisColShift = Shift( k, colAlignment, r );
        for( int i=thisColShift; i<b; i+=r )
        {
            const int sendRow = image[i];
            const int sendTo = (colAlignment+sendRow) % r;
            if( sendTo == myRow )
            {
                const int offset = recvDispls[k]+offsets[k];
                const int iLoc = (sendRow-colShift) / r;
                F* ABuffer = A.Buffer(iLoc,0);
                for( int jLoc=0; jLoc<localWidth; ++jLoc )
                    ABuffer[jLoc*ALDim] = recvData[offset+jLoc];
                offsets[k] += localWidth;
            }
        }
    }
    for( int i=0; i<b; ++i )
    {
        const int recvRow = preimage[i];
        if( recvRow >= b )
        {
            const int recvTo = (colAlignment+i) % r;
            if( recvTo == myRow )
            {
                const int recvFrom = (colAlignment+recvRow) % r; 
                const int iLoc = (i-colShift) / r;
                const int offset = recvDispls[recvFrom]+offsets[recvFrom];
                F* ABuffer = A.Buffer(iLoc,0);
                for( int jLoc=0; jLoc<localWidth; ++jLoc )
                    ABuffer[jLoc*ALDim] = recvData[offset+jLoc];
                offsets[recvFrom] += localWidth;
            }
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_APPLYROWPIVOTS_HPP
