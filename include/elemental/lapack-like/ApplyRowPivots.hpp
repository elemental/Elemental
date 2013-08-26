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
ApplyRowPivots( Matrix<F>& A, const Matrix<Int>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyRowPivots");
    if( p.Width() != 1 )
        LogicError("p must be a column vector");
    if( p.Height() > A.Height() )
        LogicError("p cannot be larger than height of A");
#endif
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = p.Height();
    const Int ldim = A.LDim();
    for( Int i=0; i<numPivots; ++i )
    {
        const Int k = p.Get(i,0);
        F* Ai = A.Buffer(i,0);
        F* Ak = A.Buffer(k,0);
        for( Int j=0; j<width; ++j )
        {
            F temp = Ai[j*ldim];
            Ai[j*ldim] = Ak[j*ldim];
            Ak[j*ldim] = temp;
        }
    }
}

template<typename F>
inline void
ApplyInverseRowPivots( Matrix<F>& A, const Matrix<Int>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyInverseRowPivots");
    if( p.Width() != 1 )
        LogicError("p must be a column vector");
    if( p.Height() > A.Height() )
        LogicError("p cannot be larger than height of A");
#endif
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = p.Height();
    const Int ldim = A.LDim();
    for( Int i=numPivots-1; i>=0; --i )
    {
        const Int k = p.Get(i,0);
        F* Ai = A.Buffer(i,0);
        F* Ak = A.Buffer(k,0);
        for( Int j=0; j<width; ++j )
        {
            F temp = Ai[j*ldim];
            Ai[j*ldim] = Ak[j*ldim];
            Ak[j*ldim] = temp;
        }
    }
}

template<typename F,Distribution U,Distribution V>
inline void
ApplyRowPivots( DistMatrix<F>& A, const DistMatrix<Int,U,V>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyRowPivots");
#endif
    DistMatrix<Int,STAR,STAR> p_STAR_STAR( p );
    ApplyRowPivots( A, p_STAR_STAR );
}

template<typename F,Distribution U,Distribution V>
inline void
ApplyInverseRowPivots
( DistMatrix<F>& A, const DistMatrix<Int,U,V>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyInverseRowPivots");
#endif
    DistMatrix<Int,STAR,STAR> p_STAR_STAR( p );
    ApplyInverseRowPivots( A, p_STAR_STAR );
}

template<typename F>
inline void
ApplyRowPivots( DistMatrix<F>& A, const DistMatrix<Int,STAR,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyRowPivots");
#endif
    std::vector<Int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyRowPivots( A, image, preimage );
}

template<typename F>
inline void
ApplyInverseRowPivots
( DistMatrix<F>& A, const DistMatrix<Int,STAR,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyInverseRowPivots");
#endif
    std::vector<Int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyRowPivots( A, preimage, image );
}

template<typename F> 
inline void
ApplyRowPivots
( Matrix<F>& A, 
  const std::vector<Int>& image,
  const std::vector<Int>& preimage )
{
    const Int b = image.size();
#ifndef RELEASE
    CallStackEntry entry("ApplyRowPivots");
    if( A.Height() < b || b != (int)preimage.size() )
        LogicError
        ("image and preimage must be vectors of equal length that are not "
         "taller than A.");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( m == 0 || n == 0 )
        return;

    // TODO: Optimize this routine

    // Make a copy of the first b rows
    auto ARowPanView = LockedView( A, 0, 0, b, n );
    auto ARowPanCopy( ARowPanView );

    // Make a copy of the preimage rows
    Matrix<F> APreimageCopy( b, n );
    for( Int i=0; i<b; ++i ) 
    {
        const Int iPre = preimage[i];
        if( iPre >= b )
            for( Int j=0; j<n; ++j )
                APreimageCopy.Set(i,j,A.Get(iPre,j));
    }

    // Apply the permutations
    for( Int i=0; i<b; ++i )
    {
        const Int iPre = preimage[i];
        const Int iPost = image[i];
        // Move row[i] into row[image[i]]
        for( Int j=0; j<n; ++j )
            A.Set(iPost,j,ARowPanCopy.Get(i,j));
        if( iPre >= b )
        {
            // Move row[preimage[i]] into row[i]
            for( Int j=0; j<n; ++j )
                A.Set(i,j,APreimageCopy.Get(i,j));
        }
    }
}

template<typename F> 
inline void
ApplyRowPivots
( DistMatrix<F>& A, 
  const std::vector<Int>& image,
  const std::vector<Int>& preimage )
{
    const Int b = image.size();
#ifndef RELEASE
    CallStackEntry entry("ApplyRowPivots");
    if( A.Height() < b || b != (int)preimage.size() )
        LogicError
        ("image and preimage must be vectors of equal length that are not "
         "taller than A.");
#endif
    const Int localWidth = A.LocalWidth();
    if( A.Height() == 0 || A.Width() == 0 )
        return;
    
    // Extract the relevant process grid information
    const Grid& g = A.Grid();
    const Int r = g.Height();
    const Int colAlignment = A.ColAlignment();
    const Int colShift = A.ColShift();
    const Int myRow = g.Row();

    // Extract the send and recv counts from the image and preimage.
    // This process's sends may be logically partitioned into two sets:
    //   (a) sends from rows [0,...,b-1]
    //   (b) sends from rows [b,...]
    // The latter is analyzed with preimage, the former deduced with image.
    std::vector<int> sendCounts(r,0), recvCounts(r,0);
    for( Int i=colShift; i<b; i+=r )
    {
        const Int sendRow = image[i];         
        const Int sendTo = (colAlignment+sendRow) % r; 
        sendCounts[sendTo] += localWidth;

        const Int recvRow = preimage[i];
        const Int recvFrom = (colAlignment+recvRow) % r;
        recvCounts[recvFrom] += localWidth;
    }
    for( Int i=0; i<b; ++i )
    {
        const Int sendRow = image[i];
        if( sendRow >= b )
        {
            const Int sendTo = (colAlignment+sendRow) % r;
            if( sendTo == myRow )
            {
                const Int sendFrom = (colAlignment+i) % r;
                recvCounts[sendFrom] += localWidth;
            }
        }

        const Int recvRow = preimage[i];
        if( recvRow >= b )
        {
            const Int recvFrom = (colAlignment+recvRow) % r;
            if( recvFrom == myRow )
            {
                const Int recvTo = (colAlignment+i) % r;
                sendCounts[recvTo] += localWidth;
            }
        }
    }

    // Construct the send and recv displacements from the counts
    std::vector<int> sendDispls(r), recvDispls(r);
    Int totalSend=0, totalRecv=0;
    for( Int i=0; i<r; ++i )
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
        LogicError( msg.str() );
    }
#endif

    // Fill vectors with the send data
    const Int ALDim = A.LDim();
    std::vector<F> sendData( mpi::Pad(totalSend) );
    std::vector<int> offsets(r,0);
    const Int localHeight = Length( b, colShift, r );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int sendRow = image[colShift+iLoc*r];
        const Int sendTo = (colAlignment+sendRow) % r;
        const Int offset = sendDispls[sendTo]+offsets[sendTo];
        const F* ABuffer = A.Buffer(iLoc,0);
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            sendData[offset+jLoc] = ABuffer[jLoc*ALDim];
        offsets[sendTo] += localWidth;
    }
    for( Int i=0; i<b; ++i )
    {
        const Int recvRow = preimage[i];
        if( recvRow >= b )
        {
            const Int recvFrom = (colAlignment+recvRow) % r; 
            if( recvFrom == myRow )
            {
                const Int recvTo = (colAlignment+i) % r;
                const Int iLoc = (recvRow-colShift) / r;
                const Int offset = sendDispls[recvTo]+offsets[recvTo];
                const F* ABuffer = A.Buffer(iLoc,0);
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                    sendData[offset+jLoc] = ABuffer[jLoc*ALDim];
                offsets[recvTo] += localWidth;
            }
        }
    }

    // Communicate all pivot rows
    std::vector<F> recvData( mpi::Pad(totalRecv) );
    mpi::AllToAll
    ( &sendData[0], &sendCounts[0], &sendDispls[0],
      &recvData[0], &recvCounts[0], &recvDispls[0], g.ColComm() );

    // Unpack the recv data
    for( Int k=0; k<r; ++k )
    {
        offsets[k] = 0;
        Int thisColShift = Shift( k, colAlignment, r );
        for( Int i=thisColShift; i<b; i+=r )
        {
            const Int sendRow = image[i];
            const Int sendTo = (colAlignment+sendRow) % r;
            if( sendTo == myRow )
            {
                const Int offset = recvDispls[k]+offsets[k];
                const Int iLoc = (sendRow-colShift) / r;
                F* ABuffer = A.Buffer(iLoc,0);
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                    ABuffer[jLoc*ALDim] = recvData[offset+jLoc];
                offsets[k] += localWidth;
            }
        }
    }
    for( Int i=0; i<b; ++i )
    {
        const Int recvRow = preimage[i];
        if( recvRow >= b )
        {
            const Int recvTo = (colAlignment+i) % r;
            if( recvTo == myRow )
            {
                const Int recvFrom = (colAlignment+recvRow) % r; 
                const Int iLoc = (i-colShift) / r;
                const Int offset = recvDispls[recvFrom]+offsets[recvFrom];
                F* ABuffer = A.Buffer(iLoc,0);
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                    ABuffer[jLoc*ALDim] = recvData[offset+jLoc];
                offsets[recvFrom] += localWidth;
            }
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_APPLYROWPIVOTS_HPP
