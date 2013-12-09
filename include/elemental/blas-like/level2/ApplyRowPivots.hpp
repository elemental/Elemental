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

#include "elemental/blas-like/level2/ComposePivots.hpp"

namespace elem {

template<typename F>
inline void
ApplyRowPivots( Matrix<F>& A, const Matrix<Int>& p )
{
#ifndef RELEASE
    CallStackEntry cse("ApplyRowPivots");
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
    CallStackEntry cse("ApplyInverseRowPivots");
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

template<typename F,Distribution U1,Distribution V1,
                    Distribution U2,Distribution V2>
inline void
ApplyRowPivots( DistMatrix<F,U1,V1>& A, const DistMatrix<Int,U2,V2>& p )
{
#ifndef RELEASE
    CallStackEntry cse("ApplyRowPivots");
#endif
    DistMatrix<Int,STAR,STAR> p_STAR_STAR( p );
    ApplyRowPivots( A, p_STAR_STAR );
}

template<typename F,Distribution U1,Distribution V1,
                    Distribution U2,Distribution V2>
inline void
ApplyInverseRowPivots
( DistMatrix<F,U1,V1>& A, const DistMatrix<Int,U2,V2>& p )
{
#ifndef RELEASE
    CallStackEntry cse("ApplyInverseRowPivots");
#endif
    DistMatrix<Int,STAR,STAR> p_STAR_STAR( p );
    ApplyInverseRowPivots( A, p_STAR_STAR );
}

template<typename F,Distribution U,Distribution V>
inline void
ApplyRowPivots( DistMatrix<F,U,V>& A, const DistMatrix<Int,STAR,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry cse("ApplyRowPivots");
#endif
    std::vector<Int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyRowPivots( A, image, preimage );
}

template<typename F,Distribution U,Distribution V>
inline void
ApplyInverseRowPivots
( DistMatrix<F,U,V>& A, const DistMatrix<Int,STAR,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry cse("ApplyInverseRowPivots");
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
    CallStackEntry cse("ApplyRowPivots");
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

template<typename F,Distribution U,Distribution V> 
inline void
ApplyRowPivots
( DistMatrix<F,U,V>& A, 
  const std::vector<Int>& image,
  const std::vector<Int>& preimage )
{
    const Int b = image.size();
#ifndef RELEASE
    CallStackEntry cse("ApplyRowPivots");
    if( A.Height() < b || b != (int)preimage.size() )
        LogicError
        ("image and preimage must be vectors of equal length that are not "
         "taller than A.");
#endif
    const Int localWidth = A.LocalWidth();
    if( A.Height() == 0 || A.Width() == 0 || !A.Participating() )
        return;
    
    const Int colStride = A.ColStride();
    const Int colAlign = A.ColAlign();
    const Int colShift = A.ColShift();
    const Int colRank = A.ColRank();
    const mpi::Comm colComm = A.ColComm();

    // Extract the send and recv counts from the image and preimage.
    // This process's sends may be logically partitioned into two sets:
    //   (a) sends from rows [0,...,b-1]
    //   (b) sends from rows [b,...]
    // The latter is analyzed with preimage, the former deduced with image.
    std::vector<int> sendCounts(colStride,0), recvCounts(colStride,0);
    for( Int i=colShift; i<b; i+=colStride )
    {
        const Int sendRow = image[i];         
        const Int sendTo = (colAlign+sendRow) % colStride; 
        sendCounts[sendTo] += localWidth;

        const Int recvRow = preimage[i];
        const Int recvFrom = (colAlign+recvRow) % colStride;
        recvCounts[recvFrom] += localWidth;
    }
    for( Int i=0; i<b; ++i )
    {
        const Int sendRow = image[i];
        if( sendRow >= b )
        {
            const Int sendTo = (colAlign+sendRow) % colStride;
            if( sendTo == colRank )
            {
                const Int sendFrom = (colAlign+i) % colStride;
                recvCounts[sendFrom] += localWidth;
            }
        }

        const Int recvRow = preimage[i];
        if( recvRow >= b )
        {
            const Int recvFrom = (colAlign+recvRow) % colStride;
            if( recvFrom == colRank )
            {
                const Int recvTo = (colAlign+i) % colStride;
                sendCounts[recvTo] += localWidth;
            }
        }
    }

    // Construct the send and recv displacements from the counts
    std::vector<int> sendDispls(colStride), recvDispls(colStride);
    Int totalSend=0, totalRecv=0;
    for( Int i=0; i<colStride; ++i )
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
    std::vector<int> offsets(colStride,0);
    const Int localHeight = Length( b, colShift, colStride );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int sendRow = image[colShift+iLoc*colStride];
        const Int sendTo = (colAlign+sendRow) % colStride;
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
            const Int recvFrom = (colAlign+recvRow) % colStride; 
            if( recvFrom == colRank )
            {
                const Int recvTo = (colAlign+i) % colStride;
                const Int iLoc = (recvRow-colShift) / colStride;
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
    ( sendData.data(), sendCounts.data(), sendDispls.data(),
      recvData.data(), recvCounts.data(), recvDispls.data(), colComm );

    // Unpack the recv data
    for( Int k=0; k<colStride; ++k )
    {
        offsets[k] = 0;
        Int thisColShift = Shift( k, colAlign, colStride );
        for( Int i=thisColShift; i<b; i+=colStride )
        {
            const Int sendRow = image[i];
            const Int sendTo = (colAlign+sendRow) % colStride;
            if( sendTo == colRank )
            {
                const Int offset = recvDispls[k]+offsets[k];
                const Int iLoc = (sendRow-colShift) / colStride;
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
            const Int recvTo = (colAlign+i) % colStride;
            if( recvTo == colRank )
            {
                const Int recvFrom = (colAlign+recvRow) % colStride; 
                const Int iLoc = (i-colShift) / colStride;
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
