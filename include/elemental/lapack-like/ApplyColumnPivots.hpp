/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_APPLYCOLUMNPIVOTS_HPP
#define ELEM_LAPACK_APPLYCOLUMNPIVOTS_HPP

#include "elemental/lapack-like/ComposePivots.hpp"

namespace elem {

template<typename F>
inline void
ApplyColumnPivots( Matrix<F>& A, const Matrix<Int>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyColumnPivots");
    if( p.Width() != 1 )
        LogicError("p must be a column vector");
    if( p.Height() > A.Width() )
        LogicError("p cannot be longer than width of A");
#endif
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = p.Height();
    for( Int j=0; j<numPivots; ++j )
    {
        const Int k = p.Get(j,0);
        F* Aj = A.Buffer(0,j);
        F* Ak = A.Buffer(0,k);
        for( Int i=0; i<height; ++i )
        {
            F temp = Aj[i];
            Aj[i] = Ak[i];
            Ak[i] = temp;
        }
    }
}

template<typename F>
inline void
ApplyInverseColumnPivots( Matrix<F>& A, const Matrix<Int>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyInverseColumnPivots");
    if( p.Width() != 1 )
        LogicError("p must be a column vector");
    if( p.Height() > A.Width() )
        LogicError("p cannot be larger than width of A");
#endif
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = p.Height();
    for( Int j=numPivots-1; j>=0; --j )
    {
        const Int k = p.Get(j,0);
        F* Aj = A.Buffer(0,j);
        F* Ak = A.Buffer(0,k);
        for( Int i=0; i<height; ++i )
        {
            F temp = Aj[i];
            Aj[i] = Ak[i];
            Ak[i] = temp;
        }
    }
}

template<typename F,Distribution U,Distribution V>
inline void
ApplyColumnPivots( DistMatrix<F>& A, const DistMatrix<Int,U,V>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyColumnPivots");
#endif
    DistMatrix<Int,STAR,STAR> p_STAR_STAR( p );
    ApplyColumnPivots( A, p_STAR_STAR );
}

template<typename F,Distribution U,Distribution V>
inline void
ApplyInverseColumnPivots
( DistMatrix<F>& A, const DistMatrix<Int,U,V>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyInverseColumnPivots");
#endif
    DistMatrix<Int,STAR,STAR> p_STAR_STAR( p );
    ApplyInverseColumnPivots( A, p_STAR_STAR );
}

template<typename F>
inline void
ApplyColumnPivots( DistMatrix<F>& A, const DistMatrix<Int,STAR,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyColumnPivots");
#endif
    std::vector<Int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyColumnPivots( A, image, preimage );
}

template<typename F>
inline void
ApplyInverseColumnPivots
( DistMatrix<F>& A, const DistMatrix<Int,STAR,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry entry("ApplyInverseColumnPivots");
#endif
    std::vector<Int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyColumnPivots( A, preimage, image );
}

template<typename F> 
inline void
ApplyColumnPivots
( Matrix<F>& A, 
  const std::vector<Int>& image,
  const std::vector<Int>& preimage )
{
    const Int b = image.size();
#ifndef RELEASE
    CallStackEntry entry("ApplyColumnPivots");
    if( A.Width() < b || b != int(preimage.size()) )
        LogicError
        ("image and preimage must be vectors of equal length that are not "
         "wider than A.");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( m == 0 || n == 0 )
        return;

    // TODO: Optimize this routine

    // Make a copy of the first b columns
    auto AColPanView = LockedView( A, 0, 0, m, b );
    auto AColPanCopy( AColPanView );

    // Make a copy of the preimage columns
    Matrix<F> APreimageCopy( m, b );
    for( Int j=0; j<b; ++j )
    {
        const Int jPre = preimage[j];
        if( jPre >= b )
            MemCopy( APreimageCopy.Buffer(0,j), A.LockedBuffer(0,jPre), m );
    }

    // Apply the permutations
    for( Int j=0; j<b; ++j )
    {
        const Int jPre = preimage[j];
        const Int jPost = image[j];
        // Move row[i] into row[image[i]]
        MemCopy( A.Buffer(0,jPost), AColPanCopy.LockedBuffer(0,j), m );
        // Move row[preimage[i]] into row[i]
        if( jPre >= b )
            MemCopy( A.Buffer(0,j), APreimageCopy.LockedBuffer(0,j), m );
    }
}

template<typename F> 
inline void
ApplyColumnPivots
( DistMatrix<F>& A, 
  const std::vector<Int>& image,
  const std::vector<Int>& preimage )
{
    const Int b = image.size();
#ifndef RELEASE
    CallStackEntry entry("ApplyColumnPivots");
    if( A.Width() < b || b != int(preimage.size()) )
        LogicError
        ("image and preimage must be vectors of equal length that are not "
         "wider than A.");
#endif
    const Int localHeight = A.LocalHeight();
    if( A.Height() == 0 || A.Width() == 0 )
        return;

    // Extract the relevant process grid information
    const Grid& g = A.Grid();
    const Int c = g.Width();
    const Int rowAlignment = A.RowAlignment();
    const Int rowShift = A.RowShift();
    const Int myCol = g.Col();

    // Extract the send and recv counts from the image and preimage.
    // This process's sends may be logically partitioned into two sets:
    //   (a) sends from rows [0,...,b-1]
    //   (b) sends from rows [b,...]
    // The latter is analyzed with preimage, the former deduced with image.
    std::vector<int> sendCounts(c,0), recvCounts(c,0);
    for( Int j=rowShift; j<b; j+=c )
    {
        const Int sendCol = image[j];         
        const Int sendTo = (rowAlignment+sendCol) % c; 
        sendCounts[sendTo] += localHeight;

        const Int recvCol = preimage[j];
        const Int recvFrom = (rowAlignment+recvCol) % c;
        recvCounts[recvFrom] += localHeight;
    }
    for( Int j=0; j<b; ++j )
    {
        const Int sendCol = image[j];
        if( sendCol >= b )
        {
            const Int sendTo = (rowAlignment+sendCol) % c;
            if( sendTo == myCol )
            {
                const Int sendFrom = (rowAlignment+j) % c;
                recvCounts[sendFrom] += localHeight;
            }
        }

        const Int recvCol = preimage[j];
        if( recvCol >= b )
        {
            const Int recvFrom = (rowAlignment+recvCol) % c;
            if( recvFrom == myCol )
            {
                const Int recvTo = (rowAlignment+j) % c;
                sendCounts[recvTo] += localHeight;
            }
        }
    }

    // Construct the send and recv displacements from the counts
    std::vector<int> sendDispls(c), recvDispls(c);
    Int totalSend=0, totalRecv=0;
    for( Int i=0; i<c; ++i )
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
    std::vector<F> sendData( mpi::Pad(totalSend) );
    std::vector<int> offsets(c,0);
    const Int localWidth = Length( b, rowShift, c );
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int sendCol = image[rowShift+jLoc*c];
        const Int sendTo = (rowAlignment+sendCol) % c;
        const Int offset = sendDispls[sendTo]+offsets[sendTo];
        MemCopy( &sendData[offset], A.Buffer(0,jLoc), localHeight );
        offsets[sendTo] += localHeight;
    }
    for( Int j=0; j<b; ++j )
    {
        const Int recvCol = preimage[j];
        if( recvCol >= b )
        {
            const Int recvFrom = (rowAlignment+recvCol) % c; 
            if( recvFrom == myCol )
            {
                const Int recvTo = (rowAlignment+j) % c;
                const Int jLoc = (recvCol-rowShift) / c;
                const Int offset = sendDispls[recvTo]+offsets[recvTo];
                MemCopy( &sendData[offset], A.Buffer(0,jLoc), localHeight );
                offsets[recvTo] += localHeight;
            }
        }
    }

    // Communicate all pivot rows
    std::vector<F> recvData( mpi::Pad(totalRecv) );
    mpi::AllToAll
    ( &sendData[0], &sendCounts[0], &sendDispls[0],
      &recvData[0], &recvCounts[0], &recvDispls[0], g.RowComm() );

    // Unpack the recv data
    for( Int k=0; k<c; ++k )
    {
        offsets[k] = 0;
        Int thisRowShift = Shift( k, rowAlignment, c );
        for( Int j=thisRowShift; j<b; j+=c )
        {
            const Int sendCol = image[j];
            const Int sendTo = (rowAlignment+sendCol) % c;
            if( sendTo == myCol )
            {
                const Int offset = recvDispls[k]+offsets[k];
                const Int jLoc = (sendCol-rowShift) / c;
                MemCopy( A.Buffer(0,jLoc), &recvData[offset], localHeight );
                offsets[k] += localHeight;
            }
        }
    }
    for( Int j=0; j<b; ++j )
    {
        const Int recvCol = preimage[j];
        if( recvCol >= b )
        {
            const Int recvTo = (rowAlignment+j) % c;
            if( recvTo == myCol )
            {
                const Int recvFrom = (rowAlignment+recvCol) % c; 
                const Int jLoc = (j-rowShift) / c;
                const Int offset = recvDispls[recvFrom]+offsets[recvFrom];
                MemCopy( A.Buffer(0,jLoc), &recvData[offset], localHeight );
                offsets[recvFrom] += localHeight;
            }
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_APPLYCOLUMNPIVOTS_HPP
