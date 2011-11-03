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

template<typename F>
inline void
elemental::advanced::internal::ApplyRowPivots
( DistMatrix<F,MC,MR>& A, const DistMatrix<int,VC,STAR>& p, int pivotOffset )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::ApplyRowPivots");
#endif
    DistMatrix<int,STAR,STAR> p_STAR_STAR( p );
    advanced::internal::ApplyRowPivots( A, p_STAR_STAR, pivotOffset );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
elemental::advanced::internal::ApplyRowPivots
( DistMatrix<F,MC,MR>& A, const DistMatrix<int,STAR,STAR>& p, int pivotOffset )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::ApplyRowPivots");
#endif
    std::vector<int> image, preimage;
    if( pivotOffset == 0 && p.Height() == A.Height() )
    {
        advanced::internal::ComposePivots( p, image, preimage );
        advanced::internal::ApplyRowPivots( A, image, preimage );
    }
    else
    {
        advanced::internal::ComposePanelPivots
        ( p, image, preimage, pivotOffset );
        advanced::internal::ApplyRowPivots
        ( A, image, preimage, pivotOffset );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // represents a real or complex number
inline void
elemental::advanced::internal::ApplyRowPivots
( DistMatrix<F,MC,MR>& A, 
  const std::vector<int>& image,
  const std::vector<int>& preimage,
  int pivotOffset )
{
    const int b = image.size();
#ifndef RELEASE
    PushCallStack("advanced::internal::ApplyRowPivots");
    if( A.Height() < b || b != preimage.size() )
        throw std::logic_error
        ("image and preimage must be vectors of equal length that are not "
         "taller than A.");
    if( pivotOffset < 0 )
        throw std::logic_error("pivot offset must be non-negative");
#endif
    const int localWidth = A.LocalWidth();
    if( A.Width() == 0 )
        return;
    
    // Extract the relevant process grid information
    const Grid& g = A.Grid();
    const int r = g.Height();
    const int colAlignment = A.ColAlignment();
    const int colShift = A.ColShift();
    const int myRank = g.MCRank();

    // Extract the send and recv counts from the image and preimage.
    // This process's sends may be logically partitioned into two sets:
    //   (a) sends from rows [0,...,b-1]
    //   (b) sends from rows [b,...]
    // The latter is analyzed with image, the former deduced with preimage.
    std::vector<int> sendCounts(r,0), recvCounts(r,0);
    for( int i=colShift; i<b; i+=r )
    {
        const int sendRow = preimage[i];         
        const int sendTo = (colAlignment+sendRow) % r; 
        sendCounts[sendTo] += localWidth;

        const int recvRow = image[i];
        const int recvFrom = (colAlignment+recvRow) % r;
        recvCounts[recvFrom] += localWidth;
    }
    for( int i=0; i<b; ++i )
    {
        const int sendRow = preimage[i];
        if( sendRow >= b )
        {
            const int sendTo = (colAlignment+sendRow) % r;
            if( sendTo == myRank )
            {
                const int sendFrom = (colAlignment+i) % r;
                recvCounts[sendFrom] += localWidth;
            }
        }

        const int recvRow = image[i];
        if( recvRow >= b )
        {
            const int recvFrom = (colAlignment+recvRow) % r;
            if( recvFrom == myRank )
            {
                const int recvTo = (colAlignment+i) % r;
                sendCounts[recvTo] += localWidth;
            }
        }
    }

    // Construct the send and recv displacements from the counts
    std::vector<int> sendDispls(r), recvDispls(r);
    int totalSend = 0;
    for( int i=0; i<r; ++i )
    {
        sendDispls[i] = totalSend;
        totalSend += sendCounts[i];
    }
    int totalRecv = 0;
    for( int i=0; i<r; ++i )
    {
        recvDispls[i] = totalRecv;
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
    std::vector<F> sendData(max(1,totalSend));
    std::vector<int> offsets(r,0);
    const int localHeight = LocalLength( b, colShift, r );
    for( int i=0; i<localHeight; ++i )
    {
        const int sendRow = preimage[colShift+i*r];
        const int sendTo = (colAlignment+sendRow) % r;
        const int offset = sendDispls[sendTo]+offsets[sendTo];
        const int ALDim = A.LocalLDim();
        const F* ABuffer = A.LocalBuffer(i,0);
        for( int j=0; j<localWidth; ++j )     
            sendData[offset+j] = ABuffer[j*ALDim];
        offsets[sendTo] += localWidth;
    }
    for( int i=0; i<b; ++i )
    {
        const int recvRow = image[i];
        if( recvRow >= b )
        {
            const int recvFrom = (colAlignment+recvRow) % r; 
            if( recvFrom == myRank )
            {
                const int recvTo = (colAlignment+i) % r;
                const int iLocal = (recvRow-colShift) / r;
                const int offset = sendDispls[recvTo]+offsets[recvTo];
                const int ALDim = A.LocalLDim();
                const F* ABuffer = A.LocalBuffer(iLocal,0);
                for( int j=0; j<localWidth; ++j )
                    sendData[offset+j] = ABuffer[j*ALDim];
                offsets[recvTo] += localWidth;
            }
        }
    }

    // Communicate all pivot rows
    std::vector<F> recvData(std::max(1,totalRecv));
    mpi::AllToAll
    ( &sendData[0], &sendCounts[0], &sendDispls[0],
      &recvData[0], &recvCounts[0], &recvDispls[0], g.MCComm() );

    // Unpack the recv data
    for( int k=0; k<r; ++k )
    {
        offsets[k] = 0;
        int thisColShift = Shift( k, colAlignment, r );
        for( int i=thisColShift; i<b; i+=r )
        {
            const int sendRow = preimage[i];
            const int sendTo = (colAlignment+sendRow) % r;
            if( sendTo == myRank )
            {
                const int offset = recvDispls[k]+offsets[k];
                const int iLocal = (sendRow-colShift) / r;
                const int ALDim = A.LocalLDim();
                F* ABuffer = A.LocalBuffer(iLocal,0);
                for( int j=0; j<localWidth; ++j )
                    ABuffer[j*ALDim] = recvData[offset+j];
                offsets[k] += localWidth;
            }
        }
    }
    for( int i=0; i<b; ++i )
    {
        const int recvRow = image[i];
        if( recvRow >= b )
        {
            const int recvTo = (colAlignment+i) % r;
            if( recvTo == myRank )
            {
                const int recvFrom = (colAlignment+recvRow) % r; 
                const int iLocal = (i-colShift) / r;
                const int offset = recvDispls[recvFrom]+offsets[recvFrom];
                const int ALDim = A.LocalLDim();
                F* ABuffer = A.LocalBuffer(iLocal,0);
                for( int j=0; j<localWidth; ++j )
                    ABuffer[j*ALDim] = recvData[offset+j];
                offsets[recvFrom] += localWidth;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
