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
elemental::advanced::ApplyRowPivots
( Matrix<F>& A, const Matrix<int>& p )
{
#ifndef RELEASE
    PushCallStack("advanced::ApplyRowPivots");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
    if( p.Height() != A.Height() )
        throw std::logic_error("p must be the same length as the height of A");
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

    const int ldim = A.LDim();
    for( int i=0; i<height; ++i )
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
elemental::advanced::ApplyInverseRowPivots
( Matrix<F>& A, const Matrix<int>& p )
{
#ifndef RELEASE
    PushCallStack("advanced::ApplyInverseRowPivots");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
    if( p.Height() != A.Height() )
        throw std::logic_error("p must be the same length as the height of A");
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

    const int ldim = A.LDim();
    for( int i=height-1; i>=0; --i )
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
elemental::advanced::ApplyRowPivots
( DistMatrix<F,MC,MR>& A, const DistMatrix<int,VC,STAR>& p )
{
#ifndef RELEASE
    PushCallStack("advanced::ApplyRowPivots");
#endif
    DistMatrix<int,STAR,STAR> p_STAR_STAR( p );
    advanced::ApplyRowPivots( A, p_STAR_STAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
elemental::advanced::ApplyInverseRowPivots
( DistMatrix<F,MC,MR>& A, const DistMatrix<int,VC,STAR>& p )
{
#ifndef RELEASE
    PushCallStack("advanced::ApplyInverseRowPivots");
#endif
    DistMatrix<int,STAR,STAR> p_STAR_STAR( p );
    advanced::ApplyInverseRowPivots( A, p_STAR_STAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
elemental::advanced::ApplyRowPivots
( DistMatrix<F,MC,MR>& A, const DistMatrix<int,STAR,STAR>& p )
{
#ifndef RELEASE
    PushCallStack("advanced::ApplyRowPivots");
#endif
    std::vector<int> image, preimage;
    advanced::ComposePivots( p, image, preimage );
    advanced::ApplyRowPivots( A, image, preimage );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
elemental::advanced::ApplyInverseRowPivots
( DistMatrix<F,MC,MR>& A, const DistMatrix<int,STAR,STAR>& p )
{
#ifndef RELEASE
    PushCallStack("advanced::ApplyInverseRowPivots");
#endif
    std::vector<int> image, preimage;
    advanced::ComposePivots( p, image, preimage );
    advanced::ApplyRowPivots( A, preimage, image );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // represents a real or complex number
inline void
elemental::advanced::ApplyRowPivots
( DistMatrix<F,MC,MR>& A, 
  const std::vector<int>& image,
  const std::vector<int>& preimage )
{
    const int b = image.size();
#ifndef RELEASE
    PushCallStack("advanced::ApplyRowPivots");
    if( A.Height() < b || b != preimage.size() )
        throw std::logic_error
        ("image and preimage must be vectors of equal length that are not "
         "taller than A.");
#endif
    const int localWidth = A.LocalWidth();
    if( A.Height() == 0 || A.Width() == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }
    
    // Extract the relevant process grid information
    const Grid& g = A.Grid();
    const int r = g.Height();
    const int colAlignment = A.ColAlignment();
    const int colShift = A.ColShift();
    const int myRow = g.MCRank();

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
            if( sendTo == myRow )
            {
                const int sendFrom = (colAlignment+i) % r;
                recvCounts[sendFrom] += localWidth;
            }
        }

        const int recvRow = image[i];
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
    const int ALDim = A.LocalLDim();
    std::vector<F> sendData(std::max(1,totalSend));
    std::vector<int> offsets(r,0);
    const int localHeight = LocalLength( b, colShift, r );
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int sendRow = preimage[colShift+iLocal*r];
        const int sendTo = (colAlignment+sendRow) % r;
        const int offset = sendDispls[sendTo]+offsets[sendTo];
        const F* ABuffer = A.LocalBuffer(iLocal,0);
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            sendData[offset+jLocal] = ABuffer[jLocal*ALDim];
        offsets[sendTo] += localWidth;
    }
    for( int i=0; i<b; ++i )
    {
        const int recvRow = image[i];
        if( recvRow >= b )
        {
            const int recvFrom = (colAlignment+recvRow) % r; 
            if( recvFrom == myRow )
            {
                const int recvTo = (colAlignment+i) % r;
                const int iLocal = (recvRow-colShift) / r;
                const int offset = sendDispls[recvTo]+offsets[recvTo];
                const F* ABuffer = A.LocalBuffer(iLocal,0);
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                    sendData[offset+jLocal] = ABuffer[jLocal*ALDim];
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
            if( sendTo == myRow )
            {
                const int offset = recvDispls[k]+offsets[k];
                const int iLocal = (sendRow-colShift) / r;
                F* ABuffer = A.LocalBuffer(iLocal,0);
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                    ABuffer[jLocal*ALDim] = recvData[offset+jLocal];
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
            if( recvTo == myRow )
            {
                const int recvFrom = (colAlignment+recvRow) % r; 
                const int iLocal = (i-colShift) / r;
                const int offset = recvDispls[recvFrom]+offsets[recvFrom];
                F* ABuffer = A.LocalBuffer(iLocal,0);
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                    ABuffer[jLocal*ALDim] = recvData[offset+jLocal];
                offsets[recvFrom] += localWidth;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
