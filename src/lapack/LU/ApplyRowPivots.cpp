/*
   Copyright (c) 2009-2010, Jack Poulson
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
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::utilities;
using namespace elemental::wrappers::mpi;

template<typename T>
void
elemental::lapack::internal::ApplyRowPivots
(       DistMatrix<T,MC,MR>& A, 
  const vector<int>& image,
  const vector<int>& preimage,
  int pivotOffset )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::ApplyRowPivots");
    if( A.Height() < (int)image.size() || image.size() != preimage.size() )
        throw logic_error
        ( "image and preimage must be vectors of equal length that are not "
          "taller than A." );
    if( pivotOffset < 0 )
        throw logic_error( "pivot offset must be non-negative." );
#endif
    if( A.Width() == 0 )
        return;

    const int b = image.size();
    const int localWidth = A.LocalWidth();
    
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
    // The former is analyzed with image, the latter deduced with preimage.
    // Similar arguments hold for the recv counts.
    vector<int> sendCounts(r,0);
    vector<int> recvCounts(r,0);
    for( int i=colShift; i<b; i+=r )
    {
        const int sendLocation = image[i];         
        const int sendTo = ( colAlignment + sendLocation ) % r; 
        sendCounts[sendTo] += localWidth;

        const int recvLocation = preimage[i];
        const int recvFrom = ( colAlignment + recvLocation ) % r;
        recvCounts[recvFrom] += localWidth;
    }
    for( int i=0; i<b; ++i )
    {
        const int sendLocation = image[i];
        if( sendLocation >= b )
        {
            const int sendTo = ( colAlignment + sendLocation ) % r;
            if( sendTo == myRank )
            {
                const int sendFrom = ( colAlignment + i ) % r;
                recvCounts[sendFrom] += localWidth;
            }
        }

        const int recvLocation = preimage[i];
        if( recvLocation >= b )
        {
            const int recvFrom = ( colAlignment + recvLocation ) % r;
            if( recvFrom == myRank )
            {
                const int recvTo = ( colAlignment + i ) % r;
                sendCounts[recvTo] += localWidth;
            }
        }
    }

    // Construct the send and recv displacements from the counts
    vector<int> sendDispls(r);
    vector<int> recvDispls(r);
    int totalSendCount = 0;
    for( int i=0; i<r; ++i )
    {
        sendDispls[i] = totalSendCount;
        totalSendCount += sendCounts[i];
    }
    int totalRecvCount = 0;
    for( int i=0; i<r; ++i )
    {
        recvDispls[i] = totalRecvCount;
        totalRecvCount += recvCounts[i];
    }
#ifndef RELEASE
    if( totalSendCount != totalRecvCount )
    {
        ostringstream msg;
        msg << "Send and recv counts do not match: (send,recv)=" 
             << totalSendCount << "," << totalRecvCount << endl;
        throw logic_error( msg.str() );
    }
#endif

    // Fill vectors with the send data
    vector<T> sendData(max(1,totalSendCount));
    vector<int> offsets(r,0);
    const int localHeight = LocalLength( b, colShift, r );
    for( int i=0; i<localHeight; ++i )
    {
        const int sendLocation = image[colShift+i*r];
        const int sendTo = ( colAlignment + sendLocation ) % r;
        const int offset = sendDispls[sendTo]+offsets[sendTo];
        for( int j=0; j<localWidth; ++j )     
            sendData[offset+j] = A.GetLocalEntry(i,j);
        offsets[sendTo] += localWidth;
    }
    for( int i=0; i<b; ++i )
    {
        const int recvLocation = preimage[i];
        if( recvLocation >= b )
        {
            const int recvFrom = ( colAlignment + recvLocation ) % r; 
            if( recvFrom == myRank )
            {
                const int recvTo = ( colAlignment + i ) % r;
                const int iLocal = ( recvLocation - colShift ) / r;
                const int offset = sendDispls[recvTo]+offsets[recvTo];
                for( int j=0; j<localWidth; ++j )
                    sendData[offset+j] = A.GetLocalEntry(iLocal,j);
                offsets[recvTo] += localWidth;
            }
        }
    }

    // Communicate all pivot rows
    vector<T> recvData(max(1,totalRecvCount));
    {
        T* sbuf = &sendData[0];
        int* scs = &sendCounts[0];
        int* sds = &sendDispls[0];
        T* rbuf = &recvData[0];
        int* rcs = &recvCounts[0];
        int* rds = &recvDispls[0];

        AllToAll( sbuf, scs, sds, rbuf, rcs, rds, g.MCComm() );
    }

    // Unpack the recv data
    for( int k=0; k<r; ++k )
    {
        offsets[k] = 0;
        int thisColShift = Shift( k, colAlignment, r );
        for( int i=thisColShift; i<b; i+=r )
        {
            const int sendLocation = image[i];
            const int sendTo = ( colAlignment + sendLocation ) % r;
            if( sendTo == myRank )
            {
                const int offset = recvDispls[k]+offsets[k];
                const int iLocal = ( sendLocation - colShift ) / r;
                for( int j=0; j<localWidth; ++j )
                    A.SetLocalEntry(iLocal,j,recvData[offset+j]);
                offsets[k] += localWidth;
            }
        }
    }
    for( int i=0; i<b; ++i )
    {
        const int recvLocation = preimage[i];
        if( recvLocation >= b )
        {
            const int recvTo = ( colAlignment + i ) % r;
            if( recvTo == myRank )
            {
                const int recvFrom = ( colAlignment + recvLocation ) % r; 
                const int iLocal = ( i - colShift ) / r;
                const int offset = recvDispls[recvFrom]+offsets[recvFrom];
                for( int j=0; j<localWidth; ++j )
                    A.SetLocalEntry(iLocal,j,recvData[offset+j]);
                offsets[recvFrom] += localWidth;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::lapack::internal::ApplyRowPivots
(       DistMatrix<float,MC,MR>& A, 
  const vector<int>& image,
  const vector<int>& preimage,
        int pivotOffset );

template void
elemental::lapack::internal::ApplyRowPivots
(       DistMatrix<double,MC,MR>& A, 
  const vector<int>& image,
  const vector<int>& preimage,
        int pivotOffset );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::internal::ApplyRowPivots
(       DistMatrix<scomplex,MC,MR>& A, 
  const vector<int>& image,
  const vector<int>& preimage,
        int pivotOffset );

template void
elemental::lapack::internal::ApplyRowPivots
(       DistMatrix<dcomplex,MC,MR>& A, 
  const vector<int>& image,
  const vector<int>& preimage,
        int pivotOffset );
#endif

