/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "ElementalLAPACK_Internal.h"
using namespace std;
using namespace Elemental;
using namespace Elemental::utilities;
using namespace Elemental::wrappers::MPI;

template<typename T>
void
Elemental::LAPACK::Internal::ApplyRowPivots
(       DistMatrix<T,MC,MR>& A, 
  const vector<int>& image,
  const vector<int>& preimage,
  const int pivotOffset         )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::ApplyRowPivots");
    if( A.Height() < (int)image.size() || image.size() != preimage.size() )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "image and preimage must be vectors of equal length that"
                 << " are not taller than A." << endl;
        }
        DumpCallStack();
        throw exception();
    }
    if( pivotOffset < 0 )
    {
        if( A.GetGrid().VCRank() == 0 )
            cerr << "The pivot offset must be non-negative." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    if( A.Width() == 0 )
        return;

    const int b = image.size();
    const int localWidth = A.LocalWidth();
    
    // Extract the relevant process grid information
    const Grid& grid = A.GetGrid();
    const int r = grid.Height();
    const int colAlignment = A.ColAlignment();
    const int colShift = A.ColShift();
    const int myRank = grid.MCRank();

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
        cerr << "Send and recv counts do not match: (send,recv)=" 
             << totalSendCount << "," << totalRecvCount << endl;
        DumpCallStack();
        throw exception();
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
            sendData[offset+j] = A.LocalEntry(i,j);
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
                    sendData[offset+j] = A.LocalEntry(iLocal,j);
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

        AllToAll( sbuf, scs, sds, rbuf, rcs, rds, grid.MCComm() );
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
                    A.LocalEntry(iLocal,j) = recvData[offset+j];
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
                    A.LocalEntry(iLocal,j) = recvData[offset+j];
                offsets[recvFrom] += localWidth;
            }
        }
    }

#ifndef RELEASE
    PopCallStack();
#endif
}

template void
Elemental::LAPACK::Internal::ApplyRowPivots
(       DistMatrix<float,MC,MR>& A, 
  const vector<int>& image,
  const vector<int>& preimage,
  const int pivotOffset             );

template void
Elemental::LAPACK::Internal::ApplyRowPivots
(       DistMatrix<double,MC,MR>& A, 
  const vector<int>& image,
  const vector<int>& preimage,
  const int pivotOffset              );

#ifndef WITHOUT_COMPLEX
template void
Elemental::LAPACK::Internal::ApplyRowPivots
(       DistMatrix<scomplex,MC,MR>& A, 
  const vector<int>& image,
  const vector<int>& preimage,
  const int pivotOffset                );

template void
Elemental::LAPACK::Internal::ApplyRowPivots
(       DistMatrix<dcomplex,MC,MR>& A, 
  const vector<int>& image,
  const vector<int>& preimage,
  const int pivotOffset                );
#endif

