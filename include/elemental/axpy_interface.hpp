/*
   Copyright (c) 2009-2011, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
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

   Authors:
   This interface is mainly due to Martin Schatz, but it was put into its
   current form by Jack Poulson.
*/
#ifndef ELEMENTAL_AXPY_INTERFACE_HPP
#define ELEMENTAL_AXPY_INTERFACE_HPP 1

/* 
   Template conventions:
     G: general datatype
  
     T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
     Z: representation of a real ring, e.g., the integers or real numbers
     std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
                      or complex numbers

     F: representation of real or complex number
     R: representation of real number
     std::complex<R>: representation of complex number
*/

namespace elemental {

enum AxpyType { LOCAL_TO_GLOBAL, GLOBAL_TO_LOCAL };

template<typename T>
class AxpyInterface
{   
    static const int DATA_TAG=1, ACK_TAG=2, EOM_TAG=3, 
                     DATA_REQUEST_TAG=4, DATA_REPLY_TAG=5;

    bool _attachedForLocalToGlobal, _attachedForGlobalToLocal;
    char _sendDummy, _recvDummy;
    DistMatrix<T,MC,MR>* _Y;
    const DistMatrix<T,MC,MR>* _X;

    std::vector<bool> _canSendTo, _haveEomFrom;
    std::vector<std::vector<char> > _sendVectors, _recvVectors;

    struct LocalToGlobalHeader
    {
        int i, j, height, width;
        T alpha;
    };

    struct GlobalToLocalRequestHeader
    {
        int i, j, height, width;
    };

    struct GlobalToLocalReplyHeader
    {
        int processRow, processCol;
    };

    // Check if we are done with this attachment's work
    bool Finished();

    // Progress functions
    void HandleAck();
    void HandleEom();
    void HandleLocalToGlobalData();
    void HandleGlobalToLocalRequest();
    void SendEoms();

    void AxpyLocalToGlobal( T alpha, const Matrix<T>& X, int i, int j );
    void AxpyGlobalToLocal( T alpha,       Matrix<T>& Y, int i, int j );

public:
    AxpyInterface( AxpyType axpyType, DistMatrix<T,MC,MR>& Z );
    void Attach( AxpyType axpyType, DistMatrix<T,MC,MR>& Z ); 
    void Axpy( T alpha, Matrix<T>& Z, int i, int j );

    // Even though axpyType == LOCAL_TO_GLOBAL is illegal for these 
    // routines, this approach was chosen to keep the calling code consistent
    // and will throw an error with the illegal option.
    AxpyInterface( AxpyType axpyType, const DistMatrix<T,MC,MR>& Z ); 
    void Attach( AxpyType axpyType, const DistMatrix<T,MC,MR>& Z ); 
    void Axpy( T alpha, const Matrix<T>& Z, int i, int j );

    AxpyInterface();
    ~AxpyInterface();
    void Detach();
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
bool
AxpyInterface<T>::Finished()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::Finished");
    if( !_attachedForLocalToGlobal && !_attachedForGlobalToLocal )
        throw std::logic_error("Not attached!");
#endif
    const Grid& g = ( _attachedForLocalToGlobal ? _Y->Grid() : _X->Grid() );
    const int p = g.Size();

    bool finished = true; 
    for( int rank=0; rank<p; ++rank )
    {
        // For every message we send, we should get an ACK back. Thus, 
        // since we start our with _canSendTo set to 'true', it should
        // end up that way as well.
        if( !_canSendTo[rank] || !_haveEomFrom[rank] )
        {
            finished = false;
            break;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return finished;
}

template<typename T>
void
AxpyInterface<T>::HandleAck()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::HandleAck");
    if( !_attachedForLocalToGlobal && !_attachedForGlobalToLocal )
        throw std::logic_error("Not attached!");
#endif
    const Grid& g = ( _attachedForLocalToGlobal ? _Y->Grid() : _X->Grid() );

    int haveAck;
    imports::mpi::Status status;
    imports::mpi::IProbe
    ( imports::mpi::ANY_SOURCE, ACK_TAG, g.VCComm(), haveAck, status );

    if( haveAck )
    {
        const int source = status.MPI_SOURCE;
        imports::mpi::Recv( &_recvDummy, 1, source, ACK_TAG, g.VCComm() );
        _canSendTo[source] = true;
        _sendVectors[source].clear();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
AxpyInterface<T>::HandleEom()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::HandleEom");
#endif
    const Grid& g = ( _attachedForLocalToGlobal ? _Y->Grid() : _X->Grid() );

    int haveEom;
    imports::mpi::Status status;
    imports::mpi::IProbe
    ( imports::mpi::ANY_SOURCE, EOM_TAG, g.VCComm(), haveEom, status );

    if( haveEom )
    {
        const int source = status.MPI_SOURCE;
        imports::mpi::Recv( &_recvDummy, 1, source, EOM_TAG, g.VCComm() );
        _haveEomFrom[source] = true;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
AxpyInterface<T>::HandleLocalToGlobalData()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::HandleLocalToGlobalData");
#endif
    using namespace elemental::imports;
    using namespace elemental::utilities;

    DistMatrix<T,MC,MR>& Y = *_Y;
    const Grid& g = Y.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int myRow = g.MCRank();
    const int myCol = g.MRRank();

    int haveData;
    mpi::Status status;
    mpi::IProbe( mpi::ANY_SOURCE, DATA_TAG, g.VCComm(), haveData, status );

    if( haveData )
    {
        // Message exists, so recv and pack    
        const int count = mpi::GetCount<char>( status );
        const int source = status.MPI_SOURCE;
        _recvVectors[source].resize( count );
        mpi::Recv
        ( &_recvVectors[source][0], count, source, DATA_TAG, g.VCComm() );

        // Extract the header
        const char* recvBuffer = &_recvVectors[source][0];
        LocalToGlobalHeader header;
        std::memcpy( &header, recvBuffer, sizeof(LocalToGlobalHeader) );
        const int i = header.i;
        const int j = header.j;
        const int height = header.height;
        const int width = header.width;

        // Update Y
        const T* XBuffer = (const T*)&recvBuffer[sizeof(LocalToGlobalHeader)];
        const int colAlignment = (Y.ColAlignment()+i) % r;
        const int rowAlignment = (Y.RowAlignment()+j) % c;
        const int colShift = Shift( myRow, colAlignment, r );
        const int rowShift = Shift( myCol, rowAlignment, c );

        const int localHeight = LocalLength( height, colShift, r );
        const int localWidth = LocalLength( width, rowShift, c );
        const int iLocalOffset = LocalLength( header.i, Y.ColShift(), r );
        const int jLocalOffset = LocalLength( header.j, Y.RowShift(), c );

        const T alpha = header.alpha;
        for( int t=0; t<localWidth; ++t )
        {
            T* YCol = Y.LocalBuffer(iLocalOffset,jLocalOffset+t);
            const T* XCol = &XBuffer[t*localHeight];
            for( int s=0; s<localHeight; ++s )
                YCol[s] += alpha*XCol[s];
        }

        // Free the memory for the recv buffer
        _recvVectors[source].clear();

        // Send an ACK to the source
        mpi::Request r;
        mpi::ISSend( &_sendDummy, 1, source, ACK_TAG, g.VCComm(), r );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
    
template<typename T>
void
AxpyInterface<T>::HandleGlobalToLocalRequest()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::HandleGlobalToLocalRequest");
#endif
    using namespace elemental::imports;
    using namespace elemental::utilities;

    const DistMatrix<T,MC,MR>& X = *_X;
    const Grid& g = X.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int myRow = g.MCRank();
    const int myCol = g.MRRank();

    int haveRequest;
    mpi::Status status;
    mpi::IProbe
    ( mpi::ANY_SOURCE, DATA_REQUEST_TAG, g.VCComm(), haveRequest, status );

    if( haveRequest )
    {
        // Request exists, so recv
        const int source = status.MPI_SOURCE;
        _recvVectors[source].resize( sizeof(GlobalToLocalRequestHeader) );
        mpi::Recv
        ( &_recvVectors[source][0], sizeof(GlobalToLocalRequestHeader),
          source, DATA_REQUEST_TAG, g.VCComm() );

        // Send ACK to let source know messages can be sent
        mpi::Request request;
        mpi::ISSend( &_sendDummy, 1, source, ACK_TAG, g.VCComm(), request );

        // Extract the header
        const char* recvBuffer = &_recvVectors[source][0];
        GlobalToLocalRequestHeader requestHeader;
        std::memcpy
        ( &requestHeader, recvBuffer, sizeof(GlobalToLocalRequestHeader) );
        const int i = requestHeader.i;
        const int j = requestHeader.j;
        const int height = requestHeader.height;
        const int width = requestHeader.width;

        const int colAlignment = (X.ColAlignment()+i) % r;
        const int rowAlignment = (X.RowAlignment()+j) % c;
        const int colShift = Shift( myRow, colAlignment, r );
        const int rowShift = Shift( myCol, rowAlignment, c );

        const int iLocalOffset = LocalLength( i, X.ColShift(), r );
        const int jLocalOffset = LocalLength( j, X.RowShift(), c );
        const int localHeight = LocalLength( height, colShift, r );
        const int localWidth = LocalLength( width, rowShift, c );
        const int numEntries = localHeight*localWidth;

        GlobalToLocalReplyHeader replyHeader;
        replyHeader.processRow = myRow;
        replyHeader.processCol = myCol;
        const int bufferSize = 
            numEntries*sizeof(T) + sizeof(GlobalToLocalReplyHeader);

        // Make sure we have a big enough buffer
        _sendVectors[source].resize( bufferSize );
        char* sendBuffer = &_sendVectors[source][0];

        // Pack the reply header
        std::memcpy
        ( sendBuffer, &replyHeader, sizeof(GlobalToLocalReplyHeader) );

        // Pack the payload
        T* sendData = (T*)&sendBuffer[sizeof(GlobalToLocalReplyHeader)];
        for( int t=0; t<localWidth; ++t )
        {
            T* sendCol = &sendData[t*localHeight];
            const T* XCol = X.LockedLocalBuffer(iLocalOffset,jLocalOffset+t);
            std::memcpy( sendCol, XCol, localHeight*sizeof(T) );
        }

        // Fire off non-blocking send
        mpi::ISSend
        ( sendBuffer, bufferSize, source, DATA_REPLY_TAG, g.VCComm(), request );

        // 'source' is no longer waiting for a message
        _canSendTo[source] = false;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
AxpyInterface<T>::SendEoms()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::SendEoms");
#endif
    const Grid& g = ( _attachedForLocalToGlobal ? _Y->Grid() : _X->Grid() );
    const int p = g.Size();

    for( int rank=0; rank<p; ++rank )
    {
        imports::mpi::Request request;
        imports::mpi::ISSend
        ( &_sendDummy, 1, rank, EOM_TAG, g.VCComm(), request );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
AxpyInterface<T>::AxpyInterface()
: _attachedForLocalToGlobal(false), _attachedForGlobalToLocal(false), 
  _Y(0), _X(0)
{ }

template<typename T>
AxpyInterface<T>::AxpyInterface( AxpyType axpyType, DistMatrix<T,MC,MR>& Z )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::AxpyInterface");
#endif
    if( axpyType == LOCAL_TO_GLOBAL )
    {
        _attachedForLocalToGlobal = true;
        _attachedForGlobalToLocal = false;
        _Y = &Z;
        _X = 0;
    }
    else
    {
        _attachedForLocalToGlobal = false;
        _attachedForGlobalToLocal = true;
        _Y = 0;
        _X = &Z;
    }

    const int p = Z.Grid().Size();
    _sendVectors.resize( p );
    _recvVectors.resize( p );
    _canSendTo.resize( p, true );
    _haveEomFrom.resize( p, false );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
AxpyInterface<T>::AxpyInterface
( AxpyType axpyType, const DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::AxpyInterface");
#endif
    if( axpyType == LOCAL_TO_GLOBAL )
    {
        throw std::logic_error("Cannot update a constant matrix");
    }
    else
    {
        _attachedForLocalToGlobal = false;
        _attachedForGlobalToLocal = true;
        _Y = 0;
        _X = &X;
    }

    const int p = X.Grid().Size();
    _sendVectors.resize( p );
    _recvVectors.resize( p );
    _canSendTo.resize( p, true );
    _haveEomFrom.resize( p, false );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
AxpyInterface<T>::~AxpyInterface()
{ 
    if( _attachedForLocalToGlobal || _attachedForGlobalToLocal )
        Detach();
}

template<typename T>
void
AxpyInterface<T>::Attach( AxpyType axpyType, DistMatrix<T,MC,MR>& Z )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::Attach");
#endif
    if( _attachedForLocalToGlobal || _attachedForGlobalToLocal )
        throw std::logic_error("Must detach before reattaching.");

    if( axpyType == LOCAL_TO_GLOBAL )
    {
        _attachedForLocalToGlobal = true;
        _Y = &Z;
    }
    else
    {
        _attachedForGlobalToLocal = true;
        _X = &Z;
    }

    const int p = Z.Grid().Size();
    _sendVectors.resize( p );
    _recvVectors.resize( p );
    _canSendTo.resize( p, true );
    _haveEomFrom.resize( p, false );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
AxpyInterface<T>::Attach( AxpyType axpyType, const DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::Attach");
#endif
    if( _attachedForLocalToGlobal || _attachedForGlobalToLocal )
        throw std::logic_error("Must detach before reattaching.");

    if( axpyType == LOCAL_TO_GLOBAL )
    {
        throw std::logic_error("Cannot update a constant matrix");
    }
    else
    {
        _attachedForGlobalToLocal = true;
        _X = &X;
    }

    const int p = X.Grid().Size();
    _sendVectors.resize( p );
    _recvVectors.resize( p );
    _canSendTo.resize( p, true );
    _haveEomFrom.resize( p, false );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void 
AxpyInterface<T>::Axpy( T alpha, Matrix<T>& Z, int i, int j )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::Axpy");
#endif
    if( _attachedForLocalToGlobal )
        AxpyLocalToGlobal( alpha, Z, i, j );
    else if( _attachedForGlobalToLocal )
        AxpyGlobalToLocal( alpha, Z, i, j );
    else
        throw std::logic_error("Cannot axpy before attaching.");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void 
AxpyInterface<T>::Axpy( T alpha, const Matrix<T>& Z, int i, int j )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::Axpy");
#endif
    if( _attachedForLocalToGlobal )
        AxpyLocalToGlobal( alpha, Z, i, j );
    else if( _attachedForGlobalToLocal )
        throw std::logic_error("Cannot update a constant matrix.");
    else
        throw std::logic_error("Cannot axpy before attaching.");
#ifndef RELEASE
    PopCallStack();
#endif
}

// Update Y(i:i+height-1,j:j+width-1) += alpha X, where X is height x width
template<typename T>
void
AxpyInterface<T>::AxpyLocalToGlobal
( T alpha, const Matrix<T>& X, int i, int j )
{
#ifndef RELEASE
    PushCallStack("axpy_interface::AxpyLocalToGlobal");
#endif
    using namespace elemental::imports;
    using namespace elemental::utilities;

    DistMatrix<T,MC,MR>& Y = *_Y;

    const Grid& g = Y.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();
    const int myProcessRow = g.MCRank();
    const int myProcessCol = g.MRRank();
    const int colAlignment = (Y.ColAlignment() + i) % r;
    const int rowAlignment = (Y.RowAlignment() + j) % c;

    const int height = X.Height();
    const int width = X.Width();

    // Fill the header
    LocalToGlobalHeader header;
    header.i = i;
    header.j = j;
    header.height = height;
    header.width = width;
    header.alpha = alpha;

    std::vector<mpi::Request> requests(p);
    int receivingRow = myProcessRow;
    int receivingCol = myProcessCol;
    for( int step=0; step<p; ++step )
    {
        const int colShift = Shift( receivingRow, colAlignment, r );
        const int rowShift = Shift( receivingCol, rowAlignment, c );
        const int localHeight = LocalLength( height, colShift, r );
        const int localWidth = LocalLength( width, rowShift, c );
        const int numEntries = localHeight*localWidth;

        if( numEntries != 0 )
        {
            const int destination = receivingRow + r*receivingCol;
            while( !_canSendTo[destination] )
            {
                HandleAck();
                HandleEom();
                HandleLocalToGlobalData();
            }

            const int bufferSize = 
                sizeof(LocalToGlobalHeader) + numEntries*sizeof(T);

            // Make sure we have a big enough buffer
            _sendVectors[destination].resize( bufferSize );

            // Pack the header
            char* sendBuffer = &_sendVectors[destination][0];
            std::memcpy( sendBuffer, &header, sizeof(LocalToGlobalHeader) );

            // Pack the payload
            T* sendData = (T*)&sendBuffer[sizeof(LocalToGlobalHeader)];
            const T* XBuffer = X.LockedBuffer();
            const int XLDim = X.LDim();
            for( int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendData[t*localHeight];
                const T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
                for( int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

            // Fire off the non-blocking send
            mpi::ISSend
            ( sendBuffer, bufferSize, destination, 
              DATA_TAG, g.VCComm(), requests[destination] );

            // 'destination' is no longer ready for a message
            _canSendTo[destination] = false;
        }

        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Update Y += alpha X(i:i+height-1,j:j+width-1), where X is the dist-matrix
template<typename T>
void
AxpyInterface<T>::AxpyGlobalToLocal
( T alpha, Matrix<T>& Y, int i, int j )
{
#ifndef RELEASE
    PushCallStack("axpy_interface::AxpyGlobalToLocal");
#endif
    using namespace elemental::imports;
    using namespace elemental::utilities;

    const DistMatrix<T,MC,MR>& X = *_X;

    const int height = Y.Height();
    const int width = Y.Width();
    if( i+height > X.Height() || j+width > X.Width() )
        throw std::logic_error("Invalid AxpyGlobalToLocal submatrix");

    const Grid& g = X.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();
    std::vector<mpi::Request> requests(p);

    // Fill the request header
    GlobalToLocalRequestHeader requestHeader;
    requestHeader.i = i;
    requestHeader.j = j;
    requestHeader.height = height;
    requestHeader.width = width;

    // Send out the requests to all processes in the grid
    for( int rank=0; rank<p; ++rank )
    {
        while( !_canSendTo[rank] )
        {
            HandleAck();
            HandleEom();
            HandleGlobalToLocalRequest();
        }

        // Resize the vector and grab out the pointer
        _sendVectors[rank].resize(sizeof(GlobalToLocalRequestHeader));
        char* sendBuffer = &_sendVectors[rank][0];

        // Copy the request header into the send buffer
        std::memcpy
        ( sendBuffer, &requestHeader, sizeof(GlobalToLocalRequestHeader) );

        // Begin the non-blocking send
        mpi::ISSend
        ( sendBuffer, sizeof(GlobalToLocalRequestHeader),
          rank, DATA_REQUEST_TAG, g.VCComm(), requests[rank] );

        // That process is no longer ready for a message
        _canSendTo[rank] = false;
    }

    // Receive all of the replies
    int numReplies = 0;
    while( numReplies < p )
    {
        HandleAck();
        HandleEom();
        HandleGlobalToLocalRequest();

        int haveReply;
        mpi::Status status;
        mpi::IProbe
        ( mpi::ANY_SOURCE, DATA_REPLY_TAG, g.VCComm(), haveReply, status );

        if( haveReply )
        {
            const int source = status.MPI_SOURCE;

            // Ensure that we have a recv buffer
            const int count = mpi::GetCount<char>( status );
            _recvVectors[source].resize( count );
            char* recvBuffer = &_recvVectors[source][0];

            // Receive the data
            mpi::Recv( recvBuffer, count, source, DATA_REPLY_TAG, g.VCComm() );

            // Send an ACK to the source
            mpi::Request request;
            mpi::ISSend( &_sendDummy, 1, source, ACK_TAG, g.VCComm(), request );

            // Unpack the reply header
            GlobalToLocalReplyHeader replyHeader;
            std::memcpy
            ( &replyHeader, recvBuffer, sizeof(GlobalToLocalReplyHeader) );
            const int row = replyHeader.processRow;
            const int col = replyHeader.processCol;
            const T* recvData = 
                (const T*)&recvBuffer[sizeof(GlobalToLocalReplyHeader)];

            // Compute the local heights and offsets
            const int colAlignment = (X.ColAlignment()+i) % r;
            const int rowAlignment = (X.RowAlignment()+j) % c;
            const int colShift = Shift( row, colAlignment, r );
            const int rowShift = Shift( col, rowAlignment, c );
            const int localHeight = LocalLength( height, colShift, r );
            const int localWidth = LocalLength( width, rowShift, c );

            // Unpack the local matrix
            for( int t=0; t<localWidth; ++t )
            {
                T* YCol = Y.Buffer(0,rowShift+t*c);
                const T* XCol = &recvData[t*localHeight];
                for( int s=0; s<localHeight; ++s )
                    YCol[colShift+s*r] += alpha*XCol[s];
            }

            ++numReplies;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
AxpyInterface<T>::Detach()
{
#ifndef RELEASE    
    PushCallStack("AxpyInterface::Detach");
#endif
    if( !_attachedForLocalToGlobal && !_attachedForGlobalToLocal )
        throw std::logic_error("Must attach before detaching.");

    SendEoms();
    while( !Finished() )
    {
        HandleAck();
        HandleEom();
        if( _attachedForLocalToGlobal )
            HandleLocalToGlobalData();
        else
            HandleGlobalToLocalRequest();
    }

    const Grid& g = ( _attachedForLocalToGlobal ? _Y->Grid() : _X->Grid() );
    imports::mpi::Barrier( g.VCComm() );

    _attachedForLocalToGlobal = false;
    _attachedForGlobalToLocal = false;
    _sendVectors.clear();
    _recvVectors.clear();
    _canSendTo.clear();
    _haveEomFrom.clear();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // elemental

#endif  /* ELEMENTAL_AXPY_INTERFACE_HPP */

