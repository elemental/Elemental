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

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

namespace elemental {

template<typename T>
class AxpyInterface
{   
    static const int DATA_TAG = 1;
    static const int ACK_TAG = 2;
    static const int EOM_TAG = 3;

    bool _attached;
    DistMatrix<T,MC,MR>* _Y;

    std::vector<bool> _canSendTo;
    std::vector<bool> _sentEomTo;
    std::vector<bool> _haveEomFrom;
    std::vector< std::vector<char> > _sendVectors;
    std::vector< std::vector<char> > _recvVectors;

    struct Header
    {
        int iStart, jStart;
        int localHeight, localWidth;
        T alpha;
    };

    // Check if we are done with this attachment's work
    bool Finished();

    // Progress functions
    void HandleAcks();
    void HandleEoms();
    void HandleData();
    void FinishSendingEoms();

public:
    AxpyInterface();
    AxpyInterface( DistMatrix<T,MC,MR>& Y );
    ~AxpyInterface();

    void Attach( DistMatrix<T,MC,MR>& Y );
    void Axpy( T alpha, const Matrix<T>& X, int i, int j );
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
#endif
    const int p = _Y->Grid().Size();

    bool finished = true; 
    for( int rank=0; rank<p; ++rank )
    {
        if( !_haveEomFrom[rank] || !_sentEomTo[rank] || !_canSendTo[rank] )
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
AxpyInterface<T>::HandleAcks()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::HandleAcks");
#endif
    const Grid& g = _Y->Grid();

    int test;
    imports::mpi::Status status;
    imports::mpi::IProbe
    ( imports::mpi::ANY_SOURCE, ACK_TAG, g.VCComm(), test, status );

    const int source = status.MPI_SOURCE;
    if( test )
    {
        int dummy;
        imports::mpi::Recv( &dummy, 1, source, ACK_TAG, g.VCComm() );
        _canSendTo[source] = true;
        _sendVectors[source].clear();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
AxpyInterface<T>::HandleEoms()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::HandleEoms");
#endif
    const Grid& g = _Y->Grid();

    int test;
    imports::mpi::Status status;
    imports::mpi::IProbe
    ( imports::mpi::ANY_SOURCE, EOM_TAG, g.VCComm(), test, status );

    const int source = status.MPI_SOURCE;
    if( test )
    {
        int dummy;
        imports::mpi::Recv( &dummy, 1, source, EOM_TAG, g.VCComm() );
        _haveEomFrom[source] = true;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
AxpyInterface<T>::HandleData()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::HandleData");
#endif
    const Grid& g = _Y->Grid();
    const int r = g.Height();
    const int c = g.Width();

    int test;
    imports::mpi::Status status;
    imports::mpi::IProbe
    ( imports::mpi::ANY_SOURCE, DATA_TAG, g.VCComm(), test, status );

    const int source = status.MPI_SOURCE;
    if( test )
    {
        // Message exists, so recv and pack    
        int count = imports::mpi::GetCount<char>( status );
        _recvVectors[source].resize( count );
        imports::mpi::Recv
        ( &_recvVectors[source][0], count, source, DATA_TAG, g.VCComm() );

        // Extract the header
        const char* recvBuffer = &_recvVectors[source][0];
        Header header;
        memcpy( &header, recvBuffer, sizeof(Header) );

        // Update Y
        const T* XBuffer = (T*)&recvBuffer[sizeof(Header)];
        const int colShift = _Y->ColShift();
        const int rowShift = _Y->RowShift();
        const int iLocalStart = (header.iStart-colShift) / r;
        const int jLocalStart = (header.jStart-rowShift) / c;
        T* YBuffer = _Y->LocalBuffer(iLocalStart,jLocalStart);
        const int YLDim = _Y->LocalLDim();
        for( int jLocal=0; jLocal<header.localWidth; ++jLocal )
        {
            T* thisYCol = &YBuffer[jLocal*YLDim];
            const T* thisXCol = &XBuffer[jLocal*header.localHeight];
            for( int iLocal=0; iLocal<header.localHeight; ++iLocal )
                thisYCol[iLocal] += header.alpha*thisXCol[iLocal];
        }

        // Free the memory for the recv buffer
        _recvVectors[source].clear();

        // Send an Ack to the source
        int dummy;
        imports::mpi::Request r;
        imports::mpi::ISSend( &dummy, 1, source, ACK_TAG, g.VCComm(), r );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
AxpyInterface<T>::FinishSendingEoms()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::FinishSendingEoms");
#endif
    const Grid& g = _Y->Grid();
    const int p = g.Size();
    const int myRank = g.VCRank();

    for( int i=0; i<p; ++i )
    {
        const int rank = (myRank+i) % p;
        if( !_sentEomTo[rank] )
        {
            int dummy;
            imports::mpi::Request r;
            imports::mpi::ISSend( &dummy, 1, rank, EOM_TAG, g.VCComm(), r );
            _sentEomTo[rank] = true;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
AxpyInterface<T>::AxpyInterface()
: _attached(false), _Y(0)
{ }

template<typename T>
AxpyInterface<T>::AxpyInterface( DistMatrix<T,MC,MR>& Y )
: _attached(true), _Y(&Y)
{
    _Y = &Y;

    const int p = Y.Grid().Size();

    _sendVectors.resize( p );
    _recvVectors.resize( p );
    _canSendTo.resize( p, true );
    _sentEomTo.resize( p, false );
    _haveEomFrom.resize( p, false );

    imports::mpi::Barrier( Y.Grid().VCComm() );
}

template<typename T>
AxpyInterface<T>::~AxpyInterface()
{ 
    if( _attached )
        Detach();
}

template<typename T>
void
AxpyInterface<T>::Attach( DistMatrix<T,MC,MR>& Y )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::Attach");
#endif
    if( _attached )
        throw std::logic_error("Must detach before reattaching.");

    _Y = &Y;

    const int p = Y.Grid().Size();

    _sendVectors.resize( p );
    _recvVectors.resize( p );
    _canSendTo.resize( p, true );
    _sentEomTo.resize( p, false );
    _haveEomFrom.resize( p, false );

    imports::mpi::Barrier( Y.Grid().VCComm() );
    _attached = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

// Update Y(i:i+height-1,j:j+width-1) += alpha X, where X is height x width
template<typename T>
void
AxpyInterface<T>::Axpy
( T alpha, const Matrix<T>& X, int i, int j )
{
#ifndef RELEASE
    PushCallStack("axpy_interface::Axpy");
#endif
    if( !_attached )
        throw std::logic_error("Must attach before axpying.");

    const Grid& g = _Y->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();
    const int myProcessRow = g.MCRank();
    const int myProcessCol = g.MRRank();
    const int colAlignment = (_Y->ColAlignment() + i) % r;
    const int rowAlignment = (_Y->RowAlignment() + j) % c;

    const int height = X.Height();
    const int width = X.Width();

    std::vector<imports::mpi::Request> requests(p);
    int receivingRow = myProcessRow;
    int receivingCol = myProcessCol;

    bool mustRetry;
    do
    {
        HandleAcks();
        HandleEoms();
        HandleData();

        const int colShift = utilities::Shift( receivingRow, colAlignment, r );
        const int rowShift = utilities::Shift( receivingCol, rowAlignment, c );

        const int localHeight = utilities::LocalLength( height, colShift, r );
        const int localWidth = utilities::LocalLength( width, rowShift, c );

        const int destination = receivingRow + r*receivingCol;

        mustRetry = false;
        if( _canSendTo[destination] )
        {
            const int numEntries = localHeight*localWidth;
            const int bufferSize = sizeof(Header) + numEntries*sizeof(T);
            if( numEntries > 0 )
            {
                // Make sure we have a big enough buffer
                _sendVectors[destination].resize( bufferSize );

                // Fill the header
                Header header;
                header.iStart = i + colShift;
                header.jStart = j + rowShift;
                header.localHeight = localHeight;
                header.localWidth = localWidth;
                header.alpha = alpha;

                // Pack the header
                char* sendBuffer = &_sendVectors[destination][0];
                memcpy( sendBuffer, &header, sizeof(Header) );

                // Pack the payload
                T* sendData = (T*)&sendBuffer[sizeof(Header)];
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
                imports::mpi::ISSend
                ( &_sendVectors[destination][0], bufferSize,
                  destination, DATA_TAG, g.VCComm(), requests[destination] );

                // Mark that we are already interacting with this process
                _canSendTo[destination] = false;
            }
        }
        else
        {
            while( !_canSendTo[destination] )
            {
                HandleAcks();
                HandleEoms();
                HandleData();
            }
            mustRetry = true;
        }

        if( !mustRetry )
        {
            receivingRow = (receivingRow + 1) % r;
            if( receivingRow == 0 )
                receivingCol = (receivingCol + 1) % c;
        }
    }
    while( mustRetry || 
           receivingRow != myProcessRow || receivingCol != myProcessCol );
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
    if( !_attached )
        throw std::logic_error("Must attach before detaching.");

    FinishSendingEoms();
    while( !Finished() )
    {
        HandleAcks();
        HandleEoms();
        HandleData();
    }

    _sendVectors.clear();
    _recvVectors.clear();
    _canSendTo.clear();
    _sentEomTo.clear();
    _haveEomFrom.clear();

    imports::mpi::Barrier( _Y->Grid().VCComm() );
    _attached = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

} // elemental

#endif  /* ELEMENTAL_AXPY_INTERFACE_HPP */

