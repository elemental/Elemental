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

namespace elemental {

namespace axpy_type_wrapper {
enum AxpyType { LOCAL_TO_GLOBAL, GLOBAL_TO_LOCAL };
}
using namespace axpy_type_wrapper;

template<typename T>
class AxpyInterface
{   
public:
    AxpyInterface();
    ~AxpyInterface();

    AxpyInterface( AxpyType type,       DistMatrix<T,MC,MR>& Z );
    AxpyInterface( AxpyType type, const DistMatrix<T,MC,MR>& Z ); 

    void Attach( AxpyType type,       DistMatrix<T,MC,MR>& Z ); 
    void Attach( AxpyType type, const DistMatrix<T,MC,MR>& Z ); 

    void Axpy( T alpha,       Matrix<T>& Z, int i, int j );
    void Axpy( T alpha, const Matrix<T>& Z, int i, int j );

    void Detach();

private:
    static const int 
        DATA_TAG        =1, 
        EOM_TAG         =2, 
        DATA_REQUEST_TAG=3, 
        DATA_REPLY_TAG  =4;

    bool attachedForLocalToGlobal_, attachedForGlobalToLocal_;
    byte sendDummy_, recvDummy_;
    DistMatrix<T,MC,MR>* localToGlobalMat_;
    const DistMatrix<T,MC,MR>* globalToLocalMat_;

    std::vector<bool> sentEomTo_, haveEomFrom_;
    std::vector<byte> recvVector_;
    std::vector<mpi::Request> eomSendRequests_;

    std::vector<std::deque<std::vector<byte> > >
        dataVectors_, requestVectors_, replyVectors_;
    std::vector<std::deque<bool> > 
        sendingData_, sendingRequest_, sendingReply_;
    std::vector<std::deque<mpi::Request> > 
        dataSendRequests_, requestSendRequests_, replySendRequests_;

    // Check if we are done with this attachment's work
    bool Finished();

    // Progress functions
    void UpdateRequestStatuses();
    void HandleEoms();
    void HandleLocalToGlobalData();
    void HandleGlobalToLocalRequest();
    void StartSendingEoms();
    void FinishSendingEoms();

    void AxpyLocalToGlobal( T alpha, const Matrix<T>& X, int i, int j );
    void AxpyGlobalToLocal( T alpha,       Matrix<T>& Y, int i, int j );

    int ReadyForSend
    ( int sendSize,
      std::deque<std::vector<byte> >& sendVectors,
      std::deque<mpi::Request>& requests, 
      std::deque<bool>& requestStatuses );
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
inline bool
AxpyInterface<T>::Finished()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::Finished");
    if( !attachedForLocalToGlobal_ && !attachedForGlobalToLocal_ )
        throw std::logic_error("Not attached!");
#endif
    const Grid& g = ( attachedForLocalToGlobal_ ? 
                      localToGlobalMat_->Grid() : 
                      globalToLocalMat_->Grid() );
    const int p = g.Size();

    bool finished = true; 
    for( int rank=0; rank<p; ++rank )
    {
        if( !sentEomTo_[rank] || !haveEomFrom_[rank] )
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
inline void
AxpyInterface<T>::HandleEoms()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::HandleEoms");
#endif
    const Grid& g = ( attachedForLocalToGlobal_ ? 
                      localToGlobalMat_->Grid() : 
                      globalToLocalMat_->Grid() );
    const int p = g.Size();

    UpdateRequestStatuses();

    // Try to progress our EOM sends
    for( int i=0; i<p; ++i )
    {
        if( !sentEomTo_[i] )
        {
            bool shouldSendEom = true;
            for( int j=0; j<sendingData_[i].size(); ++j )
            {
                if( sendingData_[i][j] )
                {
                    shouldSendEom = false;
                    break;
                }
            }
            for( int j=0; j<sendingRequest_[i].size(); ++j )
            {
                if( !shouldSendEom || sendingRequest_[i][j] )
                {
                    shouldSendEom = false; 
                    break;
                }
            }
            for( int j=0; j<sendingReply_[i].size(); ++j )
            {
                if( !shouldSendEom || sendingReply_[i][j] )
                {
                    shouldSendEom = false;
                    break;
                }
            }
            if( shouldSendEom )
            {
                mpi::Request& request = eomSendRequests_[i];
                mpi::ISSend
                ( &sendDummy_, 1, i, EOM_TAG, g.VCComm(), request );
                sentEomTo_[i] = true;
            }
        }
    }

    mpi::Status status;
    if( mpi::IProbe( mpi::ANY_SOURCE, EOM_TAG, g.VCComm(), status ) )
    {
        const int source = status.MPI_SOURCE;
        mpi::Recv( &recvDummy_, 1, source, EOM_TAG, g.VCComm() );
        haveEomFrom_[source] = true;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
AxpyInterface<T>::HandleLocalToGlobalData()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::HandleLocalToGlobalData");
#endif
    DistMatrix<T,MC,MR>& Y = *localToGlobalMat_;
    const Grid& g = Y.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int myRow = g.MCRank();
    const int myCol = g.MRRank();

    mpi::Status status;
    if( mpi::IProbe( mpi::ANY_SOURCE, DATA_TAG, g.VCComm(), status ) )
    {
        // Message exists, so recv and pack    
        const int count = mpi::GetCount<byte>( status );
#ifndef RELEASE
        if( count < 4*sizeof(int)+sizeof(T) )
            throw std::logic_error("Count was too small");
#endif
        const int source = status.MPI_SOURCE;
        recvVector_.resize( count );
        byte* recvBuffer = &recvVector_[0];
        mpi::Recv( recvBuffer, count, source, DATA_TAG, g.VCComm() );

        // Extract the header
        byte* head = recvBuffer;
        const int i = *reinterpret_cast<const int*>(head); 
        head += sizeof(int);
        const int j = *reinterpret_cast<const int*>(head); 
        head += sizeof(int);
        const int height = *reinterpret_cast<const int*>(head); 
        head += sizeof(int);
        const int width = *reinterpret_cast<const int*>(head); 
        head += sizeof(int);
        const T alpha = *reinterpret_cast<const T*>(head); 
        head += sizeof(T);
#ifndef RELEASE
        if( height < 0 || width < 0 )
        {
            std::ostringstream os;
            os << "Unpacked heights were negative:\n"
               << "  i=     " << i << std::hex << "(" << i << ")\n" << std::dec
               << "  j=     " << j << std::hex << "(" << j << ")\n" << std::dec
               << "  height=" << height << std::hex << "(" << height << ")\n"
                              << std::dec
               << "  width= " << width << std::hex << "(" << width << ")\n"
                              << std::dec
               << "  alpha= " << alpha  << std::endl;
            throw std::runtime_error( os.str().c_str() );
        }
        if( i < 0 || j < 0 )
        {
            std::ostringstream os;
            os << "Unpacked offsets were negative:\n"
               << "  i=     " << i << std::hex << "(" << i << ")\n" << std::dec
               << "  j=     " << j << std::hex << "(" << j << ")\n" << std::dec
               << "  height=" << height << std::hex << "(" << height << ")\n"
                              << std::dec
               << "  width= " << width << std::hex << "(" << width << ")\n"
                              << std::dec
               << "  alpha= " << alpha  << std::endl;
            throw std::runtime_error( os.str().c_str() );
        }
        if( i+height > Y.Height() || j+width > Y.Width() )
        {
            std::ostringstream os;
            os << "Unpacked submatrix was out of bounds:\n"
               << "  i=     " << i << std::hex << "(" << i << ")\n" << std::dec
               << "  j=     " << j << std::hex << "(" << j << ")\n" << std::dec
               << "  height=" << height << std::hex << "(" << height << ")\n"
                              << std::dec
               << "  width= " << width << std::hex << "(" << width << ")\n"
                              << std::dec
               << "  alpha= " << alpha  << std::endl;
            throw std::runtime_error( os.str().c_str() );
        }
#endif

        // Update Y
        const T* XBuffer = reinterpret_cast<const T*>(head);
        const int colAlignment = (Y.ColAlignment()+i) % r;
        const int rowAlignment = (Y.RowAlignment()+j) % c;
        const int colShift = Shift( myRow, colAlignment, r );
        const int rowShift = Shift( myCol, rowAlignment, c );

        const int localHeight = LocalLength( height, colShift, r );
        const int localWidth = LocalLength( width, rowShift, c );
        const int iLocalOffset = LocalLength( i, Y.ColShift(), r );
        const int jLocalOffset = LocalLength( j, Y.RowShift(), c );

        for( int t=0; t<localWidth; ++t )
        {
            T* YCol = Y.LocalBuffer(iLocalOffset,jLocalOffset+t);
            const T* XCol = &XBuffer[t*localHeight];
            for( int s=0; s<localHeight; ++s )
                YCol[s] += alpha*XCol[s];
        }

        // Free the memory for the recv buffer
        recvVector_.clear();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
    
template<typename T>
inline void
AxpyInterface<T>::HandleGlobalToLocalRequest()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::HandleGlobalToLocalRequest");
#endif
    const DistMatrix<T,MC,MR>& X = *globalToLocalMat_;
    const Grid& g = X.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int myRow = g.MCRank();
    const int myCol = g.MRRank();

    mpi::Status status;
    if( mpi::IProbe( mpi::ANY_SOURCE, DATA_REQUEST_TAG, g.VCComm(), status ) )
    {
        // Request exists, so recv
        const int source = status.MPI_SOURCE;
        recvVector_.resize( 4*sizeof(int) );
        byte* recvBuffer = &recvVector_[0];
        mpi::Recv
        ( recvBuffer, 4*sizeof(int), source, DATA_REQUEST_TAG, g.VCComm() );

        // Extract the header
        const byte* recvHead = recvBuffer;
        const int i = *reinterpret_cast<const int*>(recvHead); 
        recvHead += sizeof(int);
        const int j = *reinterpret_cast<const int*>(recvHead);
        recvHead += sizeof(int);
        const int height = *reinterpret_cast<const int*>(recvHead);
        recvHead += sizeof(int);
        const int width = *reinterpret_cast<const int*>(recvHead);
        recvHead += sizeof(int);

        const int colAlignment = (X.ColAlignment()+i) % r;
        const int rowAlignment = (X.RowAlignment()+j) % c;
        const int colShift = Shift( myRow, colAlignment, r );
        const int rowShift = Shift( myCol, rowAlignment, c );

        const int iLocalOffset = LocalLength( i, X.ColShift(), r );
        const int jLocalOffset = LocalLength( j, X.RowShift(), c );
        const int localHeight = LocalLength( height, colShift, r );
        const int localWidth = LocalLength( width, rowShift, c );
        const int numEntries = localHeight*localWidth;

        const int bufferSize = 2*sizeof(int) + numEntries*sizeof(T);
        const int index = 
            ReadyForSend
            ( bufferSize, replyVectors_[source], 
              replySendRequests_[source], sendingReply_[source] );

        // Pack the reply header
        byte* sendBuffer = &replyVectors_[source][index][0];
        byte* sendHead = sendBuffer;
        *reinterpret_cast<int*>(sendHead) = myRow; sendHead += sizeof(int);
        *reinterpret_cast<int*>(sendHead) = myCol; sendHead += sizeof(int);

        // Pack the payload
        T* sendData = reinterpret_cast<T*>(sendHead);
        for( int t=0; t<localWidth; ++t )
        {
            T* sendCol = &sendData[t*localHeight];
            const T* XCol = X.LockedLocalBuffer(iLocalOffset,jLocalOffset+t);
            std::memcpy( sendCol, XCol, localHeight*sizeof(T) );
        }

        // Fire off non-blocking send
        mpi::ISSend
        ( sendBuffer, bufferSize, source, DATA_REPLY_TAG, g.VCComm(), 
          replySendRequests_[source][index] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
AxpyInterface<T>::AxpyInterface()
: attachedForLocalToGlobal_(false), attachedForGlobalToLocal_(false), 
  localToGlobalMat_(0), globalToLocalMat_(0)
{ }

template<typename T>
inline
AxpyInterface<T>::AxpyInterface( AxpyType type, DistMatrix<T,MC,MR>& Z )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::AxpyInterface");
#endif
    if( type == LOCAL_TO_GLOBAL )
    {
        attachedForLocalToGlobal_ = true;
        attachedForGlobalToLocal_ = false;
        localToGlobalMat_ = &Z;
        globalToLocalMat_ = 0;
    }
    else
    {
        attachedForLocalToGlobal_ = false;
        attachedForGlobalToLocal_ = true;
        localToGlobalMat_ = 0;
        globalToLocalMat_ = &Z;
    }

    const int p = Z.Grid().Size();
    sentEomTo_.resize( p, false );
    haveEomFrom_.resize( p, false );

    sendingData_.resize( p );
    sendingRequest_.resize( p );
    sendingReply_.resize( p );

    dataVectors_.resize( p );
    requestVectors_.resize( p );
    replyVectors_.resize( p );

    dataSendRequests_.resize( p );
    requestSendRequests_.resize( p );
    replySendRequests_.resize( p );

    eomSendRequests_.resize( p );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
AxpyInterface<T>::AxpyInterface
( AxpyType type, const DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::AxpyInterface");
#endif
    if( type == LOCAL_TO_GLOBAL )
    {
        throw std::logic_error("Cannot update a constant matrix");
    }
    else
    {
        attachedForLocalToGlobal_ = false;
        attachedForGlobalToLocal_ = true;
        localToGlobalMat_ = 0;
        globalToLocalMat_ = &X;
    }

    const int p = X.Grid().Size();
    sentEomTo_.resize( p, false );
    haveEomFrom_.resize( p, false );

    sendingData_.resize( p );
    sendingRequest_.resize( p );
    sendingReply_.resize( p );

    dataVectors_.resize( p );
    requestVectors_.resize( p );
    replyVectors_.resize( p );

    dataSendRequests_.resize( p );
    requestSendRequests_.resize( p );
    replySendRequests_.resize( p );

    eomSendRequests_.resize( p );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
AxpyInterface<T>::~AxpyInterface()
{ 
    if( attachedForLocalToGlobal_ || attachedForGlobalToLocal_ )
    {
        if( std::uncaught_exception() )
        {
           const Grid& g = ( attachedForLocalToGlobal_ ? 
                             localToGlobalMat_->Grid() : 
                             globalToLocalMat_->Grid() );
           std::ostringstream os;
           os << g.VCRank()
              << "Uncaught exception detected during AxpyInterface destructor "
                 "that required a call to Detach. Instead of allowing for the "
                 "possibility of Detach throwing another exception and "
                 "resulting in a 'terminate', we instead immediately dump the "
                 "call stack (if not in RELEASE mode) since the program will "
                 "likely hang:" << std::endl;
           std::cerr << os.str();
#ifndef RELEASE
           DumpCallStack();
#endif
        }
        else
        {
            Detach(); 
        }
    }
}

template<typename T>
inline void
AxpyInterface<T>::Attach( AxpyType type, DistMatrix<T,MC,MR>& Z )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::Attach");
#endif
    if( attachedForLocalToGlobal_ || attachedForGlobalToLocal_ )
        throw std::logic_error("Must detach before reattaching.");

    if( type == LOCAL_TO_GLOBAL )
    {
        attachedForLocalToGlobal_ = true;
        localToGlobalMat_ = &Z;
    }
    else
    {
        attachedForGlobalToLocal_ = true;
        globalToLocalMat_ = &Z;
    }

    const int p = Z.Grid().Size();
    sentEomTo_.resize( p, false );
    haveEomFrom_.resize( p, false );

    sendingData_.resize( p );
    sendingRequest_.resize( p );
    sendingReply_.resize( p );

    dataVectors_.resize( p );
    requestVectors_.resize( p );
    replyVectors_.resize( p );

    dataSendRequests_.resize( p );
    requestSendRequests_.resize( p );
    replySendRequests_.resize( p );

    eomSendRequests_.resize( p );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
AxpyInterface<T>::Attach( AxpyType type, const DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::Attach");
#endif
    if( attachedForLocalToGlobal_ || attachedForGlobalToLocal_ )
        throw std::logic_error("Must detach before reattaching.");

    if( type == LOCAL_TO_GLOBAL )
    {
        throw std::logic_error("Cannot update a constant matrix");
    }
    else
    {
        attachedForGlobalToLocal_ = true;
        globalToLocalMat_ = &X;
    }

    const int p = X.Grid().Size();
    sentEomTo_.resize( p, false );
    haveEomFrom_.resize( p, false );

    sendingData_.resize( p );
    sendingRequest_.resize( p );
    sendingReply_.resize( p );

    dataVectors_.resize( p );
    requestVectors_.resize( p );
    replyVectors_.resize( p );

    dataSendRequests_.resize( p );
    requestSendRequests_.resize( p );
    replySendRequests_.resize( p );

    eomSendRequests_.resize( p );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void 
AxpyInterface<T>::Axpy( T alpha, Matrix<T>& Z, int i, int j )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::Axpy");
#endif
    if( attachedForLocalToGlobal_ )
        AxpyLocalToGlobal( alpha, Z, i, j );
    else if( attachedForGlobalToLocal_ )
        AxpyGlobalToLocal( alpha, Z, i, j );
    else
        throw std::logic_error("Cannot axpy before attaching.");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void 
AxpyInterface<T>::Axpy( T alpha, const Matrix<T>& Z, int i, int j )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::Axpy");
#endif
    if( attachedForLocalToGlobal_ )
        AxpyLocalToGlobal( alpha, Z, i, j );
    else if( attachedForGlobalToLocal_ )
        throw std::logic_error("Cannot update a constant matrix.");
    else
        throw std::logic_error("Cannot axpy before attaching.");
#ifndef RELEASE
    PopCallStack();
#endif
}

// Update Y(i:i+height-1,j:j+width-1) += alpha X, where X is height x width
template<typename T>
inline void
AxpyInterface<T>::AxpyLocalToGlobal
( T alpha, const Matrix<T>& X, int i, int j )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::AxpyLocalToGlobal");
#endif
    DistMatrix<T,MC,MR>& Y = *localToGlobalMat_;
    if( i < 0 || j < 0 )
        throw std::logic_error("Submatrix offsets must be non-negative");
    if( i+X.Height() > Y.Height() || j+X.Width() > Y.Width() )
        throw std::logic_error("Submatrix out of bounds of global matrix");

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
            const int bufferSize = 4*sizeof(int) + (numEntries+1)*sizeof(T);

            const int index = 
                ReadyForSend
                ( bufferSize, dataVectors_[destination], 
                  dataSendRequests_[destination], sendingData_[destination] );
#ifndef RELEASE
            if( dataVectors_[destination][index].size() != bufferSize )
                throw std::logic_error("Error in ReadyForSend");
#endif

            // Pack the header
            byte* sendBuffer = &dataVectors_[destination][index][0];
            byte* head = sendBuffer;
            *reinterpret_cast<int*>(head) = i; head += sizeof(int);
            *reinterpret_cast<int*>(head) = j; head += sizeof(int);
            *reinterpret_cast<int*>(head) = height; head += sizeof(int);
            *reinterpret_cast<int*>(head) = width; head += sizeof(int);
            *reinterpret_cast<T*>(head) = alpha; head += sizeof(T);

            // Pack the payload
            T* sendData = reinterpret_cast<T*>(head);
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
            ( sendBuffer, bufferSize, destination, DATA_TAG, g.VCComm(), 
              dataSendRequests_[destination][index] );
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
inline void
AxpyInterface<T>::AxpyGlobalToLocal
( T alpha, Matrix<T>& Y, int i, int j )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::AxpyGlobalToLocal");
#endif
    const DistMatrix<T,MC,MR>& X = *globalToLocalMat_;

    const int height = Y.Height();
    const int width = Y.Width();
    if( i+height > X.Height() || j+width > X.Width() )
        throw std::logic_error("Invalid AxpyGlobalToLocal submatrix");

    const Grid& g = X.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();

    // Send out the requests to all processes in the grid
    for( int rank=0; rank<p; ++rank )
    {
        const int bufferSize = 4*sizeof(int);
        const int index = 
            ReadyForSend
            ( bufferSize, requestVectors_[rank], 
              requestSendRequests_[rank], sendingRequest_[rank] );

        // Copy the request header into the send buffer
        byte* sendBuffer = &requestVectors_[rank][index][0];
        byte* head = sendBuffer;
        *reinterpret_cast<int*>(head) = i; head += sizeof(int);
        *reinterpret_cast<int*>(head) = j; head += sizeof(int);
        *reinterpret_cast<int*>(head) = height; head += sizeof(int);
        *reinterpret_cast<int*>(head) = width; head += sizeof(int);

        // Begin the non-blocking send
        mpi::ISSend
        ( sendBuffer, bufferSize, rank, DATA_REQUEST_TAG, g.VCComm(), 
          requestSendRequests_[rank][index] );
    }

    // Receive all of the replies
    int numReplies = 0;
    while( numReplies < p )
    {
        HandleGlobalToLocalRequest();

        mpi::Status status;
        if( mpi::IProbe( mpi::ANY_SOURCE, DATA_REPLY_TAG, g.VCComm(), status ) )
        {
            const int source = status.MPI_SOURCE;

            // Ensure that we have a recv buffer
            const int count = mpi::GetCount<byte>( status );
            recvVector_.resize( count );
            byte* recvBuffer = &recvVector_[0];

            // Receive the data
            mpi::Recv( recvBuffer, count, source, DATA_REPLY_TAG, g.VCComm() );

            // Unpack the reply header
            const byte* head = recvBuffer;
            const int row = *reinterpret_cast<const int*>(head); 
            head += sizeof(int);
            const int col = *reinterpret_cast<const int*>(head); 
            head += sizeof(int);
            const T* recvData = reinterpret_cast<const T*>(head);

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
inline int
AxpyInterface<T>::ReadyForSend
( int sendSize,
  std::deque<std::vector<byte> >& sendVectors,
  std::deque<mpi::Request>& requests, 
  std::deque<bool>& requestStatuses )
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::ReadyForSend");
#endif
    const int commRank = mpi::CommRank( mpi::COMM_WORLD );

    const int numCreated = sendVectors.size();
#ifndef RELEASE
    if( numCreated != requests.size() || numCreated != requestStatuses.size() )
        throw std::logic_error("size mismatch");
#endif
    for( int i=0; i<numCreated; ++i )
    {
        // If this request is still running, test to see if it finished.
        if( requestStatuses[i] )
        {
            const bool finished = mpi::Test( requests[i] );
            requestStatuses[i] = !finished;
        }

        if( !requestStatuses[i] )    
        {
            requestStatuses[i] = true;
            sendVectors[i].resize( sendSize );
#ifndef RELEASE
            PopCallStack();
#endif
            return i;
        }
    }

    sendVectors.resize( numCreated+1 );
    sendVectors[numCreated].resize( sendSize );
    requests.push_back( mpi::REQUEST_NULL );
    requestStatuses.push_back( true );

#ifndef RELEASE
    PopCallStack();
#endif
    return numCreated;
}

template<typename T>
inline void
AxpyInterface<T>::UpdateRequestStatuses()
{
#ifndef RELEASE
    PushCallStack("AxpyInterface::UpdateRequestStatuses");
#endif
    const Grid& g = ( attachedForLocalToGlobal_ ? 
                      localToGlobalMat_->Grid() : 
                      globalToLocalMat_->Grid() );
    const int p = g.Size();

    for( int i=0; i<p; ++i )
    {
        for( int j=0; j<dataSendRequests_[i].size(); ++j )
            if( sendingData_[i][j] )
                sendingData_[i][j] = 
                    !mpi::Test( dataSendRequests_[i][j] );
        for( int j=0; j<requestSendRequests_[i].size(); ++j )
            if( sendingRequest_[i][j] )
                sendingRequest_[i][j] = 
                    !mpi::Test( requestSendRequests_[i][j] );
        for( int j=0; j<replySendRequests_[i].size(); ++j )
            if( sendingReply_[i][j] )
                sendingReply_[i][j] = 
                    !mpi::Test( replySendRequests_[i][j] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
AxpyInterface<T>::Detach()
{
#ifndef RELEASE    
    PushCallStack("AxpyInterface::Detach");
#endif
    if( !attachedForLocalToGlobal_ && !attachedForGlobalToLocal_ )
        throw std::logic_error("Must attach before detaching.");

    const Grid& g = ( attachedForLocalToGlobal_ ? 
                      localToGlobalMat_->Grid() : 
                      globalToLocalMat_->Grid() );

    while( !Finished() )
    {
        if( attachedForLocalToGlobal_ )
            HandleLocalToGlobalData();
        else
            HandleGlobalToLocalRequest();
        HandleEoms();
    }

    mpi::Barrier( g.VCComm() );

    attachedForLocalToGlobal_ = false;
    attachedForGlobalToLocal_ = false;
    recvVector_.clear();
    sentEomTo_.clear();
    haveEomFrom_.clear();

    sendingData_.clear();
    sendingRequest_.clear();
    sendingReply_.clear();

    dataVectors_.clear();
    requestVectors_.clear();
    replyVectors_.clear();
    
    dataSendRequests_.clear();
    requestSendRequests_.clear();
    replySendRequests_.clear();

    eomSendRequests_.clear();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // elemental

#endif  /* ELEMENTAL_AXPY_INTERFACE_HPP */

