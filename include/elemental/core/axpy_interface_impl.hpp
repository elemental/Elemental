/*
   Copyright (c) 2009-2013, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   All rights reserved.

   Authors:
   This interface is mainly due to Martin Schatz, but it was put into its
   current form by Jack Poulson.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_AXPYINTERFACE_IMPL_HPP
#define ELEM_CORE_AXPYINTERFACE_IMPL_HPP

namespace elem {

template<typename T>
inline bool
AxpyInterface<T>::Finished()
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::Finished");
    if( !attachedForLocalToGlobal_ && !attachedForGlobalToLocal_ )
        LogicError("Not attached");
#endif
    const Grid& g = ( attachedForLocalToGlobal_ ? 
                      localToGlobalMat_->Grid() : 
                      globalToLocalMat_->Grid() );
    const Int p = g.Size();

    bool finished = true; 
    for( Int rank=0; rank<p; ++rank )
    {
        if( !sentEomTo_[rank] || !haveEomFrom_[rank] )
        {
            finished = false;
            break;
        }
    }
    return finished;
}

template<typename T>
inline void
AxpyInterface<T>::HandleEoms()
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::HandleEoms");
#endif
    const Grid& g = ( attachedForLocalToGlobal_ ? 
                      localToGlobalMat_->Grid() : 
                      globalToLocalMat_->Grid() );
    const Int p = g.Size();

    UpdateRequestStatuses();

    // Try to progress our EOM sends
    for( Int i=0; i<p; ++i )
    {
        if( !sentEomTo_[i] )
        {
            bool shouldSendEom = true;
            const Int numSends = sendingData_[i].size();
            for( Int j=0; j<numSends; ++j )
            {
                if( sendingData_[i][j] )
                {
                    shouldSendEom = false;
                    break;
                }
            }
            const Int numRequests = sendingRequest_[i].size();
            for( Int j=0; j<numRequests; ++j )
            {
                if( !shouldSendEom || sendingRequest_[i][j] )
                {
                    shouldSendEom = false; 
                    break;
                }
            }
            const Int numReplies = sendingReply_[i].size();
            for( Int j=0; j<numReplies; ++j )
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
                mpi::TaggedISSend
                ( &sendDummy_, 1, i, EOM_TAG, g.VCComm(), request );
                sentEomTo_[i] = true;
            }
        }
    }

    mpi::Status status;
    if( mpi::IProbe( mpi::ANY_SOURCE, EOM_TAG, g.VCComm(), status ) )
    {
        const Int source = status.MPI_SOURCE;
        mpi::TaggedRecv( &recvDummy_, 1, source, EOM_TAG, g.VCComm() );
        haveEomFrom_[source] = true;
    }
}

template<typename T>
inline void
AxpyInterface<T>::HandleLocalToGlobalData()
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::HandleLocalToGlobalData");
#endif
    DistMatrix<T,MC,MR>& Y = *localToGlobalMat_;
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int myRow = g.Row();
    const Int myCol = g.Col();

    mpi::Status status;
    if( mpi::IProbe( mpi::ANY_SOURCE, DATA_TAG, g.VCComm(), status ) )
    {
        // Message exists, so recv and pack    
        const Int count = mpi::GetCount<byte>( status );
#ifndef RELEASE
        if( count < Int(4*sizeof(Int)+sizeof(T)) )
            LogicError("Count was too small");
#endif
        const Int source = status.MPI_SOURCE;
        recvVector_.resize( count );
        byte* recvBuffer = recvVector_.data();
        mpi::TaggedRecv( recvBuffer, count, source, DATA_TAG, g.VCComm() );

        // Extract the header
        byte* head = recvBuffer;
        const Int i = *reinterpret_cast<const Int*>(head); 
        head += sizeof(Int);
        const Int j = *reinterpret_cast<const Int*>(head); 
        head += sizeof(Int);
        const Int height = *reinterpret_cast<const Int*>(head); 
        head += sizeof(Int);
        const Int width = *reinterpret_cast<const Int*>(head); 
        head += sizeof(Int);
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
            RuntimeError( os.str() );
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
            RuntimeError( os.str() );
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
            RuntimeError( os.str() );
        }
#endif

        // Update Y
        const T* XBuffer = reinterpret_cast<const T*>(head);
        const Int colAlignment = (Y.ColAlignment()+i) % r;
        const Int rowAlignment = (Y.RowAlignment()+j) % c;
        const Int colShift = Shift( myRow, colAlignment, r );
        const Int rowShift = Shift( myCol, rowAlignment, c );

        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int iLocalOffset = Length( i, Y.ColShift(), r );
        const Int jLocalOffset = Length( j, Y.RowShift(), c );

        for( Int t=0; t<localWidth; ++t )
        {
            T* YCol = Y.Buffer(iLocalOffset,jLocalOffset+t);
            const T* XCol = &XBuffer[t*localHeight];
            for( Int s=0; s<localHeight; ++s )
                YCol[s] += alpha*XCol[s];
        }

        // Free the memory for the recv buffer
        recvVector_.clear();
    }
}
    
template<typename T>
inline void
AxpyInterface<T>::HandleGlobalToLocalRequest()
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::HandleGlobalToLocalRequest");
#endif
    const DistMatrix<T,MC,MR>& X = *globalToLocalMat_;
    const Grid& g = X.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int myRow = g.Row();
    const Int myCol = g.Col();

    mpi::Status status;
    if( mpi::IProbe( mpi::ANY_SOURCE, DATA_REQUEST_TAG, g.VCComm(), status ) )
    {
        // Request exists, so recv
        const Int source = status.MPI_SOURCE;
        const Int recvSize = 4*sizeof(Int);
        recvVector_.resize( recvSize );
        byte* recvBuffer = recvVector_.data();
        mpi::TaggedRecv
        ( recvBuffer, recvSize, source, DATA_REQUEST_TAG, g.VCComm() );

        // Extract the header
        const byte* recvHead = recvBuffer;
        const Int i = *reinterpret_cast<const Int*>(recvHead); 
        recvHead += sizeof(Int);
        const Int j = *reinterpret_cast<const Int*>(recvHead);
        recvHead += sizeof(Int);
        const Int height = *reinterpret_cast<const Int*>(recvHead);
        recvHead += sizeof(Int);
        const Int width = *reinterpret_cast<const Int*>(recvHead);
        recvHead += sizeof(Int);

        const Int colAlignment = (X.ColAlignment()+i) % r;
        const Int rowAlignment = (X.RowAlignment()+j) % c;
        const Int colShift = Shift( myRow, colAlignment, r );
        const Int rowShift = Shift( myCol, rowAlignment, c );

        const Int iLocalOffset = Length( i, X.ColShift(), r );
        const Int jLocalOffset = Length( j, X.RowShift(), c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight*localWidth;

        const Int bufferSize = 2*sizeof(Int) + numEntries*sizeof(T);
        const Int index = 
            ReadyForSend
            ( bufferSize, replyVectors_[source], 
              replySendRequests_[source], sendingReply_[source] );

        // Pack the reply header
        byte* sendBuffer = replyVectors_[source][index].data();
        byte* sendHead = sendBuffer;
        *reinterpret_cast<Int*>(sendHead) = myRow; sendHead += sizeof(Int);
        *reinterpret_cast<Int*>(sendHead) = myCol; sendHead += sizeof(Int);

        // Pack the payload
        T* sendData = reinterpret_cast<T*>(sendHead);
        for( Int t=0; t<localWidth; ++t )
        {
            T* sendCol = &sendData[t*localHeight];
            const T* XCol = X.LockedBuffer(iLocalOffset,jLocalOffset+t);
            MemCopy( sendCol, XCol, localHeight );
        }

        // Fire off non-blocking send
        mpi::TaggedISSend
        ( sendBuffer, bufferSize, source, DATA_REPLY_TAG, g.VCComm(), 
          replySendRequests_[source][index] );
    }
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
    CallStackEntry entry("AxpyInterface::AxpyInterface");
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

    const Int p = Z.Grid().Size();
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
}

template<typename T>
inline
AxpyInterface<T>::AxpyInterface
( AxpyType type, const DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::AxpyInterface");
#endif
    if( type == LOCAL_TO_GLOBAL )
    {
        LogicError("Cannot update a constant matrix");
    }
    else
    {
        attachedForLocalToGlobal_ = false;
        attachedForGlobalToLocal_ = true;
        localToGlobalMat_ = 0;
        globalToLocalMat_ = &X;
    }

    const Int p = X.Grid().Size();
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
           os << g.Rank()
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
    CallStackEntry entry("AxpyInterface::Attach");
#endif
    if( attachedForLocalToGlobal_ || attachedForGlobalToLocal_ )
        LogicError("Must detach before reattaching.");

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

    const Int p = Z.Grid().Size();
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
}

template<typename T>
inline void
AxpyInterface<T>::Attach( AxpyType type, const DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::Attach");
#endif
    if( attachedForLocalToGlobal_ || attachedForGlobalToLocal_ )
        LogicError("Must detach before reattaching.");

    if( type == LOCAL_TO_GLOBAL )
    {
        LogicError("Cannot update a constant matrix");
    }
    else
    {
        attachedForGlobalToLocal_ = true;
        globalToLocalMat_ = &X;
    }

    const Int p = X.Grid().Size();
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
}

template<typename T>
inline void 
AxpyInterface<T>::Axpy( T alpha, Matrix<T>& Z, Int i, Int j )
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::Axpy");
#endif
    if( attachedForLocalToGlobal_ )
        AxpyLocalToGlobal( alpha, Z, i, j );
    else if( attachedForGlobalToLocal_ )
        AxpyGlobalToLocal( alpha, Z, i, j );
    else
        LogicError("Cannot axpy before attaching.");
}

template<typename T>
inline void 
AxpyInterface<T>::Axpy( T alpha, const Matrix<T>& Z, Int i, Int j )
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::Axpy");
#endif
    if( attachedForLocalToGlobal_ )
        AxpyLocalToGlobal( alpha, Z, i, j );
    else if( attachedForGlobalToLocal_ )
        LogicError("Cannot update a constant matrix.");
    else
        LogicError("Cannot axpy before attaching.");
}

// Update Y(i:i+height-1,j:j+width-1) += alpha X, where X is height x width
template<typename T>
inline void
AxpyInterface<T>::AxpyLocalToGlobal
( T alpha, const Matrix<T>& X, Int i, Int j )
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::AxpyLocalToGlobal");
#endif
    DistMatrix<T,MC,MR>& Y = *localToGlobalMat_;
    if( i < 0 || j < 0 )
        LogicError("Submatrix offsets must be non-negative");
    if( i+X.Height() > Y.Height() || j+X.Width() > Y.Width() )
        LogicError("Submatrix out of bounds of global matrix");

    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlignment = (Y.ColAlignment() + i) % r;
    const Int rowAlignment = (Y.RowAlignment() + j) % c;

    const Int height = X.Height();
    const Int width = X.Width();

    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlignment, r );
        const Int rowShift = Shift( receivingCol, rowAlignment, c );
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight*localWidth;

        if( numEntries != 0 )
        {
            const Int destination = receivingRow + r*receivingCol;
            const Int bufferSize = 4*sizeof(Int) + (numEntries+1)*sizeof(T);

            const Int index = 
                ReadyForSend
                ( bufferSize, dataVectors_[destination], 
                  dataSendRequests_[destination], sendingData_[destination] );
#ifndef RELEASE
            if( Int(dataVectors_[destination][index].size()) != bufferSize )
                LogicError("Error in ReadyForSend");
#endif

            // Pack the header
            byte* sendBuffer = dataVectors_[destination][index].data();
            byte* head = sendBuffer;
            *reinterpret_cast<Int*>(head) = i; head += sizeof(Int);
            *reinterpret_cast<Int*>(head) = j; head += sizeof(Int);
            *reinterpret_cast<Int*>(head) = height; head += sizeof(Int);
            *reinterpret_cast<Int*>(head) = width; head += sizeof(Int);
            *reinterpret_cast<T*>(head) = alpha; head += sizeof(T);

            // Pack the payload
            T* sendData = reinterpret_cast<T*>(head);
            const T* XBuffer = X.LockedBuffer();
            const Int XLDim = X.LDim();
            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendData[t*localHeight];
                const T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

            // Fire off the non-blocking send
            mpi::TaggedISSend
            ( sendBuffer, bufferSize, destination, DATA_TAG, g.VCComm(), 
              dataSendRequests_[destination][index] );
        }

        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

// Update Y += alpha X(i:i+height-1,j:j+width-1), where X is the dist-matrix
template<typename T>
inline void
AxpyInterface<T>::AxpyGlobalToLocal
( T alpha, Matrix<T>& Y, Int i, Int j )
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::AxpyGlobalToLocal");
#endif
    const DistMatrix<T,MC,MR>& X = *globalToLocalMat_;

    const Int height = Y.Height();
    const Int width = Y.Width();
    if( i+height > X.Height() || j+width > X.Width() )
        LogicError("Invalid AxpyGlobalToLocal submatrix");

    const Grid& g = X.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();

    // Send out the requests to all processes in the grid
    for( Int rank=0; rank<p; ++rank )
    {
        const Int bufferSize = 4*sizeof(Int);
        const Int index = 
            ReadyForSend
            ( bufferSize, requestVectors_[rank], 
              requestSendRequests_[rank], sendingRequest_[rank] );

        // Copy the request header into the send buffer
        byte* sendBuffer = requestVectors_[rank][index].data();
        byte* head = sendBuffer;
        *reinterpret_cast<Int*>(head) = i; head += sizeof(Int);
        *reinterpret_cast<Int*>(head) = j; head += sizeof(Int);
        *reinterpret_cast<Int*>(head) = height; head += sizeof(Int);
        *reinterpret_cast<Int*>(head) = width; head += sizeof(Int);

        // Begin the non-blocking send
        mpi::TaggedISSend
        ( sendBuffer, bufferSize, rank, DATA_REQUEST_TAG, g.VCComm(), 
          requestSendRequests_[rank][index] );
    }

    // Receive all of the replies
    Int numReplies = 0;
    while( numReplies < p )
    {
        HandleGlobalToLocalRequest();

        mpi::Status status;
        if( mpi::IProbe( mpi::ANY_SOURCE, DATA_REPLY_TAG, g.VCComm(), status ) )
        {
            const Int source = status.MPI_SOURCE;

            // Ensure that we have a recv buffer
            const Int count = mpi::GetCount<byte>( status );
            recvVector_.resize( count );
            byte* recvBuffer = recvVector_.data();

            // Receive the data
            mpi::TaggedRecv
            ( recvBuffer, count, source, DATA_REPLY_TAG, g.VCComm() );

            // Unpack the reply header
            const byte* head = recvBuffer;
            const Int row = *reinterpret_cast<const Int*>(head); 
            head += sizeof(Int);
            const Int col = *reinterpret_cast<const Int*>(head); 
            head += sizeof(Int);
            const T* recvData = reinterpret_cast<const T*>(head);

            // Compute the local heights and offsets
            const Int colAlignment = (X.ColAlignment()+i) % r;
            const Int rowAlignment = (X.RowAlignment()+j) % c;
            const Int colShift = Shift( row, colAlignment, r );
            const Int rowShift = Shift( col, rowAlignment, c );
            const Int localHeight = Length( height, colShift, r );
            const Int localWidth = Length( width, rowShift, c );

            // Unpack the local matrix
            for( Int t=0; t<localWidth; ++t )
            {
                T* YCol = Y.Buffer(0,rowShift+t*c);
                const T* XCol = &recvData[t*localHeight];
                for( Int s=0; s<localHeight; ++s )
                    YCol[colShift+s*r] += alpha*XCol[s];
            }

            ++numReplies;
        }
    }
}

template<typename T>
inline Int
AxpyInterface<T>::ReadyForSend
( Int sendSize,
  std::deque<std::vector<byte>>& sendVectors,
  std::deque<mpi::Request>& requests, 
  std::deque<bool>& requestStatuses )
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::ReadyForSend");
#endif
    const Int numCreated = sendVectors.size();
#ifndef RELEASE
    if( numCreated != Int(requests.size()) || 
        numCreated != Int(requestStatuses.size()) )
        LogicError("size mismatch");
#endif
    for( Int i=0; i<numCreated; ++i )
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
            return i;
        }
    }

    sendVectors.resize( numCreated+1 );
    sendVectors[numCreated].resize( sendSize );
    requests.push_back( mpi::REQUEST_NULL );
    requestStatuses.push_back( true );

    return numCreated;
}

template<typename T>
inline void
AxpyInterface<T>::UpdateRequestStatuses()
{
#ifndef RELEASE
    CallStackEntry entry("AxpyInterface::UpdateRequestStatuses");
#endif
    const Grid& g = ( attachedForLocalToGlobal_ ? 
                      localToGlobalMat_->Grid() : 
                      globalToLocalMat_->Grid() );
    const Int p = g.Size();

    for( Int i=0; i<p; ++i )
    {
        const Int numDataSendRequests = dataSendRequests_[i].size();
        for( Int j=0; j<numDataSendRequests; ++j )
            if( sendingData_[i][j] )
                sendingData_[i][j] = 
                    !mpi::Test( dataSendRequests_[i][j] );
        const Int numRequestSendRequests = requestSendRequests_[i].size();
        for( Int j=0; j<numRequestSendRequests; ++j )
            if( sendingRequest_[i][j] )
                sendingRequest_[i][j] = 
                    !mpi::Test( requestSendRequests_[i][j] );
        const Int numReplySendRequests = replySendRequests_[i].size();
        for( Int j=0; j<numReplySendRequests; ++j )
            if( sendingReply_[i][j] )
                sendingReply_[i][j] = 
                    !mpi::Test( replySendRequests_[i][j] );
    }
}

template<typename T>
inline void
AxpyInterface<T>::Detach()
{
#ifndef RELEASE    
    CallStackEntry entry("AxpyInterface::Detach");
#endif
    if( !attachedForLocalToGlobal_ && !attachedForGlobalToLocal_ )
        LogicError("Must attach before detaching.");

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
}

} // namespace elem

#endif // ifndef ELEM_CORE_AXPYINTERFACE_IMPL_HPP
