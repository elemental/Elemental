/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {

template<typename T>
using ADM = AbstractDistMatrix<T>;

template<typename T,Dist U,Dist V>
using GDM = GeneralDistMatrix<T,U,V>;

// NOTE: It seems that member functions cannot be defined using a 
//       fully-specified template alias, e.g., ADM<T,U,V>::AbstractDistMatrix(),
//       but DM<T> is okay if it is only partially specified, e.g., 
//       DM<T> = DistMatrix<T,MC,MR> and DM<T>::DistMatrix()

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>::GeneralDistMatrix( GDM<T,U,V>&& A )
: ADM<T>(std::move(A))
{ }

// Assignment and reconfiguration
// ==============================

template<typename T,Dist U,Dist V>
GDM<T,U,V>& 
GeneralDistMatrix<T,U,V>::operator=( GDM<T,U,V>&& A )
{
    AbstractDistMatrix<T>::operator=( std::move(A) );
    return *this;
}

// Specialized redistribution/update routines
// ------------------------------------------

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::RowSumScatterFrom( const DistMatrix<T,U,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::RowSumScatterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    this->AlignColsAndResize( A.ColAlign(), A.Height(), A.Width() );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->RowSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::ColSumScatterFrom( const DistMatrix<T,UGath,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::ColSumScatterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    this->AlignRowsAndResize( A.RowAlign(), A.Height(), A.Width() );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->ColSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SumScatterFrom( const DistMatrix<T,UGath,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::SumScatterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    this->Resize( A.Height(), A.Width() );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->SumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::PartialRowSumScatterFrom
( const DistMatrix<T,U,VPart>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialRowSumScatterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    this->AlignWith( A );
    this->Resize( A.Height(), A.Width() );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->PartialRowSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::PartialColSumScatterFrom
( const DistMatrix<T,UPart,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialColSumScatterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    this->AlignWith( A );
    this->Resize( A.Height(), A.Width() );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->PartialColSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::RowSumScatterUpdate
( T alpha, const DistMatrix<T,U,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::RowSumScatterUpdate");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
    if( !this->Participating() )
        return;

    if( this->ColAlign() == A.ColAlign() )
    {
        if( this->Width() == 1 )
        {
            const Int rowAlign = this->RowAlign();
            const Int rowRank = this->RowRank();

            const Int localHeight = this->LocalHeight();
            const Int portionSize = mpi::Pad( localHeight );
            T* buffer = this->auxMemory_.Require( 2*portionSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[portionSize];

            // Pack 
            const T* ACol = A.LockedBuffer();
            MemCopy( sendBuf, ACol, localHeight );

            // Reduce to rowAlign
            mpi::Reduce
            ( sendBuf, recvBuf, portionSize, rowAlign, this->RowComm() );

            if( rowRank == rowAlign )
            {
                T* thisCol = this->Buffer();
                FMA_PARALLEL_FOR
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    thisCol[iLoc] += alpha*recvBuf[iLoc];
            }

            this->auxMemory_.Release();
        }
        else
        {
            const Int rowStride = this->RowStride();
            const Int rowAlign = this->RowAlign();

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int maxLocalWidth = MaxLength(width,rowStride);

            const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );
            const Int sendSize = rowStride*portionSize;

            // Pack 
            const Int ALDim = A.LDim();
            const T* ABuffer = A.LockedBuffer();
            T* buffer = this->auxMemory_.Require( sendSize );
            OUTER_PARALLEL_FOR
            for( Int k=0; k<rowStride; ++k )
            {
                T* data = &buffer[k*portionSize];
                const Int thisRowShift = Shift_( k, rowAlign, rowStride );
                const Int thisLocalWidth = 
                    Length_(width,thisRowShift,rowStride);
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* ACol = 
                        &ABuffer[(thisRowShift+jLoc*rowStride)*ALDim];
                    T* dataCol = &data[jLoc*localHeight];
                    MemCopy( dataCol, ACol, localHeight );
                }
            }
            // Communicate
            mpi::ReduceScatter( buffer, portionSize, this->RowComm() );

            // Update with our received data
            T* thisBuffer = this->Buffer();
            const Int thisLDim = this->LDim();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* bufferCol = &buffer[jLoc*localHeight];
                T* thisCol = &thisBuffer[jLoc*thisLDim];
                blas::Axpy( localHeight, alpha, bufferCol, 1, thisCol, 1 );
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned RowSumScatterUpdate" << std::endl;
#endif
        if( this->Width() == 1 )
        {
            const Int colStride = this->ColStride();
            const Int rowAlign = this->RowAlign();
            const Int colRank = this->ColRank();
            const Int rowRank = this->RowRank();

            const Int height = this->Height();
            const Int localHeight = this->LocalHeight();
            const Int localHeightA = A.LocalHeight();
            const Int maxLocalHeight = MaxLength(height,colStride);
            const Int portionSize = mpi::Pad( maxLocalHeight );

            const Int colAlign = this->ColAlign();
            const Int colAlignA = A.ColAlign();
            const Int sendRow = 
                (colRank+colStride+colAlign-colAlignA) % colStride;
            const Int recvRow = 
                (colRank+colStride+colAlignA-colAlign) % colStride;

            T* buffer = this->auxMemory_.Require( 2*portionSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[portionSize];

            // Pack 
            const T* ACol = A.LockedBuffer();
            MemCopy( sendBuf, ACol, localHeightA );

            // Reduce to rowAlign
            mpi::Reduce
            ( sendBuf, recvBuf, portionSize, rowAlign, this->RowComm() );

            if( rowRank == rowAlign )
            {
                // Perform the realignment
                mpi::SendRecv
                ( recvBuf, portionSize, sendRow,
                  sendBuf, portionSize, recvRow, this->ColComm() );

                T* thisCol = this->Buffer();
                FMA_PARALLEL_FOR
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    thisCol[iLoc] += alpha*sendBuf[iLoc];
            }
            this->auxMemory_.Release();
        }
        else
        {
            const Int colStride = this->ColStride();
            const Int rowStride = this->RowStride();
            const Int colRank = this->ColRank();

            const Int colAlign = this->ColAlign();
            const Int rowAlign = this->RowAlign();
            const Int colAlignA = A.ColAlign();
            const Int sendRow = 
                (colRank+colStride+colAlign-colAlignA) % colStride;
            const Int recvRow = 
                (colRank+colStride+colAlignA-colAlign) % colStride;

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localHeightA = A.LocalHeight();
            const Int maxLocalWidth = MaxLength(width,rowStride);

            const Int recvSize_RS = mpi::Pad( localHeightA*maxLocalWidth );
            const Int sendSize_RS = rowStride * recvSize_RS;
            const Int recvSize_SR = localHeight * localWidth;

            T* buffer = this->auxMemory_.Require
                ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );
            T* firstBuf = &buffer[0];
            T* secondBuf = &buffer[recvSize_RS];

            // Pack 
            const T* ABuffer = A.LockedBuffer();
            const Int ALDim = A.LDim();
            OUTER_PARALLEL_FOR
            for( Int k=0; k<rowStride; ++k )
            {
                T* data = &secondBuf[k*recvSize_RS];
                const Int thisRowShift = Shift_( k, rowAlign, rowStride );
                const Int thisLocalWidth = 
                    Length_(width,thisRowShift,rowStride);
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* ACol = 
                        &ABuffer[(thisRowShift+jLoc*rowStride)*ALDim];
                    T* dataCol = &data[jLoc*localHeightA];
                    MemCopy( dataCol, ACol, localHeightA );
                }
            }

            // Reduce-scatter over each process row
            mpi::ReduceScatter
            ( secondBuf, firstBuf, recvSize_RS, this->RowComm() );

            // Trade reduced data with the appropriate process row
            mpi::SendRecv
            ( firstBuf,  localHeightA*localWidth, sendRow,
              secondBuf, localHeight*localWidth,  recvRow, this->ColComm() );

            // Update with our received data
            T* thisBuffer = this->Buffer();
            const Int thisLDim = this->LDim();
            FMA_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* secondBufCol = &secondBuf[jLoc*localHeight];
                T* thisCol = &thisBuffer[jLoc*thisLDim];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    thisCol[iLoc] += alpha*secondBufCol[iLoc];
            }
            this->auxMemory_.Release();
        }
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::ColSumScatterUpdate
( T alpha, const DistMatrix<T,UGath,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::ColSumScatterUpdate");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
#ifdef VECTOR_WARNINGS
    if( A.Width() == 1 && this->Grid().Rank() == 0 )
    {
        std::cerr <<
          "The vector version of ColSumScatterUpdate does not"
          " yet have a vector version implemented, but it would only "
          "require a modification of the vector version of RowSumScatterUpdate"
          << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && this->Grid().Rank() == 0 )
    {
        std::cerr <<
          "ColSumScatterUpdate potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,V] matrix instead." << std::endl;
    }
#endif
    if( !this->Participating() )
        return;

    if( this->RowAlign() == A.RowAlign() )
    {
        const Int colStride = this->ColStride();
        const Int colAlign = this->ColAlign();
        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalHeight = MaxLength(height,colStride);

        const Int recvSize = mpi::Pad( maxLocalHeight*localWidth );
        const Int sendSize = colStride*recvSize;

        // Pack 
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
        T* buffer = this->auxMemory_.Require( sendSize );
        OUTER_PARALLEL_FOR
        for( Int k=0; k<colStride; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisColShift = Shift_( k, colAlign, colStride );
            const Int thisLocalHeight = Length_(height,thisColShift,colStride);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*colStride];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, this->ColComm() );

        // Update with our received data
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        FMA_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* bufferCol = &buffer[jLoc*localHeight];
            T* thisCol = &thisBuffer[jLoc*thisLDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                thisCol[iLoc] += alpha*bufferCol[iLoc];
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned ColSumScatterUpdate" << std::endl;
#endif
        const Int colStride = this->ColStride();
        const Int rowStride = this->RowStride();
        const Int rowRank = this->RowRank();

        const Int colAlign = this->ColAlign();
        const Int rowAlign = this->RowAlign();
        const Int rowAlignA = A.RowAlign();
        const Int sendCol = (rowRank+rowStride+rowAlign-rowAlignA) % rowStride;
        const Int recvCol = (rowRank+rowStride+rowAlignA-rowAlign) % rowStride;

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localWidthA = A.LocalWidth();
        const Int maxLocalHeight = MaxLength(height,colStride);

        const Int recvSize_RS = mpi::Pad( maxLocalHeight*localWidthA );
        const Int sendSize_RS = colStride * recvSize_RS;
        const Int recvSize_SR = localHeight * localWidth;

        T* buffer = this->auxMemory_.Require
            ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[recvSize_RS];

        // Pack
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<colStride; ++k )
        {
            T* data = &secondBuf[k*recvSize_RS];
            const Int thisColShift = Shift_( k, colAlign, colStride );
            const Int thisLocalHeight = Length_(height,thisColShift,colStride);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*colStride];
            }
        }

        // Reduce-scatter over each col
        mpi::ReduceScatter( secondBuf, firstBuf, recvSize_RS, this->ColComm() );

        // Trade reduced data with the appropriate col
        mpi::SendRecv
        ( firstBuf,  localHeight*localWidthA, sendCol,
          secondBuf, localHeight*localWidth,  recvCol, this->RowComm() );

        // Update with our received data
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        FMA_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* secondBufCol = &secondBuf[jLoc*localHeight];
            T* thisCol = &thisBuffer[jLoc*thisLDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                thisCol[iLoc] += alpha*secondBufCol[iLoc];
        }
        this->auxMemory_.Release();
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SumScatterUpdate
( T alpha, const DistMatrix<T,UGath,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::SumScatterUpdate");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int colAlign = this->ColAlign();
    const Int rowAlign = this->RowAlign();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int maxLocalHeight = MaxLength(height,colStride);
    const Int maxLocalWidth = MaxLength(width,rowStride);

    const Int recvSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
    const Int sendSize = colStride*rowStride*recvSize;

    // Pack 
    const T* ABuffer = A.LockedBuffer();
    const Int ALDim = A.LDim();
    T* buffer = this->auxMemory_.Require( sendSize );
    OUTER_PARALLEL_FOR
    for( Int l=0; l<rowStride; ++l )
    {
        const Int thisRowShift = Shift_( l, rowAlign, rowStride );
        const Int thisLocalWidth = Length_( width, thisRowShift, rowStride );
        for( Int k=0; k<colStride; ++k )
        {
            T* data = &buffer[(k+l*colStride)*recvSize];
            const Int thisColShift = Shift_( k, colAlign, colStride );
            const Int thisLocalHeight = Length_(height,thisColShift,colStride);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol =
                    &ABuffer[thisColShift+(thisRowShift+jLoc*rowStride)*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*colStride];
            }
        }
    }

    // Communicate
    mpi::ReduceScatter( buffer, recvSize, this->DistComm() );

    // Unpack our received data
    T* thisBuffer = this->Buffer();
    const Int thisLDim = this->LDim();
    FMA_PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* bufferCol = &buffer[jLoc*localHeight];
        T* thisCol = &thisBuffer[jLoc*thisLDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            thisCol[iLoc] += alpha*bufferCol[iLoc];
    }
    this->auxMemory_.Release();
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::PartialRowSumScatterUpdate
( T alpha, const DistMatrix<T,U,VPart>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialRowSumScatterUpdate");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
    if( !this->Participating() )
        return;

    if( this->RowAlign() % A.RowStride() == A.RowAlign() )
    {
        const Int rowStride = this->RowStride();
        const Int rowStridePart = this->PartialRowStride();
        const Int rowStrideUnion = this->PartialUnionRowStride();
        const Int rowRankPart = this->PartialRowRank();
        const Int rowAlign = this->RowAlign();
        const Int rowShiftOfA = A.RowShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalWidth = MaxLength( width, rowStride );
        const Int recvSize = mpi::Pad( height*maxLocalWidth );
        const Int sendSize = rowStrideUnion*recvSize;

        // Pack
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        T* buffer = this->auxMemory_.Require( sendSize );
        OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = rowRankPart+k*rowStridePart;
            const Int thisRowShift = Shift_( thisRank, rowAlign, rowStride );
            const Int thisRowOffset = 
                (thisRowShift-rowShiftOfA) / rowStridePart;
            const Int thisLocalWidth = 
                Length_( width, thisRowShift, rowStride );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* ACol = 
                    &ABuf[(thisRowOffset+jLoc*rowStrideUnion)*ALDim];
                T* dataCol = &data[jLoc*height];
                MemCopy( dataCol, ACol, height );
            }
        }
    
        // Communicate
        mpi::ReduceScatter( buffer, recvSize, this->PartialUnionRowComm() );

        // Unpack our received data
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* bufferCol = &buffer[jLoc*height];
            T* thisCol = &thisBuf[jLoc*thisLDim];
            for( Int i=0; i<height; ++i )
                thisCol[i] += alpha*bufferCol[i];
        }
        this->auxMemory_.Release();
    }
    else
    {
        LogicError("Unaligned PartialRowSumScatterUpdate not implemented");
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::PartialColSumScatterUpdate
( T alpha, const DistMatrix<T,UPart,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialColSumScatterUpdate");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
    if( !this->Participating() )
        return;

#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && A.Grid().Rank() == 0 )
    {
        std::cerr <<
          "PartialColSumScatterUpdate potentially causes a large amount"
          " of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [UGath,* ] matrix instead." 
          << std::endl;
    }
#endif
    if( this->ColAlign() % A.ColStride() == A.ColAlign() )
    {
        const Int colStride = this->ColStride();
        const Int colStridePart = this->PartialColStride();
        const Int colStrideUnion = this->PartialUnionColStride();
        const Int colRankPart = this->PartialColRank();
        const Int colAlign = this->ColAlign();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int maxLocalHeight = MaxLength( height, colStride );
        const Int recvSize = mpi::Pad( maxLocalHeight*width );
        const Int sendSize = colStrideUnion*recvSize;

        T* buffer = this->auxMemory_.Require( sendSize );

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = colRankPart+k*colStridePart;
            const Int thisColShift = Shift_( thisRank, colAlign, colStride );
            const Int thisColOffset = 
                (thisColShift-colShiftOfA) / colStridePart;
            const Int thisLocalHeight = 
                Length_( height, thisColShift, colStride );
            INNER_PARALLEL_FOR
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &data[j*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+j*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*colStrideUnion];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, this->PartialUnionColComm() );

        // Unpack our received data
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* bufferCol = &buffer[j*localHeight];
            T* thisCol = &thisBuf[j*thisLDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                thisCol[iLoc] += alpha*bufferCol[iLoc];
        }
        this->auxMemory_.Release();
    }
    else
    {
        LogicError("Unaligned PartialColSumScatterUpdate not implemented");
    }
}

// Diagonal manipulation
// =====================
template<typename T,Dist U,Dist V>
bool
GeneralDistMatrix<T,U,V>::DiagonalAlignedWith
( const elem::DistData& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::DiagonalAlignedWith"))
    if( this->Grid() != *d.grid )
        return false;

    const Int diagRoot = this->DiagonalRoot(offset);
    if( diagRoot != d.root )
        return false;

    const Int diagAlign = this->DiagonalAlign(offset);
    if( d.colDist == UDiag && d.rowDist == VDiag )
        return d.colAlign == diagAlign;
    else if( d.colDist == VDiag && d.rowDist == UDiag )
        return d.rowAlign == diagAlign;
    else
        return false;
}

template<typename T,Dist U,Dist V>
Int 
GeneralDistMatrix<T,U,V>::DiagonalRoot( Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::DiagonalRoot"))
    const elem::Grid& grid = this->Grid();

    if( U == MC && V == MR )
    {
        // Result is an [MD,* ] or [* ,MD]
        Int owner;
        if( offset >= 0 )
        {
            const Int procRow = this->ColAlign();
            const Int procCol = (this->RowAlign()+offset) % this->RowStride();
            owner = procRow + this->ColStride()*procCol;
        }
        else
        {
            const Int procRow = (this->ColAlign()-offset) % this->ColStride();
            const Int procCol = this->RowAlign();
            owner = procRow + this->ColStride()*procCol;
        }
        return grid.DiagPath(owner);
    }
    else if( U == MR && V == MC )
    {
        // Result is an [MD,* ] or [* ,MD]
        Int owner;
        if( offset >= 0 )
        {
            const Int procCol = this->ColAlign();
            const Int procRow = (this->RowAlign()+offset) % this->RowStride();
            owner = procRow + this->ColStride()*procCol;
        }
        else
        {
            const Int procCol = (this->ColAlign()-offset) % this->ColStride();
            const Int procRow = this->RowAlign();
            owner = procRow + this->ColStride()*procCol;
        }
        return grid.DiagPath(owner);
    }
    else
        return this->Root();
}

template<typename T,Dist U,Dist V>
Int
GeneralDistMatrix<T,U,V>::DiagonalAlign( Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::DiagonalAlign"))
    const elem::Grid& grid = this->Grid();

    if( U == MC && V == MR )
    {
        // Result is an [MD,* ] or [* ,MD]
        Int owner;
        if( offset >= 0 )
        {
            const Int procRow = this->ColAlign();
            const Int procCol = (this->RowAlign()+offset) % this->RowStride();
            owner = procRow + this->ColStride()*procCol;
        }
        else
        {
            const Int procRow = (this->ColAlign()-offset) % this->ColStride();
            const Int procCol = this->RowAlign();
            owner = procRow + this->ColStride()*procCol;
        }
        return grid.DiagPathRank(owner);
    }
    else if( U == MR && V == MC )
    {
        // Result is an [MD,* ] or [* ,MD]
        Int owner;
        if( offset >= 0 )
        {
            const Int procCol = this->ColAlign();
            const Int procRow = (this->RowAlign()+offset) % this->RowStride();
            owner = procRow + this->ColStride()*procCol;
        }
        else
        {
            const Int procCol = (this->ColAlign()-offset) % this->ColStride();
            const Int procRow = this->RowAlign();
            owner = procRow + this->ColStride()*procCol;
        }
        return grid.DiagPathRank(owner);
    }
    else if( U == STAR )
    {
        // Result is a [V,* ] or [* ,V]
        if( offset >= 0 )
            return (this->RowAlign()+offset) % this->RowStride();
        else
            return this->RowAlign();
    }
    else
    {
        // Result is [U,V] or [V,U], where V is either STAR or CIRC
        if( offset >= 0 )
            return this->ColAlign();
        else
            return (this->ColAlign()-offset) % this->ColStride();
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::GetDiagonal
( DistMatrix<T,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::GetDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::GetRealPartOfDiagonal
( DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::GetRealPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = RealPart(beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::GetImagPartOfDiagonal
( DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::GetImagPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = ImagPart(beta); } );
}

template<typename T,Dist U,Dist V>
auto
GeneralDistMatrix<T,U,V>::GetDiagonal( Int offset ) const
-> DistMatrix<T,UDiag,VDiag>
{
    DistMatrix<T,UDiag,VDiag> d( this->Grid() );
    GetDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
GeneralDistMatrix<T,U,V>::GetRealPartOfDiagonal( Int offset ) const
-> DistMatrix<Base<T>,UDiag,VDiag>
{
    DistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
GeneralDistMatrix<T,U,V>::GetImagPartOfDiagonal( Int offset ) const
-> DistMatrix<Base<T>,UDiag,VDiag>
{
    DistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SetDiagonal
( const DistMatrix<T,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::SetDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SetRealPartOfDiagonal
( const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::SetRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetRealPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SetImagPartOfDiagonal
( const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::SetImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetImagPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::UpdateDiagonal
( T gamma, const DistMatrix<T,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::UpdateDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, [gamma]( T& alpha, T beta ) { alpha += gamma*beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::UpdateRealPartOfDiagonal
( Base<T> gamma, const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::UpdateRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      [gamma]( T& alpha, Base<T> beta ) 
      { elem::UpdateRealPart(alpha,gamma*beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::UpdateImagPartOfDiagonal
( Base<T> gamma, const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::UpdateImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      [gamma]( T& alpha, Base<T> beta ) 
      { elem::UpdateImagPart(alpha,gamma*beta); } );
}

// Private section
// ###############

// Construct using a particular process grid
// =========================================

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>::GeneralDistMatrix( const elem::Grid& grid )
: ADM<T>(grid)
{ }

// Redistribution helper functions
// ===============================

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::ColAllGather( DistMatrix<T,UGath,V>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::ColAllGather");
        this->AssertSameGrid( A.Grid() );
    )
#ifdef CACHE_WARNINGS
    if( this->Height() != 1 && this->Grid().Rank() == 0 )
    {
        std::cerr <<
          "The matrix redistribution [* ,V] <- [U,V] potentially causes a "
          "large amount of cache-thrashing. If possible, avoid it by "
          "performing the redistribution with a (conjugate)-transpose: \n"
          << "  [V,* ].(Conjugate)TransposeFrom([U,V])" << std::endl;
    }
#endif
    A.AlignRowsAndResize( this->RowAlign(), this->Height(), this->Width() );

    if( this->Participating() )
    {
        if( this->RowAlign() == A.RowAlign() )
        {
            if( this->Height() == 1 )
            {
                const Int localWidthA = A.LocalWidth();
                T* bcastBuf = A.auxMemory_.Require( localWidthA );
            
                if( this->ColRank() == this->ColAlign() )
                {
                    A.matrix_ = this->LockedMatrix();
                    // Pack
                    const T* ABuf = A.LockedBuffer(); 
                    const Int ALDim = A.LDim();
                    PARALLEL_FOR
                    for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
                        bcastBuf[jLoc] = ABuf[jLoc*ALDim];
                }

                // Broadcast within the column comm
                mpi::Broadcast
                ( bcastBuf, localWidthA, this->ColAlign(), this->ColComm() );

                // Unpack
                T* ABuf = A.Buffer();
                const Int ALDim = A.LDim(); 
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
                    ABuf[jLoc*ALDim] = bcastBuf[jLoc];

                A.auxMemory_.Release();
            }
            else
            {
                const Int colStride = this->ColStride();
                const Int height = this->Height();
                const Int localWidth = this->LocalWidth();
                const Int thisLocalHeight = this->LocalHeight();
                const Int maxLocalHeight = MaxLength(height,colStride);
                const Int portionSize = mpi::Pad( maxLocalHeight*localWidth );

                T* buffer = A.auxMemory_.Require( (colStride+1)*portionSize );
                T* sendBuf = &buffer[0];
                T* recvBuf = &buffer[portionSize];

                // Pack
                const Int ldim = this->LDim();
                const T* thisBuf = this->LockedBuffer();
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                {
                    const T* thisCol = &thisBuf[jLoc*ldim];
                    T* sendBufCol = &sendBuf[jLoc*thisLocalHeight];
                    MemCopy( sendBufCol, thisCol, thisLocalHeight );
                }

                // Communicate
                mpi::AllGather
                ( sendBuf, portionSize, recvBuf, portionSize, this->ColComm() );

                // Unpack
                T* ABuf = A.Buffer();
                const Int ALDim = A.LDim();
                const Int colAlign = this->ColAlign();
                OUTER_PARALLEL_FOR
                for( Int k=0; k<colStride; ++k )
                {
                    const T* data = &recvBuf[k*portionSize];
                    const Int colShift = Shift_( k, colAlign, colStride );
                    const Int localHeight = 
                        Length_( height, colShift, colStride );
                    INNER_PARALLEL_FOR
                    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                    {
                        T* destCol = &ABuf[colShift+jLoc*ALDim];
                        const T* sourceCol = &data[jLoc*localHeight];
                        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                            destCol[iLoc*colStride] = sourceCol[iLoc];
                    }
                }
                A.auxMemory_.Release();
            }
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( this->Grid().Rank() == 0 )
                std::cerr << "Unaligned [U,V] -> [* ,V]." << std::endl;
#endif
            const Int colStride = this->ColStride();
            const Int rowStride = this->RowStride();
            const Int rowRank = this->RowRank();

            const Int rowAlign = this->RowAlign();
            const Int rowAlignA = A.RowAlign();
            const Int sendRowRank = 
                (rowRank+rowStride+rowAlignA-rowAlign) % rowStride;
            const Int recvRowRank = 
                (rowRank+rowStride+rowAlign-rowAlignA) % rowStride;

            if( this->Height() == 1 )
            {
                const Int localWidthA = A.LocalWidth();
                T* bcastBuf;

                if( this->ColRank() == this->ColAlign() )
                {
                    const Int localWidth = this->LocalWidth(); 
                
                    T* buffer = A.auxMemory_.Require( localWidth+localWidthA );
                    T* sendBuf = &buffer[0];
                    bcastBuf   = &buffer[localWidth];

                    // Pack
                    const T* thisBuf = this->LockedBuffer();
                    const Int ldim = this->LDim();
                    PARALLEL_FOR
                    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                        sendBuf[jLoc] = thisBuf[jLoc*ldim]; 

                    // Communicate
                    mpi::SendRecv
                    ( sendBuf,  localWidth,  sendRowRank,
                      bcastBuf, localWidthA, recvRowRank, this->RowComm() ); 
                }
                else
                {
                    bcastBuf = A.auxMemory_.Require( localWidthA );
                }

                // Communicate
                mpi::Broadcast
                ( bcastBuf, localWidthA, this->ColAlign(), this->ColComm() );

                // Unpack
                T* ABuf = A.Buffer();
                const Int ALDim = A.LDim();
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
                    ABuf[jLoc*ALDim] = bcastBuf[jLoc];

                A.auxMemory_.Release();
            }
            else
            {
                const Int height = this->Height();
                const Int localWidth = this->LocalWidth();
                const Int localWidthA = A.LocalWidth();
                const Int thisLocalHeight = this->LocalHeight();
                const Int maxLocalHeight = MaxLength(height,colStride);
                const Int portionSize = mpi::Pad( maxLocalHeight*localWidthA );

                T* buffer = A.auxMemory_.Require( (colStride+1)*portionSize );
                T* firstBuf  = &buffer[0];
                T* secondBuf = &buffer[portionSize];

                // Pack
                const Int ldim = this->LDim();
                const T* thisBuf = this->LockedBuffer();
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                {
                    const T* thisCol = &thisBuf[jLoc*ldim];
                    T* secondBufCol = &secondBuf[jLoc*thisLocalHeight];
                    MemCopy( secondBufCol, thisCol, thisLocalHeight ); 
                }
 
                // Realign
                mpi::SendRecv
                ( secondBuf, portionSize, sendRowRank,
                  firstBuf,  portionSize, recvRowRank, this->RowComm() );

                // AllGather the aligned data
                mpi::AllGather
                ( firstBuf, portionSize, secondBuf, portionSize, 
                  this->ColComm() ); 

                // Unpack the contents of each member of the column team
                T* ABuf = A.Buffer();
                const Int ALDim = A.LDim();
                const Int colAlign = this->ColAlign();
                OUTER_PARALLEL_FOR
                for( Int k=0; k<colStride; ++k )
                {
                    const T* data = &secondBuf[k*portionSize];
                    const Int colShift = Shift_( k, colAlign, colStride );
                    const Int localHeight = 
                        Length_( height, colShift, colStride );
                    INNER_PARALLEL_FOR
                    for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
                    {
                        T* destCol = &ABuf[colShift+jLoc*ALDim];
                        const T* sourceCol = &data[jLoc*localHeight];
                        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                            destCol[iLoc*colStride] = sourceCol[iLoc]; 
                    }
                }

                A.auxMemory_.Release();
            }
        }
    }
    if( this->Grid().InGrid() && this->CrossComm() != mpi::COMM_SELF )
    {
        // Pack from the root
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        T* buf = A.auxMemory_.Require( localHeight*localWidth );
        if( this->CrossRank() == this->Root() )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &buf[jLoc*localHeight], A.LockedBuffer(0,jLoc), localHeight );
        }

        // Broadcast from the root
        mpi::Broadcast
        ( buf, localHeight*localWidth, this->Root(), this->CrossComm() );

        // Unpack if not the root
        if( this->CrossRank() != this->Root() )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( A.Buffer(0,jLoc), &buf[jLoc*localHeight], localHeight );
        }
        A.auxMemory_.Release();
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::RowAllGather( DistMatrix<T,U,VGath>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::RowAllGather");
        this->AssertSameGrid( A.Grid() );
    )
    A.AlignColsAndResize( this->ColAlign(), this->Height(), this->Width() );

    if( this->Participating() )
    {
        if( this->ColAlign() == A.ColAlign() )
        {
            if( this->Width() == 1 )
            {
                if( this->RowRank() == this->RowAlign() )
                    A.matrix_ = this->LockedMatrix();
                mpi::Broadcast
                ( A.matrix_.Buffer(), A.LocalHeight(), this->RowAlign(), 
                  this->RowComm() );
            }
            else
            {
                const Int rowStride = this->RowStride();
                const Int width = this->Width();
                const Int thisLocalWidth = this->LocalWidth();
                const Int localHeight = this->LocalHeight();
                const Int maxLocalWidth = MaxLength(width,rowStride);

                const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );
                T* buffer = A.auxMemory_.Require( (rowStride+1)*portionSize );
                T* sendBuf = &buffer[0];
                T* recvBuf = &buffer[portionSize];

                // Pack
                const Int ldim = this->LDim();
                const T* thisBuf = this->LockedBuffer();
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* thisCol = &thisBuf[jLoc*ldim];
                    T* sendBufCol = &sendBuf[jLoc*localHeight];
                    MemCopy( sendBufCol, thisCol, localHeight );
                }

                // Communicate
                mpi::AllGather
                ( sendBuf, portionSize, recvBuf, portionSize, this->RowComm() );

                // Unpack
                T* ABuf = A.Buffer();
                const Int ALDim = A.LDim();
                const Int rowAlign = this->RowAlign();
                OUTER_PARALLEL_FOR 
                for( Int k=0; k<rowStride; ++k )
                {
                    const T* data = &recvBuf[k*portionSize];
                    const Int rowShift = Shift_( k, rowAlign, rowStride );
                    const Int localWidth = 
                        Length_( width, rowShift, rowStride );
                    INNER_PARALLEL_FOR
                    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                    {
                        const T* dataCol = &data[jLoc*localHeight];
                        T* ACol = &ABuf[(rowShift+jLoc*rowStride)*ALDim];
                        MemCopy( ACol, dataCol, localHeight );
                    }
                }
                A.auxMemory_.Release();
            }
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( this->Grid().Rank() == 0 )
                std::cerr << "Unaligned RowAllGather." << std::endl;
#endif
            const Int colStride = this->ColStride();
            const Int rowStride = this->RowStride(); 
            const Int colRank = this->ColRank();

            const Int colAlign = this->ColAlign();
            const Int colAlignA = A.ColAlign();
            const Int sendColRank = 
                (colRank+colStride+colAlignA-colAlign) % colStride;
            const Int recvColRank = 
                (colRank+colStride+colAlign-colAlignA) % colStride;

            if( this->Width() == 1 )
            {
                const Int localHeightA = A.LocalHeight();
                if( this->RowRank() == this->RowAlign() )
                {
                    const Int localHeight = this->LocalHeight();    
                    T* buffer = A.auxMemory_.Require( localHeight );

                    // Pack
                    const T* thisCol = this->LockedBuffer();
                    MemCopy( buffer, thisCol, localHeight );

                    // Realign
                    mpi::SendRecv
                    ( buffer, localHeight, sendColRank,
                      A.matrix_.Buffer(), localHeightA, recvColRank, 
                      this->ColComm() );

                    A.auxMemory_.Release();
                }

                // Perform the row broadcast
                mpi::Broadcast
                ( A.matrix_.Buffer(), localHeightA, this->RowAlign(), 
                  this->RowComm() );
            }
            else
            {
                const Int height = this->Height();
                const Int width = this->Width();
                const Int localHeight = this->LocalHeight();
                const Int thisLocalWidth = this->LocalWidth();
                const Int localHeightA = A.LocalHeight();
                const Int maxLocalHeight = MaxLength(height,colStride);
                const Int maxLocalWidth = MaxLength(width,rowStride);

                const Int portionSize = 
                    mpi::Pad( maxLocalHeight*maxLocalWidth );
                T* buffer = A.auxMemory_.Require( (rowStride+1)*portionSize );
                T* firstBuf = &buffer[0];
                T* secondBuf = &buffer[portionSize];

                // Pack
                const Int ldim = this->LDim();
                const T* thisBuf = this->LockedBuffer();
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* thisCol = &thisBuf[jLoc*ldim];
                    T* secondBufCol = &secondBuf[jLoc*localHeight];
                    MemCopy( secondBufCol, thisCol, localHeight );
                }

                // Realign
                mpi::SendRecv
                ( secondBuf, portionSize, sendColRank,
                  firstBuf,  portionSize, recvColRank, this->ColComm() );
            
                // Perform the row AllGather
                mpi::AllGather
                ( firstBuf, portionSize, secondBuf, portionSize, 
                  this->RowComm() );

                // Unpack
                T* ABuf = A.Buffer();
                const Int ALDim = A.LDim();
                const Int rowAlign = this->RowAlign();
                OUTER_PARALLEL_FOR
                for( Int k=0; k<rowStride; ++k )
                {
                    const T* data = &secondBuf[k*portionSize];
                    const Int rowShift = Shift_( k, rowAlign, rowStride );
                    const Int localWidth = 
                        Length_( width, rowShift, rowStride );
                    INNER_PARALLEL_FOR
                    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                    {
                        const T* dataCol = &data[jLoc*localHeightA];
                        T* ACol = &ABuf[(rowShift+jLoc*rowStride)*ALDim]; 
                        MemCopy( ACol, dataCol, localHeightA );
                    }
                }
                A.auxMemory_.Release();
            }
        }
    }
    if( this->Grid().InGrid() && this->CrossComm() != mpi::COMM_SELF )
    {
        // Pack from the root
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        T* buf = A.auxMemory_.Require( localHeight*localWidth );
        if( this->CrossRank() == this->Root() )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &buf[jLoc*localHeight], A.LockedBuffer(0,jLoc), localHeight );
        }

        // Broadcast from the root
        mpi::Broadcast
        ( buf, localHeight*localWidth, this->Root(), this->CrossComm() );

        // Unpack if not the root
        if( this->CrossRank() != this->Root() )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( A.Buffer(0,jLoc), &buf[jLoc*localHeight], localHeight );
        }
        A.auxMemory_.Release();
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AllGather( DistMatrix<T,UGath,VGath>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::AllGather");
        this->AssertSameGrid( A.Grid() );
    )
    A.Resize( this->Height(), this->Width() );

    if( this->Participating() )
    {
        const Int colStride = this->ColStride(); 
        const Int rowStride = this->RowStride();
        const Int distStride = colStride*rowStride;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int thisLocalHeight = this->LocalHeight();
        const Int thisLocalWidth = this->LocalWidth();
        const Int maxLocalHeight = MaxLength(height,colStride);
        const Int maxLocalWidth = MaxLength(width,rowStride);

        const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
        T* buffer = A.auxMemory_.Require( (distStride+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        const Int ldim = this->LDim();
        const T* thisBuf = this->LockedBuffer();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            MemCopy
            ( &sendBuf[jLoc*thisLocalHeight], &thisBuf[jLoc*ldim], 
              thisLocalHeight );

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize, recvBuf, portionSize, this->DistComm() );

        // Unpack
        T* ABuf = A.Buffer();
        const Int ALDim = A.LDim();
        const Int colAlign = this->ColAlign();
        const Int rowAlign = this->RowAlign();
        OUTER_PARALLEL_FOR
        for( Int l=0; l<rowStride; ++l )
        {
            const Int rowShift = Shift_( l, rowAlign, rowStride );
            const Int localWidth = Length_( width, rowShift, rowStride );
            for( Int k=0; k<colStride; ++k )
            {
                const T* data = &recvBuf[(k+l*colStride)*portionSize];
                const Int colShift = Shift_( k, colAlign, colStride );
                const Int localHeight = Length_( height, colShift, colStride );
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                {
                    T* destCol = 
                        &ABuf[colShift+(rowShift+jLoc*rowStride)*ALDim];
                    const T* sourceCol = &data[jLoc*localHeight];
                    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                        destCol[iLoc*colStride] = sourceCol[iLoc];
                }
            }
        }
        A.auxMemory_.Release();
    }
    if( this->Grid().InGrid() && this->CrossComm() != mpi::COMM_SELF )
    {
        // Pack from the root
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        T* buf = A.auxMemory_.Require( localHeight*localWidth );
        if( this->CrossRank() == this->Root() )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &buf[jLoc*localHeight], A.LockedBuffer(0,jLoc), localHeight );
        }

        // Broadcast from the root
        mpi::Broadcast
        ( buf, localHeight*localWidth, this->Root(), this->CrossComm() );

        // Unpack if not the root
        if( this->CrossRank() != this->Root() )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( A.Buffer(0,jLoc), &buf[jLoc*localHeight], localHeight );
        }
        A.auxMemory_.Release();
    }
}

// Diagonal helper functions
// =========================
template<typename T,Dist U,Dist V>
template<typename S,class Function>
void
GeneralDistMatrix<T,U,V>::GetDiagonalHelper
( DistMatrix<S,UDiag,VDiag>& d, Int offset, Function func ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::GetDiagonalHelper"))
    d.SetGrid( this->Grid() );
    d.SetRoot( this->DiagonalRoot(offset) );
    d.AlignCols( this->DiagonalAlign(offset) );
    d.Resize( this->DiagonalLength(offset), 1 );
    if( !d.Participating() )
        return;

    const Int diagShift = d.ColShift();
    const Int diagStride = d.ColStride();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int iLocStart = (iStart-this->ColShift()) / colStride;
    const Int jLocStart = (jStart-this->RowShift()) / rowStride;

    const Int localDiagLength = d.LocalHeight();
    S* dBuf = d.Buffer();
    const T* buffer = this->LockedBuffer();
    const Int ldim = this->LDim();

    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(diagStride/colStride);
        const Int jLoc = jLocStart + k*(diagStride/rowStride);
        func( dBuf[k], buffer[iLoc+jLoc*ldim] );
    }
}

template<typename T,Dist U,Dist V>
template<typename S,class Function>
void
GeneralDistMatrix<T,U,V>::SetDiagonalHelper
( const DistMatrix<S,UDiag,VDiag>& d, Int offset, Function func ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::SetDiagonalHelper");
        if( !this->DiagonalAlignedWith( d, offset ) )
            LogicError("Invalid diagonal alignment");
    )
    if( !d.Participating() )
        return;

    const Int diagShift = d.ColShift();
    const Int diagStride = d.ColStride();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int iLocStart = (iStart-this->ColShift()) / colStride;
    const Int jLocStart = (jStart-this->RowShift()) / rowStride;

    const Int localDiagLength = d.LocalHeight();
    const S* dBuf = d.LockedBuffer();
    T* buffer = this->Buffer();
    const Int ldim = this->LDim();

    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(diagStride/colStride);
        const Int jLoc = jLocStart + k*(diagStride/rowStride);
        func( buffer[iLoc+jLoc*ldim], dBuf[k] );
    }
}

// Instantiations for {Int,Real,Complex<Real>} for each Real in {float,double}
// ###########################################################################

#define DISTPROTO(T,U,V) template class GeneralDistMatrix<T,U,V>
  
#define PROTO(T)\
  DISTPROTO(T,CIRC,CIRC);\
  DISTPROTO(T,MC,  MR  );\
  DISTPROTO(T,MC,  STAR);\
  DISTPROTO(T,MD,  STAR);\
  DISTPROTO(T,MR,  MC  );\
  DISTPROTO(T,MR,  STAR);\
  DISTPROTO(T,STAR,MC  );\
  DISTPROTO(T,STAR,MD  );\
  DISTPROTO(T,STAR,MR  );\
  DISTPROTO(T,STAR,STAR);\
  DISTPROTO(T,STAR,VC  );\
  DISTPROTO(T,STAR,VR  );\
  DISTPROTO(T,VC,  STAR);\
  DISTPROTO(T,VR,  STAR);

#ifndef DISABLE_COMPLEX
 #ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(float);
  PROTO(double);
  PROTO(Complex<float>);
  PROTO(Complex<double>);
 #else // ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(double);
  PROTO(Complex<double>);
 #endif // ifndef DISABLE_FLOAT
#else // ifndef DISABLE_COMPLEX
 #ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(float);
  PROTO(double);
 #else // ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(double);
 #endif // ifndef DISABLE_FLOAT
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
