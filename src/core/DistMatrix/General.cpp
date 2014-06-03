/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include EL_ZEROS_INC

namespace El {

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>::GeneralDistMatrix( const El::Grid& grid, Int root )
: AbstractDistMatrix<T>(grid,root)
{ }

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>::GeneralDistMatrix( GeneralDistMatrix<T,U,V>&& A ) 
EL_NOEXCEPT
: AbstractDistMatrix<T>(std::move(A))
{ }

// Assignment and reconfiguration
// ==============================

template<typename T,Dist U,Dist V>
GeneralDistMatrix<T,U,V>& 
GeneralDistMatrix<T,U,V>::operator=( GeneralDistMatrix<T,U,V>&& A )
{
    AbstractDistMatrix<T>::operator=( std::move(A) );
    return *this;
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AlignColsWith
( const El::DistData& data, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AlignColsWith")) 
    this->SetGrid( *data.grid );
    this->SetRoot( data.root );
    if( data.colDist == U || data.colDist == UPart )
        this->AlignCols( data.colAlign, constrain );
    else if( data.rowDist == U || data.rowDist == UPart )
        this->AlignCols( data.rowAlign, constrain );
    else if( data.colDist == UScat )
        this->AlignCols( data.colAlign % this->ColStride(), constrain );
    else if( data.rowDist == UScat )
        this->AlignCols( data.rowAlign % this->ColStride(), constrain );
    DEBUG_ONLY(
        else if( U != UGath && data.colDist != UGath && data.rowDist != UGath ) 
            LogicError("Nonsensical alignment");
    )
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AlignRowsWith
( const El::DistData& data, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AlignRowsWith")) 
    this->SetGrid( *data.grid );
    this->SetRoot( data.root );
    if( data.colDist == V || data.colDist == VPart )
        this->AlignRows( data.colAlign, constrain );
    else if( data.rowDist == V || data.rowDist == VPart )
        this->AlignRows( data.rowAlign, constrain );
    else if( data.colDist == VScat )
        this->AlignRows( data.colAlign % this->RowStride(), constrain );
    else if( data.rowDist == VScat )
        this->AlignRows( data.rowAlign % this->RowStride(), constrain );
    DEBUG_ONLY(
        else if( V != VGath && data.colDist != VGath && data.rowDist != VGath ) 
            LogicError("Nonsensical alignment");
    )
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::Translate( DistMatrix<T,U,V>& A ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::Translate"))
    const Grid& g = this->Grid();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int colAlign = this->ColAlign();
    const Int rowAlign = this->RowAlign();
    const Int root = this->Root();
    A.SetGrid( g );
    if( !A.RootConstrained() )
        A.SetRoot( root );
    if( !A.ColConstrained() )
        A.AlignCols( colAlign, false );
    if( !A.RowConstrained() )
        A.AlignRows( rowAlign, false );
    A.Resize( height, width );
    if( !g.InGrid() )
        return;

    const bool aligned = colAlign == A.ColAlign() && rowAlign == A.RowAlign();
    if( aligned && root == A.Root() )
    {
        A.matrix_ = this->matrix_;
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [U,V] <- [U,V]" << std::endl;
#endif
        const Int colRank = this->ColRank();
        const Int rowRank = this->RowRank();
        const Int crossRank = this->CrossRank();
        const Int colStride = this->ColStride();
        const Int rowStride = this->RowStride();
        const Int maxHeight = MaxLength( height, colStride );
        const Int maxWidth  = MaxLength( width,  rowStride );
        const Int pkgSize = mpi::Pad( maxHeight*maxWidth );
        T* buffer=0;
        if( crossRank == root || crossRank == A.Root() )
            buffer = A.auxMemory_.Require( pkgSize );

        const Int colAlignA = A.ColAlign();
        const Int rowAlignA = A.RowAlign();
        const Int localHeightA = 
            Length( height, colRank, colAlignA, colStride );
        const Int localWidthA = Length( width, rowRank, rowAlignA, rowStride );
        const Int recvSize = mpi::Pad( localHeightA*localWidthA );

        if( crossRank == root )
        {
            // Pack the local data
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &buffer[jLoc*localHeight], this->LockedBuffer(0,jLoc), 
                  localHeight );

            if( !aligned )
            {
                // If we were not aligned, then SendRecv over the DistComm
                const Int toRow = Mod(colRank+colAlignA-colAlign,colStride);
                const Int toCol = Mod(rowRank+rowAlignA-rowAlign,rowStride);
                const Int fromRow = Mod(colRank+colAlign-colAlignA,colStride);
                const Int fromCol = Mod(rowRank+rowAlign-rowAlignA,rowStride);
                const Int toRank = toRow + toCol*colStride;
                const Int fromRank = fromRow + fromCol*colStride;
                mpi::SendRecv
                ( buffer, pkgSize, toRank, fromRank, this->DistComm() );
            }
        }
        if( root != A.Root() )
        {
            // Send to the correct new root over the cross communicator
            if( crossRank == root )
                mpi::Send( buffer, recvSize, A.Root(), A.CrossComm() );
            else if( crossRank == A.Root() )
                mpi::Recv( buffer, recvSize, root, A.CrossComm() );
        }
        // Unpack
        if( crossRank == A.Root() )
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
                MemCopy
                ( A.Buffer(0,jLoc), &buffer[jLoc*localHeightA], localHeightA );
        if( crossRank == root || crossRank == A.Root() )
            A.auxMemory_.Release();
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AllGather( DistMatrix<T,UGath,VGath>& A ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AllGather"))
    const Int height = this->Height();
    const Int width = this->Width();
    A.SetGrid( this->Grid() );
    A.Resize( height, width );

    if( this->Participating() )
    {
        const Int colStride = this->ColStride(); 
        const Int rowStride = this->RowStride();
        const Int distStride = colStride*rowStride;

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
        EL_PARALLEL_FOR
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
        EL_OUTER_PARALLEL_FOR
        for( Int l=0; l<rowStride; ++l )
        {
            const Int rowShift = Shift_( l, rowAlign, rowStride );
            const Int localWidth = Length_( width, rowShift, rowStride );
            for( Int k=0; k<colStride; ++k )
            {
                const T* data = &recvBuf[(k+l*colStride)*portionSize];
                const Int colShift = Shift_( k, colAlign, colStride );
                const Int localHeight = Length_( height, colShift, colStride );
                EL_INNER_PARALLEL_FOR
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

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::ColAllGather( DistMatrix<T,UGath,V>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::ColAllGather");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = this->Height();
    const Int width = this->Width();
#ifdef EL_CACHE_WARNINGS
    if( height != 1 && this->Grid().Rank() == 0 )
    {
        std::cerr <<
          "The matrix redistribution [* ,V] <- [U,V] potentially causes a "
          "large amount of cache-thrashing. If possible, avoid it by "
          "performing the redistribution with a (conjugate)-transpose"
          << std::endl;
    }
#endif
    A.AlignRowsAndResize( this->RowAlign(), height, width, false, false );

    if( this->Participating() )
    {
        if( this->RowAlign() == A.RowAlign() )
        {
            if( height == 1 )
            {
                const Int localWidthA = A.LocalWidth();
                T* bcastBuf = A.auxMemory_.Require( localWidthA );
            
                if( this->ColRank() == this->ColAlign() )
                {
                    A.matrix_ = this->LockedMatrix();
                    // Pack
                    const T* ABuf = A.LockedBuffer(); 
                    const Int ALDim = A.LDim();
                    EL_PARALLEL_FOR
                    for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
                        bcastBuf[jLoc] = ABuf[jLoc*ALDim];
                }

                // Broadcast within the column comm
                mpi::Broadcast
                ( bcastBuf, localWidthA, this->ColAlign(), this->ColComm() );

                // Unpack
                T* ABuf = A.Buffer();
                const Int ALDim = A.LDim(); 
                EL_PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
                    ABuf[jLoc*ALDim] = bcastBuf[jLoc];
                A.auxMemory_.Release();
            }
            else
            {
                const Int colStride = this->ColStride();
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
                EL_PARALLEL_FOR
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
                EL_OUTER_PARALLEL_FOR
                for( Int k=0; k<colStride; ++k )
                {
                    const T* data = &recvBuf[k*portionSize];
                    const Int colShift = Shift_( k, colAlign, colStride );
                    const Int localHeight = 
                        Length_( height, colShift, colStride );
                    EL_INNER_PARALLEL_FOR
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
#ifdef EL_UNALIGNED_WARNINGS
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

            if( height == 1 )
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
                    EL_PARALLEL_FOR
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
                EL_PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
                    ABuf[jLoc*ALDim] = bcastBuf[jLoc];
                A.auxMemory_.Release();
            }
            else
            {
                const Int thisLocalWidth = this->LocalWidth();
                const Int localWidthA = A.LocalWidth();
                const Int thisLocalHeight = this->LocalHeight();
                const Int maxLocalHeight = MaxLength(height,colStride);
                const Int maxLocalWidth = MaxLength(width,rowStride);
                const Int portionSize = 
                    mpi::Pad( maxLocalHeight*maxLocalWidth );

                T* buffer = A.auxMemory_.Require( (colStride+1)*portionSize );
                T* firstBuf  = &buffer[0];
                T* secondBuf = &buffer[portionSize];

                // Pack
                const Int ldim = this->LDim();
                const T* thisBuf = this->LockedBuffer();
                EL_PARALLEL_FOR
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
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
                ( firstBuf, portionSize, 
                  secondBuf, portionSize, this->ColComm() ); 

                // Unpack the contents of each member of the column team
                T* ABuf = A.Buffer();
                const Int ALDim = A.LDim();
                const Int colAlign = this->ColAlign();
                EL_OUTER_PARALLEL_FOR
                for( Int k=0; k<colStride; ++k )
                {
                    const T* data = &secondBuf[k*portionSize];
                    const Int colShift = Shift_( k, colAlign, colStride );
                    const Int localHeight = 
                        Length_( height, colShift, colStride );
                    EL_INNER_PARALLEL_FOR
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
    const Int height = this->Height();
    const Int width = this->Width();
    A.AlignColsAndResize( this->ColAlign(), height, width, false, false );

    if( this->Participating() )
    {
        if( this->ColAlign() == A.ColAlign() )
        {
            if( width == 1 )
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
                EL_PARALLEL_FOR
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
                EL_OUTER_PARALLEL_FOR 
                for( Int k=0; k<rowStride; ++k )
                {
                    const T* data = &recvBuf[k*portionSize];
                    const Int rowShift = Shift_( k, rowAlign, rowStride );
                    const Int localWidth = 
                        Length_( width, rowShift, rowStride );
                    EL_INNER_PARALLEL_FOR
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
#ifdef EL_UNALIGNED_WARNINGS
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

            if( width == 1 )
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
                EL_PARALLEL_FOR
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
                ( firstBuf,  portionSize, 
                  secondBuf, portionSize, this->RowComm() );

                // Unpack
                T* ABuf = A.Buffer();
                const Int ALDim = A.LDim();
                const Int rowAlign = this->RowAlign();
                EL_OUTER_PARALLEL_FOR
                for( Int k=0; k<rowStride; ++k )
                {
                    const T* data = &secondBuf[k*portionSize];
                    const Int rowShift = Shift_( k, rowAlign, rowStride );
                    const Int localWidth = 
                        Length_( width, rowShift, rowStride );
                    EL_INNER_PARALLEL_FOR
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
GeneralDistMatrix<T,U,V>::PartialColAllGather( DistMatrix<T,UPart,V>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialColAllGather");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = this->Height();
    const Int width = this->Width();
#ifdef EL_VECTOR_WARNINGS
    if( width == 1 && this->Grid().Rank() == 0 )
    {
        std::cerr <<
          "The vector version of PartialColAllGather is not yet written but "
          "would only require modifying the vector version of "
          "PartialRowAllGather" << std::endl;
    }
#endif
#ifdef EL_CACHE_WARNINGS
    if( width && this->Grid().Rank() == 0 )
    {
        std::cerr <<
          "PartialColAllGather potentially causes a large amount of cache-"
          "thrashing. If possible, avoid it by performing the redistribution"
          "on the (conjugate-)transpose" << std::endl;
    }
#endif
    A.AlignColsAndResize
    ( this->ColAlign()%A.ColStride(), height, width, false, false );
    if( !this->Participating() )
        return;

    DEBUG_ONLY(
        if( this->LocalWidth() != this->Width() )
            LogicError("This routine assumes rows are not distributed");
    )
    const T* thisBuf = this->LockedBuffer();
    const Int ldim = this->LDim();
    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();

    const Int colAlign = this->ColAlign();
    const Int colAlignA = A.ColAlign();
    const Int colStride = this->ColStride();
    const Int colStrideUnion = this->PartialUnionColStride();
    const Int colStridePart = this->PartialColStride();
    const Int colRankPart = this->PartialColRank();
    const Int colShiftA = A.ColShift();

    const Int thisLocalHeight = this->LocalHeight();
    const Int maxLocalHeight = MaxLength(height,colStride);
    const Int portionSize = mpi::Pad( maxLocalHeight*width );
    T* buffer = A.auxMemory_.Require( (colStrideUnion+1)*portionSize );
    T* firstBuf = &buffer[0];
    T* secondBuf = &buffer[portionSize];

    if( colAlignA == colAlign % colStridePart ) 
    {
        // Pack
        EL_PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* thisCol = &thisBuf[j*ldim];
            T* firstBufCol = &firstBuf[j*thisLocalHeight];
            MemCopy( firstBufCol, thisCol, thisLocalHeight );
        }

        // Communicate
        mpi::AllGather
        ( firstBuf, portionSize, secondBuf, portionSize, 
          this->PartialUnionColComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShift = 
                Shift_( colRankPart+k*colStridePart, colAlign, colStride );
            const Int colOffset = (colShift-colShiftA) / colStridePart;
            const Int localHeight = Length_( height, colShift, colStride );
            EL_INNER_PARALLEL_FOR
            for( Int j=0; j<width; ++j )
            {
                const T* dataCol = &data[j*localHeight];
                T* ACol = &ABuf[colOffset+j*ALDim];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    ACol[iLoc*colStrideUnion] = dataCol[iLoc];
            }
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned PartialColAllGather" << std::endl;
#endif
        // Perform a SendRecv to match the row alignments
        const Int colRank = this->ColRank();
        const Int sendColRank = 
            (colRank+colStride+colAlignA-colAlign) % colStride;
        const Int recvColRank = 
            (colRank+colStride+colAlign-colAlignA) % colStride;
        EL_PARALLEL_FOR
        for( Int j=0; j<width; ++j ) 
        {
            const T* thisCol = &thisBuf[j*ldim];
            T* secondBufCol = &secondBuf[j*thisLocalHeight];
            MemCopy( secondBufCol, thisCol, thisLocalHeight );
        }
        mpi::SendRecv
        ( secondBuf, portionSize, sendColRank,
          firstBuf,  portionSize, recvColRank, this->ColComm() );

        // Use the SendRecv as an input to the partial union AllGather
        mpi::AllGather
        ( firstBuf,  portionSize, 
          secondBuf, portionSize, this->PartialUnionColComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShift = 
                Shift_( colRankPart+colStridePart*k, colAlignA, colStride );
            const Int colOffset = (colShift-colShiftA) / colStridePart;
            const Int localHeight = Length_( height, colShift, colStride );
            EL_INNER_PARALLEL_FOR
            for( Int j=0; j<width; ++j )
            {
                const T* dataCol = &data[j*localHeight];
                T* ACol = &ABuf[colOffset+j*ALDim];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    ACol[iLoc*colStrideUnion] = dataCol[iLoc];
            }
        }
    }
    A.auxMemory_.Release();
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::PartialRowAllGather( DistMatrix<T,U,VPart>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialRowAllGather");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = this->Height();
    const Int width = this->Width();
    A.AlignRowsAndResize
    ( this->RowAlign()%A.RowStride(), height, width, false, false );
    if( !this->Participating() )
        return;

    DEBUG_ONLY(
        if( this->LocalHeight() != this->Height() )
            LogicError("This routine assumes columns are not distributed");
    )
    const T* thisBuf = this->LockedBuffer();
    const Int ldim = this->LDim();
    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();

    const Int rowAlign = this->RowAlign();
    const Int rowAlignA = A.RowAlign();
    const Int rowStride = this->RowStride();
    const Int rowStrideUnion = this->PartialUnionRowStride();
    const Int rowStridePart = this->PartialRowStride();
    const Int rowRankPart = this->PartialRowRank();
    const Int rowShiftA = A.RowShift();

    const Int thisLocalWidth = this->LocalWidth();
    const Int maxLocalWidth = MaxLength(width,rowStride);
    const Int portionSize = mpi::Pad( height*maxLocalWidth );
    T* buffer = A.auxMemory_.Require( (rowStrideUnion+1)*portionSize );
    T* firstBuf = &buffer[0];
    T* secondBuf = &buffer[portionSize];

    if( rowAlignA == rowAlign % rowStridePart ) 
    {
        // Pack
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
        {
            const T* thisCol = &thisBuf[jLoc*ldim];
            T* firstBufCol = &firstBuf[jLoc*height];
            MemCopy( firstBufCol, thisCol, height );
        }

        // Communicate
        mpi::AllGather
        ( firstBuf, portionSize, secondBuf, portionSize, 
          this->PartialUnionRowComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShift = 
                Shift_( rowRankPart+k*rowStridePart, rowAlign, rowStride );
            const Int rowOffset = (rowShift-rowShiftA) / rowStridePart;
            const Int localWidth = Length_( width, rowShift, rowStride );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* ACol = &ABuf[(rowOffset+jLoc*rowStrideUnion)*ALDim];
                MemCopy( ACol, dataCol, height );
            }
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned PartialRowAllGather" << std::endl;
#endif
        // Perform a SendRecv to match the row alignments
        const Int rowRank = this->RowRank();
        const Int sendRowRank = 
            (rowRank+rowStride+rowAlignA-rowAlign) % rowStride;
        const Int recvRowRank = 
            (rowRank+rowStride+rowAlign-rowAlignA) % rowStride;
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc ) 
        {
            const T* thisCol = &thisBuf[jLoc*ldim];
            T* secondBufCol = &secondBuf[jLoc*height];
            MemCopy( secondBufCol, thisCol, height );
        }
        mpi::SendRecv
        ( secondBuf, portionSize, sendRowRank,
          firstBuf,  portionSize, recvRowRank, this->RowComm() );

        // Use the SendRecv as an input to the partial union AllGather
        mpi::AllGather
        ( firstBuf,  portionSize, 
          secondBuf, portionSize, this->PartialUnionRowComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShift = 
                Shift_( rowRankPart+rowStridePart*k, rowAlignA, rowStride );
            const Int rowOffset = (rowShift-rowShiftA) / rowStridePart;
            const Int localWidth = Length_( width, rowShift, rowStride );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* ACol = &ABuf[(rowOffset+jLoc*rowStrideUnion)*ALDim];
                MemCopy( ACol, dataCol, height );
            }
        }
    }
    A.auxMemory_.Release();
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::FilterFrom( const DistMatrix<T,UGath,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::FilterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = A.Height();
    const Int width = A.Width();
    this->Resize( height, width );
    if( !this->Participating() )
        return;

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    
    T* thisBuf = this->Buffer();
    const Int ldim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    EL_PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        T* thisCol = &thisBuf[jLoc*ldim];
        const T* ACol = &ABuf[colShift+(rowShift+jLoc*rowStride)*ALDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            thisCol[iLoc] = ACol[iLoc*colStride];
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::ColFilterFrom( const DistMatrix<T,UGath,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::ColFilterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = A.Height();
    const Int width = A.Width();
    this->AlignRowsAndResize( A.RowAlign(), height, width, false, false );
    if( !this->Participating() )
        return;

    const Int colStride = this->ColStride();
    const Int colShift = this->ColShift();
    const Int rowAlign = this->RowAlign();
    const Int rowAlignA = A.RowAlign();

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    
    T* thisBuf = this->Buffer();
    const Int ldim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    if( rowAlign == rowAlignA )
    {
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            T* thisCol = &thisBuf[jLoc*ldim];
            const T* ACol = &ABuf[colShift+jLoc*ALDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                thisCol[iLoc] = ACol[iLoc*colStride];
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned ColFilterFrom" << std::endl;
#endif
        const Int rowStride = this->RowStride();
        const Int rowRank = this->RowRank();
        const Int sendRowRank = 
            (rowRank+rowStride+rowAlign-rowAlignA) % rowStride;
        const Int recvRowRank = 
            (rowRank+rowStride+rowAlignA-rowAlign) % rowStride;
        const Int localWidthA = A.LocalWidth();
        const Int sendSize = localHeight*localWidthA;
        const Int recvSize = localHeight*localWidth;
        T* buffer = this->auxMemory_.Require( sendSize+recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];
        
        // Pack
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
        {
            T* sendCol = &sendBuf[jLoc*localHeight];
            const T* ACol = &ABuf[colShift+jLoc*ALDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                sendCol[iLoc] = ACol[iLoc*colStride];
        }

        // Realign
        mpi::SendRecv
        ( sendBuf, sendSize, sendRowRank,
          recvBuf, recvSize, recvRowRank, this->RowComm() );

        // Unpack
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuf[jLoc*ldim], &recvBuf[jLoc*localHeight], localHeight );
        this->auxMemory_.Release();
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::RowFilterFrom( const DistMatrix<T,U,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::RowFilterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = A.Height();
    const Int width = A.Width();
    this->AlignColsAndResize( A.ColAlign(), height, width, false, false );
    if( !this->Participating() )
        return;

    const Int colAlign = this->ColAlign();
    const Int colAlignA = A.ColAlign();
    const Int rowStride = this->RowStride();
    const Int rowShift = this->RowShift();

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    
    T* thisBuf = this->Buffer();
    const Int ldim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    
    if( colAlign == colAlignA )
    {
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            T* thisCol = &thisBuf[jLoc*ldim];
            const T* ACol = &ABuf[(rowShift+jLoc*rowStride)*ALDim];
            MemCopy( thisCol, ACol, localHeight );
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned RowFilterFrom" << std::endl;
#endif
        const Int colRank = this->ColRank();
        const Int colStride = this->ColStride();
        const Int sendColRank = 
            (colRank+colStride+colAlign-colAlignA) % colStride;
        const Int recvColRank = 
            (colRank+colStride+colAlignA-colAlign) % colStride;
        const Int localHeightA = A.LocalHeight();
        const Int sendSize = localHeightA*localWidth;
        const Int recvSize = localHeight *localWidth;

        T* buffer = this->auxMemory_.Require( sendSize+recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &sendBuf[jLoc*localHeightA],
              &ABuf[(rowShift+jLoc*rowStride)*ALDim], localHeightA );

        // Realign
        mpi::SendRecv
        ( sendBuf, sendSize, sendColRank, 
          recvBuf, recvSize, recvColRank, this->ColComm() );

        // Unpack
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuf[jLoc*ldim], &recvBuf[jLoc*localHeight], localHeight );
        this->auxMemory_.Release();
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::PartialColFilterFrom( const DistMatrix<T,UPart,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialColFilterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = A.Height();
    const Int width = A.Width();
    this->AlignColsAndResize( A.ColAlign(), height, width, false, false );
    if( !this->Participating() )
        return;

    const Int colAlign = this->ColAlign();
    const Int colAlignA = A.ColAlign();
    const Int colStride = this->ColStride();
    const Int colStridePart = this->PartialColStride();
    const Int colStrideUnion = this->PartialUnionColStride();
    const Int colShiftA = A.ColShift();

    const Int localHeight = this->LocalHeight();

    T* thisBuf = this->Buffer();
    const Int ldim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    if( colAlign % colStridePart == colAlignA )
    {
        const Int colShift = this->ColShift();
        const Int colOffset = (colShift-colShiftA) / colStridePart;
        EL_PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            T* thisCol = &thisBuf[j*ldim];
            const T* ACol = &ABuf[colOffset+j*ALDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                thisCol[iLoc] = ACol[iLoc*colStrideUnion];
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned PartialColFilterFrom" << std::endl;
#endif
        const Int colRankPart = this->PartialColRank();
        const Int colRankUnion = this->PartialUnionColRank();
        const Int colShiftA = A.ColShift();

        // Realign
        // -------
        const Int sendColRankPart = 
            (colRankPart+colStridePart+(colAlign%colStridePart)-colAlignA) % 
            colStridePart;
        const Int recvColRankPart =
            (colRankPart+colStridePart+colAlignA-(colAlign%colStridePart)) %
            colStridePart;
        const Int sendColRank = sendColRankPart + colStridePart*colRankUnion;
        const Int sendColShift = Shift( sendColRank, colAlign, colStride );
        const Int sendColOffset = (sendColShift-colShiftA) / colStridePart;
        const Int localHeightSend = Length( height, sendColShift, colStride );
        const Int sendSize = localHeightSend*width;
        const Int recvSize = localHeight    *width;
        T* buffer = this->auxMemory_.Require( sendSize+recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];
        // Pack
        EL_PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            T* sendCol = &sendBuf[j*localHeightSend];
            const T* ACol = &ABuf[sendColOffset+j*ALDim];
            for( Int iLoc=0; iLoc<localHeightSend; ++iLoc )
                sendCol[iLoc] = ACol[iLoc*colStrideUnion];
        }
        // Change the column alignment
        mpi::SendRecv
        ( sendBuf, sendSize, sendColRankPart,
          recvBuf, recvSize, recvColRankPart, this->PartialColComm() );

        // Unpack
        // ------
        EL_PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* recvCol = &recvBuf[j*localHeight];
            T* thisCol = &thisBuf[j*ldim];
            MemCopy( thisCol, recvCol, localHeight );
        }
        this->auxMemory_.Release();
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::PartialRowFilterFrom( const DistMatrix<T,U,VPart>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialRowFilterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = A.Height();
    const Int width = A.Width();
    this->AlignRowsAndResize( A.RowAlign(), height, width, false, false );
    if( !this->Participating() )
        return;

    const Int rowAlign = this->RowAlign();
    const Int rowAlignA = A.RowAlign();
    const Int rowStride = this->RowStride();
    const Int rowStridePart = this->PartialRowStride();
    const Int rowStrideUnion = this->PartialUnionRowStride();
    const Int rowShiftA = A.RowShift();

    const Int localWidth = this->LocalWidth();

    T* thisBuf = this->Buffer();
    const Int ldim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    if( rowAlign % rowStridePart == rowAlignA )
    {
        const Int rowShift = this->RowShift();
        const Int rowOffset = (rowShift-rowShiftA) / rowStridePart;
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            T* thisCol = &thisBuf[jLoc*ldim];
            const T* ACol = &ABuf[(rowOffset+jLoc*rowStrideUnion)*ALDim];
            MemCopy( thisCol, ACol, height );
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned PartialRowFilterFrom" << std::endl;
#endif
        const Int rowRankPart = this->PartialRowRank();
        const Int rowRankUnion = this->PartialUnionRowRank();
        const Int rowShiftA = A.RowShift();

        // Realign
        // -------
        const Int sendRowRankPart = 
            (rowRankPart+rowStridePart+(rowAlign%rowStridePart)-rowAlignA) % 
            rowStridePart;
        const Int recvRowRankPart =
            (rowRankPart+rowStridePart+rowAlignA-(rowAlign%rowStridePart)) %
            rowStridePart;
        const Int sendRowRank = sendRowRankPart + rowStridePart*rowRankUnion;
        const Int sendRowShift = Shift( sendRowRank, rowAlign, rowStride );
        const Int sendRowOffset = (sendRowShift-rowShiftA) / rowStridePart;
        const Int localWidthSend = Length( width, sendRowShift, rowStride );
        const Int sendSize = height*localWidthSend;
        const Int recvSize = height*localWidth;
        T* buffer = this->auxMemory_.Require( sendSize+recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];
        // Pack
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidthSend; ++jLoc )
        {
            T* sendCol = &sendBuf[jLoc*height];
            const T* ACol = &ABuf[(sendRowOffset+jLoc*rowStrideUnion)*ALDim];
            MemCopy( sendCol, ACol, height );
        }
        // Change the column alignment
        mpi::SendRecv
        ( sendBuf, sendSize, sendRowRankPart,
          recvBuf, recvSize, recvRowRankPart, this->PartialRowComm() );

        // Unpack
        // ------
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* recvCol = &recvBuf[jLoc*height];
            T* thisCol = &thisBuf[jLoc*ldim];
            MemCopy( thisCol, recvCol, height );
        }
        this->auxMemory_.Release();
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::PartialColAllToAllFrom
( const DistMatrix<T,UPart,VScat>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialColAllToAllFrom");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = A.Height();
    const Int width = A.Width();
    this->AlignColsAndResize( A.ColAlign(), height, width, false, false );
    if( !this->Participating() )
        return;

    const Int colAlign = this->ColAlign();
    const Int colAlignA = A.ColAlign();
    const Int rowAlignA = A.RowAlign();

    const Int colStride = this->ColStride();
    const Int colStridePart = this->PartialColStride();
    const Int colStrideUnion = this->PartialUnionColStride();
    const Int colRankPart = this->PartialColRank();

    const Int colShiftA = A.ColShift();

    const Int thisLocalHeight = this->LocalHeight();
    const Int localWidthA = A.LocalWidth();
    const Int maxLocalHeight = MaxLength(height,colStride);
    const Int maxLocalWidth = MaxLength(width,colStrideUnion);
    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

    T* thisBuf = this->Buffer();
    const Int ldim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    T* buffer = this->auxMemory_.Require( 2*colStrideUnion*portionSize );
    T* firstBuf  = &buffer[0];
    T* secondBuf = &buffer[colStrideUnion*portionSize];

    if( colAlign % colStridePart == colAlignA )
    {
        // Pack            
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            T* data = &firstBuf[k*portionSize];
            const Int colRank = colRankPart + k*colStridePart;
            const Int colShift = Shift_( colRank, colAlign, colStride );
            const Int colOffset = (colShift-colShiftA) / colStridePart;
            const Int localHeight = Length_( height, colShift, colStride );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
            {
                T* dataCol = &data[jLoc*localHeight];
                const T* ACol = &ABuf[colOffset+jLoc*ALDim];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    dataCol[iLoc] = ACol[iLoc*colStrideUnion];
            }
        }

        // Simultaneously Scatter in columns and Gather in rows
        mpi::AllToAll
        ( firstBuf,  portionSize, 
          secondBuf, portionSize, this->PartialUnionColComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShift = Shift_( k, rowAlignA, colStrideUnion );
            const Int localWidth = Length_( width, rowShift, colStrideUnion );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*thisLocalHeight];
                T* thisCol = &thisBuf[(rowShift+jLoc*colStrideUnion)*ldim]; 
                MemCopy( thisCol, dataCol, thisLocalHeight );
            }
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned PartialColAllToAllFrom" << std::endl;
#endif
        const Int sendColRankPart = 
            (colRankPart+colStridePart+(colAlign%colStridePart)-colAlignA) % 
            colStridePart;
        const Int recvColRankPart =
            (colRankPart+colStridePart+colAlignA-(colAlign%colStridePart)) %
            colStridePart; 

        // Pack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int colRank = sendColRankPart + k*colStridePart;
            const Int colShift = Shift_( colRank, colAlign, colStride );
            const Int colOffset = (colShift-colShiftA) / colStridePart;
            const Int localHeight = Length_( height, colShift, colStride );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
            {
                T* dataCol = &data[jLoc*localHeight];
                const T* ACol = &ABuf[colOffset+jLoc*ALDim];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    dataCol[iLoc] = ACol[iLoc*colStrideUnion];
            }
        }

        // Simultaneously Scatter in columns and Gather in rows
        mpi::AllToAll
        ( secondBuf, portionSize, 
          firstBuf,  portionSize, this->PartialUnionColComm() );

        // Realign the result
        mpi::SendRecv 
        ( firstBuf,  colStrideUnion*portionSize, sendColRankPart,
          secondBuf, colStrideUnion*portionSize, recvColRankPart, 
          this->PartialColComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShift = Shift_( k, rowAlignA, colStrideUnion );
            const Int localWidth = Length_( width, rowShift, colStrideUnion );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*thisLocalHeight];
                T* thisCol = &thisBuf[(rowShift+jLoc*colStrideUnion)*ldim]; 
                MemCopy( thisCol, dataCol, thisLocalHeight );
            }
        }
    }
    this->auxMemory_.Release();
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::PartialRowAllToAllFrom
( const DistMatrix<T,UScat,VPart>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialRowAllToAllFrom");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = A.Height();
    const Int width = A.Width();
    this->AlignRowsAndResize( A.RowAlign(), height, width, false, false );
    if( !this->Participating() )
        return;

    const Int rowAlign = this->RowAlign();
    const Int rowAlignA = A.RowAlign();
    const Int colAlignA = A.ColAlign();

    const Int rowStride = this->RowStride();
    const Int rowStridePart = this->PartialRowStride();
    const Int rowStrideUnion = this->PartialUnionRowStride();
    const Int rowRankPart = this->PartialRowRank();

    const Int rowShiftA = A.RowShift();

    const Int thisLocalWidth = this->LocalWidth();
    const Int localHeightA = A.LocalHeight();
    const Int maxLocalHeight = MaxLength(height,rowStrideUnion);
    const Int maxLocalWidth = MaxLength(width,rowStride);
    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

    T* thisBuf = this->Buffer();
    const Int ldim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    T* buffer = this->auxMemory_.Require( 2*rowStrideUnion*portionSize );
    T* firstBuf  = &buffer[0];
    T* secondBuf = &buffer[rowStrideUnion*portionSize];

    if( rowAlign % rowStridePart == rowAlignA )
    {
        // Pack            
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            T* data = &firstBuf[k*portionSize];
            const Int rowRank = rowRankPart + k*rowStridePart;
            const Int rowShift = Shift_( rowRank, rowAlign, rowStride );
            const Int rowOffset = (rowShift-rowShiftA) / rowStridePart;
            const Int localWidth = Length_( width, rowShift, rowStride );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* dataCol = &data[jLoc*localHeightA];
                const T* ACol = &ABuf[(rowOffset+jLoc*rowStrideUnion)*ALDim];
                MemCopy( dataCol, ACol, localHeightA );
            }
        }

        // Simultaneously Scatter in rows and Gather in columns
        mpi::AllToAll
        ( firstBuf,  portionSize, 
          secondBuf, portionSize, this->PartialUnionRowComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShift = Shift_( k, colAlignA, rowStrideUnion );
            const Int localHeight = Length_( height, colShift, rowStrideUnion );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuf[colShift+jLoc*ldim]; 
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    thisCol[iLoc*rowStrideUnion] = dataCol[iLoc];
            }
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned PartialRowAllToAllFrom" << std::endl;
#endif
        const Int sendRowRankPart = 
            (rowRankPart+rowStridePart+(rowAlign%rowStridePart)-rowAlignA) % 
            rowStridePart;
        const Int recvRowRankPart =
            (rowRankPart+rowStridePart+rowAlignA-(rowAlign%rowStridePart)) %
            rowStridePart; 

        // Pack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int rowRank = sendRowRankPart + k*rowStridePart;
            const Int rowShift = Shift_( rowRank, rowAlign, rowStride );
            const Int rowOffset = (rowShift-rowShiftA) / rowStridePart;
            const Int localWidth = Length_( width, rowShift, rowStride );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* dataCol = &data[jLoc*localHeightA];
                const T* ACol = &ABuf[(rowOffset+jLoc*rowStrideUnion)*ALDim];
                MemCopy( dataCol, ACol, localHeightA );
            }
        }

        // Simultaneously Scatter in rows and Gather in columns
        mpi::AllToAll
        ( secondBuf, portionSize, 
          firstBuf,  portionSize, this->PartialUnionRowComm() );

        // Realign the result
        mpi::SendRecv 
        ( firstBuf,  rowStrideUnion*portionSize, sendRowRankPart,
          secondBuf, rowStrideUnion*portionSize, recvRowRankPart, 
          this->PartialRowComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShift = Shift_( k, colAlignA, rowStrideUnion );
            const Int localHeight = Length_( height, colShift, rowStrideUnion );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuf[colShift+jLoc*ldim]; 
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    thisCol[iLoc*rowStrideUnion] = dataCol[iLoc];
            }
        }
    }
    this->auxMemory_.Release();
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::PartialColAllToAll
( DistMatrix<T,UPart,VScat>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialColAllToAll");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = this->Height();
    const Int width = this->Width();
    A.AlignColsAndResize
    ( this->ColAlign()%A.ColStride(), height, width, false, false );
    if( !A.Participating() )
        return;

    const Int colAlign = this->ColAlign();
    const Int colAlignA = A.ColAlign();
    const Int rowAlignA = A.RowAlign();

    const Int colStride = this->ColStride();
    const Int colStridePart = this->PartialColStride();
    const Int colStrideUnion = this->PartialUnionColStride();
    const Int colRankPart = this->PartialColRank();

    const Int colShiftA = A.ColShift();

    const Int thisLocalHeight = this->LocalHeight();
    const Int localWidthA = A.LocalWidth();
    const Int maxLocalHeight = MaxLength(height,colStride);
    const Int maxLocalWidth = MaxLength(width,colStrideUnion);
    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

    const T* thisBuf = this->LockedBuffer();
    const Int ldim = this->LDim();
    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();

    T* buffer = A.auxMemory_.Require( 2*colStrideUnion*portionSize );
    T* firstBuf  = &buffer[0];
    T* secondBuf = &buffer[colStrideUnion*portionSize];

    if( colAlignA == colAlign % colStridePart )
    {
        // Pack            
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            T* data = &firstBuf[k*portionSize];
            const Int rowShift = Shift_( k, rowAlignA, colStrideUnion );
            const Int localWidth = Length_( width, rowShift, colStrideUnion );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &data[jLoc*thisLocalHeight],
                  &thisBuf[(rowShift+jLoc*colStrideUnion)*ldim], 
                  thisLocalHeight );
        }

        // Simultaneously Gather in columns and Scatter in rows
        mpi::AllToAll
        ( firstBuf,  portionSize, 
          secondBuf, portionSize, this->PartialUnionColComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colRank = colRankPart + k*colStridePart;
            const Int colShift = Shift_( colRank, colAlign, colStride );
            const Int colOffset = (colShift-colShiftA) / colStridePart;
            const Int localHeight = Length_( height, colShift, colStride );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
            {
                T* ACol = &ABuf[colOffset+jLoc*ALDim];
                const T* dataCol = &data[jLoc*localHeight];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    ACol[iLoc*colStrideUnion] = dataCol[iLoc];
            }
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned PartialColAllToAll" << std::endl;
#endif
        const Int colAlignDiff = colAlignA - (colAlign%colStridePart);
        const Int sendColRankPart = 
            (colRankPart+colStridePart+colAlignDiff) % colStridePart;
        const Int recvColRankPart =
            (colRankPart+colStridePart-colAlignDiff) % colStridePart;

        // Pack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            T* data = &secondBuf[k*portionSize];    
            const Int rowShift = Shift_( k, rowAlignA, colStrideUnion );
            const Int localWidth = Length_( width, rowShift, colStrideUnion );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &data[jLoc*thisLocalHeight],
                  &thisBuf[(rowShift+jLoc*colStrideUnion)*ldim], 
                  thisLocalHeight );
        }

        // Realign the input
        mpi::SendRecv 
        ( secondBuf, colStrideUnion*portionSize, sendColRankPart,
          firstBuf,  colStrideUnion*portionSize, recvColRankPart, 
          this->PartialColComm() );

        // Simultaneously Scatter in columns and Gather in rows
        mpi::AllToAll
        ( firstBuf,  portionSize, 
          secondBuf, portionSize, this->PartialUnionColComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colRank = recvColRankPart + k*colStridePart;
            const Int colShift = Shift_( colRank, colAlign, colStride );
            const Int colOffset = (colShift-colShiftA) / colStridePart;
            const Int localHeight = Length_( height, colShift, colStride );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
            {
                T* ACol = &ABuf[colOffset+jLoc*ALDim];
                const T* dataCol = &data[jLoc*localHeight];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    ACol[iLoc*colStrideUnion] = dataCol[iLoc];
            }
        }
    }
    A.auxMemory_.Release();
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::PartialRowAllToAll
( DistMatrix<T,UScat,VPart>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::PartialRowAllToAll");
        this->AssertSameGrid( A.Grid() );
    )
    const Int height = this->Height();
    const Int width = this->Width();
    A.AlignRowsAndResize
    ( this->RowAlign()%A.RowStride(), height, width, false, false );
    if( !A.Participating() )
        return;

    const Int colAlignA = A.ColAlign();
    const Int rowAlign = this->RowAlign();
    const Int rowAlignA = A.RowAlign();

    const Int rowStride = this->RowStride();
    const Int rowStridePart = this->PartialRowStride();
    const Int rowStrideUnion = this->PartialUnionRowStride();
    const Int rowRankPart = this->PartialRowRank();

    const Int rowShiftA = A.RowShift();

    const Int thisLocalWidth = this->LocalWidth();
    const Int localHeightA = A.LocalHeight();
    const Int maxLocalWidth = MaxLength(width,rowStride);
    const Int maxLocalHeight = MaxLength(height,rowStrideUnion);
    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

    const T* thisBuf = this->LockedBuffer();
    const Int ldim = this->LDim();
    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();

    T* buffer = A.auxMemory_.Require( 2*rowStrideUnion*portionSize );
    T* firstBuf  = &buffer[0];
    T* secondBuf = &buffer[rowStrideUnion*portionSize];

    if( rowAlignA == rowAlign % rowStridePart )
    {
        // Pack            
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            T* data = &firstBuf[k*portionSize];
            const Int colShift = Shift_( k, colAlignA, rowStrideUnion );
            const Int localHeight = Length_( height, colShift, rowStrideUnion );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                T* dataCol = &data[jLoc*localHeight];
                const T* thisCol = &thisBuf[colShift+jLoc*ldim];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    dataCol[iLoc] = thisCol[iLoc*rowStrideUnion];
            }
        }

        // Simultaneously Gather in rows and Scatter in columns
        mpi::AllToAll
        ( firstBuf,  portionSize, 
          secondBuf, portionSize, this->PartialUnionRowComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowRank = rowRankPart + k*rowStridePart;
            const Int rowShift = Shift_( rowRank, rowAlign, rowStride );
            const Int rowOffset = (rowShift-rowShiftA) / rowStridePart;
            const Int localWidth = Length_( width, rowShift, rowStride );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &ABuf[(rowOffset+jLoc*rowStrideUnion)*ALDim],
                  &data[jLoc*localHeightA], localHeightA );
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned PartialRowAllToAll" << std::endl;
#endif
        const Int rowAlignDiff = rowAlignA - (rowAlign%rowStridePart);
        const Int sendRowRankPart = 
            (rowRankPart+rowStridePart+rowAlignDiff) % rowStridePart;
        const Int recvRowRankPart =
            (rowRankPart+rowStridePart-rowAlignDiff) % rowStridePart;

        // Pack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            T* data = &secondBuf[k*portionSize];    
            const Int colShift = Shift_( k, colAlignA, rowStrideUnion );
            const Int localHeight = Length_( height, colShift, rowStrideUnion );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                T* dataCol = &data[jLoc*localHeight];
                const T* sourceCol = &thisBuf[colShift+jLoc*ldim];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    dataCol[iLoc] = sourceCol[iLoc*rowStrideUnion]; 
            }
        }

        // Realign the input
        mpi::SendRecv 
        ( secondBuf, rowStrideUnion*portionSize, sendRowRankPart,
          firstBuf,  rowStrideUnion*portionSize, recvRowRankPart, 
          this->PartialRowComm() );

        // Simultaneously Scatter in rows and Gather in columns
        mpi::AllToAll
        ( firstBuf,  portionSize, 
          secondBuf, portionSize, this->PartialUnionRowComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowRank = recvRowRankPart + k*rowStridePart;
            const Int rowShift = Shift_( rowRank, rowAlign, rowStride );
            const Int rowOffset = (rowShift-rowShiftA) / rowStridePart;
            const Int localWidth = Length_( width, rowShift, rowStride );
            EL_INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &ABuf[(rowOffset+jLoc*rowStrideUnion)*ALDim],
                  &data[jLoc*localHeightA], localHeightA );
        }
    }
    A.auxMemory_.Release();
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::RowSumScatterFrom( const DistMatrix<T,U,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GDM::RowSumScatterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    this->AlignColsAndResize
    ( A.ColAlign(), A.Height(), A.Width(), false, false );
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
    this->AlignRowsAndResize
    ( A.RowAlign(), A.Height(), A.Width(), false, false );
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
    this->AlignAndResize
    ( A.ColAlign(), A.RowAlign(), A.Height(), A.Width(), false, false );
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
    this->AlignAndResize
    ( A.ColAlign(), A.RowAlign(), A.Height(), A.Width(), false, false );
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
                EL_FMA_PARALLEL_FOR
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
            EL_OUTER_PARALLEL_FOR
            for( Int k=0; k<rowStride; ++k )
            {
                T* data = &buffer[k*portionSize];
                const Int thisRowShift = Shift_( k, rowAlign, rowStride );
                const Int thisLocalWidth = 
                    Length_(width,thisRowShift,rowStride);
                EL_INNER_PARALLEL_FOR
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
            EL_PARALLEL_FOR
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
#ifdef EL_UNALIGNED_WARNINGS
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
                EL_FMA_PARALLEL_FOR
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
            EL_OUTER_PARALLEL_FOR
            for( Int k=0; k<rowStride; ++k )
            {
                T* data = &secondBuf[k*recvSize_RS];
                const Int thisRowShift = Shift_( k, rowAlign, rowStride );
                const Int thisLocalWidth = 
                    Length_(width,thisRowShift,rowStride);
                EL_INNER_PARALLEL_FOR
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
            EL_FMA_PARALLEL_FOR
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
#ifdef EL_VECTOR_WARNINGS
    if( A.Width() == 1 && this->Grid().Rank() == 0 )
    {
        std::cerr <<
          "The vector version of ColSumScatterUpdate does not"
          " yet have a vector version implemented, but it would only "
          "require a modification of the vector version of RowSumScatterUpdate"
          << std::endl;
    }
#endif
#ifdef EL_CACHE_WARNINGS
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
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStride; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisColShift = Shift_( k, colAlign, colStride );
            const Int thisLocalHeight = Length_(height,thisColShift,colStride);
            EL_INNER_PARALLEL_FOR
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
        EL_FMA_PARALLEL_FOR
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
#ifdef EL_UNALIGNED_WARNINGS
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
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStride; ++k )
        {
            T* data = &secondBuf[k*recvSize_RS];
            const Int thisColShift = Shift_( k, colAlign, colStride );
            const Int thisLocalHeight = Length_(height,thisColShift,colStride);
            EL_INNER_PARALLEL_FOR
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
        EL_FMA_PARALLEL_FOR
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
    EL_OUTER_PARALLEL_FOR
    for( Int l=0; l<rowStride; ++l )
    {
        const Int thisRowShift = Shift_( l, rowAlign, rowStride );
        const Int thisLocalWidth = Length_( width, thisRowShift, rowStride );
        for( Int k=0; k<colStride; ++k )
        {
            T* data = &buffer[(k+l*colStride)*recvSize];
            const Int thisColShift = Shift_( k, colAlign, colStride );
            const Int thisLocalHeight = Length_(height,thisColShift,colStride);
            EL_INNER_PARALLEL_FOR
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
    EL_FMA_PARALLEL_FOR
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
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = rowRankPart+k*rowStridePart;
            const Int thisRowShift = Shift_( thisRank, rowAlign, rowStride );
            const Int thisRowOffset = 
                (thisRowShift-rowShiftOfA) / rowStridePart;
            const Int thisLocalWidth = 
                Length_( width, thisRowShift, rowStride );
            EL_INNER_PARALLEL_FOR
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
        EL_PARALLEL_FOR
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

#ifdef EL_CACHE_WARNINGS
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
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = colRankPart+k*colStridePart;
            const Int thisColShift = Shift_( thisRank, colAlign, colStride );
            const Int thisColOffset = 
                (thisColShift-colShiftOfA) / colStridePart;
            const Int thisLocalHeight = 
                Length_( height, thisColShift, colStride );
            EL_INNER_PARALLEL_FOR
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
        EL_PARALLEL_FOR
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

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::TransposeColAllGather
( DistMatrix<T,V,UGath>& A, bool conjugate ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::TransposeColAllGather"))
    DistMatrix<T,V,U> ATrans( this->Grid() );
    ATrans.AlignWith( *this );
    ATrans.Resize( this->Width(), this->Height() );
    Transpose( this->LockedMatrix(), ATrans.Matrix(), conjugate );
    ATrans.RowAllGather( A );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::TransposePartialColAllGather
( DistMatrix<T,V,UPart>& A, bool conjugate ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::TransposePartialColAllGather"))
    DistMatrix<T,V,U> ATrans( this->Grid() );
    ATrans.AlignWith( *this );
    ATrans.Resize( this->Width(), this->Height() );
    Transpose( this->LockedMatrix(), ATrans.Matrix(), conjugate );
    ATrans.PartialRowAllGather( A );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AdjointColAllGather( DistMatrix<T,V,UGath>& A ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AdjointRowAllGather"))
    this->TransposeColAllGather( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AdjointPartialColAllGather
( DistMatrix<T,V,UPart>& A ) const
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AdjointPartialColAllGather"))
    this->TransposePartialColAllGather( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::TransposeColFilterFrom
( const DistMatrix<T,V,UGath>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::TransposeColFilterFrom"))
    DistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    AFilt.RowFilterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::TransposeRowFilterFrom
( const DistMatrix<T,VGath,U>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::TransposeRowFilterFrom"))
    DistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    AFilt.ColFilterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::TransposePartialColFilterFrom
( const DistMatrix<T,V,UPart>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::TransposePartialColFilterFrom"))
    DistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    AFilt.PartialRowFilterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::TransposePartialRowFilterFrom
( const DistMatrix<T,VPart,U>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::TransposePartialRowFilterFrom"))
    DistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    AFilt.PartialColFilterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AdjointColFilterFrom( const DistMatrix<T,V,UGath>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AdjointColFilterFrom"))
    this->TransposeColFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AdjointRowFilterFrom( const DistMatrix<T,VGath,U>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AdjointRowFilterFrom"))
    this->TransposeRowFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AdjointPartialColFilterFrom
( const DistMatrix<T,V,UPart>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AdjointPartialColFilterFrom"))
    this->TransposePartialColFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AdjointPartialRowFilterFrom
( const DistMatrix<T,VPart,U>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AdjointPartialRowFilterFrom"))
    this->TransposePartialRowFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::TransposeColSumScatterFrom
( const DistMatrix<T,V,UGath>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::TransposeColSumScatterFrom"))
    DistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.RowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( ASumFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::TransposePartialColSumScatterFrom
( const DistMatrix<T,V,UPart>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::TransposePartialColSumScatterFrom"))
    DistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.PartialRowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( ASumFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AdjointColSumScatterFrom
( const DistMatrix<T,V,UGath>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AdjointColSumScatterFrom"))
    this->TransposeColSumScatterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AdjointPartialColSumScatterFrom
( const DistMatrix<T,V,UPart>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AdjointPartialColSumScatterFrom"))
    this->TransposePartialColSumScatterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::TransposeColSumScatterUpdate
( T alpha, const DistMatrix<T,V,UGath>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::TransposeColSumScatterUpdate"))
    DistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.RowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    // ALoc += alpha ASumFiltLoc'
    El::Matrix<T>& ALoc = this->Matrix();
    const El::Matrix<T>& BLoc = ASumFilt.LockedMatrix();
    const Int localHeight = ALoc.Height();
    const Int localWidth = ALoc.Width();
    if( conjugate )
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*Conj(BLoc.Get(jLoc,iLoc)) );
    }
    else
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*BLoc.Get(jLoc,iLoc) );
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::TransposePartialColSumScatterUpdate
( T alpha, const DistMatrix<T,V,UPart>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::TransposePartialColSumScatterUpdate"))
    DistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.PartialRowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    // ALoc += alpha ASumFiltLoc'
    El::Matrix<T>& ALoc = this->Matrix();
    const El::Matrix<T>& BLoc = ASumFilt.LockedMatrix();
    const Int localHeight = ALoc.Height();
    const Int localWidth = ALoc.Width();
    if( conjugate )
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*Conj(BLoc.Get(jLoc,iLoc)) );
    }
    else
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*BLoc.Get(jLoc,iLoc) );
    }
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AdjointColSumScatterUpdate
( T alpha, const DistMatrix<T,V,UGath>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AdjointColSumScatterUpdate"))
    this->TransposeColSumScatterUpdate( alpha, A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::AdjointPartialColSumScatterUpdate
( T alpha, const DistMatrix<T,V,UPart>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::AdjointPartialColSumScatterUpdate"))
    this->TransposePartialColSumScatterUpdate( alpha, A, true );
}

// Diagonal manipulation
// =====================
template<typename T,Dist U,Dist V>
bool
GeneralDistMatrix<T,U,V>::DiagonalAlignedWith
( const El::DistData& d, Int offset ) const
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
    const El::Grid& grid = this->Grid();

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
    const El::Grid& grid = this->Grid();

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
      []( T& alpha, Base<T> beta ) { El::SetRealPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralDistMatrix<T,U,V>::SetImagPartOfDiagonal
( const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GDM::SetImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { El::SetImagPart(alpha,beta); } );
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
      { El::UpdateRealPart(alpha,gamma*beta); } );
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
      { El::UpdateImagPart(alpha,gamma*beta); } );
}

// Private section
// ###############

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
    d.AlignCols( this->DiagonalAlign(offset), false );
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

    EL_PARALLEL_FOR
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

    EL_PARALLEL_FOR
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

PROTO(Int);
#ifndef EL_DISABLE_FLOAT
PROTO(float);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef EL_DISABLE_COMPLEX
#endif // ifndef EL_DISABLE_FLOAT

PROTO(double);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El
