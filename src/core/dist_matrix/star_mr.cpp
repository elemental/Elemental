/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

template<typename T>
DistMatrix<T,STAR,MR>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->SetShifts(); }

template<typename T>
DistMatrix<T,STAR,MR>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T>
DistMatrix<T,STAR,MR>::DistMatrix
( Int height, Int width, Int rowAlignment, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Align(0,rowAlignment); this->ResizeTo(height,width); }

template<typename T>
DistMatrix<T,STAR,MR>::DistMatrix
( Int height, Int width, Int rowAlignment, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Align(0,rowAlignment); this->ResizeTo(height,width,ldim); }

template<typename T>
DistMatrix<T,STAR,MR>::DistMatrix
( Int height, Int width, Int rowAlignment, const T* buffer, Int ldim, 
  const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->LockedAttach(height,width,rowAlignment,buffer,ldim,g); }

template<typename T>
DistMatrix<T,STAR,MR>::DistMatrix
( Int height, Int width, Int rowAlignment, T* buffer, Int ldim, 
  const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Attach(height,width,rowAlignment,buffer,ldim,g); }

template<typename T>
DistMatrix<T,STAR,MR>::DistMatrix( const DistMatrix<T,STAR,MR>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,MR]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [* ,MR] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DistMatrix<T,STAR,MR>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,MR]::DistMatrix");
#endif
    this->SetShifts();
    if( STAR != U || MR != V || 
        reinterpret_cast<const DistMatrix<T,STAR,MR>*>(&A) != this )   
        *this = A;
    else
        LogicError("Tried to construct [* ,MR] with itself");
}

template<typename T>
DistMatrix<T,STAR,MR>::~DistMatrix()
{ }

template<typename T>
elem::DistData
DistMatrix<T,STAR,MR>::DistData() const
{
    elem::DistData data;
    data.colDist = STAR;
    data.rowDist = MR;
    data.colAlignment = 0;
    data.rowAlignment = this->rowAlignment_;
    data.root = 0;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T>
Int
DistMatrix<T,STAR,MR>::ColStride() const
{ return 1; }

template<typename T>
Int
DistMatrix<T,STAR,MR>::RowStride() const
{ return this->grid_->Width(); }

template<typename T>
Int
DistMatrix<T,STAR,MR>::ColRank() const
{ return 0; }

template<typename T>
Int
DistMatrix<T,STAR,MR>::RowRank() const
{ return this->grid_->Col(); }

template<typename T>
void
DistMatrix<T,STAR,MR>::AlignWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::AlignWith");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( data.colDist == MR )
        this->rowAlignment_ = data.colAlignment;
    else if( data.rowDist == MR )
        this->rowAlignment_ = data.rowAlignment;
    else if( data.colDist == VR )
        this->rowAlignment_ = data.colAlignment % this->RowStride();
    else if( data.rowDist == VR )
        this->rowAlignment_ = data.rowAlignment % this->RowStride();
#ifndef RELEASE
    else LogicError("Nonsensical alignment");
#endif
    this->constrainedRowAlignment_ = true;
    this->SetShifts();
}

template<typename T>
void
DistMatrix<T,STAR,MR>::AlignWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,STAR,MR>::AlignRowsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

template<typename T>
void
DistMatrix<T,STAR,MR>::AlignRowsWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,STAR,MR>::Attach
( Int height, Int width, Int rowAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::Attach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->viewType_ = VIEW;
    this->SetRowShift();
    if( this->Participating() )
    {
        const Int localWidth = Length(width,this->rowShift_,g.Width());
        this->matrix_.Attach_( height, localWidth, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,STAR,MR>::LockedAttach
( Int height, Int width, Int rowAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->viewType_ = LOCKED_VIEW;
    this->SetRowShift();
    if( this->Participating() )
    {
        const Int localWidth = Length(width,this->rowShift_,g.Width());
        this->matrix_.LockedAttach_( height, localWidth, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,STAR,MR>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( height, Length(width,this->RowShift(),this->Grid().Width()) );
}

template<typename T>
void
DistMatrix<T,STAR,MR>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( height, Length(width,this->RowShift(),this->Grid().Width()), ldim );
}

template<typename T>
T
DistMatrix<T,STAR,MR>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::Get");
    this->AssertValidEntry( i, j );
    if( !this->Participating() )
        LogicError("Should only be called by members of grid");
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // column within each process row
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();

    T u;
    if( g.Col() == ownerCol )
    {
        const Int jLoc = (j-this->RowShift()) / g.Width();
        u = this->GetLocal(i,jLoc);
    }
    mpi::Broadcast( u, ownerCol, g.RowComm() );
    return u;
}

template<typename T>
void
DistMatrix<T,STAR,MR>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.Col() == ownerCol )
    {
        const Int jLoc = (j-this->RowShift()) / g.Width();
        this->SetLocal(i,jLoc,u);
    }
}

template<typename T>
void
DistMatrix<T,STAR,MR>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.Col() == ownerCol )
    {
        const Int jLoc = (j-this->RowShift()) / g.Width();
        this->UpdateLocal(i,jLoc,u);
    }
}

//
// Utility functions, e.g., SumOverCol
//

template<typename T> 
void
DistMatrix<T,STAR,MR>::SumOverCol()
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::SumOverCol");
    this->AssertNotLocked();
#endif
    const Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int localSize = mpi::Pad( localHeight*localWidth );

    T* buffer = this->auxMemory_.Require( 2*localSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[localSize];

    // Pack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* thisCol = &thisBuf[jLoc*thisLDim];
        T* sendBufCol = &sendBuf[jLoc*localHeight];
        MemCopy( sendBufCol, thisCol, localHeight );
    }

    // AllReduce col
    mpi::AllReduce( sendBuf, recvBuf, localSize, g.ColComm() );

    // Unpack
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* recvBufCol = &recvBuf[jLoc*localHeight];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        MemCopy( thisCol, recvBufCol, localHeight );
    }
    this->auxMemory_.Release();
}

template<typename T>
void
DistMatrix<T,STAR,MR>::AdjointFrom( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::AdjointFrom");
#endif
    this->TransposeFrom( A, true );
}

template<typename T>
void
DistMatrix<T,STAR,MR>::TransposeFrom
( const DistMatrix<T,VR,STAR>& A, bool conjugate )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::TransposeFrom");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetRowAlignmentAndResize
    ( A.ColAlignment()%g.Width(), A.Width(), A.Height() );
    if( !this->Participating() )
        return;

    if( this->RowAlignment() == A.ColAlignment() % g.Width() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();

        const Int width = this->Width();
        const Int height = this->Height();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalHeightOfA );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        if( conjugate )
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &sendBuf[jLoc*height];
                const T* sourceCol = &ABuf[jLoc];
                for( Int i=0; i<height; ++i )
                    destCol[i] = Conj( sourceCol[i*ALDim] );
            }
        }
        else
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &sendBuf[jLoc*height];
                const T* sourceCol = &ABuf[jLoc];
                for( Int i=0; i<height; ++i )
                    destCol[i] = sourceCol[i*ALDim];
            }
        }

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int rowShift = this->RowShift();
        const Int colAlignmentOfA = A.ColAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int colShiftOfA = Shift_( col+k*c, colAlignmentOfA, p );
            const Int rowOffset = (colShiftOfA-rowShift) / c;
            const Int localWidth = Length_( width, colShiftOfA, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc ) 
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MR].TransposeFrom[VR,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();
        const Int rank = g.VRRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowShift = this->RowShift();

        const Int sendRank = (rank+p+rowAlignment-colAlignmentOfA) % p;
        const Int recvRank = (rank+p+colAlignmentOfA-rowAlignment) % p;

        const Int width = this->Width();
        const Int height = this->Height();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalHeightOfA );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        if( conjugate )
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &secondBuf[jLoc*height];
                const T* sourceCol = &ABuf[jLoc];
                for( Int i=0; i<height; ++i )
                    destCol[i] = Conj( sourceCol[i*ALDim] );
            }
        }
        else
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &secondBuf[jLoc*height];
                const T* sourceCol = &ABuf[jLoc];
                for( Int i=0; i<height; ++i )
                    destCol[i] = sourceCol[i*ALDim];
            }
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuf, portionSize, sendRank, 
          firstBuf,  portionSize, recvRank, g.VRComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShiftOfA = Shift_( col+c*k, rowAlignment, p );
            const Int rowOffset = (colShiftOfA-rowShift) / c;
            const Int localWidth = Length_( width, colShiftOfA, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "The matrix redistribution [* ,MR] <- [MC,MR] potentially causes a "
          "large amount of cache-thrashing. If possible, avoid it by "
          "performing the redistribution with a (conjugate)-transpose: \n"
          << "  [MR,* ].(Conjugate)TransposeFrom([MC,MR])" << std::endl;
    }
#endif
    this->SetRowAlignmentAndResize( A.RowAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlignment() == A.RowAlignment() )
    {
        if( A.Height() == 1 )
        {
            const Int localWidth = this->LocalWidth();
            T* bcastBuf = this->auxMemory_.Require( localWidth );

            if( g.Row() == A.ColAlignment() )
            {
                this->matrix_ = A.LockedMatrix();

                // Pack
                const T* thisBuf = this->LockedBuffer();
                const Int thisLDim = this->LDim();
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                    bcastBuf[jLoc] = thisBuf[jLoc*thisLDim];
            }

            // Communicate
            mpi::Broadcast
            ( bcastBuf, localWidth, A.ColAlignment(), g.ColComm() );

            // Unpack
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                thisBuf[jLoc*thisLDim] = bcastBuf[jLoc];

            this->auxMemory_.Release();
        }
        else
        {
            const Int r = g.Height();
            const Int height = this->Height();
            const Int localWidth = this->LocalWidth();
            const Int localHeightOfA = A.LocalHeight();
            const Int maxLocalHeight = MaxLength(height,r);
            const Int portionSize = mpi::Pad( maxLocalHeight*localWidth );

            T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[portionSize];

            // Pack 
            const Int ALDim = A.LDim();
            const T* ABuf = A.LockedBuffer();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* ACol = &ABuf[jLoc*ALDim];
                T* sendBufCol = &sendBuf[jLoc*localHeightOfA];
                MemCopy( sendBufCol, ACol, localHeightOfA );
            }

            // Communicate
            mpi::AllGather
            ( sendBuf, portionSize,
              recvBuf, portionSize, g.ColComm() );

            // Unpack
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
            const Int colAlignmentOfA = A.ColAlignment();
            OUTER_PARALLEL_FOR
            for( Int k=0; k<r; ++k )
            {
                const T* data = &recvBuf[k*portionSize];
                const Int colShift = Shift_( k, colAlignmentOfA, r );
                const Int localHeight = Length_( height, colShift, r );
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                {
                    T* destCol = &thisBuf[colShift+jLoc*thisLDim];
                    const T* sourceCol = &data[jLoc*localHeight];
                    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                        destCol[iLoc*r] = sourceCol[iLoc];
                }
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MR] <- [MC,MR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int col = g.Col();

        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const Int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

        if( A.Height() == 1 )
        {
            const Int localWidth = this->LocalWidth();
            T* bcastBuf;

            if( g.Row() == A.ColAlignment() )
            {
                const Int localWidthOfA = A.LocalWidth();

                T* buffer = this->auxMemory_.Require(localWidth+localWidthOfA);
                T* sendBuf = &buffer[0];
                bcastBuf   = &buffer[localWidthOfA];

                // Pack
                const T* ABuf = A.LockedBuffer();
                const Int ALDim = A.LDim();
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
                    sendBuf[jLoc] = ABuf[jLoc*ALDim];

                // Communicate
                mpi::SendRecv
                ( sendBuf,  localWidthOfA, sendCol, 
                  bcastBuf, localWidth,    recvCol, g.RowComm() );
            }
            else
            {
                bcastBuf = this->auxMemory_.Require( localWidth );
            }

            // Communicate
            mpi::Broadcast
            ( bcastBuf, localWidth, A.ColAlignment(), g.ColComm() );

            // Unpack
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                thisBuf[jLoc*thisLDim] = bcastBuf[jLoc];
            this->auxMemory_.Release();
        }
        else
        {
            const Int height = this->Height();
            const Int localWidth  = this->LocalWidth();
            const Int localHeightOfA = A.LocalHeight();
            const Int localWidthOfA  = A.LocalWidth();
            const Int maxLocalHeight = MaxLength(height,r);
            const Int portionSize = mpi::Pad( maxLocalHeight*localWidth );

            T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
            T* firstBuf = &buffer[0];
            T* secondBuf = &buffer[portionSize];

            // Pack
            const Int ALDim = A.LDim();
            const T* ABuf = A.LockedBuffer();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                const T* ACol = &ABuf[jLoc*ALDim];
                T* secondBufCol = &secondBuf[jLoc*localHeightOfA];
                MemCopy( secondBufCol, ACol, localHeightOfA );
            }

            // Perform the SendRecv: puts the new data into the first buffer
            mpi::SendRecv
            ( secondBuf, portionSize, sendCol,
              firstBuf,  portionSize, recvCol, g.RowComm() );

            // Use the output of the SendRecv as input to the AllGather
            mpi::AllGather
            ( firstBuf,  portionSize,
              secondBuf, portionSize, g.ColComm() );

            // Unpack the contents of each member of the process col
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
            const Int colAlignmentOfA = A.ColAlignment();
            OUTER_PARALLEL_FOR
            for( Int k=0; k<r; ++k )
            {
                const T* data = &secondBuf[k*portionSize];
                const Int colShift = Shift_( k, colAlignmentOfA, r );
                const Int localHeight = Length_( height, colShift, r );
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                {
                    T* destCol = &thisBuf[colShift+jLoc*thisLDim];
                    const T* sourceCol = &data[jLoc*localHeight];
                    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                        destCol[iLoc*r] = sourceCol[iLoc];
                }
            }
            this->auxMemory_.Release();
        }
    }
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(false,true,0,this->RowAlignment(),g);

    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetRowAlignmentAndResize( A.RowAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlignment() == A.RowAlignment() )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MR] <- [* ,MR]." << std::endl;
#endif
        const Int rank = g.Col();
        const Int c = g.Width();

        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRank = (rank+c+rowAlignment-rowAlignmentOfA) % c;
        const Int recvRank = (rank+c+rowAlignmentOfA-rowAlignment) % c;

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfA = A.LocalWidth();

        const Int sendSize = height * localWidthOfA;
        const Int recvSize = height * localWidth;

        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuf[jLoc*ALDim];
            T* sendBufCol = &sendBuf[jLoc*height];
            MemCopy( sendBufCol, ACol, height );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank,
          recvBuf, recvSize, recvRank, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* recvBufCol = &recvBuf[jLoc*height];
            T* thisCol = &thisBuf[jLoc*thisLDim];
            MemCopy( thisCol, recvBufCol, height );
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [MD,* ]");
#endif
    LogicError("[* ,MR] = [MD,* ] not yet implemented");
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,MD>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [* ,MD]");
#endif
    LogicError("[* ,MR] = [* ,MD] not yet implemented");
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [MR,MC]");
#endif
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(g) );
    *A_STAR_VC = A;

    std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(true,this->RowAlignment(),g) );
    *A_STAR_VR = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    *this = *A_STAR_VR;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [MR,* ]");
#endif
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(g) );
    *A_VR_STAR = A;

    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR
    ( new DistMatrix<T,MC,MR>(false,true,0,this->RowAlignment(),g) );
    *A_MC_MR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_MC_MR;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [* ,MC]");
#endif
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(g) );
    *A_STAR_VC = A;

    std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(true,this->RowAlignment(),g) );
    *A_STAR_VR = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_MC_MR;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [VC,* ]");
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(false,true,0,this->RowAlignment(),g);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [* ,VC]");
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR> A_STAR_VR(true,this->RowAlignment(),g);
    A_STAR_VR = A;
    *this = A_STAR_VR;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [VR,* ]");
#endif
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = A;

    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR
    ( new DistMatrix<T,MC,MR>(false,true,0,this->RowAlignment(),g) );
    *A_MC_MR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_MC_MR;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetRowAlignmentAndResize
    ( A.RowAlignment()%g.Width(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlignment() == A.RowAlignment() % g.Width() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();

        const Int width = this->Width();
        const Int height = this->Height();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidthOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalWidthOfA );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuf[jLoc*ALDim];
            T* sendBufCol = &sendBuf[jLoc*height];
            MemCopy( sendBufCol, ACol, height );
        }

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int rowShift = this->RowShift();
        const Int rowAlignmentOfA = A.RowAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int rowShiftOfA = Shift_( col+k*c, rowAlignmentOfA, p );
            const Int rowOffset = (rowShiftOfA-rowShift) / c;
            const Int localWidth = Length_( width, rowShiftOfA, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MR] <- [* ,VR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();
        const Int rank = g.VRRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int rowShift = this->RowShift();

        const Int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const Int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const Int width = this->Width();
        const Int height = this->Height();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidthOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalWidthOfA );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuf[jLoc*ALDim];
            T* secondBufCol = &secondBuf[jLoc*height];
            MemCopy( secondBufCol, ACol, height );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuf, portionSize, sendRank,
          firstBuf,  portionSize, recvRank, g.VRComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShiftOfA = Shift_( col+c*k, rowAlignment, p );
            const Int rowOffset = (rowShiftOfA-rowShift) / c;
            const Int localWidth = Length_( width, rowShiftOfA, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int c = this->Grid().Width();
    const Int rowShift = this->RowShift();

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();

    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* ACol = &ABuf[(rowShift+jLoc*c)*ALDim];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        MemCopy( thisCol, ACol, localHeight );
    }
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR] = [o ,o ]");
#endif
    DistMatrix<T,MC,MR> A_MC_MR( A );
    A_MC_MR.AlignWith( *this );
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

//
// Routines which explicitly work in the complex plane
//

template<typename T>
void
DistMatrix<T,STAR,MR>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();
    if( g.Col() == ownerCol )
    {
        const Int jLoc = (j-this->RowShift()) / g.Width();
        this->SetLocalRealPart( i, jLoc, u );
    }
}

template<typename T>
void
DistMatrix<T,STAR,MR>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();
    if( g.Col() == ownerCol )
    {
        const Int jLoc = (j-this->RowShift()) / g.Width();
        this->SetLocalImagPart( i, jLoc, u );
    }
}

template<typename T>
void
DistMatrix<T,STAR,MR>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();
    if( g.Col() == ownerCol )
    {
        const Int jLoc = (j-this->RowShift()) / g.Width();
        this->UpdateLocalRealPart( i, jLoc, u );
    }
}

template<typename T>
void
DistMatrix<T,STAR,MR>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();
    if( g.Col() == ownerCol )
    {
        const Int jLoc = (j-this->RowShift()) / g.Width();
        this->UpdateLocalImagPart( i, jLoc, u );
    }
}

#define PROTO(T) template class DistMatrix<T,STAR,MR>
#define COPY(T,CD,RD) \
  template DistMatrix<T,STAR,MR>::DistMatrix( const DistMatrix<T,CD,RD>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  MR  ); \
  COPY(T,MC,  STAR); \
  COPY(T,MD,  STAR); \
  COPY(T,MR,  MC  ); \
  COPY(T,MR,  STAR); \
  COPY(T,STAR,MC  ); \
  COPY(T,STAR,MD  ); \
  COPY(T,STAR,STAR); \
  COPY(T,STAR,VC  ); \
  COPY(T,STAR,VR  ); \
  COPY(T,VC,  STAR); \
  COPY(T,VR,  STAR);

FULL(Int);
#ifndef DISABLE_FLOAT
FULL(float);
#endif
FULL(double);

#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
FULL(Complex<float>);
#endif
FULL(Complex<double>);
#endif

} // namespace elem
