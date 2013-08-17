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
DistMatrix<T,STAR,VR>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->SetShifts(); }

template<typename T>
DistMatrix<T,STAR,VR>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T>
DistMatrix<T,STAR,VR>::DistMatrix
( Int height, Int width, Int rowAlignment, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Align(0,rowAlignment); this->ResizeTo(height,width); }

template<typename T>
DistMatrix<T,STAR,VR>::DistMatrix
( Int height, Int width, Int rowAlignment, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Align(0,rowAlignment); this->ResizeTo(height,width,ldim); }

template<typename T>
DistMatrix<T,STAR,VR>::DistMatrix
( Int height, Int width, Int rowAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->LockedAttach(height,width,rowAlignment,buffer,ldim,g); }

template<typename T>
DistMatrix<T,STAR,VR>::DistMatrix
( Int height, Int width, Int rowAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Attach(height,width,rowAlignment,buffer,ldim,g); }

template<typename T>
DistMatrix<T,STAR,VR>::DistMatrix( const DistMatrix<T,STAR,VR>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,VR]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [* ,VR] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DistMatrix<T,STAR,VR>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,VR]::DistMatrix");
#endif
    this->SetShifts();
    if( STAR != U || VR != V || 
        reinterpret_cast<const DistMatrix<T,STAR,VR>*>(&A) != this ) 
        *this = A;
    else
        LogicError("Tried to construct [* ,VR] with itself");
}

template<typename T>
DistMatrix<T,STAR,VR>::~DistMatrix()
{ }

template<typename T>
elem::DistData
DistMatrix<T,STAR,VR>::DistData() const
{
    elem::DistData data;
    data.colDist = STAR;
    data.rowDist = VR;
    data.colAlignment = 0;
    data.rowAlignment = this->rowAlignment_;
    data.root = 0;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T>
Int
DistMatrix<T,STAR,VR>::ColStride() const
{ return 1; }

template<typename T>
Int
DistMatrix<T,STAR,VR>::RowStride() const
{ return this->grid_->Size(); }

template<typename T>
Int
DistMatrix<T,STAR,VR>::ColRank() const
{ return 0; }

template<typename T>
Int
DistMatrix<T,STAR,VR>::RowRank() const
{ return this->grid_->VRRank(); }

template<typename T>
void
DistMatrix<T,STAR,VR>::AlignWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::AlignWith");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );
    
    if( data.colDist == MR || data.colDist == VR )
        this->rowAlignment_ = data.colAlignment;
    else if( data.rowDist == MR || data.rowDist == VR )
        this->rowAlignment_ = data.rowAlignment;
#ifndef RELEASE
    else LogicError("Nonsensical alignment");
#endif
    this->constrainedRowAlignment_ = true;
    this->SetShifts();
}

template<typename T>
void
DistMatrix<T,STAR,VR>::AlignWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,STAR,VR>::AlignRowsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

template<typename T>
void
DistMatrix<T,STAR,VR>::AlignRowsWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,STAR,VR>::Attach
( Int height, Int width, Int rowAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::Attach");
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
        const Int localWidth = Length(width,this->rowShift_,g.Size());
        this->matrix_.Attach_( height, localWidth, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,STAR,VR>::LockedAttach
( Int height, Int width, Int rowAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::LockedAttach");
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
        const Int localWidth = Length(width,this->rowShift_,g.Size());
        this->matrix_.LockedAttach_( height, localWidth, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,STAR,VR>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    const elem::Grid& g = this->Grid();
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( height, Length(width,this->RowShift(),g.Size()) );
}

template<typename T>
void
DistMatrix<T,STAR,VR>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    const elem::Grid& g = this->Grid();
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( height, Length(width,this->RowShift(),g.Size()), ldim );
}

template<typename T>
T
DistMatrix<T,STAR,VR>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::Get");
    this->AssertValidEntry( i, j );
    // TODO: Generalize this function to always work...
    if( !this->Participating() )
        LogicError("Should only call with processes in grid");
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    T u;
    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocal(i,jLoc);
    }
    mpi::Broadcast( u, ownerRank, g.VRComm() );
    return u;
}

template<typename T>
void
DistMatrix<T,STAR,VR>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->SetLocal(i,jLoc,u);
    }
}

template<typename T>
void
DistMatrix<T,STAR,VR>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->UpdateLocal(i,jLoc,u);
    }
}

//
// Utility functions, e.g., AdjointFrom
//

template<typename T>
void
DistMatrix<T,STAR,VR>::AdjointFrom( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[*, VR]::AdjointFrom");
#endif
    this->TransposeFrom( A, true );
}

template<typename T>
void
DistMatrix<T,STAR,VR>::TransposeFrom
( const DistMatrix<T,MR,STAR>& A, bool conjugate )
{ 
#ifndef RELEASE
    CallStackEntry entry("[*, VR]::TransposeFrom");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetRowAlignmentAndResize( A.ColAlignment(), A.Width(), A.Height() );
    if( !this->Participating() )
        return;

    if( this->RowAlignment() % g.Width() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int rowShift = this->RowShift();
        const Int colShiftOfA = A.ColShift();
        const Int rowOffset = (rowShift-colShiftOfA) / c;

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();

        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        if( conjugate )
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[jLoc*thisLDim];
                const T* sourceCol = &ABuf[rowOffset+jLoc*r];
                for( Int i=0; i<height; ++i )
                    destCol[i] = Conj( sourceCol[i*ALDim] );
            }
        }
        else
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[jLoc*thisLDim];
                const T* sourceCol = &ABuf[rowOffset+jLoc*r];
                for( Int i=0; i<height; ++i )
                    destCol[i] = sourceCol[i*ALDim];
            }
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VR]::AdjointFrom" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int col = g.Col();
        const Int colShiftOfA = A.ColShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const Int sendCol = (col+c+(rowAlignment%c)-colAlignmentOfA) % c;
        const Int recvCol = (col+c+colAlignmentOfA-(rowAlignment%c)) % c;
        const Int sendRank = sendCol + c*row;

        const Int sendRowShift = Shift( sendRank, rowAlignment, p );
        const Int sendRowOffset = (sendRowShift-colShiftOfA) / c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfSend = Length(width,sendRowShift,p);

        const Int sendSize = height * localWidthOfSend;
        const Int recvSize = height * localWidth;

        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        if( conjugate )
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthOfSend; ++jLoc )
            {
                T* destCol = &sendBuf[jLoc*height];
                const T* sourceCol = &ABuf[sendRowOffset+jLoc*r];
                for( Int i=0; i<height; ++i )
                    destCol[i] = Conj( sourceCol[i*ALDim] );
            }
        }
        else
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthOfSend; ++jLoc )
            {
                T* destCol = &sendBuf[jLoc*height];
                const T* sourceCol = &ABuf[sendRowOffset+jLoc*r];
                for( Int i=0; i<height; ++i )
                    destCol[i] = sourceCol[i*ALDim];
            }
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendCol, 
          recvBuf, recvSize, recvCol, g.RowComm() );

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
}

template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetRowAlignmentAndResize( A.RowAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int maxHeight = MaxLength(height,r);
        const Int maxWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*r*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &sendBuf[k*portionSize];
            const Int thisRank = col+k*c;
            const Int thisRowShift = Shift_(thisRank,rowAlignment,p);
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* ACol = &ABuf[(thisRowOffset+jLoc*r)*ALDim];
                T* dataCol = &data[jLoc*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuf, portionSize, recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisColShift = Shift_(k,colAlignmentOfA,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[thisColShift+jLoc*thisLDim];
                const T* sourceCol = &data[jLoc*thisLocalHeight];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc*r] = sourceCol[iLoc];
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VR] <- [MC,MR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendCol = (col+c+(rowAlignment%c)-rowAlignmentOfA) % c;
        const Int recvCol = (col+c+rowAlignmentOfA-(rowAlignment%c)) % c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int maxHeight = MaxLength(height,r);
        const Int maxWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*r*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[r*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int thisRank = sendCol+k*c;
            const Int thisRowShift = Shift_(thisRank,rowAlignment,p);
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* ACol = &ABuf[(thisRowOffset+jLoc*r)*ALDim];
                T* dataCol = &data[jLoc*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // AllToAll to gather all of the unaligned [*,VR] data into firstBuf
        mpi::AllToAll
        ( secondBuf, portionSize, firstBuf, portionSize, g.ColComm() );

        // SendRecv: properly align the [*,VR] via a trade in the column
        mpi::SendRecv
        ( firstBuf,  portionSize, sendCol, 
          secondBuf, portionSize, recvCol, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int thisColShift = Shift_(k,colAlignmentOfA,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[thisColShift+jLoc*thisLDim];
                const T* sourceCol = &data[jLoc*thisLocalHeight];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc*r] = sourceCol[iLoc];
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [MC,* ]");
#endif
    DistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetRowAlignmentAndResize( A.RowAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int rowShift = this->RowShift();
        const Int rowShiftOfA = A.RowShift();
        const Int rowOffset = (rowShift-rowShiftOfA) / c;

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();

        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* ACol = &ABuf[(rowOffset+jLoc*r)*ALDim];
            T* thisCol = &thisBuf[jLoc*thisLDim];
            MemCopy( thisCol, ACol, height );
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VR] <- [* ,MR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int col = g.Col();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const Int sendCol = (col+c+(rowAlignment%c)-rowAlignmentOfA) % c;
        const Int recvCol = (col+c+rowAlignmentOfA-(rowAlignment%c)) % c;
        const Int sendRank = sendCol + c*row;

        const Int sendRowShift = Shift( sendRank, rowAlignment, p );
        const Int sendRowOffset = (sendRowShift-rowShiftOfA) / c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfSend = Length(width,sendRowShift,p);

        const Int sendSize = height * localWidthOfSend;
        const Int recvSize = height * localWidth;

        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidthOfSend; ++jLoc )
        {
            const T* ACol = &ABuf[(sendRowOffset+jLoc*r)*ALDim];
            T* sendBufCol = &sendBuf[jLoc*height];
            MemCopy( sendBufCol, ACol, height );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendCol, 
          recvBuf, recvSize, recvCol, g.RowComm() );

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
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [MD,* ]");
#endif
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [* ,MD]");
#endif
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [MR,MC]");
#endif
    DistMatrix<T,STAR,VC> A_STAR_VC( A );
    *this = A_STAR_VC;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [MR,* ]");
#endif
    std::unique_ptr<DistMatrix<T,MR,MC>> A_MR_MC( new DistMatrix<T,MR,MC>(A) );
    std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(*A_MR_MC) );
    delete A_MR_MC.release(); // lowers memory highwater
    *this = *A_STAR_VC;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [* ,MC]");
#endif
    DistMatrix<T,STAR,VC> A_STAR_VC( A );
    *this = A_STAR_VC;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [VC,* ]");
#endif
    DistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [* ,VC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;
    
    const Int height = this->Height();
    const Int localWidth = this->LocalWidth();
    const Int localWidthOfA = A.LocalWidth();

    const Int sendSize = height * localWidthOfA;
    const Int recvSize = height * localWidth;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int rankCM = g.VCRank();
    const Int rankRM = g.VRRank();

    const Int rowShift = this->RowShift();
    const Int rowShiftOfA = A.RowShift();

    // Compute which rowmajor rank has the rowShift equal to our rowShiftOfA
    const Int sendRankRM = (rankRM+(p+rowShiftOfA-rowShift)) % p;

    // Compute which rowmajor rank has the A rowShift that we need
    const Int recvRankCM = (rankCM+(p+rowShift-rowShiftOfA)) % p;
    const Int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

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
    ( sendBuf, sendSize, sendRankRM, 
      recvBuf, recvSize, recvRankRM, g.VRComm() );

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
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [VR,* ]");
#endif
    std::unique_ptr<DistMatrix<T,MR,MC>> A_MR_MC( new DistMatrix<T,MR,MC>(A) );
    std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(*A_MR_MC) );
    delete A_MR_MC.release(); // lowers memory highwater
    *this = *A_STAR_VC;
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = this->Grid();
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
            std::cerr << "Unaligned [* ,VR] <- [* ,VR]." << std::endl;
#endif
        const Int rank = g.VRRank();
        const Int p = g.Size();

        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const Int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

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
          recvBuf, recvSize, recvRank, g.VRComm() );

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
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,VR] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int p = g.Size();
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
        const T* ACol = &ABuf[(rowShift+jLoc*p)*ALDim];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        MemCopy( thisCol, ACol, localHeight );
    }
    return *this;
}

// NOTE: This is a small modification of [MC,MR] <- [o ,o ]
template<typename T>
const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [o ,o ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = g.Size();
    this->ResizeTo( m, n );

    const Int rowAlignment = this->RowAlignment();
    const Int nLocal = this->LocalWidth();
    const Int pkgSize = mpi::Pad(m*MaxLength(n,p));
    const Int recvSize = pkgSize;
    const Int sendSize = p*pkgSize;
    T* recvBuf;
    if( A.Participating() )
    {
        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        recvBuf = &buffer[sendSize];

        // Pack the send buffer
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        for( Int t=0; t<p; ++t )
        {
            const Int tLocalWidth = Length( n, t, p );
            const Int q = (rowAlignment+t) % p;
            for( Int jLoc=0; jLoc<tLocalWidth; ++jLoc )
            {
                const Int j = t + jLoc*p;
                for( Int i=0; i<m; ++i )
                    sendBuf[q*pkgSize+i+jLoc*m] = ABuf[i+j*ALDim];
            }
        }

        // Scatter from the root
        mpi::Scatter
        ( sendBuf, pkgSize, recvBuf, pkgSize, A.Root(), g.VRComm() );
    }
    else if( this->Participating() )
    {
        recvBuf = this->auxMemory_.Require( recvSize );

        // Perform the receiving portion of the scatter from the non-root
        mpi::Scatter
        ( static_cast<T*>(0), pkgSize, 
          recvBuf,            pkgSize, A.Root(), g.VRComm() );
    }

    if( this->Participating() )
    {
        // Unpack
        const Int ldim = this->LDim();
        T* buffer = this->Buffer();
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            for( Int i=0; i<m; ++i )
                buffer[i+jLoc*ldim] = recvBuf[i+jLoc*m];
        this->auxMemory_.Release();
    }

    return *this;
}

template<typename T>
void
DistMatrix<T,STAR,VR>::SumScatterFrom( const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::SumScatterFrom( [* ,MR] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetRowAlignmentAndResize( A.RowAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return;

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myCol = g.Col();
        const Int rowAlignment = this->RowAlignment();
        const Int rowShiftOfA = A.RowShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalWidth = MaxLength( width, p );
        const Int recvSize = mpi::Pad( height*maxLocalWidth );
        const Int sendSize = r*recvSize;

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        T* buffer = this->auxMemory_.Require( sendSize );
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = myCol+k*c;
            const Int thisRowShift = Shift_( thisRank, rowAlignment, p );
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const Int thisLocalWidth = Length_( width, thisRowShift, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* ACol = &ABuf[(thisRowOffset+jLoc*r)*ALDim];
                T* dataCol = &data[jLoc*height];
                MemCopy( dataCol, ACol, height );
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.ColComm() );

        // Unpack our received data
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* bufferCol = &buffer[jLoc*height];
            T* thisCol = &thisBuf[jLoc*thisLDim];
            MemCopy( thisCol, bufferCol, height );
        }
        this->auxMemory_.Release();
    }
    else
    {
        LogicError
        ("Unaligned [* ,VR]::ReduceScatterFrom( [* ,MR] ) not implemented");
    }
}

template<typename T>
void
DistMatrix<T,STAR,VR>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::SumScatterUpdate( [* ,MR] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myCol = g.Col();
        const Int rowAlignment = this->RowAlignment();
        const Int rowShiftOfA = A.RowShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalWidth = MaxLength( width, p );
        const Int recvSize = mpi::Pad( height*maxLocalWidth );
        const Int sendSize = r*recvSize;

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        T* buffer = this->auxMemory_.Require( sendSize );
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = myCol+k*c;
            const Int thisRowShift = Shift_( thisRank, rowAlignment, p );
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const Int thisLocalWidth = Length_( width, thisRowShift, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* ACol = &ABuf[(thisRowOffset+jLoc*r)*ALDim];
                T* dataCol = &data[jLoc*height];
                MemCopy( dataCol, ACol, height );
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.ColComm() );

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
        LogicError
        ("Unaligned [* ,VR]::ReduceScatterUpdate( [* ,MR] ) not implemented");
    }
}

//
// Routines which explicitly work in the complex plane
//

template<typename T>
void
DistMatrix<T,STAR,VR>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();
    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->SetLocalRealPart(i,jLoc,u);
    }
}

template<typename T>
void
DistMatrix<T,STAR,VR>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();
    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->SetLocalImagPart(i,jLoc,u);
    }
}

template<typename T>
void
DistMatrix<T,STAR,VR>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();
    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->UpdateLocalRealPart(i,jLoc,u);
    }
}

template<typename T>
void
DistMatrix<T,STAR,VR>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,VR]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();
    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->UpdateLocalImagPart(i,jLoc,u);
    }
}

#define PROTO(T) template class DistMatrix<T,STAR,VR>
#define COPY(T,CD,RD) \
  template DistMatrix<T,STAR,VR>::DistMatrix( const DistMatrix<T,CD,RD>& A )
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
  COPY(T,STAR,MR  ); \
  COPY(T,STAR,STAR); \
  COPY(T,STAR,VC  ); \
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
