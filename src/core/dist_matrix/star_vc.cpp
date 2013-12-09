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
using ADM = AbstractDistMatrix<T>;
template<typename T>
using DM = DistMatrix<T,STAR,VC>;

template<typename T>
DM<T>::DistMatrix( const elem::Grid& g )
: ADM<T>(g)
{ this->SetShifts(); }

template<typename T>
DM<T>::DistMatrix( Int height, Int width, const elem::Grid& g )
: ADM<T>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int rowAlign, const elem::Grid& g )
: ADM<T>(g)
{ this->Align(0,rowAlign); this->ResizeTo(height,width); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int rowAlign, Int ldim, const elem::Grid& g )
: ADM<T>(g)
{ this->Align(0,rowAlign); this->ResizeTo(height,width,ldim); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int rowAlign, const T* buffer, Int ldim,
  const elem::Grid& g )
: ADM<T>(g)
{ this->LockedAttach(height,width,rowAlign,buffer,ldim,g); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int rowAlign, T* buffer, Int ldim,
  const elem::Grid& g )
: ADM<T>(g)
{ this->Attach(height,width,rowAlign,buffer,ldim,g); }

template<typename T>
DM<T>::DistMatrix( const DM<T>& A )
: ADM<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[* ,VC]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [* ,VC] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DM<T>::DistMatrix( const DistMatrix<T,U,V>& A )
: ADM<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[* ,VC]::DistMatrix");
#endif
    this->SetShifts();
    if( STAR != U || VC != V || 
        reinterpret_cast<const DM<T>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [* ,VC] with itself");
}

template<typename T>
DM<T>::DistMatrix( DM<T>&& A )
: ADM<T>(std::move(A))
{ }

template<typename T>
DM<T>&
DM<T>::operator=( DM<T>&& A )
{
    ADM<T>::operator=( std::move(A) );
    return *this;
}

template<typename T>
DM<T>::~DistMatrix()
{ }

template<typename T>
elem::DistData
DM<T>::DistData() const
{ return elem::DistData(*this); }

template<typename T>
mpi::Comm
DM<T>::DistComm() const
{ return this->grid_->VCComm(); }

template<typename T>
mpi::Comm
DM<T>::CrossComm() const
{ return mpi::COMM_SELF; }

template<typename T>
mpi::Comm
DM<T>::RedundantComm() const
{ return mpi::COMM_SELF; }

template<typename T>
mpi::Comm
DM<T>::ColComm() const
{ return mpi::COMM_SELF; }

template<typename T>
mpi::Comm
DM<T>::RowComm() const
{ return this->grid_->VCComm(); }

template<typename T>
Int
DM<T>::ColStride() const
{ return 1; }
    
template<typename T>
Int 
DM<T>::RowStride() const
{ return this->grid_->Size(); }

template<typename T>
void
DM<T>::AlignWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,VC]::AlignWith");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( data.colDist == MC || data.colDist == VC )
        this->rowAlign_ = data.colAlign;
    else if( data.rowDist == MC || data.rowDist == VC )
        this->rowAlign_ = data.rowAlign;
#ifndef RELEASE
    else LogicError("Nonsensical alignment");
#endif
    this->rowConstrained_ = true;
    this->SetShifts();
}

template<typename T>
void
DM<T>::AlignWith( const ADM<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DM<T>::AlignRowsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

template<typename T>
void
DM<T>::AlignRowsWith( const ADM<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
bool
DM<T>::AlignedWithDiagonal
( const elem::DistData& data, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[* ,VC]::AlignedWithDiagonal");
#endif
    const Grid& grid = this->Grid();
    if( grid != *data.grid )
        return false;

    bool aligned;
    if( (data.colDist == VC   && data.rowDist == STAR) ||
        (data.colDist == STAR && data.rowDist == VC  ) )
    {
        const Int alignment = ( data.colDist==VC ? data.colAlign
                                                 : data.rowAlign );
        if( offset >= 0 )
        {
            const Int proc = alignment;
            aligned = ( this->RowAlign() == proc );
        }
        else
        {
            const Int proc = (alignment-offset) % this->RowStride();
            aligned = ( this->RowAlign() == proc );
        }
    }
    else aligned = false;
    return aligned;
}

template<typename T>
bool
DM<T>::AlignedWithDiagonal
( const ADM<T>& A, Int offset ) const
{ return this->AlignedWithDiagonal( A.DistData(), offset ); }

template<typename T>
void
DM<T>::AlignWithDiagonal
( const elem::DistData& data, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,VC]::AlignWithDiagonal");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( (data.colDist == VC   && data.rowDist == STAR) ||
        (data.colDist == STAR && data.rowDist == VC  ) )
    {
        const Int alignment = ( data.colDist==VC ? data.colAlign
                                                 : data.rowAlign );
        if( offset >= 0 )
        {
            const Int proc = alignment;
            this->rowAlign_ = proc;
        }
        else
        {
            const Int proc = (alignment-offset) % this->RowStride();
            this->rowAlign_ = proc;
        }
        this->rowConstrained_ = true;
        this->SetShifts();
    }
#ifndef RELEASE
    else LogicError("Invalid diagonal alignment");
#endif
}

template<typename T>
void
DM<T>::AlignWithDiagonal
( const ADM<T>& A, Int offset )
{ this->AlignWithDiagonal( A.DistData(), offset ); }

template<typename T>
void
DM<T>::Attach
( Int height, Int width, Int rowAlign,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,VC]::Attach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->rowAlign_ = rowAlign;
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
DM<T>::LockedAttach
( Int height, Int width, Int rowAlign,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,VC]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->rowAlign_ = rowAlign;
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
DM<T>::Attach( Matrix<T>& A, Int rowAlign, const elem::Grid& g )
{ this->Attach( A.Height(), A.Width(), rowAlign, A.Buffer(), A.LDim(), g ); }

template<typename T>
void
DM<T>::LockedAttach( const Matrix<T>& A, Int rowAlign, const elem::Grid& g )
{
    this->LockedAttach
    ( A.Height(), A.Width(), rowAlign, A.LockedBuffer(), A.LDim(), g );
}

//
// Utility functions, e.g., operator=
//

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [MC,MR]");
#endif
    DistMatrix<T,STAR,VR> A_STAR_VR( A );
    *this = A_STAR_VR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [MC,* ]");
#endif
    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR( new DistMatrix<T,MC,MR>(A) );
    std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(*A_MC_MR) );
    delete A_MC_MR.release(); // lowers memory highwater
    *this = *A_STAR_VR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [* ,MR]");
#endif
    DistMatrix<T,STAR,VR> A_STAR_VR( A );
    *this = A_STAR_VR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [MD,* ]");
#endif
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [* ,MD]");
#endif
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->AlignRowsAndResize( A.RowAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlign() % g.Height() == A.RowAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlign = this->RowAlign();
        const Int colAlignOfA = A.ColAlign();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxHeight = MaxLength(height,c);
        const Int maxWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*c*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        // Pack
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuf[k*portionSize];
            const Int thisRank = row+k*r;
            const Int thisRowShift = Shift_(thisRank,rowAlign,p);
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* ACol = &ABuf[(thisRowOffset+jLoc*c)*ALDim];
                T* dataCol = &data[jLoc*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuf, portionSize, recvBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisColShift = Shift_(k,colAlignOfA,c);
            const Int thisLocalHeight = Length_(height,thisColShift,c);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[thisColShift+jLoc*thisLDim];
                const T* sourceCol = &data[jLoc*thisLocalHeight];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc*c] = sourceCol[iLoc];
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VC] <- [MR,MC]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlign = this->RowAlign();
        const Int colAlignOfA = A.ColAlign();
        const Int rowAlignOfA = A.RowAlign();

        const Int sendRow = (row+r+(rowAlign%r)-rowAlignOfA) % r;
        const Int recvRow = (row+r+rowAlignOfA-(rowAlign%r)) % r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int maxHeight = MaxLength(height,c);
        const Int maxWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*c*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[c*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int thisRank = sendRow+k*r;
            const Int thisRowShift = Shift_(thisRank,rowAlign,p);
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* ACol = &ABuf[(thisRowOffset+jLoc*c)*ALDim];
                T* dataCol = &data[jLoc*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // AllToAll to gather all of the unaligned [*,VC] data into firstBuf
        mpi::AllToAll
        ( secondBuf, portionSize, firstBuf, portionSize, g.RowComm() );

        // SendRecv: properly align the [*,VC] via a trade in the column
        mpi::SendRecv
        ( firstBuf,  portionSize, sendRow, 
          secondBuf, portionSize, recvRow, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int thisColShift = Shift_(k,colAlignOfA,c);
            const Int thisLocalHeight = Length_(height,thisColShift,c);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[thisColShift+jLoc*thisLDim];
                const T* sourceCol = &data[jLoc*thisLocalHeight];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc*c] = sourceCol[iLoc];
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [MR,* ]");
#endif
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [* ,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->AlignRowsAndResize( A.RowAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlign() % g.Height() == A.RowAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int rowShift = this->RowShift();
        const Int rowShiftOfA = A.RowShift();
        const Int rowOffset = (rowShift-rowShiftOfA) / r;

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();

        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* ACol = &ABuf[(rowOffset+jLoc*c)*ALDim];
            T* thisCol = &thisBuf[jLoc*thisLDim];
            MemCopy( thisCol, ACol, height );
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VC] <- [* ,MC]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int col = g.Col();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlign = this->RowAlign();
        const Int rowAlignOfA = A.RowAlign();

        // We will SendRecv A[*,VC] within our process column to fix alignments.
        const Int sendRow = (row+r+(rowAlign%r)-rowAlignOfA) % r;
        const Int recvRow = (row+r+rowAlignOfA-(rowAlign%r)) % r;
        const Int sendRank = sendRow + r*col;

        const Int sendRowShift = Shift( sendRank, rowAlign, p );
        const Int sendRowOffset = (sendRowShift-rowShiftOfA) / r;

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
            const T* ACol = &ABuf[(sendRowOffset+jLoc*c)*ALDim];
            T* sendBufCol = &sendBuf[jLoc*height];
            MemCopy( sendBufCol, ACol, height );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRow, 
          recvBuf, recvSize, recvRow, g.ColComm() );

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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [VC,* ]");
#endif
    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR( new DistMatrix<T,MC,MR>(A) );
    std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(*A_MC_MR) );
    delete A_MC_MR.release(); // lowers memory highwater
    *this = *A_STAR_VR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DM<T>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [* ,VC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid(); 
    this->AlignRowsAndResize( A.RowAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlign() == A.RowAlign() )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VC] <- [* ,VC]." << std::endl;
#endif
        const Int rank = g.VCRank();
        const Int p = g.Size();

        const Int rowAlign = this->RowAlign();
        const Int rowAlignOfA = A.RowAlign();

        const Int sendRank = (rank+p+rowAlign-rowAlignOfA) % p;
        const Int recvRank = (rank+p+rowAlignOfA-rowAlign) % p;

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
          recvBuf, recvSize, recvRank, g.VCComm() );

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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [VR,* ]");
#endif
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = this->Grid();
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

    // Compute which colmajor rank has the rowShift equal to our rowShiftOfA
    const Int sendRankCM = (rankCM+(p+rowShiftOfA-rowShift)) % p;

    // Compute which colmajor rank has the A rowShift that we need
    const Int recvRankRM = (rankRM+(p+rowShift-rowShiftOfA)) % p;
    const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

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
    ( sendBuf, sendSize, sendRankCM, 
      recvBuf, recvSize, recvRankCM, g.VCComm() );

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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,VC] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int p = g.Size();
    const Int rowShift = this->RowShift();

    const Int height = this->Height();
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
        MemCopy( thisCol, ACol, height );
    }
    return *this;
}

// NOTE: This is a small modification of [MC,MR] <- [o ,o ]
template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [o ,o ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = g.Size();
    this->ResizeTo( m, n );

    const Int rowAlign = this->RowAlign();
    const Int nLocal = this->LocalWidth();
    const Int pkgSize = mpi::Pad(m*MaxLength(n,p));
    const Int recvSize = pkgSize;
    const Int sendSize = p*pkgSize;
    T* recvBuf=0; // some compilers (falsely) warn otherwise
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
            const Int q = (rowAlign+t) % p;
            for( Int jLoc=0; jLoc<tLocalWidth; ++jLoc )
            {
                const Int j = t + jLoc*p;
                for( Int i=0; i<m; ++i )
                    sendBuf[q*pkgSize+i+jLoc*m] = ABuf[i+j*ALDim];
            }
        }

        // Scatter from the root
        mpi::Scatter
        ( sendBuf, pkgSize, recvBuf, pkgSize, A.Root(), g.VCComm() );
    }
    else if( this->Participating() )
    {
        recvBuf = this->auxMemory_.Require( recvSize );

        // Perform the receiving portion of the scatter from the non-root
        mpi::Scatter
        ( static_cast<T*>(0), pkgSize, 
          recvBuf,            pkgSize, A.Root(), g.VCComm() );
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
DM<T>::SumScatterFrom( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,VC]::SumScatterFrom( [* ,MC] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->AlignRowsAndResize( A.RowAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return;

    if( this->RowAlign() % g.Height() == A.RowAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int rowAlign = this->RowAlign();
        const Int rowShiftOfA = A.RowShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalWidth = MaxLength( width, p );

        const Int recvSize = mpi::Pad( height*maxLocalWidth );
        const Int sendSize = c*recvSize;

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        T* buffer = this->auxMemory_.Require( sendSize );
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = myRow+k*r;
            const Int thisRowShift = Shift_( thisRank, rowAlign, p );
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = Length_( width, thisRowShift, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* ACol = &ABuf[(thisRowOffset+jLoc*c)*ALDim];
                T* dataCol = &data[jLoc*height];
                MemCopy( dataCol, ACol, height );
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.RowComm() );

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
        ("Unaligned [* ,VC]::ReduceScatterFrom( [* ,MC] ) not yet implemented");
    }
}

template<typename T>
void
DM<T>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,VC]::SumScatterUpdate( [* ,MC] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    if( this->RowAlign() % g.Height() == A.RowAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int rowAlign = this->RowAlign();
        const Int rowShiftOfA = A.RowShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalWidth = MaxLength( width, p );
        const Int recvSize = mpi::Pad( height*maxLocalWidth );
        const Int sendSize = c*recvSize;

        // Pack
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        T* buffer = this->auxMemory_.Require( sendSize );
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = myRow+k*r;
            const Int thisRowShift = Shift_( thisRank, rowAlign, p );
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = Length_( width, thisRowShift, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* ACol = &ABuf[(thisRowOffset+jLoc*c)*ALDim];
                T* dataCol = &data[jLoc*height];
                MemCopy( dataCol, ACol, height );
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.RowComm() );

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
        ("Unaligned [* ,VC]::ReduceScatterUpdate( [* ,MC] ) not implemented");
    }
}

#define PROTO(T) template class DistMatrix<T,STAR,VC>
#define COPY(T,U,V) \
  template DistMatrix<T,STAR,VC>::DistMatrix( const DistMatrix<T,U,V>& A )
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
