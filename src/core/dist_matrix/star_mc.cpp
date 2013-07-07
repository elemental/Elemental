/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

template<typename T,typename Int>
DistMatrix<T,STAR,MC,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->SetShifts(); }

template<typename T,typename Int>
DistMatrix<T,STAR,MC,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T,typename Int>
DistMatrix<T,STAR,MC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->Align(0,rowAlignment); this->ResizeTo(height,width); }

template<typename T,typename Int>
DistMatrix<T,STAR,MC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->Align(0,rowAlignment); this->ResizeTo(height,width,ldim); }

template<typename T,typename Int>
DistMatrix<T,STAR,MC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, const T* buffer, Int ldim, 
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->LockedAttach(height,width,rowAlignment,buffer,ldim,g); }

template<typename T,typename Int>
DistMatrix<T,STAR,MC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, T* buffer, Int ldim, 
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->Attach(height,width,rowAlignment,buffer,ldim,g); }

template<typename T,typename Int>
DistMatrix<T,STAR,MC,Int>::DistMatrix( const DistMatrix<T,STAR,MC,Int>& A )
: AbstractDistMatrix<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,MC]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MC] with itself");
}

template<typename T,typename Int>
template<Distribution U,Distribution V>
DistMatrix<T,STAR,MC,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,MC]::DistMatrix");
#endif
    this->SetShifts();
    if( STAR != U || MC != V || 
        reinterpret_cast<const DistMatrix<T,STAR,MC,Int>*>(&A) != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MC] with itself");
}

template<typename T,typename Int>
DistMatrix<T,STAR,MC,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
elem::DistData<Int>
DistMatrix<T,STAR,MC,Int>::DistData() const
{
    elem::DistData<Int> data;
    data.colDist = STAR;
    data.rowDist = MC;
    data.colAlignment = 0;
    data.rowAlignment = this->rowAlignment_;
    data.root = 0;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T,typename Int>
Int
DistMatrix<T,STAR,MC,Int>::ColStride() const
{ return 1; }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,MC,Int>::RowStride() const
{ return this->grid_->Height(); }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,MC,Int>::ColRank() const
{ return 0; }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,MC,Int>::RowRank() const
{ return this->grid_->Row(); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::AlignWith( const elem::DistData<Int>& data )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::AlignWith");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( data.colDist == MC )
        this->rowAlignment_ = data.colAlignment;
    else if( data.rowDist == MC )
        this->rowAlignment_ = data.rowAlignment;
    else if( data.colDist == VC )
        this->rowAlignment_ = data.colAlignment % this->RowStride();
    else if( data.rowDist == VC )
        this->rowAlignment_ = data.rowAlignment % this->RowStride();
#ifndef RELEASE
    else throw std::logic_error("Nonsensical alignment");
#endif
    this->constrainedRowAlignment_ = true;
    this->SetShifts();
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::AlignWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::AlignRowsWith( const elem::DistData<Int>& data )
{ this->AlignWith( data ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::AlignRowsWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T,typename Int>
bool
DistMatrix<T,STAR,MC,Int>::AlignedWithDiagonal
( const elem::DistData<Int>& data, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::AlignedWithDiagonal");
#endif
    const Grid& grid = this->Grid();
    if( grid != *data.grid )
        return false;

    bool aligned;
    if( (data.colDist == MC   && data.rowDist == STAR) ||
        (data.colDist == STAR && data.rowDist == MC  ) )
    {
        const Int alignment = ( data.colDist==MC ? data.colAlignment
                                                 : data.rowAlignment );
        if( offset >= 0 )
        {
            const Int row = alignment;
            aligned = ( this->RowAlignment() == row );
        }
        else
        {
            const Int row = (alignment-offset) % this->RowStride();
            aligned = ( this->RowAlignment() == row );
        }
    }
    else aligned = false;
    return aligned;
}

template<typename T,typename Int>
bool
DistMatrix<T,STAR,MC,Int>::AlignedWithDiagonal
( const AbstractDistMatrix<T,Int>& A, Int offset ) const
{ return this->AlignedWithDiagonal( A.DistData(), offset ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::AlignWithDiagonal
( const elem::DistData<Int>& data, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::AlignWithDiagonal");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( (data.colDist == MC   && data.rowDist == STAR) ||
        (data.colDist == STAR && data.rowDist == MC  ) )
    {
        const Int alignment = ( data.colDist==MC ? data.colAlignment
                                                 : data.rowAlignment );
        if( offset >= 0 )
        {
            const Int row = alignment;
            this->rowAlignment_ = row;
        }
        else
        {
            const Int row = (alignment-offset) % this->RowStride();
            this->rowAlignment_ = row;
        }
        this->constrainedRowAlignment_ = true;
        this->SetShifts();
    }
#ifndef RELEASE
    else throw std::logic_error("Invalid diagonal alignment");
#endif
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::AlignWithDiagonal
( const AbstractDistMatrix<T,Int>& A, Int offset )
{ this->AlignWithDiagonal( A.DistData(), offset ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::Attach
( Int height, Int width, Int rowAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::Attach");
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
        const Int localWidth = Length(width,this->rowShift_,g.Height());
        this->matrix_.Attach_( height, localWidth, buffer, ldim );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::LockedAttach
( Int height, Int width, Int rowAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::LockedAttach");
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
        const Int localWidth = Length(width,this->rowShift_,g.Height());
        this->matrix_.LockedAttach_( height, localWidth, buffer, ldim );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( height, Length(width,this->RowShift(),this->Grid().Height()) );
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( height, Length(width,this->RowShift(),this->Grid().Height()), ldim );
}

template<typename T,typename Int>
T
DistMatrix<T,STAR,MC,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MR]::Get");
    this->AssertValidEntry( i, j );
    if( !this->Participating() )
        throw std::logic_error("Should only be called by processes in grid");
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that
    // row within each process column
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    T u;
    if( g.Row() == ownerRow )
    {
        const Int jLoc = (j-this->RowShift()) / g.Height();
        u = this->GetLocal(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRow, g.ColComm() );
    return u;
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->SetLocal(i,jLoc,u);
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateLocal(i,jLoc,u);
    }
}

//
// Utility functions, e.g., SumOverRow
//

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::SumOverRow()
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::SumOverRow");
    this->AssertNotLocked();
#endif
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* thisCol = &thisBuf[jLoc*thisLDim];
        T* sendBufCol = &sendBuf[jLoc*localHeight];
        MemCopy( sendBufCol, thisCol, localHeight );
    }

    // AllReduce sum
    mpi::AllReduce
    ( sendBuf, recvBuf, localSize, mpi::SUM, this->Grid().RowComm() );

    // Unpack
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* recvBufCol = &recvBuf[jLoc*localHeight];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        MemCopy( thisCol, recvBufCol, localHeight );
    }
    this->auxMemory_.Release();
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::AdjointFrom( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::AdjointFrom");
#endif
    this->TransposeFrom( A, true );
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::TransposeFrom
( const DistMatrix<T,VC,STAR,Int>& A, bool conjugate )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::TransposeFrom");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Width(), A.Height() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.ColAlignment() % g.Height();
            this->SetRowShift();
        }
        this->ResizeTo( A.Width(), A.Height() );
    }
    if( !this->Participating() )
        return;

    if( this->RowAlignment() == A.ColAlignment() % g.Height() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalHeightOfA );

        T* buffer = this->auxMemory_.Require( (c+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        if( conjugate )
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
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
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
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
          recvBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int rowShift = this->RowShift();
        const Int colAlignmentOfA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int colShiftOfA = Shift_( row+k*r, colAlignmentOfA, p );
            const Int rowOffset = (colShiftOfA-rowShift) / r;
            const Int localWidth = Length_( width, colShiftOfA, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*c)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MC]::AdjointFrom." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int rank = g.VCRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowShift = this->RowShift();

        const Int sendRank = (rank+p+rowAlignment-colAlignmentOfA) % p;
        const Int recvRank = (rank+p+colAlignmentOfA-rowAlignment) % p;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalHeightOfA );

        T* buffer = this->auxMemory_.Require( (c+1)*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        if( conjugate )
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
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
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
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
        ( secondBuf, portionSize, sendRank, 0,
          firstBuf,  portionSize, recvRank, 0, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShiftOfA = Shift_(row+r*k,rowAlignment,p);
            const Int rowOffset = (colShiftOfA-rowShift) / r;
            const Int localWidth = Length_( width, colShiftOfA, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*c)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR,Int>(g) );
    *A_STAR_VR = A;

    std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC,Int>(true,this->RowAlignment(),g) );
    *A_STAR_VC = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_STAR_VC;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MC,MR,Int> > A_MC_MR
    ( new DistMatrix<T,MC,MR,Int>(g) );
    *A_MC_MR   = A;

    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR,Int>(g) );
    *A_STAR_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC,Int>(true,this->RowAlignment(),g) );
    *A_STAR_VC = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_STAR_VC;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( A.Height() == 1 )
    {
        if( !this->Viewing() )
            this->ResizeTo( 1, A.Width() );
        if( !this->Participating() )
            return *this;

        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int myRow = g.Row();
        const Int rankCM = g.VCRank();
        const Int rankRM = g.VRRank();
        const Int rowAlignment = this->RowAlignment();
        const Int rowShift = this->RowShift();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int rowShiftOfA = A.RowShift();

        const Int width = this->Width();
        const Int maxLocalVectorWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxLocalVectorWidth );

        const Int rowShiftVC = Shift(rankCM,rowAlignment,p);
        const Int rowShiftVROfA = Shift(rankRM,rowAlignmentOfA,p);
        const Int sendRankCM = (rankCM+(p+rowShiftVROfA-rowShiftVC)) % p;
        const Int recvRankRM = (rankRM+(p+rowShiftVC-rowShiftVROfA)) % p;
        const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        T* buffer = this->auxMemory_.Require( (c+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        // A[* ,VR] <- A[* ,MR]
        {
            const Int shift = Shift(rankRM,rowAlignmentOfA,p);
            const Int offset = (shift-rowShiftOfA) / c;
            const Int thisLocalWidth = Length(width,shift,p);

            const T* ABuf = A.LockedBuffer();
            const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                sendBuf[jLoc] = ABuf[(offset+jLoc*r)*ALDim];
        }

        // A[* ,VC] <- A[* ,VR]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankCM, 0,
          recvBuf, portionSize, recvRankCM, mpi::ANY_TAG, g.VCComm() );

        // A[* ,MC] <- A[* ,VC]
        mpi::AllGather
        ( recvBuf, portionSize,
          sendBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &sendBuf[k*portionSize];
            const Int shift = Shift_(myRow+r*k,rowAlignment,p);
            const Int offset = (shift-rowShift) / r;
            const Int thisLocalWidth = Length_(width,shift,p);
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                thisBuf[(offset+jLoc*c)*thisLDim] = data[jLoc];
        }
        this->auxMemory_.Release();
    }
    else
    {
        std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
        ( new DistMatrix<T,STAR,VR,Int>(g) );
        *A_STAR_VR = A;

        std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
        ( new DistMatrix<T,STAR,VC,Int>(true,this->RowAlignment(),g) );
        *A_STAR_VC = *A_STAR_VR;
        delete A_STAR_VR.release(); // lowers memory highwater

        *this = *A_STAR_VC;
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MC] = [MD,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MC] = [MD,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Height() == 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "The vector version of [* ,MC] <- [MR,MC] is not yet written, but"
          " it would only require a modification of the vector version of "
          "[* ,MR] <- [MC,MR]." << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "The redistribution [* ,MC] <- [MR,MC] potentially causes a large"
          " amount of cache-thrashing. If possible, avoid it. "
          "Unfortunately, the following routines are not yet implemented: \n"
          << "  [MC,* ].(Conjugate)TransposeFrom([MR,MC])" << std::endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            this->SetRowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !this->Participating() )
        return *this;

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const Int c = g.Width();
        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLength(height,c);
        const Int portionSize = mpi::Pad( maxLocalHeightOfA*localWidth );

        T* buffer = this->auxMemory_.Require( (c+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* ACol = &ABuf[jLoc*ALDim];
            T* sendBufCol = &sendBuf[jLoc*localHeightOfA];
            MemCopy( sendBufCol, ACol, localHeightOfA );
        }

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.RowComm() );

        // Unpack
        const Int colAlignmentOfA = A.ColAlignment();
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int colShift = Shift_( k, colAlignmentOfA, c );
            const Int localHeight = Length_( height, colShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[colShift+jLoc*thisLDim];
                const T* sourceCol = &data[jLoc*localHeight];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc*c] = sourceCol[iLoc];
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MC] <- [MR,MC]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.Row();

        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalHeightOfA = MaxLength(height,c);
        const Int maxLocalWidth = MaxLength(width,r);
        const Int portionSize = mpi::Pad( maxLocalHeightOfA*maxLocalWidth );

        T* buffer = this->auxMemory_.Require( (c+1)*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuf[jLoc*ALDim];
            T* secondBufCol = &secondBuf[jLoc*localHeightOfA];
            MemCopy( secondBufCol, ACol, localHeightOfA );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuf, portionSize, sendRow, 0,
          firstBuf,  portionSize, recvRow, mpi::ANY_TAG, g.ColComm() );

        // Use the output of the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.RowComm() );

        // Unpack the contents of each member of the process row
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int colAlignmentOfA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShift = Shift_( k, colAlignmentOfA, c );
            const Int localHeight = Length_( height, colShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[colShift+jLoc*thisLDim];
                const T* sourceCol = &data[jLoc*localHeight];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc*c] = sourceCol[iLoc];
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MR,MC,Int> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [* ,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            if( this->Participating() )
                this->rowShift_ = A.RowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
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
            std::cerr << "Unaligned [* ,MC] <- [* ,MC]." << std::endl;
#endif
        const Int rank = g.Row();
        const Int r = g.Height();

        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRank = (rank+r+rowAlignment-rowAlignmentOfA) % r;
        const Int recvRank = (rank+r+rowAlignmentOfA-rowAlignment) % r;

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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuf[jLoc*ALDim];
            T* sendBufCol = &sendBuf[jLoc*height];
            MemCopy( sendBufCol, ACol, height );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank, 0,
          recvBuf, recvSize, recvRank, mpi::ANY_TAG, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
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

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = A;

    std::auto_ptr<DistMatrix<T,MR,MC,Int> > A_MR_MC
    ( new DistMatrix<T,MR,MC,Int>(false,true,0,this->RowAlignment(),g) );
    *A_MR_MC = *A_VR_STAR;
    delete A_VR_STAR.release();

    *this = *A_MR_MC;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [* ,VC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment() % g.Height();
            this->SetRowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !this->Participating() )
        return *this;

    if( this->RowAlignment() == A.RowAlignment() % g.Height() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidthOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalWidthOfA );

        T* buffer = this->auxMemory_.Require( (c+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuf[jLoc*ALDim];
            T* sendBufCol = &sendBuf[jLoc*height];
            MemCopy( sendBufCol, ACol, height );
        }

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int rowShift = this->RowShift();
        const Int rowAlignmentOfA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int rowShiftOfA = Shift_( row+k*r, rowAlignmentOfA, p );
            const Int rowOffset = (rowShiftOfA-rowShift) / r;
            const Int localWidth = Length_( width, rowShiftOfA, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*c)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MC] <- [* ,VC]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int rank = g.VCRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int rowShift = this->RowShift();

        const Int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const Int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidthOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalWidthOfA );

        T* buffer = this->auxMemory_.Require( (c+1)*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuf[jLoc*ALDim];
            T* secondBufCol = &secondBuf[jLoc*height];
            MemCopy( secondBufCol, ACol, height );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuf, portionSize, sendRank, 0,
          firstBuf,  portionSize, recvRank, 0, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShiftOfA = Shift_(row+r*k,rowAlignment,p);
            const Int rowOffset = (rowShiftOfA-rowShift) / r;
            const Int localWidth = Length_( width, rowShiftOfA, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*c)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MR,MC,Int> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VC,Int> A_STAR_VC(true,this->RowAlignment(),g);
    *this = A_STAR_VC = A;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int r = this->Grid().Height();
    const Int rowShift = this->RowShift();

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();

    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* ACol = &ABuf[(rowShift+jLoc*r)*ALDim];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        MemCopy( thisCol, ACol, localHeight );
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,CIRC,CIRC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,MC] = [o ,o ]");
#endif
    DistMatrix<T,MR,MC> A_MR_MC( A.Grid() );
    A_MR_MC.AlignWith( *this );
    A_MR_MC = A;
    *this = A_MR_MC;
    return *this;
}

//
// Routines which explicitly work in the complex plane
//

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->SetLocalRealPart( i, jLoc, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->SetLocalImagPart( i, jLoc, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateLocalRealPart( i, jLoc, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MC,Int>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MC]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateLocalImagPart( i, jLoc, u );
    }
}

#define PROTO(T) \
  template class DistMatrix<T,STAR,MC,int>
#define COPY(T,CD,RD) \
  template DistMatrix<T,STAR,MC,int>::DistMatrix( \
    const DistMatrix<T,CD,RD,int>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  MR  ); \
  COPY(T,MC,  STAR); \
  COPY(T,MD,  STAR); \
  COPY(T,MR,  MC  ); \
  COPY(T,MR,  STAR); \
  COPY(T,STAR,MD  ); \
  COPY(T,STAR,MR  ); \
  COPY(T,STAR,STAR); \
  COPY(T,STAR,VC  ); \
  COPY(T,STAR,VR  ); \
  COPY(T,VC,  STAR); \
  COPY(T,VR,  STAR);

FULL(int);
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
