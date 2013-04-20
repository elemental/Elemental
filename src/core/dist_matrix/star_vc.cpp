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
DistMatrix<T,STAR,VC,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   0,(g.InGrid() ? g.VCRank() : 0 ),
   0,0,g)
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,VC,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   0,(g.InGrid() ? g.VCRank() : 0),
   height,(g.InGrid() ? Length(width,g.VCRank(),0,g.Size()) : 0),
   g)
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,VC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,true,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VCRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? Length(width,g.VCRank(),rowAlignment,g.Size()) : 0),
   g)
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,VC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,true,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VCRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? Length(width,g.VCRank(),rowAlignment,g.Size()) : 0),
   ldim,g)
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,VC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VCRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? Length(width,g.VCRank(),rowAlignment,g.Size()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,VC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VCRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? Length(width,g.VCRank(),rowAlignment,g.Size()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,VC,Int>::DistMatrix( const DistMatrix<T,STAR,VC,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  0,(A.Participating() ? A.RowRank() : 0),
  0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,VC] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<Distribution U,Distribution V>
DistMatrix<T,STAR,VC,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  0,(A.Participating() ? A.RowRank() : 0),
  0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::DistMatrix");
#endif
    if( STAR != U || VC != V || 
        reinterpret_cast<const DistMatrix<T,STAR,VC,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,VC] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
DistMatrix<T,STAR,VC,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
elem::DistData<Int>
DistMatrix<T,STAR,VC,Int>::DistData() const
{
    elem::DistData<Int> data;
    data.colDist = STAR;
    data.rowDist = VC;
    data.colAlignment = 0;
    data.rowAlignment = this->rowAlignment_;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T,typename Int>
Int
DistMatrix<T,STAR,VC,Int>::ColStride() const
{ return 1; }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,VC,Int>::RowStride() const
{ return this->grid_->Size(); }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,VC,Int>::ColRank() const
{ return 0; }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,VC,Int>::RowRank() const
{ return this->grid_->VCRank(); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::AlignWith( const elem::DistData<Int>& data )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith");
    this->AssertFreeRowAlignment();
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( data.colDist == MC || data.colDist == VC )
        this->rowAlignment_ = data.colAlignment;
    else if( data.rowDist == MC || data.rowDist == VC )
        this->rowAlignment_ = data.rowAlignment;
#ifndef RELEASE
    else throw std::logic_error("Nonsensical alignment");
#endif
    this->constrainedRowAlignment_ = true;
    this->SetShifts();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::AlignWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::AlignRowsWith( const elem::DistData<Int>& data )
{ this->AlignWith( data ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::AlignRowsWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T,typename Int>
bool
DistMatrix<T,STAR,VC,Int>::AlignedWithDiagonal
( const elem::DistData<Int>& data, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignedWithDiagonal");
#endif
    const Grid& grid = this->Grid();
    if( grid != *data.grid )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return false;
    }

    bool aligned;
    if( (data.colDist == VC   && data.rowDist == STAR) ||
        (data.colDist == STAR && data.rowDist == VC  ) )
    {
        const Int alignment = ( data.colDist==VC ? data.colAlignment
                                                 : data.rowAlignment );
        if( offset >= 0 )
        {
            const Int proc = alignment;
            aligned = ( this->RowAlignment() == proc );
        }
        else
        {
            const Int proc = (alignment-offset) % this->RowStride();
            aligned = ( this->RowAlignment() == proc );
        }
    }
    else aligned = false;
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
bool
DistMatrix<T,STAR,VC,Int>::AlignedWithDiagonal
( const AbstractDistMatrix<T,Int>& A, Int offset ) const
{ return this->AlignedWithDiagonal( A.DistData(), offset ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::AlignWithDiagonal
( const elem::DistData<Int>& data, Int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWithDiagonal");
    this->AssertFreeRowAlignment();
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( (data.colDist == VC   && data.rowDist == STAR) ||
        (data.colDist == STAR && data.rowDist == VC  ) )
    {
        const Int alignment = ( data.colDist==VC ? data.colAlignment
                                                 : data.rowAlignment );
        if( offset >= 0 )
        {
            const Int proc = alignment;
            this->rowAlignment_ = proc;
        }
        else
        {
            const Int proc = (alignment-offset) % this->RowStride();
            this->rowAlignment_ = proc;
        }
        this->constrainedRowAlignment_ = true;
        this->SetShifts();
    }
#ifndef RELEASE
    else throw std::logic_error("Invalid diagonal alignment");
    PopCallStack();
#endif
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::AlignWithDiagonal
( const AbstractDistMatrix<T,Int>& A, Int offset )
{ this->AlignWithDiagonal( A.DistData(), offset ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::PrintBase");
#endif
    const elem::Grid& g = this->Grid();
    if( g.Rank() == 0 && msg != "" )
        os << msg << std::endl;

    const Int height     = this->Height();
    const Int width      = this->Width();
    const Int localWidth = this->LocalWidth();
    const Int p          = g.Size();
    const Int rowShift   = this->RowShift();

    if( height == 0 || width == 0 || !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    std::vector<T> sendBuf(height*width,0);
    const T* thisBuffer = this->LockedBuffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        T* destCol = &sendBuf[(rowShift+jLocal*p)*height];
        const T* sourceCol = &thisBuffer[jLocal*thisLDim];
        MemCopy( destCol, sourceCol, height );
    }

    // If we are the root, allocate a receive buffer
    std::vector<T> recvBuf;
    if( g.Rank() == 0 )
        recvBuf.resize( height*width );

    // Sum the contributions and send to the root
    mpi::Reduce
    ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, g.Comm() );

    if( g.Rank() == 0 )
    {
        // Print the data
        for( Int i=0; i<height; ++i )
        {
            for( Int j=0; j<width; ++j )
                os << recvBuf[i+j*height] << " ";
            os << "\n";
        }
        os << std::endl;
    }
    mpi::Barrier( g.Comm() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::Attach
( Int height, Int width, Int rowAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Attach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->viewing_ = true;
    this->SetRowShift();
    if( g.InGrid() )
    {
        const Int localWidth = Length(width,this->rowShift_,g.Size());
        this->matrix_.Attach( height, localWidth, buffer, ldim );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::LockedAttach
( Int height, Int width, Int rowAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->viewing_ = true;
    this->locked_ = true;
    this->SetRowShift();
    if( g.InGrid() )
    {
        const Int localWidth = Length(width,this->rowShift_,g.Size());
        this->matrix_.LockedAttach( height, localWidth, buffer, ldim );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    const elem::Grid& g = this->Grid();
    this->height_ = height;
    this->width_  = width;
    if( g.InGrid() )
        this->matrix_.ResizeTo
        ( height, Length( width, this->RowShift(), g.Size() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
T
DistMatrix<T,STAR,VC,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocal(i,jLoc);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->SetLocal(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->UpdateLocal(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., operator=
//

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR,Int> A_STAR_VR(g);

    A_STAR_VR = A;
    *this = A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MC,MR,Int> > A_MC_MR
    ( new DistMatrix<T,MC,MR,Int>(g) );
    *A_MC_MR = A;

    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR,Int>(g) );
    *A_STAR_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR,Int> A_STAR_VR(g);

    A_STAR_VR = A;
    *this = A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MR,MC]");
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
            this->SetRowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return *this;
    }

    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int maxHeight = MaxLength(height,c);
        const Int maxWidth = MaxLength(width,p);
        const Int portionSize = std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( 2*c*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[c*portionSize];

        // Pack
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const Int thisRank = row+k*r;
            const Int thisRowShift = Shift_(thisRank,rowAlignment,p);
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ABuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, g.RowComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const Int thisColShift = Shift_(k,colAlignmentOfA,c);
            const Int thisLocalHeight = Length_(height,thisColShift,c);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                T* destCol = &thisBuffer[thisColShift+jLocal*thisLDim];
                const T* sourceCol = &data[jLocal*thisLocalHeight];
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    destCol[iLocal*c] = sourceCol[iLocal];
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
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRow = (row+r+(rowAlignment%r)-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-(rowAlignment%r)) % r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int maxHeight = MaxLength(height,c);
        const Int maxWidth = MaxLength(width,p);
        const Int portionSize = std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( 2*c*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[c*portionSize];

        // Pack
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const Int thisRank = sendRow+k*r;
            const Int thisRowShift = Shift_(thisRank,rowAlignment,p);
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ABuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // AllToAll to gather all of the unaligned [*,VC] data into firstBuffer
        mpi::AllToAll
        ( secondBuffer, portionSize,
          firstBuffer,  portionSize, g.RowComm() );

        // SendRecv: properly align the [*,VC] via a trade in the column
        mpi::SendRecv
        ( firstBuffer,  portionSize, sendRow, 0,
          secondBuffer, portionSize, recvRow, mpi::ANY_TAG, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const Int thisColShift = Shift_(k,colAlignmentOfA,c);
            const Int thisLocalHeight = Length_(height,thisColShift,c);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                T* destCol = &thisBuffer[thisColShift+jLocal*thisLDim];
                const T* sourceCol = &data[jLocal*thisLocalHeight];
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    destCol[iLocal*c] = sourceCol[iLocal];
            }
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MR,MC,Int> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MC]");
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
            this->SetRowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return *this;
    }

    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int rowShift = this->RowShift();
        const Int rowShiftOfA = A.RowShift();
        const Int rowOffset = (rowShift-rowShiftOfA) / r;

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();

        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* ACol = &ABuffer[(rowOffset+jLocal*c)*ALDim];
            T* thisCol = &thisBuffer[jLocal*thisLDim];
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
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        // We will SendRecv A[*,VC] within our process column to fix alignments.
        const Int sendRow = (row+r+(rowAlignment%r)-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-(rowAlignment%r)) % r;
        const Int sendRank = sendRow + r*col;

        const Int sendRowShift = Shift( sendRank, rowAlignment, p );
        const Int sendRowOffset = (sendRowShift-rowShiftOfA) / r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfSend = Length(width,sendRowShift,p);

        const Int sendSize = height * localWidthOfSend;
        const Int recvSize = height * localWidth;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidthOfSend; ++jLocal )
        {
            const T* ACol = &ABuffer[(sendRowOffset+jLocal*c)*ALDim];
            T* sendBufferCol = &sendBuffer[jLocal*height];
            MemCopy( sendBufferCol, ACol, height );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendRow, 0,
          recvBuffer, recvSize, recvRow, mpi::ANY_TAG, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufferCol = &recvBuffer[jLocal*height];
            T* thisCol = &thisBuffer[jLocal*thisLDim];
            MemCopy( thisCol, recvBufferCol, height );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MC,MR,Int> > A_MC_MR
    ( new DistMatrix<T,MC,MR,Int>(g) );
    *A_MC_MR = A;

    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR,Int>(g) );
    *A_STAR_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,VC]");
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
            if( g.InGrid() )
                this->rowShift_ = A.RowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return *this;
    }

    if( this->RowAlignment() == A.RowAlignment() )
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

        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const Int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfA = A.LocalWidth();

        const Int sendSize = height * localWidthOfA;
        const Int recvSize = height * localWidth;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
        {
            const T* ACol = &ABuffer[jLocal*ALDim];
            T* sendBufferCol = &sendBuffer[jLocal*height];
            MemCopy( sendBufferCol, ACol, height );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.VCComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufferCol = &recvBuffer[jLocal*height];
            T* thisCol = &thisBuffer[jLocal*thisLDim];
            MemCopy( thisCol, recvBufferCol, height );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MR,MC,Int> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return *this;
    }
    
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

    this->auxMemory_.Require( sendSize + recvSize );

    T* buffer = this->auxMemory_.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack
    const T* ABuffer = A.LockedBuffer();
    const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
    {
        const T* ACol = &ABuffer[jLocal*ALDim];
        T* sendBufferCol = &sendBuffer[jLocal*height];
        MemCopy( sendBufferCol, ACol, height );
    }

    // Communicate
    mpi::SendRecv
    ( sendBuffer, sendSize, sendRankCM, 0,
      recvBuffer, recvSize, recvRankCM, mpi::ANY_TAG, g.VCComm() );

    // Unpack
    T* thisBuffer = this->Buffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* recvBufferCol = &recvBuffer[jLocal*height];
        T* thisCol = &thisBuffer[jLocal*thisLDim];
        MemCopy( thisCol, recvBufferCol, height );
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return *this;
    }

    const Int p = g.Size();
    const Int rowShift = this->RowShift();

    const Int height = this->Height();
    const Int localWidth = this->LocalWidth();

    T* thisBuffer = this->Buffer();
    const Int thisLDim = this->LDim();
    const T* ABuffer = A.LockedBuffer();
    const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* ACol = &ABuffer[(rowShift+jLocal*p)*ALDim];
        T* thisCol = &thisBuffer[jLocal*thisLDim];
        MemCopy( thisCol, ACol, height );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::SumScatterFrom( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SumScatterFrom( [* ,MC] )");
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
            this->SetRowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int rowAlignment = this->RowAlignment();
        const Int rowShiftOfA = A.RowShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalWidth = MaxLength( width, p );

        const Int recvSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);
        const Int sendSize = c*recvSize;

        this->auxMemory_.Require( sendSize );
        T* buffer = this->auxMemory_.Buffer();

        // Pack
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[k*recvSize];

            const Int thisRank = myRow+k*r;
            const Int thisRowShift = Shift_( thisRank, rowAlignment, p );
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = Length_( width, thisRowShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ABuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*height];
                MemCopy( dataCol, ACol, height );
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.RowComm() );

        // Unpack our received data
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* bufferCol = &buffer[jLocal*height];
            T* thisCol = &thisBuffer[jLocal*thisLDim];
            MemCopy( thisCol, bufferCol, height );
        }
        this->auxMemory_.Release();
    }
    else
    {
        throw std::logic_error
        ("Unaligned [* ,VC]::ReduceScatterFrom( [* ,MC] ) not yet implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SumScatterUpdate( [* ,MC] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int rowAlignment = this->RowAlignment();
        const Int rowShiftOfA = A.RowShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalWidth = MaxLength( width, p );

        const Int recvSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);
        const Int sendSize = c*recvSize;

        this->auxMemory_.Require( sendSize );
        T* buffer = this->auxMemory_.Buffer();

        // Pack
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[k*recvSize];

            const Int thisRank = myRow+k*r;
            const Int thisRowShift = Shift_( thisRank, rowAlignment, p );
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = Length_( width, thisRowShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ABuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*height];
                MemCopy( dataCol, ACol, height );
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.RowComm() );

        // Unpack our received data
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* bufferCol = &buffer[jLocal*height];
            T* thisCol = &thisBuffer[jLocal*thisLDim];
            for( Int i=0; i<height; ++i )
                thisCol[i] += alpha*bufferCol[i];
        }
        this->auxMemory_.Release();
    }
    else
    {
        throw std::logic_error
        ("Unaligned [* ,VC]::ReduceScatterUpdate( [* ,MC] ) not implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Routines which explicitly work in the complex plane
//

template<typename T,typename Int>
BASE(T)
DistMatrix<T,STAR,VC,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::GetRealPart");
    this->AssertValidEntry( i, j );
#endif
    typedef BASE(T) R;

    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocalRealPart(i,jLoc);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
BASE(T)
DistMatrix<T,STAR,VC,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::GetImagPart");
    this->AssertValidEntry( i, j );
#endif
    typedef BASE(T) R;

    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocalImagPart(i,jLoc);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.Size();
        this->SetLocalRealPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.Size();
        this->SetLocalImagPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.Size();
        this->UpdateLocalRealPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,VC,Int>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.Size();
        this->UpdateLocalImagPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template class DistMatrix<int,STAR,VC,int>;
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,MC,  MR,  int>& A );
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,MC,  STAR,int>& A );
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,MD,  STAR,int>& A );
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,MR,  MC,  int>& A );
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,MR,  STAR,int>& A );
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,STAR,MC,  int>& A );
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,STAR,MD,  int>& A );
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,STAR,MR,  int>& A );
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,STAR,STAR,int>& A );
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,STAR,VR,  int>& A );
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,VC,  STAR,int>& A );
template DistMatrix<int,STAR,VC,int>::DistMatrix( const DistMatrix<int,VR,  STAR,int>& A );

#ifndef DISABLE_FLOAT
template class DistMatrix<float,STAR,VC,int>;
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,MC,  MR,  int>& A );
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,MC,  STAR,int>& A );
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,MD,  STAR,int>& A );
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,MR,  MC,  int>& A );
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,MR,  STAR,int>& A );
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,STAR,MC,  int>& A );
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,STAR,MD,  int>& A );
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,STAR,MR,  int>& A );
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,STAR,STAR,int>& A );
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,STAR,VR,  int>& A );
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,VC,  STAR,int>& A );
template DistMatrix<float,STAR,VC,int>::DistMatrix( const DistMatrix<float,VR,  STAR,int>& A );
#endif // ifndef DISABLE_FLOAT

template class DistMatrix<double,STAR,VC,int>;
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,MC,  MR,  int>& A );
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,MC,  STAR,int>& A );
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,MD,  STAR,int>& A );
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,MR,  MC,  int>& A );
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,MR,  STAR,int>& A );
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,STAR,MC,  int>& A );
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,STAR,MD,  int>& A );
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,STAR,MR,  int>& A );
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,STAR,STAR,int>& A );
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,STAR,VR,  int>& A );
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,VC,  STAR,int>& A );
template DistMatrix<double,STAR,VC,int>::DistMatrix( const DistMatrix<double,VR,  STAR,int>& A );

#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
template class DistMatrix<Complex<float>,STAR,VC,int>;
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,MC,  MR,  int>& A );
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,MC,  STAR,int>& A );
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,MD,  STAR,int>& A );
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,MR,  MC,  int>& A );
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,MR,  STAR,int>& A );
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,MC,  int>& A );
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,MD,  int>& A );
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,MR,  int>& A );
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,STAR,int>& A );
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,VR,  int>& A );
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,VC,  STAR,int>& A );
template DistMatrix<Complex<float>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<float>,VR,  STAR,int>& A );
#endif // ifndef DISABLE_FLOAT
template class DistMatrix<Complex<double>,STAR,VC,int>;
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,MC,  MR,  int>& A );
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,MC,  STAR,int>& A );
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,MD,  STAR,int>& A );
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,MR,  MC,  int>& A );
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,MR,  STAR,int>& A );
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,MC,  int>& A );
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,MD,  int>& A );
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,MR,  int>& A );
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,STAR,int>& A );
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,VR,  int>& A );
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,VC,  STAR,int>& A );
template DistMatrix<Complex<double>,STAR,VC,int>::DistMatrix( const DistMatrix<Complex<double>,VR,  STAR,int>& A );
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
