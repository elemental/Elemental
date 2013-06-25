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
DistMatrix<T,STAR,STAR,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,STAR,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->ResizeTo(height,width); }

template<typename T,typename Int>
DistMatrix<T,STAR,STAR,Int>::DistMatrix
( Int height, Int width, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->ResizeTo(height,width,ldim); }

template<typename T,typename Int>
DistMatrix<T,STAR,STAR,Int>::DistMatrix
( Int height, Int width, const T* buffer, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->LockedAttach(height,width,buffer,ldim,g); }

template<typename T,typename Int>
DistMatrix<T,STAR,STAR,Int>::DistMatrix
( Int height, Int width, T* buffer, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->Attach(height,width,buffer,ldim,g); }

template<typename T,typename Int>
DistMatrix<T,STAR,STAR,Int>::DistMatrix( const DistMatrix<T,STAR,STAR,Int>& A )
: AbstractDistMatrix<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,* ] with itself");
}

template<typename T,typename Int>
template<Distribution U,Distribution V>
DistMatrix<T,STAR,STAR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,* ]::DistMatrix");
#endif
    if( STAR != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,STAR,STAR,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,* ] with itself");
}

template<typename T,typename Int>
DistMatrix<T,STAR,STAR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
elem::DistData<Int>
DistMatrix<T,STAR,STAR,Int>::DistData() const
{
    elem::DistData<Int> data;
    data.colDist = STAR;
    data.rowDist = STAR;
    data.colAlignment = 0;
    data.rowAlignment = 0;
    data.root = 0;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T,typename Int>
Int
DistMatrix<T,STAR,STAR,Int>::ColStride() const
{ return 1; }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,STAR,Int>::RowStride() const
{ return 1; }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,STAR,Int>::ColRank() const
{ return 0; }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,STAR,Int>::RowRank() const
{ return 0; }

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::Attach
( Int height, Int width, 
  T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::Attach");
#endif
    this->Empty();
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->viewType_ = VIEW;
    if( this->Participating() )
        this->matrix_.Attach_( height, width, buffer, ldim );
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::LockedAttach
( Int height, Int width, 
  const T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::LockedAttach");
#endif
    this->Empty();
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->viewType_ = LOCKED_VIEW;
    if( this->Participating() )
        this->matrix_.LockedAttach_( height, width, buffer, ldim );
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_( height, width );
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_( height, width, ldim );
}

template<typename T,typename Int>
T
DistMatrix<T,STAR,STAR,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const Int viewingSize = mpi::CommSize( g.ViewingComm() );
    const Int owningSize = mpi::GroupSize( g.OwningGroup() );
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        return this->GetLocal(i,j);
    }
    else
    {
        // Have the root broadcast its data
        T u;
        if( g.VCRank() == 0 )
            u = this->GetLocal(i,j);
        mpi::Broadcast( &u, 1, g.VCToViewingMap(0), g.ViewingComm() );
        return u;
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->SetLocal(i,j,u);
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->UpdateLocal(i,j,u);
}

//
// Utility functions, e.g., operator=
//

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int r = g.Height();
    const Int c = g.Width(); 
    const Int p = g.Size();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalHeight = MaxLength(height,r);
    const Int maxLocalWidth = MaxLength(width,c);

    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
    T* buffer = this->auxMemory_.Require( (p+1)*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        MemCopy
        ( &sendBuf[jLoc*localHeightOfA], 
          &ABuf[jLoc*ALDim], localHeightOfA );

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.VCComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const Int colAlignmentOfA = A.ColAlignment();
    const Int rowAlignmentOfA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int l=0; l<c; ++l )
    {
        const Int rowShift = Shift_( l, rowAlignmentOfA, c );
        const Int localWidth = Length_( width, rowShift, c );
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[(k+l*r)*portionSize];
            const Int colShift = Shift_( k, colAlignmentOfA, r );
            const Int localHeight = Length_( height, colShift, r );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[colShift+(rowShift+jLoc*c)*thisLDim];
                const T* sourceCol = &data[jLoc*localHeight];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc*r] = sourceCol[iLoc];
            }
        }
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int r = g.Height();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int maxLocalHeight = MaxLength(height,r);

    const Int portionSize = mpi::Pad( maxLocalHeight*width );
    T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
        MemCopy
        ( &sendBuf[j*localHeightOfA], &ABuf[j*ALDim], localHeightOfA );

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.ColComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const Int colAlignmentOfA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<r; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int colShift = Shift_( k, colAlignmentOfA, r );
        const Int localHeight = Length_( height, colShift, r );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisBuf[colShift+j*thisLDim];
            const T* sourceCol = &data[j*localHeight];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc*r] = sourceCol[iLoc];
        }
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int c = g.Width();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalWidth = MaxLength(width,c);

    const Int portionSize = mpi::Pad( height*maxLocalWidth );
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
        MemCopy( &sendBuf[jLoc*height], &ABuf[jLoc*ALDim], height );

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.RowComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const Int rowAlignmentOfA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<c; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int rowShift = Shift_( k, rowAlignmentOfA, c );
        const Int localWidth = Length_( width, rowShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuf[(rowShift+jLoc*c)*thisLDim],
              &data[jLoc*height], height );
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int p = g.Size();
    const Int lcm = g.LCM();
    const Int ownerPath = A.diagPath_;
    const Int ownerPathRank = A.colAlignment_;

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = A.LocalHeight();
    const Int maxLocalHeight = MaxLength( height, lcm );
    const Int portionSize = mpi::Pad( maxLocalHeight*width );

    // Since a MD communicator has not been implemented, we will take
    // the suboptimal route of 'rounding up' everyone's contribution over 
    // the VC communicator.
    T* buffer = this->auxMemory_.Require( (p+1)*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    if( A.Participating() )
    {
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
            MemCopy( &sendBuf[j*localHeight], &ABuf[j*ALDim], localHeight );
    }

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.VCComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        if( g.DiagPath( k ) == ownerPath )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisPathRank = g.DiagPathRank( k );
            const Int thisColShift = Shift_( thisPathRank, ownerPathRank, lcm );
            const Int thisLocalHeight = Length_( height, thisColShift, lcm );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &thisBuf[thisColShift+j*thisLDim];
                const T* sourceCol = &data[j*thisLocalHeight];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc*lcm] = sourceCol[iLoc];
            }
        }
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int p = g.Size();
    const Int lcm = g.LCM();
    const Int ownerPath = A.diagPath_;
    const Int ownerPathRank = A.rowAlignment_;

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localWidth = A.LocalWidth();
    const Int maxLocalWidth = MaxLength( width, lcm );
    const Int portionSize = mpi::Pad( height*maxLocalWidth );

    // Since a MD communicator has not been implemented, we will take
    // the suboptimal route of 'rounding up' everyone's contribution over 
    // the VC communicator.
    T* buffer = this->auxMemory_.Require( (p+1)*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    if( A.Participating() )
    {
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy( &sendBuf[jLoc*height], &ABuf[jLoc*ALDim], height );
    }

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.VCComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        if( g.DiagPath( k ) == ownerPath )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisPathRank = g.DiagPathRank( k );
            const Int thisRowShift = Shift_( thisPathRank, ownerPathRank, lcm );
            const Int thisLocalWidth = Length_( width, thisRowShift, lcm );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                MemCopy
                ( &thisBuf[(thisRowShift+jLoc*lcm)*thisLDim], 
                  &data[jLoc*height], height );
        }
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalHeight = MaxLength(height,c);
    const Int maxLocalWidth = MaxLength(width,r);

    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
    T* buffer = this->auxMemory_.Require( (p+1)*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        MemCopy
        ( &sendBuf[jLoc*localHeightOfA], 
          &ABuf[jLoc*ALDim], localHeightOfA );

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.VRComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const Int colAlignmentOfA = A.ColAlignment();
    const Int rowAlignmentOfA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int l=0; l<r; ++l )
    {
        const Int rowShift = Shift_( l, rowAlignmentOfA, r );
        const Int localWidth = Length_( width, rowShift, r );
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[(k+l*c)*portionSize];
            const Int colShift = Shift_( k, colAlignmentOfA, c );
            const Int localHeight = Length_( height, colShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[colShift+(rowShift+jLoc*r)*thisLDim];
                const T* sourceCol = &data[jLoc*localHeight];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc*c] = sourceCol[iLoc];
            }
        }
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int c = g.Width();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int maxLocalHeight = MaxLength(height,c);

    const Int portionSize = mpi::Pad( maxLocalHeight*width );
    T* buffer = this->auxMemory_.Require( (c+1)*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
        MemCopy
        ( &sendBuf[j*localHeightOfA], &ABuf[j*ALDim], localHeightOfA );

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.RowComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const Int colAlignmentOfA = A.ColAlignment();
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
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisBuf[colShift+j*thisLDim];
            const T* sourceCol = &data[j*localHeight];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc*c] = sourceCol[iLoc];
        }
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [* ,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int r = g.Height();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalWidth = MaxLength(width,r);

    const Int portionSize = mpi::Pad( height*maxLocalWidth );
    T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        MemCopy( &sendBuf[jLoc*height], &ABuf[jLoc*ALDim], height );

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.ColComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const Int rowAlignmentOfA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<r; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int rowShift = Shift_( k, rowAlignmentOfA, r );
        const Int localWidth = Length_( width, rowShift, r );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuf[(rowShift+jLoc*r)*thisLDim], 
              &data[jLoc*height], height );
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int p = g.Size();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int maxLocalHeight = MaxLength(height,p);

    const Int portionSize = mpi::Pad( maxLocalHeight*width );
    T* buffer = this->auxMemory_.Require( (p+1)*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
        MemCopy
        ( &sendBuf[j*localHeightOfA], &ABuf[j*ALDim], localHeightOfA );

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.VCComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const Int colAlignmentOfA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int colShift = Shift_( k, colAlignmentOfA, p );
        const Int localHeight = Length_( height, colShift, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisBuf[colShift+j*thisLDim];
            const T* sourceCol = &data[j*localHeight];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc*p] = sourceCol[iLoc];
        }
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int p = g.Size();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalWidth = MaxLength(width,p);

    const Int portionSize = mpi::Pad( height*maxLocalWidth );
    T* buffer = this->auxMemory_.Require( (p+1)*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        MemCopy( &sendBuf[jLoc*height], &ABuf[jLoc*ALDim], height );

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.VCComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const Int rowAlignmentOfA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int rowShift = Shift_( k, rowAlignmentOfA, p );
        const Int localWidth = Length_( width, rowShift, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuf[(rowShift+jLoc*p)*thisLDim], 
              &data[jLoc*height], height );
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int p = g.Size();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int maxLocalHeight = MaxLength(height,p);

    const Int portionSize = mpi::Pad( maxLocalHeight*width );
    T* buffer = this->auxMemory_.Require( (p+1)*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
        MemCopy
        ( &sendBuf[j*localHeightOfA], &ABuf[j*ALDim], localHeightOfA );

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.VRComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const Int colAlignmentOfA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int colShift = Shift_( k, colAlignmentOfA, p );
        const Int localHeight = Length_( height, colShift, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisBuf[colShift+j*thisLDim];
            const T* sourceCol = &data[j*localHeight];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc*p] = sourceCol[iLoc];
        }
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int p = g.Size();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalWidth = MaxLength(width,p);

    const Int portionSize = mpi::Pad( height*maxLocalWidth );
    T* buffer = this->auxMemory_.Require( (p+1)*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        MemCopy( &sendBuf[jLoc*height], &ABuf[jLoc*ALDim], height );

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.VRComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const Int rowAlignmentOfA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int rowShift = Shift_( k, rowAlignmentOfA, p );
        const Int localWidth = Length_( width, rowShift, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuf[(rowShift+jLoc*p)*thisLDim], 
              &data[jLoc*height], height );
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [* ,* ]");
    this->AssertNotLocked();
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    if( this->Grid() == A.Grid() )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
        // TODO: Remember why I wrote this...
        if( !mpi::CongruentComms( A.Grid().ViewingComm(),
                                  this->Grid().ViewingComm() ) )
            throw std::logic_error
            ("Redistributing between nonmatching grids currently requires"
             " the viewing communicators to match.");

        // Compute and allocate the amount of required memory
        Int requiredMemory = 0;
        if( A.Grid().VCRank() == 0 )
            requiredMemory += A.Height()*A.Width();
        if( this->Participating() )
            requiredMemory += A.Height()*A.Width();
        T* buffer = this->auxMemory_.Require( requiredMemory );
        Int offset = 0;
        T* sendBuf = &buffer[offset];
        if( A.Grid().VCRank() == 0 )
            offset += A.Height()*A.Width();
        T* bcastBuffer = &buffer[offset];

        // Send from the root of A to the root of this matrix's grid
        mpi::Request sendRequest;
        if( A.Grid().VCRank() == 0 )
        {
            for( Int j=0; j<A.Width(); ++j ) 
                for( Int i=0; i<A.Height(); ++i )
                    sendBuf[i+j*A.Height()] = A.GetLocal(i,j);
            const Int recvViewingRank = this->Grid().VCToViewingMap(0);
            mpi::ISend
            ( sendBuf, A.Height()*A.Width(), recvViewingRank, 0,
              this->Grid().ViewingComm(), sendRequest );
        }

        // Receive on the root of this matrix's grid and then broadcast
        // over this matrix's owning communicator
        if( this->Participating() )
        {
            if( this->Grid().VCRank() == 0 )
            {
                const Int sendViewingRank = A.Grid().VCToViewingMap(0);
                mpi::Recv
                ( bcastBuffer, A.Height()*A.Width(), sendViewingRank, 0,
                  this->Grid().ViewingComm() );
            }

            mpi::Broadcast
            ( bcastBuffer, A.Height()*A.Width(), 0, this->Grid().VCComm() );

            for( Int j=0; j<A.Width(); ++j )
                for( Int i=0; i<A.Height(); ++i )
                    this->SetLocal(i,j,bcastBuffer[i+j*A.Height()]);
        }

        if( A.Grid().VCRank() == 0 )
            mpi::Wait( sendRequest );
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,CIRC,CIRC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ] = [o ,o ]");
    this->AssertNotLocked();
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Grid& g = A.Grid();
    const int m = A.Height(); 
    const int n = A.Width();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    if( this->Participating() )
    {
        const int pkgSize = mpi::Pad( m*n );
        T* commBuffer = this->auxMemory_.Require( pkgSize );

        if( A.Participating() )
        {
            // Pack            
            const int ALDim = A.LDim();
            const T* ABuf = A.LockedBuffer();
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    commBuffer[i+j*m] = ABuf[i+j*ALDim];
        }

        // Broadcast from the process that packed
        mpi::Broadcast( commBuffer, pkgSize, A.Root(), g.VCComm() );

        // Unpack
        T* buffer = this->Buffer();
        const int ldim = this->LDim();
        for( int j=0; j<n; ++j )
            for( int i=0; i<m; ++i )
                buffer[i+j*ldim] = commBuffer[i+j*m];        
    }

    return *this;
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::SumOverCol()
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::SumOverCol");
    this->AssertNotLocked();
#endif
    const Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int localSize = localHeight*localWidth;
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
        MemCopy
        ( &sendBuf[jLoc*localHeight], 
          &thisBuf[jLoc*thisLDim], localHeight );

    // Sum
    mpi::AllReduce( sendBuf, recvBuf, localSize, mpi::SUM, g.ColComm() );

    // Unpack
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        MemCopy
        ( &thisBuf[jLoc*thisLDim], 
          &recvBuf[jLoc*localHeight], localHeight );
    this->auxMemory_.Release();
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::SumOverRow()
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::SumOverRow");
    this->AssertNotLocked();
#endif
    const Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int localSize = localHeight*localWidth;
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
        MemCopy
        ( &sendBuf[jLoc*localHeight], 
          &thisBuf[jLoc*thisLDim], localHeight );

    // Sum
    mpi::AllReduce( sendBuf, recvBuf, localSize, mpi::SUM, g.RowComm() );

    // Unpack
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        MemCopy
        ( &thisBuf[jLoc*thisLDim], 
          &recvBuf[jLoc*localHeight], localHeight );
    this->auxMemory_.Release();
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::SumOverGrid()
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::SumOverGrid");
    this->AssertNotLocked();
#endif
    const Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int localSize = localHeight*localWidth;
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
        MemCopy
        ( &sendBuf[jLoc*localHeight], 
          &thisBuf[jLoc*thisLDim], localHeight );

    // Sum
    mpi::AllReduce( sendBuf, recvBuf, localSize, mpi::SUM, g.VCComm() );

    // Unpack
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        MemCopy
        ( &thisBuf[jLoc*thisLDim], &recvBuf[jLoc*localHeight], localHeight );
    this->auxMemory_.Release();
}

//
// Routines which explicitly work in the complex plane
//

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->SetLocalRealPart(i,j,u);
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");
    if( this->Participating() )
        this->SetLocalImagPart(i,j,u);
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->UpdateLocalRealPart(i,j,u);
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,STAR,Int>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,* ]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");
    if( this->Participating() )
        this->UpdateLocalImagPart(i,j,u);
}

#define PROTO(T) \
  template class DistMatrix<T,STAR,STAR,int>
#define COPY(T,CD,RD) \
  template DistMatrix<T,STAR,STAR,int>::DistMatrix( \
    const DistMatrix<T,CD,RD,int>& A )
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
