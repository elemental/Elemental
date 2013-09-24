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
DistMatrix<T,STAR,STAR>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ }

template<typename T>
DistMatrix<T,STAR,STAR>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->ResizeTo(height,width); }

template<typename T>
DistMatrix<T,STAR,STAR>::DistMatrix
( Int height, Int width, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->ResizeTo(height,width,ldim); }

template<typename T>
DistMatrix<T,STAR,STAR>::DistMatrix
( Int height, Int width, const T* buffer, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->LockedAttach(height,width,buffer,ldim,g); }

template<typename T>
DistMatrix<T,STAR,STAR>::DistMatrix
( Int height, Int width, T* buffer, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Attach(height,width,buffer,ldim,g); }

template<typename T>
DistMatrix<T,STAR,STAR>::DistMatrix( const DistMatrix<T,STAR,STAR>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[* ,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [* ,* ] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DistMatrix<T,STAR,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[* ,* ]::DistMatrix");
#endif
    if( STAR != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,STAR,STAR>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [* ,* ] with itself");
}

template<typename T>
DistMatrix<T,STAR,STAR>::DistMatrix( DistMatrix<T,STAR,STAR>&& A )
: AbstractDistMatrix<T>(std::move(A))
{ }

template<typename T>
DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( DistMatrix<T,STAR,STAR>&& A )
{
    AbstractDistMatrix<T>::operator=( std::move(A) );
    return *this;
}

template<typename T>
DistMatrix<T,STAR,STAR>::~DistMatrix()
{ }

template<typename T>
elem::DistData
DistMatrix<T,STAR,STAR>::DistData() const
{
    elem::DistData data;
    data.colDist = STAR;
    data.rowDist = STAR;
    data.colAlignment = 0;
    data.rowAlignment = 0;
    data.root = 0;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T>
Int
DistMatrix<T,STAR,STAR>::ColStride() const
{ return 1; }

template<typename T>
Int
DistMatrix<T,STAR,STAR>::RowStride() const
{ return 1; }

template<typename T>
Int
DistMatrix<T,STAR,STAR>::ColRank() const
{ return 0; }

template<typename T>
Int
DistMatrix<T,STAR,STAR>::RowRank() const
{ return 0; }

template<typename T>
void
DistMatrix<T,STAR,STAR>::Attach
( Int height, Int width, 
  T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::Attach");
#endif
    this->Empty();
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->viewType_ = VIEW;
    if( this->Participating() )
        this->matrix_.Attach_( height, width, buffer, ldim );
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::LockedAttach
( Int height, Int width, 
  const T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::LockedAttach");
#endif
    this->Empty();
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->viewType_ = LOCKED_VIEW;
    if( this->Participating() )
        this->matrix_.LockedAttach_( height, width, buffer, ldim );
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_( height, width );
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_( height, width, ldim );
}

template<typename T>
T
DistMatrix<T,STAR,STAR>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::Get");
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
        mpi::Broadcast( u, g.VCToViewingMap(0), g.ViewingComm() );
        return u;
    }
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->SetLocal(i,j,u);
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->SetLocalRealPart(i,j,u);
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    if( this->Participating() )
        this->SetLocalImagPart(i,j,u);
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->UpdateLocal(i,j,u);
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->UpdateLocalRealPart(i,j,u);
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    if( this->Participating() )
        this->UpdateLocalImagPart(i,j,u);
}

template<typename T>
template<typename S,class Function>
void
DistMatrix<T,STAR,STAR>::GetDiagonalHelper
( DistMatrix<S,STAR,STAR>& d, Int offset, Function func ) const
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::GetDiagonalHelper");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
        d.SetGrid( g );
    const Int diagLength = this->DiagonalLength(offset);
    d.ResizeTo( diagLength, 1 );
    if( !d.Participating() )
        return;

    const Int iStart = ( offset>=0 ? 0      : -offset );
    const Int jStart = ( offset>=0 ? offset : 0       );
    S* dBuf = d.Buffer();  
    const T* buffer = this->LockedBuffer();
    const Int ldim = this->LDim();
    for( Int k=0; k<diagLength; ++k )
        func( dBuf[k], buffer[(iStart)+k+(jStart+k)*ldim] );
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::GetDiagonal
( DistMatrix<T,STAR,STAR>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::GetDiagonal");
#endif
    this->GetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::GetRealPartOfDiagonal
( DistMatrix<BASE(T),STAR,STAR>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::GetRealPartOfDiagonal");
#endif
    this->GetDiagonalHelper
    ( d, offset, []( BASE(T)& alpha, T beta ) { alpha = RealPart(beta); } );
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::GetImagPartOfDiagonal
( DistMatrix<BASE(T),STAR,STAR>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::GetImagPartOfDiagonal");
#endif
    this->GetDiagonalHelper
    ( d, offset, []( BASE(T)& alpha, T beta ) { alpha = ImagPart(beta); } );
}

template<typename T>
DistMatrix<T,STAR,STAR>
DistMatrix<T,STAR,STAR>::GetDiagonal( Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::GetDiagonal");
#endif
    DistMatrix<T,STAR,STAR> d( this->Grid() );
    GetDiagonal( d, offset );
    return d;
}

template<typename T>
DistMatrix<BASE(T),STAR,STAR>
DistMatrix<T,STAR,STAR>::GetRealPartOfDiagonal( Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::GetRealPartOfDiagonal");
#endif
    DistMatrix<BASE(T),STAR,STAR> d( this->Grid() );
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
DistMatrix<BASE(T),STAR,STAR>
DistMatrix<T,STAR,STAR>::GetImagPartOfDiagonal( Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::GetImagPartOfDiagonal");
#endif
    DistMatrix<BASE(T),STAR,STAR> d( this->Grid() );
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
template<typename S,class Function>
void
DistMatrix<T,STAR,STAR>::SetDiagonalHelper
( const DistMatrix<S,STAR,STAR>& d, Int offset, Function func )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::SetDiagonalHelper");
    this->AssertSameGrid( d.Grid() );
    if( d.Height() != 1 && d.Width() != 1 )
        LogicError("d must be a vector");
#endif
    const Int diagLength = this->DiagonalLength(offset);
    Int dStride;
    if( d.Width() == 1 )
    {
#ifndef RELEASE
        if( d.Height() != diagLength )
            LogicError("d is not the right length");
#endif
        dStride = 1;
    }
    else
    {
#ifndef RELEASE
        if( d.Width() != diagLength )
            LogicError("d is not the right length");
#endif
        dStride = d.LDim();
    }
    const Int iStart = ( offset>=0 ? 0      : -offset );
    const Int jStart = ( offset>=0 ? offset : 0       );
    const S* dBuf = d.LockedBuffer();
    T* buffer = this->Buffer();
    const Int ldim = this->LDim();
    for( Int k=0; k<diagLength; ++k )
        func( buffer[(iStart+k)+(jStart+k)*ldim], dBuf[k*dStride] );
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::SetDiagonal
( const DistMatrix<T,STAR,STAR>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::SetDiagonal");
#endif
    this->SetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::SetRealPartOfDiagonal
( const DistMatrix<BASE(T),STAR,STAR>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::SetRealPartOfDiagonal");
#endif
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, BASE(T) beta ) { elem::SetRealPart(alpha,beta); } );
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::SetImagPartOfDiagonal
( const DistMatrix<BASE(T),STAR,STAR>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::SetImagPartOfDiagonal");
#endif
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, BASE(T) beta ) { elem::SetImagPart(alpha,beta); } );
}

//
// Utility functions, e.g., operator=
//

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
    PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int l=0; l<c; ++l )
    {
        const Int rowShift = Shift_( l, rowAlignmentOfA, c );
        const Int localWidth = Length_( width, rowShift, c );
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[(k+l*r)*portionSize];
            const Int colShift = Shift_( k, colAlignmentOfA, r );
            const Int localHeight = Length_( height, colShift, r );
            INNER_PARALLEL_FOR
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

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
    PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int k=0; k<r; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int colShift = Shift_( k, colAlignmentOfA, r );
        const Int localHeight = Length_( height, colShift, r );
        INNER_PARALLEL_FOR
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

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
    PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int k=0; k<c; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int rowShift = Shift_( k, rowAlignmentOfA, c );
        const Int localWidth = Length_( width, rowShift, c );
        INNER_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuf[(rowShift+jLoc*c)*thisLDim],
              &data[jLoc*height], height );
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
        PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int k=0; k<p; ++k )
    {
        if( g.DiagPath( k ) == ownerPath )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisPathRank = g.DiagPathRank( k );
            const Int thisColShift = Shift_( thisPathRank, ownerPathRank, lcm );
            const Int thisLocalHeight = Length_( height, thisColShift, lcm );
            INNER_PARALLEL_FOR
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

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
        PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int k=0; k<p; ++k )
    {
        if( g.DiagPath( k ) == ownerPath )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisPathRank = g.DiagPathRank( k );
            const Int thisRowShift = Shift_( thisPathRank, ownerPathRank, lcm );
            const Int thisLocalWidth = Length_( width, thisRowShift, lcm );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                MemCopy
                ( &thisBuf[(thisRowShift+jLoc*lcm)*thisLDim], 
                  &data[jLoc*height], height );
        }
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
    PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int l=0; l<r; ++l )
    {
        const Int rowShift = Shift_( l, rowAlignmentOfA, r );
        const Int localWidth = Length_( width, rowShift, r );
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[(k+l*c)*portionSize];
            const Int colShift = Shift_( k, colAlignmentOfA, c );
            const Int localHeight = Length_( height, colShift, c );
            INNER_PARALLEL_FOR
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

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
    PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int k=0; k<c; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int colShift = Shift_( k, colAlignmentOfA, c );
        const Int localHeight = Length_( height, colShift, c );
        INNER_PARALLEL_FOR
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

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [* ,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
    PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int k=0; k<r; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int rowShift = Shift_( k, rowAlignmentOfA, r );
        const Int localWidth = Length_( width, rowShift, r );
        INNER_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuf[(rowShift+jLoc*r)*thisLDim], 
              &data[jLoc*height], height );
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,VC,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
    PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int k=0; k<p; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int colShift = Shift_( k, colAlignmentOfA, p );
        const Int localHeight = Length_( height, colShift, p );
        INNER_PARALLEL_FOR
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

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,VC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
    PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int k=0; k<p; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int rowShift = Shift_( k, rowAlignmentOfA, p );
        const Int localWidth = Length_( width, rowShift, p );
        INNER_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuf[(rowShift+jLoc*p)*thisLDim], 
              &data[jLoc*height], height );
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,VR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
    PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int k=0; k<p; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int colShift = Shift_( k, colAlignmentOfA, p );
        const Int localHeight = Length_( height, colShift, p );
        INNER_PARALLEL_FOR
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

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,VR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
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
    PARALLEL_FOR
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
    OUTER_PARALLEL_FOR
    for( Int k=0; k<p; ++k )
    {
        const T* data = &recvBuf[k*portionSize];
        const Int rowShift = Shift_( k, rowAlignmentOfA, p );
        const Int localWidth = Length_( width, rowShift, p );
        INNER_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuf[(rowShift+jLoc*p)*thisLDim], 
              &data[jLoc*height], height );
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [* ,* ]");
    this->AssertNotLocked();
#endif
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
            LogicError
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
            ( sendBuf, A.Height()*A.Width(), recvViewingRank,
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
                ( bcastBuffer, A.Height()*A.Width(), sendViewingRank,
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

template<typename T>
const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ] = [o ,o ]");
    this->AssertNotLocked();
#endif
    const Grid& g = A.Grid();
    const Int m = A.Height(); 
    const Int n = A.Width();
    this->ResizeTo( A.Height(), A.Width() );

    if( this->Participating() )
    {
        const Int pkgSize = mpi::Pad( m*n );
        T* commBuffer = this->auxMemory_.Require( pkgSize );

        if( A.Participating() )
        {
            // Pack            
            const Int ALDim = A.LDim();
            const T* ABuf = A.LockedBuffer();
            for( Int j=0; j<n; ++j )
                for( Int i=0; i<m; ++i )
                    commBuffer[i+j*m] = ABuf[i+j*ALDim];
        }

        // Broadcast from the process that packed
        mpi::Broadcast( commBuffer, pkgSize, A.Root(), g.VCComm() );

        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                buffer[i+j*ldim] = commBuffer[i+j*m];        
    }

    return *this;
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::SumOverCol()
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::SumOverCol");
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
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        MemCopy
        ( &sendBuf[jLoc*localHeight], 
          &thisBuf[jLoc*thisLDim], localHeight );

    // Sum
    mpi::AllReduce( sendBuf, recvBuf, localSize, g.ColComm() );

    // Unpack
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        MemCopy
        ( &thisBuf[jLoc*thisLDim], 
          &recvBuf[jLoc*localHeight], localHeight );
    this->auxMemory_.Release();
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::SumOverRow()
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::SumOverRow");
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
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        MemCopy
        ( &sendBuf[jLoc*localHeight], 
          &thisBuf[jLoc*thisLDim], localHeight );

    // Sum
    mpi::AllReduce( sendBuf, recvBuf, localSize, g.RowComm() );

    // Unpack
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        MemCopy
        ( &thisBuf[jLoc*thisLDim], 
          &recvBuf[jLoc*localHeight], localHeight );
    this->auxMemory_.Release();
}

template<typename T>
void
DistMatrix<T,STAR,STAR>::SumOverGrid()
{
#ifndef RELEASE
    CallStackEntry cse("[* ,* ]::SumOverGrid");
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
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        MemCopy
        ( &sendBuf[jLoc*localHeight], 
          &thisBuf[jLoc*thisLDim], localHeight );

    // Sum
    mpi::AllReduce( sendBuf, recvBuf, localSize, g.VCComm() );

    // Unpack
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        MemCopy
        ( &thisBuf[jLoc*thisLDim], &recvBuf[jLoc*localHeight], localHeight );
    this->auxMemory_.Release();
}

#define PROTO(T) template class DistMatrix<T,STAR,STAR>
#define COPY(T,CD,RD) \
  template DistMatrix<T,STAR,STAR>::DistMatrix( const DistMatrix<T,CD,RD>& A )
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
