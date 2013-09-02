/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Transpose.hpp"

namespace elem {

template<typename T>
DistMatrix<T,MC,MR>::DistMatrix( const elem::Grid& grid )
: AbstractDistMatrix<T>(grid)
{ this->SetShifts(); }

template<typename T>
DistMatrix<T,MC,MR>::DistMatrix
( Int height, Int width, const elem::Grid& grid )
: AbstractDistMatrix<T>(grid)
{ this->SetShifts(); this->ResizeTo( height, width ); }

template<typename T>
DistMatrix<T,MC,MR>::DistMatrix
( Int height, Int width, 
  Int colAlignment, Int rowAlignment, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ 
    this->Align( colAlignment, rowAlignment );
    this->ResizeTo( height, width );
}

template<typename T>
DistMatrix<T,MC,MR>::DistMatrix
( Int height, Int width,
  Int colAlignment, Int rowAlignment, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ 
    this->Align( colAlignment, rowAlignment );
    this->ResizeTo( height, width, ldim );
}

template<typename T>
DistMatrix<T,MC,MR>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ 
    this->LockedAttach
    ( height, width, colAlignment, rowAlignment, buffer, ldim, g ); 
}

template<typename T>
DistMatrix<T,MC,MR>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ 
    this->Attach
    ( height, width, colAlignment, rowAlignment, buffer, ldim, g );
}

template<typename T>
DistMatrix<T,MC,MR>::DistMatrix( const DistMatrix<T,MC,MR>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MC,MR]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [MC,MR] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DistMatrix<T,MC,MR>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MC,MR]::DistMatrix");
#endif
    this->SetShifts();
    if( MC != U || MR != V ||
        reinterpret_cast<const DistMatrix<T,MC,MR>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [MC,MR] with itself");
}

template<typename T>
DistMatrix<T,MC,MR>::DistMatrix( DistMatrix<T,MC,MR>&& A )
: AbstractDistMatrix<T>(std::move(A))
{ }

template<typename T>
DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( DistMatrix<T,MC,MR>&& A )
{
    AbstractDistMatrix<T>::operator=( std::move(A) );
    return *this;
}

template<typename T>
DistMatrix<T,MC,MR>::~DistMatrix()
{ }

template<typename T>
elem::DistData
DistMatrix<T,MC,MR>::DistData() const
{
    elem::DistData data;
    data.colDist = MC;
    data.rowDist = MR;
    data.colAlignment = this->colAlignment_;
    data.rowAlignment = this->rowAlignment_;
    data.root = 0;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T>
Int
DistMatrix<T,MC,MR>::ColStride() const
{ return this->grid_->Height(); }

template<typename T>
Int
DistMatrix<T,MC,MR>::RowStride() const
{ return this->grid_->Width(); }

template<typename T>
Int
DistMatrix<T,MC,MR>::ColRank() const
{ return this->grid_->Row(); }

template<typename T>
Int
DistMatrix<T,MC,MR>::RowRank() const
{ return this->grid_->Col(); }

template<typename T>
void
DistMatrix<T,MC,MR>::AlignWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AlignWith");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );
    if( data.colDist == MC && data.rowDist == MR )
    {
        this->colAlignment_ = data.colAlignment;
        this->rowAlignment_ = data.rowAlignment;
        this->constrainedColAlignment_ = true;
        this->constrainedRowAlignment_ = true;
    }
    else if( data.colDist == MC && data.rowDist == STAR )
    {
        this->colAlignment_ = data.colAlignment;
        this->constrainedColAlignment_ = true; 
    }
    else if( data.colDist == MR && data.rowDist == MC )
    {
        this->colAlignment_ = data.rowAlignment;
        this->rowAlignment_ = data.colAlignment;
        this->constrainedColAlignment_ = true;
        this->constrainedRowAlignment_ = true;
    }
    else if( data.colDist == MR && data.rowDist == STAR )
    {
        this->rowAlignment_ = data.colAlignment;
        this->constrainedRowAlignment_ = true;
    }
    else if( data.colDist == STAR && data.rowDist == MC )
    {
        this->colAlignment_ = data.rowAlignment;
        this->constrainedColAlignment_ = true;
    }
    else if( data.colDist == STAR && data.rowDist == MR )
    {
        this->rowAlignment_ = data.rowAlignment;
        this->constrainedRowAlignment_ = true;
    }
    else if( data.colDist == STAR && data.rowDist == VC )
    {
        this->colAlignment_ = data.rowAlignment % this->ColStride();
        this->constrainedColAlignment_ = true;
    }
    else if( data.colDist == STAR && data.rowDist == VR )
    {
        this->rowAlignment_ = data.rowAlignment % this->RowStride();
        this->constrainedRowAlignment_ = true;
    }
    else if( data.colDist == VC && data.rowDist == STAR )
    {
        this->colAlignment_ = data.colAlignment % this->ColStride();
        this->constrainedColAlignment_ = true;
    }
    else if( data.colDist == VR && data.rowDist == STAR )
    {
        this->rowAlignment_ = data.colAlignment % this->RowStride();
        this->constrainedRowAlignment_ = true;
    }
#ifndef RELEASE
    else LogicError("Nonsensical alignment");
#endif
    this->SetShifts();
}

template<typename T>
void
DistMatrix<T,MC,MR>::AlignWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,MC,MR>::AlignColsWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AlignColsWith");
    if( *this->grid_ != *data.grid )
        LogicError("Grids do not match");
#endif
    if( data.colDist == MC )
        this->colAlignment_ = data.colAlignment;
    else if( data.rowDist == MC )
        this->colAlignment_ = data.rowAlignment;
    else if( data.colDist == VC )
        this->colAlignment_ = data.colAlignment % this->ColStride();
    else if( data.rowDist == VC )
        this->colAlignment_ = data.rowAlignment % this->ColStride();
#ifndef RELEASE
    else LogicError("Nonsensical alignment");
#endif
    this->constrainedColAlignment_ = true;
    this->SetShifts();
}

template<typename T>
void
DistMatrix<T,MC,MR>::AlignColsWith( const AbstractDistMatrix<T>& A )
{ this->AlignColsWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,MC,MR>::AlignRowsWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AlignRowsWith");
    if( *this->grid_ != *data.grid )
        LogicError("Grids do not match");
#endif
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
DistMatrix<T,MC,MR>::AlignRowsWith( const AbstractDistMatrix<T>& A )
{ this->AlignRowsWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,MC,MR>::Attach
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::Attach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->rowAlignment_ = rowAlignment;
    this->viewType_ = VIEW;
    this->SetShifts();
    if( this->Participating() )
    {
        Int localHeight = Length(height,this->colShift_,this->ColStride());
        Int localWidth = Length(width,this->rowShift_,this->RowStride());
        this->matrix_.Attach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::LockedAttach
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->rowAlignment_ = rowAlignment;
    this->viewType_ = LOCKED_VIEW;
    this->SetShifts();
    if( this->Participating() )
    {
        Int localHeight = Length(height,this->colShift_,this->ColStride());
        Int localWidth = Length(width,this->rowShift_,this->RowStride());
        this->matrix_.LockedAttach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::ResizeTo");
    this->AssertNotLocked();
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),this->ColStride()),
          Length(width, this->RowShift(),this->RowStride()) );
}

template<typename T>
void
DistMatrix<T,MC,MR>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::ResizeTo");
    this->AssertNotLocked();
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),this->ColStride()),
          Length(width, this->RowShift(),this->RowStride()), ldim );
}

template<typename T>
T
DistMatrix<T,MC,MR>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const Int ownerRow = (i + this->ColAlignment()) % this->ColStride();
    const Int ownerCol = (j + this->RowAlignment()) % this->RowStride();
    const Int ownerRank = ownerRow + ownerCol*this->ColStride();

    T u;
    const elem::Grid& g = this->Grid();
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / this->ColStride();
        const Int jLoc = (j-this->RowShift()) / this->RowStride();
        u = this->GetLocal(iLoc,jLoc);
    }
    mpi::Broadcast( u, g.VCToViewingMap(ownerRank), g.ViewingComm() );
    return u;
}

template<typename T>
void
DistMatrix<T,MC,MR>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::Set");
    this->AssertValidEntry( i, j );
#endif
    const Int ownerRow = (i + this->ColAlignment()) % this->ColStride();
    const Int ownerCol = (j + this->RowAlignment()) % this->RowStride();
    const Int ownerRank = ownerRow + ownerCol*this->ColStride();
    if( this->Grid().VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / this->ColStride();
        const Int jLoc = (j-this->RowShift()) / this->RowStride();
        this->SetLocal(iLoc,jLoc,u);
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::Update");
    this->AssertValidEntry( i, j );
#endif
    const Int ownerRow = (i + this->ColAlignment()) % this->ColStride();
    const Int ownerCol = (j + this->RowAlignment()) % this->RowStride();
    const Int ownerRank = ownerRow + ownerCol*this->ColStride();
    if( this->Grid().VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / this->ColStride();
        const Int jLoc = (j-this->RowShift()) / this->RowStride();
        this->UpdateLocal(iLoc,jLoc,u);
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::GetDiagonal
( DistMatrix<T,MD,STAR>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
    {
        std::ostringstream os;
        os << mpi::WorldRank() << "\n"
           << "offset:         " << offset << "\n"
           << "colAlignment:   " << this->colAlignment_ << "\n"
           << "rowAlignment:   " << this->rowAlignment_ << "\n"
           << "d.diagPath:     " << d.diagPath_ << "\n"
           << "d.colAlignment: " << d.colAlignment_ << std::endl;
        std::cerr << os.str();
        LogicError("d must be aligned with the 'offset' diagonal");
    }
#endif
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( *this, offset );
    }
    const Int diagLength = this->DiagonalLength(offset);
    d.ResizeTo( diagLength, 1 );
    if( !d.Participating() )
        return;

    Int iStart, jStart;
    const Int diagShift = d.ColShift();
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int iLocStart = (iStart-colShift) / colStride;
    const Int jLocStart = (jStart-rowShift) / rowStride;

    const Int lcm = g.LCM();
    const Int localDiagLength = d.LocalHeight();
    T* dBuf = d.Buffer();
    const T* buffer = this->LockedBuffer();
    const Int ldim = this->LDim();

    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/colStride);
        const Int jLoc = jLocStart + k*(lcm/rowStride);
        dBuf[k] = buffer[iLoc+jLoc*ldim];
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::GetDiagonal
( DistMatrix<T,STAR,MD>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        LogicError("d must be aligned with the 'offset' diagonal");
#endif
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( *this, offset );
    }
    const Int diagLength = this->DiagonalLength(offset);
    d.ResizeTo( 1, diagLength );
    if( !d.Participating() )
        return;

    Int iStart, jStart;
    const Int diagShift = d.RowShift();
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int iLocStart = (iStart-colShift) / colStride;
    const Int jLocStart = (jStart-rowShift) / rowStride;

    const Int localDiagLength = d.LocalWidth();
    T* dBuf = d.Buffer();
    const Int dLDim = d.LDim();
    const T* buffer = this->LockedBuffer();
    const Int ldim = this->LDim();
    const Int lcm = g.LCM();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/colStride);
        const Int jLoc = jLocStart + k*(lcm/rowStride);
        dBuf[k*dLDim] = buffer[iLoc+jLoc*ldim];
    }
}

template<typename T>
DistMatrix<T,MD,STAR>
DistMatrix<T,MC,MR>::GetDiagonal( Int offset ) const
{
    DistMatrix<T,MD,STAR> d( this->Grid() );
    GetDiagonal( d, offset );
    return d;
}

template<typename T>
void
DistMatrix<T,MC,MR>::SetDiagonal
( const DistMatrix<T,MD,STAR>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    const Int diagLength = this->DiagonalLength(offset);
    if( diagLength != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        LogicError( msg.str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        LogicError("d must be aligned with the 'offset' diagonal");
#endif
    if( !d.Participating() )
        return;

    Int iStart,jStart;
    const Int diagShift = d.ColShift();
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int iLocStart = (iStart-colShift) / colStride;
    const Int jLocStart = (jStart-rowShift) / rowStride;

    const Int localDiagLength = d.LocalHeight();
    const T* dBuf = d.LockedBuffer();
    T* buffer = this->Buffer();
    const Int ldim = this->LDim();
    const Int lcm = this->Grid().LCM();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/colStride);
        const Int jLoc = jLocStart + k*(lcm/rowStride);
        buffer[iLoc+jLoc*ldim] = dBuf[k];
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::SetDiagonal
( const DistMatrix<T,STAR,MD>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Height() != 1 )
        LogicError("d must be a row vector");
    const Int diagLength = this->DiagonalLength(offset);
    if( diagLength != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        LogicError( msg.str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        LogicError("d must be aligned with the 'offset' diagonal");
#endif
    if( !d.Participating() )
        return;

    Int iStart,jStart;
    const Int diagShift = d.RowShift();
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int iLocStart = (iStart-colShift) / colStride;
    const Int jLocStart = (jStart-rowShift) / rowStride;

    const Int localDiagLength = d.LocalWidth();
    const T* dBuf = d.LockedBuffer();
    T* buffer = this->Buffer();
    const Int dLDim = d.LDim();
    const Int ldim = this->LDim();
    const Int lcm = this->Grid().LCM();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/colStride);
        const Int jLoc = jLocStart + k*(lcm/rowStride);
        buffer[iLoc+jLoc*ldim] = dBuf[k*dLDim];
    }
}

//
// Utility functions, e.g., TransposeFrom
//

template<typename T>
void
DistMatrix<T,MC,MR>::AdjointFrom( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AdjointFrom");
#endif
    this->TransposeFrom( A, true );
}

template<typename T>
void
DistMatrix<T,MC,MR>::TransposeFrom
( const DistMatrix<T,STAR,MC>& A, bool conjugate )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::TransposeFrom");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    this->SetColAlignmentAndResize( A.RowAlignment(), A.Width(), A.Height() );
    if( !this->Participating() )
        return;

    if( this->ColAlignment() == A.RowAlignment() )
    {
        const Int rowStride = this->RowStride();
        const Int rowShift = this->RowShift();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            T* destCol = &buffer[jLoc*ldim];
            const T* sourceCol = &ABuffer[rowShift+jLoc*rowStride];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc] = ( conjugate ? Conj(sourceCol[iLoc*ALDim])
                                            : sourceCol[iLoc*ALDim] );
        }
    }
    else
    {
        const Grid& g = this->Grid();
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MC,MR]::TransposeFrom." << std::endl;
#endif
        const Int colStride = this->ColStride();
        const Int rowStride = this->RowStride();
        const Int colRank = this->ColRank();
        const Int rowShift = this->RowShift();
        const Int colAlign = this->ColAlignment();
        const Int rowAlignA = A.RowAlignment();
        const Int sendRank = (colRank+colStride+colAlign-rowAlignA) % colStride;
        const Int recvRank = (colRank+colStride+rowAlignA-colAlign) % colStride;

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localWidthA = A.LocalWidth();
        const Int sendSize = localWidthA*localWidth;
        const Int recvSize = localHeight*localWidth;
        T* auxBuf = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &auxBuf[0];
        T* recvBuf = &auxBuf[sendSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            T* destCol = &sendBuf[jLoc*localWidth];
            const T* sourceCol = &ABuffer[rowShift+jLoc*rowStride];
            for( Int iLoc=0; iLoc<localWidthA; ++iLoc )
                destCol[iLoc] = ( conjugate ? Conj(sourceCol[iLoc*ALDim]) 
                                            : sourceCol[iLoc*ALDim] );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank, 
          recvBuf, recvSize, recvRank, g.ColComm() );

        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &buffer[jLoc*ldim], &recvBuf[jLoc*localHeight], localHeight );
        this->auxMemory_.Release();
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::AdjointFrom( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AdjointFrom");
#endif
    this->TransposeFrom( A, true );
}

template<typename T>
void
DistMatrix<T,MC,MR>::TransposeFrom
( const DistMatrix<T,MR,STAR>& A, bool conjugate )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::TransposeFrom");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    this->ForceRowAlignmentAndResize( A.ColAlignment(), A.Width(), A.Height() );
    if( !this->Participating() )
        return;

    const Int colStride = this->ColStride();
    const Int colShift = this->ColShift();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const T* ABuffer = A.LockedBuffer();
    const Int ALDim = A.LDim();
    T* buffer = this->Buffer();
    const Int ldim = this->LDim();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        T* destCol = &buffer[jLoc*ldim];
        const T* sourceCol = &ABuffer[jLoc+colShift*ALDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            destCol[iLoc] = 
                ( conjugate ? Conj(sourceCol[iLoc*colStride*ALDim])
                            : sourceCol[iLoc*colStride*ALDim] );
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::AdjointSumScatterFrom
( const DistMatrix<T,MR,STAR>& AAdj_MR_STAR )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AdjointSumScatterFrom");
#endif
    this->TransposeSumScatterFrom( AAdj_MR_STAR, true );
}

template<typename T>
void
DistMatrix<T,MC,MR>::TransposeSumScatterFrom
( const DistMatrix<T,MR,STAR>& ATrans_MR_STAR, bool conjugate )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::TransposeSumScatterFrom");
#endif
    const Grid& g = ATrans_MR_STAR.Grid();
    DistMatrix<T,MR,MC> ATrans( g );
    if( this->Viewing() )
        ATrans.AlignWith( *this );
    ATrans.SumScatterFrom( ATrans_MR_STAR );
    Transpose( ATrans, *this, conjugate );
}

template<typename T>
void
DistMatrix<T,MC,MR>::AdjointSumScatterUpdate
( T alpha, const DistMatrix<T,MR,STAR>& AAdj_MR_STAR )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AdjointSumScatterUpdate");
#endif
    this->TransposeSumScatterUpdate( alpha, AAdj_MR_STAR, true );
}

template<typename T>
void
DistMatrix<T,MC,MR>::TransposeSumScatterUpdate
( T alpha, const DistMatrix<T,MR,STAR>& ATrans_MR_STAR, bool conjugate )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::TransposeSumScatterUpdate");
#endif
    const Grid& g = ATrans_MR_STAR.Grid();
    DistMatrix<T,MR,MC> ATrans( g );
    ATrans.SumScatterFrom( ATrans_MR_STAR );
    DistMatrix<T,MC,MR> A( g );
    if( this->Viewing() )
        A.AlignWith( *this );
    Transpose( ATrans, A, conjugate );
    Axpy( alpha, A, *this );
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [MC,MR]");
    this->AssertNotLocked();
#endif
    if( this->Grid() == A.Grid() )
    {
        this->SetAlignmentsAndResize
        ( A.ColAlignment(), A.RowAlignment(), A.Height(), A.Width() );
        if( !this->Participating() && !A.Participating() )
            return *this;
        if( this->ColAlignment() == A.ColAlignment() &&
            this->RowAlignment() == A.RowAlignment() )
        {
            this->matrix_ = A.LockedMatrix();
        }
        else
        {
            const elem::Grid& g = this->Grid();
#ifdef UNALIGNED_WARNINGS
            if( g.Rank() == 0 )
                std::cerr << "Unaligned [MC,MR] <- [MC,MR]." << std::endl;
#endif
            const Int colRank = this->ColRank();
            const Int rowRank = this->RowRank();
            const Int colStride = this->ColStride();
            const Int rowStride = this->RowStride();
            const Int colAlignment = this->ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int colAlignmentA = A.ColAlignment();
            const Int rowAlignmentA = A.RowAlignment();
            const Int colDiff = colAlignment - colAlignmentA;
            const Int rowDiff = rowAlignment - rowAlignmentA;
            const Int sendRow = (colRank + colStride + colDiff) % colStride;
            const Int recvRow = (colRank + colStride - colDiff) % colStride;
            const Int sendCol = (rowRank + rowStride + rowDiff) % rowStride;
            const Int recvCol = (rowRank + rowStride - rowDiff) % rowStride;
            const Int sendRank = sendRow + sendCol*colStride;
            const Int recvRank = recvRow + recvCol*colStride;

            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localHeightA = A.LocalHeight();
            const Int localWidthA = A.LocalWidth();
            const Int sendSize = localHeightA*localWidthA;
            const Int recvSize = localHeight*localWidth;
            T* auxBuf = this->auxMemory_.Require( sendSize + recvSize );
            T* sendBuf = &auxBuf[0];
            T* recvBuf = &auxBuf[sendSize];

            // Pack
            const Int ALDim = A.LDim();
            const T* ABuffer = A.LockedBuffer();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
                MemCopy
                ( &sendBuf[jLoc*localHeightA], 
                  &ABuffer[jLoc*ALDim], localHeightA );

            // Communicate
            mpi::SendRecv
            ( sendBuf, sendSize, sendRank, 
              recvBuf, recvSize, recvRank, g.VCComm() );

            // Unpack
            T* buffer = this->Buffer();
            const Int ldim = this->LDim();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &buffer[jLoc*ldim], 
                  &recvBuf[jLoc*localHeight], localHeight );
            this->auxMemory_.Release();
        }
    }
    else // the grids don't match
    {
        CopyFromDifferentGrid( A );
    }
    return *this;
}

template<typename T>
void DistMatrix<T,MC,MR>::CopyFromDifferentGrid( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MC,MR]::CopyFromDifferentGrid");
#endif
    this->ResizeTo( A.Height(), A.Width() ); 
    // Just need to ensure that each viewing comm contains the other team's
    // owning comm. Congruence is too strong.

    // Compute the number of process rows and columns that each process 
    // needs to send to.
    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int colRank = this->ColRank();
    const Int rowRank = this->RowRank();
    const Int colStrideA = A.ColStride();
    const Int rowStrideA = A.RowStride();
    const Int colRankA = A.ColRank();
    const Int rowRankA = A.RowRank();
    const Int colGCD = GCD( colStride, colStrideA );
    const Int rowGCD = GCD( rowStride, rowStrideA );
    const Int colLCM = colStride*colStrideA / colGCD;
    const Int rowLCM = rowStride*rowStrideA / rowGCD;
    const Int numColSends = colStride / colGCD;
    const Int numRowSends = rowStride / rowGCD;
    const Int localColStride = colLCM / colStride;
    const Int localRowStride = rowLCM / rowStride;
    const Int localColStrideA = numColSends;
    const Int localRowStrideA = numRowSends;

    const Int colAlign = this->ColAlignment();
    const Int rowAlign = this->RowAlignment();
    const Int colAlignA = A.ColAlignment();
    const Int rowAlignA = A.RowAlignment();

    const bool inThisGrid = this->Participating();
    const bool inAGrid = A.Participating();
    if( !inThisGrid && !inAGrid )
        return;

    const Int maxSendSize = 
        (A.Height()/(colStrideA*localColStrideA)+1) * 
        (A.Width()/(rowStrideA*localRowStrideA)+1);

    // Translate the ranks from A's VC communicator to this's viewing so that
    // we can match send/recv communicators
    const int sizeA = A.Grid().Size();
    std::vector<int> rankMap(sizeA), ranks(sizeA);
    for( int j=0; j<sizeA; ++j )
        ranks[j] = j;
    mpi::Group viewingGroup;
    mpi::CommGroup( this->Grid().ViewingComm(), viewingGroup );
    mpi::GroupTranslateRanks
    ( A.Grid().OwningGroup(), sizeA, &ranks[0], viewingGroup, &rankMap[0] );

    // Have each member of A's grid individually send to all numRow x numCol
    // processes in order, while the members of this grid receive from all 
    // necessary processes at each step.
    Int requiredMemory = 0;
    if( inAGrid )
        requiredMemory += maxSendSize;
    if( inThisGrid )
        requiredMemory += maxSendSize;
    T* auxBuf = this->auxMemory_.Require( requiredMemory );
    Int offset = 0;
    T* sendBuf = &auxBuf[offset];
    if( inAGrid )
        offset += maxSendSize;
    T* recvBuf = &auxBuf[offset];

    Int recvRow = 0; // avoid compiler warnings...
    if( inAGrid )
        recvRow = (((colRankA+colStrideA-colAlignA)%colStrideA)+colAlign) % 
                  colStride;
    for( Int colSend=0; colSend<numColSends; ++colSend )
    {
        Int recvCol = 0; // avoid compiler warnings...
        if( inAGrid )
            recvCol = (((rowRankA+rowStrideA-rowAlignA)%rowStrideA)+rowAlign) % 
                      rowStride;
        for( Int rowSend=0; rowSend<numRowSends; ++rowSend )
        {
            mpi::Request sendRequest;
            // Fire off this round of non-blocking sends
            if( inAGrid )
            {
                // Pack the data
                Int sendHeight = Length(A.LocalHeight(),colSend,numColSends);
                Int sendWidth = Length(A.LocalWidth(),rowSend,numRowSends);
                const T* ABuffer = A.LockedBuffer();
                const Int ALDim = A.LDim();
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<sendWidth; ++jLoc )
                {
                    const Int j = rowSend+jLoc*localRowStrideA;
                    for( Int iLoc=0; iLoc<sendHeight; ++iLoc )
                    {
                        const Int i = colSend+iLoc*localColStrideA;
                        sendBuf[iLoc+jLoc*sendHeight] = ABuffer[i+j*ALDim];
                    }
                }
                // Send data
                const Int recvVCRank = recvRow + recvCol*colStride;
                const Int recvViewingRank = 
                    this->Grid().VCToViewingMap( recvVCRank );
                mpi::ISend
                ( sendBuf, sendHeight*sendWidth, recvViewingRank,
                  this->Grid().ViewingComm(), sendRequest );
            }
            // Perform this round of recv's
            if( inThisGrid )
            {
                const Int sendColOffset = (colSend*colStrideA+colAlignA) % colStrideA;
                const Int recvColOffset = (colSend*colStrideA+colAlign) % colStride;
                const Int sendRowOffset = (rowSend*rowStrideA+rowAlignA) % rowStrideA;
                const Int recvRowOffset = (rowSend*rowStrideA+rowAlign) % rowStride;

                const Int firstSendRow = (((colRank+colStride-recvColOffset)%colStride)+sendColOffset)%colStrideA;
                const Int firstSendCol = (((rowRank+rowStride-recvRowOffset)%rowStride)+sendRowOffset)%rowStrideA;

                const Int colShift = (colRank+colStride-recvColOffset)%colStride;
                const Int rowShift = (rowRank+rowStride-recvRowOffset)%rowStride;
                const Int numColRecvs = Length( colStrideA, colShift, colStride ); 
                const Int numRowRecvs = Length( rowStrideA, rowShift, rowStride );

                // Recv data
                // For now, simply receive sequentially. Until we switch to 
                // nonblocking recv's, we won't be using much of the 
                // recvBuf
                Int sendRow = firstSendRow;
                for( Int colRecv=0; colRecv<numColRecvs; ++colRecv )
                {
                    const Int sendColShift = Shift( sendRow, colAlignA, colStrideA ) + colSend*colStrideA;
                    const Int sendHeight = Length( A.Height(), sendColShift, colLCM );
                    const Int localColOffset = (sendColShift-this->ColShift()) / colStride;

                    Int sendCol = firstSendCol;
                    for( Int rowRecv=0; rowRecv<numRowRecvs; ++rowRecv )
                    {
                        const Int sendRowShift = Shift( sendCol, rowAlignA, rowStrideA ) + rowSend*rowStrideA;
                        const Int sendWidth = Length( A.Width(), sendRowShift, rowLCM );
                        const Int localRowOffset = (sendRowShift-this->RowShift()) / rowStride;

                        const Int sendVCRank = sendRow+sendCol*colStrideA;
                        mpi::Recv
                        ( recvBuf, sendHeight*sendWidth, rankMap[sendVCRank],
                          this->Grid().ViewingComm() );
                        
                        // Unpack the data
                        T* buffer = this->Buffer();
                        const Int ldim = this->LDim();
                        PARALLEL_FOR
                        for( Int jLoc=0; jLoc<sendWidth; ++jLoc )
                        {
                            const Int j = localRowOffset+jLoc*localRowStride;
                            for( Int iLoc=0; iLoc<sendHeight; ++iLoc )
                            {
                                const Int i = localColOffset+iLoc*localColStride;
                                buffer[i+j*ldim] = recvBuf[iLoc+jLoc*sendHeight];
                            }
                        }
                        // Set up the next send col
                        sendCol = (sendCol + rowStride) % rowStrideA;
                    }
                    // Set up the next send row
                    sendRow = (sendRow + colStride) % colStrideA;
                }
            }
            // Ensure that this round of non-blocking sends completes
            if( inAGrid )
            {
                mpi::Wait( sendRequest );
                recvCol = (recvCol + rowStrideA) % rowStride;
            }
        }
        if( inAGrid )
            recvRow = (recvRow + colStrideA) % colStride;
    }
    this->auxMemory_.Release();
}

// PAUSED PASS HERE

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetColAlignmentAndResize( A.ColAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlignment() == A.ColAlignment() )
    {
        const Int c = g.Width();
        const Int rowShift = this->RowShift();

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();

        const Int ALDim = A.LDim();
        const Int thisLDim = this->LDim();
        T* thisBuffer = this->Buffer();
        const T* ABuffer = A.LockedBuffer();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuffer[jLoc*thisLDim], 
              &ABuffer[(rowShift+jLoc*c)*ALDim], localHeight );
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MC,MR] <- [MC,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int rank = g.Row();
        const Int rowShift = this->RowShift();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentA = A.ColAlignment();

        const Int sendRank = (rank+r+colAlignment-colAlignmentA) % r;
        const Int recvRank = (rank+r+colAlignmentA-colAlignment) % r;

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localHeightA = A.LocalHeight();

        const Int sendSize = localHeightA*localWidth;
        const Int recvSize = localHeight*localWidth;

        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &sendBuf[jLoc*localHeightA], 
              &ABuffer[(rowShift+jLoc*c)*ALDim], localHeightA );

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank, 
          recvBuf, recvSize, recvRank, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuffer[jLoc*thisLDim], 
              &recvBuf[jLoc*localHeight], localHeight );
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetRowAlignmentAndResize( A.RowAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int colShift = this->ColShift();

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();

        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            T* destCol = &thisBuffer[jLoc*thisLDim];
            const T* sourceCol = &ABuffer[colShift+jLoc*ALDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*r];
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MC,MR] <- [* ,MR]." << std::endl;
#endif
        const Int r = g.Height(); 
        const Int c = g.Width();
        const Int col = g.Col();
        const Int colShift = this->ColShift();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentA = A.RowAlignment();

        const Int sendCol = (col+c+rowAlignment-rowAlignmentA) % c;
        const Int recvCol = (col+c+rowAlignmentA-rowAlignment) % c;

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localWidthA = A.LocalWidth();

        const Int sendSize = localHeight*localWidthA;
        const Int recvSize = localHeight*localWidth;

        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
        {
            T* destCol = &sendBuf[jLoc*localHeight];
            const T* sourceCol = &ABuffer[colShift+jLoc*ALDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*r];
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendCol, 
          recvBuf, recvSize, recvCol, g.RowComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuffer[jLoc*thisLDim], 
              &recvBuf[jLoc*localHeight], localHeight );
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [MD,* ]");
#endif
    LogicError("[MC,MR] = [MD,* ] not yet implemented");
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,MD>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,MD]");
#endif
    LogicError("[MC,MR] = [* ,MD] not yet implemented");
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( A.Width() == 1 )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int myRow = g.Row();
        const Int myCol = g.Col();
        const Int rankCM = g.VCRank();
        const Int rankRM = g.VRRank();
        const Int ownerCol = this->RowAlignment();
        const Int ownerRow = A.RowAlignment();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentA = A.ColAlignment();
        const Int colShift = this->ColShift();
        const Int colShiftA = A.ColShift();

        const Int height = A.Height();
        const Int maxLocalHeight = MaxLength(height,p);
        const Int portionSize = mpi::Pad( maxLocalHeight );

        const Int colShiftVC = Shift(rankCM,colAlignment,p);
        const Int colShiftVRA = Shift(rankRM,colAlignmentA,p);
        const Int sendRankCM = (rankCM+(p+colShiftVRA-colShiftVC)) % p;
        const Int recvRankRM = (rankRM+(p+colShiftVC-colShiftVRA)) % p;
        const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        T* buffer = this->auxMemory_.Require( (r+c)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        if( myRow == ownerRow )
        {
            // Pack
            const T* ABuffer = A.LockedBuffer();
            PARALLEL_FOR
            for( Int k=0; k<r; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const Int shift = Shift_(myCol+c*k,colAlignmentA,p);
                const Int offset = (shift-colShiftA) / c;
                const Int thisLocalHeight = Length_(height,shift,p);

                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    data[iLoc] = ABuffer[offset+iLoc*r];
            }
        }

        // A[VR,* ] <- A[MR,MC]
        mpi::Scatter
        ( recvBuf, portionSize, sendBuf, portionSize, ownerRow, g.ColComm() );

        // A[VC,* ] <- A[VR,* ]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankCM,
          recvBuf, portionSize, recvRankCM, g.VCComm() );

        // A[MC,MR] <- A[VC,* ]
        mpi::Gather
        ( recvBuf, portionSize, sendBuf, portionSize, ownerCol, g.RowComm() );

        if( myCol == ownerCol )
        {
            // Unpack
            T* thisBuffer = this->Buffer();
            PARALLEL_FOR
            for( Int k=0; k<c; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const Int shift = Shift_(myRow+r*k,colAlignment,p);
                const Int offset = (shift-colShift) / r;
                const Int thisLocalHeight = Length_(height,shift,p);

                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    thisBuffer[offset+iLoc*c] = data[iLoc];
            }
        }
        this->auxMemory_.Release();
    }
    else if( A.Height() == 1 )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int myRow = g.Row();
        const Int myCol = g.Col();
        const Int rankCM = g.VCRank();
        const Int rankRM = g.VRRank();
        const Int ownerRow = this->ColAlignment();
        const Int ownerCol = A.ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentA = A.RowAlignment();
        const Int rowShift = this->RowShift();
        const Int rowShiftA = A.RowShift();

        const Int width = A.Width();
        const Int maxLocalWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxLocalWidth );

        const Int rowShiftVR = Shift(rankRM,rowAlignment,p);
        const Int rowShiftVCA = Shift(rankCM,rowAlignmentA,p);
        const Int sendRankRM = (rankRM+(p+rowShiftVCA-rowShiftVR)) % p;
        const Int recvRankCM = (rankCM+(p+rowShiftVR-rowShiftVCA)) % p;
        const Int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

        T* buffer = this->auxMemory_.Require( (r+c)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        if( myCol == ownerCol )
        {
            // Pack
            const T* ABuffer = A.LockedBuffer();
            const Int ALDim = A.LDim();
            PARALLEL_FOR
            for( Int k=0; k<c; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const Int shift = Shift_(myRow+r*k,rowAlignmentA,p);
                const Int offset = (shift-rowShiftA) / r;
                const Int thisLocalWidth = Length_(width,shift,p);

                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                    data[jLoc] = ABuffer[(offset+jLoc*c)*ALDim];
            }
        }

        // A[* ,VC] <- A[MR,MC]
        mpi::Scatter
        ( recvBuf, portionSize, sendBuf, portionSize, ownerCol, g.RowComm() );

        // A[* ,VR] <- A[* ,VC]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankRM,
          recvBuf, portionSize, recvRankRM, g.VRComm() );

        // A[MC,MR] <- A[* ,VR]
        mpi::Gather
        ( recvBuf, portionSize, sendBuf, portionSize, ownerRow, g.ColComm() );
    
        if( myRow == ownerRow )
        {
            // Unpack
            T* thisBuffer = this->Buffer();
            const Int thisLDim = this->LDim();
            PARALLEL_FOR
            for( Int k=0; k<r; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const Int shift = Shift_(myCol+c*k,rowAlignment,p);
                const Int offset = (shift-rowShift) / c;
                const Int thisLocalWidth = Length_(width,shift,p);

                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                    thisBuffer[(offset+jLoc*r)*thisLDim] = data[jLoc];
            }
        }

        this->auxMemory_.Release();
    }
    else
    {
        if( A.Height() >= A.Width() )
        {
            std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
            ( new DistMatrix<T,VR,STAR>(g) );

            *A_VR_STAR = A;

            std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
            ( new DistMatrix<T,VC,STAR>(true,this->ColAlignment(),g) );
            *A_VC_STAR = *A_VR_STAR;
            delete A_VR_STAR.release(); // lowers memory highwater

            *this = *A_VC_STAR;
        }
        else
        {
            std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
            ( new DistMatrix<T,STAR,VC>(g) );
            *A_STAR_VC = A;

            std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
            ( new DistMatrix<T,STAR,VR>(true,this->RowAlignment(),g) );
            *A_STAR_VR = *A_STAR_VC;
            delete A_STAR_VC.release(); // lowers memory highwater

            *this = *A_STAR_VR;
            this->ResizeTo( A_STAR_VR->Height(), A_STAR_VR->Width() );
        }
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [MR,* ]");
#endif
    const Grid& g = A.Grid();
    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(A) );
    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,MC]");
#endif
    const Grid& g = A.Grid();
    std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(A) );
    std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(true,this->RowAlignment(),g) );
    *A_STAR_VR = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater
    *this = *A_STAR_VR;
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetColAlignmentAndResize
    ( A.ColAlignment()%g.Height(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlignment() == A.ColAlignment() % g.Height() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int row = g.Row();
        const Int colShift = this->ColShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentA = A.ColAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightA = A.LocalHeight();

        const Int maxHeight = MaxLength(height,p);
        const Int maxWidth = MaxLength(width,c);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*c*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        // Pack
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuf[k*portionSize];
            const Int thisRowShift = Shift_(k,rowAlignment,c);
            const Int thisLocalWidth = Length_(width,thisRowShift,c);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                MemCopy
                ( &data[jLoc*localHeightA], 
                  &ABuffer[(thisRowShift+jLoc*c)*ALDim], localHeightA );
        }

        // Communicate
        mpi::AllToAll
        ( sendBuf, portionSize, recvBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisRank = row+k*r;
            const Int thisColShift = Shift_(thisRank,colAlignmentA,p);
            const Int thisColOffset = (thisColShift-colShift) / r;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuffer[thisColOffset+jLoc*thisLDim];
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
            std::cerr << "Unaligned [MC,MR] <- [VC,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int row = g.Row();
        const Int colShift = this->ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentA = A.ColAlignment();

        const Int sendRow = (row+r+colAlignment-(colAlignmentA%r)) % r;
        const Int recvRow = (row+r+(colAlignmentA%r)-colAlignment) % r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightA = A.LocalHeight();

        const Int maxHeight = MaxLength(height,p);
        const Int maxWidth = MaxLength(width,c);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*c*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[c*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int thisRowShift = Shift_(k,rowAlignment,c);
            const Int thisLocalWidth = Length_(width,thisRowShift,c);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                MemCopy
                ( &data[jLoc*localHeightA], 
                  &ABuffer[(thisRowShift+jLoc*c)*ALDim], localHeightA );
        }

        // SendRecv: properly align A[VC,*] via a trade in the column
        mpi::SendRecv
        ( secondBuf, c*portionSize, sendRow,
          firstBuf,  c*portionSize, recvRow, g.ColComm() );

        // AllToAll to gather all of the aligned A[VC,*] data into 
        // secondBuff.
        mpi::AllToAll
        ( firstBuf,  portionSize, secondBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int thisRank = recvRow+k*r;
            const Int thisColShift = Shift_(thisRank,colAlignmentA,p);
            const Int thisColOffset = (thisColShift-colShift) / r;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuffer[thisColOffset+jLoc*thisLDim];
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
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,VC]");
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR> A_STAR_VR(true,this->RowAlignment(),g);
    A_STAR_VR = A;
    *this = A_STAR_VR;
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [VR,* ]");
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR> A_VC_STAR(true,this->ColAlignment(),g);
    A_VC_STAR = A;
    *this = A_VC_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,VR]");
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
        const Int p = r * c;
        const Int col = g.Col();
        const Int rowShift = this->RowShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignmentA = A.RowAlignment();
    
        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthA = A.LocalWidth();

        const Int maxHeight = MaxLength(height,r);
        const Int maxWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*r*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &sendBuf[k*portionSize];
            const Int thisColShift = Shift_(k,colAlignment,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuf, portionSize, recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisRank = col+k*c;
            const Int thisRowShift = Shift_(thisRank,rowAlignmentA,p);
            const Int thisRowOffset = (thisRowShift-rowShift) / c;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                MemCopy
                ( &thisBuffer[(thisRowOffset+jLoc*r)*thisLDim], 
                  &data[jLoc*localHeight], localHeight );
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MC,MR] <- [* ,VR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int col = g.Col();
        const Int rowShift = this->RowShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentA = A.RowAlignment();

        const Int sendCol = (col+c+rowAlignment-(rowAlignmentA%c)) % c;
        const Int recvCol = (col+c+(rowAlignmentA%c)-rowAlignment) % c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthA = A.LocalWidth();
            
        const Int maxHeight = MaxLength(height,r);
        const Int maxWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*r*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[r*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int thisColShift = Shift_(k,colAlignment,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // SendRecv: properly align A[*,VR] via a trade in the column
        mpi::SendRecv
        ( secondBuf, r*portionSize, sendCol,
          firstBuf,  r*portionSize, recvCol, g.RowComm() );

        // AllToAll to gather all of the aligned [*,VR] data into 
        // secondBuf
        mpi::AllToAll
        ( firstBuf, portionSize, secondBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int thisRank = recvCol+k*c;
            const Int thisRowShift = Shift_(thisRank,rowAlignmentA,p);
            const Int thisRowOffset = (thisRowShift-rowShift) / c;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                MemCopy
                ( &thisBuffer[(thisRowOffset+jLoc*r)*thisLDim], 
                  &data[jLoc*localHeight], localHeight );
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int r = this->Grid().Height();
    const Int c = this->Grid().Width();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();

    const Int ALDim = A.LDim();
    const Int thisLDim = this->LDim();
    T* thisBuffer = this->Buffer();
    const T* ABuffer = A.LockedBuffer();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        T* destCol = &thisBuffer[jLoc*thisLDim];
        const T* sourceCol = &ABuffer[colShift+(rowShift+jLoc*c)*ALDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            destCol[iLoc] = sourceCol[iLoc*r];
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [o ,o ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int p = g.Size();
    this->ResizeTo( m, n );

    const Int colAlignment = this->ColAlignment();
    const Int rowAlignment = this->RowAlignment();
    const Int mLocal = this->LocalHeight();
    const Int nLocal = this->LocalWidth();
    const Int pkgSize = mpi::Pad(MaxLength(m,colStride)*MaxLength(n,rowStride));
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
        const T* ABuffer = A.LockedBuffer();
        for( Int t=0; t<rowStride; ++t )
        {
            const Int tLocalWidth = Length( n, t, rowStride );
            const Int col = (rowAlignment+t) % rowStride;
            for( Int s=0; s<colStride; ++s )
            {
                const Int sLocalHeight = Length( m, s, colStride );
                const Int row = (colAlignment+s) % colStride;
                const Int q = row + col*colStride;
                for( Int jLoc=0; jLoc<tLocalWidth; ++jLoc ) 
                {
                    const Int j = t + jLoc*rowStride;
                    for( Int iLoc=0; iLoc<sLocalHeight; ++iLoc )
                    {
                        const Int i = s + iLoc*colStride;
                        sendBuf[q*pkgSize+iLoc+jLoc*sLocalHeight] = 
                            ABuffer[i+j*ALDim];
                    }
                }
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
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                buffer[iLoc+jLoc*ldim] = recvBuf[iLoc+jLoc*mLocal];     
        this->auxMemory_.Release();
    }

    return *this;
}

template<typename T>
void
DistMatrix<T,MC,MR>::SumScatterFrom( const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SumScatterFrom([MC,* ])");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetColAlignmentAndResize( A.ColAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return;

    if( this->ColAlignment() == A.ColAlignment() )
    {
        if( this->Width() == 1 )
        {
            const Int rowAlignment = this->RowAlignment();
            const Int myCol = g.Col();

            const Int localHeight = this->LocalHeight();
            const Int recvSize = mpi::Pad( localHeight );
            const Int sendSize = recvSize;

            T* buffer = this->auxMemory_.Require( sendSize + recvSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[sendSize];

            // Pack 
            const T* ACol = A.LockedBuffer();
            MemCopy( sendBuf, ACol, localHeight );

            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuf, recvBuf, sendSize, rowAlignment, g.RowComm() );

            if( myCol == rowAlignment )
                MemCopy( this->Buffer(), recvBuf, localHeight );
            this->auxMemory_.Release();
        }
        else
        {
            const Int c = g.Width();
            const Int rowAlignment = this->RowAlignment();
            
            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int maxLocalWidth = MaxLength(width,c);
            const Int recvSize = mpi::Pad( localHeight*maxLocalWidth );
            const Int sendSize = c * recvSize;

            // Pack 
            const Int ALDim = A.LDim();
            const T* ABuffer = A.LockedBuffer();
            T* buffer = this->auxMemory_.Require( sendSize );
            OUTER_PARALLEL_FOR
            for( Int k=0; k<c; ++k )
            {
                T* data = &buffer[k*recvSize];
                const Int thisRowShift = Shift_( k, rowAlignment, c );
                const Int thisLocalWidth = Length_(width,thisRowShift,c);
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                    MemCopy
                    ( &data[jLoc*localHeight], 
                      &ABuffer[(thisRowShift+jLoc*c)*ALDim], localHeight );
            }

            // Communicate
            mpi::ReduceScatter( buffer, recvSize, g.RowComm() );

            // Unpack our received data
            T* thisBuffer = this->Buffer();
            const Int thisLDim = this->LDim();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &thisBuffer[jLoc*thisLDim], 
                  &buffer[jLoc*localHeight], localHeight );
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned SumScatterFrom [MC,MR] <- [MC,* ]." 
                      << std::endl;
#endif
        if( this->Width() == 1 )
        {
            const Int r = g.Height();
            const Int rowAlignment = this->RowAlignment();
            const Int myRow = g.Row();
            const Int myCol = g.Col();

            const Int height = this->Height();
            const Int localHeight = this->LocalHeight();
            const Int localHeightA = A.LocalHeight();
            const Int maxLocalHeight = MaxLength(height,r);
            const Int portionSize = mpi::Pad( maxLocalHeight );

            const Int colAlignment = this->ColAlignment();
            const Int colAlignmentA = A.ColAlignment();
            const Int sendRow = (myRow+r+colAlignment-colAlignmentA) % r;
            const Int recvRow = (myRow+r+colAlignmentA-colAlignment) % r;

            T* buffer = this->auxMemory_.Require( 2*portionSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[portionSize];

            // Pack 
            const T* ACol = A.LockedBuffer();
            MemCopy( sendBuf, ACol, localHeightA );
            
            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuf, recvBuf, portionSize, rowAlignment, g.RowComm() );

            if( myCol == rowAlignment )
            {
                // Perform the realignment
                mpi::SendRecv
                ( recvBuf, portionSize, sendRow,
                  sendBuf, portionSize, recvRow, g.ColComm() );
                MemCopy( this->Buffer(), sendBuf, localHeight );
            }
            this->auxMemory_.Release();
        }
        else
        {
            const Int r = g.Height();
            const Int c = g.Width();
            const Int row = g.Row();

            const Int colAlignment = this->ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int colAlignmentA = A.ColAlignment();
            const Int sendRow = (row+r+colAlignment-colAlignmentA) % r;
            const Int recvRow = (row+r+colAlignmentA-colAlignment) % r;

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localHeightA = A.LocalHeight();
            const Int maxLocalWidth = MaxLength(width,c);

            const Int recvSize_RS = mpi::Pad( localHeightA*maxLocalWidth );
            const Int sendSize_RS = c * recvSize_RS;
            const Int recvSize_SR = localHeight * localWidth;

            T* buffer = this->auxMemory_.Require
            ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );
            T* firstBuf = &buffer[0];
            T* secondBuf = &buffer[recvSize_RS];

            // Pack 
            // TODO: Stick an optional outer parallelization here?
            const Int ALDim = A.LDim();
            const T* ABuffer = A.LockedBuffer();
            for( Int k=0; k<c; ++k )
            {
                T* data = &secondBuf[k*recvSize_RS];
                const Int thisRowShift = Shift_( k, rowAlignment, c );
                const Int thisLocalWidth = Length_(width,thisRowShift,c);
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                    MemCopy
                    ( &data[jLoc*localHeightA], 
                      &ABuffer[(thisRowShift+jLoc*c)*ALDim], localHeightA );
            }

            // Reduce-scatter over each process row
            mpi::ReduceScatter( secondBuf, firstBuf, recvSize_RS, g.RowComm() );

            // Trade reduced data with the appropriate process row
            mpi::SendRecv
            ( firstBuf,  localHeightA*localWidth, sendRow,
              secondBuf, localHeight*localWidth,  recvRow, g.ColComm() );

            // Unpack the received data
            T* thisBuffer = this->Buffer();
            const Int thisLDim = this->LDim();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &thisBuffer[jLoc*thisLDim], 
                  &secondBuf[jLoc*localHeight], localHeight );
            this->auxMemory_.Release();
        }
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::SumScatterFrom( const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SumScatterFrom([* ,MR])");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Width() == 1 && g.Rank() == 0 )
    {
        std::cerr <<
          "The vector version of [MC,MR].SumScatterFrom([* ,MR]) does not "
          "yet have a vector version implemented, but it would only require"
          " a modification of the vector version of "
          "[MC,MR].SumScatterFrom([MC,* ])" << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "[MC,MR]::SumScatterFrom([* ,MR]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MR] matrix instead." << std::endl;
    }
#endif
    this->SetRowAlignmentAndResize( A.RowAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return;

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int colAlignment = this->ColAlignment();

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalHeight = MaxLength(height,r);

        const Int recvSize = mpi::Pad( maxLocalHeight*localWidth );
        const Int sendSize = r * recvSize;

        // Pack 
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
        T* buffer = this->auxMemory_.Require( sendSize );
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisColShift = Shift_(k,colAlignment,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.ColComm() );

        // Unpack our received data
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuffer[jLoc*thisLDim], 
              &buffer[jLoc*localHeight], localHeight );
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned SumScatterFrom [MC,MR] <- [* ,MR]." 
                      << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int col = g.Col();

        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentA = A.RowAlignment();
        const Int sendCol = (col+c+rowAlignment-rowAlignmentA) % c;
        const Int recvCol = (col+c+rowAlignmentA-rowAlignment) % c;

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localWidthA = A.LocalWidth();
        const Int maxLocalHeight = MaxLength(height,r);

        const Int recvSize_RS = mpi::Pad( maxLocalHeight*localWidthA );
        const Int sendSize_RS = r * recvSize_RS;
        const Int recvSize_SR = localHeight * localWidth;

        T* buffer = this->auxMemory_.Require
        ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[recvSize_RS];

        // Pack 
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &secondBuf[k*recvSize_RS];
            const Int thisColShift = Shift_(k,colAlignment,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Reduce-scatter over each process col
        mpi::ReduceScatter( secondBuf, firstBuf, recvSize_RS, g.ColComm() );

        // Trade reduced data with the appropriate process col
        mpi::SendRecv
        ( firstBuf,  localHeight*localWidthA, sendCol,
          secondBuf, localHeight*localWidth,  recvCol, g.RowComm() );

        // Unpack the received data
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuffer[jLoc*thisLDim], 
              &secondBuf[jLoc*localHeight], localHeight );
        this->auxMemory_.Release();
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::SumScatterFrom( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SumScatterFrom([* ,* ])");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int colAlignment = this->ColAlignment();
    const Int rowAlignment = this->RowAlignment();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int maxLocalHeight = MaxLength(height,r);
    const Int maxLocalWidth = MaxLength(width,c);

    const Int recvSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
    const Int sendSize = r*c*recvSize;

    // Pack 
    const Int ALDim = A.LDim();
    const T* ABuffer = A.LockedBuffer();
    T* buffer = this->auxMemory_.Require( sendSize );
    OUTER_PARALLEL_FOR
    for( Int l=0; l<c; ++l )
    {
        const Int thisRowShift = Shift_( l, rowAlignment, c );
        const Int thisLocalWidth = Length_( width, thisRowShift, c );

        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[(k+l*r)*recvSize];
            const Int thisColShift = Shift_(k,colAlignment,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = 
                    &ABuffer[thisColShift+(thisRowShift+jLoc*c)*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }
    }

    // Communicate
    mpi::ReduceScatter( buffer, recvSize, g.VCComm() );

    // Unpack our received data
    T* thisBuffer = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        MemCopy
        ( &thisBuffer[jLoc*thisLDim], 
          &buffer[jLoc*localHeight], localHeight );
    this->auxMemory_.Release();
}

template<typename T>
void
DistMatrix<T,MC,MR>::SumScatterUpdate
( T alpha, const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SumScatterUpdate([MC,* ])");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    if( this->ColAlignment() == A.ColAlignment() )
    {
        if( this->Width() == 1 )
        {
            const Int rowAlignment = this->RowAlignment();
            const Int myCol = g.Col();

            const Int localHeight = this->LocalHeight();
            const Int portionSize = mpi::Pad( localHeight );
            T* buffer = this->auxMemory_.Require( 2*portionSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[portionSize];

            // Pack 
            const T* ACol = A.LockedBuffer(0,0);
            MemCopy( sendBuf, ACol, localHeight );
            
            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuf, recvBuf, portionSize, rowAlignment, g.RowComm() );

            if( myCol == rowAlignment )
            {
                T* thisCol = this->Buffer(0,0);
                FMA_PARALLEL_FOR
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    thisCol[iLoc] += alpha*recvBuf[iLoc];
            }

            this->auxMemory_.Release();
        }
        else
        {
            const Int c = g.Width();
            const Int rowAlignment = this->RowAlignment();

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int maxLocalWidth = MaxLength(width,c);

            const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );
            const Int sendSize = c*portionSize;

            // Pack 
            const Int ALDim = A.LDim();
            const T* ABuffer = A.LockedBuffer();
            T* buffer = this->auxMemory_.Require( sendSize );
            OUTER_PARALLEL_FOR
            for( Int k=0; k<c; ++k )
            {
                T* data = &buffer[k*portionSize];
                const Int thisRowShift = Shift_( k, rowAlignment, c );
                const Int thisLocalWidth = Length_(width,thisRowShift,c);
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* ACol = &ABuffer[(thisRowShift+jLoc*c)*ALDim];
                    T* dataCol = &data[jLoc*localHeight];
                    MemCopy( dataCol, ACol, localHeight );
                }
            }
            
            // Communicate
            mpi::ReduceScatter( buffer, portionSize, g.RowComm() );

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
        if( g.Rank() == 0 )
        {
            std::cerr << "Unaligned SumScatterUpdate [MC,MR] <- [MC,* ]." 
                      << std::endl;
        }
#endif
        if( this->Width() == 1 )
        {
            const Int r = g.Height();
            const Int rowAlignment = this->RowAlignment();
            const Int myRow = g.Row();
            const Int myCol = g.Col();

            const Int height = this->Height();
            const Int localHeight = this->LocalHeight();
            const Int localHeightA = A.LocalHeight();
            const Int maxLocalHeight = MaxLength(height,r);
            const Int portionSize = mpi::Pad( maxLocalHeight );

            const Int colAlignment = this->ColAlignment();
            const Int colAlignmentA = A.ColAlignment();
            const Int sendRow = (myRow+r+colAlignment-colAlignmentA) % r;
            const Int recvRow = (myRow+r+colAlignmentA-colAlignment) % r;

            T* buffer = this->auxMemory_.Require( 2*portionSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[portionSize];

            // Pack 
            const T* ACol = A.LockedBuffer();
            MemCopy( sendBuf, ACol, localHeightA );
            
            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuf, recvBuf, portionSize, rowAlignment, g.RowComm() );

            if( myCol == rowAlignment )
            {
                // Perform the realignment
                mpi::SendRecv
                ( recvBuf, portionSize, sendRow,
                  sendBuf, portionSize, recvRow, g.ColComm() );

                T* thisCol = this->Buffer(0,0);
                FMA_PARALLEL_FOR
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    thisCol[iLoc] += alpha*sendBuf[iLoc];
            }
            this->auxMemory_.Release();
        }
        else
        {
            const Int r = g.Height();
            const Int c = g.Width();
            const Int row = g.Row();

            const Int colAlignment = this->ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int colAlignmentA = A.ColAlignment();
            const Int sendRow = (row+r+colAlignment-colAlignmentA) % r;
            const Int recvRow = (row+r+colAlignmentA-colAlignment) % r;

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localHeightA = A.LocalHeight();
            const Int maxLocalWidth = MaxLength(width,c);

            const Int recvSize_RS = mpi::Pad( localHeightA*maxLocalWidth );
            const Int sendSize_RS = c * recvSize_RS;
            const Int recvSize_SR = localHeight * localWidth;

            T* buffer = this->auxMemory_.Require
            ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );
            T* firstBuf = &buffer[0];
            T* secondBuf = &buffer[recvSize_RS];

            // Pack 
            const T* ABuffer = A.LockedBuffer();
            const Int ALDim = A.LDim();
            OUTER_PARALLEL_FOR
            for( Int k=0; k<c; ++k )
            {
                T* data = &secondBuf[k*recvSize_RS];
                const Int thisRowShift = Shift_( k, rowAlignment, c );
                const Int thisLocalWidth = Length_(width,thisRowShift,c);
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* ACol = &ABuffer[(thisRowShift+jLoc*c)*ALDim];
                    T* dataCol = &data[jLoc*localHeightA];
                    MemCopy( dataCol, ACol, localHeightA );
                }
            }

            // Reduce-scatter over each process row
            mpi::ReduceScatter( secondBuf, firstBuf, recvSize_RS, g.RowComm() );

            // Trade reduced data with the appropriate process row
            mpi::SendRecv
            ( firstBuf,  localHeightA*localWidth, sendRow,
              secondBuf, localHeight*localWidth,  recvRow, g.ColComm() );

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

template<typename T>
void
DistMatrix<T,MC,MR>::SumScatterUpdate( T alpha, const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SumScatterUpdate([* ,MR])");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Width() == 1 && g.Rank() == 0 )
    {
        std::cerr <<
          "The vector version of [MC,MR].SumScatterUpdate([* ,MR]) does not"
          " yet have a vector version implemented, but it would only "
          "require a modification of the vector version of "
          "[MC,MR].SumScatterUpdate([MC,* ])" << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "[MC,MR]::SumScatterUpdate([* ,MR]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MR] matrix instead." << std::endl;
    }
#endif
    if( !this->Participating() )
        return;

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int colAlignment = this->ColAlignment();

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalHeight = MaxLength(height,r);

        const Int recvSize = mpi::Pad( maxLocalHeight*localWidth );
        const Int sendSize = r*recvSize;

        // Pack 
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
        T* buffer = this->auxMemory_.Require( sendSize );
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisColShift = Shift_( k, colAlignment, r );
            const Int thisLocalHeight = Length_(height,thisColShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.ColComm() );

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
        if( g.Rank() == 0 )
        {
            std::cerr << "Unaligned SumScatterUpdate [MC,MR] <- [* ,MR]." 
                      << std::endl;
        }
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int col = g.Col();

        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentA = A.RowAlignment();
        const Int sendCol = (col+c+rowAlignment-rowAlignmentA) % c;
        const Int recvCol = (col+c+rowAlignmentA-rowAlignment) % c;

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localWidthA = A.LocalWidth();
        const Int maxLocalHeight = MaxLength(height,r);

        const Int recvSize_RS = mpi::Pad( maxLocalHeight*localWidthA );
        const Int sendSize_RS = r * recvSize_RS;
        const Int recvSize_SR = localHeight * localWidth;

        T* buffer = this->auxMemory_.Require
        ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[recvSize_RS];

        // Pack
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &secondBuf[k*recvSize_RS];
            const Int thisColShift = Shift_( k, colAlignment, r );
            const Int thisLocalHeight = Length_(height,thisColShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Reduce-scatter over each process col
        mpi::ReduceScatter( secondBuf, firstBuf, recvSize_RS, g.ColComm() );

        // Trade reduced data with the appropriate process col
        mpi::SendRecv
        ( firstBuf,  localHeight*localWidthA, sendCol,
          secondBuf, localHeight*localWidth,  recvCol, g.RowComm() );

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

template<typename T>
void
DistMatrix<T,MC,MR>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SumScatterUpdate([* ,* ])");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int colAlignment = this->ColAlignment();
    const Int rowAlignment = this->RowAlignment();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int maxLocalHeight = MaxLength(height,r);
    const Int maxLocalWidth = MaxLength(width,c);

    const Int recvSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
    const Int sendSize = r * c * recvSize;

    // Pack 
    const T* ABuffer = A.LockedBuffer();
    const Int ALDim = A.LDim();
    T* buffer = this->auxMemory_.Require( sendSize );
    OUTER_PARALLEL_FOR
    for( Int l=0; l<c; ++l )
    {
        const Int thisRowShift = Shift_( l, rowAlignment, c );
        const Int thisLocalWidth = Length_( width, thisRowShift, c );
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[(k+l*r)*recvSize];
            const Int thisColShift = Shift_( k, colAlignment, r );
            const Int thisLocalHeight = Length_(height,thisColShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = 
                    &ABuffer[thisColShift+(thisRowShift+jLoc*c)*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }
    }

    // Communicate
    mpi::ReduceScatter( buffer, recvSize, g.VCComm() );

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

//
// Functions which explicitly work in the complex plane
//

template<typename T>
void
DistMatrix<T,MC,MR>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid(); 
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol*g.Height();
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        const Int jLoc = (j-this->RowShift()) / g.Width();
        this->SetLocalRealPart( iLoc, jLoc, u );
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid(); 
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol*g.Height();
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        const Int jLoc = (j-this->RowShift()) / g.Width();
        this->SetLocalImagPart( iLoc, jLoc, u );
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid(); 
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol*g.Height();
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        const Int jLoc = (j-this->RowShift()) / g.Width();
        this->UpdateLocalRealPart( iLoc, jLoc, u );
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid(); 
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol*g.Height();
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        const Int jLoc = (j-this->RowShift()) / g.Width();
        this->UpdateLocalImagPart( iLoc, jLoc, u );
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::GetRealPartOfDiagonal
( DistMatrix<BASE(T),MD,STAR>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetRealPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    typedef BASE(T) R;
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( this->DistData(), offset );
    }
    const Int length = this->DiagonalLength( offset );
    d.ResizeTo( length, 1 );
    if( !d.Participating() )
        return;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int diagShift = d.ColShift();

    Int iStart, jStart;
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int iLocStart = (iStart-colShift) / r;
    const Int jLocStart = (jStart-rowShift) / c;

    const Int localDiagLength = d.LocalHeight();

    const T* thisBuffer = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    R* dBuf = d.Buffer();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        dBuf[k] = RealPart(thisBuffer[iLoc+jLoc*thisLDim]);
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::GetImagPartOfDiagonal
( DistMatrix<BASE(T),MD,STAR>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetImagPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    typedef BASE(T) R;
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( this->DistData(), offset );
    }
    const Int length = this->DiagonalLength( offset );
    d.ResizeTo( length, 1 );
    if( !d.Participating() )
        return;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int diagShift = d.ColShift();

    Int iStart, jStart;
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int iLocStart = (iStart-colShift) / r;
    const Int jLocStart = (jStart-rowShift) / c;

    const Int localDiagLength = d.LocalHeight();

    const T* thisBuffer = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    R* dBuf = d.Buffer();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        dBuf[k] = ImagPart(thisBuffer[iLoc+jLoc*thisLDim]);
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::GetRealPartOfDiagonal
( DistMatrix<BASE(T),STAR,MD>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetRealPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    typedef BASE(T) R;
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( this->DistData(), offset );
    }
    const Int length = this->DiagonalLength( offset );
    d.ResizeTo( 1, length );
    if( !d.Participating() )
        return;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int diagShift = d.RowShift();

    Int iStart, jStart;
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int iLocStart = (iStart-colShift) / r;
    const Int jLocStart = (jStart-rowShift) / c;

    const Int localDiagLength = d.LocalWidth();

    const T* thisBuffer = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    R* dBuf = d.Buffer();
    const Int dLDim = d.LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        dBuf[k*dLDim] = RealPart(thisBuffer[iLoc+jLoc*thisLDim]);
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::GetImagPartOfDiagonal
( DistMatrix<BASE(T),STAR,MD>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetImagPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    typedef BASE(T) R;
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( this->DistData(), offset );
    }
    const Int length = this->DiagonalLength( offset );
    d.ResizeTo( 1, length );
    if( !d.Participating() )
        return;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int diagShift = d.RowShift();

    Int iStart, jStart;
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int iLocStart = (iStart-colShift) / r;
    const Int jLocStart = (jStart-rowShift) / c;

    const Int localDiagLength = d.LocalWidth();

    const T* thisBuffer = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    R* dBuf = d.Buffer();
    const Int dLDim = d.LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        dBuf[k*dLDim] = ImagPart(thisBuffer[iLoc+jLoc*thisLDim]);
    }
}

template<typename T>
DistMatrix<BASE(T),MD,STAR>
DistMatrix<T,MC,MR>::GetRealPartOfDiagonal( Int offset ) const
{
    DistMatrix<BASE(T),MD,STAR> d( this->Grid() );
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
DistMatrix<BASE(T),MD,STAR>
DistMatrix<T,MC,MR>::GetImagPartOfDiagonal( Int offset ) const
{
    DistMatrix<BASE(T),MD,STAR> d( this->Grid() );
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
void
DistMatrix<T,MC,MR>::SetRealPartOfDiagonal
( const DistMatrix<BASE(T),MD,STAR>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetRealPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    const Int length = this->DiagonalLength( offset );
    if( length != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        LogicError( msg.str() );
    }
#endif
    typedef BASE(T) R;
    if( !d.Participating() )
        return;

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int diagShift = d.ColShift();

    Int iStart,jStart;
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int iLocStart = (iStart-colShift) / r;
    const Int jLocStart = (jStart-rowShift) / c;

    const Int localDiagLength = d.LocalHeight();
    const R* dBuf = d.LockedBuffer();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        this->SetLocalRealPart( iLoc, jLoc, dBuf[k] );
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::SetImagPartOfDiagonal
( const DistMatrix<BASE(T),MD,STAR>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetImagPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    const Int length = this->DiagonalLength( offset );
    if( length != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        LogicError( msg.str() );
    }
#endif
    this->ComplainIfReal();
    typedef BASE(T) R;
    if( !d.Participating() )
        return;

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int diagShift = d.ColShift();

    Int iStart,jStart;
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int iLocStart = (iStart-colShift) / r;
    const Int jLocStart = (jStart-rowShift) / c;

    const Int localDiagLength = d.LocalHeight();
    const R* dBuf = d.LockedBuffer();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        this->SetLocalImagPart( iLoc, jLoc, dBuf[k] );
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::SetRealPartOfDiagonal
( const DistMatrix<BASE(T),STAR,MD>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetRealPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Height() != 1 )
        LogicError("d must be a row vector");
    const Int length = this->DiagonalLength( offset );
    if( length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        LogicError( msg.str() );
    }
#endif
    typedef BASE(T) R;
    if( !d.Participating() )
        return;

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int diagShift = d.RowShift();

    Int iStart,jStart;
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int iLocStart = (iStart-colShift) / r;
    const Int jLocStart = (jStart-rowShift) / c;

    const Int localDiagLength = d.LocalWidth();

    const R* dBuf = d.LockedBuffer();
    const Int dLDim = d.LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        this->SetLocalRealPart( iLoc, jLoc, dBuf[k*dLDim] );
    }
}

template<typename T>
void
DistMatrix<T,MC,MR>::SetImagPartOfDiagonal
( const DistMatrix<BASE(T),STAR,MD>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetImagPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Height() != 1 )
        LogicError("d must be a row vector");
    const Int length = this->DiagonalLength( offset );
    if( length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        LogicError( msg.str() );
    }
#endif
    this->ComplainIfReal();
    typedef BASE(T) R;
    if( !d.Participating() )
        return;

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();
    const Int diagShift = d.RowShift();

    Int iStart,jStart;
    if( offset >= 0 )
    {
        iStart = diagShift;
        jStart = diagShift+offset;
    }
    else
    {
        iStart = diagShift-offset;
        jStart = diagShift;
    }

    const Int iLocStart = (iStart-colShift) / r;
    const Int jLocStart = (jStart-rowShift) / c;

    const Int localDiagLength = d.LocalWidth();
    const R* dBuf = d.LockedBuffer();
    const Int dLDim = d.LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        this->SetLocalImagPart( iLoc, jLoc, dBuf[k*dLDim] );
    }
}

#define PROTO(T) template class DistMatrix<T,MC,MR>
#define COPY(T,CD,RD) \
  template DistMatrix<T,MC,MR>::DistMatrix( const DistMatrix<T,CD,RD>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  STAR); \
  COPY(T,MD,  STAR); \
  COPY(T,MR,  MC  ); \
  COPY(T,MR,  STAR); \
  COPY(T,STAR,MC  ); \
  COPY(T,STAR,MD  ); \
  COPY(T,STAR,MR  ); \
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
