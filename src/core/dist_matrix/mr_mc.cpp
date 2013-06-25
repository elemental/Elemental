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
DistMatrix<T,MR,MC,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->SetShifts(); }

template<typename T,typename Int>
DistMatrix<T,MR,MC,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T,typename Int>
DistMatrix<T,MR,MC,Int>::DistMatrix
( Int height, Int width,
  Int colAlignment, Int rowAlignment, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->Align(colAlignment,rowAlignment); this->ResizeTo(height,width); }

template<typename T,typename Int>
DistMatrix<T,MR,MC,Int>::DistMatrix
( Int height, Int width,
  Int colAlignment, Int rowAlignment, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->Align(colAlignment,rowAlignment); this->ResizeTo(height,width,ldim); }

template<typename T,typename Int>
DistMatrix<T,MR,MC,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  const T* buffer, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->LockedAttach(height,width,colAlignment,rowAlignment,buffer,ldim,g); }

template<typename T,typename Int>
DistMatrix<T,MR,MC,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  T* buffer, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->Attach(height,width,colAlignment,rowAlignment,buffer,ldim,g); }

template<typename T,typename Int>
DistMatrix<T,MR,MC,Int>::DistMatrix( const DistMatrix<T,MR,MC,Int>& A )
: AbstractDistMatrix<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MR,MC]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [MR,MC] with itself");
}

template<typename T,typename Int>
template<Distribution U,Distribution V>
DistMatrix<T,MR,MC,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MR,MC]::DistMatrix");
#endif
    this->SetShifts();
    if( MR != U || MC != V || 
        reinterpret_cast<const DistMatrix<T,MR,MC,Int>*>(&A) != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [MR,MC] with itself");
}

template<typename T,typename Int>
DistMatrix<T,MR,MC,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
elem::DistData<Int>
DistMatrix<T,MR,MC,Int>::DistData() const
{
    elem::DistData<Int> data;
    data.colDist = MR;
    data.rowDist = MC;
    data.colAlignment = this->colAlignment_;
    data.rowAlignment = this->rowAlignment_;
    data.root = 0;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T,typename Int>
Int
DistMatrix<T,MR,MC,Int>::ColStride() const
{ return this->grid_->Width(); }

template<typename T,typename Int>
Int
DistMatrix<T,MR,MC,Int>::RowStride() const
{ return this->grid_->Height(); }

template<typename T,typename Int>
Int
DistMatrix<T,MR,MC,Int>::ColRank() const
{ return this->grid_->Col(); }

template<typename T,typename Int>
Int
DistMatrix<T,MR,MC,Int>::RowRank() const
{ return this->grid_->Row(); }

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::AlignWith( const elem::DistData<Int>& data )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::AlignWith");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( data.colDist == MC && data.rowDist == MR )
    {
        this->colAlignment_ = data.rowAlignment;
        this->rowAlignment_ = data.colAlignment;
        this->constrainedColAlignment_ = true;
        this->constrainedRowAlignment_ = true;
    }
    else if( data.colDist == MC && data.rowDist == STAR )
    {
        this->rowAlignment_ = data.colAlignment;
        this->constrainedRowAlignment_ = true;
    }
    else if( data.colDist == MR && data.rowDist == MC )
    {
        this->colAlignment_ = data.colAlignment;
        this->rowAlignment_ = data.rowAlignment;
        this->constrainedColAlignment_ = true;
        this->constrainedRowAlignment_ = true;
    }
    else if( data.colDist == MR && data.rowDist == STAR )
    {
        this->colAlignment_ = data.colAlignment;
        this->constrainedColAlignment_ = true;
    }
    else if( data.colDist == STAR && data.rowDist == MC )
    {
        this->rowAlignment_ = data.rowAlignment;
        this->constrainedRowAlignment_ = true;
    }
    else if( data.colDist == STAR && data.rowDist == MR )
    {
        this->colAlignment_ = data.rowAlignment;
        this->constrainedColAlignment_ = true;
    }
    else if( data.colDist == STAR && data.rowDist == VC )
    {
        this->rowAlignment_ = data.rowAlignment % this->RowStride();
        this->constrainedRowAlignment_ = true;
    }
    else if( data.colDist == STAR && data.rowDist == VR )
    {
        this->colAlignment_ = data.rowAlignment % this->ColStride();
        this->constrainedColAlignment_ = true;
    }
    else if( data.colDist == VC && data.rowDist == STAR )
    {
        this->rowAlignment_ = data.colAlignment % this->RowStride();
        this->constrainedRowAlignment_ = true;
    }
    else if( data.colDist == VR && data.rowDist == STAR )
    {
        this->colAlignment_ = data.colAlignment % this->ColStride();
        this->constrainedColAlignment_ = true;
    }
#ifndef RELEASE
    else throw std::logic_error("Nonsensical alignment");
#endif
    this->SetShifts();
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::AlignWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::AlignColsWith( const elem::DistData<Int>& data )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::AlignColsWith");
    this->AssertFreeColAlignment();
    if( *this->grid_ != *data.grid )
        throw std::logic_error("Grids do not match");
#endif
    if( data.colDist == MR )
        this->colAlignment_ = data.colAlignment;
    else if( data.rowDist == MR )
        this->colAlignment_ = data.rowAlignment;
    else if( data.colDist == VR )
        this->colAlignment_ = data.colAlignment % this->ColStride();
    else if( data.rowDist == VR )
        this->colAlignment_ = data.rowAlignment % this->ColStride();
#ifndef RELEASE
    else throw std::logic_error("Nonsensical alignment");
#endif
    this->constrainedColAlignment_ = true;
    this->SetShifts();
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::AlignColsWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignColsWith( A.DistData() ); }

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::AlignRowsWith( const elem::DistData<Int>& data )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::AlignRowsWith");
    this->AssertFreeRowAlignment();
    if( *this->grid_ != *data.grid )
        throw std::logic_error("Grids do not match");
#endif
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
DistMatrix<T,MR,MC,Int>::AlignRowsWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignRowsWith( A.DistData() ); }

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::Attach
( Int height, Int width, Int colAlignment, Int rowAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::Attach");
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
        const Int localHeight = Length(height,this->colShift_,g.Width());
        const Int localWidth = Length(width,this->rowShift_,g.Height());
        this->LocalAttach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::LockedAttach
( Int height, Int width, Int colAlignment, Int rowAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::LockedAttach");
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
        const Int localHeight = Length(height,this->colShift_,g.Width());
        const Int localWidth = Length(width,this->rowShift_,g.Height());
        this->LocalLockedAttach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    const elem::Grid& g = this->Grid();
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->LocalResize_
        ( Length( height, this->ColShift(), g.Width() ),
          Length( width,  this->RowShift(), g.Height() ) );
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    const elem::Grid& g = this->Grid();
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->LocalResize_
        ( Length( height, this->ColShift(), g.Width() ),
          Length( width,  this->RowShift(), g.Height() ), ldim );
}

template<typename T,typename Int>
T
DistMatrix<T,MR,MC,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol * g.Height();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        const Int jLoc = (j-this->RowShift()) / g.Height();
        u = this->GetLocal(iLoc,jLoc);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
    return u;
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->SetLocal(iLoc,jLoc,u);
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateLocal(iLoc,jLoc,u);
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::GetDiagonal
( DistMatrix<T,MD,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::GetDiagonal([MD,* ])");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Height() || d.Width() != 1) )
        throw std::logic_error("d is not of the correct dimensions");
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( diagLength, 1 );
    }

    if( d.Participating() )
    {
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

        const Int iLocStart = (iStart-colShift) / c;
        const Int jLocStart = (jStart-rowShift) / r;

        const Int localDiagLength = d.LocalHeight();

        const T* thisBuf = this->LockedBuffer();
        const Int thisLDim = this->LDim();
        T* dBuf = d.Buffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLoc = iLocStart + k*(lcm/c);
            const Int jLoc = jLocStart + k*(lcm/r);
            dBuf[k] = thisBuf[iLoc+jLoc*thisLDim]; 
        }
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::GetDiagonal
( DistMatrix<T,STAR,MD,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::GetDiagonal([* ,MD])");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Width() || d.Height() != 1) )
        throw std::logic_error("d is not of the correct dimensions");
    if( ( d.Viewing() && d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( 1, diagLength );
    }

    if( d.Participating() )
    {
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

        const Int iLocStart = (iStart-colShift) / c;
        const Int jLocStart = (jStart-rowShift) / r;

        const Int localDiagLength = d.LocalWidth();

        const T* thisBuf = this->LockedBuffer();
        const Int thisLDim = this->LDim();
        T* dBuf = d.Buffer();
        const Int dLDim = d.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLoc = iLocStart + k*(lcm/c);
            const Int jLoc = jLocStart + k*(lcm/r);
            dBuf[k*dLDim] = thisBuf[iLoc+jLoc*thisLDim];
        }
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SetDiagonal
( const DistMatrix<T,MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SetDiagonal([MD,* ])");
    this->AssertSameGrid( d.Grid() );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const Int diagLength = this->DiagonalLength(offset);
    if( diagLength != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    if( d.Participating() )
    {
        const elem::Grid& g = this->Grid();
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

        const Int iLocStart = (iStart-colShift) / c;
        const Int jLocStart = (jStart-rowShift) / r;

        const Int localDiagLength = d.LocalHeight();

        const T* dBuf = d.LockedBuffer();
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLoc = iLocStart + k*(lcm/c);
            const Int jLoc = jLocStart + k*(lcm/r);
            thisBuf[iLoc+jLoc*thisLDim] = dBuf[k];
        }
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SetDiagonal
( const DistMatrix<T,STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SetDiagonal([* ,MD])");
    this->AssertSameGrid( d.Grid() );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const Int diagLength = this->DiagonalLength(offset);
    if( diagLength != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    if( d.Participating() )
    {
        const elem::Grid& g = this->Grid();
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

        const Int iLocStart = (iStart-colShift) / c;
        const Int jLocStart = (jStart-rowShift) / r;

        const Int localDiagLength = d.LocalWidth();

        const T* dBuf = d.LockedBuffer();
        const Int dLDim = d.LDim();
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLoc = iLocStart + k*(lcm/c);
            const Int jLoc = jLocStart + k*(lcm/r);
            thisBuf[iLoc+jLoc*thisLDim] = dBuf[k*dLDim];
        }
    }
}

//
// Utility functions, e.g., operator=
//

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        ResizeTo( A.Height(), A.Width() );

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
        const Int ownerRow = this->RowAlignment();
        const Int ownerCol = A.RowAlignment();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int height = A.Height();
        const Int maxLocalHeight = MaxLength(height,p);
        const Int portionSize = mpi::Pad( maxLocalHeight );

        const Int colShiftVR = Shift(rankRM,colAlignment,p);
        const Int colShiftVCOfA = Shift(rankCM,colAlignmentOfA,p);
        const Int sendRankRM = (rankRM+(p+colShiftVCOfA-colShiftVR)) % p;
        const Int recvRankCM = (rankCM+(p+colShiftVR-colShiftVCOfA)) % p;
        const Int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

        T* buffer = this->auxMemory_.Require( (r+c)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        if( myCol == ownerCol )
        {
            // Pack
            const Int AColShift = A.ColShift();
            const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const Int shift = Shift_(myRow+r*k,colAlignmentOfA,p);
                const Int offset = (shift-AColShift) / r;
                const Int thisLocalHeight = Length_(height,shift,p);

                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    data[iLoc] = ABuf[offset+iLoc*c];
            }
        }

        // A[VC,* ] <- A[MC,MR]
        mpi::Scatter
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerCol, g.RowComm() );

        // A[VR,* ] <- A[VC,* ]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankRM, 0,
          recvBuf, portionSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

        // A[MR,MC] <- A[VR,* ]
        mpi::Gather
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerRow, g.ColComm() );

        if( myRow == ownerRow )
        {
            // Unpack
            const Int thisColShift = this->ColShift();
            T* thisBuf = this->Buffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const Int shift = Shift_(myCol+c*k,colAlignment,p);
                const Int offset = (shift-thisColShift) / c;
                const Int thisLocalHeight = Length_(height,shift,p);

                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    thisBuf[offset+iLoc*r] = data[iLoc];
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
        const Int ownerCol = this->ColAlignment();
        const Int ownerRow = A.ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int width = A.Width();
        const Int maxLocalWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxLocalWidth );

        const Int rowShiftVC = Shift(rankCM,rowAlignment,p);
        const Int rowShiftVROfA = Shift(rankRM,rowAlignmentOfA,p);
        const Int sendRankCM = (rankCM+(p+rowShiftVROfA-rowShiftVC)) % p;
        const Int recvRankRM = (rankRM+(p+rowShiftVC-rowShiftVROfA)) % p;
        const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        T* buffer = this->auxMemory_.Require( (r+c)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        if( myRow == ownerRow )
        {
            // Pack
            const Int ARowShift = A.RowShift();
            const T* ABuf = A.LockedBuffer();
            const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const Int shift = Shift_(myCol+c*k,rowAlignmentOfA,p);
                const Int offset = (shift-ARowShift) / c;
                const Int thisLocalWidth = Length_(width,shift,p);

                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                    data[jLoc] = ABuf[(offset+jLoc*r)*ALDim];
            }
        }

        // A[* ,VR] <- A[MC,MR]
        mpi::Scatter
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerRow, g.ColComm() );

        // A[* ,VC] <- A[* ,VR]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankCM, 0,
          recvBuf, portionSize, recvRankCM, mpi::ANY_TAG, g.VCComm() );

        // A[MR,MC] <- A[* ,VC]
        mpi::Gather
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerCol, g.RowComm() );

        if( myCol == ownerCol )
        {
            // Unpack
            const Int thisRowShift = this->RowShift();
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const Int shift = Shift_(myRow+r*k,rowAlignment,p);
                const Int offset = (shift-thisRowShift) / r;
                const Int thisLocalWidth = Length_(width,shift,p);

                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                    thisBuf[(offset+jLoc*c)*thisLDim] = data[jLoc];
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
        if( A.Height() >= A.Width() )
        {
            std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
            ( new DistMatrix<T,VC,STAR,Int>(g) );
            *A_VC_STAR = A;

            std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
            ( new DistMatrix<T,VR,STAR,Int>(true,this->ColAlignment(),g) );
            *A_VR_STAR = *A_VC_STAR;
            delete A_VC_STAR.release(); // lowers memory highwater

            *this = *A_VR_STAR;
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
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(g) );
    *A_VC_STAR = A;

    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VR_STAR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_VR_STAR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [* ,MR]");
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
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MR,MC] = [MD,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MR,MC] = [* ,MD] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            if( this->Participating() )
                this->colShift_ = A.ColShift();
        }
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

    if( this->ColAlignment() == A.ColAlignment() &&
        this->RowAlignment() == A.RowAlignment() )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MR,MC] <- [MR,MC]" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.Row();
        const Int col = g.Col();

        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const Int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
        const Int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;
        const Int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;
        const Int sendRank = sendRow + sendCol*r;
        const Int recvRank = recvRow + recvCol*r;

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int sendSize = localHeightOfA * localWidthOfA;
        const Int recvSize = localHeight * localWidth;

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
            T* sendBufCol = &sendBuf[jLoc*localHeightOfA];
            MemCopy( sendBufCol, ACol, localHeightOfA );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank, 0,
          recvBuf, recvSize, recvRank, mpi::ANY_TAG, g.VCComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
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
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            this->SetColShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !this->Participating() )
        return *this;

    if( this->ColAlignment() == A.ColAlignment() )
    {
        const Int r = g.Height();
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
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MR,MC] <- [MR,* ]" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int col = g.Col();

        const Int rowShift = this->RowShift();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int sendRank = (col+c+colAlignment-colAlignmentOfA) % c;
        const Int recvRank = (col+c+colAlignmentOfA-colAlignment) % c;

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int sendSize = localHeightOfA * localWidth;
        const Int recvSize = localHeight * localWidth;

        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* ACol = &ABuf[(rowShift+jLoc*r)*ALDim];
            T* sendBufCol = &sendBuf[jLoc*localHeightOfA];
            MemCopy( sendBufCol, ACol, localHeightOfA );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank, 0,
          recvBuf, recvSize, recvRank, mpi::ANY_TAG, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
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
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [* ,MC]");
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
    if( !this->Participating() )
        return *this;

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const Int c = g.Width();
        const Int colShift = this->ColShift();

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
            T* destCol = &thisBuf[jLoc*thisLDim];
            const T* sourceCol = &ABuf[colShift+jLoc*ALDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*c];
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MR,MC] <- [* ,MC]" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.Row(); 

        const Int colShift = this->ColShift();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfA = A.LocalWidth();

        const Int sendSize = localHeight * localWidthOfA;
        const Int recvSize = localHeight * localWidth;

        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            T* destCol = &sendBuf[jLoc*localHeight];
            const T* sourceCol = &ABuf[colShift+jLoc];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*c];
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRow, 0,
          recvBuf, recvSize, recvRow, mpi::ANY_TAG, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
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
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VR,STAR,Int> A_VR_STAR(g);

    A_VR_STAR = A;
    *this = A_VR_STAR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [* ,VC]");
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
        const Int p = r * c;
        const Int row = g.Row();

        const Int rowShift = this->RowShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLength(height,c);
        const Int maxWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*c*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuf[k*portionSize];
            const Int thisColShift = Shift_(k,colAlignment,c);
            const Int thisLocalHeight = Length_(height,thisColShift,c);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[k*portionSize];

            const Int thisRank = row+k*r;
            const Int thisRowShift = Shift_(thisRank,rowAlignmentOfA,p);
            const Int thisRowOffset = (thisRowShift-rowShift) / r;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuf[(thisRowOffset+jLoc*c)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MR,MC] <- [* ,VC]" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int row = g.Row();

        const Int rowShift = this->RowShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRow = (row+r+rowAlignment-(rowAlignmentOfA%r)) % r;
        const Int recvRow = (row+r+(rowAlignmentOfA%r)-rowAlignment) % r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLength(height,c);
        const Int maxWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*c*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[c*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int thisColShift = Shift_(k,colAlignment,c);
            const Int thisLocalHeight = Length_(height,thisColShift,c);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // SendRecv to align A[* ,VC] via a trade in the column
        mpi::SendRecv
        ( secondBuf, c*portionSize, sendRow, 0,
          firstBuf,  c*portionSize, recvRow, mpi::ANY_TAG, g.ColComm() );

        // AllToAll to gather all of the aligned [* ,VC] into secondBuf
        mpi::AllToAll
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
            const Int thisRank = recvRow+k*r;
            const Int thisRowShift = Shift_(thisRank,rowAlignmentOfA,p);
            const Int thisRowOffset = (thisRowShift-rowShift) / r;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuf[(thisRowOffset+jLoc*c)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment() % g.Width();
            this->SetColShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !this->Participating() )
        return *this;

    if( this->ColAlignment() == A.ColAlignment() % g.Width() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int col = g.Col();

        const Int colShift = this->ColShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int maxHeight = MaxLength(height,p);
        const Int maxWidth = MaxLength(width,r);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*r*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &sendBuf[k*portionSize];
            const Int thisRowShift = Shift_(k,rowAlignment,r);
            const Int thisLocalWidth = Length_(width,thisRowShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* ACol = &ABuf[(thisRowShift+jLoc*r)*ALDim];
                T* dataCol = &data[jLoc*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisRank = col+k*c;
            const Int thisColShift = Shift_(thisRank,colAlignmentOfA,p);
            const Int thisColOffset = (thisColShift-colShift) / c;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[thisColOffset+jLoc*thisLDim];
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
            std::cerr << "Unaligned [MR,MC] <- [* ,VC]" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int col = g.Col();

        const Int colShift = this->ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int sendCol = (col+c+colAlignment-(colAlignmentOfA%c)) % c;
        const Int recvCol = (col+c+(colAlignmentOfA%c)-colAlignment) % c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int maxHeight = MaxLength(height,p);
        const Int maxWidth = MaxLength(width,r);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*r*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[r*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int thisRowShift = Shift_(k,rowAlignment,r);
            const Int thisLocalWidth = Length_(width,thisRowShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* ACol = &ABuf[(thisRowShift+jLoc*r)*ALDim];
                T* dataCol = &data[jLoc*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // SendRecv to align A[VR,* ] via a trade in the row
        mpi::SendRecv
        ( secondBuf, r*portionSize, sendCol, 0,
          firstBuf,  r*portionSize, recvCol, mpi::ANY_TAG, g.RowComm() );

        // AllToAll to gather all of the aligned [VR,* ] data into secondBuf
        mpi::AllToAll
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int thisRank = recvCol+k*c;
            const Int thisColShift = Shift_(thisRank,colAlignmentOfA,p);
            const Int thisColOffset = (thisColShift-colShift) / c;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &thisBuf[thisColOffset+jLoc*thisLDim];
                const T* sourceCol = &data[jLoc*thisLocalHeight];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc*r] = sourceCol[iLoc];
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VC,Int> A_STAR_VC(g);

    A_STAR_VC = A;
    *this = A_STAR_VC;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();

    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        T* destCol = &thisBuf[jLoc*thisLDim];
        const T* sourceCol = &ABuf[colShift+(rowShift+jLoc*r)*ALDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            destCol[iLoc] = sourceCol[iLoc*c];
    }
    return *this;
}

// NOTE: This is almost an exact duplicate of [MC,MR] <- [o, o ]
template<typename T,typename Int>
const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,CIRC,CIRC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC] = [o ,o ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const int colStride = this->ColStride();
    const int rowStride = this->RowStride();
    const int p = g.Size();
    if( !this->Viewing() )
        this->ResizeTo( m, n );

    const int colAlignment = this->ColAlignment();
    const int rowAlignment = this->RowAlignment();
    const int mLocal = this->LocalHeight();
    const int nLocal = this->LocalWidth();
    const int pkgSize = mpi::Pad(MaxLength(m,colStride)*MaxLength(n,rowStride));
    const int recvSize = pkgSize;
    const int sendSize = p*pkgSize;
    T* recvBuf;
    if( A.Participating() )
    {
        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        recvBuf = &buffer[sendSize];

        // Pack the send buffer
        const int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        for( int t=0; t<rowStride; ++t )
        {
            const int tLocalWidth = Length( n, t, rowStride );
            // NOTE: switched vs. [MC,MR] variant of [o, o] redist
            const int row = (rowAlignment+t) % rowStride;
            for( int s=0; s<colStride; ++s )
            {
                const int sLocalHeight = Length( m, s, colStride );
                // NOTE: switched vs. [MC,MR] variant of [o, o] redist
                const int col = (colAlignment+s) % colStride;
                const int q = row + col*colStride;
                for( int jLoc=0; jLoc<tLocalWidth; ++jLoc )
                {
                    const int j = t + jLoc*rowStride;
                    for( int iLoc=0; iLoc<sLocalHeight; ++iLoc )
                    {
                        const int i = s + iLoc*colStride;
                        sendBuf[q*pkgSize+iLoc+jLoc*sLocalHeight] =
                            ABuf[i+j*ALDim];
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
          recvBuf,         pkgSize, A.Root(), g.VCComm() );
    }

    if( this->Participating() )
    {
        // Unpack
        const int ldim = this->LDim();
        T* buffer = this->Buffer();
        for( int jLoc=0; jLoc<nLocal; ++jLoc )
            for( int iLoc=0; iLoc<mLocal; ++iLoc )
                buffer[iLoc+jLoc*ldim] = recvBuf[iLoc+jLoc*mLocal];
        this->auxMemory_.Release();
    }

    return *this;
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SumScatterFrom( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SumScatterFrom([MR,* ])");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            this->SetColShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !this->Participating() )
        return;

    if( this->ColAlignment() == A.ColAlignment() )
    {
        if( this->Width() == 1 )
        {
            const Int rowAlignment = this->RowAlignment();
            const Int myRow = g.Row();

            const Int localHeight = this->LocalHeight();
            const Int portionSize = mpi::Pad( localHeight );

            T* buffer = this->auxMemory_.Require( 2*portionSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[portionSize];

            // Pack 
            const T* ACol = A.LockedBuffer();
            MemCopy( sendBuf, ACol, localHeight );

            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuf, recvBuf, portionSize,
              mpi::SUM, rowAlignment, g.ColComm() );

            if( myRow == rowAlignment )
            {
                T* thisCol = this->Buffer();
                MemCopy( thisCol, recvBuf, localHeight );
            }

            this->auxMemory_.Release();
        }
        else
        {
            const Int r = g.Height();
            const Int rowAlignment = this->RowAlignment();

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int maxLocalWidth = MaxLength(width,r);
            const Int recvSize = mpi::Pad( localHeight*maxLocalWidth );
            const Int sendSize = r * recvSize;

            T* buffer = this->auxMemory_.Require( sendSize );

            // Pack 
            const Int ALDim = A.LDim();
            const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &buffer[k*recvSize];
                const Int thisRowShift = Shift_( k, rowAlignment, r );
                const Int thisLocalWidth = Length_( width, thisRowShift, r );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* ACol = &ABuf[(thisRowShift+jLoc*r)*ALDim];
                    T* dataCol = &data[jLoc*localHeight];
                    MemCopy( dataCol, ACol, localHeight );
                }
            }

            // Communicate
            mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.ColComm() );

            // Unpack our received data
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int j=0; j<localWidth; ++j )
            {
                const T* bufferCol = &buffer[j*localHeight];
                T* thisCol = &thisBuf[j*thisLDim];
                MemCopy( thisCol, bufferCol, localHeight );
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned SumScatterFrom [MR,MC] <- [MR,* ]" 
                      << std::endl;
#endif
        if( this->Width() == 1 )
        {
            const Int c = g.Width();
            const Int rowAlignment = this->RowAlignment();
            const Int myRow = g.Row();
            const Int myCol = g.Col();

            const Int height = this->Height();
            const Int localHeight = this->LocalHeight();
            const Int localHeightOfA = A.LocalHeight();
            const Int maxLocalHeight = MaxLength(height,c);
            const Int portionSize = mpi::Pad( maxLocalHeight );

            const Int colAlignment = this->ColAlignment();
            const Int colAlignmentOfA = A.ColAlignment();
            const Int sendCol = (myCol+c+colAlignment-colAlignmentOfA) % c;
            const Int recvCol = (myCol+c+colAlignmentOfA-colAlignment) % c;

            T* buffer = this->auxMemory_.Require( 2*portionSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[portionSize];

            // Pack
            const T* ACol = A.LockedBuffer();
            MemCopy( sendBuf, ACol, localHeightOfA );

            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuf, recvBuf, portionSize,
              mpi::SUM, rowAlignment, g.ColComm() );

            if( myRow == rowAlignment )
            {
                // Perform the realignment
                mpi::SendRecv
                ( recvBuf, portionSize, sendCol, 0,
                  sendBuf, portionSize, recvCol, 0, g.RowComm() );

                T* thisCol = this->Buffer();
                MemCopy( thisCol, sendBuf, localHeight );
            }

            this->auxMemory_.Release();
        }
        else
        {
            const Int r = g.Height();
            const Int c = g.Width();
            const Int col = g.Col();

            const Int colAlignment = this->ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int colAlignmentOfA = A.ColAlignment();
            const Int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
            const Int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localHeightOfA = A.LocalHeight();
            const Int maxLocalWidth = MaxLength(width,r);

            const Int recvSize_RS = mpi::Pad( localHeightOfA*maxLocalWidth );
            const Int sendSize_RS = r * recvSize_RS;
            const Int recvSize_SR = localHeight * localWidth;

            T* buffer = this->auxMemory_.Require
            ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );
            T* firstBuf = &buffer[0];
            T* secondBuf = &buffer[recvSize_RS];

            // Pack
            const T* ABuf = A.LockedBuffer();
            const Int ALDim = A.LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &secondBuf[k*recvSize_RS];

                const Int thisRowShift = Shift_( k, rowAlignment, r );
                const Int thisLocalWidth = Length_(width,thisRowShift,r);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* ACol = &ABuf[(thisRowShift+jLoc*r)*ALDim];
                    T* dataCol = &data[jLoc*localHeightOfA];
                    MemCopy( dataCol, ACol, localHeightOfA );
                }
            }

            // Reduce-scatter over each process col
            mpi::ReduceScatter
            ( secondBuf, firstBuf, recvSize_RS, mpi::SUM, g.ColComm() );

            // Trade reduced data with the appropriate process col
            mpi::SendRecv
            ( firstBuf,  localHeightOfA*localWidth, sendCol, 0,
              secondBuf, localHeight*localWidth,    recvCol, mpi::ANY_TAG,
              g.RowComm() );

            // Unpack the received data
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* secondBufCol = &secondBuf[jLoc*localHeight];
                T* thisCol = &thisBuf[jLoc*thisLDim];
                MemCopy( thisCol, secondBufCol, localHeight );
            }
            this->auxMemory_.Release();
        }
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SumScatterFrom( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SumScatterFrom([* ,MC])");
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
          "The vector version of [MR,MC].SumScatterFrom([* ,MC]) is not yet"
          " written, but it only requires a modification of the vector "
          "version of [MR,MC].SumScatterFrom([MR,* ])" << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "[MR,MC]::SumScatterFrom([* ,MC]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MC] matrix instead." << std::endl;
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
        return;

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const Int c = g.Width();
        const Int colAlignment = this->ColAlignment();

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalHeight = MaxLength(height,c);
        const Int recvSize = mpi::Pad( maxLocalHeight*localWidth );
        const Int sendSize = c * recvSize;

        T* buffer = this->auxMemory_.Require( sendSize );

        // Pack 
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisColShift = Shift_( k, colAlignment, c );
            const Int thisLocalHeight = Length_( height, thisColShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.RowComm() );

        // Unpack our received data
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* bufferCol = &buffer[jLoc*localHeight];
            T* thisCol = &thisBuf[jLoc*thisLDim];
            MemCopy( thisCol, bufferCol, localHeight );
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned SumScatterFrom [MR,MC] <- [* ,MC]" 
                      << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.Row();

        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalHeight = MaxLength(height,c);
        const Int recvSize_RS = mpi::Pad( maxLocalHeight*localWidthOfA );
        const Int sendSize_RS = c* recvSize_RS;
        const Int recvSize_SR = localHeight * localWidth;

        T* buffer = this->auxMemory_.Require
        ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[recvSize_RS];

        // Pack 
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &secondBuf[k*recvSize_RS];
            const Int thisColShift = Shift_( k, colAlignment, c );
            const Int thisLocalHeight = Length_( height, thisColShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // Reduce-scatter over each process row
        mpi::ReduceScatter
        ( secondBuf, firstBuf, recvSize_RS, mpi::SUM, g.RowComm() );

        // Trade reduced data with the appropriate process row
        mpi::SendRecv
        ( firstBuf,  localHeight*localWidthOfA, sendRow, 0,
          secondBuf, localHeight*localWidth,    recvRow, mpi::ANY_TAG, 
          g.ColComm() );

        // Unpack the received data
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* secondBufCol = &secondBuf[jLoc*localHeight];
            T* thisCol = &thisBuf[jLoc*thisLDim];
            MemCopy( thisCol, secondBufCol, localHeight );
        }
        this->auxMemory_.Release();
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SumScatterFrom( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SumScatterFrom([* ,* ])");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

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
    const Int maxLocalHeight = MaxLength(height,c);
    const Int maxLocalWidth = MaxLength(width,r);
    const Int recvSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
    const Int sendSize = r * c * recvSize;

    T* buffer = this->auxMemory_.Require( sendSize );

    // Pack 
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int l=0; l<r; ++l )
    {
        const Int thisRowShift = Shift_( l, rowAlignment, r );
        const Int thisLocalWidth = Length_( width, thisRowShift, r );
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[(k+l*c)*recvSize];
            const Int thisColShift = Shift_( k, colAlignment, c );
            const Int thisLocalHeight = Length_( height, thisColShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = 
                    &ABuf[thisColShift+(thisRowShift+jLoc*r)*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }
    }

    // Communicate
    mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.VRComm() );

    // Unpack our received data
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* bufferCol = &buffer[jLoc*localHeight];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        MemCopy( thisCol, bufferCol, localHeight );
    }
    this->auxMemory_.Release();
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SumScatterUpdate([MR,* ])");
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
            const Int myRow = g.Row();

            const Int localHeight = this->LocalHeight();
            const Int portionSize = mpi::Pad( localHeight );

            T* buffer = this->auxMemory_.Require( 2*portionSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[portionSize];

            // Pack
            const T* ACol = A.LockedBuffer();
            MemCopy( sendBuf, ACol, localHeight );

            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuf, recvBuf, portionSize,
              mpi::SUM, rowAlignment, g.ColComm() );

            if( myRow == rowAlignment )
            {
                T* thisCol = this->Buffer();
#if defined(HAVE_OPENMP) && !defined(AVOID_OMP_FMA)
#pragma omp parallel for
#endif
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    thisCol[iLoc] += alpha*recvBuf[iLoc];
            }
            this->auxMemory_.Release();
        }
        else
        {
            const Int r = g.Height();
            const Int rowAlignment = this->RowAlignment();

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int maxLocalWidth = MaxLength(width,r);
            const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );
            const Int sendSize = r*portionSize;

            T* buffer = this->auxMemory_.Require( sendSize );

            // Pack 
            const Int ALDim = A.LDim();
            const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &buffer[k*portionSize];
                const Int thisRowShift = Shift_( k, rowAlignment, r );
                const Int thisLocalWidth = Length_(width,thisRowShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* ACol = &ABuf[(thisRowShift+jLoc*r)*ALDim];
                    T* dataCol = &data[jLoc*localHeight];
                    MemCopy( dataCol, ACol, localHeight );
                }
            }

            // Communicate
            mpi::ReduceScatter( buffer, portionSize, mpi::SUM, g.ColComm() );

            // Update with our received data
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(AVOID_OMP_FMA)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* bufferCol = &buffer[jLoc*localHeight];
                T* thisCol = &thisBuf[jLoc*thisLDim];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    thisCol[iLoc] += alpha*bufferCol[iLoc];
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned SumScatterUpdate [MR,MC] <- [MR,* ]" 
                      << std::endl;
#endif
        if( this->Width() == 1 )
        {
            const Int c = g.Width();
            const Int rowAlignment = this->RowAlignment();
            const Int myRow = g.Row();
            const Int myCol = g.Col();

            const Int height = this->Height();
            const Int localHeight = this->LocalHeight();
            const Int localHeightOfA = A.LocalHeight();
            const Int maxLocalHeight = MaxLength(height,c);
            const Int portionSize = mpi::Pad( maxLocalHeight );

            const Int colAlignment = this->ColAlignment();
            const Int colAlignmentOfA = A.ColAlignment();
            const Int sendCol = (myCol+c+colAlignment-colAlignmentOfA) % c;
            const Int recvCol = (myCol+c+colAlignmentOfA-colAlignment) % c;

            T* buffer = this->auxMemory_.Require( 2*portionSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[portionSize];

            // Pack
            const T* ACol = A.LockedBuffer();
            MemCopy( sendBuf, ACol, localHeightOfA );

            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuf, recvBuf, portionSize,
              mpi::SUM, rowAlignment, g.ColComm() );

            if( myRow == rowAlignment )
            {
                // Perform the realignment
                mpi::SendRecv
                ( recvBuf, portionSize, sendCol, 0,
                  sendBuf, portionSize, recvCol, 0, g.RowComm() );

                T* thisCol = this->Buffer();
#if defined(HAVE_OPENMP) && !defined(AVOID_OMP_FMA)
#pragma omp parallel for
#endif
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    thisCol[iLoc] += alpha*sendBuf[iLoc];
            }
            this->auxMemory_.Release();
        }
        else
        {
            const Int r = g.Height();
            const Int c = g.Width();
            const Int col = g.Col();

            const Int colAlignment = this->ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int colAlignmentOfA = A.ColAlignment();
            const Int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
            const Int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localHeightOfA = A.LocalHeight();
            const Int maxLocalWidth = MaxLength(width,r);
            const Int recvSize_RS = mpi::Pad( localHeightOfA*maxLocalWidth );
            const Int sendSize_RS = r * recvSize_RS;
            const Int recvSize_SR = localHeight * localWidth;

            T* buffer = this->auxMemory_.Require
            ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );
            T* firstBuf = &buffer[0];
            T* secondBuf = &buffer[recvSize_RS];

            // Pack 
            const Int ALDim = A.LDim();
            const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &secondBuf[k*recvSize_RS];
                const Int thisRowShift = Shift_( k, rowAlignment, r );
                const Int thisLocalWidth = Length_(width,thisRowShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* ACol = &ABuf[(thisRowShift+jLoc*r)*ALDim];
                    T* dataCol = &data[jLoc*localHeightOfA];
                    MemCopy( dataCol, ACol, localHeightOfA );
                }
            }

            // Reduce-scatter over each process col
            mpi::ReduceScatter
            ( secondBuf, firstBuf, recvSize_RS, mpi::SUM, g.ColComm() );

            // Trade reduced data with the appropriate process col
            mpi::SendRecv
            ( firstBuf,  localHeightOfA*localWidth, sendCol, 0, 
              secondBuf, localHeight*localWidth,    recvCol, mpi::ANY_TAG,
              g.RowComm() );

            // Update with our received data
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(AVOID_OMP_FMA)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* secondBufCol = &secondBuf[jLoc*localHeight];
                T* thisCol = &thisBuf[jLoc*thisLDim];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    thisCol[iLoc] += alpha*secondBufCol[iLoc];
            }
            this->auxMemory_.Release();
        }
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SumScatterUpdate([* ,MC])");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

#ifdef VECTOR_WARNINGS
    if( A.Height() == 1 && g.Rank() == 0 )
    {
        std::cerr <<    
          "The vector version of [MR,MC].SumScatterUpdate([* ,MC]) is not "
          "yet written, but it only requires a modification of the vector "
          "version of [MR,MC].SumScatterUpdate([MR,* ])." << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.Rank() == 0 )
    {
        std::cerr <<
          "[MR,MC]::SumScatterUpdate([* ,MC]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MC] matrix instead." << std::endl;
    }
#endif
    if( this->RowAlignment() == A.RowAlignment() )
    {
        const Int c = g.Width();
        const Int colAlignment = this->ColAlignment();

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalHeight = MaxLength(height,c);
        const Int recvSize = mpi::Pad( maxLocalHeight*localWidth );
        const Int sendSize = c * recvSize;

        T* buffer = this->auxMemory_.Require( sendSize );

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisColShift = Shift_( k, colAlignment, c );
            const Int thisLocalHeight = Length_( height, thisColShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.RowComm() );

        // Update with our received data
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(AVOID_OMP_FMA)
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* bufferCol = &buffer[jLoc*localHeight];
            T* thisCol = &thisBuf[jLoc*thisLDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                thisCol[iLoc] += alpha*bufferCol[iLoc];
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned SumScatterUpdate [MR,MC] <- [* ,MC]" 
                      << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.Row();

        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalHeight = MaxLength(height,c);
        const Int recvSize_RS = mpi::Pad( maxLocalHeight*localWidthOfA );
        const Int sendSize_RS = c * recvSize_RS;
        const Int recvSize_SR = localHeight * localWidth;

        T* buffer = this->auxMemory_.Require
        ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[recvSize_RS];

        // Pack 
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &secondBuf[k*recvSize_RS];
            const Int thisColShift = Shift_( k, colAlignment, c );
            const Int thisLocalHeight = Length_( height, thisColShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // Reduce-scatter over each process row
        mpi::ReduceScatter
        ( secondBuf, firstBuf, recvSize_RS, mpi::SUM, g.RowComm() );

        // Trade reduced data with the appropriate process row
        mpi::SendRecv
        ( firstBuf,  localHeight*localWidthOfA, sendRow, 0,
          secondBuf, localHeight*localWidth,    recvRow, mpi::ANY_TAG,
          g.RowComm() );

        // Update with our received data
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(AVOID_OMP_FMA)
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* secondBufCol = &secondBuf[jLoc*localHeight];
            T* thisCol = &thisBuf[jLoc*thisLDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                thisCol[iLoc] += alpha*secondBufCol[iLoc];
        }
        this->auxMemory_.Release();
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SumScatterUpdate([* ,* ])");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

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
    const Int maxLocalHeight = MaxLength(height,c);
    const Int maxLocalWidth = MaxLength(width,r);
    const Int recvSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
    const Int sendSize = r * c * recvSize;

    T* buffer = this->auxMemory_.Require( sendSize );

    // Pack 
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int l=0; l<r; ++l )
    {
        const Int thisRowShift = Shift_( l, rowAlignment, r );
        const Int thisLocalWidth = Length_( width, thisRowShift, r );
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[(k+l*c)*recvSize];
            const Int thisColShift = Shift_( k, colAlignment, c );
            const Int thisLocalHeight = Length_( height, thisColShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = 
                    &ABuf[thisColShift+(thisRowShift+jLoc*r)*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }
    }

    // Communicate
    mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.VRComm() );

    // Unpack our received data
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* bufferCol = &buffer[jLoc*localHeight];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        blas::Axpy( localHeight, alpha, bufferCol, 1, thisCol, 1 );
    }
    this->auxMemory_.Release();
}

//
// Routines which explicitly work in the complex plane
//

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol*g.Height();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->SetLocalRealPart( iLoc, jLoc, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol*g.Height();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->SetLocalImagPart( iLoc, jLoc, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol*g.Height();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateLocalRealPart( iLoc, jLoc, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol*g.Height();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateLocalImagPart( iLoc, jLoc, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::GetRealPartOfDiagonal
( DistMatrix<BASE(T),MD,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::GetRealPartOfDiagonal([MD,* ])");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (length != d.Height() || d.Width() != 1) )
        throw std::logic_error("d is not of the correct dimensions");
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( this->DistData(), offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef BASE(T) R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( this->DistData(), offset );
        d.ResizeTo( length, 1 );
    }
    if( !d.Participating() )
        return;

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

    const Int iLocStart = (iStart-colShift) / c;
    const Int jLocStart = (jStart-rowShift) / r;

    const Int localDiagLength = d.LocalHeight();

    R* dBuf = d.Buffer();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/c);
        const Int jLoc = jLocStart + k*(lcm/r);
        dBuf[k] = RealPart( thisBuf[iLoc+jLoc*thisLDim] );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::GetImagPartOfDiagonal
( DistMatrix<BASE(T),MD,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::GetImagPartOfDiagonal([MD,* ])");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (length != d.Height() || d.Width() != 1) )
        throw std::logic_error("d is not of the correct dimensions");
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( this->DistData(), offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef BASE(T) R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( this->DistData(), offset );
        d.ResizeTo( length, 1 );
    }
    if( !d.Participating() )
        return;

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

    const Int iLocStart = (iStart-colShift) / c;
    const Int jLocStart = (jStart-rowShift) / r;

    const Int localDiagLength = d.LocalHeight();

    R* dBuf = d.Buffer();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/c);
        const Int jLoc = jLocStart + k*(lcm/r);
        dBuf[k] = ImagPart( thisBuf[iLoc+jLoc*thisLDim] );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::GetRealPartOfDiagonal
( DistMatrix<BASE(T),STAR,MD,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::GetRealPartOfDiagonal([* ,MD])");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (length != d.Width() || d.Height() != 1) )
        throw std::logic_error("d is not of the correct dimensions");
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( this->DistData(), offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef BASE(T) R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( this->DistData(), offset );
        d.ResizeTo( 1, length );
    }
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

    const Int iLocStart = (iStart-colShift) / c;
    const Int jLocStart = (jStart-rowShift) / r;

    const Int localDiagLength = d.LocalWidth();

    R* dBuf = d.Buffer();
    const Int dLDim = d.LDim();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/c);
        const Int jLoc = jLocStart + k*(lcm/r);
        dBuf[k*dLDim] = RealPart( thisBuf[iLoc+jLoc*thisLDim] );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::GetImagPartOfDiagonal
( DistMatrix<BASE(T),STAR,MD,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::GetImagPartOfDiagonal([* ,MD])");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (length != d.Width() || d.Height() != 1) )
        throw std::logic_error("d is not of the correct dimensions");
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( this->DistData(), offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef BASE(T) R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( this->DistData(), offset );
        d.ResizeTo( 1, length );
    }
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

    const Int iLocStart = (iStart-colShift) / c;
    const Int jLocStart = (jStart-rowShift) / r;

    const Int localDiagLength = d.LocalWidth();

    R* dBuf = d.Buffer();
    const Int dLDim = d.LDim();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/c);
        const Int jLoc = jLocStart + k*(lcm/r);
        dBuf[k*dLDim] = ImagPart( thisBuf[iLoc+jLoc*thisLDim] );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SetRealPartOfDiagonal
( const DistMatrix<BASE(T),MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SetRealPartOfDiagonal([MD,* ])");
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    if( !d.AlignedWithDiagonal( this->DistData(), offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
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
    const Int iLocStart = (iStart-colShift) / c;
    const Int jLocStart = (jStart-rowShift) / r;

    const Int localDiagLength = d.LocalHeight();
    const R* dBuf = d.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/c);
        const Int jLoc = jLocStart + k*(lcm/r);
        this->SetLocalRealPart( iLoc, jLoc, dBuf[k] );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SetImagPartOfDiagonal
( const DistMatrix<BASE(T),MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SetImagPartOfDiagonal([MD,* ])");
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    if( !d.AlignedWithDiagonal( this->DistData(), offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    typedef BASE(T) R;
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");
    if( !d.Participating() )
        return;

    const elem::Grid& g = this->Grid();
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
    const Int iLocStart = (iStart-colShift) / c;
    const Int jLocStart = (jStart-rowShift) / r;

    const Int localDiagLength = d.LocalHeight();
    const R* dBuf = d.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/c);
        const Int jLoc = jLocStart + k*(lcm/r);
        this->SetLocalImagPart( iLoc, jLoc, dBuf[k] );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SetRealPartOfDiagonal
( const DistMatrix<BASE(T),STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SetRealPartOfDiagonal([* ,MD])");
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    if( !d.AlignedWithDiagonal( this->DistData(), offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
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

    const Int iLocStart = (iStart-colShift) / c;
    const Int jLocStart = (jStart-rowShift) / r;

    const Int localDiagLength = d.LocalWidth();
    const R* dBuf = d.LockedBuffer();
    const Int dLDim = d.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/c);
        const Int jLoc = jLocStart + k*(lcm/r);
        this->SetLocalRealPart( iLoc, jLoc, dBuf[k*dLDim] );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,MC,Int>::SetImagPartOfDiagonal
( const DistMatrix<BASE(T),STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,MC]::SetImagPartOfDiagonal([* ,MD])");
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    if( !d.AlignedWithDiagonal( this->DistData(), offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
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

    const Int iLocStart = (iStart-colShift) / c;
    const Int jLocStart = (jStart-rowShift) / r;

    const Int localDiagLength = d.LocalWidth();
    const R* dBuf = d.LockedBuffer();
    const Int dLDim = d.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/c);
        const Int jLoc = jLocStart + k*(lcm/r);
        this->SetLocalImagPart( iLoc, jLoc, dBuf[k*dLDim] );
    }
}

template class DistMatrix<int,MR,MC,int>;
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,CIRC,CIRC,int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,MC,  MR,  int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,MC,  STAR,int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,MD,  STAR,int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,MR,  STAR,int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,STAR,MC,  int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,STAR,MD,  int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,STAR,MR,  int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,STAR,STAR,int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,STAR,VC,  int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,STAR,VR,  int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,VC,  STAR,int>& A );
template DistMatrix<int,MR,MC,int>::DistMatrix( const DistMatrix<int,VR,  STAR,int>& A );

#ifndef DISABLE_FLOAT
template class DistMatrix<float,MR,MC,int>;
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,CIRC,CIRC,int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,MC,  MR,  int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,MC,  STAR,int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,MD,  STAR,int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,MR,  STAR,int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,STAR,MC,  int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,STAR,MD,  int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,STAR,MR,  int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,STAR,STAR,int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,STAR,VC,  int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,STAR,VR,  int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,VC,  STAR,int>& A );
template DistMatrix<float,MR,MC,int>::DistMatrix( const DistMatrix<float,VR,  STAR,int>& A );
#endif // ifndef DISABLE_FLOAT

template class DistMatrix<double,MR,MC,int>;
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,CIRC,CIRC,int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,MC,  MR,  int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,MC,  STAR,int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,MD,  STAR,int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,MR,  STAR,int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,STAR,MC,  int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,STAR,MD,  int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,STAR,MR,  int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,STAR,STAR,int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,STAR,VC,  int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,STAR,VR,  int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,VC,  STAR,int>& A );
template DistMatrix<double,MR,MC,int>::DistMatrix( const DistMatrix<double,VR,  STAR,int>& A );

#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
template class DistMatrix<Complex<float>,MR,MC,int>;
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,CIRC,CIRC,int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,MC,  MR,  int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,MC,  STAR,int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,MD,  STAR,int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,MR,  STAR,int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,MC,  int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,MD,  int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,MR,  int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,STAR,int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,VC,  int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,VR,  int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,VC,  STAR,int>& A );
template DistMatrix<Complex<float>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<float>,VR,  STAR,int>& A );
#endif // ifndef DISABLE_FLOAT
template class DistMatrix<Complex<double>,MR,MC,int>;
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,CIRC,CIRC,int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,MC,  MR,  int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,MC,  STAR,int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,MD,  STAR,int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,MR,  STAR,int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,MC,  int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,MD,  int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,MR,  int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,STAR,int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,VC,  int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,VR,  int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,VC,  STAR,int>& A );
template DistMatrix<Complex<double>,MR,MC,int>::DistMatrix( const DistMatrix<Complex<double>,VR,  STAR,int>& A );
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
