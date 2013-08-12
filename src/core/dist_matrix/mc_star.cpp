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
DistMatrix<T,MC,STAR>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->SetShifts(); }

template<typename T>
DistMatrix<T,MC,STAR>::DistMatrix( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T>
DistMatrix<T,MC,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ 
    this->Align( colAlignment, 0 );
    this->ResizeTo( height, width );
}

template<typename T>
DistMatrix<T,MC,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ 
    this->Align( colAlignment, 0 );
    this->ResizeTo( height, width, ldim );
}

template<typename T>
DistMatrix<T,MC,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->LockedAttach( height, width, colAlignment, buffer, ldim, g ); }

template<typename T>
DistMatrix<T,MC,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Attach( height, width, colAlignment, buffer, ldim, g ); }

template<typename T>
DistMatrix<T,MC,STAR>::DistMatrix( const DistMatrix<T,MC,STAR>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MC,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this ) 
        *this = A;
    else
        LogicError("Tried to construct [MC,* ] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DistMatrix<T,MC,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MC,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( MC != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,MC,STAR>*>(&A) != this ) 
        *this = A;
    else
        LogicError("Tried to construct [MC,* ] with itself");
}

template<typename T>
DistMatrix<T,MC,STAR>::~DistMatrix()
{ }

template<typename T>
elem::DistData
DistMatrix<T,MC,STAR>::DistData() const
{
    elem::DistData data;
    data.colDist = MC;
    data.rowDist = STAR;
    data.colAlignment = this->colAlignment_;
    data.rowAlignment = 0;
    data.root = 0;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T>
Int
DistMatrix<T,MC,STAR>::ColStride() const
{ return this->grid_->Height(); }

template<typename T>
Int
DistMatrix<T,MC,STAR>::RowStride() const
{ return 1; }

template<typename T>
Int
DistMatrix<T,MC,STAR>::ColRank() const
{ return this->grid_->Row(); }

template<typename T>
Int
DistMatrix<T,MC,STAR>::RowRank() const
{ return 0; }

template<typename T>
void
DistMatrix<T,MC,STAR>::AlignWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,STAR]::AlignWith");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

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
DistMatrix<T,MC,STAR>::AlignWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,MC,STAR>::AlignColsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

template<typename T>
void
DistMatrix<T,MC,STAR>::AlignColsWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
bool
DistMatrix<T,MC,STAR>::AlignedWithDiagonal
( const elem::DistData& data, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::AlignedWithDiagonal");
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
            aligned = ( this->ColAlignment() == row );
        }
        else
        {
            const Int row = (alignment-offset) % this->ColStride();
            aligned = ( this->ColAlignment() == row );
        }
    }
    else aligned = false;
    return aligned;
}

template<typename T>
bool
DistMatrix<T,MC,STAR>::AlignedWithDiagonal
( const AbstractDistMatrix<T>& A, Int offset ) const
{ return this->AlignedWithDiagonal( A.DistData(), offset ); }

template<typename T>
void
DistMatrix<T,MC,STAR>::AlignWithDiagonal
( const elem::DistData& data, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::AlignWithDiagonal");
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
            this->colAlignment_ = row;
        }
        else 
        {
            const Int row = (alignment-offset) % this->ColStride();
            this->colAlignment_ = row;
        }
        this->constrainedColAlignment_ = true;
        this->SetShifts();
    }
#ifndef RELEASE
    else LogicError("Invalid diagonal alignment");
#endif
}

template<typename T>
void
DistMatrix<T,MC,STAR>::AlignWithDiagonal
( const AbstractDistMatrix<T>& A, Int offset )
{ this->AlignWithDiagonal( A.DistData(), offset ); }

template<typename T>
void
DistMatrix<T,MC,STAR>::Attach
( Int height, Int width, Int colAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::Attach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->viewType_ = VIEW;
    this->SetColShift();
    if( this->Participating() )
    {
        const Int localHeight = Length(height,this->colShift_,g.Height());
        this->matrix_.Attach_( localHeight, width, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::LockedAttach
( Int height, Int width, Int colAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->viewType_ = LOCKED_VIEW;
    this->SetColShift();
    if( this->Participating() )
    {
        const Int localHeight = Length(height,this->colShift_,g.Height());
        this->matrix_.LockedAttach_( localHeight, width, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),this->Grid().Height()), width );
}

template<typename T>
void
DistMatrix<T,MC,STAR>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),this->Grid().Height()), width, ldim );
}

template<typename T>
T DistMatrix<T,MC,STAR>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::Get");
    this->AssertValidEntry( i, j );
    if( !this->Participating() )
        LogicError("Should only be called by grid members");
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that
    // row within each process column
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();

    T u;
    if( g.Row() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        u = this->GetLocal( iLoc, j );
    }
    mpi::Broadcast( u, ownerRow, g.ColComm() );
    return u;
}

template<typename T>
void DistMatrix<T,MC,STAR>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        this->SetLocal(iLoc,j,u);
    }
}

template<typename T>
void DistMatrix<T,MC,STAR>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        this->UpdateLocal(iLoc,j,u);
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::GetDiagonal
( DistMatrix<T,MC,STAR>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
    if( (d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        LogicError("d must be aligned with the 'offset' diagonal");
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
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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
    const Int localDiagLength = d.LocalHeight();
    T* dBuf = d.Buffer();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        dBuf[k] = thisBuf[iLoc+jLoc*thisLDim];
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::GetDiagonal
( DistMatrix<T,STAR,MC>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetDiagonal");
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
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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
    const Int localDiagLength = d.LocalWidth();
    T* dBuf = d.Buffer();
    const Int dLDim = d.LDim();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        dBuf[k*dLDim] = thisBuf[iLoc+jLoc*thisLDim];
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::SetDiagonal
( const DistMatrix<T,MC,STAR>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetDiagonal");
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
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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

    const Int iLocStart = (iStart-colShift)/r;
    const Int localDiagLength = d.LocalHeight();
    const T* dBuf = d.LockedBuffer();
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        thisBuf[iLoc+jLoc*thisLDim] = dBuf[k];
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::SetDiagonal
( const DistMatrix<T,STAR,MC>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetDiagonal");
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
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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

    const Int iLocStart = (iStart-colShift)/r;
    const Int localDiagLength = d.LocalWidth();
    const T* dBuf = d.LockedBuffer();
    T* thisBuf = this->Buffer();
    const Int dLDim = d.LDim();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        thisBuf[iLoc+jLoc*thisLDim] = dBuf[k*dLDim];
    }
}

//
// Utility functions, e.g., SumOverRow
//

template<typename T>
void DistMatrix<T,MC,STAR>::SumOverRow()
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SumOverRow");
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
    mpi::AllReduce( sendBuf, recvBuf, localSize, this->Grid().RowComm() );

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

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetColAlignmentAndResize( A.ColAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlignment() == A.ColAlignment() )
    {
        if( A.Width() == 1 )
        {
            if( g.Col() == A.RowAlignment() )
                this->matrix_ = A.LockedMatrix();

            // Communicate
            mpi::Broadcast
            ( this->matrix_.Buffer(), this->LocalHeight(), A.RowAlignment(), 
              g.RowComm() );
        }
        else
        {
            const Int c = g.Width();

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidthOfA = A.LocalWidth();
            const Int maxLocalWidth = MaxLength(width,c);

            const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );
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
                T* sendBufCol = &sendBuf[jLoc*localHeight];
                MemCopy( sendBufCol, ACol, localHeight );
            }

            // Communicate
            mpi::AllGather
            ( sendBuf, portionSize, recvBuf, portionSize, g.RowComm() );

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
                {
                    const T* dataCol = &data[jLoc*localHeight];
                    T* thisCol = &thisBuf[(rowShift+jLoc*c)*thisLDim];
                    MemCopy( thisCol, dataCol, localHeight );
                }
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MC,* ] <- [MC,MR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.Row();

        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
        const Int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

        if( A.Width()==1 )
        {
            const Int localHeight = this->LocalHeight();

            if( this->grid_->Col() == A.RowAlignment() )
            {
                const Int localHeightOfA = A.LocalHeight();
                T* buffer = this->auxMemory_.Require( localHeightOfA );

                // Pack
                const T* ACol = A.LockedBuffer(0,0);
                MemCopy( buffer, ACol, localHeightOfA );

                // Communicate
                mpi::SendRecv
                ( buffer, localHeightOfA, sendRow,
                  this->matrix_.Buffer(), localHeight, recvRow, g.ColComm() );

                this->auxMemory_.Release();
            }

            // Communicate
            mpi::Broadcast
            ( this->matrix_.Buffer(), localHeight, A.RowAlignment(),
              g.RowComm() );
        }
        else
        {
            const Int height = this->Height();
            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localHeightOfA = A.LocalHeight();
            const Int localWidthOfA  = A.LocalWidth();
            const Int maxLocalHeight = MaxLength(height,r);
            const Int maxLocalWidth  = MaxLength(width,c);

            const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
            T* buffer = this->auxMemory_.Require( (c+1)*portionSize );
            T* firstBuffer = &buffer[0];
            T* secondBuf = &buffer[portionSize];

            // Pack the currently owned local data of A into the second 
            // buffer
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
            ( secondBuf,   portionSize, sendRow,
              firstBuffer, portionSize, recvRow, g.ColComm() );

            // Use the output of the SendRecv as the input to the AllGather
            mpi::AllGather
            ( firstBuffer, portionSize, secondBuf, portionSize, g.RowComm() );

            // Unpack the contents of each member of the process row
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
            const Int rowAlignmentOfA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                const T* data = &secondBuf[k*portionSize];
                const Int rowShift = Shift_( k, rowAlignmentOfA, c ); 
                const Int localWidth = Length_( width, rowShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                {
                    const T* dataCol = &data[jLoc*localHeight];
                    T* thisCol = &thisBuf[(rowShift+jLoc*c)*thisLDim];
                    MemCopy( thisCol, dataCol, localHeight );
                }
            }
            this->auxMemory_.Release();
        }
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetColAlignmentAndResize( A.ColAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlignment() == A.ColAlignment() )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MC,* ] <- [MC,* ]." << std::endl;
#endif
        const Int rank = g.Row();
        const Int r = g.Height();

        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int sendRank = (rank+r+colAlignment-colAlignmentOfA) % r;
        const Int recvRank = (rank+r+colAlignmentOfA-colAlignment) % r;

        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localHeightOfA = A.LocalHeight();

        const Int sendSize = localHeightOfA * width;
        const Int recvSize = localHeight * width;

        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ABuf[j*ALDim];
            T* sendBufCol = &sendBuf[j*localHeightOfA];
            MemCopy( sendBufCol, ACol, localHeightOfA );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank,
          recvBuf, recvSize, recvRank, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* recvBufCol = &recvBuf[j*localHeight];
            T* thisCol = &thisBuf[j*thisLDim];
            MemCopy( thisCol, recvBufCol, localHeight );
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,MR]");
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(true,false,this->ColAlignment(),0,g);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [MD,* ]");
#endif
    LogicError("[MC,* ] = [MD,* ] not yet implemented");
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,STAR,MD>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,MD]");
#endif
    LogicError("[MC,* ] = [* ,MD] not yet implemented");
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [MR,MC]");
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,VR,STAR> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(A) );
    std::auto_ptr<DistMatrix<T,VC,STAR> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [MR,* ]");
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,VR,STAR> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(A) );
    std::auto_ptr<DistMatrix<T,VC,STAR> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,MC]");
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MR,MC> > 
        A_MR_MC( new DistMatrix<T,MR,MC>(A) );
    std::auto_ptr<DistMatrix<T,VR,STAR> > 
        A_VR_STAR( new DistMatrix<T,VR,STAR>(g) );
    *A_VR_STAR = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    std::auto_ptr<DistMatrix<T,VC,STAR> > 
        A_VC_STAR( new DistMatrix<T,VC,STAR>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,VC,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Width() == 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "The vector version of [MC,* ] <- [VC,* ] is not yet written, but"
          " it only requires a modification of the vector version of "
          "[* ,MR] <- [* ,VR]" << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "[MC,* ] <- [VC,* ] potentially causes a large amount of cache-"
          "thrashing. If possible avoid it by performing the redistribution"
          " with a (conjugate)-transpose: \n"
          "  [* ,MC].TransposeFrom([VC,* ])" << std::endl;
    }
#endif
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

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLength(height,p);

        const Int portionSize = mpi::Pad( maxLocalHeightOfA*width );
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
        {
            const T* ACol = &ABuf[j*ALDim];
            T* sendBufCol = &sendBuf[j*localHeightOfA];
            MemCopy( sendBufCol, ACol, localHeightOfA );
        }

        // Communicate 
        mpi::AllGather
        ( sendBuf, portionSize, recvBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int colShift = this->ColShift();
        const Int colAlignmentOfA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[k*portionSize];    
            const Int colShiftOfA = Shift_( row+r*k, colAlignmentOfA, p );
            const Int colOffset = (colShiftOfA-colShift) / r;
            const Int localHeight = Length_( height, colShiftOfA, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &thisBuf[colOffset+j*thisLDim];
                const T* sourceCol = &data[j*localHeight];
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
            std::cerr << "Unaligned [MC,* ] <- [VC,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int row = g.Row();
        const Int rank = g.VCRank();

        // Perform the SendRecv to make A have the same colAlignment
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int colShift = this->ColShift();

        const Int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const Int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLength(height,p);

        const Int portionSize = mpi::Pad( maxLocalHeightOfA*width );
        T* buffer = this->auxMemory_.Require( (c+1)*portionSize );
        T* firstBuffer = &buffer[0];
        T* secondBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ABuf[j*ALDim];
            T* secondBufCol = &secondBuf[j*localHeightOfA];
            MemCopy( secondBufCol, ACol, localHeightOfA );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuf,   portionSize, sendRank,
          firstBuffer, portionSize, recvRank, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer, portionSize, secondBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShiftOfA = Shift_( row+r*k, colAlignment, p );
            const Int colOffset = (colShiftOfA-colShift) / r;
            const Int localHeight = Length_( height, colShiftOfA, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &thisBuf[colOffset+j*thisLDim];
                const T* sourceCol = &data[j*localHeight];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc*c] = sourceCol[iLoc];
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,STAR,VC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,VC]");
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,STAR,VR> > 
        A_STAR_VR( new DistMatrix<T,STAR,VR>(A) );
    std::auto_ptr<DistMatrix<T,MC,MR> > 
        A_MC_MR
        ( new DistMatrix<T,MC,MR>(true,false,this->ColAlignment(),0,g) );
    *A_MC_MR = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater
    *this = *A_MC_MR;
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,VR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [VR,* ]");
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR> A_VC_STAR(true,this->ColAlignment(),g);
    A_VC_STAR = A;
    *this = A_VC_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,STAR,VR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,VR]");
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(true,false,this->ColAlignment(),0,g);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    this->ResizeTo( A.Height(), A.Width() );

    const Int r = this->Grid().Height(); 
    const Int colShift = this->ColShift();

    const Int localHeight = this->LocalHeight();
    const Int width = this->Width();

    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
    {
        T* destCol = &thisBuf[j*thisLDim];
        const T* sourceCol = &ABuf[colShift+j*ALDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            destCol[iLoc] = sourceCol[iLoc*r];
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MC,STAR>&
DistMatrix<T,MC,STAR>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [o ,o ]");
#endif
    DistMatrix<T,MC,MR> A_MC_MR( A.Grid() );
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
DistMatrix<T,MC,STAR>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetRealPart");
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        this->SetLocalRealPart( iLoc, j, u );
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetImagPart");
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        this->SetLocalImagPart( iLoc, j, u );
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::UpdateRealPart");
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        this->UpdateLocalRealPart( iLoc, j, u );
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::UpdateImagPart");
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        this->UpdateLocalImagPart( iLoc, j, u );
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::GetRealPartOfDiagonal
( DistMatrix<BASE(T),MC,STAR>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetRealPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( this->DistData(), offset ) )
        LogicError("d must be aligned with the offset diag");
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
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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
    const Int localDiagLength = d.LocalHeight();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    R* dBuf = d.Buffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        dBuf[k] = RealPart(thisBuf[iLoc+jLoc*thisLDim]);
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::GetImagPartOfDiagonal
( DistMatrix<BASE(T),MC,STAR>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetImagPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( this->DistData(), offset ) )
        LogicError("d must be aligned with the offset diag");
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
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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
    const Int localDiagLength = d.LocalHeight();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    R* dBuf = d.Buffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        dBuf[k] = ImagPart(thisBuf[iLoc+jLoc*thisLDim]);
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::GetRealPartOfDiagonal
( DistMatrix<BASE(T),STAR,MC>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetRealPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( this->DistData(), offset ) )
        LogicError("d must be aligned with the offset diag");
#endif
    typedef BASE(T) R;
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( this->DistData(), offset );
    }
    const Int length = this->DiagonalLength(offset);
    d.ResizeTo( 1, length );
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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
    const Int localDiagLength = d.LocalWidth();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    R* dBuf = d.Buffer();
    const Int dLDim = d.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        dBuf[k*dLDim] = RealPart(thisBuf[iLoc+jLoc*thisLDim]);
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::GetImagPartOfDiagonal
( DistMatrix<BASE(T),STAR,MC>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetImagPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( this->DistData(), offset ) )
        LogicError("d must be aligned with the offset diag");
#endif
    typedef BASE(T) R;
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( this->DistData(), offset );
    }
    const Int length = this->DiagonalLength(offset);
    d.ResizeTo( 1, length );
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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
    const Int localDiagLength = d.LocalWidth();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    R* dBuf = d.Buffer();
    const Int dLDim = d.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        dBuf[k*dLDim] = ImagPart(thisBuf[iLoc+jLoc*thisLDim]);
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::SetRealPartOfDiagonal
( const DistMatrix<BASE(T),MC,STAR>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetRealPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        LogicError( msg.str() );
    }
    if( !d.AlignedWithDiagonal( this->DistData(), offset ) )
        LogicError("d must be aligned with the 'offset' diagonal");
#endif
    typedef BASE(T) R;
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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

    const Int iLocStart = (iStart-colShift)/r;
    const Int localDiagLength = d.LocalHeight();
    const R* dBuf = d.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        this->SetLocalRealPart( iLoc, jLoc, dBuf[k] );
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::SetImagPartOfDiagonal
( const DistMatrix<BASE(T),MC,STAR>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetImagPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        LogicError( msg.str() );
    }
    if( !d.AlignedWithDiagonal( this->DistData(), offset ) )
        LogicError("d must be aligned with the 'offset' diagonal");
#endif
    this->ComplainIfReal();
    typedef BASE(T) R;
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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

    const Int iLocStart = (iStart-colShift)/r;
    const Int localDiagLength = d.LocalHeight();
    const R* dBuf = d.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        this->SetLocalImagPart( iLoc, jLoc, dBuf[k] );
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::SetRealPartOfDiagonal
( const DistMatrix<BASE(T),STAR,MC>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetRealPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Height() != 1 )
        LogicError("d must be a row vector");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        LogicError( msg.str() );
    }
    if( !d.AlignedWithDiagonal( this->DistData(), offset ) )
        LogicError("d must be aligned with the 'offset' diagonal");
#endif
    typedef BASE(T) R;
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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

    const Int iLocStart = (iStart-colShift)/r;
    const Int localDiagLength = d.LocalWidth();

    const R* dBuf = d.LockedBuffer();
    const Int dLDim = d.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        this->SetLocalRealPart( iLoc, jLoc, dBuf[k*dLDim] );
    }
}

template<typename T>
void
DistMatrix<T,MC,STAR>::SetImagPartOfDiagonal
( const DistMatrix<BASE(T),STAR,MC>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetImagPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Height() != 1 )
        LogicError("d must be a row vector");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        LogicError( msg.str() );
    }
    if( !d.AlignedWithDiagonal( this->DistData(), offset ) )
        LogicError("d must be aligned with the 'offset' diagonal");
#endif
    this->ComplainIfReal();
    typedef BASE(T) R;
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
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

    const Int iLocStart = (iStart-colShift)/r;
    const Int localDiagLength = d.LocalWidth();

    const R* dBuf = d.LockedBuffer();
    const Int dLDim = d.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        this->SetLocalImagPart( iLoc, jLoc, dBuf[k*dLDim] );
    }
}

#define PROTO(T) template class DistMatrix<T,MC,STAR>
#define COPY(T,CD,RD) \
  template DistMatrix<T,MC,STAR>::DistMatrix( const DistMatrix<T,CD,RD>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  MR); \
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
