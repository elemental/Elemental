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
using DM = DistMatrix<T,VC,STAR>;

template<typename T>
DM<T>::DistMatrix( const elem::Grid& g )
: ADM<T>(g)
{ this->SetShifts(); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: ADM<T>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, const elem::Grid& g )
: ADM<T>(g)
{ this->Align(colAlign,0); this->ResizeTo(height,width); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, Int ldim, const elem::Grid& g )
: ADM<T>(g)
{ this->Align(colAlign,0); this->ResizeTo(height,width,ldim); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, const T* buffer, Int ldim,
  const elem::Grid& g )
: ADM<T>(g)
{ this->LockedAttach(height,width,colAlign,buffer,ldim,g); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, T* buffer, Int ldim,
  const elem::Grid& g )
: ADM<T>(g)
{ this->Attach(height,width,colAlign,buffer,ldim,g); }

template<typename T>
DM<T>::DistMatrix( const DM<T>& A )
: ADM<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[VC,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [VC,* ] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DM<T>::DistMatrix( const DistMatrix<T,U,V>& A )
: ADM<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[VC,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( VC != U || STAR != V || 
        reinterpret_cast<const DM<T>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [VC,* ] with itself");
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
{ return this->grid_->VCComm(); }

template<typename T>
mpi::Comm
DM<T>::RowComm() const
{ return mpi::COMM_SELF; }

template<typename T>
Int
DM<T>::ColStride() const
{ return this->grid_->Size(); }

template<typename T>
Int
DM<T>::RowStride() const
{ return 1; }

template<typename T>
void
DM<T>::AlignWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::AlignWith");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );
    
    if( data.colDist == MC || data.colDist == VC )
        this->colAlign_ = data.colAlign;
    else if( data.rowDist == MC || data.rowDist == VC )
        this->colAlign_ = data.rowAlign;
#ifndef RELEASE
    else LogicError("Nonsensical alignment");
#endif
    this->colConstrained_ = true;
    this->SetShifts();
}

template<typename T>
void
DM<T>::AlignWith( const ADM<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DM<T>::AlignColsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

template<typename T>
void
DM<T>::AlignColsWith( const ADM<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
bool
DM<T>::AlignedWithDiagonal( const elem::DistData& data, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::AlignedWithDiagonal");
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
            aligned = ( this->ColAlign() == proc );
        }
        else
        {
            const Int proc = (alignment-offset) % this->ColStride();
            aligned = ( this->ColAlign() == proc );
        }
    }
    else aligned = false;
    return aligned;
}

template<typename T>
bool
DM<T>::AlignedWithDiagonal( const ADM<T>& A, Int offset ) const
{ return this->AlignedWithDiagonal( A.DistData(), offset ); }

template<typename T>
void
DM<T>::AlignWithDiagonal( const elem::DistData& data, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::AlignWithDiagonal");
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
            this->colAlign_ = proc;
        }
        else
        {
            const Int proc = (alignment-offset) % this->ColStride();
            this->colAlign_ = proc;
        }
        this->colConstrained_ = true;
        this->SetShifts();
    }
#ifndef RELEASE
    else LogicError("Invalid diagonal alignment");
#endif
}

template<typename T>
void
DM<T>::AlignWithDiagonal( const ADM<T>& A, Int offset )
{ this->AlignWithDiagonal( A.DistData(), offset ); }

template<typename T>
void
DM<T>::Attach
( Int height, Int width, Int colAlign,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::Attach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlign_ = colAlign;
    this->viewType_ = VIEW;
    this->SetColShift();
    if( this->Participating() )
    {
        const Int localHeight = Length(height,this->colShift_,g.Size());
        this->matrix_.Attach_( localHeight, width, buffer, ldim );
    }
}

template<typename T>
void
DM<T>::LockedAttach
( Int height, Int width, Int colAlign,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlign_ = colAlign;
    this->viewType_ = LOCKED_VIEW;
    this->SetColShift();
    if( this->Participating() )
    {
        const Int localHeight = Length(height,this->colShift_,g.Size());
        this->matrix_.LockedAttach_( localHeight, width, buffer, ldim );
    }
}

template<typename T>
void
DM<T>::Attach( Matrix<T>& A, Int colAlign, const elem::Grid& g )
{ this->Attach( A.Height(), A.Width(), colAlign, A.Buffer(), A.LDim(), g ); }

template<typename T>
void
DM<T>::LockedAttach( const Matrix<T>& A, Int colAlign, const elem::Grid& g )
{
    this->LockedAttach
    ( A.Height(), A.Width(), colAlign, A.LockedBuffer(), A.LDim(), g );
}

template<typename T>
template<typename S,class Function>
void
DM<T>::GetDiagonalHelper
( DistMatrix<S,VC,STAR>& d, Int offset, Function func ) const
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::GetDiagonalHelper");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
    if( (d.Viewing() || d.ColConstrained() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        LogicError("d must be aligned with the 'offset' diagonal");
#endif
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ColConstrained() )
            d.AlignWithDiagonal( *this, offset );
    }
    d.ResizeTo( this->DiagonalLength(offset), 1 );
    if( !this->Participating() )
        return;

    const Int diagShift = d.ColShift();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int p = g.Size();
    const Int iLocStart = (iStart-this->ColShift()) / p;

    S* dBuf = d.Buffer();
    const Int localDiagLength = d.LocalHeight();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*p;
        func( dBuf[k], thisBuf[iLoc+jLoc*thisLDim] );
    }
}

template<typename T>
template<typename S,class Function>
void
DM<T>::GetDiagonalHelper
( DistMatrix<S,STAR,VC>& d, Int offset, Function func ) const
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::GetDiagonalHelper");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
    if( ( d.Viewing() || d.RowConstrained() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        LogicError("d must be aligned with the 'offset' diagonal");
#endif
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.RowConstrained() )
            d.AlignWithDiagonal( *this, offset );
    }
    d.ResizeTo( 1, this->DiagonalLength(offset) );
    if( !this->Participating() )
        return;

    const Int diagShift = d.RowShift();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int p = g.Size();
    const Int iLocStart = (iStart-this->ColShift()) / p;
    const Int localDiagLength = d.LocalWidth();
    S* dBuf = d.Buffer();
    const Int dLDim = d.LDim();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*p;
        func( dBuf[k*dLDim], thisBuf[iLoc+jLoc*thisLDim] );
    }
}

template<typename T>
void
DM<T>::GetDiagonal( DM<T>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::GetDiagonal");
#endif
    this->GetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T>
void
DM<T>::GetDiagonal( DistMatrix<T,STAR,VC>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::GetDiagonal");
#endif
    this->GetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T>
void
DM<T>::GetRealPartOfDiagonal( DistMatrix<Base<T>,VC,STAR>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::GetRealPartOfDiagonal");
#endif
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = RealPart(beta); } );
}

template<typename T>
void
DM<T>::GetRealPartOfDiagonal( DistMatrix<Base<T>,STAR,VC>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::GetRealPartOfDiagonal");
#endif
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = RealPart(beta); } );
}

template<typename T>
void
DM<T>::GetImagPartOfDiagonal( DistMatrix<Base<T>,VC,STAR>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::GetImagPartOfDiagonal");
#endif
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = ImagPart(beta); } );
}

template<typename T>
void
DM<T>::GetImagPartOfDiagonal( DistMatrix<Base<T>,STAR,VC>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::GetImagPartOfDiagonal");
#endif
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = ImagPart(beta); } );
}

template<typename T>
DM<T>
DM<T>::GetDiagonal( Int offset ) const
{
    DM<T> d( this->Grid() );
    GetDiagonal( d, offset );
    return d;
}

template<typename T>
DistMatrix<Base<T>,VC,STAR>
DM<T>::GetRealPartOfDiagonal( Int offset ) const
{
    DistMatrix<Base<T>,VC,STAR> d( this->Grid() );
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
DistMatrix<Base<T>,VC,STAR>
DM<T>::GetImagPartOfDiagonal( Int offset ) const
{
    DistMatrix<Base<T>,VC,STAR> d( this->Grid() );
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
template<typename S,class Function>
void
DM<T>::SetDiagonalHelper
( const DistMatrix<S,VC,STAR>& d, Int offset, Function func )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SetDiagonalHelper");
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

    const Int diagShift = d.ColShift();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int p = g.Size();
    const Int iLocStart = (iStart-this->ColShift())/p;

    const S* dBuf = d.LockedBuffer();
    const Int localDiagLength = d.LocalHeight();
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*p;
        func( thisBuf[iLoc+jLoc*thisLDim], dBuf[k] );
    }
}

template<typename T>
template<typename S,class Function>
void
DM<T>::SetDiagonalHelper
( const DistMatrix<S,STAR,VC>& d, Int offset, Function func )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SetDiagonalHelper");
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

    const Int diagShift = d.RowShift();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int p = g.Size();
    const Int iLocStart = (iStart-this->ColShift())/p;

    const S* dBuf = d.LockedBuffer();
    const Int localDiagLength = d.LocalWidth();
    T* thisBuf = this->Buffer();
    const Int dLDim = d.LDim();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*p;
        func( thisBuf[iLoc+jLoc*thisLDim], dBuf[k*dLDim] );
    }
}

template<typename T>
void
DM<T>::SetDiagonal( const DM<T>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SetDiagonal");
#endif
    this->SetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T>
void
DM<T>::SetDiagonal( const DistMatrix<T,STAR,VC>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SetDiagonal");
#endif
    this->SetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T>
void
DM<T>::SetRealPartOfDiagonal( const DistMatrix<Base<T>,VC,STAR>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SetRealPartOfDiagonal");
#endif
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetRealPart(alpha,beta); } );
}

template<typename T>
void
DM<T>::SetRealPartOfDiagonal( const DistMatrix<Base<T>,STAR,VC>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SetRealPartOfDiagonal");
#endif
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetRealPart(alpha,beta); } );
}

template<typename T>
void
DM<T>::SetImagPartOfDiagonal( const DistMatrix<Base<T>,VC,STAR>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SetImagPartOfDiagonal");
#endif
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetImagPart(alpha,beta); } );
}

template<typename T>
void
DM<T>::SetImagPartOfDiagonal( const DistMatrix<Base<T>,STAR,VC>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SetImagPartOfDiagonal");
#endif
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetImagPart(alpha,beta); } );
}

//
// Utility functions, e.g., operator=
//

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->AlignColsAndResize( A.ColAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlign() % g.Height() == A.ColAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int colShiftOfA = A.ColShift();
        const Int colAlign = this->ColAlign();
        const Int rowAlignOfA = A.RowAlign();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLength(height,p);
        const Int maxWidth = MaxLength(width,c);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*c*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuf[k*portionSize];
            const Int thisRank = row+k*r;
            const Int thisColShift = Shift_(thisRank,colAlign,p); 
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
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
            const Int thisRowShift = Shift_(k,rowAlignOfA,c);
            const Int thisLocalWidth = Length_(width,thisRowShift,c);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuf[(thisRowShift+jLoc*c)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [VC,* ] <- [MC,MR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int colShiftOfA = A.ColShift();
        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();
        const Int rowAlignOfA = A.RowAlign();
        
        const Int sendRow = (row+r+(colAlign%r)-colAlignOfA) % r;
        const Int recvRow = (row+r+colAlignOfA-(colAlign%r)) % r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLength(height,p);
        const Int maxWidth = MaxLength(width,c);
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
            const Int thisColShift = Shift_(thisRank,colAlign,p);
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // AllToAll to gather all of the unaligned [VC,*] data into firstBuf
        mpi::AllToAll
        ( secondBuf, portionSize, firstBuf, portionSize, g.RowComm() );

        // SendRecv: properly align the [VC,*] via a trade in the column
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
            const Int thisRowShift = Shift_(k,rowAlignOfA,c);
            const Int thisLocalWidth = Length_(width,thisRowShift,c);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuf[(thisRowShift+jLoc*c)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->AlignColsAndResize( A.ColAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlign() % g.Height() == A.ColAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int colShift = this->ColShift();
        const Int colShiftOfA = A.ColShift();
        const Int colOffset = (colShift-colShiftOfA) / r;
        
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();

        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisBuf[j*thisLDim];
            const T* sourceCol = &ABuf[colOffset+j*ALDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*c];
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [VC,* ] <- [MC,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int col = g.Col();
        const Int colShiftOfA = A.ColShift();
        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();

        // We will SendRecv A[VC,*] within our process column to fix alignments.
        const Int sendRow = (row+r+(colAlign%r)-colAlignOfA) % r;
        const Int recvRow = (row+r+colAlignOfA-(colAlign%r)) % r;
        const Int sendRank = sendRow + r*col;

        const Int sendColShift = Shift( sendRank, colAlign, p );
        const Int sendColOffset = (sendColShift-colShiftOfA) / r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localHeightOfSend = Length(height,sendColShift,p);

        const Int sendSize = localHeightOfSend * width;
        const Int recvSize = localHeight * width;

        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &sendBuf[j*localHeightOfSend];
            const T* sourceCol = &ABuf[sendColOffset+j*ALDim];
            for( Int iLoc=0; iLoc<localHeightOfSend; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*c];
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRow,
          recvBuf, recvSize, recvRow, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [* ,MR]");
#endif
    DistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [MD,* ]");
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
    CallStackEntry cse("[VC,* ] = [* ,MD]");
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
    CallStackEntry cse("[VC,* ] = [MR,MC]");
#endif
    DistMatrix<T,VR,STAR> A_VR_STAR( A );
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [MR,* ]");
#endif
    DistMatrix<T,VR,STAR> A_VR_STAR( A );
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [* ,MC]");
#endif
    std::unique_ptr<DistMatrix<T,MR,MC>> A_MR_MC( new DistMatrix<T,MR,MC>(A) );
    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(*A_MR_MC) );
    delete A_MR_MC.release(); // lowers memory highwater
    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DM<T>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = this->Grid();
    this->AlignColsAndResize( A.ColAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlign() == A.ColAlign() )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [VC,* ] <- [VC,* ]." << std::endl;
#endif
        const Int rank = g.VCRank();
        const Int p = g.Size();

        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();

        const Int sendRank = (rank+p+colAlign-colAlignOfA) % p;
        const Int recvRank = (rank+p+colAlignOfA-colAlign) % p;

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
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ABuf[j*ALDim];
            T* sendBufCol = &sendBuf[j*localHeightOfA];
            MemCopy( sendBufCol, ACol, localHeightOfA );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank,
          recvBuf, recvSize, recvRank, g.VCComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [* ,VC]");
#endif
    std::unique_ptr<DistMatrix<T,MR,MC>> A_MR_MC( new DistMatrix<T,MR,MC>(A) );
    std::unique_ptr<DM<T>> A_VR_STAR
    ( new DM<T>(*A_MR_MC) );
    delete A_MR_MC.release(); // lowers memory highwater
    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;
    
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localHeightOfA = A.LocalHeight();

    const Int sendSize = localHeightOfA * width;
    const Int recvSize = localHeight * width;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int rankCM = g.VCRank();
    const Int rankRM = g.VRRank();

    const Int colShift = this->ColShift();
    const Int colShiftOfA = A.ColShift();

    // Compute which colmajor rank has the colShift equal to our colShiftOfA
    const Int sendRankCM = (rankCM+(p+colShiftOfA-colShift)) % p;

    // Compute which colmajor rank has the A colShift that we need
    const Int recvRankRM = (rankRM+(p+colShift-colShiftOfA)) % p;
    const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

    T* buffer = this->auxMemory_.Require( sendSize + recvSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[sendSize];

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        const T* ACol = &ABuf[j*ALDim];
        T* sendBufCol = &sendBuf[j*localHeightOfA];
        MemCopy( sendBufCol, ACol, localHeightOfA );
    }

    // Communicate
    mpi::SendRecv
    ( sendBuf, sendSize, sendRankCM,
      recvBuf, recvSize, recvRankCM, g.VCComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        const T* recvBufCol = &recvBuf[j*localHeight];
        T* thisCol = &thisBuf[j*thisLDim];
        MemCopy( thisCol, recvBufCol, localHeight );
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [* ,VR]");
#endif
    DistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[VC,* ] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int p = g.Size();
    const Int colShift = this->ColShift();

    const Int localHeight = this->LocalHeight();
    const Int width = this->Width();

    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        T* destCol = &thisBuf[j*thisLDim];
        const T* sourceCol = &ABuf[colShift+j*ALDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            destCol[iLoc] = sourceCol[iLoc*p];
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

    const Int colAlign = this->ColAlign();
    const Int mLocal = this->LocalHeight();
    const Int pkgSize = mpi::Pad(MaxLength(m,p)*n);
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
        for( Int s=0; s<p; ++s )
        {
            const Int sLocalHeight = Length( m, s, p );
            const Int q = (colAlign+s) % p;
            for( Int j=0; j<n; ++j )
            {
                for( Int iLoc=0; iLoc<sLocalHeight; ++iLoc )
                {
                    const Int i = s + iLoc*p;
                    sendBuf[q*pkgSize+iLoc+j*sLocalHeight] =
                        ABuf[i+j*ALDim];
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
        for( Int j=0; j<n; ++j )
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                buffer[iLoc+j*ldim] = recvBuf[iLoc+j*mLocal];
        this->auxMemory_.Release();
    }

    return *this;
}

template<typename T>
void
DM<T>::SumScatterFrom( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SumScatterFrom( [MC,* ] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && A.Grid().Rank() == 0 )
    {
        std::cerr <<
          "[VC,* ]::SumScatterFrom([MC,* ]) potentially causes a large amount "
          "of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MC,* ] matrix instead." << std::endl;
    }
#endif
    this->AlignColsAndResize( A.ColAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return;

    if( this->ColAlign() % g.Height() == A.ColAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int colAlign = this->ColAlign();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int maxLocalHeight = MaxLength( height, p );
        const Int recvSize = mpi::Pad( maxLocalHeight*width );
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
            const Int thisColShift = Shift_( thisRank, colAlign, p );
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = Length_( height, thisColShift, p );
            INNER_PARALLEL_FOR
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &data[j*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+j*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.RowComm() );

        // Unpack our received data
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* bufferCol = &buffer[j*localHeight];
            T* thisCol = &thisBuf[j*thisLDim];
            MemCopy( thisCol, bufferCol, localHeight );
        }
        this->auxMemory_.Release();
    }
    else
    {
        LogicError
        ("Unaligned [VC,* ]::ReduceScatterFrom( [MC,* ] ) not implemented");
    }
}

template<typename T>
void
DM<T>::SumScatterFrom( const DistMatrix<T,STAR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SumScatterFrom( [* ,* ] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return;

    const Int p = g.Size();
    const Int colAlign = this->ColAlign();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int maxLocalHeight = MaxLength( height, p );
    const Int recvSize = mpi::Pad( maxLocalHeight*width );
    const Int sendSize = p*recvSize;

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    T* buffer = this->auxMemory_.Require( sendSize );
    OUTER_PARALLEL_FOR
    for( Int k=0; k<p; ++k )
    {
        T* data = &buffer[k*recvSize];
        const Int thisColShift = Shift_( k, colAlign, p );
        const Int thisLocalHeight = Length_( height, thisColShift, p );
        INNER_PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &data[j*thisLocalHeight];
            const T* sourceCol = &ABuf[thisColShift+j*ALDim];
            for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*p];
        }
    }

    // Communicate
    mpi::ReduceScatter( buffer, recvSize, g.VCComm() );

    // Unpack our received data
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        const T* bufferCol = &buffer[j*localHeight];
        T* thisCol = &thisBuf[j*thisLDim];
        MemCopy( thisCol, bufferCol, localHeight );
    }
    this->auxMemory_.Release();
}

template<typename T>
void
DM<T>::SumScatterUpdate
( T alpha, const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SumScatterUpdate( [MC,* ] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && A.Grid().Rank() == 0 )
    {
        std::cerr <<
          "[VC,* ]::SumScatterUpdate([MC,* ]) potentially causes a large amount"
          " of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MC,* ] matrix instead." << std::endl;
    }
#endif
    if( this->ColAlign() % g.Height() == A.ColAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int colAlign = this->ColAlign();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int maxLocalHeight = MaxLength( height, p );
        const Int recvSize = mpi::Pad( maxLocalHeight*width );
        const Int sendSize = c*recvSize;

        T* buffer = this->auxMemory_.Require( sendSize );

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = myRow+k*r;
            const Int thisColShift = Shift_( thisRank, colAlign, p );
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = Length_( height, thisColShift, p );
            INNER_PARALLEL_FOR
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &data[j*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+j*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.RowComm() );

        // Unpack our received data
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* bufferCol = &buffer[j*localHeight];
            T* thisCol = &thisBuf[j*thisLDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                thisCol[iLoc] += alpha*bufferCol[iLoc];
        }
        this->auxMemory_.Release();
    }
    else
    {
        LogicError
        ("Unaligned [VC,* ]::ReduceScatterUpdate( [MC,* ] ) not implemented");
    }
}

template<typename T>
void
DM<T>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[VC,* ]::SumScatterUpdate( [* ,* ] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int p = g.Size();
    const Int colAlign = this->ColAlign();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int maxLocalHeight = MaxLength( height, p );
    const Int recvSize = mpi::Pad( maxLocalHeight*width );
    const Int sendSize = p*recvSize;

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    T* buffer = this->auxMemory_.Require( sendSize );
    OUTER_PARALLEL_FOR
    for( Int k=0; k<p; ++k )
    {
        T* data = &buffer[k*recvSize];
        const Int thisColShift = Shift_( k, colAlign, p );
        const Int thisLocalHeight = Length_( height, thisColShift, p );
        INNER_PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &data[j*thisLocalHeight];
            const T* sourceCol = &ABuf[thisColShift+j*ALDim];
            for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*p];
        }
    }

    // Communicate
    mpi::ReduceScatter( buffer, recvSize, g.VCComm() );

    // Unpack our received data
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        const T* bufferCol = &buffer[j*localHeight];
        T* thisCol = &thisBuf[j*thisLDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            thisCol[iLoc] += alpha*bufferCol[iLoc];
    }
    this->auxMemory_.Release();
}

#define PROTO(T) template class DistMatrix<T,VC,STAR>
#define COPY(T,U,V) \
  template DistMatrix<T,VC,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
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
  COPY(T,STAR,VR  ); \
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
