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
using DM = DistMatrix<T,MC,STAR>;

template<typename T>
DM<T>::DistMatrix( const elem::Grid& g )
: ADM<T>(g)
{ this->SetShifts(); }

template<typename T>
DM<T>::DistMatrix( Int height, Int width, const elem::Grid& g )
: ADM<T>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T>
DM<T>::DistMatrix( Int height, Int width, Int colAlign, const elem::Grid& g )
: ADM<T>(g)
{ 
    this->Align( colAlign, 0 );
    this->ResizeTo( height, width );
}

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, Int ldim, const elem::Grid& g )
: ADM<T>(g)
{ 
    this->Align( colAlign, 0 );
    this->ResizeTo( height, width, ldim );
}

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, const T* buffer, Int ldim,
  const elem::Grid& g )
: ADM<T>(g)
{ this->LockedAttach( height, width, colAlign, buffer, ldim, g ); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, T* buffer, Int ldim,
  const elem::Grid& g )
: ADM<T>(g)
{ this->Attach( height, width, colAlign, buffer, ldim, g ); }

template<typename T>
DM<T>::DistMatrix( const DM<T>& A )
: ADM<T>(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::DistMatrix"))
    this->SetShifts();
    if( &A != this ) 
        *this = A;
    else
        LogicError("Tried to construct [MC,* ] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DM<T>::DistMatrix( const DistMatrix<T,U,V>& A )
: ADM<T>(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::DistMatrix"))
    this->SetShifts();
    if( MC != U || STAR != V || 
        reinterpret_cast<const DM<T>*>(&A) != this ) 
        *this = A;
    else
        LogicError("Tried to construct [MC,* ] with itself");
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
{ return this->grid_->ColComm(); }

template<typename T>
mpi::Comm
DM<T>::RedundantComm() const
{ return this->grid_->RowComm(); }

template<typename T>
mpi::Comm
DM<T>::CrossComm() const
{ return mpi::COMM_SELF; }

template<typename T>
mpi::Comm
DM<T>::ColComm() const
{ return this->grid_->ColComm(); }

template<typename T>
mpi::Comm
DM<T>::RowComm() const
{ return mpi::COMM_SELF; }

template<typename T>
Int
DM<T>::ColStride() const
{ return this->grid_->Height(); }
    
template<typename T>
Int 
DM<T>::RowStride() const
{ return 1; }

template<typename T>
void
DM<T>::AlignWith( const elem::DistData& data )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR]::AlignWith"))
    this->SetGrid( *data.grid );
    if( data.colDist == MC )
        this->AlignCols( data.colAlign );
    else if( data.rowDist == MC )
        this->AlignCols( data.rowAlign );
    else if( data.colDist == VC )
        this->AlignCols( data.colAlign % this->ColStride() );
    else if( data.rowDist == VC )
        this->AlignCols( data.rowAlign % this->ColStride() );
    DEBUG_ONLY(else LogicError("Nonsensical alignment"))
}

template<typename T>
void
DM<T>::AlignColsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

template<typename T>
bool
DM<T>::AlignedWithDiagonal( const elem::DistData& data, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::AlignedWithDiagonal"))
    const Grid& grid = this->Grid();
    if( grid != *data.grid )
        return false;

    bool aligned;
    if( (data.colDist == MC   && data.rowDist == STAR) ||
        (data.colDist == STAR && data.rowDist == MC  ) )
    {
        const Int alignment = ( data.colDist==MC ? data.colAlign 
                                                 : data.rowAlign );
        if( offset >= 0 )
        {
            const Int row = alignment;
            aligned = ( this->ColAlign() == row );
        }
        else
        {
            const Int row = (alignment-offset) % this->ColStride();
            aligned = ( this->ColAlign() == row );
        }
    }
    else aligned = false;
    return aligned;
}

template<typename T>
void
DM<T>::AlignWithDiagonal( const elem::DistData& data, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::AlignWithDiagonal"))
    this->SetGrid( *data.grid );
    if( (data.colDist == MC   && data.rowDist == STAR) ||
        (data.colDist == STAR && data.rowDist == MC  ) )
    {
        const Int alignment = ( data.colDist==MC ? data.colAlign
                                                 : data.rowAlign );
        if( offset >= 0 )
            this->AlignCols( alignment );
        else 
            this->AlignCols( (alignment-offset) % this->ColStride() );
    }
    DEBUG_ONLY(else LogicError("Invalid diagonal alignment"))
}

template<typename T>
void
DM<T>::Attach
( Int height, Int width, Int colAlign,
  T* buffer, Int ldim, const elem::Grid& g )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::Attach"))
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlign_ = colAlign;
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
DM<T>::LockedAttach
( Int height, Int width, Int colAlign,
  const T* buffer, Int ldim, const elem::Grid& g )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::LockedAttach"))
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlign_ = colAlign;
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
( DistMatrix<S,MC,STAR>& d, Int offset, Function func ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,* ]::GetDiagonalHelper");
        if( d.Viewing() )
            this->AssertSameGrid( d.Grid() );
        if( (d.Viewing() || d.ColConstrained() ) &&
            !d.AlignedWithDiagonal( *this, offset ) )
            LogicError("d must be aligned with the 'offset' diagonal");
    )
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ColConstrained() )
            d.AlignWithDiagonal( *this, offset );
    }
    const Int diagLength = this->DiagonalLength(offset);
    d.ResizeTo( diagLength, 1 );
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
    const Int diagShift = d.ColShift();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int iLocStart = (iStart-colShift) / r;
    const Int localDiagLength = d.LocalHeight();
    S* dBuf = d.Buffer();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        func( dBuf[k], thisBuf[iLoc+jLoc*thisLDim] );
    }
}

template<typename T>
template<typename S,class Function>
void
DM<T>::GetDiagonalHelper
( DistMatrix<S,STAR,MC>& d, Int offset, Function func ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,* ]::GetDiagonalHelper");
        if( d.Viewing() )
            this->AssertSameGrid( d.Grid() );
        if( ( d.Viewing() || d.RowConstrained() ) &&
            !d.AlignedWithDiagonal( *this, offset ) )
            LogicError("d must be aligned with the 'offset' diagonal");
    )
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.RowConstrained() )
            d.AlignWithDiagonal( *this, offset );
    }
    const Int diagLength = this->DiagonalLength(offset);
    d.ResizeTo( 1, diagLength );
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
    const Int diagShift = d.RowShift();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int iLocStart = (iStart-colShift) / r;
    const Int localDiagLength = d.LocalWidth();
    S* dBuf = d.Buffer();
    const Int dLDim = d.LDim();
    const T* thisBuf = this->LockedBuffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        func( dBuf[k*dLDim], thisBuf[iLoc+jLoc*thisLDim] );
    }
}

template<typename T>
void
DM<T>::GetDiagonal( DM<T>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::GetDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T>
void
DM<T>::GetDiagonal( DistMatrix<T,STAR,MC>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::GetDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T>
void
DM<T>::GetRealPartOfDiagonal( DistMatrix<Base<T>,MC,STAR>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::GetRealPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = RealPart(beta); } );
}

template<typename T>
void
DM<T>::GetRealPartOfDiagonal( DistMatrix<Base<T>,STAR,MC>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::GetRealPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = RealPart(beta); } );
}

template<typename T>
void
DM<T>::GetImagPartOfDiagonal( DistMatrix<Base<T>,MC,STAR>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::GetImagPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = ImagPart(beta); } );
}

template<typename T>
void
DM<T>::GetImagPartOfDiagonal( DistMatrix<Base<T>,STAR,MC>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::GetImagPartOfDiagonal"))
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
DistMatrix<Base<T>,MC,STAR>
DM<T>::GetRealPartOfDiagonal( Int offset ) const
{
    DistMatrix<Base<T>,MC,STAR> d( this->Grid() );
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
DistMatrix<Base<T>,MC,STAR>
DM<T>::GetImagPartOfDiagonal( Int offset ) const
{
    DistMatrix<Base<T>,MC,STAR> d( this->Grid() );
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
template<typename S,class Function>
void
DM<T>::SetDiagonalHelper
( const DistMatrix<S,MC,STAR>& d, Int offset, Function func )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,* ]::SetDiagonalHelper");
        this->AssertSameGrid( d.Grid() );
        if( d.Width() != 1 )
            LogicError("d must be a column vector");
        const Int diagLength = this->DiagonalLength(offset);
        if( diagLength != d.Height() )
            LogicError
            ("d is not of the same length as the diagonal:\n",
             DimsString(*this,"A"),"\n",
             DimsString(d,"d"),"\n",
             "  A diag length: ",diagLength);
        if( !d.AlignedWithDiagonal( *this, offset ) )
            LogicError("d must be aligned with the 'offset' diagonal");
    )
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
    const Int diagShift = d.ColShift();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int iLocStart = (iStart-colShift)/r;
    const Int localDiagLength = d.LocalHeight();
    const S* dBuf = d.LockedBuffer();
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        func( thisBuf[iLoc+jLoc*thisLDim], dBuf[k] );
    }
}

template<typename T>
template<typename S,class Function>
void
DM<T>::SetDiagonalHelper
( const DistMatrix<S,STAR,MC>& d, Int offset, Function func )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,* ]::SetDiagonalHelper");
        this->AssertSameGrid( d.Grid() );
        if( d.Height() != 1 )
            LogicError("d must be a row vector");
        const Int diagLength = this->DiagonalLength(offset);
        if( diagLength != d.Width() )
            LogicError
            ("d is not of the same length as the diagonal:\n",
             DimsString(*this,"A"),"\n",
             DimsString(d,"d"),"\n",
             "  A diag length: ",diagLength);
        if( !d.AlignedWithDiagonal( *this, offset ) )
            LogicError("d must be aligned with the 'offset' diagonal");
    )
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int r = g.Height();
    const Int colShift = this->ColShift();
    const Int diagShift = d.RowShift();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int iLocStart = (iStart-colShift)/r;
    const Int localDiagLength = d.LocalWidth();
    const S* dBuf = d.LockedBuffer();
    T* thisBuf = this->Buffer();
    const Int dLDim = d.LDim();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart+k;
        const Int jLoc = jStart+k*r;
        func( thisBuf[iLoc+jLoc*thisLDim], dBuf[k*dLDim] );
    }
}

template<typename T>
void
DM<T>::SetDiagonal( const DM<T>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::SetDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T>
void
DM<T>::SetDiagonal( const DistMatrix<T,STAR,MC>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::SetDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T>
void
DM<T>::SetRealPartOfDiagonal( const DistMatrix<Base<T>,MC,STAR>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::SetRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetRealPart(alpha,beta); } );
}

template<typename T>
void
DM<T>::SetRealPartOfDiagonal( const DistMatrix<Base<T>,STAR,MC>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::SetRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetRealPart(alpha,beta); } );
}

template<typename T>
void
DM<T>::SetImagPartOfDiagonal( const DistMatrix<Base<T>,MC,STAR>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::SetImagPartOfDiagonal"))
    this->ComplainIfReal();
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetImagPart(alpha,beta); } );
}

template<typename T>
void
DM<T>::SetImagPartOfDiagonal( const DistMatrix<Base<T>,STAR,MC>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::SetImagPartOfDiagonal"))
    this->ComplainIfReal();
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetImagPart(alpha,beta); } );
}

//
// Utility functions, e.g., SumOverRow
//

template<typename T>
void DM<T>::SumOverRow()
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,* ]::SumOverRow");
        this->AssertNotLocked();
    )
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
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* thisCol = &thisBuf[jLoc*thisLDim];
        T* sendBufCol = &sendBuf[jLoc*localHeight];
        MemCopy( sendBufCol, thisCol, localHeight );
    }

    // AllReduce sum
    mpi::AllReduce( sendBuf, recvBuf, localSize, this->Grid().RowComm() );

    // Unpack
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* recvBufCol = &recvBuf[jLoc*localHeight];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        MemCopy( thisCol, recvBufCol, localHeight );
    }
    this->auxMemory_.Release();
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,* ] = [MC,MR]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
    this->AlignColsAndResize( A.ColAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlign() == A.ColAlign() )
    {
        if( A.Width() == 1 )
        {
            if( g.Col() == A.RowAlign() )
                this->matrix_ = A.LockedMatrix();

            // Communicate
            mpi::Broadcast
            ( this->matrix_.Buffer(), this->LocalHeight(), A.RowAlign(), 
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
            PARALLEL_FOR
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
            const Int rowAlignOfA = A.RowAlign();
            OUTER_PARALLEL_FOR
            for( Int k=0; k<c; ++k )
            {
                const T* data = &recvBuf[k*portionSize];
                const Int rowShift = Shift_( k, rowAlignOfA, c );
                const Int localWidth = Length_( width, rowShift, c );
                INNER_PARALLEL_FOR
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

        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();
        const Int sendRow = (row+r+colAlign-colAlignOfA) % r;
        const Int recvRow = (row+r+colAlignOfA-colAlign) % r;

        if( A.Width()==1 )
        {
            const Int localHeight = this->LocalHeight();

            if( this->grid_->Col() == A.RowAlign() )
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
            ( this->matrix_.Buffer(), localHeight, A.RowAlign(),
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
            PARALLEL_FOR
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
            const Int rowAlignOfA = A.RowAlign();
            OUTER_PARALLEL_FOR
            for( Int k=0; k<c; ++k )
            {
                const T* data = &secondBuf[k*portionSize];
                const Int rowShift = Shift_( k, rowAlignOfA, c ); 
                const Int localWidth = Length_( width, rowShift, c );
                INNER_PARALLEL_FOR
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
const DM<T>&
DM<T>::operator=( const DM<T>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,* ] = [MC,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
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
            std::cerr << "Unaligned [MC,* ] <- [MC,* ]." << std::endl;
#endif
        const Int rank = g.Row();
        const Int r = g.Height();

        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();

        const Int sendRank = (rank+r+colAlign-colAlignOfA) % r;
        const Int recvRank = (rank+r+colAlignOfA-colAlign) % r;

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
          recvBuf, recvSize, recvRank, g.ColComm() );

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
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [* ,MR]"))
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(true,false,this->ColAlign(),0,g);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [MD,* ]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [* ,MD]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [MR,MC]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(A) );
    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(true,this->ColAlign(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [MR,* ]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(A) );
    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(true,this->ColAlign(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [* ,MC]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,MR,MC>> 
        A_MR_MC( new DistMatrix<T,MR,MC>(A) );
    std::unique_ptr<DistMatrix<T,VR,STAR>> 
        A_VR_STAR( new DistMatrix<T,VR,STAR>(g) );
    *A_VR_STAR = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    std::unique_ptr<DistMatrix<T,VC,STAR>> 
        A_VC_STAR( new DistMatrix<T,VC,STAR>(true,this->ColAlign(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VC,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,* ] = [VC,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
    this->AlignColsAndResize( A.ColAlign()%g.Height(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlign() == A.ColAlign() % g.Height() )
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
        PARALLEL_FOR
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
        const Int colAlignOfA = A.ColAlign();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[k*portionSize];    
            const Int colShiftOfA = Shift_( row+r*k, colAlignOfA, p );
            const Int colOffset = (colShiftOfA-colShift) / r;
            const Int localHeight = Length_( height, colShiftOfA, p );
            INNER_PARALLEL_FOR
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

        // Perform the SendRecv to make A have the same colAlign
        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();
        const Int colShift = this->ColShift();

        const Int sendRank = (rank+p+colAlign-colAlignOfA) % p;
        const Int recvRank = (rank+p+colAlignOfA-colAlign) % p;

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
        PARALLEL_FOR
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
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShiftOfA = Shift_( row+r*k, colAlign, p );
            const Int colOffset = (colShiftOfA-colShift) / r;
            const Int localHeight = Length_( height, colShiftOfA, p );
            INNER_PARALLEL_FOR
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [* ,VC]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,STAR,VR>> 
        A_STAR_VR( new DistMatrix<T,STAR,VR>(A) );
    std::unique_ptr<DistMatrix<T,MC,MR>> 
        A_MC_MR
        ( new DistMatrix<T,MC,MR>(true,false,this->ColAlign(),0,g) );
    *A_MC_MR = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater
    *this = *A_MC_MR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [VR,* ]"))
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR> A_VC_STAR(true,this->ColAlign(),g);
    A_VC_STAR = A;
    *this = A_VC_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [* ,VR]"))
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(true,false,this->ColAlign(),0,g);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,* ] = [* ,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    this->ResizeTo( A.Height(), A.Width() );

    const Int r = this->Grid().Height(); 
    const Int colShift = this->ColShift();

    const Int localHeight = this->LocalHeight();
    const Int width = this->Width();

    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [o ,o ]"))
    DistMatrix<T,MC,MR> A_MC_MR( A.Grid() );
    A_MC_MR.AlignWith( *this );
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

#define PROTO(T) template class DistMatrix<T,MC,STAR>
#define COPY(T,U,V) \
  template DistMatrix<T,MC,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
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
