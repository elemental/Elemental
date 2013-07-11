/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

/*
 * DistMatrix_Dist
 */

template<typename Int>
DistMatrix_Dist<MC,STAR,Int>::DistMatrix_Dist( const elem::Grid& g )
: DistMatrix_Base<Int>(g)
{ this->SetShifts(); }

template<typename Int>
DistMatrix_Dist<MC,STAR,Int>::DistMatrix_Dist( const elem::Grid& g, Int colAlignment )
: DistMatrix_Base<Int>(g)
{ this->AlignCols( colAlignment ); }

template <typename Int>
elem::Distribution
DistMatrix_Dist<MC,STAR,Int>::ColDist() const { return MC; }

template <typename Int>
elem::Distribution
DistMatrix_Dist<MC,STAR,Int>::RowDist() const { return STAR; }

template<typename Int>
Int
DistMatrix_Dist<MC,STAR,Int>::ColRank() const
{ return this->grid_->Row(); }

template<typename Int>
Int
DistMatrix_Dist<MC,STAR,Int>::RowRank() const
{ return 0; }

template<typename Int>
Int
DistMatrix_Dist<MC,STAR,Int>::ColStride() const
{ return this->grid_->Height(); }

template<typename Int>
Int
DistMatrix_Dist<MC,STAR,Int>::RowStride() const
{ return 1; }

template <typename Int>
void 
DistMatrix_Dist<MC,STAR,Int>::Attach
( Int height, Int width, Int colAlignment, void* buffer, Int ldim, const elem::Grid& g )
{ DistMatrix_Base<Int>::Attach( height, width, colAlignment, 0, buffer, ldim, g ); }

template <typename Int>
void 
DistMatrix_Dist<MC,STAR,Int>::LockedAttach
( Int height, Int width, Int colAlignment, const void* buffer, Int ldim, const elem::Grid& g )
{ DistMatrix_Base<Int>::LockedAttach( height, width, colAlignment, 0, buffer, ldim, g ); }

template<typename Int>
void
DistMatrix_Dist<MC,STAR,Int>::AlignWith( const DistMatrix_Base<Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,STAR]::AlignWith");
#endif
    this->SetGrid( A.Grid() );
    elem::Distribution CD = A.ColDist(), RD = A.RowDist();
    if( CD == MC )
        this->colAlignment_ = A.colAlignment_;
    else if( RD == MC )
        this->colAlignment_ = A.rowAlignment_;
    else if( CD == VC )
        this->colAlignment_ = A.colAlignment_ % this->ColStride();
    else if( RD == VC )
        this->colAlignment_ = A.rowAlignment_ % this->ColStride();
#ifndef RELEASE
    else throw std::logic_error("Nonsensical alignment");
#endif
    this->constrainedColAlignment_ = true;
    this->SetShifts();
}

template<typename Int>
void
DistMatrix_Dist<MC,STAR,Int>::AlignColsWith( const DistMatrix_Base<Int>& A )
{ this->AlignWith( A ); }

template<typename Int>
bool
DistMatrix_Dist<MC,STAR,Int>::AlignedWithDiagonal
( const DistMatrix_Base<Int>& A, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::AlignedWithDiagonal");
#endif
    if( this->Grid() != A.Grid() )
        return false;
    elem::Distribution CD = A.ColDist(), RD = A.RowDist();

    bool aligned;
    if( (CD == MC   && RD == STAR) ||
        (CD == STAR && RD == MC  ) )
    {
        const Int alignment = CD == MC ? A.colAlignment_ : A.rowAlignment_;
        if( offset >= 0 )
        {
            const Int row = alignment;
            aligned = ( this->colAlignment_ == row );
        }
        else
        {
            const Int row = (alignment-offset) % this->ColStride();
            aligned = ( this->colAlignment_ == row );
        }
    }
    else aligned = false;
    return aligned;
}

template<typename Int>
void
DistMatrix_Dist<MC,STAR,Int>::AlignWithDiagonal
( const DistMatrix_Base<Int>& A, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::AlignWithDiagonal");
#endif
    this->SetGrid( A.Grid() );
    elem::Distribution CD = A.ColDist(), RD = A.RowDist();

    if( (CD == MC   && RD == STAR) ||
        (CD == STAR && RD == MC  ) )
    {
        const Int alignment = CD == MC ? A.colAlignment_ : A.rowAlignment_;
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
    else throw std::logic_error("Invalid diagonal alignment");
#endif
}

template<typename Int>
bool
DistMatrix_Dist<MC,STAR,Int>::Index( Int i, Int j, Int& iLocal, Int& jLocal, int& mpiSrc, mpi::Comm& mpiDst ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::Get");
    this->AssertValidEntry( i, j );
    if( !this->Participating() )
        throw std::logic_error("Should only be called by grid members");
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that
    // row within each process column
    const elem::Grid& g = this->Grid();
    mpiSrc = (i+this->ColAlignment()) % g.Height();
    mpiDst = g.ColComm();
    if ( g.Row() != mpiSrc ) return false;
    iLocal = (i-this->ColShift()) / g.Height();
    jLocal = j;
    return true;
}

/*
 * DistMatrix
 */

template<typename T,typename Int>
DistMatrix<T,MC,STAR,Int>::DistMatrix( const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,STAR>(g), DistMatrix_Type<T,Int>(g)
{ }

template<typename T,typename Int>
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,STAR>(g), DistMatrix_Type<T,Int>(g)
{ this->ResizeTo( height, width ); }

template<typename T,typename Int>
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,STAR>(g,colAlignment), DistMatrix_Type<T,Int>(g)
{ this->ResizeTo( height, width ); }

template<typename T,typename Int>
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int ldim, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,STAR>(g,colAlignment), DistMatrix_Type<T,Int>(g)
{ this->ResizeTo( height, width, ldim ); }

template<typename T,typename Int>
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,STAR>(g,colAlignment), DistMatrix_Type<T,Int>(g)
{ this->LockedAttach( height, width, colAlignment, buffer, ldim, g ); }

template<typename T,typename Int>
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,STAR>(g,colAlignment), DistMatrix_Type<T,Int>(g)
{ this->Attach( height, width, colAlignment, buffer, ldim, g ); }

template<typename T,typename Int>
DistMatrix<T,MC,STAR,Int>::DistMatrix( const DistMatrix<T,MC,STAR,Int>& A )
: DistMatrix_Base<Int>(A.Grid()), DistMatrix_Dist<MC,STAR>(A.Grid()), DistMatrix_Type<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MC,* ]::DistMatrix");
#endif
    if( &A != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [MC,* ] with itself");
}

template<typename T,typename Int>
template<Distribution U,Distribution V>
DistMatrix<T,MC,STAR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: DistMatrix_Base<Int>(A.Grid()), DistMatrix_Dist<MC,STAR>(A.Grid()), DistMatrix_Type<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MC,* ]::DistMatrix");
#endif
    if( MC != U || STAR != V || reinterpret_cast<const DistMatrix_Base<Int>*>(&A) != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [MC,* ] with itself");
}

template<typename T,typename Int>
DistMatrix<T,MC,STAR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::GetDiagonal
( DistMatrix<T,MC,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Height() || d.Width() != 1) )
    {
        std::ostringstream msg;
        msg << "d is not a column vec of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( (d.Viewing() || d.ConstrainedColAlignment() ) &&
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

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::GetDiagonal
( DistMatrix<T,STAR,MC,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Width() || d.Height() != 1) )
    {
        std::ostringstream msg;
        msg << "d is not a row vec of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
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

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::SetDiagonal
( const DistMatrix<T,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetDiagonal");
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

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::SetDiagonal
( const DistMatrix<T,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetDiagonal");
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

template<typename T,typename Int>
void DistMatrix<T,MC,STAR,Int>::SumOverRow()
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
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [MC,MR]");
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
                ( buffer, localHeightOfA, sendRow, 0,
                  this->matrix_.Buffer(), localHeight, recvRow, mpi::ANY_TAG,
                  g.ColComm() );

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
            ( secondBuf, portionSize, sendRow, 0,
              firstBuffer,  portionSize, recvRow, mpi::ANY_TAG, 
              g.ColComm() );

            // Use the output of the SendRecv as the input to the AllGather
            mpi::AllGather
            ( firstBuffer,  portionSize, 
              secondBuf, portionSize, g.RowComm() );

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

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [MC,* ]");
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
        this->ResizeTo( A.Height(), A.Width() );
    }
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
        ( sendBuf, sendSize, sendRank, 0,
          recvBuf, recvSize, recvRank, mpi::ANY_TAG, g.ColComm() );

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

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR,Int> A_MC_MR(true,false,this->ColAlignment(),0,g);

    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MC,* ] = [MD,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MC,* ] = [* ,MD] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = A;

    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = A;

    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MR,MC,Int> > 
        A_MR_MC( new DistMatrix<T,MR,MC,Int>(g) );
    *A_MR_MC = A;

    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > 
        A_VR_STAR( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > 
        A_VC_STAR( new DistMatrix<T,VC,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
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
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment() % g.Height();
            this->SetColShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
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
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.RowComm() );

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
        ( secondBuf, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, mpi::ANY_TAG, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
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

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,VC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > 
        A_STAR_VR( new DistMatrix<T,STAR,VR,Int>(g) );
    *A_STAR_VR = A;

    std::auto_ptr<DistMatrix<T,MC,MR,Int> > 
        A_MC_MR
        ( new DistMatrix<T,MC,MR,Int>(true,false,this->ColAlignment(),0,g) );
    *A_MC_MR = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_MC_MR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR,Int> A_VC_STAR(true,this->ColAlignment(),g);

    A_VC_STAR = A;
    *this = A_VC_STAR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR,Int> A_MC_MR(true,false,this->ColAlignment(),0,g);

    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
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

template<typename T,typename Int>
const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,CIRC,CIRC,Int>& A )
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

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::GetRealPartOfDiagonal
( DistMatrix<BASE(T),MC,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetRealPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength( offset );
#ifndef RELEASE
    if( d.Viewing() && (length != d.Height() || d.Width() != 1) )
    {
        std::ostringstream msg;
        msg << "d is not a column vec of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef BASE(T) R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( length, 1 );
    }
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

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::GetImagPartOfDiagonal
( DistMatrix<BASE(T),MC,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetImagPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength( offset );
#ifndef RELEASE
    if( d.Viewing() && (length != d.Height() || d.Width() != 1) )
    {
        std::ostringstream msg;
        msg << "d is not a column vec of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef BASE(T) R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( length, 1 );
    }
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

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::GetRealPartOfDiagonal
( DistMatrix<BASE(T),STAR,MC,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetRealPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef BASE(T) R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( 1, length );
    }
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

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::GetImagPartOfDiagonal
( DistMatrix<BASE(T),STAR,MC,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::GetImagPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef BASE(T) R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( 1, length );
    }
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

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::SetRealPartOfDiagonal
( const DistMatrix<BASE(T),MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetRealPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
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
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
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

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::SetImagPartOfDiagonal
( const DistMatrix<BASE(T),MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetImagPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
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
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
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

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::SetRealPartOfDiagonal
( const DistMatrix<BASE(T),STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetRealPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
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
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
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

template<typename T,typename Int>
void
DistMatrix<T,MC,STAR,Int>::SetImagPartOfDiagonal
( const DistMatrix<BASE(T),STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,* ]::SetImagPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
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
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
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

template class DistMatrix_Dist<MC,STAR,int>;

#define PROTO(T) \
  template class DistMatrix<T,MC,STAR,int>
#define COPY(T,CD,RD) \
  template DistMatrix<T,MC,STAR,int>::DistMatrix( \
    const DistMatrix<T,CD,RD,int>& A )
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
