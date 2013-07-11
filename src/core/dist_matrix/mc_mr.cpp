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

/*
 * DistMatrix_Dist
 */

template <typename Int>
DistMatrix_Dist<MC,MR,Int>::DistMatrix_Dist( const elem::Grid& g )
: DistMatrix_Base<Int>(g)
{ this->SetShifts(); }

template <typename Int>
DistMatrix_Dist<MC,MR,Int>::DistMatrix_Dist( const elem::Grid& g, Int rowAlignment, Int colAlignment )
: DistMatrix_Base<Int>(g)
{ this->Align( rowAlignment, colAlignment ); }

template <typename Int>
elem::Distribution
DistMatrix_Dist<MC,MR,Int>::ColDist() const { return MC; }

template <typename Int>
elem::Distribution
DistMatrix_Dist<MC,MR,Int>::RowDist() const { return MR; }

template <typename Int>
Int
DistMatrix_Dist<MC,MR,Int>::ColRank() const
{ return this->grid_->Row(); }

template <typename Int>
Int
DistMatrix_Dist<MC,MR,Int>::RowRank() const
{ return this->grid_->Col(); }

template <typename Int>
Int
DistMatrix_Dist<MC,MR,Int>::ColStride() const
{ return this->grid_->Height(); }

template <typename Int>
Int
DistMatrix_Dist<MC,MR,Int>::RowStride() const
{ return this->grid_->Width(); }

template <typename Int>
void
DistMatrix_Dist<MC,MR,Int>::AlignWith( const DistMatrix_Base<Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AlignWith");
#endif
    this->SetGrid( A.Grid() );
    elem::Distribution CD = A.ColDist(), RD = A.RowDist();
    if( CD == MC && RD == MR )
    {
        this->colAlignment_ = A.colAlignment_;
        this->rowAlignment_ = A.rowAlignment_;
        this->constrainedColAlignment_ = true;
        this->constrainedRowAlignment_ = true;
    }
    else if( CD == MC && RD == STAR )
    {
        this->colAlignment_ = A.colAlignment_;
        this->constrainedColAlignment_ = true; 
    }
    else if( CD == MR && RD == MC )
    {
        this->colAlignment_ = A.rowAlignment_;
        this->rowAlignment_ = A.colAlignment_;
        this->constrainedColAlignment_ = true;
        this->constrainedRowAlignment_ = true;
    }
    else if( CD == MR && RD == STAR )
    {
        this->rowAlignment_ = A.colAlignment_;
        this->constrainedRowAlignment_ = true;
    }
    else if( CD == STAR && RD == MC )
    {
        this->colAlignment_ = A.rowAlignment_;
        this->constrainedColAlignment_ = true;
    }
    else if( CD == STAR && RD == MR )
    {
        this->rowAlignment_ = A.rowAlignment_;
        this->constrainedRowAlignment_ = true;
    }
    else if( CD == STAR && RD == VC )
    {
        this->colAlignment_ = A.rowAlignment_ % this->ColStride();
        this->constrainedColAlignment_ = true;
    }
    else if( CD == STAR && RD == VR )
    {
        this->rowAlignment_ = A.rowAlignment_ % this->RowStride();
        this->constrainedRowAlignment_ = true;
    }
    else if( CD == VC && RD == STAR )
    {
        this->colAlignment_ = A.colAlignment_ % this->ColStride();
        this->constrainedColAlignment_ = true;
    }
    else if( CD == VR && RD == STAR )
    {
        this->rowAlignment_ = A.colAlignment_ % this->RowStride();
        this->constrainedRowAlignment_ = true;
    }
#ifndef RELEASE
    else throw std::logic_error("Nonsensical alignment");
#endif
    this->SetShifts();
}

template <typename Int>
void
DistMatrix_Dist<MC,MR,Int>::AlignColsWith( const DistMatrix_Base<Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AlignColsWith");
    if( this->grid_ != A.grid_ )
        throw std::logic_error("Grids do not match");
#endif
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

template <typename Int>
void
DistMatrix_Dist<MC,MR,Int>::AlignRowsWith( const DistMatrix_Base<Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AlignRowsWith");
    if( this->grid_ != A.grid_ )
        throw std::logic_error("Grids do not match");
#endif
    elem::Distribution CD = A.ColDist(), RD = A.RowDist();
    if( CD == MR )
        this->rowAlignment_ = A.colAlignment_;
    else if( RD == MR )
        this->rowAlignment_ = A.rowAlignment_;
    else if( CD == VR )
        this->rowAlignment_ = A.colAlignment_ % this->RowStride();
    else if( RD == VR )
        this->rowAlignment_ = A.rowAlignment_ % this->RowStride();
#ifndef RELEASE
    else throw std::logic_error("Nonsensical alignment");
#endif
    this->constrainedRowAlignment_ = true;
    this->SetShifts();
}

template <typename Int>
bool
DistMatrix_Dist<MC,MR,Int>::Index( Int i, Int j, Int& iLocal, Int& jLocal, int& mpiSrc, mpi::Comm& mpiDst ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::Index");
    this->AssertValidEntry( i, j );
#endif
    const Int ownerRow = (i + this->ColAlignment()) % this->ColStride();
    const Int ownerCol = (j + this->RowAlignment()) % this->RowStride();
    const Int ownerRank = ownerRow + ownerCol*this->ColStride();
    const elem::Grid& g = this->Grid();
    mpiSrc = g.VCToViewingMap(ownerRank);
    mpiDst = g.ViewingComm();
    if ( g.VCRank() != ownerRank ) return false;
    iLocal = (i-this->ColShift()) / this->ColStride();
    jLocal = (j-this->RowShift()) / this->RowStride();
    return true;
}

/*
 * DistMatrix
 */

template<typename T,typename Int>
DistMatrix<T,MC,MR,Int>::DistMatrix( const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,MR,Int>(g), DistMatrix_Type<T,Int>(g)
{ }

template<typename T,typename Int>
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,MR,Int>(g), DistMatrix_Type<T,Int>(g)
{ this->ResizeTo( height, width );  }

template<typename T,typename Int>
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,MR,Int>(g,colAlignment,rowAlignment), DistMatrix_Type<T,Int>(g)
{ this->ResizeTo( height, width );  }

template<typename T,typename Int>
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment, Int ldim, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,MR,Int>(g,colAlignment,rowAlignment), DistMatrix_Type<T,Int>(g)
{ this->ResizeTo( height, width, ldim );  }

template<typename T,typename Int>
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,MR,Int>(g), DistMatrix_Type<T,Int>(g)
{ this->LockedAttach( height, width, colAlignment, rowAlignment, buffer, ldim, g );  }

template<typename T,typename Int>
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<MC,MR,Int>(g), DistMatrix_Type<T,Int>(g)
{ this->Attach( height, width, colAlignment, rowAlignment, buffer, ldim, g ); }

template<typename T,typename Int>
DistMatrix<T,MC,MR,Int>::DistMatrix( const DistMatrix<T,MC,MR,Int>& A )
: DistMatrix_Base<Int>(A.Grid()), DistMatrix_Dist<MC,MR>(A.Grid()), DistMatrix_Type<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MC,MR]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [MC,MR] with itself");
}

template<typename T,typename Int>
template<Distribution U,Distribution V>
DistMatrix<T,MC,MR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: DistMatrix_Base<Int>(A.Grid()), DistMatrix_Dist<MC,MR>(A.Grid()), DistMatrix_Type<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MC,MR]::DistMatrix");
#endif
    if( MC != U || MR != V || reinterpret_cast<const DistMatrix_Base<Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [MC,MR] with itself");
}

template<typename T,typename Int>
DistMatrix<T,MC,MR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::GetDiagonal
( DistMatrix<T,MD,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && diagLength != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
    {
        std::ostringstream os;
        os << "offset:         " << offset << "\n"
           << "colAlignment:   " << this->colAlignment_ << "\n"
           << "rowAlignment:   " << this->rowAlignment_ << "\n"
           << "d.diagPath:     " << d.diagPath_ << "\n"
           << "d.colAlignment: " << d.colAlignment_ << std::endl;
        std::cerr << os.str();
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
    }
#endif
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( diagLength, 1 );
    }
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/colStride);
        const Int jLoc = jLocStart + k*(lcm/rowStride);
        dBuf[k] = buffer[iLoc+jLoc*ldim];
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::GetDiagonal
( DistMatrix<T,STAR,MD,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && diagLength != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/colStride);
        const Int jLoc = jLocStart + k*(lcm/rowStride);
        dBuf[k*dLDim] = buffer[iLoc+jLoc*ldim];
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::SetDiagonal
( const DistMatrix<T,MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetDiagonal");
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/colStride);
        const Int jLoc = jLocStart + k*(lcm/rowStride);
        buffer[iLoc+jLoc*ldim] = dBuf[k];
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::SetDiagonal
( const DistMatrix<T,STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetDiagonal");
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
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

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::AdjointFrom( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AdjointFrom");
#endif
    this->TransposeFrom( A, true );
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::TransposeFrom
( const DistMatrix<T,STAR,MC,Int>& A, bool conjugate )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::TransposeFrom");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Width(), A.Height() );
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.RowAlignment();
            this->SetColShift();
        }
        this->ResizeTo( A.Width(), A.Height() );
    }
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
        if( conjugate )
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &buffer[jLoc*ldim];
                const T* sourceCol = &ABuffer[rowShift+jLoc*rowStride];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc] = Conj( sourceCol[iLoc*ALDim] );
            }
        }
        else
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &buffer[jLoc*ldim];
                const T* sourceCol = &ABuffer[rowShift+jLoc*rowStride];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*ALDim];
            }
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
        if( conjugate )
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &sendBuf[jLoc*localWidth];
                const T* sourceCol = &ABuffer[rowShift+jLoc*rowStride];
                for( Int iLoc=0; iLoc<localWidthA; ++iLoc )
                    destCol[iLoc] = Conj( sourceCol[iLoc*ALDim] );
            }
        }
        else
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &sendBuf[jLoc*localWidth];
                const T* sourceCol = &ABuffer[rowShift+jLoc*rowStride];
                for( Int iLoc=0; iLoc<localWidthA; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*ALDim];
            }
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank, 0,
          recvBuf, recvSize, recvRank, mpi::ANY_TAG, g.ColComm() );

        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &buffer[jLoc*ldim], &recvBuf[jLoc*localHeight], localHeight );
        this->auxMemory_.Release();
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::AdjointFrom
( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AdjointFrom");
#endif
    this->TransposeFrom( A, true );
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::TransposeFrom
( const DistMatrix<T,MR,STAR,Int>& A, bool conjugate )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::TransposeFrom");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Width(), A.Height() );
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.ColAlignment();
            this->SetRowShift();
        }
        this->ResizeTo( A.Width(), A.Height() );
    }
    if( this->rowAlignment_ != A.ColAlignment() )
        throw std::logic_error("Unaligned TransposeFrom");

    if( this->Participating() ) 
    { 
        const Int colStride = this->ColStride();
        const Int colShift = this->ColShift();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        if( conjugate )
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &buffer[jLoc*ldim];
                const T* sourceCol = &ABuffer[jLoc+colShift*ALDim];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc] = Conj( sourceCol[iLoc*colStride*ALDim] );
            }
        }
        else
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &buffer[jLoc*ldim];
                const T* sourceCol = &ABuffer[jLoc+colShift*ALDim];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*colStride*ALDim];
            }
        }
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::AdjointSumScatterFrom
( const DistMatrix<T,MR,STAR,Int>& AAdj_MR_STAR )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AdjointSumScatterFrom");
#endif
    this->TransposeSumScatterFrom( AAdj_MR_STAR, true );
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::TransposeSumScatterFrom
( const DistMatrix<T,MR,STAR,Int>& ATrans_MR_STAR, bool conjugate )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::TransposeSumScatterFrom");
    this->AssertNotLocked();
    this->AssertSameGrid( ATrans_MR_STAR.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( ATrans_MR_STAR.Width(), ATrans_MR_STAR.Height() );
#endif
    const Grid& g = ATrans_MR_STAR.Grid();
    DistMatrix<T,MR,MC,Int> ATrans( g );
    if( this->Viewing() )
        ATrans.AlignWith( *this );
    ATrans.SumScatterFrom( ATrans_MR_STAR );
    Transpose( ATrans, *this, conjugate );
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::AdjointSumScatterUpdate
( T alpha, const DistMatrix<T,MR,STAR,Int>& AAdj_MR_STAR )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::AdjointSumScatterUpdate");
#endif
    this->TransposeSumScatterUpdate( alpha, AAdj_MR_STAR, true );
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::TransposeSumScatterUpdate
( T alpha, const DistMatrix<T,MR,STAR,Int>& ATrans_MR_STAR, bool conjugate )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::TransposeSumScatterUpdate");
    this->AssertNotLocked();
    this->AssertSameGrid( ATrans_MR_STAR.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( ATrans_MR_STAR.Width(), ATrans_MR_STAR.Height() );
#endif
    const Grid& g = ATrans_MR_STAR.Grid();
    DistMatrix<T,MR,MC,Int> ATrans( g );
    ATrans.SumScatterFrom( ATrans_MR_STAR );
    DistMatrix<T,MC,MR,Int> A( g );
    if( this->Viewing() )
        A.AlignWith( *this );
    Transpose( ATrans, A, conjugate );
    Axpy( alpha, A, *this );
}

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [MC,MR]");
    this->AssertNotLocked();
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
    {
        if( this->Grid() == A.Grid() )
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
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !this->Participating() && !A.Participating() )
        return *this;

    if( this->Grid() == A.Grid() )
    {
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
                MemCopy
                ( &sendBuf[jLoc*localHeightA], 
                  &ABuffer[jLoc*ALDim], localHeightA );

            // Communicate
            mpi::SendRecv
            ( sendBuf, sendSize, sendRank, 0,
              recvBuf, recvSize, recvRank, mpi::ANY_TAG, g.VCComm() );

            // Unpack
            T* buffer = this->Buffer();
            const Int ldim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &buffer[jLoc*ldim], 
                  &recvBuf[jLoc*localHeight], localHeight );
            this->auxMemory_.Release();
        }
    }
    else // the grids don't match
    {
        if( !mpi::CongruentComms( A.Grid().ViewingComm(), 
                                  this->Grid().ViewingComm() ) )
            throw std::logic_error
            ("Redistributing between nonmatching grids currently requires"
             " the viewing communicators to match.");

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

        const Int maxSendSize = 
            (A.Height()/(colStrideA*localColStrideA)+1) * 
            (A.Width()/(rowStrideA*localRowStrideA)+1);

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
            recvRow = (((colRankA+colStrideA-colAlignA) % colStrideA) + colAlign) % colStride;
        for( Int colSend=0; colSend<numColSends; ++colSend )
        {
            Int recvCol = 0; // avoid compiler warnings...
            if( inAGrid )
                recvCol = (((rowRankA+rowStrideA-rowAlignA) % rowStrideA) + rowAlign) % rowStride;
            for( Int rowSend=0; rowSend<numRowSends; ++rowSend )
            {
                mpi::Request sendRequest;
                // Fire off this round of non-blocking sends
                if( inAGrid )
                {
                    // Pack the data
                    Int sendHeight = Length( A.LocalHeight(), colSend, numColSends );
                    Int sendWidth = Length( A.LocalWidth(), rowSend, numRowSends );
                    const T* ABuffer = A.LockedBuffer();
                    const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
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
                    const Int recvViewingRank = this->Grid().VCToViewingMap( recvVCRank );
                    mpi::ISend
                    ( sendBuf, sendHeight*sendWidth, recvViewingRank,
                      0, this->Grid().ViewingComm(), sendRequest );
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
                            const Int sendViewingRank = A.Grid().VCToViewingMap( sendVCRank );

                            mpi::Recv
                            ( recvBuf, sendHeight*sendWidth, sendViewingRank, 0, this->Grid().ViewingComm() );
                            
                            // Unpack the data
                            T* buffer = this->Buffer();
                            const Int ldim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
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
    return *this;
}

// PAUSED PASS HERE

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [MC,* ]");
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
        const Int c = g.Width();
        const Int rowShift = this->RowShift();

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();

        const Int ALDim = A.LDim();
        const Int thisLDim = this->LDim();
        T* thisBuffer = this->Buffer();
        const T* ABuffer = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &sendBuf[jLoc*localHeightA], 
              &ABuffer[(rowShift+jLoc*c)*ALDim], localHeightA );

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank, 0,
          recvBuf, recvSize, recvRank, mpi::ANY_TAG, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuffer[jLoc*thisLDim], 
              &recvBuf[jLoc*localHeight], localHeight );
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,MR]");
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
        const Int r = g.Height();
        const Int colShift = this->ColShift();

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();

        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
        {
            T* destCol = &sendBuf[jLoc*localHeight];
            const T* sourceCol = &ABuffer[colShift+jLoc*ALDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*r];
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendCol, 0,
          recvBuf, recvSize, recvCol, mpi::ANY_TAG, g.RowComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for  
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuffer[jLoc*thisLDim], 
              &recvBuf[jLoc*localHeight], localHeight );
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MC,MR] = [MD,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MC,MR] = [* ,MD] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [MR,MC]");
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
#ifdef HAVE_OPENMP
#pragma omp parallel for  
#endif
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
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerRow, g.ColComm() );

        // A[VC,* ] <- A[VR,* ]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankCM, 0,
          recvBuf, portionSize, recvRankCM, mpi::ANY_TAG, g.VCComm() );

        // A[MC,MR] <- A[VC,* ]
        mpi::Gather
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerCol, g.RowComm() );

        if( myCol == ownerCol )
        {
            // Unpack
            T* thisBuffer = this->Buffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for  
#endif
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
#ifdef HAVE_OPENMP
#pragma omp parallel for  
#endif
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
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerCol, g.RowComm() );

        // A[* ,VR] <- A[* ,VC]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankRM, 0,
          recvBuf, portionSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

        // A[MC,MR] <- A[* ,VR]
        mpi::Gather
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerRow, g.ColComm() );
    
        if( myRow == ownerRow )
        {
            // Unpack
            T* thisBuffer = this->Buffer();
            const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for  
#endif
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
            std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
            ( new DistMatrix<T,VR,STAR,Int>(g) );

            *A_VR_STAR = A;

            std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
            ( new DistMatrix<T,VC,STAR,Int>(true,this->ColAlignment(),g) );
            *A_VC_STAR = *A_VR_STAR;
            delete A_VR_STAR.release(); // lowers memory highwater

            *this = *A_VC_STAR;
        }
        else
        {
            std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
            ( new DistMatrix<T,STAR,VC,Int>(g) );
            *A_STAR_VC = A;

            std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
            ( new DistMatrix<T,STAR,VR,Int>(true,this->RowAlignment(),g) );
            *A_STAR_VR = *A_STAR_VC;
            delete A_STAR_VC.release(); // lowers memory highwater

            *this = *A_STAR_VR;
            this->ResizeTo( A_STAR_VR->Height(), A_STAR_VR->Width() );
        }
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [MR,* ]");
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
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();

    std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC,Int>(g) );
    *A_STAR_VC = A;

    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR,Int>(true,this->RowAlignment(),g) );
    *A_STAR_VR = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    *this = *A_STAR_VR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [VC,* ]");
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
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for  
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuf[k*portionSize];
            const Int thisRowShift = Shift_(k,rowAlignment,c);
            const Int thisLocalWidth = Length_(width,thisRowShift,c);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                MemCopy
                ( &data[jLoc*localHeightA], 
                  &ABuffer[(thisRowShift+jLoc*c)*ALDim], localHeightA );
        }

        // Communicate
        mpi::AllToAll
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisRank = row+k*r;
            const Int thisColShift = Shift_(thisRank,colAlignmentA,p);
            const Int thisColOffset = (thisColShift-colShift) / r;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int thisRowShift = Shift_(k,rowAlignment,c);
            const Int thisLocalWidth = Length_(width,thisRowShift,c);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                MemCopy
                ( &data[jLoc*localHeightA], 
                  &ABuffer[(thisRowShift+jLoc*c)*ALDim], localHeightA );
        }

        // SendRecv: properly align A[VC,*] via a trade in the column
        mpi::SendRecv
        ( secondBuf, c*portionSize, sendRow, 0,
          firstBuf,  c*portionSize, recvRow, 0, g.ColComm() );

        // AllToAll to gather all of the aligned A[VC,*] data into 
        // secondBuff.
        mpi::AllToAll
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int thisRank = recvRow+k*r;
            const Int thisColShift = Shift_(thisRank,colAlignmentA,p);
            const Int thisColOffset = (thisColShift-colShift) / r;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,VC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR,Int> A_STAR_VR(true,this->RowAlignment(),g);

    A_STAR_VR = A;
    *this = A_STAR_VR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [VR,* ]");
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
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,VR]");
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
            this->rowAlignment_ = A.RowAlignment() % g.Width();
            this->SetRowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
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
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &sendBuf[k*portionSize];
            const Int thisColShift = Shift_(k,colAlignment,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisRank = col+k*c;
            const Int thisRowShift = Shift_(thisRank,rowAlignmentA,p);
            const Int thisRowOffset = (thisRowShift-rowShift) / c;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int thisColShift = Shift_(k,colAlignment,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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
        ( secondBuf, r*portionSize, sendCol, 0,
          firstBuf,  r*portionSize, recvCol, 0, g.RowComm() );

        // AllToAll to gather all of the aligned [*,VR] data into 
        // secondBuf
        mpi::AllToAll
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int thisRank = recvCol+k*c;
            const Int thisRowShift = Shift_(thisRank,rowAlignmentA,p);
            const Int thisRowOffset = (thisRowShift-rowShift) / c;
            const Int thisLocalWidth = Length_(width,thisRowShift,p);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                MemCopy
                ( &thisBuffer[(thisRowOffset+jLoc*r)*thisLDim], 
                  &data[jLoc*localHeight], localHeight );
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [* ,* ]");
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
    const Int c = this->Grid().Width();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();

    const Int ALDim = A.LDim();
    const Int thisLDim = this->LDim();
    T* thisBuffer = this->Buffer();
    const T* ABuffer = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        T* destCol = &thisBuffer[jLoc*thisLDim];
        const T* sourceCol = &ABuffer[colShift+(rowShift+jLoc*c)*ALDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            destCol[iLoc] = sourceCol[iLoc*r];
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,CIRC,CIRC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR] = [o ,o ]");
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
        const T* ABuffer = A.LockedBuffer();
        for( int t=0; t<rowStride; ++t )
        {
            const int tLocalWidth = Length( n, t, rowStride );
            const int col = (rowAlignment+t) % rowStride;
            for( int s=0; s<colStride; ++s )
            {
                const int sLocalHeight = Length( m, s, colStride );
                const int row = (colAlignment+s) % colStride;
                const int q = row + col*colStride;
                for( int jLoc=0; jLoc<tLocalWidth; ++jLoc ) 
                {
                    const int j = t + jLoc*rowStride;
                    for( int iLoc=0; iLoc<sLocalHeight; ++iLoc )
                    {
                        const int i = s + iLoc*colStride;
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
DistMatrix<T,MC,MR,Int>::SumScatterFrom( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SumScatterFrom([MC,* ])");
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
            ( sendBuf, recvBuf, sendSize, 
              mpi::SUM, rowAlignment, g.RowComm() );

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

            T* buffer = this->auxMemory_.Require( sendSize );
            
            // Pack 
            const Int ALDim = A.LDim();
            const T* ABuffer = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                T* data = &buffer[k*recvSize];
                const Int thisRowShift = Shift_( k, rowAlignment, c );
                const Int thisLocalWidth = Length_(width,thisRowShift,c);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                    MemCopy
                    ( &data[jLoc*localHeight], 
                      &ABuffer[(thisRowShift+jLoc*c)*ALDim], localHeight );
            }

            // Communicate
            mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.RowComm() );

            // Unpack our received data
            T* thisBuffer = this->Buffer();
            const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
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
            ( sendBuf, recvBuf, portionSize, 
              mpi::SUM, rowAlignment, g.RowComm() );

            if( myCol == rowAlignment )
            {
                // Perform the realignment
                mpi::SendRecv
                ( recvBuf, portionSize, sendRow, 0,
                  sendBuf, portionSize, recvRow, 0, g.ColComm() );
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                    MemCopy
                    ( &data[jLoc*localHeightA], 
                      &ABuffer[(thisRowShift+jLoc*c)*ALDim], localHeightA );
            }

            // Reduce-scatter over each process row
            mpi::ReduceScatter
            ( secondBuf, firstBuf, recvSize_RS, mpi::SUM, g.RowComm() );

            // Trade reduced data with the appropriate process row
            mpi::SendRecv
            ( firstBuf,  localHeightA*localWidth, sendRow, 0,
              secondBuf, localHeight*localWidth,    recvRow, 0, 
              g.ColComm() );

            // Unpack the received data
            T* thisBuffer = this->Buffer();
            const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                MemCopy
                ( &thisBuffer[jLoc*thisLDim], 
                  &secondBuf[jLoc*localHeight], localHeight );
            this->auxMemory_.Release();
        }
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::SumScatterFrom( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SumScatterFrom([* ,MR])");
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
        const Int r = g.Height();
        const Int colAlignment = this->ColAlignment();

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalHeight = MaxLength(height,r);

        const Int recvSize = mpi::Pad( maxLocalHeight*localWidth );
        const Int sendSize = r * recvSize;

        T* buffer = this->auxMemory_.Require( sendSize );

        // Pack 
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisColShift = Shift_(k,colAlignment,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.ColComm() );

        // Unpack our received data
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
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
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &secondBuf[k*recvSize_RS];
            const Int thisColShift = Shift_(k,colAlignment,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Reduce-scatter over each process col
        mpi::ReduceScatter
        ( secondBuf, firstBuf, recvSize_RS, mpi::SUM, g.ColComm() );

        // Trade reduced data with the appropriate process col
        mpi::SendRecv
        ( firstBuf,  localHeight*localWidthA, sendCol, 0,
          secondBuf, localHeight*localWidth,    recvCol, 0, g.RowComm() );

        // Unpack the received data
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( &thisBuffer[jLoc*thisLDim], 
              &secondBuf[jLoc*localHeight], localHeight );
        this->auxMemory_.Release();
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::SumScatterFrom( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SumScatterFrom([* ,* ])");
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
    const Int maxLocalHeight = MaxLength(height,r);
    const Int maxLocalWidth = MaxLength(width,c);

    const Int recvSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
    const Int sendSize = r*c*recvSize;

    T* buffer = this->auxMemory_.Require( sendSize );

    // Pack 
    const Int ALDim = A.LDim();
    const T* ABuffer = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int l=0; l<c; ++l )
    {
        const Int thisRowShift = Shift_( l, rowAlignment, c );
        const Int thisLocalWidth = Length_( width, thisRowShift, c );

        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[(k+l*r)*recvSize];
            const Int thisColShift = Shift_(k,colAlignment,r);
            const Int thisLocalHeight = Length_(height,thisColShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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
    mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.VCComm() );

    // Unpack our received data
    T* thisBuffer = this->Buffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        MemCopy
        ( &thisBuffer[jLoc*thisLDim], 
          &buffer[jLoc*localHeight], localHeight );
    this->auxMemory_.Release();
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,MC,STAR,Int>& A )
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
            ( sendBuf, recvBuf, portionSize, 
              mpi::SUM, rowAlignment, g.RowComm() );

            if( myCol == rowAlignment )
            {
                T* thisCol = this->Buffer(0,0);
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
            const Int c = g.Width();
            const Int rowAlignment = this->RowAlignment();

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int maxLocalWidth = MaxLength(width,c);

            const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );
            const Int sendSize = c*portionSize;
            T* buffer = this->auxMemory_.Require( sendSize );

            // Pack 
            const Int ALDim = A.LDim();
            const T* ABuffer = A.LockedBuffer();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                T* data = &buffer[k*portionSize];
                const Int thisRowShift = Shift_( k, rowAlignment, c );
                const Int thisLocalWidth = Length_(width,thisRowShift,c);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* ACol = &ABuffer[(thisRowShift+jLoc*c)*ALDim];
                    T* dataCol = &data[jLoc*localHeight];
                    MemCopy( dataCol, ACol, localHeight );
                }
            }
            
            // Communicate
            mpi::ReduceScatter( buffer, portionSize, mpi::SUM, g.RowComm() );

            // Update with our received data
            T* thisBuffer = this->Buffer();
            const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
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
            ( sendBuf, recvBuf, portionSize, 
              mpi::SUM, rowAlignment, g.RowComm() );

            if( myCol == rowAlignment )
            {
                // Perform the realignment
                mpi::SendRecv
                ( recvBuf, portionSize, sendRow, 0,
                  sendBuf, portionSize, recvRow, 0, g.ColComm() );

                T* thisCol = this->Buffer(0,0);
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
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                T* data = &secondBuf[k*recvSize_RS];
                const Int thisRowShift = Shift_( k, rowAlignment, c );
                const Int thisLocalWidth = Length_(width,thisRowShift,c);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                {
                    const T* ACol = &ABuffer[(thisRowShift+jLoc*c)*ALDim];
                    T* dataCol = &data[jLoc*localHeightA];
                    MemCopy( dataCol, ACol, localHeightA );
                }
            }

            // Reduce-scatter over each process row
            mpi::ReduceScatter
            ( secondBuf, firstBuf, recvSize_RS, mpi::SUM, g.RowComm() );

            // Trade reduced data with the appropriate process row
            mpi::SendRecv
            ( firstBuf,  localHeightA*localWidth, sendRow, 0,
              secondBuf, localHeight*localWidth,    recvRow, 0, 
              g.ColComm() );

            // Update with our received data
            T* thisBuffer = this->Buffer();
            const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(AVOID_OMP_FMA)
#pragma omp parallel for
#endif
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

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MR,Int>& A )
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
        T* buffer = this->auxMemory_.Require( sendSize );

        // Pack 
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisColShift = Shift_( k, colAlignment, r );
            const Int thisLocalHeight = Length_(height,thisColShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.ColComm() );

        // Update with our received data
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(AVOID_OMP_FMA)
#pragma omp parallel for
#endif
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
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &secondBuf[k*recvSize_RS];
            const Int thisColShift = Shift_( k, colAlignment, r );
            const Int thisLocalHeight = Length_(height,thisColShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuffer[thisColShift+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Reduce-scatter over each process col
        mpi::ReduceScatter
        ( secondBuf, firstBuf, recvSize_RS, mpi::SUM, g.ColComm() );

        // Trade reduced data with the appropriate process col
        mpi::SendRecv
        ( firstBuf,  localHeight*localWidthA, sendCol, 0,
          secondBuf, localHeight*localWidth,  recvCol, mpi::ANY_TAG,
          g.RowComm() );

        // Update with our received data
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(AVOID_OMP_FMA)
#pragma omp parallel for
#endif
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

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SumScatterUpdate([* ,* ])");
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
    const Int maxLocalHeight = MaxLength(height,r);
    const Int maxLocalWidth = MaxLength(width,c);

    const Int recvSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
    const Int sendSize = r * c * recvSize;
    T* buffer = this->auxMemory_.Require( sendSize );

    // Pack 
    const T* ABuffer = A.LockedBuffer();
    const Int ALDim = A.LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int l=0; l<c; ++l )
    {
        const Int thisRowShift = Shift_( l, rowAlignment, c );
        const Int thisLocalWidth = Length_( width, thisRowShift, c );
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[(k+l*r)*recvSize];
            const Int thisColShift = Shift_( k, colAlignment, r );
            const Int thisLocalHeight = Length_(height,thisColShift,r);
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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
    mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.VCComm() );

    // Unpack our received data
    T* thisBuffer = this->Buffer();
    const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(AVOID_OMP_FMA)
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* bufferCol = &buffer[jLoc*localHeight];
        T* thisCol = &thisBuffer[jLoc*thisLDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            thisCol[iLoc] += alpha*bufferCol[iLoc];
    }
    this->auxMemory_.Release();
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::GetRealPartOfDiagonal
( DistMatrix<BASE(T),MD,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetRealPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength( offset );
#ifndef RELEASE
    if( d.Viewing() && length != d.Height() )
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

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( *this, offset );
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        dBuf[k] = RealPart(thisBuffer[iLoc+jLoc*thisLDim]);
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::GetImagPartOfDiagonal
( DistMatrix<BASE(T),MD,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetImagPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength( offset );
#ifndef RELEASE
    if( d.Viewing() && length != d.Height() )
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

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( *this, offset );
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        dBuf[k] = ImagPart(thisBuffer[iLoc+jLoc*thisLDim]);
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::GetRealPartOfDiagonal
( DistMatrix<BASE(T),STAR,MD,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetRealPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength( offset );
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        dBuf[k*dLDim] = RealPart(thisBuffer[iLoc+jLoc*thisLDim]);
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::GetImagPartOfDiagonal
( DistMatrix<BASE(T),STAR,MD,Int>& d, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::GetImagPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d.Grid() );
#endif
    const Int length = this->DiagonalLength( offset );
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        dBuf[k*dLDim] = ImagPart(thisBuffer[iLoc+jLoc*thisLDim]);
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::SetRealPartOfDiagonal
( const DistMatrix<BASE(T),MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetRealPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const Int length = this->DiagonalLength( offset );
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        this->SetLocalRealPart( iLoc, jLoc, dBuf[k] );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::SetImagPartOfDiagonal
( const DistMatrix<BASE(T),MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetImagPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const Int length = this->DiagonalLength( offset );
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        this->SetLocalImagPart( iLoc, jLoc, dBuf[k] );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::SetRealPartOfDiagonal
( const DistMatrix<BASE(T),STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetRealPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const Int length = this->DiagonalLength( offset );
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        this->SetLocalRealPart( iLoc, jLoc, dBuf[k*dLDim] );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MC,MR,Int>::SetImagPartOfDiagonal
( const DistMatrix<BASE(T),STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MC,MR]::SetImagPartOfDiagonal");
    this->AssertSameGrid( d.Grid() );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const Int length = this->DiagonalLength( offset );
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(lcm/r);
        const Int jLoc = jLocStart + k*(lcm/c);
        this->SetLocalImagPart( iLoc, jLoc, dBuf[k*dLDim] );
    }
}

template class DistMatrix_Dist<MC,MR,int>;

#define PROTO(T) \
  template class DistMatrix<T,MC,MR,int>
#define COPY(T,CD,RD) \
  template DistMatrix<T,MC,MR,int>::DistMatrix( \
    const DistMatrix<T,CD,RD,int>& A )
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
