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
DistMatrix_Dist<STAR,MD,Int>::DistMatrix_Dist( const elem::Grid& g )
: DistMatrix_Base<Int>(g), diagPath_(0)
{ this->SetShifts(); }

template<typename Int>
DistMatrix_Dist<STAR,MD,Int>::DistMatrix_Dist( const elem::Grid& g, Int rowAlignmentVC )
: DistMatrix_Base<Int>(g), diagPath_(g.DiagPathRank(rowAlignmentVC))
{ this->Align( 0, diagPath_ ); }

template <typename Int>
elem::Distribution
DistMatrix_Dist<STAR,MD,Int>::ColDist() const { return STAR; }

template <typename Int>
elem::Distribution
DistMatrix_Dist<STAR,MD,Int>::RowDist() const { return MD; }

template<typename Int>
Int
DistMatrix_Dist<STAR,MD,Int>::ColStride() const
{ return 1; }

template<typename Int>
Int
DistMatrix_Dist<STAR,MD,Int>::RowStride() const
{ return this->grid_->LCM(); }

template<typename Int>
Int
DistMatrix_Dist<STAR,MD,Int>::ColRank() const
{ return 0; }

template<typename Int>
Int
DistMatrix_Dist<STAR,MD,Int>::RowRank() const
{ return this->grid_->DiagPathRank(); }

template<typename Int>
bool
DistMatrix_Dist<STAR,MD,Int>::Participating() const
{
    const Grid& g = this->Grid();
    return ( g.InGrid() && g.DiagPath()==diagPath_ );
}

template<typename Int>
Int
DistMatrix_Dist<STAR,MD,Int>::DiagPath() const
{ return diagPath_; }

template <typename Int>
void 
DistMatrix_Dist<STAR,MD,Int>::Attach
( Int height, Int width, Int rowAlignmentVC, void* buffer, Int ldim, const elem::Grid& g )
{ 
    diagPath_ = this->Grid().DiagPath( rowAlignmentVC );
    DistMatrix_Base<Int>::Attach( height, width, 0, diagPath_, buffer, ldim, g ); 
}

template <typename Int>
void 
DistMatrix_Dist<STAR,MD,Int>::LockedAttach
( Int height, Int width, Int rowAlignmentVC, const void* buffer, Int ldim, const elem::Grid& g )
{ 
    diagPath_ = this->Grid().DiagPath( rowAlignmentVC );
    DistMatrix_Base<Int>::LockedAttach( height, width, 0, diagPath_, buffer, ldim, g ); 
}

template<typename Int>
void
DistMatrix_Dist<STAR,MD,Int>::AlignWith( const DistMatrix_Base<Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::AlignWith");
#endif
    this->SetGrid( A.Grid() );
    elem::Distribution CD = A.ColDist(), RD = A.RowDist();

    if( CD == MD && RD == STAR )
    {
        this->rowAlignment_ = A.colAlignment_;
        diagPath_ = A.DiagPath();
    }
    else if( CD == STAR && RD == MD )
    {
        this->rowAlignment_ = A.rowAlignment_;
        diagPath_ = A.DiagPath();
    }
#ifndef RELEASE
    else throw std::logic_error("Invalid alignment");
#endif
    this->constrainedRowAlignment_ = true;
    this->SetShifts();
}

template<typename Int>
void
DistMatrix_Dist<STAR,MD,Int>::AlignRowsWith( const DistMatrix_Base<Int>& A )
{ this->AlignWith( A ); }

template<typename Int>
bool
DistMatrix_Dist<STAR,MD,Int>::AlignedWithDiagonal
( const DistMatrix_Base<Int>& A, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::AlignedWithDiagonal");
#endif
    const elem::Grid& grid = this->Grid();
    if( grid != A.Grid() )
        return false;
    elem::Distribution CD = A.ColDist(), RD = A.RowDist();

    bool aligned;
    const Int r = grid.Height();
    const Int c = grid.Width();
    const Int firstDiagRow = 0;
    const Int firstDiagCol = diagPath_;
    const Int diagRow = (firstDiagRow+this->RowAlignment()) % r;
    const Int diagCol = (firstDiagCol+this->RowAlignment()) % c;
    if( CD == MC && RD == MR )
    {
        if( offset >= 0 )
        {
            const Int ownerRow = A.colAlignment_;
            const Int ownerCol = (A.rowAlignment_ + offset) % c;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
        else
        {
            const Int ownerRow = (A.colAlignment_-offset) % r;
            const Int ownerCol = A.rowAlignment_;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
    }
    else if( CD == MR && RD == MC )
    {
        if( offset >= 0 )
        {
            const Int ownerCol = A.colAlignment_;
            const Int ownerRow = (A.rowAlignment_ + offset) % r;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
        else
        {
            const Int ownerCol = (A.colAlignment_-offset) % c;
            const Int ownerRow = A.rowAlignment_;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
    }
    else if( CD == MD && RD == STAR )
    {
        aligned = ( diagPath_==A.DiagPath() &&
                    this->rowAlignment_==A.colAlignment_ );
    }
    else if( CD == STAR && RD == MD )
    {
        aligned = ( diagPath_==A.DiagPath() &&
                    this->rowAlignment_==A.rowAlignment_ );
    }
    else aligned = false;
    return aligned;
}

template<typename Int>
void
DistMatrix_Dist<STAR,MD,Int>::AlignWithDiagonal
( const DistMatrix_Base<Int>& A, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::AlignWithDiagonal");
#endif
    const elem::Grid& grid = this->Grid();
    this->SetGrid( grid );
    elem::Distribution CD = A.ColDist(), RD = A.RowDist();

    const Int r = grid.Height();
    const Int c = grid.Width();
    if( CD == MC && RD == MR )
    {
        Int owner;
        if( offset >= 0 )
        {
            const Int ownerRow = A.colAlignment_;
            const Int ownerCol = (A.rowAlignment_ + offset) % c;
            owner = ownerRow + r*ownerCol;
        }
        else
        {
            const Int ownerRow = (A.colAlignment_-offset) % r;
            const Int ownerCol = A.rowAlignment_;
            owner = ownerRow + r*ownerCol;
        }
        diagPath_ = grid.DiagPath(owner);
        this->rowAlignment_ = grid.DiagPathRank(owner);
    }
    else if( CD == MR && RD == MC )
    {
        Int owner;
        if( offset >= 0 )
        {
            const Int ownerCol = A.colAlignment_;
            const Int ownerRow = (A.rowAlignment_ + offset) % r;
            owner = ownerRow + r*ownerCol;
        }
        else
        {
            const Int ownerCol = (A.colAlignment_-offset) % c;
            const Int ownerRow = A.rowAlignment_;
            owner = ownerRow + r*ownerCol;
        }
        diagPath_ = grid.DiagPath(owner);
        this->rowAlignment_ = grid.DiagPathRank(owner);
    }
    else if( CD == MD && RD == STAR )
    {
        diagPath_ = A.DiagPath();
        this->rowAlignment_ = A.colAlignment_;
    }
    else if( CD == STAR && RD == MD )
    {
        diagPath_ = A.DiagPath();
        this->rowAlignment_ = A.rowAlignment_;
    }
#ifndef RELEASE
    else throw std::logic_error("Nonsensical AlignWithDiagonal");
#endif
    this->constrainedRowAlignment_ = true;
    this->SetShifts();
}

template<typename Int>
bool
DistMatrix_Dist<STAR,MD,Int>::Index( Int i, Int j, Int& iLocal, Int& jLocal, int& mpiSrc, mpi::Comm& mpiDst ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::Index");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (j + this->rowAlignment_) % r;
    const Int ownerCol = (j + this->rowAlignment_ + diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;
    mpiSrc = g.VCToViewingMap(ownerRank);
    mpiDst = g.ViewingComm();
    if( g.VCRank() != ownerRank ) return false;
    iLocal = i;
    jLocal = (j-this->RowShift()) / g.LCM();
    return true;
}

template<typename Int>
void
DistMatrix_Dist<STAR,MD,Int>::MakeConsistent()
{
#ifndef RELEASE
    CallStackEntry cse("[* ,MD]::MakeConsistent");
#endif
    const elem::Grid& g = this->Grid();
    const int root = g.VCToViewingMap(0);
    int message[6];
    if( g.ViewingRank() == root )
    {   
        message[0] = this->viewType_;
        message[1] = this->height_;
        message[2] = this->width_;
        message[3] = this->constrainedRowAlignment_;
        message[4] = this->rowAlignment_;
        message[5] = this->diagPath_;
    }
    mpi::Broadcast( message, 6, root, g.ViewingComm() );
    const ViewType newViewType = static_cast<ViewType>(message[0]);
    const int newHeight = message[1];
    const int newWidth = message[2];
    const bool newConstrainedRow = message[3];
    const int newRowAlignment = message[4];
    const int newDiagPath = message[5];
    if( !this->Participating() )
    {
        this->viewType_ = newViewType;
        this->height_ = newHeight;
        this->width_ = newWidth;
        this->constrainedRowAlignment_ = newConstrainedRow;
        this->rowAlignment_ = newRowAlignment;
        this->diagPath_ = newDiagPath;
        this->constrainedColAlignment_ = 0;
        this->colAlignment_ = 0;
        this->colShift_ = 0;
        this->rowShift_ = 0;
    }
#ifndef RELEASE
    else
    {
        if( this->viewType_ != newViewType )
            throw std::logic_error("Inconsistent ViewType");
        if( this->height_ != newHeight )
            throw std::logic_error("Inconsistent height");
        if( this->width_ != newWidth )
            throw std::logic_error("Inconsistent width");
        if( this->constrainedRowAlignment_ != newConstrainedRow ||
            this->rowAlignment_ != newRowAlignment )
            throw std::logic_error("Inconsistent row constraint");
        if( this->diagPath_ != newDiagPath )
            throw std::logic_error("Inconsistent diagonal path");
    }
#endif
}

/*
 * DistMatrix
 */

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix( const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<STAR,MD,Int>(g), DistMatrix_Type<T,Int>(g)
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<STAR,MD,Int>(g), DistMatrix_Type<T,Int>(g)
{ this->ResizeTo(height,width); }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignmentVC, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<STAR,MD,Int>(g,rowAlignmentVC), DistMatrix_Type<T,Int>(g)
{ this->ResizeTo( height, width ); }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignmentVC, Int ldim, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<STAR,MD,Int>(g,rowAlignmentVC), DistMatrix_Type<T,Int>(g)
{ this->ResizeTo( height, width, ldim ); }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignmentVC, const T* buffer, Int ldim, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<STAR,MD,Int>(g,rowAlignmentVC), DistMatrix_Type<T,Int>(g)
{ this->LockedAttach(height,width,rowAlignmentVC,buffer,ldim,g); }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignmentVC, T* buffer, Int ldim, const elem::Grid& g )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<STAR,MD,Int>(g,rowAlignmentVC), DistMatrix_Type<T,Int>(g)
{ this->Attach(height,width,rowAlignmentVC,buffer,ldim,g); }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix( const DistMatrix<T,STAR,MD,Int>& A )
: DistMatrix_Base<Int>(A.Grid()), DistMatrix_Dist<STAR,MD,Int>(A.Grid()), DistMatrix_Type<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,MD]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MD] with itself");
}

template<typename T,typename Int>
template<Distribution U,Distribution V>
DistMatrix<T,STAR,MD,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: DistMatrix_Base<Int>(A.Grid()), DistMatrix_Dist<STAR,MD,Int>(A.Grid()), DistMatrix_Type<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,MD]::DistMatrix");
#endif
    this->SetShifts();
    if( STAR != U || MD != V || reinterpret_cast<const DistMatrix_Base<Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MD] with itself");
}

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::~DistMatrix()
{ }

//
// Utility functions, e.g., operator=
//

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MD] = [MC,MR] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MD] = [MC,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MD] = [* ,MR] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MD] = [MD,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->diagPath_ = A.DiagPath();
            this->rowAlignment_ = A.rowAlignment_;
            if( this->Participating() )
                this->rowShift_ = A.RowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->diagPath_ == A.DiagPath() && 
        this->rowAlignment_ == A.rowAlignment_ )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned [* ,MD] <- [* ,MD]." << std::endl;
#endif
        throw std::logic_error
        ("Unaligned [* ,MD] = [* ,MD] not yet implemented");
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MD] = [MR,MC] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MD] = [MR,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [* ,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MD] = [* ,MC] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MD] = [VC,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [* ,VC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MD] = [* ,VC] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MD] = [VR,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[* ,MD] = [* ,VR] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int lcm = this->Grid().LCM();
    const Int rowShift = this->RowShift();

    const Int height = this->Height();
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
        const T* ACol = &ABuf[(rowShift+jLoc*lcm)*ALDim];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        MemCopy( thisCol, ACol, height );
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,CIRC,CIRC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD] = [o ,o ]");
#endif
    DistMatrix<T,MC,MR> A_MC_MR( A.Grid() );
    A_MC_MR.AlignWith( *this );
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}


template class DistMatrix_Dist<STAR,MD,int>;

#define PROTO(T) \
  template class DistMatrix<T,STAR,MD,int>
#define COPY(T,CD,RD) \
  template DistMatrix<T,STAR,MD,int>::DistMatrix( \
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
