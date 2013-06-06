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
DistMatrix<T,STAR,MD,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   0,(g.InGrid() && g.DiagPath()==0 ? g.DiagPathRank() : 0),
   0,0,g),
  diagPath_(0)
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   0,
   (g.InGrid() && g.DiagPath()==0 ? g.DiagPathRank() : 0),height,
   (g.InGrid() && g.DiagPath()==0 ? 
    Length(width,g.DiagPathRank(),0,g.LCM()) : 0),g),
  diagPath_(0)
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignmentVC, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,true,0,g.DiagPathRank(rowAlignmentVC),
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignmentVC) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignmentVC),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignmentVC) ?
    Length(width,g.DiagPathRank(),g.DiagPathRank(rowAlignmentVC),g.LCM()) :
    0),g),
  diagPath_(g.DiagPath(rowAlignmentVC))
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignmentVC, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,true,0,g.DiagPathRank(rowAlignmentVC),
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignmentVC) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignmentVC),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignmentVC) ?
    Length(width,g.DiagPathRank(),g.DiagPathRank(rowAlignmentVC),g.LCM()) :
    0),ldim,g),
  diagPath_(g.DiagPath(rowAlignmentVC))
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignmentVC, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,g.DiagPathRank(rowAlignmentVC),
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignmentVC) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignmentVC),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignmentVC) ?
    Length(width,g.DiagPathRank(),g.DiagPathRank(rowAlignmentVC),g.LCM()) :
    0),buffer,ldim,g),
  diagPath_(g.DiagPath(rowAlignmentVC))
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignmentVC, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,g.DiagPathRank(rowAlignmentVC),
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignmentVC) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignmentVC),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignmentVC) ?
    Length(width,g.DiagPathRank(),g.DiagPathRank(rowAlignmentVC),g.LCM()) :
    0),buffer,ldim,g),
  diagPath_(g.DiagPath(rowAlignmentVC))
{ }

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::DistMatrix( const DistMatrix<T,STAR,MD,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  0,(A.Participating() ? A.RowRank() : 0),
  0,0,A.Grid()),
  diagPath_(0)
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,MD]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MD] with itself");
}

template<typename T,typename Int>
template<Distribution U,Distribution V>
DistMatrix<T,STAR,MD,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  0,(A.Participating() ? A.RowRank() : 0),
  0,0,A.Grid()),
  diagPath_(0)
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[* ,MD]::DistMatrix");
#endif
    if( STAR != U || MD != V || 
        reinterpret_cast<const DistMatrix<T,STAR,MD,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MD] with itself");
}

template<typename T,typename Int>
DistMatrix<T,STAR,MD,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
elem::DistData<Int>
DistMatrix<T,STAR,MD,Int>::DistData() const
{
    elem::DistData<Int> data;
    data.colDist = STAR;
    data.rowDist = MD;
    data.colAlignment = 0;
    data.rowAlignment = this->rowAlignment_;
    data.diagPath = this->diagPath_;
    data.grid = this->grid_;
    return data;
}

template<typename T,typename Int>
Int
DistMatrix<T,STAR,MD,Int>::ColStride() const
{ return 1; }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,MD,Int>::RowStride() const
{ return this->grid_->LCM(); }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,MD,Int>::ColRank() const
{ return 0; }

template<typename T,typename Int>
Int
DistMatrix<T,STAR,MD,Int>::RowRank() const
{ return this->grid_->DiagPathRank(); }

template<typename T,typename Int>
bool
DistMatrix<T,STAR,MD,Int>::Participating() const
{
    const Grid& g = this->Grid();
    return ( g.InGrid() && g.DiagPath()==this->diagPath_ );
}

template<typename T,typename Int>
Int
DistMatrix<T,STAR,MD,Int>::DiagPath() const
{ return this->diagPath_; }

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::AlignWith( const elem::DistData<Int>& data )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::AlignWith");
    this->AssertFreeRowAlignment();
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( data.colDist == MD && data.rowDist == STAR )
    {
        this->rowAlignment_ = data.colAlignment;
        this->diagPath_ = data.diagPath;
    }
    else if( data.colDist == STAR && data.rowDist == MD )
    {
        this->rowAlignment_ = data.rowAlignment;
        this->diagPath_ = data.diagPath;
    }
#ifndef RELEASE
    else throw std::logic_error("Invalid alignment");
#endif
    this->constrainedRowAlignment_ = true;
    this->SetShifts();
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::AlignWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::AlignRowsWith( const elem::DistData<Int>& data )
{ this->AlignWith( data ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::AlignRowsWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T,typename Int>
bool
DistMatrix<T,STAR,MD,Int>::AlignedWithDiagonal
( const elem::DistData<Int>& data, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::AlignedWithDiagonal");
#endif
    const Grid& grid = this->Grid();
    if( grid != *data.grid )
        return false;

    bool aligned;
    const Int r = grid.Height();
    const Int c = grid.Width();
    const Int firstDiagRow = 0;
    const Int firstDiagCol = this->diagPath_;
    const Int diagRow = (firstDiagRow+this->RowAlignment()) % r;
    const Int diagCol = (firstDiagCol+this->RowAlignment()) % c;
    if( data.colDist == MC && data.rowDist == MR )
    {
        if( offset >= 0 )
        {
            const Int ownerRow = data.colAlignment;
            const Int ownerCol = (data.rowAlignment + offset) % c;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
        else
        {
            const Int ownerRow = (data.colAlignment-offset) % r;
            const Int ownerCol = data.rowAlignment;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
    }
    else if( data.colDist == MR && data.rowDist == MC )
    {
        if( offset >= 0 )
        {
            const Int ownerCol = data.colAlignment;
            const Int ownerRow = (data.rowAlignment + offset) % r;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
        else
        {
            const Int ownerCol = (data.colAlignment-offset) % c;
            const Int ownerRow = data.rowAlignment;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
    }
    else if( data.colDist == MD && data.rowDist == STAR )
    {
        aligned = ( this->diagPath_==data.diagPath &&
                    this->rowAlignment_==data.colAlignment );
    }
    else if( data.colDist == STAR && data.rowDist == MD )
    {
        aligned = ( this->diagPath_==data.diagPath &&
                    this->rowAlignment_==data.rowAlignment );
    }
    else aligned = false;
    return aligned;
}

template<typename T,typename Int>
bool
DistMatrix<T,STAR,MD,Int>::AlignedWithDiagonal
( const AbstractDistMatrix<T,Int>& A, Int offset ) const
{ return this->AlignedWithDiagonal( A.DistData(), offset ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::AlignWithDiagonal
( const elem::DistData<Int>& data, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::AlignWithDiagonal");
    this->AssertFreeRowAlignment();
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    const Int r = grid.Height();
    const Int c = grid.Width();
    if( data.colDist == MC && data.rowDist == MR )
    {
        Int owner;
        if( offset >= 0 )
        {
            const Int ownerRow = data.colAlignment;
            const Int ownerCol = (data.rowAlignment + offset) % c;
            owner = ownerRow + r*ownerCol;
        }
        else
        {
            const Int ownerRow = (data.colAlignment-offset) % r;
            const Int ownerCol = data.rowAlignment;
            owner = ownerRow + r*ownerCol;
        }
        this->diagPath_ = grid.DiagPath(owner);
        this->rowAlignment_ = grid.DiagPathRank(owner);
    }
    else if( data.colDist == MR && data.rowDist == MC )
    {
        Int owner;
        if( offset >= 0 )
        {
            const Int ownerCol = data.colAlignment;
            const Int ownerRow = (data.rowAlignment + offset) % r;
            owner = ownerRow + r*ownerCol;
        }
        else
        {
            const Int ownerCol = (data.colAlignment-offset) % c;
            const Int ownerRow = data.rowAlignment;
            owner = ownerRow + r*ownerCol;
        }
        this->diagPath_ = grid.DiagPath(owner);
        this->rowAlignment_ = grid.DiagPathRank(owner);
    }
    else if( data.colDist == MD && data.rowDist == STAR )
    {
        this->diagPath_ = data.diagPath;
        this->rowAlignment_ = data.colAlignment;
    }
    else if( data.colDist == STAR && data.rowDist == MD )
    {
        this->diagPath_ = data.diagPath;
        this->rowAlignment_ = data.rowAlignment;
    }
#ifndef RELEASE
    else throw std::logic_error("Nonsensical AlignWithDiagonal");
#endif
    this->constrainedRowAlignment_ = true;
    this->SetShifts();
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::AlignWithDiagonal
( const AbstractDistMatrix<T,Int>& A, Int offset )
{ this->AlignWithDiagonal( A.DistData(), offset ); }

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::PrintBase");
#endif
    if( this->Grid().Rank() == 0 && msg != "" )
        os << msg << std::endl;
        
    const Int height     = this->Height();
    const Int width      = this->Width();
    const Int localWidth = this->LocalWidth();
    const Int lcm        = this->Grid().LCM();

    if( height == 0 || width == 0 || !this->Grid().InGrid() )
        return;

    std::vector<T> sendBuf(height*width,0);
    if( this->Participating() )
    {
        const Int colShift = this->ColShift();
        const T* thisBuffer = this->LockedBuffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            T* destCol = &sendBuf[colShift+jLoc*lcm*height];
            const T* sourceCol = &thisBuffer[jLoc*thisLDim];
            for( Int i=0; i<height; ++i )
                destCol[i] = sourceCol[i];
        }
    }

    // If we are the root, allocate a receive buffer
    std::vector<T> recvBuf;
    if( this->Grid().Rank() == 0 )
        recvBuf.resize( height*width );

    // Sum the contributions and send to the root
    mpi::Reduce
    ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, 
      this->Grid().Comm() );

    if( this->Grid().Rank() == 0 )
    {
        // Print the data
        for( Int i=0; i<height; ++i )
        {
            for( Int j=0; j<width; ++j )
                os << recvBuf[i+j*height] << " ";
            os << "\n";
        }
        os << std::endl;
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::Attach
( Int height, Int width, Int rowAlignmentVC,
  T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::Attach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->diagPath_ = grid.DiagPath(rowAlignmentVC);
    this->rowAlignment_ = grid.DiagPathRank(rowAlignmentVC);
    this->viewing_ = true;
    this->SetRowShift();
    if( this->Participating() )
    {
        const Int localWidth = Length(width,this->rowShift_,grid.LCM());
        this->matrix_.Attach( height, localWidth, buffer, ldim );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::LockedAttach
( Int height, Int width, Int rowAlignmentVC,
  const T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->diagPath_ = grid.DiagPath(rowAlignmentVC);
    this->rowAlignment_ = grid.DiagPathRank(rowAlignmentVC);
    this->viewing_ = true;
    this->locked_ = true;
    this->SetRowShift();
    if( this->Participating() )
    {
        const Int localWidth = Length(width,this->rowShift_,grid.LCM());
        this->matrix_.LockedAttach( height, localWidth, buffer, ldim );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
    {
        const Int lcm = this->Grid().LCM();
        this->matrix_.ResizeTo( height, Length(width,this->RowShift(),lcm) );
    }
}

template<typename T,typename Int>
T
DistMatrix<T,STAR,MD,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (j + this->rowAlignment_) % r;
    const Int ownerCol = (j + this->rowAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    T u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
        u = this->GetLocal(i,jLoc);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
    return u;
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (j + this->rowAlignment_) % r;
    const Int ownerCol = (j + this->rowAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
        this->SetLocal(i,jLoc,u);
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (j + this->rowAlignment_) % r;
    const Int ownerCol = (j + this->rowAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
        this->UpdateLocal(i,jLoc,u);
    }
}

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
            this->diagPath_ = A.diagPath_;
            this->rowAlignment_ = A.rowAlignment_;
            if( this->Participating() )
                this->rowShift_ = A.RowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->diagPath_ == A.diagPath_ && 
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

    if( this->Participating() )
    {
        const Int lcm = this->Grid().LCM();
        const Int rowShift = this->RowShift();

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();

        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        const T* ABuffer = A.LockedBuffer();
        const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* ACol = &ABuffer[(rowShift+jLoc*lcm)*ALDim];
            T* thisCol = &thisBuffer[jLoc*thisLDim];
            MemCopy( thisCol, ACol, height );
        }
    }
    return *this;
}

//
// Routines which explicitly work in the complex plane
//

template<typename T,typename Int>
BASE(T)
DistMatrix<T,STAR,MD,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::GetRealPart");
    this->AssertValidEntry( i, j );
#endif
    typedef BASE(T) R;

    // We will determine the owner of entry (i,j) and broadcast from it
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (j + this->rowAlignment_) % r;
    const Int ownerCol = (j + this->rowAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
        u = this->GetLocalRealPart( i, jLoc );
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
    return u;
}

template<typename T,typename Int>
BASE(T)
DistMatrix<T,STAR,MD,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::GetImagPart");
    this->AssertValidEntry( i, j );
#endif
    typedef BASE(T) R;

    // We will determine the owner of entry (i,j) and broadcast from it
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (j + this->rowAlignment_) % r;
    const Int ownerCol = (j + this->rowAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
        u = this->GetLocalImagPart( i, jLoc );
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
    return u;
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (j + this->rowAlignment_) % r;
    const Int ownerCol = (j + this->rowAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
        this->SetLocalRealPart( i, jLoc, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (j + this->rowAlignment_) % r;
    const Int ownerCol = (j + this->rowAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
        this->SetLocalImagPart( i, jLoc, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (j + this->rowAlignment_) % r;
    const Int ownerCol = (j + this->rowAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
        this->UpdateLocalRealPart( i, jLoc, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,STAR,MD,Int>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[* ,MD]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (j + this->rowAlignment_) % r;
    const Int ownerCol = (j + this->rowAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
        this->UpdateLocalImagPart( i, jLoc, u );
    }
}

template class DistMatrix<int,STAR,MD,int>;
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,MC,  MR,  int>& A );
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,MC,  STAR,int>& A );
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,MD,  STAR,int>& A );
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,MR,  MC,  int>& A );
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,MR,  STAR,int>& A );
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,STAR,MC,  int>& A );
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,STAR,MR,  int>& A );
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,STAR,STAR,int>& A );
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,STAR,VC,  int>& A );
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,STAR,VR,  int>& A );
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,VC,  STAR,int>& A );
template DistMatrix<int,STAR,MD,int>::DistMatrix( const DistMatrix<int,VR,  STAR,int>& A );

#ifndef DISABLE_FLOAT
template class DistMatrix<float,STAR,MD,int>;
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,MC,  MR,  int>& A );
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,MC,  STAR,int>& A );
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,MD,  STAR,int>& A );
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,MR,  MC,  int>& A );
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,MR,  STAR,int>& A );
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,STAR,MC,  int>& A );
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,STAR,MR,  int>& A );
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,STAR,STAR,int>& A );
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,STAR,VC,  int>& A );
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,STAR,VR,  int>& A );
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,VC,  STAR,int>& A );
template DistMatrix<float,STAR,MD,int>::DistMatrix( const DistMatrix<float,VR,  STAR,int>& A );
#endif // ifndef DISABLE_FLOAT

template class DistMatrix<double,STAR,MD,int>;
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,MC,  MR,  int>& A );
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,MC,  STAR,int>& A );
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,MD,  STAR,int>& A );
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,MR,  MC,  int>& A );
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,MR,  STAR,int>& A );
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,STAR,MC,  int>& A );
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,STAR,MR,  int>& A );
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,STAR,STAR,int>& A );
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,STAR,VC,  int>& A );
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,STAR,VR,  int>& A );
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,VC,  STAR,int>& A );
template DistMatrix<double,STAR,MD,int>::DistMatrix( const DistMatrix<double,VR,  STAR,int>& A );

#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
template class DistMatrix<Complex<float>,STAR,MD,int>;
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,MC,  MR,  int>& A );
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,MC,  STAR,int>& A );
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,MD,  STAR,int>& A );
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,MR,  MC,  int>& A );
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,MR,  STAR,int>& A );
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,MC,  int>& A );
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,MR,  int>& A );
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,STAR,int>& A );
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,VC,  int>& A );
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,VR,  int>& A );
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,VC,  STAR,int>& A );
template DistMatrix<Complex<float>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<float>,VR,  STAR,int>& A );
#endif // ifndef DISABLE_FLOAT
template class DistMatrix<Complex<double>,STAR,MD,int>;
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,MC,  MR,  int>& A );
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,MC,  STAR,int>& A );
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,MD,  STAR,int>& A );
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,MR,  MC,  int>& A );
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,MR,  STAR,int>& A );
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,MC,  int>& A );
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,MR,  int>& A );
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,STAR,int>& A );
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,VC,  int>& A );
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,VR,  int>& A );
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,VC,  STAR,int>& A );
template DistMatrix<Complex<double>,STAR,MD,int>::DistMatrix( const DistMatrix<Complex<double>,VR,  STAR,int>& A );
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
