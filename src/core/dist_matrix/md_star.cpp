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
DistMatrix<T,MD,STAR,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g), diagPath_(0)
{ this->SetShifts(); } 

template<typename T,typename Int>
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g), diagPath_(0)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T,typename Int>
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignmentVC, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g), diagPath_(g.DiagPath(colAlignmentVC))
{ 
    this->Align( g.DiagPathRank(colAlignmentVC), 0 );
    this->ResizeTo( height, width );
}

template<typename T,typename Int>
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignmentVC, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g), diagPath_(g.DiagPath(colAlignmentVC))
{ 
    this->Align( g.DiagPathRank(colAlignmentVC), 0 );
    this->ResizeTo( height, width, ldim );
}

template<typename T,typename Int>
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignmentVC, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g), diagPath_(g.DiagPath(colAlignmentVC))
{ this->LockedAttach(height,width,colAlignmentVC,buffer,ldim,g); }

template<typename T,typename Int>
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignmentVC, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g), diagPath_(g.DiagPath(colAlignmentVC))
{ this->Attach(height,width,colAlignmentVC,buffer,ldim,g); }

template<typename T,typename Int>
DistMatrix<T,MD,STAR,Int>::DistMatrix( const DistMatrix<T,MD,STAR,Int>& A )
: AbstractDistMatrix<T,Int>(A.Grid()), diagPath_(0)
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MD,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [MD,* ] with itself");
}

template<typename T,typename Int>
template<Distribution U,Distribution V>
DistMatrix<T,MD,STAR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(A.Grid()), diagPath_(0)
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MD,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( MD != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,MD,STAR,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [MD,* ] with itself");
}

template<typename T,typename Int>
DistMatrix<T,MD,STAR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
elem::DistData<Int>
DistMatrix<T,MD,STAR,Int>::DistData() const
{
    elem::DistData<Int> data;
    data.colDist = MD;
    data.rowDist = STAR;
    data.colAlignment = this->colAlignment_;
    data.rowAlignment = 0;
    data.root = 0;
    data.diagPath = this->diagPath_;
    data.grid = this->grid_;
    return data;
}

template<typename T,typename Int>
Int
DistMatrix<T,MD,STAR,Int>::ColStride() const
{ return this->grid_->LCM(); }

template<typename T,typename Int>
Int
DistMatrix<T,MD,STAR,Int>::RowStride() const
{ return 1; }

template<typename T,typename Int>
Int
DistMatrix<T,MD,STAR,Int>::ColRank() const
{ return this->grid_->DiagPathRank(); }

template<typename T,typename Int>
Int
DistMatrix<T,MD,STAR,Int>::RowRank() const
{ return 0; }

template<typename T,typename Int>
bool
DistMatrix<T,MD,STAR,Int>::Participating() const
{
    const Grid& g = this->Grid();
    return ( g.InGrid() && g.DiagPath()==this->diagPath_ );
}

template<typename T,typename Int>
Int
DistMatrix<T,MD,STAR,Int>::DiagPath() const
{ return this->diagPath_; }

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::AlignWith( const elem::DistData<Int>& data )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::AlignWith");
    this->AssertFreeColAlignment();
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( data.colDist == MD && data.rowDist == STAR )
    {
        this->colAlignment_ = data.colAlignment;
        this->diagPath_ = data.diagPath;
    }
    else if( data.colDist == STAR && data.rowDist == MD )
    {
        this->colAlignment_ = data.rowAlignment;
        this->diagPath_ = data.diagPath;
    }
#ifndef RELEASE
    else throw std::logic_error("Invalid alignment");
#endif
    this->constrainedColAlignment_ = true;
    this->SetShifts();
}

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::AlignWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::AlignColsWith( const elem::DistData<Int>& data )
{ this->AlignWith( data ); }

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::AlignColsWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T,typename Int>
bool
DistMatrix<T,MD,STAR,Int>::AlignedWithDiagonal
( const elem::DistData<Int>& data, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::AlignedWithDiagonal");
#endif
    const Grid& grid = this->Grid();
    if( grid != *data.grid )
        return false;

    bool aligned;
    const Int r = grid.Height();
    const Int c = grid.Width();
    const Int firstDiagRow = 0;
    const Int firstDiagCol = this->diagPath_;
    const Int diagRow = (firstDiagRow+this->ColAlignment()) % r;
    const Int diagCol = (firstDiagCol+this->ColAlignment()) % c;
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
                    this->colAlignment_==data.colAlignment );
    }
    else if( data.colDist == STAR && data.rowDist == MD )
    {
        aligned = ( this->diagPath_==data.diagPath && 
                    this->colAlignment_==data.rowAlignment );
    }
    else aligned = false;
    return aligned;
}

template<typename T,typename Int>
bool
DistMatrix<T,MD,STAR,Int>::AlignedWithDiagonal
( const AbstractDistMatrix<T,Int>& A, Int offset ) const
{ return this->AlignedWithDiagonal( A.DistData(), offset ); }

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::AlignWithDiagonal
( const elem::DistData<Int>& data, Int offset )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::AlignWithDiagonal");
    this->AssertFreeColAlignment();
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
        this->colAlignment_ = grid.DiagPathRank(owner);
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
        this->colAlignment_ = grid.DiagPathRank(owner);
    }
    else if( data.colDist == MD && data.rowDist == STAR )
    {
        this->diagPath_ = data.diagPath;
        this->colAlignment_ = data.colAlignment;
    }
    else if( data.colDist == STAR && data.rowDist == MD )
    {
        this->diagPath_ = data.diagPath;
        this->colAlignment_ = data.rowAlignment;
    }
#ifndef RELEASE
    else throw std::logic_error("Nonsensical AlignWithDiagonal");
#endif
    this->constrainedColAlignment_ = true;
    this->SetShifts();
}

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::AlignWithDiagonal
( const AbstractDistMatrix<T,Int>& A, Int offset )
{ this->AlignWithDiagonal( A.DistData(), offset ); }

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::Attach
( Int height, Int width, Int colAlignmentVC,
  T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::Attach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->diagPath_ = grid.DiagPath(colAlignmentVC);
    this->colAlignment_ = grid.DiagPathRank(colAlignmentVC);
    this->viewing_ = true;
    this->SetColShift();
    if( this->Participating() )
    {
        const Int localHeight = Length(height,this->colShift_,grid.LCM());
        this->matrix_.Attach( localHeight, width, buffer, ldim );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::LockedAttach
( Int height, Int width, Int colAlignmentVC,
  const T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->diagPath_ = grid.DiagPath(colAlignmentVC);
    this->colAlignment_ = grid.DiagPathRank(colAlignmentVC);
    this->viewing_ = true;
    this->locked_ = true;
    this->SetColShift();
    if( this->Participating() )
    {
        const Int localHeight = Length(height,this->colShift_,grid.LCM());
        this->matrix_.LockedAttach( localHeight, width, buffer, ldim );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo
        ( Length(height,this->ColShift(),this->Grid().LCM()), width );
}

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo
        ( Length(height,this->ColShift(),this->Grid().LCM()), width, ldim );
}

template<typename T,typename Int>
T
DistMatrix<T,MD,STAR,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    T u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.LCM();
        u = this->GetLocal(iLoc,j);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
    return u;
}

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.LCM();
        this->SetLocal(iLoc,j,u);
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.LCM();
        this->UpdateLocal(iLoc,j,u);
    }
}

//
// Utility functions, e.g., operator=
//

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MD,* ] = [MC,MR] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MD,* ] = [MC,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MD,* ] = [* ,MR] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->diagPath_ = A.diagPath_;
            this->colAlignment_ = A.colAlignment_;
            if( this->Participating() )
                this->colShift_ = A.ColShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->diagPath_ == A.diagPath_ && 
        this->colAlignment_ == A.colAlignment_ )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned [MD,* ] <- [MD,* ]." << std::endl;
#endif
        throw std::logic_error
        ("Unaligned [MD,* ] = [MD,* ] not yet implemented");
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MD,* ] = [* ,MD] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() ); 
#endif
    throw std::logic_error("[MD,* ] = [MR,MC] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MD,* ] = [MR,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [* ,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MD,* ] = [* ,MC] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MD,* ] = [VC,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [* ,VC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MD,* ] = [* ,VC] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MD,* ] = [VR,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MD,* ] = [* ,VR] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    if( this->Participating() )
    {
        const Int lcm = this->grid_->LCM();
        const Int colShift = this->ColShift();

        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();

        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for 
#endif
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisBuffer[j*thisLDim];
            const T* sourceCol = &ABuf[colShift+j*ALDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*lcm];
        }
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,CIRC,CIRC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ] = [o ,o ]");
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

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    typedef BASE(T) R;

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.LCM();
        this->SetLocalRealPart( iLoc, j, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    typedef BASE(T) R;
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.LCM();
        this->SetLocalImagPart( iLoc, j, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    typedef BASE(T) R;

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.LCM();
        this->UpdateLocalRealPart( iLoc, j, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MD,STAR,Int>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MD,* ]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    typedef BASE(T) R;
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.LCM();
        this->UpdateLocalImagPart( iLoc, j, u );
    }
}

template class DistMatrix<int,MD,STAR,int>;
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,CIRC,CIRC,int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,MC,  MR,  int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,MC,  STAR,int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,MR,  MC,  int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,MR,  STAR,int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,STAR,MC,  int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,STAR,MD,  int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,STAR,MR,  int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,STAR,STAR,int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,STAR,VC,  int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,STAR,VR,  int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,VC,  STAR,int>& A );
template DistMatrix<int,MD,STAR,int>::DistMatrix( const DistMatrix<int,VR,  STAR,int>& A );

#ifndef DISABLE_FLOAT
template class DistMatrix<float,MD,STAR,int>;
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,CIRC,CIRC,int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,MC,  MR,  int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,MC,  STAR,int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,MR,  MC,  int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,MR,  STAR,int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,STAR,MC,  int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,STAR,MD,  int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,STAR,MR,  int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,STAR,STAR,int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,STAR,VC,  int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,STAR,VR,  int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,VC,  STAR,int>& A );
template DistMatrix<float,MD,STAR,int>::DistMatrix( const DistMatrix<float,VR,  STAR,int>& A );
#endif // ifndef DISABLE_FLOAT

template class DistMatrix<double,MD,STAR,int>;
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,CIRC,CIRC,int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,MC,  MR,  int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,MC,  STAR,int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,MR,  MC,  int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,MR,  STAR,int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,STAR,MC,  int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,STAR,MD,  int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,STAR,MR,  int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,STAR,STAR,int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,STAR,VC,  int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,STAR,VR,  int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,VC,  STAR,int>& A );
template DistMatrix<double,MD,STAR,int>::DistMatrix( const DistMatrix<double,VR,  STAR,int>& A );

#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
template class DistMatrix<Complex<float>,MD,STAR,int>;
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,CIRC,CIRC,int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,MC,  MR,  int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,MC,  STAR,int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,MR,  MC,  int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,MR,  STAR,int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,MC,  int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,MD,  int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,MR,  int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,STAR,int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,VC,  int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,STAR,VR,  int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,VC,  STAR,int>& A );
template DistMatrix<Complex<float>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<float>,VR,  STAR,int>& A );
#endif // ifndef DISABLE_FLOAT
template class DistMatrix<Complex<double>,MD,STAR,int>;
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,CIRC,CIRC,int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,MC,  MR,  int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,MC,  STAR,int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,MR,  MC,  int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,MR,  STAR,int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,MC,  int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,MD,  int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,MR,  int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,STAR,int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,VC,  int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,STAR,VR,  int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,VC,  STAR,int>& A );
template DistMatrix<Complex<double>,MD,STAR,int>::DistMatrix( const DistMatrix<Complex<double>,VR,  STAR,int>& A );
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
