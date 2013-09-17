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
DistMatrix<T,MD,STAR>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T>(g), diagPath_(0)
{ this->SetShifts(); } 

template<typename T>
DistMatrix<T,MD,STAR>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T>(g), diagPath_(0)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T>
DistMatrix<T,MD,STAR>::DistMatrix
( Int height, Int width, Int colAlignmentVC, const elem::Grid& g )
: AbstractDistMatrix<T>(g), diagPath_(g.DiagPath(colAlignmentVC))
{ 
    this->Align( g.DiagPathRank(colAlignmentVC), 0 );
    this->ResizeTo( height, width );
}

template<typename T>
DistMatrix<T,MD,STAR>::DistMatrix
( Int height, Int width, Int colAlignmentVC, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g), diagPath_(g.DiagPath(colAlignmentVC))
{ 
    this->Align( g.DiagPathRank(colAlignmentVC), 0 );
    this->ResizeTo( height, width, ldim );
}

template<typename T>
DistMatrix<T,MD,STAR>::DistMatrix
( Int height, Int width, Int colAlignmentVC, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T>(g), diagPath_(g.DiagPath(colAlignmentVC))
{ this->LockedAttach(height,width,colAlignmentVC,buffer,ldim,g); }

template<typename T>
DistMatrix<T,MD,STAR>::DistMatrix
( Int height, Int width, Int colAlignmentVC, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T>(g), diagPath_(g.DiagPath(colAlignmentVC))
{ this->Attach(height,width,colAlignmentVC,buffer,ldim,g); }

template<typename T>
DistMatrix<T,MD,STAR>::DistMatrix( const DistMatrix<T,MD,STAR>& A )
: AbstractDistMatrix<T>(A.Grid()), diagPath_(0)
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[MD,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [MD,* ] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DistMatrix<T,MD,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(A.Grid()), diagPath_(0)
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[MD,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( MD != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,MD,STAR>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [MD,* ] with itself");
}

template<typename T>
DistMatrix<T,MD,STAR>::DistMatrix( DistMatrix<T,MD,STAR>&& A )
: AbstractDistMatrix<T>(std::move(A)), diagPath_(A.diagPath_)
{ }

template<typename T>
DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( DistMatrix<T,MD,STAR>&& A )
{
    AbstractDistMatrix<T>::operator=( std::move(A) );
    diagPath_ = A.diagPath_;
    return *this;
}

template<typename T>
DistMatrix<T,MD,STAR>::~DistMatrix()
{ }

template<typename T>
void
DistMatrix<T,MD,STAR>::Swap( DistMatrix<T,MD,STAR>& A )
{
    AbstractDistMatrix<T>::Swap( A );
    std::swap( diagPath_, A.diagPath_ );
}

template<typename T>
elem::DistData
DistMatrix<T,MD,STAR>::DistData() const
{
    elem::DistData data;
    data.colDist = MD;
    data.rowDist = STAR;
    data.colAlignment = this->colAlignment_;
    data.rowAlignment = 0;
    data.root = 0;
    data.diagPath = this->diagPath_;
    data.grid = this->grid_;
    return data;
}

template<typename T>
Int
DistMatrix<T,MD,STAR>::ColStride() const
{ return this->grid_->LCM(); }

template<typename T>
Int
DistMatrix<T,MD,STAR>::RowStride() const
{ return 1; }

template<typename T>
Int
DistMatrix<T,MD,STAR>::ColRank() const
{ return this->grid_->DiagPathRank(); }

template<typename T>
Int
DistMatrix<T,MD,STAR>::RowRank() const
{ return 0; }

template<typename T>
bool
DistMatrix<T,MD,STAR>::Participating() const
{
    const Grid& g = this->Grid();
    return ( g.InGrid() && g.DiagPath()==this->diagPath_ );
}

template<typename T>
Int
DistMatrix<T,MD,STAR>::DiagPath() const
{ return this->diagPath_; }

template<typename T>
void
DistMatrix<T,MD,STAR>::AlignWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::AlignWith");
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
    else LogicError("Invalid alignment");
#endif
    this->constrainedColAlignment_ = true;
    this->SetShifts();
}

template<typename T>
void
DistMatrix<T,MD,STAR>::AlignWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,MD,STAR>::AlignColsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

template<typename T>
void
DistMatrix<T,MD,STAR>::AlignColsWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
bool
DistMatrix<T,MD,STAR>::AlignedWithDiagonal
( const elem::DistData& data, Int offset ) const
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::AlignedWithDiagonal");
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

template<typename T>
bool
DistMatrix<T,MD,STAR>::AlignedWithDiagonal
( const AbstractDistMatrix<T>& A, Int offset ) const
{ return this->AlignedWithDiagonal( A.DistData(), offset ); }

template<typename T>
void
DistMatrix<T,MD,STAR>::AlignWithDiagonal
( const elem::DistData& data, Int offset )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::AlignWithDiagonal");
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
    else LogicError("Nonsensical AlignWithDiagonal");
#endif
    this->constrainedColAlignment_ = true;
    this->SetShifts();
}

template<typename T>
void
DistMatrix<T,MD,STAR>::AlignWithDiagonal
( const AbstractDistMatrix<T>& A, Int offset )
{ this->AlignWithDiagonal( A.DistData(), offset ); }

template<typename T>
void
DistMatrix<T,MD,STAR>::Attach
( Int height, Int width, Int colAlignmentVC,
  T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::Attach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->diagPath_ = grid.DiagPath(colAlignmentVC);
    this->colAlignment_ = grid.DiagPathRank(colAlignmentVC);
    this->viewType_ = VIEW;
    this->SetColShift();
    if( this->Participating() )
    {
        const Int localHeight = Length(height,this->colShift_,grid.LCM());
        this->matrix_.Attach_( localHeight, width, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,MD,STAR>::LockedAttach
( Int height, Int width, Int colAlignmentVC,
  const T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->diagPath_ = grid.DiagPath(colAlignmentVC);
    this->colAlignment_ = grid.DiagPathRank(colAlignmentVC);
    this->viewType_ = LOCKED_VIEW;
    this->SetColShift();
    if( this->Participating() )
    {
        const Int localHeight = Length(height,this->colShift_,grid.LCM());
        this->matrix_.LockedAttach_( localHeight, width, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,MD,STAR>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),this->Grid().LCM()), width );
}

template<typename T>
void
DistMatrix<T,MD,STAR>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),this->Grid().LCM()), width, ldim );
}

template<typename T>
T
DistMatrix<T,MD,STAR>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::Get");
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
    mpi::Broadcast( u, g.VCToViewingMap(ownerRank), g.ViewingComm() );
    return u;
}

template<typename T>
void
DistMatrix<T,MD,STAR>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::Set");
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

template<typename T>
void
DistMatrix<T,MD,STAR>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::SetRealPart");
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
        this->SetLocalRealPart( iLoc, j, u );
    }
}

template<typename T>
void
DistMatrix<T,MD,STAR>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
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

template<typename T>
void
DistMatrix<T,MD,STAR>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::Update");
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

template<typename T>
void
DistMatrix<T,MD,STAR>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::UpdateRealPart");
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
        this->UpdateLocalRealPart( iLoc, j, u );
    }
}

template<typename T>
void
DistMatrix<T,MD,STAR>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
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

//
// Utility functions, e.g., operator=
//

template<typename T>
void
DistMatrix<T,MD,STAR>::MakeConsistent()
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::MakeConsistent");
#endif
    const elem::Grid& g = this->Grid();
    const Int root = g.VCToViewingMap(0);
    Int message[6];
    if( g.ViewingRank() == root )
    {
        message[0] = this->viewType_;
        message[1] = this->height_;
        message[2] = this->width_;
        message[3] = this->constrainedColAlignment_;
        message[4] = this->colAlignment_;
        message[5] = this->diagPath_;
    }
    mpi::Broadcast( message, 6, root, g.ViewingComm() );
    const ViewType newViewType = static_cast<ViewType>(message[0]);
    const Int newHeight = message[1];
    const Int newWidth = message[2];
    const bool newConstrainedCol = message[3];
    const Int newColAlignment = message[4];
    const Int newDiagPath = message[5];
    if( !this->Participating() )
    {
        this->viewType_ = newViewType;
        this->height_ = newHeight;
        this->width_ = newWidth;
        this->constrainedColAlignment_ = newConstrainedCol;
        this->colAlignment_ = newColAlignment;
        this->diagPath_ = newDiagPath;
        this->constrainedRowAlignment_ = 0;
        this->rowAlignment_ = 0;
        this->colShift_ = 0;
        this->rowShift_ = 0;
    }
#ifndef RELEASE
    else
    {
        if( this->viewType_ != newViewType )
            LogicError("Inconsistent ViewType");
        if( this->height_ != newHeight )
            LogicError("Inconsistent height");
        if( this->width_ != newWidth )
            LogicError("Inconsistent width");
        if( this->constrainedColAlignment_ != newConstrainedCol ||
            this->colAlignment_ != newColAlignment )
            LogicError("Inconsistent column constraint");
        if( this->diagPath_ != newDiagPath )
            LogicError("Inconsistent diagonal path");
    }
#endif
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [MC,MR]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [MC,* ]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [* ,MR]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    if( !this->Viewing() && !this->ConstrainedColAlignment() )
    {
        this->diagPath_ = A.diagPath_;
        this->colAlignment_ = A.colAlignment_;
        if( this->Participating() )
            this->colShift_ = A.ColShift();
    }
    this->ResizeTo( A.Height(), A.Width() );

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
        // TODO: More efficient implementation?
        DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
        *this = A_STAR_STAR;
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,MD>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [* ,MD]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [MR,MC]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [MR,* ]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [* ,MC]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,VC,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [VC,* ]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,VC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [* ,VC]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,VR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [VR,* ]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR(A);
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,VR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [* ,VR]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int lcm = this->grid_->LCM();
    const Int colShift = this->ColShift();

    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();

    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    T* thisBuffer = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        T* destCol = &thisBuffer[j*thisLDim];
        const T* sourceCol = &ABuf[colShift+j*ALDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            destCol[iLoc] = sourceCol[iLoc*lcm];
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [o ,o ]");
#endif
    DistMatrix<T,MC,MR> A_MC_MR( A.Grid() );
    A_MC_MR.AlignWith( *this );
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

#define PROTO(T) template class DistMatrix<T,MD,STAR>
#define COPY(T,CD,RD) \
  template DistMatrix<T,MD,STAR>::DistMatrix( const DistMatrix<T,CD,RD>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  MR); \
  COPY(T,MC,  STAR); \
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
