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
using DM = DistMatrix<T,MD,STAR>;

template<typename T>
DM<T>::DistMatrix( const elem::Grid& g )
: ADM<T>(g)
{ this->SetShifts(); } 

template<typename T>
DM<T>::DistMatrix( Int height, Int width, const elem::Grid& g )
: ADM<T>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, Int root, const elem::Grid& g )
: ADM<T>(g)
{ 
    this->root_ = root;
    this->Align( colAlign, 0 );
    this->ResizeTo( height, width );
}

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, Int root, Int ldim, const elem::Grid& g )
: ADM<T>(g)
{ 
    this->root_ = root;
    this->Align( colAlign, 0 );
    this->ResizeTo( height, width, ldim );
}

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, Int root, const T* buffer, Int ldim,
  const elem::Grid& g )
: ADM<T>(g)
{ this->LockedAttach(height,width,colAlign,root,buffer,ldim,g); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, Int root, T* buffer, Int ldim,
  const elem::Grid& g )
: ADM<T>(g)
{ this->Attach(height,width,colAlign,root,buffer,ldim,g); }

template<typename T>
DM<T>::DistMatrix( const DM<T>& A )
: ADM<T>(A.Grid())
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
DM<T>::DistMatrix( const DistMatrix<T,U,V>& A )
: ADM<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[MD,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( MD != U || STAR != V || 
        reinterpret_cast<const DM<T>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [MD,* ] with itself");
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
{ return this->grid_->MDComm(); }

template<typename T>
mpi::Comm
DM<T>::CrossComm() const
{ return this->grid_->MDPerpComm(); }

template<typename T>
mpi::Comm
DM<T>::RedundantComm() const
{ return mpi::COMM_SELF; }

template<typename T>
mpi::Comm
DM<T>::ColComm() const
{ return this->grid_->MDComm(); }

template<typename T>
mpi::Comm
DM<T>::RowComm() const
{ return mpi::COMM_SELF; }

template<typename T>
Int
DM<T>::ColStride() const
{ return this->grid_->LCM(); }
    
template<typename T>
Int 
DM<T>::RowStride() const
{ return 1; }

template<typename T>
void
DM<T>::ShallowSwap( DM<T>& A )
{ ADM<T>::ShallowSwap( A ); }

template<typename T>
void
DM<T>::AlignWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::AlignWith");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( data.colDist == MD && data.rowDist == STAR )
    {
        this->colAlign_ = data.colAlign;
        this->root_ = data.root;
    }
    else if( data.colDist == STAR && data.rowDist == MD )
    {
        this->colAlign_ = data.rowAlign;
        this->root_ = data.root;
    }
#ifndef RELEASE
    else LogicError("Invalid alignment");
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
    CallStackEntry cse("[MD,* ]::AlignedWithDiagonal");
#endif
    const Grid& grid = this->Grid();
    if( grid != *data.grid )
        return false;

    bool aligned;
    const Int r = grid.Height();
    const Int c = grid.Width();
    const Int firstDiagRow = 0;
    const Int firstDiagCol = this->root_;
    const Int diagRow = (firstDiagRow+this->ColAlign()) % r;
    const Int diagCol = (firstDiagCol+this->ColAlign()) % c;
    if( data.colDist == MC && data.rowDist == MR )
    {
        if( offset >= 0 )
        {
            const Int ownerRow = data.colAlign;
            const Int ownerCol = (data.rowAlign + offset) % c;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
        else
        {
            const Int ownerRow = (data.colAlign-offset) % r;
            const Int ownerCol = data.rowAlign;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
    }
    else if( data.colDist == MR && data.rowDist == MC )
    {
        if( offset >= 0 )
        {
            const Int ownerCol = data.colAlign;
            const Int ownerRow = (data.rowAlign + offset) % r;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
        else
        {
            const Int ownerCol = (data.colAlign-offset) % c;
            const Int ownerRow = data.rowAlign;
            aligned = ( ownerRow==diagRow && ownerCol==diagCol );
        }
    }
    else if( data.colDist == MD && data.rowDist == STAR )
    {
        aligned = ( this->root_==data.root && 
                    this->colAlign_==data.colAlign );
    }
    else if( data.colDist == STAR && data.rowDist == MD )
    {
        aligned = ( this->root_==data.root && 
                    this->colAlign_==data.rowAlign );
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
            const Int ownerRow = data.colAlign;
            const Int ownerCol = (data.rowAlign + offset) % c;
            owner = ownerRow + r*ownerCol;
        }
        else
        {
            const Int ownerRow = (data.colAlign-offset) % r;
            const Int ownerCol = data.rowAlign;
            owner = ownerRow + r*ownerCol;
        }
        this->root_ = grid.DiagPath(owner);
        this->colAlign_ = grid.DiagPathRank(owner);
    }
    else if( data.colDist == MR && data.rowDist == MC )
    {
        Int owner;
        if( offset >= 0 )
        {
            const Int ownerCol = data.colAlign;
            const Int ownerRow = (data.rowAlign + offset) % r;
            owner = ownerRow + r*ownerCol;
        }
        else
        {
            const Int ownerCol = (data.colAlign-offset) % c;
            const Int ownerRow = data.rowAlign;
            owner = ownerRow + r*ownerCol;
        }
        this->root_ = grid.DiagPath(owner);
        this->colAlign_ = grid.DiagPathRank(owner);
    }
    else if( data.colDist == MD && data.rowDist == STAR )
    {
        this->root_ = data.root;
        this->colAlign_ = data.colAlign;
    }
    else if( data.colDist == STAR && data.rowDist == MD )
    {
        this->root_ = data.root;
        this->colAlign_ = data.rowAlign;
    }
#ifndef RELEASE
    else LogicError("Nonsensical AlignWithDiagonal");
#endif
    this->colConstrained_ = true;
    this->SetShifts();
}

template<typename T>
void
DM<T>::AlignWithDiagonal( const ADM<T>& A, Int offset )
{ this->AlignWithDiagonal( A.DistData(), offset ); }

template<typename T>
void
DM<T>::Attach
( Int height, Int width, Int colAlign, Int root,
  T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::Attach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->root_ = root;
    this->colAlign_ = colAlign;
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
DM<T>::LockedAttach
( Int height, Int width, Int colAlign, Int root,
  const T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->root_ = root;
    this->colAlign_ = colAlign;
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
DM<T>::Attach
( Matrix<T>& A, Int colAlign, Int root, const elem::Grid& g )
{
    this->Attach
    ( A.Height(), A.Width(), colAlign, root, A.Buffer(), A.LDim(), g );
}

template<typename T>
void
DM<T>::LockedAttach
( const Matrix<T>& A, Int colAlign, Int root, const elem::Grid& g )
{
    this->LockedAttach
    ( A.Height(), A.Width(), colAlign, root, A.LockedBuffer(), A.LDim(), g );
}

//
// Utility functions, e.g., operator=
//

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,MR>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,STAR>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MR>& A )
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
const DM<T>&
DM<T>::operator=( const DM<T>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MD,* ] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    if( !this->Viewing() && !this->ColConstrained() )
    {
        this->root_ = A.root_;
        this->colAlign_ = A.colAlign_;
        if( this->Participating() )
            this->colShift_ = A.ColShift();
    }
    this->ResizeTo( A.Height(), A.Width() );

    if( this->root_ == A.root_ && this->colAlign_ == A.colAlign_ )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MD>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,MC>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,STAR>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MC>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VC,STAR>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VC>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VR,STAR>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VR>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,STAR>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
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
#define COPY(T,U,V) \
  template DistMatrix<T,MD,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
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
