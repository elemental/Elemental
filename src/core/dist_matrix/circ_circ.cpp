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
DistMatrix<T,CIRC,CIRC>::DistMatrix( const elem::Grid& g, Int root )
: AbstractDistMatrix<T>(g)
{ 
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::DistMatrix");
    if( root < 0 || root >= this->grid_->Size() )
        LogicError("Invalid root");
#endif
    this->root_ = root; 
}

template<typename T>
DistMatrix<T,CIRC,CIRC>::DistMatrix
( Int height, Int width, const elem::Grid& g, Int root )
: AbstractDistMatrix<T>(g)
{ 
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::DistMatrix");
    if( root < 0 || root >= this->grid_->Size() )
        LogicError("Invalid root");
#endif
    this->root_ = root;
    this->ResizeTo( height, width );
}

template<typename T>
DistMatrix<T,CIRC,CIRC>::DistMatrix
( Int height, Int width, Int ldim, const elem::Grid& g, Int root )
: AbstractDistMatrix<T>(g)
{ 
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::DistMatrix");
    if( root < 0 || root >= this->grid_->Size() )
        LogicError("Invalid root");
#endif
    this->root_ = root;
    this->ResizeTo( height, width, ldim );
}

template<typename T>
DistMatrix<T,CIRC,CIRC>::DistMatrix
( Int height, Int width, const T* buffer, Int ldim, const elem::Grid& g, 
  Int root )
: AbstractDistMatrix<T>(g)
{ 
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::DistMatrix");
    if( root < 0 || root >= this->grid_->Size() )
        LogicError("Invalid root");
#endif
    this->root_ = root;
    this->LockedAttach( height, width, buffer, ldim, g, root );
}

template<typename T>
DistMatrix<T,CIRC,CIRC>::DistMatrix
( Int height, Int width, T* buffer, Int ldim, const elem::Grid& g, Int root )
: AbstractDistMatrix<T>(g)
{ 
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::DistMatrix");
    if( root < 0 || root >= this->grid_->Size() )
        LogicError("Invalid root");
#endif
    this->root_ = root;
    this->Attach( height, width, buffer, ldim, g, root );
}

template<typename T>
DistMatrix<T,CIRC,CIRC>::DistMatrix( const DistMatrix<T,CIRC,CIRC>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[o ,o ]::DistMatrix");
#endif
    this->root_ = A.Root();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [o ,o ] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DistMatrix<T,CIRC,CIRC>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[o ,o ]::DistMatrix");
#endif
    this->root_ = 0;
    if( CIRC != U || CIRC != V || 
        reinterpret_cast<const DistMatrix<T,CIRC,CIRC>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [o ,o ] with itself");
}

template<typename T>
DistMatrix<T,CIRC,CIRC>::DistMatrix( DistMatrix<T,CIRC,CIRC>&& A )
: AbstractDistMatrix<T>(std::move(A)), root_(A.root_)
{ }

template<typename T>
DistMatrix<T,CIRC,CIRC>& 
DistMatrix<T,CIRC,CIRC>::operator=( DistMatrix<T,CIRC,CIRC>&& A )
{
    AbstractDistMatrix<T>::operator=( std::move(A) );
    root_ = A.root_;
    return *this;
}

template<typename T>
DistMatrix<T,CIRC,CIRC>::~DistMatrix()
{ }

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::Swap( DistMatrix<T,CIRC,CIRC>& A )
{
    AbstractDistMatrix<T>::Swap( A );        
    std::swap( root_, A.root_ );
}

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::SetRoot( Int root )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::SetRoot");
    if( root < 0 || root >= this->grid_->Size() )
        LogicError("Invalid root");
#endif
    if( root != this->root_ )
        this->Empty();
    this->root_ = root;
}

template<typename T>
Int
DistMatrix<T,CIRC,CIRC>::Root() const
{ return this->root_; }

template<typename T>
elem::DistData
DistMatrix<T,CIRC,CIRC>::DistData() const
{
    elem::DistData data;
    data.colDist = CIRC;
    data.rowDist = CIRC;
    data.colAlignment = 0;
    data.rowAlignment = 0;
    data.root = 0;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T>
Int
DistMatrix<T,CIRC,CIRC>::ColStride() const
{ return 1; }

template<typename T>
Int
DistMatrix<T,CIRC,CIRC>::RowStride() const
{ return 1; }

template<typename T>
Int
DistMatrix<T,CIRC,CIRC>::ColRank() const
{ return 0; }

template<typename T>
Int
DistMatrix<T,CIRC,CIRC>::RowRank() const
{ return 0; }

template<typename T>
bool
DistMatrix<T,CIRC,CIRC>::Participating() const
{ return ( this->Grid().Rank() == this->root_ ); }

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::Attach
( Int height, Int width, 
  T* buffer, Int ldim, const elem::Grid& grid, Int root )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::Attach");
#endif
    this->grid_ = &grid;
    this->SetRoot( root );
    this->height_ = height;
    this->width_ = width;
    this->viewType_ = VIEW;
    if( this->Participating() )
        this->matrix_.Attach_( height, width, buffer, ldim );
}

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::LockedAttach
( Int height, Int width, 
  const T* buffer, Int ldim, const elem::Grid& grid, Int root )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::LockedAttach");
#endif
    this->grid_ = &grid;
    this->SetRoot( root );
    this->height_ = height;
    this->width_ = width;
    this->viewType_ = LOCKED_VIEW;
    if( this->Participating() )
        this->matrix_.LockedAttach_( height, width, buffer, ldim );
}

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_( height, width );
}

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_( height, width, ldim );
}

template<typename T>
T
DistMatrix<T,CIRC,CIRC>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::Get");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    T u;
    if( this->Participating() )
        u = this->GetLocal( i, j );
    mpi::Broadcast( u, g.VCToViewingMap(this->root_), g.ViewingComm() );
    return u;
}

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::Set");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->SetLocal(i,j,u);
}

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::Update");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->UpdateLocal(i,j,u);
}

//
// Utility functions, e.g., operator=
//

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::MakeConsistent()
{
#ifndef RELEASE
    CallStackEntry cse("[o ,o ]::MakeConsistent");
#endif
    const elem::Grid& g = this->Grid();
    const Int root = g.VCToViewingMap(0);
    Int message[4];
    if( g.ViewingRank() == root )
    {
        message[0] = this->viewType_;
        message[1] = this->height_;
        message[2] = this->width_;
        message[3] = this->root_;
    }
    mpi::Broadcast( message, 4, root, g.ViewingComm() );
    const ViewType newViewType = static_cast<ViewType>(message[0]);
    const Int newHeight = message[1];
    const Int newWidth = message[2];
    const Int newRoot = message[3];
    if( !this->Participating() )
    {
        this->viewType_ = newViewType;
        this->height_ = newHeight;
        this->width_ = newWidth;
        this->root_ = newRoot;
        this->constrainedColAlignment_ = false;
        this->constrainedColAlignment_ = false;
        this->colAlignment_ = 0;
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
        if( this->root_ != newRoot )
            LogicError("Inconsistent root");
    }
#endif
}

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::CopyFromRoot( const Matrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[o ,o ]::CopyFromRoot");
#endif
    const Grid& grid = this->Grid();
    if( grid.VCRank() != this->Root() )
        LogicError("Called CopyFromRoot from non-root");

    Int dims[2];
    dims[0] = A.Height();
    dims[1] = A.Width();
    mpi::Broadcast( dims, 2, this->Root(), grid.VCComm() );

    this->ResizeTo( dims[0], dims[1] );
    this->matrix_ = A;
}

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::CopyFromNonRoot()
{
#ifndef RELEASE
    CallStackEntry cse("[o ,o ]::CopyFromNonRoot");
#endif
    const Grid& grid = this->Grid();
    if( grid.VCRank() == this->Root() )
        LogicError("Called CopyFromNonRoot from root");

    Int dims[2];
    mpi::Broadcast( dims, 2, this->Root(), grid.VCComm() );

    this->ResizeTo( dims[0], dims[1] );
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int mLocalA = A.LocalHeight();
    const Int nLocalA = A.LocalWidth();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    const Int mLocalMax = MaxLength(m,colStride);
    const Int nLocalMax = MaxLength(n,rowStride);

    const Int pkgSize = mpi::Pad( mLocalMax*nLocalMax );
    const Int p = g.Size();
    const Int root = this->Root();
    T *sendBuf, *recvBuf;
    if( g.VCRank() == root )
    {
        T* buffer = this->auxMemory_.Require( (p+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
        recvBuf = 0;
    }

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*mLocalA], &ABuf[jLoc*ALDim], mLocalA );

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
        const Int rowAlignA = A.RowAlignment();
        OUTER_PARALLEL_FOR
        for( Int l=0; l<rowStride; ++l )
        {
            const Int rowShift = Shift_( l, rowAlignA, rowStride );
            const Int nLocal = Length_( n, rowShift, rowStride );
            for( Int k=0; k<colStride; ++k )
            {
                const T* data = &recvBuf[(k+l*colStride)*pkgSize];
                const Int colShift = Shift_( k, colAlignA, colStride );
                const Int mLocal = Length_( m, colShift, colStride );
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                {
                    T* destCol = 
                      &buffer[colShift+(rowShift+jLoc*rowStride)*ldim];
                    const T* sourceCol = &data[jLoc*mLocal];
                    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                        destCol[iLoc*colStride] = sourceCol[iLoc];
                }
            }
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );

    const Int root = this->Root();
    const elem::Grid& g = this->Grid();
    const Int owningRow = root % g.Height();
    const Int owningCol = root / g.Height();
    if( !g.InGrid() || g.Col() != owningCol )
        return *this;

    const Int colStride = A.ColStride();
    const Int mLocalA = A.LocalHeight();
    const Int mLocalMax = MaxLength(m,colStride);

    const Int pkgSize = mpi::Pad( mLocalMax*n );
    T *sendBuf, *recvBuf=0; // some compilers (falsely) warn otherwise
    if( g.Row() == owningRow )
    {
        T* buffer = this->auxMemory_.Require( (colStride+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
    }

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    PARALLEL_FOR
    for( Int j=0; j<n; ++j )
        MemCopy( &sendBuf[j*mLocalA], &ABuf[j*ALDim], mLocalA );

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, owningRow, g.ColComm() );

    if( g.Row() == owningRow )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<colStride; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int colShift = Shift_( k, colAlignA, colStride );
            const Int mLocal = Length_( m, colShift, colStride );
            INNER_PARALLEL_FOR
            for( Int j=0; j<n; ++j )
            {
                T* destCol = &buffer[colShift+j*ldim];
                const T* sourceCol = &data[j*mLocal];
                for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                    destCol[iLoc*colStride] = sourceCol[iLoc];
            }
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );

    const Int root = this->Root();
    const elem::Grid& g = this->Grid();
    const Int owningRow = root % g.Height();
    const Int owningCol = root / g.Height();
    if( !g.InGrid() || g.Row() != owningRow )
        return *this;

    const Int rowStride = A.RowStride();
    const Int nLocalA = A.LocalWidth();
    const Int nLocalMax = MaxLength(n,rowStride);

    const Int pkgSize = mpi::Pad( m*nLocalMax );
    T *sendBuf, *recvBuf;
    if( g.Col() == owningCol )
    {
        T* buffer = this->auxMemory_.Require( (rowStride+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
        recvBuf = 0;
    }

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*m], &ABuf[jLoc*ALDim], m );

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, owningCol, g.RowComm() );

    if( g.Col() == owningCol )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int rowAlignA = A.RowAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStride; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int rowShift = Shift_( k, rowAlignA, rowStride );
            const Int nLocal = Length_( n, rowShift, rowStride );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                MemCopy
                ( &buffer[(rowShift+jLoc*rowStride)*ldim], &data[jLoc*m], m );
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int lcm = g.LCM();
    const Int ownerPath = A.diagPath_;
    const Int ownerPathRank = A.colAlignment_;

    const Int mLocalA = A.LocalHeight();
    const Int mLocalMax = MaxLength( m, lcm );
    const Int pkgSize = mpi::Pad( mLocalMax*n );

    // Since a MD communicator has not been implemented, we will take
    // the suboptimal route of 'rounding up' everyone's contribution over 
    // the VC communicator.
    const Int root = this->Root();
    T *sendBuf, *recvBuf;
    if( g.VCRank() == root )
    {
        T* buffer = this->auxMemory_.Require( (p+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
        recvBuf = 0;
    }

    if( A.Participating() )
    {
        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        PARALLEL_FOR
        for( Int j=0; j<n; ++j )
            MemCopy( &sendBuf[j*mLocalA], &ABuf[j*ALDim], mLocalA );
    }

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<p; ++k )
        {
            if( g.DiagPath( k ) == ownerPath )
            {
                const T* data = &recvBuf[k*pkgSize];
                const Int pathRank = g.DiagPathRank( k );
                const Int colShift = Shift_( pathRank, ownerPathRank, lcm );
                const Int mLocal = Length_( m, colShift, lcm );
                INNER_PARALLEL_FOR
                for( Int j=0; j<n; ++j )
                {
                    T* destCol = &buffer[colShift+j*ldim];
                    const T* sourceCol = &data[j*mLocal];
                    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                        destCol[iLoc*lcm] = sourceCol[iLoc];
                }
            }
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int lcm = g.LCM();
    const Int ownerPath = A.diagPath_;
    const Int ownerPathRank = A.rowAlignment_;

    const Int nLocalA = A.LocalWidth();
    const Int nLocalMax = MaxLength( n, lcm );
    const Int pkgSize = mpi::Pad( m*nLocalMax );

    // Since a MD communicator has not been implemented, we will take
    // the suboptimal route of 'rounding up' everyone's contribution over 
    // the VC communicator.
    const Int root = this->Root();
    T *sendBuf, *recvBuf;
    if( g.VCRank() == root )
    {
        T* buffer = this->auxMemory_.Require( (p+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
        recvBuf = 0;
    }

    if( A.Participating() )
    {
        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
            MemCopy( &sendBuf[jLoc*m], &ABuf[jLoc*ALDim], m );
    }

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<p; ++k )
        {
            if( g.DiagPath( k ) == ownerPath )
            {
                const T* data = &recvBuf[k*pkgSize];
                const Int pathRank = g.DiagPathRank( k );
                const Int rowShift = Shift_( pathRank, ownerPathRank, lcm );
                const Int nLocal = Length_( n, rowShift, lcm );
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                    MemCopy
                    ( &buffer[(rowShift+jLoc*lcm)*ldim], &data[jLoc*m], m );
            }
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int mLocalA = A.LocalHeight();
    const Int nLocalA = A.LocalWidth();
    const Int rowStride = A.RowStride();
    const Int colStride = A.ColStride();
    const Int mLocalMax = MaxLength( m, colStride );
    const Int nLocalMax = MaxLength( n, rowStride );

    const Int pkgSize = mpi::Pad( mLocalMax*nLocalMax );
    const Int p = g.Size();
    const Int root = this->Root();
    T *sendBuf, *recvBuf;
    if( g.VCRank() == root )
    {
        T* buffer = this->auxMemory_.Require( (p+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
        recvBuf = 0;
    }

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*mLocalA], &ABuf[jLoc*ALDim], mLocalA );

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
        const Int rowAlignA = A.RowAlignment();
        OUTER_PARALLEL_FOR
        for( Int l=0; l<rowStride; ++l )
        {
            const Int rowShift = Shift_( l, rowAlignA, rowStride );
            const Int nLocal = Length_( n, rowShift, rowStride );
            for( Int k=0; k<colStride; ++k )
            {
                const T* data = &recvBuf[(l+k*rowStride)*pkgSize];
                const Int colShift = Shift_( k, colAlignA, colStride );
                const Int mLocal = Length_( m, colShift, colStride );
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                {
                    T* destCol = 
                      &buffer[colShift+(rowShift+jLoc*rowStride)*ldim];
                    const T* sourceCol = &data[jLoc*mLocal];
                    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                        destCol[iLoc*colStride] = sourceCol[iLoc];
                }
            }
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );

    const Int root = this->Root();
    const elem::Grid& g = this->Grid();
    const Int owningRow = root % g.Height();
    const Int owningCol = root / g.Height();
    if( !g.InGrid() || g.Row() != owningRow )
        return *this;

    const Int colStride = A.ColStride();
    const Int mLocalA = A.LocalHeight();
    const Int mLocalMax = MaxLength(m,colStride);

    const Int pkgSize = mpi::Pad( mLocalMax*n );
    T *sendBuf, *recvBuf;
    if( g.Col() == owningCol )
    {
        T* buffer = this->auxMemory_.Require( (colStride+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
        recvBuf = 0; 
    }

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    PARALLEL_FOR
    for( Int j=0; j<n; ++j )
        MemCopy( &sendBuf[j*mLocalA], &ABuf[j*ALDim], mLocalA );

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, owningCol, g.RowComm() );

    if( g.Col() == owningCol )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<colStride; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int colShift = Shift_( k, colAlignA, colStride );
            const Int mLocal = Length_( m, colShift, colStride );
            INNER_PARALLEL_FOR
            for( Int j=0; j<n; ++j )
            {
                T* destCol = &buffer[colShift+j*ldim];
                const T* sourceCol = &data[j*mLocal];
                for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                    destCol[iLoc*colStride] = sourceCol[iLoc];
            }
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [* ,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );

    const Int root = this->Root();
    const elem::Grid& g = this->Grid();
    const Int owningRow = root % g.Height();
    const Int owningCol = root / g.Height();
    if( !g.InGrid() || g.Col() != owningCol )
        return *this;

    const Int rowStride = A.RowStride();
    const Int nLocalA = A.LocalWidth();
    const Int nLocalMax = MaxLength(n,rowStride);

    const Int pkgSize = mpi::Pad( m*nLocalMax );
    T *sendBuf, *recvBuf;
    if( g.Row() == owningRow )
    {
        T* buffer = this->auxMemory_.Require( (rowStride+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
        recvBuf = 0;
    }

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*m], &ABuf[jLoc*ALDim], m );

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, owningRow, g.ColComm() );

    if( g.Row() == owningRow )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int rowAlignA = A.RowAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStride; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int rowShift = Shift_( k, rowAlignA, rowStride );
            const Int nLocal = Length_( n, rowShift, rowStride );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                MemCopy
                ( &buffer[(rowShift+jLoc*rowStride)*ldim], &data[jLoc*m], m );
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,VC,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int mLocalA = A.LocalHeight();
    const Int mLocalMax = MaxLength(m,p);

    const Int pkgSize = mpi::Pad( mLocalMax*n );
    const Int root = this->Root();
    T *sendBuf, *recvBuf;
    if( g.VCRank() == root )
    {
        T* buffer = this->auxMemory_.Require( (p+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
        recvBuf = 0;
    }

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    PARALLEL_FOR
    for( Int j=0; j<n; ++j )
        MemCopy( &sendBuf[j*mLocalA], &ABuf[j*ALDim], mLocalA );

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<p; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int colShift = Shift_( k, colAlignA, p );
            const Int mLocal = Length_( m, colShift, p );
            INNER_PARALLEL_FOR
            for( Int j=0; j<n; ++j )
            {
                T* destCol = &buffer[colShift+j*ldim];
                const T* sourceCol = &data[j*mLocal];
                for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                    destCol[iLoc*p] = sourceCol[iLoc];
            }
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,STAR,VC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [o ,o ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( A.Height(), A.Width() );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int nLocalA = A.LocalWidth();
    const Int nLocalMax = MaxLength(n,p);

    const Int pkgSize = mpi::Pad( m*nLocalMax );
    const Int root = this->Root();
    T *sendBuf, *recvBuf;
    if( g.VCRank() == root )
    {
        T* buffer = this->auxMemory_.Require( (p+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
        recvBuf = 0;
    }

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*m], &ABuf[jLoc*ALDim], m );

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int rowAlignA = A.RowAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<p; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int rowShift = Shift_( k, rowAlignA, p );
            const Int nLocal = Length_( n, rowShift, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                MemCopy( &buffer[(rowShift+jLoc*p)*ldim], &data[jLoc*m], m );
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,VR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int mLocalA = A.LocalHeight();
    const Int mLocalMax = MaxLength(m,p);

    const Int pkgSize = mpi::Pad( mLocalMax*n );
    const Int root = this->Root();
    T *sendBuf, *recvBuf;
    if( g.VCRank() == root )
    {
        T* buffer = this->auxMemory_.Require( (p+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
        recvBuf = 0;
    }

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    PARALLEL_FOR
    for( Int j=0; j<n; ++j )
        MemCopy( &sendBuf[j*mLocalA], &ABuf[j*ALDim], mLocalA );

    // Communicate
    const Int rootRow = root % g.Height();
    const Int rootCol = root / g.Height();
    const Int rootVR = rootCol + rootRow*g.Width();
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, rootVR, g.VRComm() );

    if( g.VRRank() == rootVR )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<p; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int colShift = Shift_( k, colAlignA, p );
            const Int mLocal = Length_( m, colShift, p );
            INNER_PARALLEL_FOR
            for( Int j=0; j<n; ++j )
            {
                T* destCol = &buffer[colShift+j*ldim];
                const T* sourceCol = &data[j*mLocal];
                for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                    destCol[iLoc*p] = sourceCol[iLoc];
            }
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,STAR,VR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int nLocalA = A.LocalWidth();
    const Int nLocalMax = MaxLength(n,p);

    const Int pkgSize = mpi::Pad( m*nLocalMax );
    const Int root = this->Root();
    T *sendBuf, *recvBuf;
    if( g.VCRank() == root )
    {
        T* buffer = this->auxMemory_.Require( (p+1)*pkgSize );
        sendBuf = &buffer[0];
        recvBuf = &buffer[pkgSize];
    }
    else
    {
        sendBuf = this->auxMemory_.Require( pkgSize );
        recvBuf = 0;
    }

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*m], &ABuf[jLoc*ALDim], m );

    // Communicate
    const Int rootRow = root % g.Height();
    const Int rootCol = root / g.Height();
    const Int rootVR = rootCol + rootRow*g.Width();
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, rootVR, g.VRComm() );

    if( g.VRRank() == rootVR )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int rowAlignA = A.RowAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<p; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int rowShift = Shift_( k, rowAlignA, p );
            const Int nLocal = Length_( n, rowShift, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                MemCopy( &buffer[(rowShift+jLoc*p)*ldim], &data[jLoc*m], m );
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [* ,* ]");
#endif
    this->ResizeTo( A.Height(), A.Width() );
    if( A.Grid().VCRank() == this->Root() )
        this->matrix_ = A.LockedMatrix();
    return *this;
}

template<typename T>
const DistMatrix<T,CIRC,CIRC>&
DistMatrix<T,CIRC,CIRC>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [o ,o ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );

    const Grid& g = A.Grid();
    if( this->Root() == A.Root() )
    {
        if( g.VCRank() == A.Root() )
            this->matrix_ = A.matrix_;
    }
    else
    {
        if( g.VCRank() == A.Root() )
        {
            T* sendBuf = this->auxMemory_.Require( m*n );
            // Pack
            const Int ALDim = A.LDim();
            const T* ABuf = A.LockedBuffer();
            for( Int j=0; j<n; ++j )
                for( Int i=0; i<m; ++i )
                    sendBuf[i+j*m] = ABuf[i+j*ALDim];
            // Send
            mpi::Send( sendBuf, m*n, this->Root(), g.VCComm() );
        }
        else if( g.VCRank() == this->Root() )
        {
            // Recv
            T* recvBuf = this->auxMemory_.Require( m*n );
            mpi::Recv( recvBuf, m*n, A.Root(), g.VCComm() );
            // Unpack
            const Int ldim = this->LDim();
            T* buffer = this->Buffer();
            for( Int j=0; j<n; ++j )
                for( Int i=0; i<m; ++i )
                    buffer[i+j*ldim] = recvBuf[i+j*m];
        }
        this->auxMemory_.Release();
    }

    return *this;
}

//
// Routines which explicitly work in the complex plane
//

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->SetLocalRealPart(i,j,u);
}

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    if( this->Participating() )
        this->SetLocalImagPart(i,j,u);
}

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->UpdateLocalRealPart(i,j,u);
}

template<typename T>
void
DistMatrix<T,CIRC,CIRC>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    if( this->Participating() )
        this->UpdateLocalImagPart(i,j,u);
}

#define PROTO(T) template class DistMatrix<T,CIRC,CIRC>
#define COPY(T,CD,RD) \
  template DistMatrix<T,CIRC,CIRC>::DistMatrix( const DistMatrix<T,CD,RD>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,MC,  MR); \
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
