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
 * DistMatrix_Dist<CIRC,CIRC,Int>
 */

template<typename Int>
void
DistMatrix_Dist<CIRC,CIRC,Int>::AssertValidRoot( Int root )
{
    if( root < 0 || root >= this->Grid().Size() )
        throw std::logic_error("Invalid root");
}

template<typename Int>
DistMatrix_Dist<CIRC,CIRC,Int>::DistMatrix_Dist( const elem::Grid& g, Int root )
: DistMatrix_Base<Int>( g ), root_( root )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::DistMatrix");
    AssertValidRoot( root );
#endif
    this->SetShifts();
}

template<typename Int>
elem::Distribution
DistMatrix_Dist<CIRC,CIRC,Int>::ColDist() const { return CIRC; }

template<typename Int>
elem::Distribution
DistMatrix_Dist<CIRC,CIRC,Int>::RowDist() const { return CIRC; }

template<typename Int>
Int
DistMatrix_Dist<CIRC,CIRC,Int>::ColRank() const
{ return 0; }

template<typename Int>
Int
DistMatrix_Dist<CIRC,CIRC,Int>::RowRank() const
{ return 0; }

template<typename Int>
Int
DistMatrix_Dist<CIRC,CIRC,Int>::ColStride() const
{ return 1; }

template<typename Int>
Int
DistMatrix_Dist<CIRC,CIRC,Int>::RowStride() const
{ return 1; }

template<typename Int>
bool
DistMatrix_Dist<CIRC,CIRC,Int>::Participating() const
{ return this->Grid().Rank() == this->root_; }

template<typename Int>
Int
DistMatrix_Dist<CIRC,CIRC,Int>::Root() const
{ return root_; }

template<typename Int>
void
DistMatrix_Dist<CIRC,CIRC,Int>::SetRoot( Int root )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::SetRoot");
    AssertValidRoot( root );
#endif
    if( root != this->root_ )
        this->Empty();
    root_ = root;
}

template <typename Int>
void 
DistMatrix_Dist<CIRC,CIRC,Int>::Attach
( Int height, Int width, void* buffer, Int ldim, const elem::Grid& g, Int root )
{
    this->root_ = root;
    DistMatrix_Base<Int>::Attach( height, width, 0, 0, buffer, ldim, g );
}

template <typename Int>
void 
DistMatrix_Dist<CIRC,CIRC,Int>::LockedAttach
( Int height, Int width, const void* buffer, Int ldim, const elem::Grid& g, Int root )
{
    this->root_ = root;
    DistMatrix_Base<Int>::LockedAttach( height, width, 0, 0, buffer, ldim, g );
}

template<typename Int>
bool
DistMatrix_Dist<CIRC,CIRC,Int>::Index( Int i, Int j, Int& iLocal, Int& jLocal, int& mpiSrc, mpi::Comm& mpiDst ) const
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ]::Index");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    mpiSrc = g.VCToViewingMap(root_);
    mpiDst = g.ViewingComm();
    if ( g.Rank() != root_ ) return false;
    iLocal = i; 
    jLocal = j;
    return true;
}

template<typename Int>
void
DistMatrix_Dist<CIRC,CIRC,Int>::MakeConsistent()
{
#ifndef RELEASE
    CallStackEntry cse("[o ,o ]::MakeConsistent");
#endif
    const elem::Grid& g = this->Grid();
    const int root = g.VCToViewingMap(0);
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
            throw std::logic_error("Inconsistent ViewType");
        if( this->height_ != newHeight )
            throw std::logic_error("Inconsistent height");
        if( this->width_ != newWidth )
            throw std::logic_error("Inconsistent width");
        if( this->root_ != newRoot )
            throw std::logic_error("Inconsistent root");
    }
#endif
}

/*
 * DistMatrix<T,CIRC,CIRC,Int>
 */

template<typename T,typename Int>
DistMatrix<T,CIRC,CIRC,Int>::DistMatrix( const elem::Grid& g, Int root )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<CIRC,CIRC,Int>(g,root), DistMatrix_Type<T,Int>(g)
{}

template<typename T,typename Int>
DistMatrix<T,CIRC,CIRC,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g, Int root )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<CIRC,CIRC,Int>(g,root), DistMatrix_Type<T,Int>(g)
{ this->ResizeTo( height, width ); }

template<typename T,typename Int>
DistMatrix<T,CIRC,CIRC,Int>::DistMatrix
( Int height, Int width, Int ldim, const elem::Grid& g, Int root )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<CIRC,CIRC,Int>(g,root), DistMatrix_Type<T,Int>(g)
{ this->ResizeTo( height, width, ldim ); }

template<typename T,typename Int>
DistMatrix<T,CIRC,CIRC,Int>::DistMatrix
( Int height, Int width, const T* buffer, Int ldim, const elem::Grid& g, Int root )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<CIRC,CIRC,Int>(g,root), DistMatrix_Type<T,Int>(g)
{ this->LockedAttach( height, width, buffer, ldim, g, root ); }

template<typename T,typename Int>
DistMatrix<T,CIRC,CIRC,Int>::DistMatrix
( Int height, Int width, T* buffer, Int ldim, const elem::Grid& g, Int root )
: DistMatrix_Base<Int>(g), DistMatrix_Dist<CIRC,CIRC,Int>(g,root), DistMatrix_Type<T,Int>(g)
{ this->Attach( height, width, buffer, ldim, g, root ); }

template<typename T,typename Int>
DistMatrix<T,CIRC,CIRC,Int>::DistMatrix( const DistMatrix<T,CIRC,CIRC,Int>& A )
: DistMatrix_Base<Int>(A.Grid()), DistMatrix_Dist<CIRC,CIRC,Int>(A.Grid(),A.Root()), DistMatrix_Type<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[o ,o ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [o ,o ] with itself");
}

template<typename T,typename Int>
template<Distribution U,Distribution V>
DistMatrix<T,CIRC,CIRC,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: DistMatrix_Base<Int>(A.Grid()), DistMatrix_Dist<CIRC,CIRC,Int>(A.Grid(),0), DistMatrix_Type<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[o ,o ]::DistMatrix");
#endif
    if( CIRC != U || CIRC != V || reinterpret_cast<const DistMatrix_Base<Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [CIRC,CIRC] with itself");
}

template<typename T,typename Int>
DistMatrix<T,CIRC,CIRC,Int>::~DistMatrix()
{ }

//
// Utility functions, e.g., operator=
//

template<typename T,typename Int>
void
DistMatrix<T,CIRC,CIRC,Int>::CopyFromRoot( const Matrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[o ,o ]::CopyFromRoot");
#endif
    const Grid& grid = this->Grid();
    if( grid.VCRank() != this->Root() )
        throw std::logic_error("Called CopyFromRoot from non-root");

    int dims[2];
    dims[0] = A.Height();
    dims[1] = A.Width();
    mpi::Broadcast( dims, 2, this->Root(), grid.VCComm() );

    this->ResizeTo( dims[0], dims[1] );
    this->matrix_ = A;
}

template<typename T,typename Int>
void
DistMatrix<T,CIRC,CIRC,Int>::CopyFromNonRoot()
{
#ifndef RELEASE
    CallStackEntry cse("[o ,o ]::CopyFromNonRoot");
#endif
    const Grid& grid = this->Grid();
    if( grid.VCRank() == this->Root() )
        throw std::logic_error("Called CopyFromNonRoot from root");

    int dims[2];
    mpi::Broadcast( dims, 2, this->Root(), grid.VCComm() );

    this->ResizeTo( dims[0], dims[1] );
}

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
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
    const int root = this->Root();
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*mLocalA], &ABuf[jLoc*ALDim], mLocalA );

    // Communicate
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
        const Int rowAlignA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int l=0; l<rowStride; ++l )
        {
            const Int rowShift = Shift_( l, rowAlignA, rowStride );
            const Int nLocal = Length_( n, rowShift, rowStride );
            for( Int k=0; k<colStride; ++k )
            {
                const T* data = &recvBuf[(k+l*colStride)*pkgSize];
                const Int colShift = Shift_( k, colAlignA, colStride );
                const Int mLocal = Length_( m, colShift, colStride );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
        this->ResizeTo( m, n );

    const int root = this->Root();
    const elem::Grid& g = this->Grid();
    const int owningRow = root % g.Height();
    const int owningCol = root / g.Height();
    if( !g.InGrid() || g.Col() != owningCol )
        return *this;

    const int colStride = A.ColStride();
    const Int mLocalA = A.LocalHeight();
    const Int mLocalMax = MaxLength(m,colStride);

    const Int pkgSize = mpi::Pad( mLocalMax*n );
    T *sendBuf, *recvBuf;
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<n; ++j )
        MemCopy( &sendBuf[j*mLocalA], &ABuf[j*ALDim], mLocalA );

    // Communicate
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, owningRow, g.ColComm() );

    if( g.Row() == owningRow )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<colStride; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int colShift = Shift_( k, colAlignA, colStride );
            const Int mLocal = Length_( m, colShift, colStride );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
        this->ResizeTo( m, n );

    const int root = this->Root();
    const elem::Grid& g = this->Grid();
    const int owningRow = root % g.Height();
    const int owningCol = root / g.Height();
    if( !g.InGrid() || g.Row() != owningRow )
        return *this;

    const int rowStride = A.RowStride();
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*m], &ABuf[jLoc*ALDim], m );

    // Communicate
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, owningCol, g.RowComm() );

    if( g.Col() == owningCol )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int rowAlignA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<rowStride; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int rowShift = Shift_( k, rowAlignA, rowStride );
            const Int nLocal = Length_( n, rowShift, rowStride );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                MemCopy
                ( &buffer[(rowShift+jLoc*rowStride)*ldim], &data[jLoc*m], m );
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
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
    const int root = this->Root();
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int j=0; j<n; ++j )
            MemCopy( &sendBuf[j*mLocalA], &ABuf[j*ALDim], mLocalA );
    }

    // Communicate
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<p; ++k )
        {
            if( g.DiagPath( k ) == ownerPath )
            {
                const T* data = &recvBuf[k*pkgSize];
                const Int pathRank = g.DiagPathRank( k );
                const Int colShift = Shift_( pathRank, ownerPathRank, lcm );
                const Int mLocal = Length_( m, colShift, lcm );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
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
    const int root = this->Root();
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
            MemCopy( &sendBuf[jLoc*m], &ABuf[jLoc*ALDim], m );
    }

    // Communicate
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<p; ++k )
        {
            if( g.DiagPath( k ) == ownerPath )
            {
                const T* data = &recvBuf[k*pkgSize];
                const Int pathRank = g.DiagPathRank( k );
                const Int rowShift = Shift_( pathRank, ownerPathRank, lcm );
                const Int nLocal = Length_( n, rowShift, lcm );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
                for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                    MemCopy
                    ( &buffer[(rowShift+jLoc*lcm)*ldim], &data[jLoc*m], m );
            }
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
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
    const int root = this->Root();
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*mLocalA], &ABuf[jLoc*ALDim], mLocalA );

    // Communicate
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
        const Int rowAlignA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int l=0; l<rowStride; ++l )
        {
            const Int rowShift = Shift_( l, rowAlignA, rowStride );
            const Int nLocal = Length_( n, rowShift, rowStride );
            for( Int k=0; k<colStride; ++k )
            {
                const T* data = &recvBuf[(l+k*rowStride)*pkgSize];
                const Int colShift = Shift_( k, colAlignA, colStride );
                const Int mLocal = Length_( m, colShift, colStride );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
        this->ResizeTo( m, n );

    const int root = this->Root();
    const elem::Grid& g = this->Grid();
    const int owningRow = root % g.Height();
    const int owningCol = root / g.Height();
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<n; ++j )
        MemCopy( &sendBuf[j*mLocalA], &ABuf[j*ALDim], mLocalA );

    // Communicate
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, owningCol, g.RowComm() );

    if( g.Col() == owningCol )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<colStride; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int colShift = Shift_( k, colAlignA, colStride );
            const Int mLocal = Length_( m, colShift, colStride );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
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

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [* ,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
        this->ResizeTo( m, n );

    const int root = this->Root();
    const elem::Grid& g = this->Grid();
    const int owningRow = root % g.Height();
    const int owningCol = root / g.Height();
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*m], &ABuf[jLoc*ALDim], m );

    // Communicate
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, owningRow, g.ColComm() );

    if( g.Row() == owningRow )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int rowAlignA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<rowStride; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int rowShift = Shift_( k, rowAlignA, rowStride );
            const Int nLocal = Length_( n, rowShift, rowStride );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                MemCopy
                ( &buffer[(rowShift+jLoc*rowStride)*ldim], &data[jLoc*m], m );
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
        this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int mLocalA = A.LocalHeight();
    const Int mLocalMax = MaxLength(m,p);

    const Int pkgSize = mpi::Pad( mLocalMax*n );
    const int root = this->Root();
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<n; ++j )
        MemCopy( &sendBuf[j*mLocalA], &ABuf[j*ALDim], mLocalA );

    // Communicate
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<p; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int colShift = Shift_( k, colAlignA, p );
            const Int mLocal = Length_( m, colShift, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
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

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [o ,o ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int nLocalA = A.LocalWidth();
    const Int nLocalMax = MaxLength(n,p);

    const Int pkgSize = mpi::Pad( m*nLocalMax );
    const int root = this->Root();
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*m], &ABuf[jLoc*ALDim], m );

    // Communicate
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, root, g.VCComm() );

    if( g.VCRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int rowAlignA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<p; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int rowShift = Shift_( k, rowAlignA, p );
            const Int nLocal = Length_( n, rowShift, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                MemCopy( &buffer[(rowShift+jLoc*p)*ldim], &data[jLoc*m], m );
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
        this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int mLocalA = A.LocalHeight();
    const Int mLocalMax = MaxLength(m,p);

    const Int pkgSize = mpi::Pad( mLocalMax*n );
    const int root = this->Root();
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<n; ++j )
        MemCopy( &sendBuf[j*mLocalA], &ABuf[j*ALDim], mLocalA );

    // Communicate
    const int rootRow = root % g.Height();
    const int rootCol = root / g.Height();
    const int rootVR = rootCol + rootRow*g.Width();
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, rootVR, g.VRComm() );

    if( g.VRRank() == rootVR )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<p; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int colShift = Shift_( k, colAlignA, p );
            const Int mLocal = Length_( m, colShift, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
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

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
        this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int nLocalA = A.LocalWidth();
    const Int nLocalMax = MaxLength(n,p);

    const Int pkgSize = mpi::Pad( m*nLocalMax );
    const int root = this->Root();
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
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*m], &ABuf[jLoc*ALDim], m );

    // Communicate
    const int rootRow = root % g.Height();
    const int rootCol = root / g.Height();
    const int rootVR = rootCol + rootRow*g.Width();
    mpi::Gather
    ( sendBuf, pkgSize,
      recvBuf, pkgSize, rootVR, g.VRComm() );

    if( g.VRRank() == rootVR )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int rowAlignA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<p; ++k )
        {
            const T* data = &recvBuf[k*pkgSize];
            const Int rowShift = Shift_( k, rowAlignA, p );
            const Int nLocal = Length_( n, rowShift, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                MemCopy( &buffer[(rowShift+jLoc*p)*ldim], &data[jLoc*m], m );
        }
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    if( A.Grid().VCRank() == this->Root() )
        this->matrix_ = A.LockedMatrix();

    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,CIRC,CIRC,Int>&
DistMatrix<T,CIRC,CIRC,Int>::operator=( const DistMatrix<T,CIRC,CIRC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[o ,o ] = [o ,o ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( !this->Viewing() )
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
            mpi::Send( sendBuf, m*n, this->Root(), 0, g.VCComm() );
        }
        else if( g.VCRank() == this->Root() )
        {
            // Recv
            T* recvBuf = this->auxMemory_.Require( m*n );
            mpi::Recv( recvBuf, m*n, A.Root(), 0, g.VCComm() );
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

template class DistMatrix_Dist<CIRC,CIRC,int>;

#define PROTO(T) \
  template class DistMatrix<T,CIRC,CIRC,int>
#define COPY(T,CD,RD) \
  template DistMatrix<T,CIRC,CIRC,int>::DistMatrix( \
    const DistMatrix<T,CD,RD,int>& A )
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
