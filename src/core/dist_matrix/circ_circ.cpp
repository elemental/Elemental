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
using DM = DistMatrix<T,CIRC,CIRC>;

template<typename T>
DM<T>::DistMatrix( const elem::Grid& g, Int root )
: ADM<T>(g)
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ]::DistMatrix");
        if( root < 0 || root >= this->CrossSize() )
            LogicError("Invalid root");
    )
    this->root_ = root; 
}

template<typename T>
DM<T>::DistMatrix( Int height, Int width, const elem::Grid& g, Int root )
: ADM<T>(g)
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ]::DistMatrix");
        if( root < 0 || root >= this->CrossSize() )
            LogicError("Invalid root");
    )
    this->root_ = root;
    this->ResizeTo( height, width );
}

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int ldim, const elem::Grid& g, Int root )
: ADM<T>(g)
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ]::DistMatrix");
        if( root < 0 || root >= this->CrossSize() )
            LogicError("Invalid root");
    )
    this->root_ = root;
    this->ResizeTo( height, width, ldim );
}

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, const T* buffer, Int ldim, const elem::Grid& g, 
  Int root )
: ADM<T>(g)
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ]::DistMatrix");
        if( root < 0 || root >= this->CrossSize() )
            LogicError("Invalid root");
    )
    this->root_ = root;
    this->LockedAttach( height, width, root, buffer, ldim, g );
}

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, T* buffer, Int ldim, const elem::Grid& g, Int root )
: ADM<T>(g)
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ]::DistMatrix");
        if( root < 0 || root >= this->CrossSize() )
            LogicError("Invalid root");
    )
    this->root_ = root;
    this->Attach( height, width, root, buffer, ldim, g );
}

template<typename T>
DM<T>::DistMatrix( const DM<T>& A )
: ADM<T>(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[o ,o ]::DistMatrix"))
    this->root_ = A.Root();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [o ,o ] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DM<T>::DistMatrix( const DistMatrix<T,U,V>& A )
: ADM<T>(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[o ,o ]::DistMatrix"))
    this->root_ = 0;
    if( CIRC != U || CIRC != V || 
        reinterpret_cast<const DM<T>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [o ,o ] with itself");
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
void
DM<T>::ShallowSwap( DM<T>& A )
{
    ADM<T>::ShallowSwap( A );        
}

template<typename T>
elem::DistData
DM<T>::DistData() const
{ return elem::DistData(*this); }

template<typename T>
mpi::Comm
DM<T>::DistComm() const
{ return mpi::COMM_SELF; }

template<typename T>
mpi::Comm
DM<T>::RedundantComm() const
{ return mpi::COMM_SELF; }

template<typename T>
mpi::Comm
DM<T>::CrossComm() const
{ return this->grid_->VCComm(); }

template<typename T>
mpi::Comm
DM<T>::ColComm() const
{ return mpi::COMM_SELF; }

template<typename T>
mpi::Comm
DM<T>::RowComm() const
{ return mpi::COMM_SELF; }

template<typename T>
Int
DM<T>::ColStride() const
{ return 1; }

template<typename T>
Int
DM<T>::RowStride() const
{ return 1; }

template<typename T>
void
DM<T>::Attach
( Int height, Int width, Int root,
  T* buffer, Int ldim, const elem::Grid& grid )
{
    DEBUG_ONLY(CallStackEntry cse("[o ,o ]::Attach"))
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
DM<T>::LockedAttach
( Int height, Int width, Int root,
  const T* buffer, Int ldim, const elem::Grid& grid )
{
    DEBUG_ONLY(CallStackEntry cse("[o ,o ]::LockedAttach"))
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
DM<T>::Attach( Matrix<T>& A, Int root, const elem::Grid& g )
{ this->Attach( A.Height(), A.Width(), root, A.Buffer(), A.LDim(), g ); }

template<typename T>
void
DM<T>::LockedAttach( const Matrix<T>& A, Int root, const elem::Grid& g )
{
    this->LockedAttach
    ( A.Height(), A.Width(), root, A.LockedBuffer(), A.LDim(), g );
}

//
// Utility functions, e.g., operator=
//

template<typename T>
void
DM<T>::CopyFromRoot( const Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[o ,o ]::CopyFromRoot"))
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
DM<T>::CopyFromNonRoot()
{
    DEBUG_ONLY(CallStackEntry cse("[o ,o ]::CopyFromNonRoot"))
    const Grid& grid = this->Grid();
    if( grid.VCRank() == this->Root() )
        LogicError("Called CopyFromNonRoot from root");

    Int dims[2];
    mpi::Broadcast( dims, 2, this->Root(), grid.VCComm() );

    this->ResizeTo( dims[0], dims[1] );
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [MC,MR]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
        const Int colAlignA = A.ColAlign();
        const Int rowAlignA = A.RowAlign();
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [MC,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
        const Int colAlignA = A.ColAlign();
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [* ,MR]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
        const Int rowAlignA = A.RowAlign();
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [MD,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int lcm = g.LCM();
    const Int ownerPath = A.root_;
    const Int ownerPathRank = A.colAlign_;

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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [* ,MD]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const Int m = A.Height();
    const Int n = A.Width();
    this->ResizeTo( m, n );
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
        return *this;

    const Int p = g.Size();
    const Int lcm = g.LCM();
    const Int ownerPath = A.root_;
    const Int ownerPathRank = A.rowAlign_;

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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,MC>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [MR,MC]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
        const Int colAlignA = A.ColAlign();
        const Int rowAlignA = A.RowAlign();
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [MR,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
        const Int colAlignA = A.ColAlign();
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MC>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [* ,MC]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
        const Int rowAlignA = A.RowAlign();
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VC,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [VC,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
        const Int colAlignA = A.ColAlign();
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VC>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [o ,o ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
        const Int rowAlignA = A.RowAlign();
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VR,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [VR,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
        const Int colAlignA = A.ColAlign();
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [* ,VR]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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
        const Int rowAlignA = A.RowAlign();
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[o ,o ] = [* ,* ]"))
    this->ResizeTo( A.Height(), A.Width() );
    if( A.Grid().VCRank() == this->Root() )
        this->matrix_ = A.LockedMatrix();
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DM<T>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[o ,o ] = [o ,o ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
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

#define PROTO(T) template class DistMatrix<T,CIRC,CIRC>
#define COPY(T,U,V) \
  template DistMatrix<T,CIRC,CIRC>::DistMatrix( const DistMatrix<T,U,V>& A )
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
