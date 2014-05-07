/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

#define ColDist CIRC
#define RowDist CIRC

#include "./setup.hpp"

namespace elem {

// Public section
// ##############

// Assignment and reconfiguration
// ==============================

template<typename T>
DM&
DM::operator=( const DM& A )
{
    DEBUG_ONLY(CallStackEntry cse("DM[U,V] = DM[U,V]"))
    A.Translate( *this );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MC,MR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MC,STAR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,MR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[CIRC,CIRC] = [MD,STAR]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const Int m = A.Height();
    const Int n = A.Width();
    this->Resize( m, n );
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
        ELEM_PARALLEL_FOR
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
        ELEM_OUTER_PARALLEL_FOR
        for( Int k=0; k<p; ++k )
        {
            if( g.DiagPath( k ) == ownerPath )
            {
                const T* data = &recvBuf[k*pkgSize];
                const Int pathRank = g.DiagPathRank( k );
                const Int colShift = Shift_( pathRank, ownerPathRank, lcm );
                const Int mLocal = Length_( m, colShift, lcm );
                ELEM_INNER_PARALLEL_FOR
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
DM&
DM::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[CIRC,CIRC] = [STAR,MD]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const Int m = A.Height();
    const Int n = A.Width();
    this->Resize( m, n );
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
        ELEM_PARALLEL_FOR
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
        ELEM_OUTER_PARALLEL_FOR
        for( Int k=0; k<p; ++k )
        {
            if( g.DiagPath( k ) == ownerPath )
            {
                const T* data = &recvBuf[k*pkgSize];
                const Int pathRank = g.DiagPathRank( k );
                const Int rowShift = Shift_( pathRank, ownerPathRank, lcm );
                const Int nLocal = Length_( n, rowShift, lcm );
                ELEM_INNER_PARALLEL_FOR
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
DM&
DM::operator=( const DistMatrix<T,MR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MR,MC]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MR,STAR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,MC]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,VC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [VC,STAR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,VC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,VC]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,VR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [VR,STAR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,VR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,VR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,STAR]"))
    this->Resize( A.Height(), A.Width() );
    if( A.Grid().VCRank() == this->Root() )
        this->matrix_ = A.LockedMatrix();
    return *this;
}

template<typename T>
void
DM::CopyFromRoot( const Matrix<T>& A, bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::CopyFromRoot"))
    if( this->CrossRank() != this->Root() )
        LogicError("Called CopyFromRoot from non-root");

    this->Resize( A.Height(), A.Width() );
    this->MakeSizeConsistent( includingViewers );

    this->matrix_ = A;
}

template<typename T>
void
DM::CopyFromNonRoot( bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::CopyFromNonRoot"))
    if( this->CrossRank() == this->Root() )
        LogicError("Called CopyFromNonRoot from root");

    this->MakeSizeConsistent( includingViewers );
}

// Basic queries
// =============

template<typename T>
mpi::Comm DM::DistComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::CrossComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm DM::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RowComm() const { return mpi::COMM_SELF; }

template<typename T>
Int DM::ColStride() const { return 1; }
template<typename T>
Int DM::RowStride() const { return 1; }
template<typename T>
Int DM::DistSize() const { return 1; }
template<typename T>
Int DM::CrossSize() const { return this->grid_->VCSize(); }
template<typename T>
Int DM::RedundantSize() const { return 1; }

// Private section
// ###############

template<typename T>
template<Dist U,Dist V>
void
DM::CollectFrom( const DistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::CollectFrom"))
    const Int m = A.Height();
    const Int n = A.Width();
    this->Resize( m, n ); 
    if( A.CrossSize() != 1 )
        LogicError("This routine does not yet support non-trivial cross-teams");
    if( !A.Grid().InGrid() )
        return;

    const Int root = this->Root();
    // Translate the root into our DistComm (if possible)
    const Int target = mpi::Translate( this->CrossComm(), root, A.DistComm() );
    if( target == mpi::UNDEFINED )
        return;

    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    const Int mLocalA = A.LocalHeight();
    const Int nLocalA = A.LocalWidth();
    const Int mLocalMax = MaxLength(m,colStride);
    const Int nLocalMax = MaxLength(n,rowStride);
    const Int pkgSize = mpi::Pad( mLocalMax*nLocalMax );
    const Int numDist = A.DistSize();

    T *sendBuf, *recvBuf;
    if( this->CrossRank() == root )
    {
        T* buffer = this->auxMemory_.Require( (numDist+1)*pkgSize );
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
    ELEM_PARALLEL_FOR
    for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        MemCopy( &sendBuf[jLoc*mLocalA], &ABuf[jLoc*ALDim], mLocalA );

    // Communicate
    mpi::Gather( sendBuf, pkgSize, recvBuf, pkgSize, target, A.DistComm() );

    if( this->CrossRank() == root )
    {
        // Unpack
        T* buffer = this->Buffer();
        const Int ldim = this->LDim();
        const Int colAlignA = A.ColAlign();
        const Int rowAlignA = A.RowAlign();
        ELEM_OUTER_PARALLEL_FOR
        for( Int l=0; l<rowStride; ++l )
        {
            const Int rowShift = Shift_( l, rowAlignA, rowStride );
            const Int nLocal = Length_( n, rowShift, rowStride );
            for( Int k=0; k<colStride; ++k )
            {
                const T* data = &recvBuf[(k+l*colStride)*pkgSize];
                const Int colShift = Shift_( k, colAlignA, colStride );
                const Int mLocal = Length_( m, colShift, colStride );
                ELEM_INNER_PARALLEL_FOR
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
}

template<typename T>
template<Dist U,Dist V>
void
DM::Scatter( DistMatrix<T,U,V>& A ) const
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::Scatter"))
    if( A.CrossSize() != 1 )
        LogicError("This routine does not yet support non-trivial cross-teams");
    LogicError("This routine is not yet written");
}

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class DistMatrix<T,ColDist,RowDist>
#define SELF(T,U,V) \
  template DistMatrix<T,ColDist,RowDist>::DistMatrix \
  ( const DistMatrix<T,U,V>& A );
#define OTHER(T,U,V) \
  template DistMatrix<T,ColDist,RowDist>::DistMatrix \
  ( const BlockDistMatrix<T,U,V>& A ); \
  template DistMatrix<T,ColDist,RowDist>& \
           DistMatrix<T,ColDist,RowDist>::operator= \
           ( const BlockDistMatrix<T,U,V>& A )
#define BOTH(T,U,V) \
  SELF(T,U,V); \
  OTHER(T,U,V)
#define FULL(T) \
  PROTO(T); \
  OTHER(T,CIRC,CIRC); \
  BOTH( T,MC,  MR  ); \
  BOTH( T,MC,  STAR); \
  BOTH( T,MD,  STAR); \
  BOTH( T,MR,  MC  ); \
  BOTH( T,MR,  STAR); \
  BOTH( T,STAR,MC  ); \
  BOTH( T,STAR,MD  ); \
  BOTH( T,STAR,MR  ); \
  BOTH( T,STAR,STAR); \
  BOTH( T,STAR,VC  ); \
  BOTH( T,STAR,VR  ); \
  BOTH( T,VC,  STAR); \
  BOTH( T,VR,  STAR);

FULL(Int);
#ifndef ELEM_DISABLE_FLOAT
FULL(float);
#endif
FULL(double);

#ifndef ELEM_DISABLE_COMPLEX
#ifndef ELEM_DISABLE_FLOAT
FULL(Complex<float>);
#endif
FULL(Complex<double>);
#endif

} // namespace elem
