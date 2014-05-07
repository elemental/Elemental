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
BDM&
BDM::operator=( const BlockDistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MC,MR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MC,STAR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,MR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MD,STAR]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,MD]"))
    LogicError("This routine is not yet written");
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MR,MC]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [MR,STAR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,MC]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [VC,STAR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,VC]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [VR,STAR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,VR]"))
    this->CollectFrom( A );
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BlockDistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [STAR,STAR]"))
    this->Resize( A.Height(), A.Width() );
    if( A.Grid().VCRank() == this->Root() )
        this->matrix_ = A.LockedMatrix();
    return *this;
}

template<typename T>
BDM&
BDM::operator=( const BDM& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC] = [CIRC,CIRC]"))
    A.Translate( *this );
    return *this;
}

template<typename T>
void
BDM::CopyFromRoot( const Matrix<T>& A, bool includingViewers )
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
BDM::CopyFromNonRoot( bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::CopyFromNonRoot"))
    if( this->CrossRank() == this->Root() )
        LogicError("Called CopyFromNonRoot from root");

    this->MakeSizeConsistent( includingViewers );
}

// Basic queries
// =============

template<typename T>
mpi::Comm BDM::DistComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::CrossComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm BDM::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm BDM::RowComm() const { return mpi::COMM_SELF; }

template<typename T>
Int BDM::ColStride() const { return 1; }
template<typename T>
Int BDM::RowStride() const { return 1; }
template<typename T>
Int BDM::DistSize() const { return 1; }
template<typename T>
Int BDM::CrossSize() const { return this->grid_->VCSize(); }
template<typename T>
Int BDM::RedundantSize() const { return 1; }

// Private section
// ###############

template<typename T>
template<Dist U,Dist V>
void
BDM::CollectFrom( const BlockDistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::CollectFrom"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mb = A.BlockHeight();
    const Int nb = A.BlockWidth();
    const Int colCut = A.ColCut();
    const Int rowCut = A.RowCut();
    this->Align( mb, nb, 0, 0, false );
    this->Resize( m, n );
    if( A.RedundantSize() != 1 )
        LogicError("This routine does not yet support non-trivial redundancy");
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
    const Int mLocalMax = MaxBlockedLength(m,mb,colCut,colStride);
    const Int nLocalMax = MaxBlockedLength(n,nb,rowCut,rowStride);
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
        const Int colAlignA = A.ColAlign();
        const Int rowAlignA = A.RowAlign();
        ELEM_OUTER_PARALLEL_FOR
        for( Int l=0; l<rowStride; ++l )
        {
            const Int rowShift = Shift_( l, rowAlignA, rowStride );
            const Int nLocal = 
                BlockedLength_( n, rowShift, nb, rowCut, rowStride );
            for( Int k=0; k<colStride; ++k )
            {
                const T* data = &recvBuf[(k+l*colStride)*pkgSize];
                const Int colShift = Shift_( k, colAlignA, colStride );
                const Int mLocal = 
                    BlockedLength_( m, colShift, mb, colCut, colStride );
                ELEM_INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                {
                    const Int jBefore = rowShift*nb - rowCut;
                    const Int jLocAdj = ( rowShift==0 ? jLoc+rowCut : jLoc );
                    const Int numFilledLocalBlocks = jLocAdj / nb;
                    const Int jMid = numFilledLocalBlocks*nb*rowStride;
                    const Int jPost = jLocAdj-numFilledLocalBlocks*nb;
                    const Int j = jBefore + jMid + jPost;
                    const T* sourceCol = &data[jLoc*mLocal];
                    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                    {
                        const Int iBefore = colShift*mb - colCut;
                        const Int iLocAdj = (colShift==0 ? iLoc+colCut : iLoc);
                        const Int numFilledLocalBlocks = iLocAdj / mb;
                        const Int iMid = numFilledLocalBlocks*mb*colStride;
                        const Int iPost = iLocAdj-numFilledLocalBlocks*mb;
                        const Int i = iBefore + iMid + iPost;
                        this->SetLocal(i,j,sourceCol[iLoc]);
                    }
                }
            }
        }
    }
    this->auxMemory_.Release();
}

template<typename T>
template<Dist U,Dist V>
void
BDM::Scatter( BlockDistMatrix<T,U,V>& A ) const
{
    DEBUG_ONLY(CallStackEntry cse("[CIRC,CIRC]::Scatter"))
    if( A.CrossSize() != 1 )
        LogicError("This routine does not yet support non-trivial cross-teams");
    LogicError("This routine is not yet written");
}

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class BlockDistMatrix<T,ColDist,RowDist>
#define SELF(T,U,V) \
  template BlockDistMatrix<T,ColDist,RowDist>::BlockDistMatrix \
  ( const BlockDistMatrix<T,U,V>& A );
#define OTHER(T,U,V) \
  template BlockDistMatrix<T,ColDist,RowDist>::BlockDistMatrix \
  ( const DistMatrix<T,U,V>& A ); \
  template BlockDistMatrix<T,ColDist,RowDist>& \
           BlockDistMatrix<T,ColDist,RowDist>::operator= \
           ( const DistMatrix<T,U,V>& A )
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
