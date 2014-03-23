/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

#define ColDist VC
#define RowDist STAR

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
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [MC,MR]"))
    this->PartialColAllToAllFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [MC,STAR]"))
    this->PartialColFilterFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,MR]"))
    DistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [MD,STAR]"))
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,MD]"))
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [MR,MC]"))
    DistMatrix<T,VR,STAR> A_VR_STAR( A );
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [MR,STAR]"))
    DistMatrix<T,VR,STAR> A_VR_STAR( A );
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,MC]"))
    std::unique_ptr<DistMatrix<T,MR,MC>> A_MR_MC( new DistMatrix<T,MR,MC>(A) );
    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(*A_MR_MC) );
    delete A_MR_MC.release(); // lowers memory highwater
    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,VC]"))
    std::unique_ptr<DistMatrix<T,MR,MC>> A_MR_MC( new DistMatrix<T,MR,MC>(A) );
    std::unique_ptr<DistMatrix<T,VR,STAR>> 
        A_VR_STAR( new DistMatrix<T,VR,STAR>(*A_MR_MC) );
    delete A_MR_MC.release(); // lowers memory highwater
    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[VC,STAR] = [VR,STAR]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const Grid& g = this->Grid();
    this->Resize( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;
    
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localHeightOfA = A.LocalHeight();

    const Int sendSize = localHeightOfA * width;
    const Int recvSize = localHeight * width;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int rankCM = g.VCRank();
    const Int rankRM = g.VRRank();

    const Int colShift = this->ColShift();
    const Int colShiftOfA = A.ColShift();

    // Compute which colmajor rank has the colShift equal to our colShiftOfA
    const Int sendRankCM = (rankCM+(p+colShiftOfA-colShift)) % p;

    // Compute which colmajor rank has the A colShift that we need
    const Int recvRankRM = (rankRM+(p+colShift-colShiftOfA)) % p;
    const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

    T* buffer = this->auxMemory_.Require( sendSize + recvSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[sendSize];

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    ELEM_PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        const T* ACol = &ABuf[j*ALDim];
        T* sendBufCol = &sendBuf[j*localHeightOfA];
        MemCopy( sendBufCol, ACol, localHeightOfA );
    }

    // Communicate
    mpi::SendRecv
    ( sendBuf, sendSize, sendRankCM,
      recvBuf, recvSize, recvRankCM, g.VCComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    ELEM_PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        const T* recvBufCol = &recvBuf[j*localHeight];
        T* thisCol = &thisBuf[j*thisLDim];
        MemCopy( thisCol, recvBufCol, localHeight );
    }

    this->auxMemory_.Release();
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,VR]"))
    DistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[VC,STAR] = [STAR,STAR]"))
    this->ColFilterFrom( A );
    return *this;
}

// NOTE: This is a small modification of [MC,MR] <- [CIRC,CIRC]
template<typename T>
DM&
DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[VC,STAR] = [CIRC,CIRC]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = g.Size();
    this->Resize( m, n );

    const Int colAlign = this->ColAlign();
    const Int mLocal = this->LocalHeight();
    const Int pkgSize = mpi::Pad(MaxLength(m,p)*n);
    const Int recvSize = pkgSize;
    const Int sendSize = p*pkgSize;
    T* recvBuf=0; // some compilers (falsely) warn otherwise
    if( A.Participating() )
    {
        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        recvBuf = &buffer[sendSize];

        // Pack the send buffer
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        for( Int s=0; s<p; ++s )
        {
            const Int sLocalHeight = Length( m, s, p );
            const Int q = (colAlign+s) % p;
            for( Int j=0; j<n; ++j )
            {
                for( Int iLoc=0; iLoc<sLocalHeight; ++iLoc )
                {
                    const Int i = s + iLoc*p;
                    sendBuf[q*pkgSize+iLoc+j*sLocalHeight] =
                        ABuf[i+j*ALDim];
                }
            }
        }

        // Scatter from the root
        mpi::Scatter
        ( sendBuf, pkgSize, recvBuf, pkgSize, A.Root(), g.VCComm() );
    }
    else if( this->Participating() )
    {
        recvBuf = this->auxMemory_.Require( recvSize );

        // Perform the receiving portion of the scatter from the non-root
        mpi::Scatter
        ( static_cast<T*>(0), pkgSize, 
          recvBuf,            pkgSize, A.Root(), g.VCComm() );
    }

    if( this->Participating() )
    {
        // Unpack
        const Int ldim = this->LDim();
        T* buffer = this->Buffer();
        for( Int j=0; j<n; ++j )
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                buffer[iLoc+j*ldim] = recvBuf[iLoc+j*mLocal];
        this->auxMemory_.Release();
    }

    return *this;
}

// Basic queries
// =============

template<typename T>
mpi::Comm DM::DistComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm DM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::ColComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm DM::RowComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::PartialColComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm DM::PartialUnionColComm() const { return this->grid_->MRComm(); }

template<typename T>
Int DM::ColStride() const { return this->grid_->VCSize(); }
template<typename T>
Int DM::RowStride() const { return 1; }
template<typename T>
Int DM::PartialColStride() const { return this->grid_->MCSize(); }
template<typename T>
Int DM::PartialUnionColStride() const { return this->grid_->MRSize(); }
template<typename T>
Int DM::DistSize() const { return this->grid_->VCSize(); }
template<typename T>
Int DM::CrossSize() const { return 1; }
template<typename T>
Int DM::RedundantSize() const { return 1; }

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
  BOTH( T,CIRC,CIRC); \
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
  OTHER(T,VC,  STAR); \
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
