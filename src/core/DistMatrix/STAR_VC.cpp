/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#define ColDist STAR
#define RowDist VC

#include "./setup.hpp"

namespace El {

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
    DEBUG_ONLY(CallStackEntry cse("[STAR,VC] = [MC,MR]"))
    DistMatrix<T,STAR,VR> A_STAR_VR( A );
    *this = A_STAR_VR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VC] = [MC,STAR]"))
    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR( new DistMatrix<T,MC,MR>(A) );
    std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(*A_MC_MR) );
    delete A_MC_MR.release(); // lowers memory highwater
    *this = *A_STAR_VR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VC] = [STAR,MR]"))
    DistMatrix<T,STAR,VR> A_STAR_VR( A );
    *this = A_STAR_VR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,VC] = [MD,STAR]"))
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VC] = [STAR,MD]"))
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VC] = [MR,MC]"))
    this->PartialRowAllToAllFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VC] = [MR,STAR]"))
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VC] = [STAR,MC]"))
    this->PartialRowFilterFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VC] = [VC,STAR]"))
    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR( new DistMatrix<T,MC,MR>(A) );
    std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(*A_MC_MR) );
    delete A_MC_MR.release(); // lowers memory highwater
    *this = *A_STAR_VR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VC] = [VR,STAR]"))
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[STAR,VC] = [STAR,VR]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const Grid& g = this->Grid();
    this->Resize( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;
    
    const Int height = this->Height();
    const Int localWidth = this->LocalWidth();
    const Int localWidthOfA = A.LocalWidth();
    
    const Int sendSize = height * localWidthOfA;
    const Int recvSize = height * localWidth;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int rankCM = g.VCRank();
    const Int rankRM = g.VRRank(); 

    const Int rowShift = this->RowShift();
    const Int rowShiftOfA = A.RowShift();

    // Compute which colmajor rank has the rowShift equal to our rowShiftOfA
    const Int sendRankCM = (rankCM+(p+rowShiftOfA-rowShift)) % p;

    // Compute which colmajor rank has the A rowShift that we need
    const Int recvRankRM = (rankRM+(p+rowShift-rowShiftOfA)) % p;
    const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

    T* buffer = this->auxMemory_.Require( sendSize + recvSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[sendSize];

    // Pack
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    EL_PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
    {
        const T* ACol = &ABuf[jLoc*ALDim];
        T* sendBufCol = &sendBuf[jLoc*height];
        MemCopy( sendBufCol, ACol, height );
    }

    // Communicate
    mpi::SendRecv
    ( sendBuf, sendSize, sendRankCM, 
      recvBuf, recvSize, recvRankCM, g.VCComm() );

    // Unpack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    EL_PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* recvBufCol = &recvBuf[jLoc*height];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        MemCopy( thisCol, recvBufCol, height );
    }
    this->auxMemory_.Release();
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,VC] = [STAR,STAR]"))
    this->RowFilterFrom( A );
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

    const Int rowAlign = this->RowAlign();
    const Int nLocal = this->LocalWidth();
    const Int pkgSize = mpi::Pad(m*MaxLength(n,p));
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
        for( Int t=0; t<p; ++t )
        {
            const Int tLocalWidth = Length( n, t, p );
            const Int q = (rowAlign+t) % p;
            for( Int jLoc=0; jLoc<tLocalWidth; ++jLoc )
            {
                const Int j = t + jLoc*p;
                for( Int i=0; i<m; ++i )
                    sendBuf[q*pkgSize+i+jLoc*m] = ABuf[i+j*ALDim];
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
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            for( Int i=0; i<m; ++i )
                buffer[i+jLoc*ldim] = recvBuf[i+jLoc*m];
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
mpi::Comm DM::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RowComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm DM::PartialRowComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm DM::PartialUnionRowComm() const { return this->grid_->MRComm(); }

template<typename T>
Int DM::ColStride() const { return 1; }
template<typename T>
Int DM::RowStride() const { return this->grid_->VCSize(); }
template<typename T>
Int DM::PartialRowStride() const { return this->grid_->MCSize(); }
template<typename T>
Int DM::PartialUnionRowStride() const { return this->grid_->MRSize(); }
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
  OTHER(T,STAR,VC  ); \
  BOTH( T,STAR,VR  ); \
  BOTH( T,VC,  STAR); \
  BOTH( T,VR,  STAR);

FULL(Int);
#ifndef EL_DISABLE_FLOAT
FULL(float);
#endif
FULL(double);

#ifndef EL_DISABLE_COMPLEX
#ifndef EL_DISABLE_FLOAT
FULL(Complex<float>);
#endif
FULL(Complex<double>);
#endif

} // namespace El
