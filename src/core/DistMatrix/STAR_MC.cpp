/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#define ColDist STAR
#define RowDist MC

#include "./setup.hpp"

namespace El {

// Public section
// ##############

// Assignment and reconfiguration
// ==============================

// Return a view
// -------------
template<typename T>
DM DM::operator()( Range<Int> indVert, Range<Int> indHorz )
{
    DEBUG_ONLY(CallStackEntry cse("DM[STAR,MC]( ind, ind )"))
    if( this->Locked() )
        return LockedView( *this, indVert, indHorz );
    else
        return View( *this, indVert, indHorz );
}

template<typename T>
const DM DM::operator()( Range<Int> indVert, Range<Int> indHorz ) const
{
    DEBUG_ONLY(CallStackEntry cse("DM[STAR,MC]( ind, ind )"))
    return LockedView( *this, indVert, indHorz );
}

// Make a copy
// -----------
template<typename T>
DM& DM::operator=( const DM& A )
{
    DEBUG_ONLY(CallStackEntry cse("DM[STAR,MC] = DM[STAR,MC]"))
    copy::Translate( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MC,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [MC,MR]"))
    auto A_STAR_VR = MakeUnique<DistMatrix<T,STAR,VR>>( A );
    auto A_STAR_VC = MakeUnique<DistMatrix<T,STAR,VC>>( this->Grid() );
    A_STAR_VC->AlignRowsWith(*this);
    *A_STAR_VC = *A_STAR_VR;
    A_STAR_VR.reset();

    *this = *A_STAR_VC;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [MC,STAR]"))
    auto A_MC_MR = MakeUnique<DistMatrix<T,MC,MR>>( A );
    auto A_STAR_VR = MakeUnique<DistMatrix<T,STAR,VR>>( *A_MC_MR );
    A_MC_MR.reset(); 

    auto A_STAR_VC = MakeUnique<DistMatrix<T,STAR,VC>>( this->Grid() );
    A_STAR_VC->AlignRowsWith(*this);
    *A_STAR_VC = *A_STAR_VR;
    A_STAR_VR.reset();

    *this = *A_STAR_VC;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[STAR,MC] = [STAR,MR]");
        AssertSameGrids( *this, A );
        this->AssertNotLocked();
    )
    const El::Grid& g = this->Grid();
    if( A.Height() == 1 )
    {
        this->Resize( 1, A.Width() );
        if( !this->Participating() )
            return *this;

        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int myRow = g.Row();
        const Int rankCM = g.VCRank();
        const Int rankRM = g.VRRank();
        const Int rowAlign = this->RowAlign();
        const Int rowShift = this->RowShift();
        const Int rowAlignOfA = A.RowAlign();
        const Int rowShiftOfA = A.RowShift();

        const Int width = this->Width();
        const Int maxLocalVectorWidth = MaxLength(width,p);
        const Int portionSize = mpi::Pad( maxLocalVectorWidth );

        const Int rowShiftVC = Shift(rankCM,rowAlign,p);
        const Int rowShiftVROfA = Shift(rankRM,rowAlignOfA,p);
        const Int sendRankCM = (rankCM+(p+rowShiftVROfA-rowShiftVC)) % p;
        const Int recvRankRM = (rankRM+(p+rowShiftVC-rowShiftVROfA)) % p;
        const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        T* buffer = this->auxMemory_.Require( (c+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        // A[STAR,VR] <- A[STAR,MR]
        {
            const Int shift = Shift(rankRM,rowAlignOfA,p);
            const Int offset = (shift-rowShiftOfA) / c;
            const Int thisLocalWidth = Length(width,shift,p);

            const T* ABuf = A.LockedBuffer();
            const Int ALDim = A.LDim();
            EL_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                sendBuf[jLoc] = ABuf[(offset+jLoc*r)*ALDim];
        }

        // A[STAR,VC] <- A[STAR,VR]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankCM,
          recvBuf, portionSize, recvRankCM, g.VCComm() );

        // A[STAR,MC] <- A[STAR,VC]
        mpi::AllGather
        ( recvBuf, portionSize,
          sendBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        EL_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            const T* data = &sendBuf[k*portionSize];
            const Int shift = Shift_(myRow+r*k,rowAlign,p);
            const Int offset = (shift-rowShift) / r;
            const Int thisLocalWidth = Length_(width,shift,p);
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                thisBuf[(offset+jLoc*c)*thisLDim] = data[jLoc];
        }
        this->auxMemory_.Release();
    }
    else
    {
        auto A_STAR_VR = MakeUnique<DistMatrix<T,STAR,VR>>( A );
        auto A_STAR_VC = MakeUnique<DistMatrix<T,STAR,VC>>( g );
        A_STAR_VC->AlignRowsWith(*this);
        *A_STAR_VC = *A_STAR_VR;
        A_STAR_VR.reset();
        *this = *A_STAR_VC;
    }
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [MD,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [STAR,MD]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [MR,MC]"))
    copy::ColAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [MR,STAR]"))
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [VC,STAR]"))
    auto A_VR_STAR = MakeUnique<DistMatrix<T,VR,STAR>>( A );
    auto A_MR_MC = MakeUnique<DistMatrix<T,MR,MC>>( this->Grid() );
    A_MR_MC->AlignRowsWith(*this);
    *A_MR_MC = *A_VR_STAR;
    A_VR_STAR.reset();

    *this = *A_MR_MC;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [STAR,VC]"))
    copy::PartialRowAllGather( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [VR,STAR]"))
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [STAR,VR]"))
    DistMatrix<T,STAR,VC> A_STAR_VC(this->Grid());
    A_STAR_VC.AlignRowsWith(*this);
    *this = A_STAR_VC = A;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [STAR,STAR]"))
    copy::RowFilter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[STAR,MC] = [CIRC,CIRC]"))
    DistMatrix<T,MR,MC> A_MR_MC( A.Grid() );
    A_MR_MC.AlignWith( *this );
    A_MR_MC = A;
    *this = A_MR_MC;
    return *this;
}

// Basic queries
// =============

template<typename T>
mpi::Comm DM::DistComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm DM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RedundantComm() const { return this->grid_->MRComm(); }
template<typename T>
mpi::Comm DM::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RowComm() const { return this->grid_->MCComm(); }

template<typename T>
Int DM::ColStride() const { return 1; }
template<typename T>
Int DM::RowStride() const { return this->grid_->MCSize(); }
template<typename T>
Int DM::DistSize() const { return this->grid_->MCSize(); }
template<typename T>
Int DM::CrossSize() const { return 1; }
template<typename T>
Int DM::RedundantSize() const { return this->grid_->MRSize(); }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

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
#define PROTO(T) \
  template class DistMatrix<T,ColDist,RowDist>; \
  BOTH( T,CIRC,CIRC); \
  BOTH( T,MC,  MR  ); \
  BOTH( T,MC,  STAR); \
  BOTH( T,MD,  STAR); \
  BOTH( T,MR,  MC  ); \
  BOTH( T,MR,  STAR); \
  OTHER(T,STAR,MC  ); \
  BOTH( T,STAR,MD  ); \
  BOTH( T,STAR,MR  ); \
  BOTH( T,STAR,STAR); \
  BOTH( T,STAR,VC  ); \
  BOTH( T,STAR,VR  ); \
  BOTH( T,VC,  STAR); \
  BOTH( T,VR,  STAR);

#include "El/macros/Instantiate.h"

} // namespace El
