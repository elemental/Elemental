/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

#define DM DistMatrix<T,STAR,MC>
#define GDM GeneralDistMatrix<T,STAR,MC>

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
DM::DistMatrix( const elem::Grid& g )
: GDM(g)
{ this->SetShifts(); }

template<typename T>
DM::DistMatrix( Int height, Int width, const elem::Grid& g )
: GDM(g)
{ this->SetShifts(); this->Resize(height,width); }

template<typename T>
DM::DistMatrix( Int height, Int width, Int rowAlign, const elem::Grid& g )
: GDM(g)
{ this->Align(0,rowAlign); this->Resize(height,width); }

template<typename T>
DM::DistMatrix
( Int height, Int width, Int rowAlign, Int ldim, const elem::Grid& g )
: GDM(g)
{ this->Align(0,rowAlign); this->Resize(height,width,ldim); }

template<typename T>
DM::DistMatrix
( Int height, Int width, Int rowAlign, const T* buffer, Int ldim, 
  const elem::Grid& g )
: GDM(g)
{ this->LockedAttach(height,width,0,rowAlign,buffer,ldim,g); }

template<typename T>
DM::DistMatrix
( Int height, Int width, Int rowAlign, T* buffer, Int ldim, 
  const elem::Grid& g )
: GDM(g)
{ this->Attach(height,width,0,rowAlign,buffer,ldim,g); }

template<typename T>
DM::DistMatrix( const DM& A )
: GDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MC]::DistMatrix"))
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [* ,MC] with itself");
}

template<typename T>
template<Dist U,Dist V>
DM::DistMatrix( const DistMatrix<T,U,V>& A )
: GDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MC]::DistMatrix"))
    this->SetShifts();
    if( STAR != U || MC != V || 
        reinterpret_cast<const DM*>(&A) != this ) 
        *this = A;
    else
        LogicError("Tried to construct [* ,MC] with itself");
}

template<typename T> DM::DistMatrix( DM&& A ) : GDM(std::move(A)) { }

template<typename T> DM::~DistMatrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MC,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [MC,MR]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(A) );

    std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(true,this->RowAlign(),g) );
    *A_STAR_VC = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_STAR_VC;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [MC,* ]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR( new DistMatrix<T,MC,MR>(A) );

    std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(*A_MC_MR) );
    delete A_MC_MR.release(); // lowers memory highwater

    std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(true,this->RowAlign(),g) );
    *A_STAR_VC = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_STAR_VC;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[* ,MC] = [* ,MR]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
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

        // A[* ,VR] <- A[* ,MR]
        {
            const Int shift = Shift(rankRM,rowAlignOfA,p);
            const Int offset = (shift-rowShiftOfA) / c;
            const Int thisLocalWidth = Length(width,shift,p);

            const T* ABuf = A.LockedBuffer();
            const Int ALDim = A.LDim();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                sendBuf[jLoc] = ABuf[(offset+jLoc*r)*ALDim];
        }

        // A[* ,VC] <- A[* ,VR]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankCM,
          recvBuf, portionSize, recvRankCM, g.VCComm() );

        // A[* ,MC] <- A[* ,VC]
        mpi::AllGather
        ( recvBuf, portionSize,
          sendBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
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
        std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
        ( new DistMatrix<T,STAR,VR>(g) );
        *A_STAR_VR = A;

        std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
        ( new DistMatrix<T,STAR,VC>(true,this->RowAlign(),g) );
        *A_STAR_VC = *A_STAR_VR;
        delete A_STAR_VR.release(); // lowers memory highwater

        *this = *A_STAR_VC;
    }
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [MD,* ]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [* ,MD]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [MR,MC]"))
    A.ColAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [MR,* ]"))
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DM& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[* ,MC] = [* ,MC]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
    this->AlignRowsAndResize( A.RowAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlign() == A.RowAlign() )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MC] <- [* ,MC]." << std::endl;
#endif
        const Int rank = g.Row();
        const Int r = g.Height();

        const Int rowAlign = this->RowAlign();
        const Int rowAlignOfA = A.RowAlign();

        const Int sendRank = (rank+r+rowAlign-rowAlignOfA) % r;
        const Int recvRank = (rank+r+rowAlignOfA-rowAlign) % r;

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfA = A.LocalWidth();

        const Int sendSize = height * localWidthOfA;
        const Int recvSize = height * localWidth;

        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuf[jLoc*ALDim];
            T* sendBufCol = &sendBuf[jLoc*height];
            MemCopy( sendBufCol, ACol, height );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank, 
          recvBuf, recvSize, recvRank, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* recvBufCol = &recvBuf[jLoc*height];
            T* thisCol = &thisBuf[jLoc*thisLDim];
            MemCopy( thisCol, recvBufCol, height );
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [VC,* ]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(A) );

    std::unique_ptr<DistMatrix<T,MR,MC>> A_MR_MC
    ( new DistMatrix<T,MR,MC>(false,true,0,this->RowAlign(),g) );
    *A_MR_MC = *A_VR_STAR;
    delete A_VR_STAR.release();

    *this = *A_MR_MC;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [* ,VC]"))
    A.PartialRowAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [VR,* ]"))
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [* ,VR]"))
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VC> A_STAR_VC(true,this->RowAlign(),g);
    *this = A_STAR_VC = A;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [* ,* ]"))
    this->RowFilterFrom( A );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MC] = [o ,o ]"))
    DistMatrix<T,MR,MC> A_MR_MC( A.Grid() );
    A_MR_MC.AlignWith( *this );
    A_MR_MC = A;
    *this = A_MR_MC;
    return *this;
}

template<typename T>
DM&
DM::operator=( DM&& A )
{
    GDM::operator=( std::move(A) );
    return *this;
}

// Realignment
// -----------

template<typename T>
void
DM::AlignWith( const elem::DistData& data )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MC]::AlignWith"))
    this->SetGrid( *data.grid );

    if( data.colDist == MC )
        this->AlignRows( data.colAlign );
    else if( data.rowDist == MC )
        this->AlignRows( data.rowAlign );
    else if( data.colDist == VC )
        this->AlignRows( data.colAlign % this->RowStride() );
    else if( data.rowDist == VC )
        this->AlignRows( data.rowAlign % this->RowStride() );
    DEBUG_ONLY(else LogicError("Nonsensical alignment"))
}

template<typename T>
void
DM::AlignRowsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

// Basic queries
// =============

template<typename T>
elem::DistData DM::DistData() const { return elem::DistData(*this); }

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
Int DM::RowStride() const { return this->grid_->Height(); }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class DistMatrix<T,STAR,MC>
#define COPY(T,U,V) \
  template DistMatrix<T,STAR,MC>::DistMatrix( const DistMatrix<T,U,V>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  MR  ); \
  COPY(T,MC,  STAR); \
  COPY(T,MD,  STAR); \
  COPY(T,MR,  MC  ); \
  COPY(T,MR,  STAR); \
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
