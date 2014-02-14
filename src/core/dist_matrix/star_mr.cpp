/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

#define DM DistMatrix<T,STAR,MR>
#define GDM GeneralDistMatrix<T,STAR,MR>

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
DM::DistMatrix
( Int height, Int width, Int rowAlign, const elem::Grid& g )
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
    DEBUG_ONLY(CallStackEntry cse("[* ,MR]::DistMatrix"))
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [* ,MR] with itself");
}

template<typename T>
template<Dist U,Dist V>
DM::DistMatrix( const DistMatrix<T,U,V>& A )
: GDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MR]::DistMatrix"))
    this->SetShifts();
    if( STAR != U || MR != V || 
        reinterpret_cast<const DM*>(&A) != this )   
        *this = A;
    else
        LogicError("Tried to construct [* ,MR] with itself");
}

template<typename T>
DM::DistMatrix( DM&& A ) : GDM(std::move(A)) { }

template<typename T> DM::~DistMatrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MC,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [MC,MR]"))
    A.ColAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[* ,MR] = [MC,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        if( this->Viewing() )
            this->AssertSameSize( A.Height(), A.Width() );
    )
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(false,true,0,this->RowAlign(),g);

    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DM& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[* ,MR] = [* ,MR]");
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
            std::cerr << "Unaligned [* ,MR] <- [* ,MR]." << std::endl;
#endif
        const Int rank = g.Col();
        const Int c = g.Width();

        const Int rowAlign = this->RowAlign();
        const Int rowAlignOfA = A.RowAlign();

        const Int sendRank = (rank+c+rowAlign-rowAlignOfA) % c;
        const Int recvRank = (rank+c+rowAlignOfA-rowAlign) % c;

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
          recvBuf, recvSize, recvRank, g.RowComm() );

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
DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [MD,* ]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [* ,MD]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [MR,MC]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(g) );
    *A_STAR_VC = A;

    std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(true,this->RowAlign(),g) );
    *A_STAR_VR = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    *this = *A_STAR_VR;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [MR,* ]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(g) );
    *A_VR_STAR = A;

    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR
    ( new DistMatrix<T,MC,MR>(false,true,0,this->RowAlign(),g) );
    *A_MC_MR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_MC_MR;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [* ,MC]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(g) );
    *A_STAR_VC = A;

    std::unique_ptr<DistMatrix<T,STAR,VR>> A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(true,this->RowAlign(),g) );
    *A_STAR_VR = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_MC_MR;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [VC,* ]"))
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(false,true,0,this->RowAlign(),g);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [* ,VC]"))
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR> A_STAR_VR(true,this->RowAlign(),g);
    A_STAR_VR = A;
    *this = A_STAR_VR;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [VR,* ]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = A;

    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR
    ( new DistMatrix<T,MC,MR>(false,true,0,this->RowAlign(),g) );
    *A_MC_MR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_MC_MR;
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [* ,VR]"))
    A.PartialRowAllGather( *this );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [* ,* ]"))
    this->RowFilterFrom( A );
    return *this;
}

template<typename T>
const DM&
DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [o ,o ]"))
    DistMatrix<T,MC,MR> A_MC_MR( A );
    A_MC_MR.AlignWith( *this );
    A_MC_MR = A;
    *this = A_MC_MR;
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
    DEBUG_ONLY(CallStackEntry cse("[* ,MR]::AlignWith"))
    this->SetGrid( *data.grid );

    if( data.colDist == MR )
        this->AlignRows( data.colAlign );
    else if( data.rowDist == MR )
        this->AlignRows( data.rowAlign );
    else if( data.colDist == VR )
        this->AlignRows( data.colAlign % this->RowStride() );
    else if( data.rowDist == VR )
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
mpi::Comm DM::DistComm() const { return this->grid_->MRComm(); }
template<typename T>
mpi::Comm DM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RedundantComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm DM::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RowComm() const { return this->grid_->MRComm(); }

template<typename T>
Int DM::ColStride() const { return 1; }
template<typename T>
Int DM::RowStride() const { return this->grid_->Width(); }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class DistMatrix<T,STAR,MR>
#define COPY(T,U,V) \
  template DistMatrix<T,STAR,MR>::DistMatrix( const DistMatrix<T,U,V>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  MR  ); \
  COPY(T,MC,  STAR); \
  COPY(T,MD,  STAR); \
  COPY(T,MR,  MC  ); \
  COPY(T,MR,  STAR); \
  COPY(T,STAR,MC  ); \
  COPY(T,STAR,MD  ); \
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
