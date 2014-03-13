/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

#define DM DistMatrix<T,MC,STAR>
#define GDM GeneralDistMatrix<T,MC,STAR>

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
DM::DistMatrix( const elem::Grid& grid, Int root )
: GDM(grid,root)
{ this->SetShifts(); }

template<typename T>
DM::DistMatrix( Int height, Int width, const elem::Grid& grid, Int root )
: GDM(grid,root)
{ this->SetShifts(); this->Resize(height,width); }

template<typename T>
DM::DistMatrix
( Int height, Int width, Int colAlign, Int rowAlign, const elem::Grid& grid,
  Int root )
: GDM(grid,root)
{
    this->SetShifts();
    this->Align(colAlign,rowAlign);
    this->Resize(height,width);
}

template<typename T>
DM::DistMatrix
( Int height, Int width, Int colAlign, Int rowAlign, Int ldim,
  const elem::Grid& grid, Int root )
: GDM(grid,root)
{
    this->SetShifts();
    this->Align(colAlign,rowAlign);
    this->Resize(height,width,ldim);
}

template<typename T>
DM::DistMatrix
( Int height, Int width, Int colAlign, Int rowAlign,
  const T* buffer, Int ldim, const elem::Grid& grid, Int root )
: GDM(grid,root)
{ this->LockedAttach(height,width,colAlign,rowAlign,buffer,ldim,grid,root); }

template<typename T>
DM::DistMatrix
( Int height, Int width, Int colAlign, Int rowAlign,
  T* buffer, Int ldim, const elem::Grid& grid, Int root )
: GDM(grid,root)
{ this->Attach(height,width,colAlign,rowAlign,buffer,ldim,grid,root); }

template<typename T>
DM::DistMatrix( const DM& A )
: GDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::DistMatrix"))
    this->SetShifts();
    if( &A != this ) 
        *this = A;
    else
        LogicError("Tried to construct [MC,* ] with itself");
}

template<typename T>
template<Dist U,Dist V>
DM::DistMatrix( const DistMatrix<T,U,V>& A )
: GDM(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ]::DistMatrix"))
    this->SetShifts();
    if( MC != U || STAR != V || 
        reinterpret_cast<const DM*>(&A) != this ) 
        *this = A;
    else
        LogicError("Tried to construct [MC,* ] with itself");
}

template<typename T>
DM::DistMatrix( DM&& A ) noexcept : GDM(std::move(A)) { }

template<typename T> DM::~DistMatrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MC,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [MC,MR]"))
    A.RowAllGather( *this );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DM& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[MC,* ] = [MC,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
    this->AlignColsAndResize( A.ColAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlign() == A.ColAlign() )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MC,* ] <- [MC,* ]." << std::endl;
#endif
        const Int rank = g.Row();
        const Int r = g.Height();

        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();

        const Int sendRank = (rank+r+colAlign-colAlignOfA) % r;
        const Int recvRank = (rank+r+colAlignOfA-colAlign) % r;

        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localHeightOfA = A.LocalHeight();

        const Int sendSize = localHeightOfA * width;
        const Int recvSize = localHeight * width;

        T* buffer = this->auxMemory_.Require( sendSize + recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ABuf[j*ALDim];
            T* sendBufCol = &sendBuf[j*localHeightOfA];
            MemCopy( sendBufCol, ACol, localHeightOfA );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank,
          recvBuf, recvSize, recvRank, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* recvBufCol = &recvBuf[j*localHeight];
            T* thisCol = &thisBuf[j*thisLDim];
            MemCopy( thisCol, recvBufCol, localHeight );
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [* ,MR]"))
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(true,false,this->ColAlign(),0,g);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [MD,* ]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [* ,MD]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [MR,MC]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(A) );
    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(true,this->ColAlign(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,MR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [MR,* ]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(A) );
    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(true,this->ColAlign(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,MC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [* ,MC]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,MR,MC>> 
        A_MR_MC( new DistMatrix<T,MR,MC>(A) );
    std::unique_ptr<DistMatrix<T,VR,STAR>> 
        A_VR_STAR( new DistMatrix<T,VR,STAR>(g) );
    *A_VR_STAR = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    std::unique_ptr<DistMatrix<T,VC,STAR>> 
        A_VC_STAR( new DistMatrix<T,VC,STAR>(true,this->ColAlign(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,VC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [VC,* ]"))
    A.PartialColAllGather( *this );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,VC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [* ,VC]"))
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,STAR,VR>> 
        A_STAR_VR( new DistMatrix<T,STAR,VR>(A) );
    std::unique_ptr<DistMatrix<T,MC,MR>> 
        A_MC_MR
        ( new DistMatrix<T,MC,MR>(true,false,this->ColAlign(),0,g) );
    *A_MC_MR = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater
    *this = *A_MC_MR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,VR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [VR,* ]"))
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR> A_VC_STAR(true,this->ColAlign(),g);
    A_VC_STAR = A;
    *this = A_VC_STAR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,VR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [* ,VR]"))
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(true,false,this->ColAlign(),0,g);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [* ,* ]"))
    this->ColFilterFrom( A );
    return *this;
}

template<typename T>
DM&
DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,* ] = [o ,o ]"))
    DistMatrix<T,MC,MR> A_MC_MR( A.Grid() );
    A_MC_MR.AlignWith( *this );
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
DM&
DM::operator=( DM&& A )
{
    if( this->Viewing() && !A.Viewing() )
    {
        const DM& AConst = A;
        this->operator=( AConst );
    }
    else
    {
        GDM::operator=( std::move(A) );
    }
    return *this;
}

// Realignment
// -----------

template<typename T>
void
DM::AlignWith( const elem::DistData& data )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,STAR]::AlignWith"))
    this->SetGrid( *data.grid );
    if( data.colDist == MC )
        this->AlignCols( data.colAlign );
    else if( data.rowDist == MC )
        this->AlignCols( data.rowAlign );
    else if( data.colDist == VC )
        this->AlignCols( data.colAlign % this->ColStride() );
    else if( data.rowDist == VC )
        this->AlignCols( data.rowAlign % this->ColStride() );
    DEBUG_ONLY(else LogicError("Nonsensical alignment"))
}

template<typename T>
void
DM::AlignColsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

// Basic queries
// =============

template<typename T>
elem::DistData DM::DistData() const { return elem::DistData(*this); }

template<typename T>
mpi::Comm DM::DistComm() const { return this->grid_->ColComm(); }
template<typename T>
mpi::Comm DM::RedundantComm() const { return this->grid_->RowComm(); }
template<typename T>
mpi::Comm DM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::ColComm() const { return this->grid_->ColComm(); }
template<typename T>
mpi::Comm DM::RowComm() const { return mpi::COMM_SELF; }

template<typename T>
Int DM::ColStride() const { return this->grid_->Height(); }
template<typename T>
Int DM::RowStride() const { return 1; }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class DistMatrix<T,MC,STAR>
#define COPY(T,U,V) \
  template DistMatrix<T,MC,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  MR); \
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
