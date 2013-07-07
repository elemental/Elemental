/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

template<typename T,typename Int>
DistMatrix<T,MR,STAR,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->SetShifts(); } 

template<typename T,typename Int>
DistMatrix<T,MR,STAR,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T,typename Int>
DistMatrix<T,MR,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->Align(colAlignment,0); this->ResizeTo(height,width); }

template<typename T,typename Int>
DistMatrix<T,MR,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->Align(colAlignment,0); this->ResizeTo(height,width,ldim); }

template<typename T,typename Int>
DistMatrix<T,MR,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->LockedAttach( height, width, colAlignment, buffer, ldim, g ); }

template<typename T,typename Int>
DistMatrix<T,MR,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>(g)
{ this->Attach( height, width, colAlignment, buffer, ldim, g ); }

template<typename T,typename Int>
DistMatrix<T,MR,STAR,Int>::DistMatrix( const DistMatrix<T,MR,STAR,Int>& A )
: AbstractDistMatrix<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MR,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [MR,* ] with itself");
}

template<typename T,typename Int>
template<Distribution U,Distribution V>
DistMatrix<T,MR,STAR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[MR,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( MR != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,MR,STAR,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [MR,* ] with itself");
}

template<typename T,typename Int>
DistMatrix<T,MR,STAR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
elem::DistData<Int>
DistMatrix<T,MR,STAR,Int>::DistData() const
{
    elem::DistData<Int> data;
    data.colDist = MR;
    data.rowDist = STAR;
    data.colAlignment = this->colAlignment_;
    data.rowAlignment = 0;
    data.root = 0;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T,typename Int>
Int
DistMatrix<T,MR,STAR,Int>::ColStride() const
{ return this->grid_->Width(); }

template<typename T,typename Int>
Int
DistMatrix<T,MR,STAR,Int>::RowStride() const
{ return 1; }

template<typename T,typename Int>
Int
DistMatrix<T,MR,STAR,Int>::ColRank() const
{ return this->grid_->Col(); }

template<typename T,typename Int>
Int
DistMatrix<T,MR,STAR,Int>::RowRank() const
{ return 0; }

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::AlignWith( const elem::DistData<Int>& data )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::AlignWith");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( data.colDist == MR )
        this->colAlignment_ = data.colAlignment;
    else if( data.rowDist == MR )
        this->colAlignment_ = data.rowAlignment;
    else if( data.colDist == VR )
        this->colAlignment_ = data.colAlignment % this->ColStride();
    else if( data.rowDist == VR )
        this->colAlignment_ = data.rowAlignment % this->ColStride();
#ifndef RELEASE
    else throw std::logic_error("Nonsensical alignment");
#endif
    this->constrainedColAlignment_ = true;
    this->SetShifts();
}

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::AlignWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::AlignColsWith( const elem::DistData<Int>& data )
{ this->AlignWith( data ); }

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::AlignColsWith( const AbstractDistMatrix<T,Int>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::Attach
( Int height, Int width, Int colAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::Attach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->viewType_ = VIEW;
    this->SetColShift();
    if( this->Participating() )
    {
        const Int localHeight = Length(height,this->colShift_,g.Width());
        this->matrix_.Attach_( localHeight, width, buffer, ldim );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::LockedAttach
( Int height, Int width, Int colAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->viewType_ = LOCKED_VIEW;
    this->SetColShift();
    if( this->Participating() )
    {
        const Int localHeight = Length(height,this->colShift_,g.Width());
        this->matrix_.LockedAttach_( localHeight, width, buffer, ldim );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),this->Grid().Width()), width );
}

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),this->Grid().Width()), width, ldim );
}

template<typename T,typename Int>
T
DistMatrix<T,MR,STAR,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::Get");
    this->AssertValidEntry( i, j );
    if( !this->Participating() )
        throw std::logic_error("Should only be called by grid members");
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // columns within each process row
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();

    T u;
    if( g.Col() == ownerCol )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        u = this->GetLocal(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerCol, g.RowComm() );
    return u;
}

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();

    if( g.Col() == ownerCol )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        this->SetLocal(iLoc,j,u);
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();

    if( g.Col() == ownerCol )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        this->UpdateLocal(iLoc,j,u);
    }
}

//
// Utility functions, e.g., SumOverCol
//

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::SumOverCol()
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::SumOverCol");
    this->AssertNotLocked();
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;
    
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localSize = mpi::Pad( localHeight*width );

    T* buffer = this->auxMemory_.Require( 2*localSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[localSize];

    // Pack
    T* thisBuffer = this->Buffer();
    const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
    {
        const T* thisCol = &thisBuffer[j*thisLDim];
        T* sendBufCol = &sendBuf[j*localHeight];
        MemCopy( sendBufCol, thisCol, localHeight );
    }

    // AllReduce sum
    mpi::AllReduce
    ( sendBuf, recvBuf, localSize, mpi::SUM, g.ColComm() );

    // Unpack
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
    {
        const T* recvBufCol = &recvBuf[j*localHeight];
        T* thisCol = &thisBuffer[j*thisLDim];
        MemCopy( thisCol, recvBufCol, localHeight );
    }
    this->auxMemory_.Release();
}

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::AdjointFrom( const DistMatrix<T,MC,MR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::AdjointFrom");
#endif
    this->TransposeFrom( A, true );
}

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::TransposeFrom
( const DistMatrix<T,MC,MR,Int>& A, bool conjugate )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::TransposeFrom");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Width(), A.Height() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.RowAlignment();
            this->SetColShift();
        }
        this->ResizeTo( A.Width(), A.Height() );
    }
    if( !this->Participating() )
        return;

    if( this->ColAlignment() == A.RowAlignment() )
    {
        const Int r = g.Height();

        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalWidth = MaxLength(width,r);
        const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack 
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
        if( conjugate )
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &sendBuf[jLoc*localHeight];
                const T* sourceCol = &ABuffer[jLoc];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc] = Conj( sourceCol[iLoc*ALDim] );
            }
        }
        else
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &sendBuf[jLoc*localHeight];
                const T* sourceCol = &ABuffer[jLoc];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*ALDim];
            }
        }

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int colAlignmentOfA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int rowShift = Shift_( k, colAlignmentOfA, r );
            const Int localWidth = Length_( width, rowShift, r );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuffer[(rowShift+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MR,* ]::TransposeFrom" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int col = g.Col();

        const Int colAlignment = this->ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int sendCol = (col+c+colAlignment-rowAlignmentOfA) % c;
        const Int recvCol = (col+c+rowAlignmentOfA-colAlignment) % c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localHeightOfA = A.LocalHeight();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalHeight = MaxLength(height,c);
        const Int maxLocalWidth = MaxLength(width,r);
        const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[portionSize];

        // Pack the currently owned local data of A into the second buffer
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
        if( conjugate )
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &secondBuf[jLoc*localWidthOfA];
                const T* sourceCol = &ABuffer[jLoc];
                for( Int iLoc=0; iLoc<localWidthOfA; ++iLoc )
                    destCol[iLoc] = Conj( sourceCol[iLoc*ALDim] );
            }
        }
        else
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &secondBuf[jLoc*localWidthOfA];
                const T* sourceCol = &ABuffer[jLoc];
                for( Int iLoc=0; iLoc<localWidthOfA; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*ALDim];
            }
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuf, portionSize, sendCol, 0,
          firstBuf,  portionSize, recvCol, mpi::ANY_TAG, g.RowComm() );

        // Use the output of the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.ColComm() );

        // Unpack the contents of each member of the process col
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int colAlignmentOfA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShift = Shift_( k, colAlignmentOfA, r );
            const Int localWidth = Length_( width, rowShift, r );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuffer[(rowShift+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [MC,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(g) );
    *A_VC_STAR = A;

    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VR_STAR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_VR_STAR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
    {
        if( !this->Viewing() )
            this->ResizeTo( A.Height(), A.Width() );
        return *this;
    }

    if( A.Width() == 1 )
    {
        if( !this->Viewing() )
            this->ResizeTo( A.Height(), 1 );
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int myCol = g.Col();
        const Int rankCM = g.VCRank();
        const Int rankRM = g.VRRank();
        const Int colAlignment = this->ColAlignment();
        const Int colShift = this->ColShift();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int maxLocalVectorHeight = MaxLength(height,p);
        const Int portionSize = mpi::Pad( maxLocalVectorHeight );

        const Int colShiftVR = Shift(rankRM,colAlignment,p);
        const Int colShiftVCOfA = Shift(rankCM,colAlignmentOfA,p);
        const Int sendRankRM = (rankRM+(p+colShiftVCOfA-colShiftVR)) % p;
        const Int recvRankCM = (rankCM+(p+colShiftVR-colShiftVCOfA)) % p;
        const Int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        // A[VC,* ] <- A[MC,* ]
        {
            const Int shift = Shift(rankCM,colAlignmentOfA,p);
            const Int offset = (shift-colShiftOfA) / r;
            const Int thisLocalHeight = Length(height,shift,p);

            const T* ABuffer = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                sendBuf[iLoc] = ABuffer[offset+iLoc*c];
        }

        // A[VR,* ] <- A[VC,* ]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankRM, 0,
          recvBuf, portionSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

        // A[MR,* ] <- A[VR,* ]
        mpi::AllGather
        ( recvBuf, portionSize,
          sendBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &sendBuf[k*portionSize];
            const Int shift = Shift_(myCol+c*k,colAlignment,p);
            const Int offset = (shift-colShift) / c;
            const Int thisLocalHeight = Length_(height,shift,p);
            for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                thisBuffer[offset+iLoc*r] = data[iLoc];
        }
        this->auxMemory_.Release();
    }
    else
    {
        std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
        ( new DistMatrix<T,VC,STAR,Int>(g) );
        *A_VC_STAR = A;

        std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
        ( new DistMatrix<T,VR,STAR,Int>(true,this->ColAlignment(),g) );
        *A_VR_STAR = *A_VC_STAR;
        delete A_VC_STAR.release(); // lowers memory highwater

        *this = *A_VR_STAR;
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [* ,MR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MC,MR,Int> > A_MC_MR
    ( new DistMatrix<T,MC,MR,Int>(g) );
    *A_MC_MR   = A;

    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(g) );
    *A_VC_STAR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VR_STAR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_VR_STAR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [MD,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MR,* ] = [MD,* ] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [* ,MD]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    throw std::logic_error("[MR,* ] = [* ,MD] not yet implemented");
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            this->SetColShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !this->Participating() )
        return *this;

    if( this->ColAlignment() == A.ColAlignment() )
    {
        const Int r = g.Height();

        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidth = MaxLength(width,r);
        const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack 
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuffer[jLoc*ALDim];
            T* sendBufCol = &sendBuf[jLoc*localHeight];
            MemCopy( sendBufCol, ACol, localHeight );
        }

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int rowAlignmentOfA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int rowShift = Shift_( k, rowAlignmentOfA, r );
            const Int localWidth = Length_( width, rowShift, r );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuffer[(rowShift+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MR,* ] <- [MR,MC]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int col = g.Col();

        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
        const Int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localHeightOfA = A.LocalHeight();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalHeight = MaxLength(height,c);
        const Int maxLocalWidth = MaxLength(width,r);
        const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[portionSize];

        // Pack the currently owned local data of A into the second buffer
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuffer[jLoc*ALDim];
            T* secondBufCol = &secondBuf[jLoc*localHeightOfA];
            MemCopy( secondBufCol, ACol, localHeightOfA );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuf, portionSize, sendCol, 0,
          firstBuf,  portionSize, recvCol, mpi::ANY_TAG, g.RowComm() );

        // Use the output of the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.ColComm() );

        // Unpack the contents of each member of the process col
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int rowAlignmentOfA = A.RowAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShift = Shift_( k, rowAlignmentOfA, r );
            const Int localWidth = Length_( width, rowShift, r );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuffer[(rowShift+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            if( this->Participating() )
                this->colShift_ = A.ColShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !this->Participating() )
        return *this;

    if( this->ColAlignment() == A.ColAlignment() )
    {
        this->matrix_ = A.LockedMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MR,* ] <- [MR,* ]." << std::endl;
#endif
        const Int rank = g.Col();
        const Int c = g.Width();

        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int sendRank = (rank+c+colAlignment-colAlignmentOfA) % c;
        const Int recvRank = (rank+c+colAlignmentOfA-colAlignment) % c;

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
        const T* ABuffer = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ABuffer[j*ALDim];
            T* sendBufCol = &sendBuf[j*localHeightOfA];
            MemCopy( sendBufCol, ACol, localHeightOfA );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank, 0,
          recvBuf, recvSize, recvRank, mpi::ANY_TAG, g.RowComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* recvBufCol = &recvBuf[j*localHeight];
            T* thisCol = &thisBuffer[j*thisLDim];
            MemCopy( thisCol, recvBufCol, localHeight );
        }

        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [* ,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MR,MC,Int> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VR,STAR,Int> A_VR_STAR(true,this->ColAlignment(),g);

    A_VR_STAR = A;
    *this = A_VR_STAR;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [* ,VC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MR,MC,Int> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "[MR,* ] <- [VR,* ] potentially causes a large amount of cache-"
          "thrashing. If possible avoid it by performing the redistribution "
          "with a (conjugate)-transpose: \n" << 
          "  [* ,MR].(Conjugate)TransposeFrom([VR,* ])" << std::endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment() % g.Width();
            this->SetColShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !this->Participating() )
        return *this;

    if( this->ColAlignment() == A.ColAlignment() % g.Width() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int col = g.Col();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLength(height,p);
        const Int portionSize = mpi::Pad( maxLocalHeightOfA*width );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ABuffer[j*ALDim];
            T* sendBufCol = &sendBuf[j*localHeightOfA];
            MemCopy( sendBufCol, ACol, localHeightOfA );
        }

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int colShift = this->ColShift();
        const Int colAlignmentOfA = A.ColAlignment();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int colShiftOfA = Shift_( col+c*k, colAlignmentOfA, p );
            const Int colOffset = (colShiftOfA-colShift) / c;
            const Int localHeight = Length_( height, colShiftOfA, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &thisBuffer[colOffset+j*thisLDim];
                const T* sourceCol = &data[j*localHeight];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc*r] = sourceCol[iLoc];
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MR,* ] <- [VR,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();
        const Int rank = g.VRRank();

        // Perform the SendRecv to make A have the same colAlignment
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int colShift = this->ColShift();

        const Int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const Int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLength(height,p);
        const Int portionSize = mpi::Pad( maxLocalHeightOfA*width );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuffer = A.LockedBuffer();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ABuffer[j*ALDim];
            T* secondBufCol = &secondBuf[j*localHeightOfA];
            MemCopy( secondBufCol, ACol, localHeightOfA );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuf, portionSize, sendRank, 0,
          firstBuf,  portionSize, recvRank, mpi::ANY_TAG, g.VRComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShiftOfA = Shift_( col+c*k, colAlignment, p );
            const Int colOffset = (colShiftOfA-colShift) / c;
            const Int localHeight = Length_( height, colShiftOfA, p );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for 
#endif
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &thisBuffer[colOffset+j*thisLDim];
                const T* sourceCol = &data[j*localHeight];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    destCol[iLoc*r] = sourceCol[iLoc];
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [* ,VR]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC,Int>(g) );
    *A_STAR_VC = A;

    std::auto_ptr<DistMatrix<T,MR,MC,Int> > A_MR_MC
    ( new DistMatrix<T,MR,MC,Int>(true,false,this->ColAlignment(),0,g) );
    *A_MR_MC = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    *this = *A_MR_MC;
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    if( this->Viewing() )
        this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Int c = g.Width();
    const Int colShift = this->ColShift();

    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();

    T* thisBuffer = this->Buffer();
    const Int thisLDim = this->LDim();
    const T* ABuffer = A.LockedBuffer();
    const Int ALDim = A.LDim();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
    {
        T* destCol = &thisBuffer[j*thisLDim];
        const T* sourceCol = &ABuffer[colShift+j*ALDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            destCol[iLoc] = sourceCol[iLoc*c];
    }
    return *this;
}

template<typename T,typename Int>
const DistMatrix<T,MR,STAR,Int>&
DistMatrix<T,MR,STAR,Int>::operator=( const DistMatrix<T,CIRC,CIRC,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ] = [o ,o ]");
#endif
    DistMatrix<T,MR,MC> A_MR_MC( A.Grid() );
    A_MR_MC.AlignWith( *this );
    A_MR_MC = A;
    *this = A_MR_MC;
    return *this;
}

//
// Routines which explicitly work in the complex plane
//

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    if( g.Col() == ownerCol )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        this->SetLocalRealPart( iLoc, j, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    if( g.Col() == ownerCol )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        this->SetLocalImagPart( iLoc, j, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    if( g.Col() == ownerCol )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        this->UpdateLocalRealPart( iLoc, j, u );
    }
}

template<typename T,typename Int>
void
DistMatrix<T,MR,STAR,Int>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[MR,* ]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    if( g.Col() == ownerCol )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        this->UpdateLocalImagPart( iLoc, j, u );
    }
}

#define PROTO(T) \
  template class DistMatrix<T,MR,STAR,int>
#define COPY(T,CD,RD) \
  template DistMatrix<T,MR,STAR,int>::DistMatrix( \
    const DistMatrix<T,CD,RD,int>& A )
#define FULL(T) \
  PROTO(T); \
  COPY(T,CIRC,CIRC); \
  COPY(T,MC,  MR); \
  COPY(T,MC,  STAR); \
  COPY(T,MD,  STAR); \
  COPY(T,MR,  MC  ); \
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
