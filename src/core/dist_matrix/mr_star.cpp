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
DistMatrix<T,MR,STAR>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->SetShifts(); } 

template<typename T>
DistMatrix<T,MR,STAR>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T>
DistMatrix<T,MR,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Align(colAlignment,0); this->ResizeTo(height,width); }

template<typename T>
DistMatrix<T,MR,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Align(colAlignment,0); this->ResizeTo(height,width,ldim); }

template<typename T>
DistMatrix<T,MR,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->LockedAttach( height, width, colAlignment, buffer, ldim, g ); }

template<typename T>
DistMatrix<T,MR,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Attach( height, width, colAlignment, buffer, ldim, g ); }

template<typename T>
DistMatrix<T,MR,STAR>::DistMatrix( const DistMatrix<T,MR,STAR>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[MR,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [MR,* ] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DistMatrix<T,MR,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[MR,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( MR != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,MR,STAR>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [MR,* ] with itself");
}

template<typename T>
DistMatrix<T,MR,STAR>::DistMatrix( DistMatrix<T,MR,STAR>&& A )
: AbstractDistMatrix<T>(std::move(A))
{ }

template<typename T>
DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( DistMatrix<T,MR,STAR>&& A )
{
    AbstractDistMatrix<T>::operator=( std::move(A) );
    return *this;
}

template<typename T>
DistMatrix<T,MR,STAR>::~DistMatrix()
{ }

template<typename T>
elem::DistData
DistMatrix<T,MR,STAR>::DistData() const
{
    elem::DistData data;
    data.colDist = MR;
    data.rowDist = STAR;
    data.colAlignment = this->colAlignment_;
    data.rowAlignment = 0;
    data.root = 0;
    data.diagPath = 0;
    data.grid = this->grid_;
    return data;
}

template<typename T>
Int
DistMatrix<T,MR,STAR>::ColStride() const
{ return this->grid_->Width(); }

template<typename T>
Int
DistMatrix<T,MR,STAR>::RowStride() const
{ return 1; }

template<typename T>
Int
DistMatrix<T,MR,STAR>::ColRank() const
{ return this->grid_->Col(); }

template<typename T>
Int
DistMatrix<T,MR,STAR>::RowRank() const
{ return 0; }

template<typename T>
void
DistMatrix<T,MR,STAR>::AlignWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::AlignWith");
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
    else LogicError("Nonsensical alignment");
#endif
    this->constrainedColAlignment_ = true;
    this->SetShifts();
}

template<typename T>
void
DistMatrix<T,MR,STAR>::AlignWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,MR,STAR>::AlignColsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

template<typename T>
void
DistMatrix<T,MR,STAR>::AlignColsWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,MR,STAR>::Attach
( Int height, Int width, Int colAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::Attach");
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

template<typename T>
void
DistMatrix<T,MR,STAR>::LockedAttach
( Int height, Int width, Int colAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::LockedAttach");
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

template<typename T>
void
DistMatrix<T,MR,STAR>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),this->Grid().Width()), width );
}

template<typename T>
void
DistMatrix<T,MR,STAR>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),this->Grid().Width()), width, ldim );
}

template<typename T>
T
DistMatrix<T,MR,STAR>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::Get");
    this->AssertValidEntry( i, j );
    if( !this->Participating() )
        LogicError("Should only be called by grid members");
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
    mpi::Broadcast( u, ownerCol, g.RowComm() );
    return u;
}

template<typename T>
void
DistMatrix<T,MR,STAR>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::Set");
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

template<typename T>
void
DistMatrix<T,MR,STAR>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::SetRealPart");
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

template<typename T>
void
DistMatrix<T,MR,STAR>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::SetImagPart");
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

template<typename T>
void
DistMatrix<T,MR,STAR>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::Update");
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

template<typename T>
void
DistMatrix<T,MR,STAR>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::UpdateRealPart");
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

template<typename T>
void
DistMatrix<T,MR,STAR>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::UpdateImagPart");
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

//
// Utility functions, e.g., SumOverCol
//

template<typename T>
void
DistMatrix<T,MR,STAR>::SumOverCol()
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::SumOverCol");
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
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        const T* thisCol = &thisBuffer[j*thisLDim];
        T* sendBufCol = &sendBuf[j*localHeight];
        MemCopy( sendBufCol, thisCol, localHeight );
    }

    // AllReduce sum
    mpi::AllReduce( sendBuf, recvBuf, localSize, g.ColComm() );

    // Unpack
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        const T* recvBufCol = &recvBuf[j*localHeight];
        T* thisCol = &thisBuffer[j*thisLDim];
        MemCopy( thisCol, recvBufCol, localHeight );
    }
    this->auxMemory_.Release();
}

template<typename T>
void
DistMatrix<T,MR,STAR>::AdjointFrom( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::AdjointFrom");
#endif
    this->TransposeFrom( A, true );
}

template<typename T>
void
DistMatrix<T,MR,STAR>::TransposeFrom
( const DistMatrix<T,MC,MR>& A, bool conjugate )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::TransposeFrom");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetColAlignmentAndResize( A.RowAlignment(), A.Width(), A.Height() );
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
            PARALLEL_FOR
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
            PARALLEL_FOR
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
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int rowShift = Shift_( k, colAlignmentOfA, r );
            const Int localWidth = Length_( width, rowShift, r );
            INNER_PARALLEL_FOR
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
            PARALLEL_FOR
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
            PARALLEL_FOR
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
        ( secondBuf, portionSize, sendCol,
          firstBuf,  portionSize, recvCol, g.RowComm() );

        // Use the output of the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.ColComm() );

        // Unpack the contents of each member of the process col
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int colAlignmentOfA = A.ColAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShift = Shift_( k, colAlignmentOfA, r );
            const Int localWidth = Length_( width, rowShift, r );
            INNER_PARALLEL_FOR
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

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [MC,MR]");
#endif
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = A;

    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(true,this->ColAlignment(),g) );
    *A_VR_STAR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [MC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
    {
        this->ResizeTo( A.Height(), A.Width() );
        return *this;
    }

    if( A.Width() == 1 )
    {
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
            PARALLEL_FOR
            for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                sendBuf[iLoc] = ABuffer[offset+iLoc*c];
        }

        // A[VR,* ] <- A[VC,* ]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankRM,
          recvBuf, portionSize, recvRankRM, g.VRComm() );

        // A[MR,* ] <- A[VR,* ]
        mpi::AllGather
        ( recvBuf, portionSize,
          sendBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        PARALLEL_FOR
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
        std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
        ( new DistMatrix<T,VC,STAR>(g) );
        *A_VC_STAR = A;

        std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
        ( new DistMatrix<T,VR,STAR>(true,this->ColAlignment(),g) );
        *A_VR_STAR = *A_VC_STAR;
        delete A_VC_STAR.release(); // lowers memory highwater

        *this = *A_VR_STAR;
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [* ,MR]");
#endif
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR( new DistMatrix<T,MC,MR>(A) );

    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(true,this->ColAlignment(),g) );
    *A_VR_STAR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [MD,* ]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [* ,MD]");
#endif
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetColAlignmentAndResize( A.ColAlignment(), A.Height(), A.Width() );
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
        PARALLEL_FOR
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
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int rowShift = Shift_( k, rowAlignmentOfA, r );
            const Int localWidth = Length_( width, rowShift, r );
            INNER_PARALLEL_FOR
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
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuffer[jLoc*ALDim];
            T* secondBufCol = &secondBuf[jLoc*localHeightOfA];
            MemCopy( secondBufCol, ACol, localHeightOfA );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuf, portionSize, sendCol,
          firstBuf,  portionSize, recvCol, g.RowComm() );

        // Use the output of the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.ColComm() );

        // Unpack the contents of each member of the process col
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int rowAlignmentOfA = A.RowAlignment();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShift = Shift_( k, rowAlignmentOfA, r );
            const Int localWidth = Length_( width, rowShift, r );
            INNER_PARALLEL_FOR
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

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetColAlignmentAndResize( A.ColAlignment(), A.Height(), A.Width() );
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
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ABuffer[j*ALDim];
            T* sendBufCol = &sendBuf[j*localHeightOfA];
            MemCopy( sendBufCol, ACol, localHeightOfA );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRank,
          recvBuf, recvSize, recvRank, g.RowComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
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

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [* ,MC]");
#endif
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [VC,* ]");
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VR,STAR> A_VR_STAR(true,this->ColAlignment(),g);

    A_VR_STAR = A;
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [* ,VC]");
#endif
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [VR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
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
    this->SetColAlignmentAndResize
    ( A.ColAlignment()%g.Width(), A.Height(), A.Width() );
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
        PARALLEL_FOR
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
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int colShiftOfA = Shift_( col+c*k, colAlignmentOfA, p );
            const Int colOffset = (colShiftOfA-colShift) / c;
            const Int localHeight = Length_( height, colShiftOfA, p );
            INNER_PARALLEL_FOR
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
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ABuffer[j*ALDim];
            T* secondBufCol = &secondBuf[j*localHeightOfA];
            MemCopy( secondBufCol, ACol, localHeightOfA );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuf, portionSize, sendRank,
          firstBuf,  portionSize, recvRank, g.VRComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuffer = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShiftOfA = Shift_( col+c*k, colAlignment, p );
            const Int colOffset = (colShiftOfA-colShift) / c;
            const Int localHeight = Length_( height, colShiftOfA, p );
            INNER_PARALLEL_FOR
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

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [* ,VR]");
#endif
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(A) );

    std::unique_ptr<DistMatrix<T,MR,MC>> A_MR_MC
    ( new DistMatrix<T,MR,MC>(true,false,this->ColAlignment(),0,g) );
    *A_MR_MC = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    *this = *A_MR_MC;
    return *this;
}

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );

    const Int c = g.Width();
    const Int colShift = this->ColShift();

    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();

    T* thisBuffer = this->Buffer();
    const Int thisLDim = this->LDim();
    const T* ABuffer = A.LockedBuffer();
    const Int ALDim = A.LDim();
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        T* destCol = &thisBuffer[j*thisLDim];
        const T* sourceCol = &ABuffer[colShift+j*ALDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            destCol[iLoc] = sourceCol[iLoc*c];
    }
    return *this;
}

template<typename T>
const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [o ,o ]");
#endif
    DistMatrix<T,MR,MC> A_MR_MC( A.Grid() );
    A_MR_MC.AlignWith( *this );
    A_MR_MC = A;
    *this = A_MR_MC;
    return *this;
}

#define PROTO(T) template class DistMatrix<T,MR,STAR>
#define COPY(T,CD,RD) \
  template DistMatrix<T,MR,STAR>::DistMatrix( const DistMatrix<T,CD,RD>& A )
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
