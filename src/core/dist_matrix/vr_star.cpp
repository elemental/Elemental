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
DistMatrix<T,VR,STAR>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->SetShifts(); }

template<typename T>
DistMatrix<T,VR,STAR>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T>
DistMatrix<T,VR,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Align(colAlignment,0); this->ResizeTo(height,width); }

template<typename T>
DistMatrix<T,VR,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Align(colAlignment,0); this->ResizeTo(height,width,ldim); }

template<typename T>
DistMatrix<T,VR,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->LockedAttach(height,width,colAlignment,buffer,ldim,g); }

template<typename T>
DistMatrix<T,VR,STAR>::DistMatrix
( Int height, Int width, Int colAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T>(g)
{ this->Attach(height,width,colAlignment,buffer,ldim,g); }

template<typename T>
DistMatrix<T,VR,STAR>::DistMatrix( const DistMatrix<T,VR,STAR>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[VR,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [VR,* ] with itself");
}

template<typename T>
template<Distribution U,Distribution V>
DistMatrix<T,VR,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix[VR,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( VR != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,VR,STAR>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [VR,* ] with itself");
}

template<typename T>
DistMatrix<T,VR,STAR>::~DistMatrix()
{ }

template<typename T>
elem::DistData
DistMatrix<T,VR,STAR>::DistData() const
{
    elem::DistData data;
    data.colDist = VR;
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
DistMatrix<T,VR,STAR>::ColStride() const
{ return this->grid_->Size(); }

template<typename T>
Int
DistMatrix<T,VR,STAR>::RowStride() const
{ return 1; }

template<typename T>
Int
DistMatrix<T,VR,STAR>::ColRank() const
{ return this->grid_->VRRank(); }

template<typename T>
Int
DistMatrix<T,VR,STAR>::RowRank() const
{ return 0; }

template<typename T>
void
DistMatrix<T,VR,STAR>::AlignWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::AlignWith");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( data.colDist == MR || data.colDist == VR )
        this->colAlignment_ = data.colAlignment;
    else if( data.rowDist == MR || data.rowDist == VR )
        this->colAlignment_ = data.rowAlignment;
#ifndef RELEASE
    else LogicError("Nonsensical alignment");
#endif
    this->constrainedColAlignment_ = true;
    this->SetShifts();
}

template<typename T>
void
DistMatrix<T,VR,STAR>::AlignWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,VR,STAR>::AlignColsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

template<typename T>
void
DistMatrix<T,VR,STAR>::AlignColsWith( const AbstractDistMatrix<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DistMatrix<T,VR,STAR>::Attach
( Int height, Int width, Int colAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::Attach");
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
        const Int localHeight = Length(height,this->colShift_,g.Size());
        this->matrix_.Attach_( localHeight, width, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,VR,STAR>::LockedAttach
( Int height, Int width, Int colAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::LockedAttach");
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
        const Int localHeight = Length(height,this->colShift_,g.Size());
        this->matrix_.LockedAttach_( localHeight, width, buffer, ldim );
    }
}

template<typename T>
void
DistMatrix<T,VR,STAR>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    const elem::Grid& g = this->Grid();
    this->height_ = height;
    this->width_  = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),g.Size()) ,width );
}

template<typename T>
void
DistMatrix<T,VR,STAR>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::ResizeTo");
    this->AssertNotLocked();
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
#endif
    const elem::Grid& g = this->Grid();
    this->height_ = height;
    this->width_  = width;
    if( this->Participating() )
        this->matrix_.ResizeTo_
        ( Length(height,this->ColShift(),g.Size()), width, ldim );
}

template<typename T>
T
DistMatrix<T,VR,STAR>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::Get");
    this->AssertValidEntry( i, j );
    // TODO: Generalize this function to always work...
    if( !this->Participating() )
        LogicError("Should only call with processes in grid");
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    T u;
    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        u = this->GetLocal(iLoc,j);
    }
    mpi::Broadcast( u, ownerRank, g.VRComm() );
    return u;
}

template<typename T>
void
DistMatrix<T,VR,STAR>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->SetLocal(iLoc,j,u);
    }
}

template<typename T>
void
DistMatrix<T,VR,STAR>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();
    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->UpdateLocal(iLoc,j,u);
    }
}

//
// Utility functions, e.g., operator=
//

template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [MC,MR]");
#endif
    DistMatrix<T,VC,STAR> A_VC_STAR( A );
    *this = A_VC_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [MC,* ]");
#endif
    DistMatrix<T,VC,STAR> A_VC_STAR( A );
    *this = A_VC_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [* ,MR]");
#endif
    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR( new DistMatrix<T,MC,MR>(A) );
    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(*A_MC_MR) );
    delete A_MC_MR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [MD,* ]");
#endif
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [* ,MD]");
#endif
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetColAlignmentAndResize( A.ColAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();
        const Int colShiftOfA = A.ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLength(height,p);
        const Int maxWidth = MaxLength(width,r);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*r*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &sendBuf[k*portionSize];
            const Int thisRank = col+k*c;
            const Int thisColShift = Shift_(thisRank,colAlignment,p);
            const Int thisColOffset = (thisColShift-colShiftOfA) / c;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuf, portionSize, recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisRowShift = Shift_(k,rowAlignmentOfA,r);
            const Int thisLocalWidth = Length_(width,thisRowShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuf[(thisRowShift+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [VR,* ] <- [MR,MC]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();
        const Int colShiftOfA = A.ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendCol = (col+c+(colAlignment%c)-colAlignmentOfA) % c;
        const Int recvCol = (col+c+colAlignmentOfA-(colAlignment%c)) % c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLength(height,p);
        const Int maxWidth = MaxLength(width,r);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*r*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[r*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int thisRank = sendCol+k*c;
            const Int thisColShift = Shift_(thisRank,colAlignment,p);
            const Int thisColOffset = (thisColShift-colShiftOfA) / c;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // AllToAll to gather all of the unaligned [VR,*] data into firstBuf
        mpi::AllToAll
        ( secondBuf, portionSize, firstBuf, portionSize, g.ColComm() );

        // SendRecv: properly align the [VR,*] via a trade in the row
        mpi::SendRecv
        ( firstBuf,  portionSize, sendCol, 
          secondBuf, portionSize, recvCol, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int thisRowShift = Shift_(k,rowAlignmentOfA,r);
            const Int thisLocalWidth = Length_(width,thisRowShift,r);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuf[(thisRowShift+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->SetColAlignmentAndResize( A.ColAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int colShift = this->ColShift();
        const Int colShiftOfA = A.ColShift();
        const Int colOffset = (colShift-colShiftOfA) / c;

        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();

        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisBuf[j*thisLDim];
            const T* sourceCol = &ABuf[colOffset+j*ALDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*r];
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [VR,* ] <- [MR,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int col = g.Col();
        const Int colShiftOfA = A.ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[VR,*] within our process row to fix alignments.
        const Int sendCol = (col+c+(colAlignment%c)-colAlignmentOfA) % c;
        const Int recvCol = (col+c+colAlignmentOfA-(colAlignment%c)) % c;
        const Int sendRank = sendCol + c*row;

        const Int sendColShift = Shift( sendRank, colAlignment, p );
        const Int sendColOffset = (sendColShift-colShiftOfA) / c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localHeightOfSend = Length(height,sendColShift,p);
        const Int maxLocalHeight = MaxLength(height,p);

        const Int portionSize = maxLocalHeight * width;

        T* buffer = this->auxMemory_.Require( 2*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &sendBuf[j*localHeightOfSend];
            const T* sourceCol = &ABuf[sendColOffset+j*ALDim];
            for( Int iLoc=0; iLoc<localHeightOfSend; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*r];
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, portionSize, sendCol, 
          recvBuf, portionSize, recvCol, g.RowComm() );

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
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [* ,MC]");
#endif
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [VC,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int rankCM = g.VCRank();
    const Int rankRM = g.VRRank();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localHeightOfA = A.LocalHeight();
    const Int maxLocalHeight = MaxLength(height,p);

    const Int portionSize = maxLocalHeight * width;

    const Int colShift = this->ColShift();
    const Int colShiftOfA = A.ColShift();

    // Compute which rowmajor rank has the colShift equal to our colShiftOfA
    const Int sendRankRM = (rankRM+(p+colShiftOfA-colShift)) % p;

    // Compute which rowmajor rank has the A colShift that we need
    const Int recvRankCM = (rankCM+(p+colShift-colShiftOfA)) % p;
    const Int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

    T* buffer = this->auxMemory_.Require( 2*portionSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

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
    ( sendBuf, portionSize, sendRankRM, 
      recvBuf, portionSize, recvRankRM, g.VRComm() );

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
    return *this;
}

template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [* ,VC]");
#endif
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [VR,* ]");
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
            std::cerr << "Unaligned [VR,* ] <- [VR,* ]." << std::endl;
#endif
        const Int rank = g.VRRank();
        const Int p = g.Size();

        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const Int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

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
          recvBuf, recvSize, recvRank, g.VRComm() );

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
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [* ,VR]");
#endif
    std::unique_ptr<DistMatrix<T,MC,MR>> A_MC_MR( new DistMatrix<T,MC,MR>(A) );
    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(*A_MC_MR) );
    delete A_MC_MR.release(); // lowers memory highwater
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [* ,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int p = g.Size();
    const Int colShift = this->ColShift();

    const Int localHeight = this->LocalHeight();
    const Int width = this->Width();

    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        T* destCol = &thisBuf[j*thisLDim];
        const T* sourceCol = &ABuf[colShift+j*ALDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            destCol[iLoc] = sourceCol[iLoc*p];
    }
    return *this;
}

// NOTE: This is a small modification of [MC,MR] <- [o ,o ]
template<typename T>
const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ] = [o ,o ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = g.Size();
    this->ResizeTo( m, n );

    const Int colAlignment = this->ColAlignment();
    const Int mLocal = this->LocalHeight();
    const Int pkgSize = mpi::Pad(MaxLength(m,p)*n);
    const Int recvSize = pkgSize;
    const Int sendSize = p*pkgSize;
    T* recvBuf;
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
            const Int q = (colAlignment+s) % p;
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
        ( sendBuf, pkgSize, recvBuf, pkgSize, A.Root(), g.VRComm() );
    }
    else if( this->Participating() )
    {
        recvBuf = this->auxMemory_.Require( recvSize );

        // Perform the receiving portion of the scatter from the non-root
        mpi::Scatter
        ( static_cast<T*>(0), pkgSize, 
          recvBuf,            pkgSize, A.Root(), g.VRComm() );
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

template<typename T>
void
DistMatrix<T,VR,STAR>::SumScatterFrom( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::SumScatterFrom( [MR,* ] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.Rank() == 0 )
    {
        std::cerr <<
          "[VR,* ]::SumScatterFrom([MR,* ]) potentially causes a large amount "
          "of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MR,* ] matrix instead." << std::endl;
    }
#endif
    this->SetColAlignmentAndResize( A.ColAlignment(), A.Height(), A.Width() );
    if( !this->Participating() )
        return;

    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myCol = g.Col();
        const Int colAlignment = this->ColAlignment();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int maxLocalHeight = MaxLength( height, p );
        const Int recvSize = mpi::Pad( maxLocalHeight*width );
        const Int sendSize = r*recvSize;

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        T* buffer = this->auxMemory_.Require( sendSize );
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = myCol+k*c;
            const Int thisColShift = Shift_( thisRank, colAlignment, p );
            const Int thisColOffset = (thisColShift-colShiftOfA) / c;
            const Int thisLocalHeight = Length_( height, thisColShift, p );
            INNER_PARALLEL_FOR
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &data[j*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+j*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.ColComm() );

        // Unpack our received data
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* bufferCol = &buffer[j*localHeight];
            T* thisCol = &thisBuf[j*thisLDim];
            MemCopy( thisCol, bufferCol, localHeight );
        }
        this->auxMemory_.Release();
    }
    else
    {
        LogicError
        ("Unaligned [VR,* ]::ReduceScatterFrom( [MR,* ] ) not implemented");
    }
}

template<typename T>
void
DistMatrix<T,VR,STAR>::SumScatterFrom( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::SumScatterFrom( [* ,* ] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->ResizeTo( A.Height(), A.Width() );
    if( !this->Participating() )
        return;

    const Int p = g.Size();
    const Int colAlignment = this->ColAlignment();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int maxLocalHeight = MaxLength( height, p );

    const Int recvSize = mpi::Pad( maxLocalHeight*width );
    const Int sendSize = p*recvSize;

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    T* buffer = this->auxMemory_.Require( sendSize );
    OUTER_PARALLEL_FOR
    for( Int k=0; k<p; ++k )
    {
        T* data = &buffer[k*recvSize];
        const Int thisColShift = Shift_( k, colAlignment, p );
        const Int thisLocalHeight = Length_( height, thisColShift, p );
        INNER_PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &data[j*thisLocalHeight];
            const T* sourceCol = &ABuf[thisColShift+j*ALDim];
            for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*p];
        }
    }

    // Communicate
    mpi::ReduceScatter( buffer, recvSize, g.VRComm() );

    // Unpack our received data
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        const T* bufferCol = &buffer[j*localHeight];
        T* thisCol = &thisBuf[j*thisLDim];
        MemCopy( thisCol, bufferCol, localHeight );
    }
    this->auxMemory_.Release();
}

template<typename T>
void
DistMatrix<T,VR,STAR>::SumScatterUpdate
( T alpha, const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::SumScatterUpdate( [MR,* ] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.Rank() == 0 )
    {
        std::cerr <<
          "[VR,* ]::SumScatterUpdate([MR,* ]) potentially causes a large amount"
          " of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MR,* ] matrix instead." << std::endl;
    }
#endif
    if( !this->Participating() )
        return;

    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myCol = g.Col();
        const Int colAlignment = this->ColAlignment();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int maxLocalHeight = MaxLength( height, p );

        const Int recvSize = mpi::Pad( maxLocalHeight*width );
        const Int sendSize = r*recvSize;

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        T* buffer = this->auxMemory_.Require( sendSize );
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = myCol+k*c;
            const Int thisColShift = Shift_( thisRank, colAlignment, p );
            const Int thisColOffset = (thisColShift-colShiftOfA) / c;
            const Int thisLocalHeight = Length_( height, thisColShift, p );
            INNER_PARALLEL_FOR
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &data[j*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+j*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*r];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.ColComm() );

        // Unpack our received data
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const T* bufferCol = &buffer[j*localHeight];
            T* thisCol = &thisBuf[j*thisLDim];
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                thisCol[iLoc] += alpha*bufferCol[iLoc];
        }
        this->auxMemory_.Release();
    }
    else
    {
        LogicError
        ("Unaligned [VR,* ]::ReduceScatterUpdate( [MR,* ] ) not implemented");
    }
}

template<typename T>
void
DistMatrix<T,VR,STAR>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::SumScatterUpdate( [* ,* ] )");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
    this->AssertSameSize( A.Height(), A.Width() );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int p = g.Size();
    const Int colAlignment = this->ColAlignment();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int maxLocalHeight = MaxLength( height, p );

    const Int recvSize = mpi::Pad( maxLocalHeight*width );
    const Int sendSize = p*recvSize;

    // Pack
    const Int ALDim = A.LDim();
    const T* ABuf = A.LockedBuffer();
    T* buffer = this->auxMemory_.Require( sendSize );
    OUTER_PARALLEL_FOR
    for( Int k=0; k<p; ++k )
    {
        T* data = &buffer[k*recvSize];
        const Int thisColShift = Shift_( k, colAlignment, p );
        const Int thisLocalHeight = Length_( height, thisColShift, p );
        INNER_PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &data[j*thisLocalHeight];
            const T* sourceCol = &ABuf[thisColShift+j*ALDim];
            for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*p];
        }
    }

    // Communicate
    mpi::ReduceScatter( buffer, recvSize, g.VRComm() );

    // Unpack our received data
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int j=0; j<width; ++j )
    {
        const T* bufferCol = &buffer[j*localHeight];
        T* thisCol = &thisBuf[j*thisLDim];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            thisCol[iLoc] += alpha*bufferCol[iLoc];
    }
    this->auxMemory_.Release();
}

//
// Routines which explicitly work in the complex plane
//

template<typename T>
void
DistMatrix<T,VR,STAR>::SetRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();
    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->SetLocalRealPart(iLoc,j,u);
    }
}

template<typename T>
void
DistMatrix<T,VR,STAR>::SetImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();
    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->SetLocalImagPart(iLoc,j,u);
    }
}

template<typename T>
void
DistMatrix<T,VR,STAR>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();
    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->UpdateLocalRealPart(iLoc,j,u);
    }
}

template<typename T>
void
DistMatrix<T,VR,STAR>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
#ifndef RELEASE
    CallStackEntry entry("[VR,* ]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    this->ComplainIfReal();
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();
    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->UpdateLocalImagPart(iLoc,j,u);
    }
}

#define PROTO(T) \
  template class DistMatrix<T,VR,STAR>
#define COPY(T,CD,RD) \
  template DistMatrix<T,VR,STAR>::DistMatrix( const DistMatrix<T,CD,RD>& A )
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
  COPY(T,STAR,MR  ); \
  COPY(T,STAR,STAR); \
  COPY(T,STAR,VC  ); \
  COPY(T,STAR,VR  ); \
  COPY(T,VC,  STAR); 

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
