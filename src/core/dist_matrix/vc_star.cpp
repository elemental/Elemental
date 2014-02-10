/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

template<typename T>
using GDM = GeneralDistMatrix<T,VC,STAR>;
template<typename T>
using DM = DistMatrix<T,VC,STAR>;

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
DM<T>::DistMatrix( const elem::Grid& g )
: GDM<T>(g)
{ this->SetShifts(); }

template<typename T>
DM<T>::DistMatrix( Int height, Int width, const elem::Grid& g )
: GDM<T>(g)
{ this->SetShifts(); this->Resize(height,width); }

template<typename T>
DM<T>::DistMatrix( Int height, Int width, Int colAlign, const elem::Grid& g )
: GDM<T>(g)
{ this->Align(colAlign,0); this->Resize(height,width); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, Int ldim, const elem::Grid& g )
: GDM<T>(g)
{ this->Align(colAlign,0); this->Resize(height,width,ldim); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, const T* buffer, Int ldim,
  const elem::Grid& g )
: GDM<T>(g)
{ this->LockedAttach(height,width,colAlign,0,buffer,ldim,g); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, T* buffer, Int ldim,
  const elem::Grid& g )
: GDM<T>(g)
{ this->Attach(height,width,colAlign,0,buffer,ldim,g); }

template<typename T>
DM<T>::DistMatrix( const DM<T>& A )
: GDM<T>(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[VC,* ]::DistMatrix"))
    this->SetShifts();
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct [VC,* ] with itself");
}

template<typename T>
template<Dist U,Dist V>
DM<T>::DistMatrix( const DistMatrix<T,U,V>& A )
: GDM<T>(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[VC,* ]::DistMatrix"))
    this->SetShifts();
    if( VC != U || STAR != V || 
        reinterpret_cast<const DM<T>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [VC,* ] with itself");
}

template<typename T>
DM<T>::DistMatrix( DM<T>&& A ) : GDM<T>(std::move(A)) { }

template<typename T> DM<T>::~DistMatrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[VC,* ] = [MC,MR]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
    this->AlignColsAndResize( A.ColAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlign() % g.Height() == A.ColAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int colShiftOfA = A.ColShift();
        const Int colAlign = this->ColAlign();
        const Int rowAlignOfA = A.RowAlign();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLength(height,p);
        const Int maxWidth = MaxLength(width,c);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*c*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuf[k*portionSize];
            const Int thisRank = row+k*r;
            const Int thisColShift = Shift_(thisRank,colAlign,p); 
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuf, portionSize, recvBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int thisRowShift = Shift_(k,rowAlignOfA,c);
            const Int thisLocalWidth = Length_(width,thisRowShift,c);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuf[(thisRowShift+jLoc*c)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [VC,* ] <- [MC,MR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int colShiftOfA = A.ColShift();
        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();
        const Int rowAlignOfA = A.RowAlign();
        
        const Int sendRow = (row+r+(colAlign%r)-colAlignOfA) % r;
        const Int recvRow = (row+r+colAlignOfA-(colAlign%r)) % r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLength(height,p);
        const Int maxWidth = MaxLength(width,c);
        const Int portionSize = mpi::Pad( maxHeight*maxWidth );

        T* buffer = this->auxMemory_.Require( 2*c*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[c*portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &secondBuf[k*portionSize];
            const Int thisRank = sendRow+k*r;
            const Int thisColShift = Shift_(thisRank,colAlign,p);
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = Length_(height,thisColShift,p);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                T* destCol = &data[jLoc*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+jLoc*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // AllToAll to gather all of the unaligned [VC,*] data into firstBuf
        mpi::AllToAll
        ( secondBuf, portionSize, firstBuf, portionSize, g.RowComm() );

        // SendRecv: properly align the [VC,*] via a trade in the column
        mpi::SendRecv
        ( firstBuf,  portionSize, sendRow,
          secondBuf, portionSize, recvRow, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int thisRowShift = Shift_(k,rowAlignOfA,c);
            const Int thisLocalWidth = Length_(width,thisRowShift,c);
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*localHeight];
                T* thisCol = &thisBuf[(thisRowShift+jLoc*c)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[VC,* ] = [MC,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
    this->AlignColsAndResize( A.ColAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlign() % g.Height() == A.ColAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int colShift = this->ColShift();
        const Int colShiftOfA = A.ColShift();
        const Int colOffset = (colShift-colShiftOfA) / r;
        
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
                destCol[iLoc] = sourceCol[iLoc*c];
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [VC,* ] <- [MC,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int col = g.Col();
        const Int colShiftOfA = A.ColShift();
        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();

        // We will SendRecv A[VC,*] within our process column to fix alignments.
        const Int sendRow = (row+r+(colAlign%r)-colAlignOfA) % r;
        const Int recvRow = (row+r+colAlignOfA-(colAlign%r)) % r;
        const Int sendRank = sendRow + r*col;

        const Int sendColShift = Shift( sendRank, colAlign, p );
        const Int sendColOffset = (sendColShift-colShiftOfA) / r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localHeightOfSend = Length(height,sendColShift,p);

        const Int sendSize = localHeightOfSend * width;
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
            T* destCol = &sendBuf[j*localHeightOfSend];
            const T* sourceCol = &ABuf[sendColOffset+j*ALDim];
            for( Int iLoc=0; iLoc<localHeightOfSend; ++iLoc )
                destCol[iLoc] = sourceCol[iLoc*c];
        }

        // Communicate
        mpi::SendRecv
        ( sendBuf, sendSize, sendRow,
          recvBuf, recvSize, recvRow, g.ColComm() );

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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,* ] = [* ,MR]"))
    DistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[VC,* ] = [MD,* ]"))
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,* ] = [* ,MD]"))
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,* ] = [MR,MC]"))
    DistMatrix<T,VR,STAR> A_VR_STAR( A );
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,* ] = [MR,* ]"))
    DistMatrix<T,VR,STAR> A_VR_STAR( A );
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,* ] = [* ,MC]"))
    std::unique_ptr<DistMatrix<T,MR,MC>> A_MR_MC( new DistMatrix<T,MR,MC>(A) );
    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(*A_MR_MC) );
    delete A_MR_MC.release(); // lowers memory highwater
    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DM<T>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[VC,* ] = [VC,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const Grid& g = this->Grid();
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
            std::cerr << "Unaligned [VC,* ] <- [VC,* ]." << std::endl;
#endif
        const Int rank = g.VCRank();
        const Int p = g.Size();

        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();

        const Int sendRank = (rank+p+colAlign-colAlignOfA) % p;
        const Int recvRank = (rank+p+colAlignOfA-colAlign) % p;

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
          recvBuf, recvSize, recvRank, g.VCComm() );

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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,* ] = [* ,VC]"))
    std::unique_ptr<DistMatrix<T,MR,MC>> A_MR_MC( new DistMatrix<T,MR,MC>(A) );
    std::unique_ptr<DM<T>> A_VR_STAR
    ( new DM<T>(*A_MR_MC) );
    delete A_MR_MC.release(); // lowers memory highwater
    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[VC,* ] = [VR,* ]");
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
    PARALLEL_FOR
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[VC,* ] = [* ,VR]"))
    DistMatrix<T,MC,MR> A_MC_MR( A );
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[VC,* ] = [* ,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const Grid& g = this->Grid();
    this->Resize( A.Height(), A.Width() );
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[VC,* ] = [o ,o ]");
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

template<typename T>
DM<T>&
DM<T>::operator=( DM<T>&& A )
{
    GDM<T>::operator=( std::move(A) );
    return *this;
}

// Realignment
// -----------
template<typename T>
void
DM<T>::AlignWith( const elem::DistData& data )
{
    DEBUG_ONLY(CallStackEntry cse("[VC,* ]::AlignWith"))
    this->SetGrid( *data.grid );
    
    if( data.colDist == MC || data.colDist == VC )
        this->AlignCols( data.colAlign );
    else if( data.rowDist == MC || data.rowDist == VC )
        this->AlignCols( data.rowAlign );
    DEBUG_ONLY(else LogicError("Nonsensical alignment"))
}

template<typename T>
void
DM<T>::AlignColsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

// Specialized redistributions
// ---------------------------

template<typename T>
void
DM<T>::SumScatterFrom( const DistMatrix<T,MC,STAR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[VC,* ]::SumScatterFrom( [MC,* ] )");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && A.Grid().Rank() == 0 )
    {
        std::cerr <<
          "[VC,* ]::SumScatterFrom([MC,* ]) potentially causes a large amount "
          "of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MC,* ] matrix instead." << std::endl;
    }
#endif
    this->AlignColsAndResize( A.ColAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return;

    if( this->ColAlign() % g.Height() == A.ColAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int colAlign = this->ColAlign();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int maxLocalHeight = MaxLength( height, p );
        const Int recvSize = mpi::Pad( maxLocalHeight*width );
        const Int sendSize = c*recvSize;

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        T* buffer = this->auxMemory_.Require( sendSize );
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = myRow+k*r;
            const Int thisColShift = Shift_( thisRank, colAlign, p );
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = Length_( height, thisColShift, p );
            INNER_PARALLEL_FOR
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &data[j*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+j*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.RowComm() );

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
        ("Unaligned [VC,* ]::ReduceScatterFrom( [MC,* ] ) not implemented");
    }
}

template<typename T>
void
DM<T>::SumScatterFrom( const DistMatrix<T,STAR,STAR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[VC,* ]::SumScatterFrom( [* ,* ] )");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
    this->Resize( A.Height(), A.Width() );
    if( !this->Participating() )
        return;

    const Int p = g.Size();
    const Int colAlign = this->ColAlign();

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
        const Int thisColShift = Shift_( k, colAlign, p );
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
    mpi::ReduceScatter( buffer, recvSize, g.VCComm() );

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
DM<T>::SumScatterUpdate( T alpha, const DistMatrix<T,MC,STAR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[VC,* ]::SumScatterUpdate( [MC,* ] )");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && A.Grid().Rank() == 0 )
    {
        std::cerr <<
          "[VC,* ]::SumScatterUpdate([MC,* ]) potentially causes a large amount"
          " of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MC,* ] matrix instead." << std::endl;
    }
#endif
    if( this->ColAlign() % g.Height() == A.ColAlign() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int colAlign = this->ColAlign();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int maxLocalHeight = MaxLength( height, p );
        const Int recvSize = mpi::Pad( maxLocalHeight*width );
        const Int sendSize = c*recvSize;

        T* buffer = this->auxMemory_.Require( sendSize );

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[k*recvSize];
            const Int thisRank = myRow+k*r;
            const Int thisColShift = Shift_( thisRank, colAlign, p );
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = Length_( height, thisColShift, p );
            INNER_PARALLEL_FOR
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &data[j*thisLocalHeight];
                const T* sourceCol = &ABuf[thisColOffset+j*ALDim];
                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    destCol[iLoc] = sourceCol[iLoc*c];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, g.RowComm() );

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
        ("Unaligned [VC,* ]::ReduceScatterUpdate( [MC,* ] ) not implemented");
    }
}

template<typename T>
void
DM<T>::SumScatterUpdate( T alpha, const DistMatrix<T,STAR,STAR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[VC,* ]::SumScatterUpdate( [* ,* ] )");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
    const elem::Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int p = g.Size();
    const Int colAlign = this->ColAlign();

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
        const Int thisColShift = Shift_( k, colAlign, p );
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
    mpi::ReduceScatter( buffer, recvSize, g.VCComm() );

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

// Basic queries
// =============

template<typename T>
elem::DistData DM<T>::DistData() const { return elem::DistData(*this); }

template<typename T>
mpi::Comm DM<T>::DistComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm DM<T>::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM<T>::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM<T>::ColComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm DM<T>::RowComm() const { return mpi::COMM_SELF; }

template<typename T>
Int DM<T>::ColStride() const { return this->grid_->Size(); }
template<typename T>
Int DM<T>::RowStride() const { return 1; }

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define PROTO(T) template class DistMatrix<T,VC,STAR>
#define COPY(T,U,V) \
  template DistMatrix<T,VC,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
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
