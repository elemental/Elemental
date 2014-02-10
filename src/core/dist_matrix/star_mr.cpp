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
using GDM = GeneralDistMatrix<T,STAR,MR>;
template<typename T>
using DM = DistMatrix<T,STAR,MR>;

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
DM<T>::DistMatrix
( Int height, Int width, Int rowAlign, const elem::Grid& g )
: GDM<T>(g)
{ this->Align(0,rowAlign); this->Resize(height,width); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int rowAlign, Int ldim, const elem::Grid& g )
: GDM<T>(g)
{ this->Align(0,rowAlign); this->Resize(height,width,ldim); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int rowAlign, const T* buffer, Int ldim, 
  const elem::Grid& g )
: GDM<T>(g)
{ this->LockedAttach(height,width,0,rowAlign,buffer,ldim,g); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int rowAlign, T* buffer, Int ldim, 
  const elem::Grid& g )
: GDM<T>(g)
{ this->Attach(height,width,0,rowAlign,buffer,ldim,g); }

template<typename T>
DM<T>::DistMatrix( const DM<T>& A )
: GDM<T>(A.Grid())
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
DM<T>::DistMatrix( const DistMatrix<T,U,V>& A )
: GDM<T>(A.Grid())
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MR]::DistMatrix"))
    this->SetShifts();
    if( STAR != U || MR != V || 
        reinterpret_cast<const DM<T>*>(&A) != this )   
        *this = A;
    else
        LogicError("Tried to construct [* ,MR] with itself");
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
        CallStackEntry cse("[* ,MR] = [MC,MR]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "The matrix redistribution [* ,MR] <- [MC,MR] potentially causes a "
          "large amount of cache-thrashing. If possible, avoid it by "
          "performing the redistribution with a (conjugate)-transpose: \n"
          << "  [MR,* ].(Conjugate)TransposeFrom([MC,MR])" << std::endl;
    }
#endif
    this->AlignRowsAndResize( A.RowAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlign() == A.RowAlign() )
    {
        if( A.Height() == 1 )
        {
            const Int localWidth = this->LocalWidth();
            T* bcastBuf = this->auxMemory_.Require( localWidth );

            if( g.Row() == A.ColAlign() )
            {
                this->matrix_ = A.LockedMatrix();

                // Pack
                const T* thisBuf = this->LockedBuffer();
                const Int thisLDim = this->LDim();
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                    bcastBuf[jLoc] = thisBuf[jLoc*thisLDim];
            }

            // Communicate
            mpi::Broadcast
            ( bcastBuf, localWidth, A.ColAlign(), g.ColComm() );

            // Unpack
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                thisBuf[jLoc*thisLDim] = bcastBuf[jLoc];

            this->auxMemory_.Release();
        }
        else
        {
            const Int r = g.Height();
            const Int height = this->Height();
            const Int localWidth = this->LocalWidth();
            const Int localHeightOfA = A.LocalHeight();
            const Int maxLocalHeight = MaxLength(height,r);
            const Int portionSize = mpi::Pad( maxLocalHeight*localWidth );

            T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[portionSize];

            // Pack 
            const Int ALDim = A.LDim();
            const T* ABuf = A.LockedBuffer();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* ACol = &ABuf[jLoc*ALDim];
                T* sendBufCol = &sendBuf[jLoc*localHeightOfA];
                MemCopy( sendBufCol, ACol, localHeightOfA );
            }

            // Communicate
            mpi::AllGather
            ( sendBuf, portionSize,
              recvBuf, portionSize, g.ColComm() );

            // Unpack
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
            const Int colAlignOfA = A.ColAlign();
            OUTER_PARALLEL_FOR
            for( Int k=0; k<r; ++k )
            {
                const T* data = &recvBuf[k*portionSize];
                const Int colShift = Shift_( k, colAlignOfA, r );
                const Int localHeight = Length_( height, colShift, r );
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                {
                    T* destCol = &thisBuf[colShift+jLoc*thisLDim];
                    const T* sourceCol = &data[jLoc*localHeight];
                    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                        destCol[iLoc*r] = sourceCol[iLoc];
                }
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MR] <- [MC,MR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int col = g.Col();

        const Int rowAlign = this->RowAlign();
        const Int rowAlignOfA = A.RowAlign();
        const Int sendCol = (col+c+rowAlign-rowAlignOfA) % c;
        const Int recvCol = (col+c+rowAlignOfA-rowAlign) % c;

        if( A.Height() == 1 )
        {
            const Int localWidth = this->LocalWidth();
            T* bcastBuf;

            if( g.Row() == A.ColAlign() )
            {
                const Int localWidthOfA = A.LocalWidth();

                T* buffer = this->auxMemory_.Require(localWidth+localWidthOfA);
                T* sendBuf = &buffer[0];
                bcastBuf   = &buffer[localWidthOfA];

                // Pack
                const T* ABuf = A.LockedBuffer();
                const Int ALDim = A.LDim();
                PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
                    sendBuf[jLoc] = ABuf[jLoc*ALDim];

                // Communicate
                mpi::SendRecv
                ( sendBuf,  localWidthOfA, sendCol, 
                  bcastBuf, localWidth,    recvCol, g.RowComm() );
            }
            else
            {
                bcastBuf = this->auxMemory_.Require( localWidth );
            }

            // Communicate
            mpi::Broadcast
            ( bcastBuf, localWidth, A.ColAlign(), g.ColComm() );

            // Unpack
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                thisBuf[jLoc*thisLDim] = bcastBuf[jLoc];
            this->auxMemory_.Release();
        }
        else
        {
            const Int height = this->Height();
            const Int localWidth  = this->LocalWidth();
            const Int localHeightOfA = A.LocalHeight();
            const Int localWidthOfA  = A.LocalWidth();
            const Int maxLocalHeight = MaxLength(height,r);
            const Int portionSize = mpi::Pad( maxLocalHeight*localWidth );

            T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
            T* firstBuf = &buffer[0];
            T* secondBuf = &buffer[portionSize];

            // Pack
            const Int ALDim = A.LDim();
            const T* ABuf = A.LockedBuffer();
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
            {
                const T* ACol = &ABuf[jLoc*ALDim];
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
            T* thisBuf = this->Buffer();
            const Int thisLDim = this->LDim();
            const Int colAlignOfA = A.ColAlign();
            OUTER_PARALLEL_FOR
            for( Int k=0; k<r; ++k )
            {
                const T* data = &secondBuf[k*portionSize];
                const Int colShift = Shift_( k, colAlignOfA, r );
                const Int localHeight = Length_( height, colShift, r );
                INNER_PARALLEL_FOR
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                {
                    T* destCol = &thisBuf[colShift+jLoc*thisLDim];
                    const T* sourceCol = &data[jLoc*localHeight];
                    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                        destCol[iLoc*r] = sourceCol[iLoc];
                }
            }
            this->auxMemory_.Release();
        }
    }
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,STAR>& A )
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
const DM<T>&
DM<T>::operator=( const DM<T>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [MD,* ]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [* ,MD]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,MC>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,STAR>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MC>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [VC,* ]"))
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(false,true,0,this->RowAlign(),g);
    A_MC_MR = A;
    *this = A_MC_MR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [* ,VC]"))
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR> A_STAR_VR(true,this->RowAlign(),g);
    A_STAR_VR = A;
    *this = A_STAR_VR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VR,STAR>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[* ,MR] = [* ,VR]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
    this->AlignRowsAndResize( A.RowAlign()%g.Width(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->RowAlign() == A.RowAlign() % g.Width() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();

        const Int width = this->Width();
        const Int height = this->Height();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidthOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalWidthOfA );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuf[jLoc*ALDim];
            T* sendBufCol = &sendBuf[jLoc*height];
            MemCopy( sendBufCol, ACol, height );
        }

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int rowShift = this->RowShift();
        const Int rowAlignOfA = A.RowAlign();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int rowShiftOfA = Shift_( col+k*c, rowAlignOfA, p );
            const Int rowOffset = (rowShiftOfA-rowShift) / c;
            const Int localWidth = Length_( width, rowShiftOfA, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MR] <- [* ,VR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();
        const Int rank = g.VRRank();

        // Perform the SendRecv to make A have the same rowAlign
        const Int rowAlign = this->RowAlign();
        const Int rowAlignOfA = A.RowAlign();
        const Int rowShift = this->RowShift();

        const Int sendRank = (rank+p+rowAlign-rowAlignOfA) % p;
        const Int recvRank = (rank+p+rowAlignOfA-rowAlign) % p;

        const Int width = this->Width();
        const Int height = this->Height();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidthOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalWidthOfA );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidthOfA; ++jLoc )
        {
            const T* ACol = &ABuf[jLoc*ALDim];
            T* secondBufCol = &secondBuf[jLoc*height];
            MemCopy( secondBufCol, ACol, height );
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
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShiftOfA = Shift_( col+c*k, rowAlign, p );
            const Int rowOffset = (rowShiftOfA-rowShift) / c;
            const Int localWidth = Length_( width, rowShiftOfA, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("[* ,MR] = [* ,* ]");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    this->Resize( A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    const Int c = this->Grid().Width();
    const Int rowShift = this->RowShift();

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();

    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* ACol = &ABuf[(rowShift+jLoc*c)*ALDim];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        MemCopy( thisCol, ACol, localHeight );
    }
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[* ,MR] = [o ,o ]"))
    DistMatrix<T,MC,MR> A_MC_MR( A );
    A_MC_MR.AlignWith( *this );
    A_MC_MR = A;
    *this = A_MC_MR;
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
DM<T>::AlignRowsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

// Specialized redistributions
// ---------------------------

template<typename T> 
void
DM<T>::SumOverCol()
{
    DEBUG_ONLY(
        CallStackEntry cse("[* ,MR]::SumOverCol");
        this->AssertNotLocked();
    )
    const Grid& g = this->Grid();
    if( !this->Participating() )
        return;

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int localSize = mpi::Pad( localHeight*localWidth );

    T* buffer = this->auxMemory_.Require( 2*localSize );
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[localSize];

    // Pack
    T* thisBuf = this->Buffer();
    const Int thisLDim = this->LDim();
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* thisCol = &thisBuf[jLoc*thisLDim];
        T* sendBufCol = &sendBuf[jLoc*localHeight];
        MemCopy( sendBufCol, thisCol, localHeight );
    }

    // AllReduce col
    mpi::AllReduce( sendBuf, recvBuf, localSize, g.ColComm() );

    // Unpack
    PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* recvBufCol = &recvBuf[jLoc*localHeight];
        T* thisCol = &thisBuf[jLoc*thisLDim];
        MemCopy( thisCol, recvBufCol, localHeight );
    }
    this->auxMemory_.Release();
}

template<typename T>
void
DM<T>::TransposeFrom( const DistMatrix<T,VR,STAR>& A, bool conjugate )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("[* ,MR]::TransposeFrom");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
    )
    const elem::Grid& g = this->Grid();
    this->AlignRowsAndResize( A.ColAlign()%g.Width(), A.Width(), A.Height() );
    if( !this->Participating() )
        return;

    if( this->RowAlign() == A.ColAlign() % g.Width() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();

        const Int width = this->Width();
        const Int height = this->Height();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalHeightOfA );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        if( conjugate )
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &sendBuf[jLoc*height];
                const T* sourceCol = &ABuf[jLoc];
                for( Int i=0; i<height; ++i )
                    destCol[i] = Conj( sourceCol[i*ALDim] );
            }
        }
        else
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &sendBuf[jLoc*height];
                const T* sourceCol = &ABuf[jLoc];
                for( Int i=0; i<height; ++i )
                    destCol[i] = sourceCol[i*ALDim];
            }
        }

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.ColComm() );

        // Unpack
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        const Int rowShift = this->RowShift();
        const Int colAlignOfA = A.ColAlign();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int colShiftOfA = Shift_( col+k*c, colAlignOfA, p );
            const Int rowOffset = (colShiftOfA-rowShift) / c;
            const Int localWidth = Length_( width, colShiftOfA, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc ) 
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MR].TransposeFrom[VR,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();
        const Int rank = g.VRRank();

        // Perform the SendRecv to make A have the same rowAlign
        const Int rowAlign = this->RowAlign();
        const Int colAlignOfA = A.ColAlign();
        const Int rowShift = this->RowShift();

        const Int sendRank = (rank+p+rowAlign-colAlignOfA) % p;
        const Int recvRank = (rank+p+colAlignOfA-rowAlign) % p;

        const Int width = this->Width();
        const Int height = this->Height();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLength(width,p);
        const Int portionSize = mpi::Pad( height*maxLocalHeightOfA );

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[portionSize];

        // Pack
        const Int ALDim = A.LDim();
        const T* ABuf = A.LockedBuffer();
        if( conjugate )
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &secondBuf[jLoc*height];
                const T* sourceCol = &ABuf[jLoc];
                for( Int i=0; i<height; ++i )
                    destCol[i] = Conj( sourceCol[i*ALDim] );
            }
        }
        else
        {
            PARALLEL_FOR
            for( Int jLoc=0; jLoc<localHeightOfA; ++jLoc )
            {
                T* destCol = &secondBuf[jLoc*height];
                const T* sourceCol = &ABuf[jLoc];
                for( Int i=0; i<height; ++i )
                    destCol[i] = sourceCol[i*ALDim];
            }
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
        T* thisBuf = this->Buffer();
        const Int thisLDim = this->LDim();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int colShiftOfA = Shift_( col+c*k, rowAlign, p );
            const Int rowOffset = (colShiftOfA-rowShift) / c;
            const Int localWidth = Length_( width, colShiftOfA, p );
            INNER_PARALLEL_FOR
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const T* dataCol = &data[jLoc*height];
                T* thisCol = &thisBuf[(rowOffset+jLoc*r)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
}

template<typename T>
void
DM<T>::AdjointFrom( const DistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[* ,MR]::AdjointFrom"))
    this->TransposeFrom( A, true );
}

// Basic queries
// =============

template<typename T>
elem::DistData DM<T>::DistData() const { return elem::DistData(*this); }

template<typename T>
mpi::Comm DM<T>::DistComm() const { return this->grid_->MRComm(); }
template<typename T>
mpi::Comm DM<T>::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM<T>::RedundantComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm DM<T>::ColComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM<T>::RowComm() const { return this->grid_->MRComm(); }

template<typename T>
Int DM<T>::ColStride() const { return 1; }
template<typename T>
Int DM<T>::RowStride() const { return this->grid_->Width(); }

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
