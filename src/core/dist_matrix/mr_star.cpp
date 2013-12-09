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
using ADM = AbstractDistMatrix<T>;
template<typename T>
using DM = DistMatrix<T,MR,STAR>;

template<typename T>
DM<T>::DistMatrix( const elem::Grid& g )
: ADM<T>(g)
{ this->SetShifts(); } 

template<typename T>
DM<T>::DistMatrix( Int height, Int width, const elem::Grid& g )
: ADM<T>(g)
{ this->SetShifts(); this->ResizeTo(height,width); }

template<typename T>
DM<T>::DistMatrix( Int height, Int width, Int colAlign, const elem::Grid& g )
: ADM<T>(g)
{ this->Align(colAlign,0); this->ResizeTo(height,width); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, Int ldim, const elem::Grid& g )
: ADM<T>(g)
{ this->Align(colAlign,0); this->ResizeTo(height,width,ldim); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, const T* buffer, Int ldim,
  const elem::Grid& g )
: ADM<T>(g)
{ this->LockedAttach( height, width, colAlign, buffer, ldim, g ); }

template<typename T>
DM<T>::DistMatrix
( Int height, Int width, Int colAlign, T* buffer, Int ldim,
  const elem::Grid& g )
: ADM<T>(g)
{ this->Attach( height, width, colAlign, buffer, ldim, g ); }

template<typename T>
DM<T>::DistMatrix( const DM<T>& A )
: ADM<T>(A.Grid())
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
DM<T>::DistMatrix( const DistMatrix<T,U,V>& A )
: ADM<T>(A.Grid())
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix[MR,* ]::DistMatrix");
#endif
    this->SetShifts();
    if( MR != U || STAR != V || 
        reinterpret_cast<const DM<T>*>(&A) != this )
        *this = A;
    else
        LogicError("Tried to construct [MR,* ] with itself");
}

template<typename T>
DM<T>::DistMatrix( DM<T>&& A )
: ADM<T>(std::move(A))
{ }

template<typename T>
DM<T>&
DM<T>::operator=( DM<T>&& A )
{
    ADM<T>::operator=( std::move(A) );
    return *this;
}

template<typename T>
DM<T>::~DistMatrix()
{ }

template<typename T>
elem::DistData
DM<T>::DistData() const
{ return elem::DistData(*this); }

template<typename T>
mpi::Comm
DM<T>::DistComm() const
{ return this->grid_->MRComm(); }

template<typename T>
mpi::Comm
DM<T>::CrossComm() const
{ return mpi::COMM_SELF; }

template<typename T>
mpi::Comm
DM<T>::RedundantComm() const
{ return this->grid_->MCComm(); }

template<typename T>
mpi::Comm
DM<T>::ColComm() const
{ return this->grid_->MRComm(); }

template<typename T>
mpi::Comm
DM<T>::RowComm() const
{ return mpi::COMM_SELF; }

template<typename T>
Int
DM<T>::ColStride() const
{ return this->grid_->Width(); }
    
template<typename T>
Int 
DM<T>::RowStride() const
{ return 1; }

template<typename T>
void
DM<T>::AlignWith( const elem::DistData& data )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::AlignWith");
#endif
    const Grid& grid = *data.grid;
    this->SetGrid( grid );

    if( data.colDist == MR )
        this->colAlign_ = data.colAlign;
    else if( data.rowDist == MR )
        this->colAlign_ = data.rowAlign;
    else if( data.colDist == VR )
        this->colAlign_ = data.colAlign % this->ColStride();
    else if( data.rowDist == VR )
        this->colAlign_ = data.rowAlign % this->ColStride();
#ifndef RELEASE
    else LogicError("Nonsensical alignment");
#endif
    this->colConstrained_ = true;
    this->SetShifts();
}

template<typename T>
void
DM<T>::AlignWith( const ADM<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DM<T>::AlignColsWith( const elem::DistData& data )
{ this->AlignWith( data ); }

template<typename T>
void
DM<T>::AlignColsWith( const ADM<T>& A )
{ this->AlignWith( A.DistData() ); }

template<typename T>
void
DM<T>::Attach
( Int height, Int width, Int colAlign,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::Attach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlign_ = colAlign;
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
DM<T>::LockedAttach
( Int height, Int width, Int colAlign,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlign_ = colAlign;
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
DM<T>::Attach( Matrix<T>& A, Int colAlign, const elem::Grid& g )
{ this->Attach( A.Height(), A.Width(), colAlign, A.Buffer(), A.LDim(), g ); }

template<typename T>
void
DM<T>::LockedAttach( const Matrix<T>& A, Int colAlign, const elem::Grid& g )
{
    this->LockedAttach
    ( A.Height(), A.Width(), colAlign, A.LockedBuffer(), A.LDim(), g );
}


//
// Utility functions, e.g., SumOverCol
//

template<typename T>
void
DM<T>::SumOverCol()
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
DM<T>::AdjointFrom( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::AdjointFrom");
#endif
    this->TransposeFrom( A, true );
}

template<typename T>
void
DM<T>::TransposeFrom( const DistMatrix<T,MC,MR>& A, bool conjugate )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ]::TransposeFrom");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->AlignColsAndResize( A.RowAlign(), A.Width(), A.Height() );
    if( !this->Participating() )
        return;

    if( this->ColAlign() == A.RowAlign() )
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
        const Int colAlignOfA = A.ColAlign();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int rowShift = Shift_( k, colAlignOfA, r );
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

        const Int colAlign = this->ColAlign();
        const Int rowAlignOfA = A.RowAlign();
        const Int sendCol = (col+c+colAlign-rowAlignOfA) % c;
        const Int recvCol = (col+c+rowAlignOfA-colAlign) % c;

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
        const Int colAlignOfA = A.ColAlign();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShift = Shift_( k, colAlignOfA, r );
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [MC,MR]");
#endif
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,VC,STAR>> A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = A;

    std::unique_ptr<DistMatrix<T,VR,STAR>> A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(true,this->ColAlign(),g) );
    *A_VR_STAR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MC,STAR>& A )
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
        const Int colAlign = this->ColAlign();
        const Int colShift = this->ColShift();
        const Int colAlignOfA = A.ColAlign();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int maxLocalVectorHeight = MaxLength(height,p);
        const Int portionSize = mpi::Pad( maxLocalVectorHeight );

        const Int colShiftVR = Shift(rankRM,colAlign,p);
        const Int colShiftVCOfA = Shift(rankCM,colAlignOfA,p);
        const Int sendRankRM = (rankRM+(p+colShiftVCOfA-colShiftVR)) % p;
        const Int recvRankCM = (rankCM+(p+colShiftVR-colShiftVCOfA)) % p;
        const Int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

        T* buffer = this->auxMemory_.Require( (r+1)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        // A[VC,* ] <- A[MC,* ]
        {
            const Int shift = Shift(rankCM,colAlignOfA,p);
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
            const Int shift = Shift_(myCol+c*k,colAlign,p);
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
        ( new DistMatrix<T,VR,STAR>(true,this->ColAlign(),g) );
        *A_VR_STAR = *A_VC_STAR;
        delete A_VC_STAR.release(); // lowers memory highwater

        *this = *A_VR_STAR;
    }
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MR>& A )
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
    ( new DistMatrix<T,VR,STAR>(true,this->ColAlign(),g) );
    *A_VR_STAR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_VR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MD,STAR>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MD>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [MR,MC]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
    const elem::Grid& g = this->Grid();
    this->AlignColsAndResize( A.ColAlign(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlign() == A.ColAlign() )
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
        const Int rowAlignOfA = A.RowAlign();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int rowShift = Shift_( k, rowAlignOfA, r );
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

        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();
        const Int sendCol = (col+c+colAlign-colAlignOfA) % c;
        const Int recvCol = (col+c+colAlignOfA-colAlign) % c;

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
        const Int rowAlignOfA = A.RowAlign();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuf[k*portionSize];
            const Int rowShift = Shift_( k, rowAlignOfA, r );
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
const DM<T>&
DM<T>::operator=( const DM<T>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [MR,* ]");
    this->AssertNotLocked();
    this->AssertSameGrid( A.Grid() );
#endif
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
            std::cerr << "Unaligned [MR,* ] <- [MR,* ]." << std::endl;
#endif
        const Int rank = g.Col();
        const Int c = g.Width();

        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();

        const Int sendRank = (rank+c+colAlign-colAlignOfA) % c;
        const Int recvRank = (rank+c+colAlignOfA-colAlign) % c;

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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [* ,MC]");
#endif
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [VC,* ]");
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VR,STAR> A_VR_STAR(true,this->ColAlign(),g);

    A_VR_STAR = A;
    *this = A_VR_STAR;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [* ,VC]");
#endif
    DistMatrix<T,MR,MC> A_MR_MC( A );
    *this = A_MR_MC;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,VR,STAR>& A )
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
    this->AlignColsAndResize
    ( A.ColAlign()%g.Width(), A.Height(), A.Width() );
    if( !this->Participating() )
        return *this;

    if( this->ColAlign() == A.ColAlign() % g.Width() )
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
        const Int colAlignOfA = A.ColAlign();
        OUTER_PARALLEL_FOR
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuf[k*portionSize];
            const Int colShiftOfA = Shift_( col+c*k, colAlignOfA, p );
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

        // Perform the SendRecv to make A have the same colAlign
        const Int colAlign = this->ColAlign();
        const Int colAlignOfA = A.ColAlign();
        const Int colShift = this->ColShift();

        const Int sendRank = (rank+p+colAlign-colAlignOfA) % p;
        const Int recvRank = (rank+p+colAlignOfA-colAlign) % p;

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
            const Int colShiftOfA = Shift_( col+c*k, colAlign, p );
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    CallStackEntry cse("[MR,* ] = [* ,VR]");
#endif
    const elem::Grid& g = this->Grid();
    std::unique_ptr<DistMatrix<T,STAR,VC>> A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(A) );

    std::unique_ptr<DistMatrix<T,MR,MC>> A_MR_MC
    ( new DistMatrix<T,MR,MC>(true,false,this->ColAlign(),0,g) );
    *A_MR_MC = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    *this = *A_MR_MC;
    return *this;
}

template<typename T>
const DM<T>&
DM<T>::operator=( const DistMatrix<T,STAR,STAR>& A )
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
const DM<T>&
DM<T>::operator=( const DistMatrix<T,CIRC,CIRC>& A )
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
#define COPY(T,U,V) \
  template DistMatrix<T,MR,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
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
