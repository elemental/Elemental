/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elemental {

template<typename T>
inline void
DistMatrix<T,MR,STAR>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;

    const int height      = this->Height();
    const int width       = this->Width();
    const int localHeight = this->LocalHeight();
    const int c           = g.Width();
    const int colShift    = this->ColShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Only one process row needs to participate
    if( g.MCRank() == 0 )
    {
        std::vector<T> sendBuf(height*width,0);
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            for( int j=0; j<width; ++j )
                sendBuf[(colShift+iLocal*c)+j*height] = 
                    thisLocalBuffer[iLocal+j*thisLDim];

        // If we are the root, allocate a receive buffer
        std::vector<T> recvBuf;
        if( g.MRRank() == 0 )
            recvBuf.resize( height*width );

        // Sum the contributions and send to the root
        mpi::Reduce
        ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, g.MRComm() );

        if( g.MRRank() == 0 )
        {
            // Print the data
            for( int i=0; i<height; ++i )
            {
                for( int j=0; j<width; ++j )
                    os << WrapScalar(recvBuf[i+j*height]) << " ";
                os << "\n";
            }
            os << std::endl;
        }
    }
    mpi::Barrier( g.VCComm() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::Align( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::Align");
    this->AssertFreeColAlignment();
#endif
    this->AlignCols( colAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::AlignCols( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Width() )
        throw std::runtime_error("Invalid column alignment for [MR,* ]");
#endif
    this->colAlignment_ = colAlignment;
    this->colShift_ = Shift( g.MRRank(), colAlignment, g.Width() );
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::View( DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ = A.ColShift();
    this->localMatrix_.View( A.LocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::View
( int height, int width, int colAlignment,
  T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->colShift_ = Shift(grid.MRRank(),colAlignment,grid.Width());
    const int localHeight = LocalLength(height,this->colShift_,grid.Width());
    this->localMatrix_.View( localHeight, width, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::LockedView( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ = A.ColShift();
    this->localMatrix_.LockedView( A.LockedLocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::LockedView
( int height, int width, int colAlignment,
  const T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->colShift_ = Shift(grid.MRRank(),colAlignment,grid.Width());
    const int localHeight = LocalLength(height,this->colShift_,grid.Width());
    this->localMatrix_.LockedView( localHeight, width, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::View
( DistMatrix<T,MR,STAR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    {
        const elemental::Grid& g = this->Grid();
        const int c   = g.Width();
        const int col = g.MRRank();

        this->colAlignment_ = (A.ColAlignment()+i) % c;
        this->colShift_ = Shift( col, this->ColAlignment(), c );

        const int localHeightBefore = LocalLength( i, A.ColShift(), c );
        const int localHeight = LocalLength( height, this->ColShift(), c );

        this->localMatrix_.View
        ( A.LocalMatrix(), localHeightBefore, j, localHeight, width );
    }
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::LockedView
( const DistMatrix<T,MR,STAR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    {
        const elemental::Grid& g = this->Grid();
        const int c   = g.Width();
        const int col = g.MRRank();

        this->colAlignment_ = (A.ColAlignment()+i) % c;
        this->colShift_ = Shift( col, this->ColAlignment(), c );

        const int localHeightBefore = LocalLength( i, A.ColShift(), c );
        const int localHeight = LocalLength( height, this->ColShift(), c );

        this->localMatrix_.LockedView
        ( A.LockedLocalMatrix(), localHeightBefore, j, localHeight, width );
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::View1x2
( DistMatrix<T,MR,STAR>& AL, DistMatrix<T,MR,STAR>& AR )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::View1x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->colShift_ = AL.ColShift();
    this->localMatrix_.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::LockedView1x2
( const DistMatrix<T,MR,STAR>& AL, const DistMatrix<T,MR,STAR>& AR )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::LockedView1x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->colShift_ = AL.ColShift();
    this->localMatrix_.LockedView1x2
    ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::View2x1
( DistMatrix<T,MR,STAR>& AT,
  DistMatrix<T,MR,STAR>& AB )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::View2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->colShift_ = AT.ColShift();
    this->localMatrix_.View2x1
    ( AT.LocalMatrix(),
      AB.LocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::LockedView2x1
( const DistMatrix<T,MR,STAR>& AT,
  const DistMatrix<T,MR,STAR>& AB )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::LockedView2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->colShift_ = AT.ColShift();
    this->localMatrix_.LockedView2x1
    ( AT.LockedLocalMatrix(),
      AB.LockedLocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::View2x2
( DistMatrix<T,MR,STAR>& ATL, DistMatrix<T,MR,STAR>& ATR,
  DistMatrix<T,MR,STAR>& ABL, DistMatrix<T,MR,STAR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::View2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->colAlignment_ = ATL.ColAlignment();
    this->colShift_ = ATL.ColShift();
    this->localMatrix_.View2x2
    ( ATL.LocalMatrix(), ATR.LocalMatrix(),
      ABL.LocalMatrix(), ABR.LocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::LockedView2x2
( const DistMatrix<T,MR,STAR>& ATL, const DistMatrix<T,MR,STAR>& ATR,
  const DistMatrix<T,MR,STAR>& ABL, const DistMatrix<T,MR,STAR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::LockedView2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->colAlignment_ = ATL.ColAlignment();
    this->colShift_ = ATL.ColShift();
    this->localMatrix_.LockedView2x2
    ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
      ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    this->localMatrix_.ResizeTo
    ( LocalLength(height,this->ColShift(),this->Grid().Width()), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
DistMatrix<T,MR,STAR>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // columns within each process row
    const elemental::Grid& g = this->Grid();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();

    T u;
    if( g.MRRank() == ownerCol )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        u = this->GetLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        this->SetLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        this->UpdateLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., SetToIdentity and MakeTrapezoidal
//

template<typename T>
inline void
DistMatrix<T,MR,STAR>::MakeTrapezoidal
( Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();    
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int c = this->Grid().Width();
    const int colShift = this->ColShift();

    if( uplo == LOWER )
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            int lastZeroRow = ( side==LEFT ? j-offset-1
                                           : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                int boundary = std::min( lastZeroRow+1, height );
                int numZeroRows = RawLocalLength( boundary, colShift, c );
                T* thisCol = &thisLocalBuffer[j*thisLDim];
                memset( thisCol, 0, numZeroRows*sizeof(T) );
            }
        }
    }
    else
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            int firstZeroRow = 
                ( side==LEFT ? std::max(j-offset+1,0)
                             : std::max(j-offset+height-width+1,0) );
            int numNonzeroRows = RawLocalLength(firstZeroRow,colShift,c);
            if( numNonzeroRows < localHeight )
            {
                T* thisCol = &thisLocalBuffer[numNonzeroRows+j*thisLDim];
                memset( thisCol, 0, (localHeight-numNonzeroRows)*sizeof(T) );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::ScaleTrapezoid
( T alpha, Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::ScaleTrapezoid");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();    
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int c = this->Grid().Width();
    const int colShift = this->ColShift();

    if( uplo == UPPER )
    {
        T* thisLocalBuffer =this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            int boundary = std::min( lastRow+1, height );
            int numRows = RawLocalLength( boundary, colShift, c );
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            for( int iLocal=0; iLocal<numRows; ++iLocal )
                thisCol[iLocal] *= alpha;
        }
    }
    else
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            int firstRow = ( side==LEFT ? std::max(j-offset,0)
                                        : std::max(j-offset+height-width,0) );
            int numZeroRows = RawLocalLength( firstRow, colShift, c );
            T* thisCol = &thisLocalBuffer[numZeroRows+j*thisLDim];
            for( int iLocal=0; iLocal<(localHeight-numZeroRows); ++iLocal )
                thisCol[iLocal] *= alpha;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int c = this->Grid().Width();
    const int colShift = this->ColShift();

    this->SetToZero();

    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*c;
        if( i < width )
            thisLocalBuffer[iLocal+i*thisLDim] = 1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const elemental::Grid& g = this->Grid();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int bufSize = localHeight*width;

    this->auxMemory_.Require( bufSize );

    // Create random matrix on process row 0, then broadcast
    T* buffer = this->auxMemory_.Buffer();
    if( g.MCRank() == 0 )
    {
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                buffer[i+j*localHeight] = SampleUnitBall<T>();
    }
    mpi::Broadcast( buffer, bufSize, 0, g.MCComm() );

    // Unpack
    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int j=0; j<width; ++j )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            thisLocalBuffer[iLocal+j*thisLDim] = buffer[iLocal+j*localHeight];

    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::SumOverCol()
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SumOverCol");
    this->AssertNotLockedView();
#endif
    const elemental::Grid& g = this->Grid();
    
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localSize = std::max( localHeight*width, mpi::MIN_COLL_MSG );

    this->auxMemory_.Require( 2*localSize );
    T* buffer = this->auxMemory_.Buffer();
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[localSize];

    // Pack
    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<width; ++j )
    {
        const T* thisCol = &thisLocalBuffer[j*thisLDim];
        T* sendBufCol = &sendBuf[j*localHeight];
        memcpy( sendBufCol, thisCol, localHeight*sizeof(T) );
    }

    // AllReduce sum
    mpi::AllReduce
    ( sendBuf, recvBuf, localSize, mpi::SUM, g.MCComm() );

    // Unpack
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<width; ++j )
    {
        const T* recvBufCol = &recvBuf[j*localHeight];
        T* thisCol = &thisLocalBuffer[j*thisLDim];
        memcpy( thisCol, recvBufCol, localHeight*sizeof(T) );
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::AdjointFrom( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ]::AdjointFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.RowAlignment();
            this->colShift_ = 
                Shift( g.MRRank(), this->ColAlignment(), g.Width() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->ColAlignment() == A.RowAlignment() )
    {
        const int r = g.Height();

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            std::max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack 
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localHeightOfA; ++jLocal )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                originalData[iLocal+jLocal*localHeight] = 
                    Conj( ALocalBuffer[jLocal+iLocal*ALDim] );

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShift = RawShift( k, colAlignmentOfA, r );
            const int localWidth = RawLocalLength( width, rowShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(rowShift+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [MR,* ]::AdjointFrom" << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendCol = (col+c+colAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-colAlignment) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack the currently owned local data of A into the second buffer
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localHeightOfA; ++jLocal )
            for( int iLocal=0; iLocal<localWidthOfA; ++iLocal )
                secondBuffer[iLocal+jLocal*localWidthOfA] = 
                    Conj( ALocalBuffer[jLocal+iLocal*ALDim] );

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuffer, portionSize, sendCol, 0,
          firstBuffer,  portionSize, recvCol, mpi::ANY_TAG, g.MRComm() );

        // Use the output of the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

        // Unpack the contents of each member of the process col
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int rowShift = RawShift( k, colAlignmentOfA, r );
            const int localWidth = RawLocalLength( width, rowShift, r );
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(rowShift+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,STAR>::TransposeFrom( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ]::TransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.RowAlignment();
            this->colShift_ = 
                Shift( g.MRRank(), this->ColAlignment(), g.Width() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->ColAlignment() == A.RowAlignment() )
    {
        const int r = g.Height();

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            std::max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack 
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localHeightOfA; ++jLocal )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                originalData[iLocal+jLocal*localHeight] = 
                    ALocalBuffer[jLocal+iLocal*ALDim];

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShift = RawShift( k, colAlignmentOfA, r );
            const int localWidth = RawLocalLength( width, rowShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(rowShift+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [MR,* ]::TransposeFrom" << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendCol = (col+c+colAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-colAlignment) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack the currently owned local data of A into the second buffer
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localHeightOfA; ++jLocal )
            for( int iLocal=0; iLocal<localWidthOfA; ++iLocal )
                secondBuffer[iLocal+jLocal*localWidthOfA] = 
                    ALocalBuffer[jLocal+iLocal*ALDim];

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuffer, portionSize, sendCol, 0,
          firstBuffer,  portionSize, recvCol, mpi::ANY_TAG, g.MRComm() );

        // Use the output of the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

        // Unpack the contents of each member of the process col
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int rowShift = RawShift( k, colAlignmentOfA, r );
            const int localWidth = RawLocalLength( width, rowShift, r );
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(rowShift+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,VC,STAR> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = A;

    std::auto_ptr< DistMatrix<T,VR,STAR> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(true,this->ColAlignment(),g) );
    *A_VR_STAR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( A.Width() == 1 )
    {
        if( !this->Viewing() )
            this->ResizeTo( A.Height(), 1 );

        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int myCol = g.MRRank();
        const int rankCM = g.VCRank();
        const int rankRM = g.VRRank();
        const int colAlignment = this->ColAlignment();
        const int colShift = this->ColShift();
        const int colAlignmentOfA = A.ColAlignment();
        const int colShiftOfA = A.ColShift();

        const int height = this->Height();
        const int maxLocalVectorHeight = MaxLocalLength(height,p);
        const int portionSize = 
            std::max(maxLocalVectorHeight,mpi::MIN_COLL_MSG);

        const int colShiftVR = Shift(rankRM,colAlignment,p);
        const int colShiftVCOfA = Shift(rankCM,colAlignmentOfA,p);
        const int sendRankRM = (rankRM+(p+colShiftVCOfA-colShiftVR)) % p;
        const int recvRankCM = (rankCM+(p+colShiftVR-colShiftVCOfA)) % p;
        const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

        this->auxMemory_.Require( (r+1)*portionSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        // A[VC,* ] <- A[MC,* ]
        {
            const int shift = Shift(rankCM,colAlignmentOfA,p);
            const int offset = (shift-colShiftOfA) / r;
            const int thisLocalHeight = LocalLength(height,shift,p);

            const T* ALocalBuffer = A.LockedLocalBuffer();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                sendBuf[iLocal] = ALocalBuffer[offset+iLocal*c];
        }

        // A[VR,* ] <- A[VC,* ]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankRM, 0,
          recvBuf, portionSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

        // A[MR,* ] <- A[VR,* ]
        mpi::AllGather
        ( recvBuf, portionSize,
          sendBuf, portionSize, g.MCComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &sendBuf[k*portionSize];

            const int shift = RawShift(myCol+c*k,colAlignment,p);
            const int offset = (shift-colShift) / c;
            const int thisLocalHeight = RawLocalLength(height,shift,p);

            for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                thisLocalBuffer[offset+iLocal*r] = data[iLocal];
        }
        this->auxMemory_.Release();
    }
    else
    {
        std::auto_ptr< DistMatrix<T,VC,STAR> > A_VC_STAR
        ( new DistMatrix<T,VC,STAR>(g) );
        *A_VC_STAR = A;

        std::auto_ptr< DistMatrix<T,VR,STAR> > A_VR_STAR
        ( new DistMatrix<T,VR,STAR>(true,this->ColAlignment(),g) );
        *A_VR_STAR = *A_VC_STAR;
        delete A_VC_STAR.release(); // lowers memory highwater

        *this = *A_VR_STAR;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR   = A;

    std::auto_ptr< DistMatrix<T,VC,STAR> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    std::auto_ptr< DistMatrix<T,VR,STAR> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(true,this->ColAlignment(),g) );
    *A_VR_STAR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MR,* ] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MR,* ] = [* ,MD] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            this->colShift_ = 
                Shift( g.MRRank(), this->ColAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() )
    {
        const int r = g.Height();

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            std::max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack 
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[jLocal*ALDim];
            T* originalDataCol = &originalData[jLocal*localHeight];
            memcpy( originalDataCol, ACol, localHeight*sizeof(T) );
        }

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowAlignmentOfA = A.RowAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShift = RawShift( k, rowAlignmentOfA, r );
            const int localWidth = RawLocalLength( width, rowShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(rowShift+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [MR,* ] <- [MR,MC]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack the currently owned local data of A into the second buffer
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[jLocal*ALDim];
            T* secondBufferCol = &secondBuffer[jLocal*localHeightOfA];
            memcpy( secondBufferCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuffer, portionSize, sendCol, 0,
          firstBuffer,  portionSize, recvCol, mpi::ANY_TAG, g.MRComm() );

        // Use the output of the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

        // Unpack the contents of each member of the process col
        const int rowAlignmentOfA = A.RowAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int rowShift = RawShift( k, rowAlignmentOfA, r );
            const int localWidth = RawLocalLength( width, rowShift, r );
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(rowShift+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            this->colShift_ = A.ColShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() )
    {
        this->localMatrix_ = A.LockedLocalMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [MR,* ] <- [MR,* ]." << std::endl;
#endif
        const int rank = g.MRRank();
        const int c = g.Width();

        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (rank+c+colAlignment-colAlignmentOfA) % c;
        const int recvRank = (rank+c+colAlignmentOfA-colAlignment) % c;

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();

        const int sendSize = localHeightOfA * width;
        const int recvSize = localHeight * width;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            const T* ACol = &ALocalBuffer[j*ALDim];
            T* sendBufferCol = &sendBuffer[j*localHeightOfA];
            memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.MRComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*localHeight];
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }

        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MR,MC> A_MR_MC(g);

    A_MR_MC = A;
    *this   = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,VR,STAR> A_VR_STAR(true,this->ColAlignment(),g);

    A_VR_STAR = A;
    *this     = A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MR,MC> A_MR_MC(g);

    A_MR_MC = A;
    *this   = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.VCRank() == 0 )
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
            this->colShift_ = 
                Shift( g.MRRank(), this->ColAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() % g.Width() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(height,p);

        const int portionSize = 
            std::max(maxLocalHeightOfA*width,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            const T* ACol = &ALocalBuffer[j*ALDim];
            T* originalDataCol = &originalData[j*localHeightOfA];
            memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int colShift = this->ColShift();
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShiftOfA = RawShift( col+c*k, colAlignmentOfA, p );
            const int colOffset = (colShiftOfA-colShift) / c;
            const int localHeight = RawLocalLength( height, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colOffset+iLocal*r)+j*thisLDim] = 
                        data[iLocal+j*localHeight];
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [MR,* ] <- [VR,* ]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rank = g.VRRank();

        // Perform the SendRecv to make A have the same colAlignment
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int colShift = this->ColShift();

        const int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(height,p);

        const int portionSize = 
            std::max(maxLocalHeightOfA*width,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            const T* ACol = &ALocalBuffer[j*ALDim];
            T* secondBufferCol = &secondBuffer[j*localHeightOfA];
            memcpy( secondBufferCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, mpi::ANY_TAG, g.VRComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int colShiftOfA = RawShift( col+c*k, colAlignment, p );
            const int colOffset = (colShiftOfA-colShift) / c;
            const int localHeight = RawLocalLength( height, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colOffset+iLocal*r)+j*thisLDim] = 
                        data[iLocal+j*localHeight];
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(g) );
    *A_STAR_VC = A;

    std::auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
    ( new DistMatrix<T,MR,MC>(true,false,this->ColAlignment(),0,g) );
    *A_MR_MC = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    *this = *A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,STAR>&
DistMatrix<T,MR,STAR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int c = g.Width();
    const int colShift = this->ColShift();

    const int width = this->Width();
    const int localHeight = this->LocalHeight();

    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const int ALDim = A.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int j=0; j<width; ++j )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            thisLocalBuffer[iLocal+j*thisLDim] = 
                ALocalBuffer[(colShift+iLocal*c)+j*ALDim];
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

} // namespace elemental
