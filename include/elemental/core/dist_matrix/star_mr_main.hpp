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
DistMatrix<T,STAR,MR>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;

    const int height     = this->Height();
    const int width      = this->Width();
    const int localWidth = this->LocalWidth();
    const int c          = g.Width();
    const int rowShift   = this->RowShift();

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
        for( int i=0; i<height; ++i )
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                sendBuf[i+(rowShift+jLocal*c)*height] = 
                    thisLocalBuffer[i+jLocal*thisLDim];

        // If we are the root, allocate the receive buffer
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
DistMatrix<T,STAR,MR>::Align( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Align");
    this->AssertFreeRowAlignment();
#endif
    this->AlignRows( rowAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::AlignRows( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Width() )
        throw std::runtime_error("Invalid row alignment for [* ,MR]");
#endif
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift( g.MRRank(), rowAlignment, g.Width() );
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::View( DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->rowAlignment_ = A.RowAlignment();
    this->rowShift_ = A.RowShift();
    this->localMatrix_.View( A.LocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::View
( int height, int width, int rowAlignment,
  T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift(grid.MRRank(),rowAlignment,grid.Width());
    const int localWidth = LocalLength(width,this->rowShift_,grid.Width());
    this->localMatrix_.View( height, localWidth, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::LockedView( const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[*, MR]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->rowAlignment_ = A.RowAlignment();
    this->rowShift_ = A.RowShift();
    this->localMatrix_.LockedView( A.LockedLocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::LockedView
( int height, int width, int rowAlignment,
  const T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift(grid.MRRank(),rowAlignment,grid.Width());
    const int localWidth = LocalLength(width,this->rowShift_,grid.Width());
    this->localMatrix_.LockedView( height, localWidth, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::View
( DistMatrix<T,STAR,MR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::View");
    this->AssertFreeRowAlignment();
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

        this->rowAlignment_ = (A.RowAlignment()+j) % c;
        this->rowShift_ = Shift( col, this->RowAlignment(), c );

        const int localWidthBefore = LocalLength( j, A.RowShift(), c );
        const int localWidth = LocalLength( width, this->RowShift(), c );

        this->localMatrix_.View
        ( A.LocalMatrix(), i, localWidthBefore, height, localWidth );
    }
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::LockedView
( const DistMatrix<T,STAR,MR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    {
        const elemental::Grid& g = this->Grid();
        const int c = g.Width();
        const int col = g.MRRank();

        this->rowAlignment_ = (A.RowAlignment()+j) % c;
        this->rowShift_ = Shift( col, this->RowAlignment(), c );

        const int localWidthBefore = LocalLength( j, A.RowShift(), c );
        const int localWidth = LocalLength( width, this->RowShift(), c );

        this->localMatrix_.LockedView
        ( A.LockedLocalMatrix(), i, localWidthBefore, height, localWidth );
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::View1x2
( DistMatrix<T,STAR,MR>& AL, DistMatrix<T,STAR,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::View1x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->rowAlignment_ = AL.RowAlignment();
    this->rowShift_ = AL.RowShift();
    this->localMatrix_.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::LockedView1x2
( const DistMatrix<T,STAR,MR>& AL, const DistMatrix<T,STAR,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::LockedView1x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->rowAlignment_ = AL.RowAlignment();
    this->rowShift_ = AL.RowShift();
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
DistMatrix<T,STAR,MR>::View2x1
( DistMatrix<T,STAR,MR>& AT,
  DistMatrix<T,STAR,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::View2x1");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->rowAlignment_ = AT.RowAlignment();
    this->rowShift_ = AT.RowShift();
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
DistMatrix<T,STAR,MR>::LockedView2x1
( const DistMatrix<T,STAR,MR>& AT,
  const DistMatrix<T,STAR,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::LockedView2x1");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->rowAlignment_ = AT.RowAlignment();
    this->rowShift_ = AT.RowShift();
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
DistMatrix<T,STAR,MR>::View2x2
( DistMatrix<T,STAR,MR>& ATL, DistMatrix<T,STAR,MR>& ATR,
  DistMatrix<T,STAR,MR>& ABL, DistMatrix<T,STAR,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::View2x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->rowAlignment_ = ATL.RowAlignment();
    this->rowShift_ = ATL.RowShift();
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
DistMatrix<T,STAR,MR>::LockedView2x2
( const DistMatrix<T,STAR,MR>& ATL, const DistMatrix<T,STAR,MR>& ATR,
  const DistMatrix<T,STAR,MR>& ABL, const DistMatrix<T,STAR,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::LockedView2x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->rowAlignment_ = ATL.RowAlignment();
    this->rowShift_ = ATL.RowShift();
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
DistMatrix<T,STAR,MR>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    this->localMatrix_.ResizeTo
    ( height, LocalLength(width,this->RowShift(),this->Grid().Width()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
DistMatrix<T,STAR,MR>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // column within each process row
    const elemental::Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    T u;
    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = this->GetLocalEntry(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        this->SetLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        this->UpdateLocalEntry(i,jLoc,u);
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
DistMatrix<T,STAR,MR>::MakeTrapezoidal
( Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const elemental::Grid& g = this->Grid();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int c = g.Width();
    const int rowShift = this->RowShift();

    if( uplo == LOWER )
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int lastZeroRow = ( side==LEFT ? j-offset-1
                                           : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                int boundary = std::min( lastZeroRow+1, height );
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                memset( thisCol, 0, boundary*sizeof(T) );
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
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int firstZeroRow = 
                ( side==LEFT ? std::max(j-offset+1,0)
                             : std::max(j-offset+height-width+1,0) );
            if( firstZeroRow < height )
            {
                T* thisCol = &thisLocalBuffer[firstZeroRow+jLocal*thisLDim];
                memset( thisCol, 0, (height-firstZeroRow)*sizeof(T) );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::ScaleTrapezoidal
( T alpha, Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const elemental::Grid& g = this->Grid();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int c = g.Width();
    const int rowShift = this->RowShift();

    if( uplo == UPPER )
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            int boundary = std::min( lastRow+1, height );
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            for( int i=0; i<boundary; ++i )
                thisCol[i] *= alpha;
        }
    }
    else
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int firstRow = ( side==LEFT ? std::max(j-offset,0)
                                        : std::max(j-offset+height-width,0) );
            T* thisCol = &thisLocalBuffer[firstRow+jLocal*thisLDim];
            for( int i=0; i<(height-firstRow); ++i )
                thisCol[i] *= alpha;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int c = this->Grid().Width();
    const int rowShift = this->RowShift();

    this->SetToZero();

    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*c;
        if( j < height )
            thisLocalBuffer[j+jLocal*thisLDim] = 1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const elemental::Grid& g = this->Grid();
    const int height     = this->Height();
    const int localWidth = this->LocalWidth();
    const int bufSize    = height*localWidth;

    this->auxMemory_.Require( bufSize );

    // Create random matrix on process row 0, then broadcast
    T* buffer = this->auxMemory_.Buffer();
    if( g.MCRank() == 0 )
    {
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                buffer[i+j*height] = SampleUnitBall<T>();
    }
    mpi::Broadcast( buffer, bufSize, 0, g.MCComm() );

    // Unpack
    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* bufferCol = &buffer[jLocal*height];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        memcpy( thisCol, bufferCol, height*sizeof(T) );
    }
    this->auxMemory_.Release();

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> 
inline void
DistMatrix<T,STAR,MR>::SumOverCol()
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SumOverCol");
    this->AssertNotLockedView();
#endif
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int localSize = std::max( localHeight*localWidth, mpi::MIN_COLL_MSG );

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
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        T* sendBufCol = &sendBuf[jLocal*localHeight];
        memcpy( sendBufCol, thisCol, localHeight*sizeof(T) );
    }

    // AllReduce col
    mpi::AllReduce
    ( sendBuf, recvBuf, localSize, mpi::SUM, this->Grid().MCComm() );

    // Unpack
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* recvBufCol = &recvBuf[jLocal*localHeight];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        memcpy( thisCol, recvBufCol, localHeight*sizeof(T) );
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::AdjointFrom( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR]::AdjointFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.ColAlignment() % g.Width();
            this->rowShift_ = 
                Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->RowAlignment() == A.ColAlignment() % g.Width() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();

        const int width = this->Width();
        const int height = this->Height();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

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
            for( int i=0; i<height; ++i )
                originalData[i+jLocal*height] = 
                    Conj( ALocalBuffer[jLocal+i*ALDim] );

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShiftOfA = RawShift( col+k*c, colAlignmentOfA, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = RawLocalLength( width, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal ) 
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [* ,MR].AdjointFrom[VR,* ]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rank = g.VRRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowShift = this->RowShift();

        const int sendRank = (rank+p+rowAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-rowAlignment) % p;

        const int width = this->Width();
        const int height = this->Height();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localHeightOfA; ++jLocal )
            for( int i=0; i<height; ++i )
                secondBuffer[i+jLocal*height] = 
                    Conj( ALocalBuffer[jLocal+i*ALDim] );

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

            const int colShiftOfA = RawShift( col+c*k, rowAlignment, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = RawLocalLength( width, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, height*sizeof(T) );
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
DistMatrix<T,STAR,MR>::TransposeFrom
( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR]::TransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.ColAlignment() % g.Width();
            this->rowShift_ = 
                Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->RowAlignment() == A.ColAlignment() % g.Width() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();

        const int width = this->Width();
        const int height = this->Height();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

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
            for( int i=0; i<height; ++i )
                originalData[i+jLocal*height] = ALocalBuffer[jLocal+i*ALDim];

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShiftOfA = RawShift( col+k*c, colAlignmentOfA, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = RawLocalLength( width, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [* ,MR].TransposeFrom[VR,* ]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rank = g.VRRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowShift = this->RowShift();

        const int sendRank = (rank+p+rowAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-rowAlignment) % p;

        const int width = this->Width();
        const int height = this->Height();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localHeightOfA; ++jLocal )
            for( int i=0; i<height; ++i )
                secondBuffer[i+jLocal*height] = ALocalBuffer[jLocal+i*ALDim];

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

            const int colShiftOfA = RawShift( col+c*k, rowAlignment, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = RawLocalLength( width, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.VCRank() == 0 )
    {
        std::cerr << 
          "The matrix redistribution [* ,MR] <- [MC,MR] potentially causes a "
          "large amount of cache-thrashing. If possible, avoid it by "
          "performing the redistribution with a (conjugate)-transpose: \n"
          << "  [MR,* ].(Conjugate)TransposeFrom([MC,MR])" << std::endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            this->rowShift_ = 
                Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        if( A.Height() == 1 )
        {
            const int localWidth = this->LocalWidth();

            this->auxMemory_.Require( localWidth );
            T* bcastBuf = this->auxMemory_.Buffer();

            if( g.MCRank() == A.ColAlignment() )
            {
                this->localMatrix_ = A.LockedLocalMatrix();

                // Pack
                const T* thisLocalBuffer = this->LockedLocalBuffer();
                const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                    bcastBuf[jLocal] = thisLocalBuffer[jLocal*thisLDim];
            }

            // Communicate
            mpi::Broadcast
            ( bcastBuf, localWidth, A.ColAlignment(), g.MCComm() );

            // Unpack
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                thisLocalBuffer[jLocal*thisLDim] = bcastBuf[jLocal];

            this->auxMemory_.Release();
        }
        else
        {
            const int r = g.Height();
            const int height = this->Height();
            const int localWidth = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int portionSize = 
                std::max(maxLocalHeight*localWidth,mpi::MIN_COLL_MSG);

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
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[jLocal*ALDim];
                T* originalDataCol = &originalData[jLocal*localHeightOfA];
                memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
            }

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

                const int colShift = RawShift( k, colAlignmentOfA, r );
                const int localHeight = RawLocalLength( height, colShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                    for( int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisLocalBuffer[(colShift+iLocal*r)+jLocal*thisLDim] =
                            data[iLocal+jLocal*localHeight];
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [* ,MR] <- [MC,MR]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

        if( A.Height() == 1 )
        {
            const int localWidth = this->LocalWidth();
            T* bcastBuf;

            if( g.MCRank() == A.ColAlignment() )
            {
                const int localWidthOfA = A.LocalWidth();

                this->auxMemory_.Require( localWidth+localWidthOfA );
                T* buffer = this->auxMemory_.Buffer();
                T* sendBuf = &buffer[0];
                bcastBuf   = &buffer[localWidthOfA];

                // Pack
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const int ALDim = A.LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                    sendBuf[jLocal] = ALocalBuffer[jLocal*ALDim];

                // Communicate
                mpi::SendRecv
                ( sendBuf,  localWidthOfA, sendCol, 0,
                  bcastBuf, localWidth,    recvCol, mpi::ANY_TAG,
                  g.MRComm() );
            }
            else
            {
                this->auxMemory_.Require( localWidth );
                bcastBuf = this->auxMemory_.Buffer();
            }

            // Communicate
            mpi::Broadcast
            ( bcastBuf, localWidth, A.ColAlignment(), g.MCComm() );

            // Unpack
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                thisLocalBuffer[jLocal*thisLDim] = bcastBuf[jLocal];
            this->auxMemory_.Release();
        }
        else
        {
            const int height = this->Height();
            const int localWidth  = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int localWidthOfA  = A.LocalWidth();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int portionSize = 
                std::max(maxLocalHeight*localWidth,mpi::MIN_COLL_MSG);

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
            const int colAlignmentOfA = A.ColAlignment();
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                const T* data = &secondBuffer[k*portionSize];

                const int colShift = RawShift( k, colAlignmentOfA, r );
                const int localHeight = RawLocalLength( height, colShift, r );
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                    for( int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisLocalBuffer[(colShift+iLocal*r)+jLocal*thisLDim] =
                            data[iLocal+jLocal*localHeight];
            }
            this->auxMemory_.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(false,true,0,this->RowAlignment(),g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            this->rowShift_ = A.RowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        this->localMatrix_ = A.LockedLocalMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [* ,MR] <- [* ,MR]." << std::endl;
#endif
        const int rank = g.MRRank();
        const int c = g.Width();

        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvRank = (rank+c+rowAlignmentOfA-rowAlignment) % c;

        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = height * localWidthOfA;
        const int recvSize = height * localWidth;

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
        for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[jLocal*ALDim];
            T* sendBufferCol = &sendBuffer[jLocal*height];
            memcpy( sendBufferCol, ACol, height*sizeof(T) );
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
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufferCol = &recvBuffer[jLocal*height];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            memcpy( thisCol, recvBufferCol, height*sizeof(T) );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MR] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MR] = [* ,MD] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(g) );
    *A_STAR_VC = A;

    std::auto_ptr< DistMatrix<T,STAR,VR> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(true,this->RowAlignment(),g) );
    *A_STAR_VR = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    *this = *A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,VR,STAR> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(g) );
    *A_VR_STAR = A;

    std::auto_ptr< DistMatrix<T,VC,STAR> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    std::auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(false,true,0,this->RowAlignment(),g) );
    *A_MC_MR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(g) );
    *A_STAR_VC = A;

    std::auto_ptr< DistMatrix<T,STAR,VR> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(true,this->RowAlignment(),g) );
    *A_STAR_VR = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    std::auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(false,true,0,this->RowAlignment(),g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR> A_STAR_VR(true,this->RowAlignment(),g);

    A_STAR_VR = A;
    *this = A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,VC,STAR> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = A;

    std::auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(false,true,0,this->RowAlignment(),g) );
    *A_MC_MR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment() % g.Width();
            this->rowShift_ = 
                Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() % g.Width() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();

        const int width = this->Width();
        const int height = this->Height();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidthOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalWidthOfA,mpi::MIN_COLL_MSG);

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
            T* originalDataCol = &originalData[jLocal*height];
            memcpy( originalDataCol, ACol, height*sizeof(T) );
        }

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int rowAlignmentOfA = A.RowAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShiftOfA = RawShift( col+k*c, rowAlignmentOfA, p );
            const int rowOffset = (rowShiftOfA-rowShift) / c;
            const int localWidth = RawLocalLength( width, rowShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [* ,MR] <- [* ,VR]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rank = g.VRRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int rowShift = this->RowShift();

        const int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const int width = this->Width();
        const int height = this->Height();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidthOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalWidthOfA,mpi::MIN_COLL_MSG);

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
        for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[jLocal*ALDim];
            T* secondBufferCol = &secondBuffer[jLocal*height];
            memcpy( secondBufferCol, ACol, height*sizeof(T) );
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

            const int rowShiftOfA = RawShift( col+c*k, rowAlignment, p );
            const int rowOffset = (rowShiftOfA-rowShift) / c;
            const int localWidth = RawLocalLength( width, rowShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, height*sizeof(T) );
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
inline const DistMatrix<T,STAR,MR>&
DistMatrix<T,STAR,MR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int c = this->Grid().Width();
    const int rowShift = this->RowShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();

    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const int ALDim = A.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* ACol = &ALocalBuffer[(rowShift+jLocal*c)*ALDim];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        memcpy( thisCol, ACol, localHeight*sizeof(T) );
    }

#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

} // namespace elemental
