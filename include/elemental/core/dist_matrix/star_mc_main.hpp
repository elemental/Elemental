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
DistMatrix<T,STAR,MC>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;

    const int height     = this->Height();
    const int width      = this->Width();
    const int localWidth = this->LocalWidth();
    const int r          = g.Height();
    const int rowShift   = this->RowShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Only one process col needs to participate
    if( g.MRRank() == 0 )
    {
        std::vector<T> sendBuf(height*width,0);
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int i=0; i<height; ++i )
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                sendBuf[i+(rowShift+jLocal*r)*height] = 
                    thisLocalBuffer[i+jLocal*thisLDim];

        // If we are the root, allocate a receive buffer
        std::vector<T> recvBuf;
        if( g.MCRank() == 0 )
            recvBuf.resize( height*width );

        // Sum the contributions and send to the root
        mpi::Reduce
        ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, g.MCComm() );

        if( g.MCRank() == 0 )
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
DistMatrix<T,STAR,MC>::Align( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::Align");
    this->AssertFreeRowAlignment();
#endif
    this->AlignRows( rowAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MC>::AlignRows( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Height() )
        throw std::runtime_error("Invalid row alignment for [* ,MC]");
#endif
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift( g.MCRank(), rowAlignment, g.Height() );
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
DistMatrix<T,STAR,MC>::View( DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View");
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
DistMatrix<T,STAR,MC>::View
( int height, int width, int rowAlignment,
  T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift(grid.MCRank(),rowAlignment,grid.Height());
    const int localWidth = LocalLength(width,this->rowShift_,grid.Height());
    this->localMatrix_.View( height, localWidth, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MC>::LockedView( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE 
    PushCallStack("[* ,MC]::LockedView");
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
DistMatrix<T,STAR,MC>::LockedView
( int height, int width, int rowAlignment,
  const T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift(grid.MCRank(),rowAlignment,grid.Height());
    const int localWidth = LocalLength(width,this->rowShift_,grid.Height());
    this->localMatrix_.LockedView( height, localWidth, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MC>::View
( DistMatrix<T,STAR,MC>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    {
        const elemental::Grid& g = this->Grid();
        const int r   = g.Height();
        const int row = g.MCRank();

        this->rowAlignment_ = (A.RowAlignment()+j) % r;
        this->rowShift_ = Shift( row, this->RowAlignment(), r ); 

        const int localWidthBefore = LocalLength( j, A.RowShift(), r );
        const int localWidth = LocalLength( width, this->RowShift(), r );

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
DistMatrix<T,STAR,MC>::LockedView
( const DistMatrix<T,STAR,MC>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    {
        const elemental::Grid& g = this->Grid();
        const int r   = g.Height();
        const int row = g.MCRank();

        this->rowAlignment_ = (A.RowAlignment()+j) % r;
        this->rowShift_ = Shift( row, this->RowAlignment(), r );

        const int localWidthBefore = LocalLength( j, A.RowShift(), r );
        const int localWidth = LocalLength( width, this->RowShift(), r );

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
DistMatrix<T,STAR,MC>::View1x2
( DistMatrix<T,STAR,MC>& AL, DistMatrix<T,STAR,MC>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View1x2");
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
DistMatrix<T,STAR,MC>::LockedView1x2
( const DistMatrix<T,STAR,MC>& AL, const DistMatrix<T,STAR,MC>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::LockedView1x2");
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
DistMatrix<T,STAR,MC>::View2x1
( DistMatrix<T,STAR,MC>& AT,
  DistMatrix<T,STAR,MC>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View2x1");
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
DistMatrix<T,STAR,MC>::LockedView2x1
( const DistMatrix<T,STAR,MC>& AT,
  const DistMatrix<T,STAR,MC>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::LockedView2x1");
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
DistMatrix<T,STAR,MC>::View2x2
( DistMatrix<T,STAR,MC>& ATL, DistMatrix<T,STAR,MC>& ATR,
  DistMatrix<T,STAR,MC>& ABL, DistMatrix<T,STAR,MC>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View2x2");
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
DistMatrix<T,STAR,MC>::LockedView2x2
( const DistMatrix<T,STAR,MC>& ATL, const DistMatrix<T,STAR,MC>& ATR,
  const DistMatrix<T,STAR,MC>& ABL, const DistMatrix<T,STAR,MC>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::LockedView2x2");
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
DistMatrix<T,STAR,MC>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    this->localMatrix_.ResizeTo
    ( height, LocalLength(width,this->RowShift(),this->Grid().Height()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
DistMatrix<T,STAR,MC>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that
    // row within each process column
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();

    T u;
    if( g.MCRank() == ownerRow )
    {
        const int jLoc = (j-this->RowShift()) / g.Height();
        u = this->GetLocalEntry(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRow, g.MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,STAR,MC>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const int jLoc = (j-this->RowShift()) / g.Height();
        this->SetLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MC>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const int jLoc = (j-this->RowShift()) / g.Height();
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
DistMatrix<T,STAR,MC>::MakeTrapezoidal
( Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int r = this->Grid().Height();
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
            int j = rowShift + jLocal*r;
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
            const int j = rowShift + jLocal*r;
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
DistMatrix<T,STAR,MC>::ScaleTrapezoidal
( T alpha, Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int r = this->Grid().Height();
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
            int j = rowShift + jLocal*r;
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
            int j = rowShift + jLocal*r;
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
DistMatrix<T,STAR,MC>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int r = this->Grid().Height();
    const int rowShift = this->RowShift();

    this->SetToZero();

    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*r;
        if( j < height )
            thisLocalBuffer[j+jLocal*thisLDim] = 1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MC>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const elemental::Grid& g = this->Grid();
    const int height     = this->Height();    
    const int localWidth = this->LocalWidth();
    const int bufSize    = height*localWidth;

    this->auxMemory_.Require( bufSize );

    // Create a random matrix on process column 0, then broadcast
    T* buffer = this->auxMemory_.Buffer();
    if( g.MRRank() == 0 )
    {
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            for( int i=0; i<height; ++i )
                buffer[i+jLocal*height] = SampleUnitBall<T>();
    }
    mpi::Broadcast( buffer, bufSize, 0, g.MRComm() );

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
DistMatrix<T,STAR,MC>::SumOverRow()
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::SumOverRow");
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

    // AllReduce sum
    mpi::AllReduce
    ( sendBuf, recvBuf, localSize, mpi::SUM, this->Grid().MRComm() );

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
DistMatrix<T,STAR,MC>::AdjointFrom( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC]::AdjointFrom");
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
            this->rowAlignment_ = A.ColAlignment() % g.Height();
            this->rowShift_ = 
                Shift( g.MCRank(), this->RowAlignment(), g.Height() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->RowAlignment() == A.ColAlignment() % g.Height() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

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
          gatheredData, portionSize, g.MRComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShiftOfA = RawShift( row+k*r, colAlignmentOfA, p );
            const int rowOffset = (colShiftOfA-rowShift) / r;
            const int localWidth = RawLocalLength( width, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [* ,MC]::AdjointFrom." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int rank = g.VCRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowShift = this->RowShift();

        const int sendRank = (rank+p+rowAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-rowAlignment) % p;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

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
          firstBuffer,  portionSize, recvRank, 0, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MRComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int colShiftOfA = RawShift(row+r*k,rowAlignment,p);
            const int rowOffset = (colShiftOfA-rowShift) / r;
            const int localWidth = RawLocalLength( width, colShiftOfA, p );
            
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
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
DistMatrix<T,STAR,MC>::TransposeFrom
( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC]::TransposeFrom");
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
            this->rowAlignment_ = A.ColAlignment() % g.Height();
            this->rowShift_ = 
                Shift( g.MCRank(), this->RowAlignment(), g.Height() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->RowAlignment() == A.ColAlignment() % g.Height() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

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
          gatheredData, portionSize, g.MRComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShiftOfA = RawShift( row+k*r, colAlignmentOfA, p );
            const int rowOffset = (colShiftOfA-rowShift) / r;
            const int localWidth = RawLocalLength( width, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [* ,MC]::TransposeFrom." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int rank = g.VCRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowShift = this->RowShift();

        const int sendRank = (rank+p+rowAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-rowAlignment) % p;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

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
          firstBuffer,  portionSize, recvRank, 0, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MRComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int colShiftOfA = RawShift(row+r*k,rowAlignment,p);
            const int rowOffset = (colShiftOfA-rowShift) / r;
            const int localWidth = RawLocalLength( width, colShiftOfA, p );
            
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
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
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,STAR,VR> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(g) );
    *A_STAR_VR = A;

    std::auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(true,this->RowAlignment(),g) );
    *A_STAR_VC = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR   = A;

    std::auto_ptr< DistMatrix<T,STAR,VR> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(g) );
    *A_STAR_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    std::auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(true,this->RowAlignment(),g) );
    *A_STAR_VC = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( A.Height() == 1 )
    {
        if( !this->Viewing() )
            this->ResizeTo( 1, A.Width() );

        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int myRow = g.MCRank();
        const int rankCM = g.VCRank();
        const int rankRM = g.VRRank();
        const int rowAlignment = this->RowAlignment();
        const int rowShift = this->RowShift();
        const int rowAlignmentOfA = A.RowAlignment();
        const int rowShiftOfA = A.RowShift();

        const int width = this->Width();
        const int maxLocalVectorWidth = MaxLocalLength(width,p);
        const int portionSize = std::max(maxLocalVectorWidth,mpi::MIN_COLL_MSG);

        const int rowShiftVC = Shift(rankCM,rowAlignment,p);
        const int rowShiftVROfA = Shift(rankRM,rowAlignmentOfA,p);
        const int sendRankCM = (rankCM+(p+rowShiftVROfA-rowShiftVC)) % p;
        const int recvRankRM = (rankRM+(p+rowShiftVC-rowShiftVROfA)) % p;
        const int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        this->auxMemory_.Require( (c+1)*portionSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        // A[* ,VR] <- A[* ,MR]
        {
            const int shift = Shift(rankRM,rowAlignmentOfA,p);
            const int offset = (shift-rowShiftOfA) / c;
            const int thisLocalWidth = LocalLength(width,shift,p);

            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                sendBuf[jLocal] = ALocalBuffer[(offset+jLocal*r)*ALDim];
        }

        // A[* ,VC] <- A[* ,VR]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankCM, 0,
          recvBuf, portionSize, recvRankCM, mpi::ANY_TAG, g.VCComm() );

        // A[* ,MC] <- A[* ,VC]
        mpi::AllGather
        ( recvBuf, portionSize,
          sendBuf, portionSize, g.MRComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &sendBuf[k*portionSize];

            const int shift = RawShift(myRow+r*k,rowAlignment,p);
            const int offset = (shift-rowShift) / r;
            const int thisLocalWidth = RawLocalLength(width,shift,p);

            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                thisLocalBuffer[(offset+jLocal*c)*thisLDim] = data[jLocal];
        }
        this->auxMemory_.Release();
    }
    else
    {
        std::auto_ptr< DistMatrix<T,STAR,VR> > A_STAR_VR
        ( new DistMatrix<T,STAR,VR>(g) );
        *A_STAR_VR = A;

        std::auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
        ( new DistMatrix<T,STAR,VC>(true,this->RowAlignment(),g) );
        *A_STAR_VC = *A_STAR_VR;
        delete A_STAR_VR.release(); // lowers memory highwater

        *this = *A_STAR_VC;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MC] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MC] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Height() == 1 && g.VCRank() == 0 )
    {
        std::cerr << 
          "The vector version of [* ,MC] <- [MR,MC] is not yet written, but"
          " it would only require a modification of the vector version of "
          "[* ,MR] <- [MC,MR]." << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.VCRank() == 0 )
    {
        std::cerr << 
          "The redistribution [* ,MC] <- [MR,MC] potentially causes a large"
          " amount of cache-thrashing. If possible, avoid it. "
          "Unfortunately, the following routines are not yet implemented: \n"
          << "  [MC,* ].(Conjugate)TransposeFrom([MR,MC])" << std::endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            this->rowShift_ = 
                Shift( g.MCRank(), this->RowAlignment(), g.Height() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const int c = g.Width();
        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(height,c);

        const int portionSize = 
            std::max(maxLocalHeightOfA*localWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

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
          gatheredData, portionSize, g.MRComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShift = RawShift( k, colAlignmentOfA, c );
            const int localHeight = RawLocalLength( height, colShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colShift+iLocal*c)+jLocal*thisLDim] =
                        data[iLocal+jLocal*localHeight];
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [* ,MC] <- [MR,MC]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int row = g.MCRank();

        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeightOfA = MaxLocalLength(height,c);
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            std::max(maxLocalHeightOfA*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

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
        ( secondBuffer, portionSize, sendRow, 0,
          firstBuffer,  portionSize, recvRow, mpi::ANY_TAG, g.MCComm() );

        // Use the output of the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MRComm() );

        // Unpack the contents of each member of the process row
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int colShift = RawShift( k, colAlignmentOfA, c );
            const int localHeight = RawLocalLength( height, colShift, c );
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colShift+iLocal*c)+jLocal*thisLDim] =
                        data[iLocal+jLocal*localHeight];
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MR,MC> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,MC]");
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
            std::cerr << "Unaligned [* ,MC] <- [* ,MC]." << std::endl;
#endif
        const int rank = g.MCRank();
        const int r = g.Height();

        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRank = (rank+r+rowAlignmentOfA-rowAlignment) % r;

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
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.MCComm() );

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
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,VR,STAR> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(g) );
    *A_VR_STAR = A;

    std::auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
    ( new DistMatrix<T,MR,MC>(false,true,0,this->RowAlignment(),g) );
    *A_MR_MC = *A_VR_STAR;
    delete A_VR_STAR.release();

    *this = *A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,VC]");
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
            this->rowAlignment_ = A.RowAlignment() % g.Height();
            this->rowShift_ = 
                Shift( g.MCRank(), this->RowAlignment(), g.Height() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() % g.Height() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidthOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalWidthOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

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
          gatheredData, portionSize, g.MRComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int rowAlignmentOfA = A.RowAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShiftOfA = RawShift( row+k*r, rowAlignmentOfA, p );
            const int rowOffset = (rowShiftOfA-rowShift) / r;
            const int localWidth = RawLocalLength( width, rowShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [* ,MC] <- [* ,VC]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int rank = g.VCRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int rowShift = this->RowShift();

        const int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidthOfA = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalWidthOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

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
          firstBuffer,  portionSize, recvRank, 0, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MRComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int rowShiftOfA = RawShift(row+r*k,rowAlignment,p);
            const int rowOffset = (rowShiftOfA-rowShift) / r;
            const int localWidth = RawLocalLength( width, rowShiftOfA, p );
            
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
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
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MR,MC> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,STAR,VC> A_STAR_VC(true,this->RowAlignment(),g);
    *this = A_STAR_VC = A;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MC>&
DistMatrix<T,STAR,MC>::operator=( const DistMatrix<T,STAR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int r = this->Grid().Height();
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
        const T* ACol = &ALocalBuffer[(rowShift+jLocal*r)*ALDim];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        memcpy( thisCol, ACol, localHeight*sizeof(T) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

} // namespace elemental
