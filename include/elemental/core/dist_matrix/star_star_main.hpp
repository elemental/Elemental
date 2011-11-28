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
DistMatrix<T,STAR,STAR>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;

    const int height = this->Height();
    const int width  = this->Width();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( g.InGrid() )
    {
        if( g.VCRank() == 0 )
        {
            for( int i=0; i<height; ++i )
            {
                for( int j=0; j<width; ++j )
                    os << WrapScalar(this->GetLocalEntry(i,j)) << " ";
                os << "\n";
            }
            os << std::endl;
        }
        mpi::Barrier( g.VCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::View( DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View");
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
        this->localMatrix_.View( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::View
( int height, int width, 
  T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View");
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
        this->localMatrix_.View( height, width, buffer, ldim );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::LockedView
( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView");
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
        this->localMatrix_.LockedView( A.LockedLocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::LockedView
( int height, int width, 
  const T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView");
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
        this->localMatrix_.LockedView( height, width, buffer, ldim );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::View
( DistMatrix<T,STAR,STAR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View");
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
        this->localMatrix_.View( A.LocalMatrix(), i, j, height, width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::LockedView
( const DistMatrix<T,STAR,STAR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView");
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.LockedView
        ( A.LockedLocalMatrix(), i, j, height, width );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::View1x2
( DistMatrix<T,STAR,STAR>& AL, DistMatrix<T,STAR,STAR>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View1x2");
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
        this->localMatrix_.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::LockedView1x2
( const DistMatrix<T,STAR,STAR>& AL, const DistMatrix<T,STAR,STAR>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView1x2");
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.LockedView1x2
        ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::View2x1
( DistMatrix<T,STAR,STAR>& AT,
  DistMatrix<T,STAR,STAR>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View2x1");
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.View2x1
        ( AT.LocalMatrix(),
          AB.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::LockedView2x1
( const DistMatrix<T,STAR,STAR>& AT,
  const DistMatrix<T,STAR,STAR>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView2x1");
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.LockedView2x1
        ( AT.LockedLocalMatrix(),
          AB.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::View2x2
( DistMatrix<T,STAR,STAR>& ATL, DistMatrix<T,STAR,STAR>& ATR,
  DistMatrix<T,STAR,STAR>& ABL, DistMatrix<T,STAR,STAR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View2x2");
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.View2x2
        ( ATL.LocalMatrix(), ATR.LocalMatrix(),
          ABL.LocalMatrix(), ABR.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::LockedView2x2
( const DistMatrix<T,STAR,STAR>& ATL, const DistMatrix<T,STAR,STAR>& ATR,
  const DistMatrix<T,STAR,STAR>& ABL, const DistMatrix<T,STAR,STAR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView2x2");
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.LockedView2x2
        ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
          ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Grid().InGrid() )
        this->localMatrix_.ResizeTo( height, width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
DistMatrix<T,STAR,STAR>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    const int viewingSize = mpi::CommSize( this->Grid().ViewingComm() );
    const int owningSize = mpi::GroupSize( this->Grid().OwningGroup() );
    T u;
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        u = this->GetLocalEntry(i,j);
    }
    else
    {
        // Have the root broadcast its data
        if( this->Grid().VCRank() == 0 )
            u = this->GetLocalEntry(i,j);
        mpi::Broadcast
        ( &u, 1, this->Grid().VCToViewingMap(0), 
          this->Grid().ViewingComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    if( this->Grid().InGrid() )
        this->SetLocalEntry(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    if( this->Grid().InGrid() )
    {
        this->UpdateLocalEntry(i,j,u);
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
DistMatrix<T,STAR,STAR>::MakeTrapezoidal
( Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();

    if( this->Grid().InGrid() )
    {
        if( uplo == LOWER )
        {
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int lastZeroRow = 
                    ( side==LEFT ? j-offset-1 : j-offset+height-width-1 );
                if( lastZeroRow >= 0 )
                {
                    int boundary = std::min( lastZeroRow+1, height );
                    T* thisCol = &thisLocalBuffer[j*thisLDim];
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
            for( int j=0; j<width; ++j )
            {
                int firstZeroRow = 
                    ( side==LEFT ? std::max(j-offset+1,0)
                                 : std::max(j-offset+height-width+1,0) );
                if( firstZeroRow < height )
                {
                    T* thisCol = &thisLocalBuffer[firstZeroRow+j*thisLDim];
                    memset( thisCol, 0, (height-firstZeroRow)*sizeof(T) );
                }
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::ScaleTrapezoid
( T alpha, Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::ScaleTrapezoid");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();

    if( this->Grid().InGrid() )
    {
        if( uplo == UPPER )
        {
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int lastRow = 
                    ( side==LEFT ? j-offset : j-offset+height-width );
                int boundary = std::min( lastRow+1, height );
                T* thisCol = &thisLocalBuffer[j*thisLDim];
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
            for( int j=0; j<width; ++j )
            {
                int firstRow = 
                    ( side==LEFT ? std::max(j-offset,0)
                                 : std::max(j-offset+height-width,0) );
                T* thisCol = &thisLocalBuffer[firstRow+j*thisLDim];
                for( int i=0; i<(height-firstRow); ++i )
                    thisCol[i] *= alpha;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();

    if( this->Grid().InGrid() )
    {
        this->SetToZero();

        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<std::min(height,width); ++j )
            thisLocalBuffer[j+j*thisLDim] = 1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    // Create random matrix on process 0 and then broadcast
    const elemental::Grid& g = this->Grid();
    const int height = this->Height();
    const int width = this->Width();
    const int bufSize = height*width;

    if( g.InGrid() )
    {
        this->auxMemory_.Require( bufSize );

        T* buffer = this->auxMemory_.Buffer();
        if( g.VCRank() == 0 )
        {
            for( int j=0; j<width; ++j )
                for( int i=0; i<height; ++i )
                    buffer[i+j*height] = SampleUnitBall<T>();
        }
        mpi::Broadcast( buffer, bufSize, 0, g.VCComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            const T* bufferCol = &buffer[j*height];
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            memcpy( thisCol, bufferCol, height*sizeof(T) );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int r = g.Height();
        const int c = g.Width(); 
        const int p = g.Size();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);
        const int maxLocalWidth = MaxLocalLength(width,c);

        const int portionSize = 
            std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (p+1)*portionSize );

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
            T* originalDataCol = &originalData[jLocal*localHeightOfA];
            memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.VCComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int l=0; l<c; ++l )
        {
            const int rowShift = RawShift( l, rowAlignmentOfA, c );
            const int localWidth = RawLocalLength( width, rowShift, c );

            for( int k=0; k<r; ++k )
            {
                const T* data = &gatheredData[(k+l*r)*portionSize];

                const int colShift = RawShift( k, colAlignmentOfA, r );
                const int localHeight = RawLocalLength( height, colShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                    for( int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisLocalBuffer[(colShift+iLocal*r)+
                                        (rowShift+jLocal*c)*thisLDim] = 
                            data[iLocal+jLocal*localHeight];
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int r = g.Height();
        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeight = MaxLocalLength(height,r);

        const int portionSize = 
            std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

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
            for( int j=0; j<width; ++j )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colShift+iLocal*r)+j*thisLDim] =
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int c = g.Width();
        const int height = this->Height();
        const int width = this->Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength(width,c);

        const int portionSize = 
            std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

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
        const int rowAlignmentOfA = A.RowAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShift = RawShift( k, rowAlignmentOfA, c );
            const int localWidth = RawLocalLength( width, rowShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowShift+jLocal*c)*thisLDim];
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int lcm = g.LCM();
        const int ownerPath = g.DiagPath( A.ColAlignment() );
        const int ownerPathRank = g.DiagPathRank( A.ColAlignment() );

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = A.LocalHeight();
        const int maxLocalHeight = MaxLocalLength( height, lcm );
        const int portionSize = 
            std::max( maxLocalHeight*width, mpi::MIN_COLL_MSG );

        // Since a MD communicator has not been implemented, we will take
        // the suboptimal route of 'rounding up' everyone's contribution over 
        // the VC communicator.
        this->auxMemory_.Require( (p+1)*portionSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        if( A.InDiagonal() )
        {
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                const T* ACol = &ALocalBuffer[j*ALDim];
                T* sendBufCol = &sendBuf[j*localHeight];
                memcpy( sendBufCol, ACol, localHeight*sizeof(T) );
            }
        }

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.VCComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            if( g.DiagPath( k ) == ownerPath )
            {
                const T* data = &recvBuf[k*portionSize];

                const int thisPathRank = g.DiagPathRank( k );
                const int thisColShift = 
                    RawShift( thisPathRank, ownerPathRank, lcm );
                const int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, lcm );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int j=0; j<width; ++j )
                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        thisLocalBuffer[(thisColShift+iLocal*lcm)+j*thisLDim] =
                            data[iLocal+j*thisLocalHeight];
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int lcm = g.LCM();
        const int ownerPath = g.DiagPath( A.RowAlignment() );
        const int ownerPathRank = g.DiagPathRank( A.RowAlignment() );

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength( width, lcm );
        const int portionSize = 
            std::max( height*maxLocalWidth, mpi::MIN_COLL_MSG );

        // Since a MD communicator has not been implemented, we will take
        // the suboptimal route of 'rounding up' everyone's contribution over 
        // the VC communicator.
        this->auxMemory_.Require( (p+1)*portionSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        if( A.InDiagonal() )
        {
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[jLocal*ALDim];
                T* sendBufCol = &sendBuf[jLocal*height];
                memcpy( sendBufCol, ACol, height*sizeof(T) );
            }
        }

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.VCComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            if( g.DiagPath( k ) == ownerPath )
            {
                const T* data = &recvBuf[k*portionSize];

                const int thisPathRank = g.DiagPathRank( k );
                const int thisRowShift = 
                    RawShift( thisPathRank, ownerPathRank, lcm );
                const int thisLocalWidth = 
                    RawLocalLength( width, thisRowShift, lcm );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* dataCol = &data[jLocal*height];
                    T* thisCol = 
                        &thisLocalBuffer[(thisRowShift+jLocal*lcm)*thisLDim];
                    memcpy( thisCol, dataCol, height*sizeof(T) );
                }
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (p+1)*portionSize );

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
            T* originalDataCol = &originalData[jLocal*localHeightOfA];
            memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.VRComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int l=0; l<r; ++l )
        {
            const int rowShift = RawShift( l, rowAlignmentOfA, r );
            const int localWidth = RawLocalLength( width, rowShift, r );

            for( int k=0; k<c; ++k )
            {
                const T* data = &gatheredData[(k+l*c)*portionSize];

                const int colShift = RawShift( k, colAlignmentOfA, c );
                const int localHeight = RawLocalLength( height, colShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                    for( int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisLocalBuffer[(colShift+iLocal*c)+
                                        (rowShift+jLocal*r)*thisLDim] = 
                            data[iLocal+jLocal*localHeight];
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int c = g.Width();
        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeight = MaxLocalLength(height,c);

        const int portionSize = 
            std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

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
        for( int j=0; j<width; ++j )
        {
            const T* ACol = &ALocalBuffer[j*ALDim];
            T* originalDataCol = &originalData[j*localHeightOfA];
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
            for( int j=0; j<width; ++j )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colShift+iLocal*c)+j*thisLDim] = 
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int r = g.Height();
        const int height = this->Height();
        const int width = this->Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

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
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowShift+jLocal*r)*thisLDim];
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,VC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = 
            std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (p+1)*portionSize );

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
          gatheredData, portionSize, g.VCComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShift = RawShift( k, colAlignmentOfA, p );
            const int localHeight = RawLocalLength( height, colShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colShift+iLocal*p)+j*thisLDim] = 
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int height = this->Height();
        const int width = this->Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (p+1)*portionSize );

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
          gatheredData, portionSize, g.VCComm() );

        // Unpack
        const int rowAlignmentOfA = A.RowAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShift = RawShift( k, rowAlignmentOfA, p );
            const int localWidth = RawLocalLength( width, rowShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowShift+jLocal*p)*thisLDim];
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,VR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = 
            std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (p+1)*portionSize );

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
          gatheredData, portionSize, g.VRComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShift = RawShift( k, colAlignmentOfA, p );
            const int localHeight = RawLocalLength( height, colShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colShift+iLocal*p)+j*thisLDim] = 
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int height = this->Height();
        const int width = this->Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength(width,p);

        const int portionSize = 
            std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (p+1)*portionSize );

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
          gatheredData, portionSize, g.VRComm() );

        // Unpack
        const int rowAlignmentOfA = A.RowAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShift = RawShift( k, rowAlignmentOfA, p );
            const int localWidth = RawLocalLength( width, rowShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowShift+jLocal*p)*thisLDim];
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
inline const DistMatrix<T,STAR,STAR>&
DistMatrix<T,STAR,STAR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,* ]");
    this->AssertNotLockedView();
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( this->Grid() == A.Grid() )
    {
        this->localMatrix_ = A.LockedLocalMatrix();
    }
    else
    {
        if( !mpi::CongruentComms( A.Grid().ViewingComm(),
                                  this->Grid().ViewingComm() ) )
            throw std::logic_error
                ("Redistributing between nonmatching grids currently requires"
                 " the viewing communicators to match.");

        // Compute and allocate the amount of required memory
        int requiredMemory = 0;
        if( A.Grid().VCRank() == 0 )
            requiredMemory += A.Height()*A.Width();
        if( this->Grid().InGrid() )
            requiredMemory += A.Height()*A.Width();
        this->auxMemory_.Require( requiredMemory );
        T* buffer = this->auxMemory_.Buffer();
        int offset = 0;
        T* sendBuffer = &buffer[offset];
        if( A.Grid().VCRank() == 0 )
            offset += A.Height()*A.Width();
        T* bcastBuffer = &buffer[offset];

        // Send from the root of A to the root of this matrix's grid
        mpi::Request sendRequest;
        if( A.Grid().VCRank() == 0 )
        {
            for( int j=0; j<A.Width(); ++j ) 
                for( int i=0; i<A.Height(); ++i )
                    sendBuffer[i+j*A.Height()] = A.GetLocalEntry(i,j);
            const int recvViewingRank = this->Grid().VCToViewingMap(0);
            mpi::ISend
            ( sendBuffer, A.Height()*A.Width(), recvViewingRank, 0,
              this->Grid().ViewingComm(), sendRequest );
        }

        // Receive on the root of this matrix's grid and then broadcast
        // over this matrix's owning communicator
        if( this->Grid().InGrid() )
        {
            if( this->Grid().VCRank() == 0 )
            {
                const int sendViewingRank = A.Grid().VCToViewingMap(0);
                mpi::Recv
                ( bcastBuffer, A.Height()*A.Width(), sendViewingRank, 0,
                  this->Grid().ViewingComm() );
            }

            mpi::Broadcast
            ( bcastBuffer, A.Height()*A.Width(), 0, this->Grid().VCComm() );

            for( int j=0; j<A.Width(); ++j )
                for( int i=0; i<A.Height(); ++i )
                    this->SetLocalEntry(i,j,bcastBuffer[i+j*A.Height()]);
        }

        if( A.Grid().VCRank() == 0 )
            mpi::Wait( sendRequest );
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::SumOverCol()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SumOverCol");
    this->AssertNotLockedView();
#endif
    if( this->Grid().InGrid() )
    {
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localSize = localHeight*localWidth;
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

        // Sum
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
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::SumOverRow()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SumOverRow");
    this->AssertNotLockedView();
#endif
    if( this->Grid().InGrid() )
    {
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localSize = localHeight*localWidth;
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

        // Sum
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
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::SumOverGrid()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SumOverGrid");
    this->AssertNotLockedView();
#endif
    if( this->Grid().InGrid() )
    {
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localSize = localHeight*localWidth;
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

        // Sum
        mpi::AllReduce
        ( sendBuf, recvBuf, localSize, mpi::SUM, this->Grid().VCComm() );

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
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
