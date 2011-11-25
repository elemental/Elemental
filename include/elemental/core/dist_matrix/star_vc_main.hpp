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
DistMatrix<T,STAR,VC>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;

    const int height     = this->Height();
    const int width      = this->Width();
    const int localWidth = this->LocalWidth();
    const int p          = g.Size();
    const int rowShift   = this->RowShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    std::vector<T> sendBuf(height*width,0);
    const T* thisLocalBuffer = this->LockedLocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int i=0; i<height; ++i )
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            sendBuf[i+(rowShift+jLocal*p)*height] = 
                thisLocalBuffer[i+jLocal*thisLDim];

    // If we are the root, allocate a receive buffer
    std::vector<T> recvBuf;
    if( g.VCRank() == 0 )
        recvBuf.resize( height*width );

    // Sum the contributions and send to the root
    mpi::Reduce
    ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, g.VCComm() );

    if( g.VCRank() == 0 )
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
    mpi::Barrier( g.VCComm() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,VC>::Align( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Align");
    this->AssertFreeRowAlignment();
#endif
    this->AlignRows( rowAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,VC>::AlignRows( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Size() )
        throw std::runtime_error("Invalid row alignment for [* ,VC]");
#endif
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift( g.VCRank(), rowAlignment, g.Size() );
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
DistMatrix<T,STAR,VC>::View( DistMatrix<T,STAR,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View");
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
DistMatrix<T,STAR,VC>::View
( int height, int width, int rowAlignment,
  T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift(grid.VCRank(),rowAlignment,grid.Size());
    const int localWidth = LocalLength(width,this->rowShift_,grid.Size());
    this->localMatrix_.View( height, localWidth, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,VC>::LockedView( const DistMatrix<T,STAR,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView");
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
DistMatrix<T,STAR,VC>::LockedView
( int height, int width, int rowAlignment,
  const T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift(grid.VCRank(),rowAlignment,grid.Size());
    const int localWidth = LocalLength(width,this->rowShift_,grid.Size());
    this->localMatrix_.LockedView( height, localWidth, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}
template<typename T>
inline void
DistMatrix<T,STAR,VC>::View
( DistMatrix<T,STAR,VC>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    {
        const elemental::Grid& g = this->Grid();
        const int colMajorRank = g.VCRank();
        const int size = g.Size();

        this->rowAlignment_ = (A.RowAlignment()+j) % size;
        this->rowShift_ = Shift( colMajorRank, this->RowAlignment(), size );

        const int localWidthBefore = LocalLength( j, A.RowShift(), size );
        const int localWidth = LocalLength( width, this->RowShift(), size );

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
DistMatrix<T,STAR,VC>::LockedView
( const DistMatrix<T,STAR,VC>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    {
        const elemental::Grid& g = this->Grid();
        const int colMajorRank = g.VCRank();
        const int size = g.Size();

        this->rowAlignment_ = (A.RowAlignment()+j) % size;
        this->rowShift_ = Shift( colMajorRank, this->RowAlignment(), size );

        const int localWidth = LocalLength( width, this->RowShift(), size );

        this->localMatrix_.LockedView
        ( A.LockedLocalMatrix(), i, A.RowShift(), height, localWidth );
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,VC>::View1x2
( DistMatrix<T,STAR,VC>& AL, DistMatrix<T,STAR,VC>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View1x2");
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
DistMatrix<T,STAR,VC>::LockedView1x2
( const DistMatrix<T,STAR,VC>& AL, const DistMatrix<T,STAR,VC>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView1x2");
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
DistMatrix<T,STAR,VC>::View2x1
( DistMatrix<T,STAR,VC>& AT,
  DistMatrix<T,STAR,VC>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View2x1");
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
    this->localMatrix_.View2x1( AT.LocalMatrix(), AB.LocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,VC>::LockedView2x1
( const DistMatrix<T,STAR,VC>& AT,
  const DistMatrix<T,STAR,VC>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView2x1");
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
DistMatrix<T,STAR,VC>::View2x2
( DistMatrix<T,STAR,VC>& ATL, DistMatrix<T,STAR,VC>& ATR,
  DistMatrix<T,STAR,VC>& ABL, DistMatrix<T,STAR,VC>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View2x2");
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
DistMatrix<T,STAR,VC>::LockedView2x2
( const DistMatrix<T,STAR,VC>& ATL, const DistMatrix<T,STAR,VC>& ATR,
  const DistMatrix<T,STAR,VC>& ABL, const DistMatrix<T,STAR,VC>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView2x2");
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
DistMatrix<T,STAR,VC>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    const elemental::Grid& g = this->Grid();
    this->height_ = height;
    this->width_  = width;
    this->localMatrix_.ResizeTo
    ( height, LocalLength( width, this->RowShift(), g.Size() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
DistMatrix<T,STAR,VC>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocalEntry(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );
    
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,STAR,VC>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.Size();
        this->SetLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,VC>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.Size();
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
DistMatrix<T,STAR,VC>::MakeTrapezoidal
( Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int p = this->Grid().Size();
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
            int j = rowShift + jLocal*p;
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
            int j = rowShift + jLocal*p;
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
DistMatrix<T,STAR,VC>::ScaleTrapezoidal
( T alpha, Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int p = this->Grid().Size();
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
            int j = rowShift + jLocal*p;
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
            int j = rowShift + jLocal*p;
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
DistMatrix<T,STAR,VC>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int p = this->Grid().Size();
    const int rowShift = this->RowShift();

    this->SetToZero();

    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*p;
        if( j < height )
            thisLocalBuffer[j+jLocal*thisLDim] = 1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,VC>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<height; ++i )
            this->SetLocalEntry(i,j,SampleUnitBall<T>());
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR> A_STAR_VR(g);

    A_STAR_VR = A;
    *this = A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR = A;

    std::auto_ptr< DistMatrix<T,STAR,VR> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(g) );
    *A_STAR_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR> A_STAR_VR(g);

    A_STAR_VR = A;
    *this = A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,VC] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,VC] = [* ,MD] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MR,MC]");
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
            this->rowShift_ = 
                Shift( g.VCRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,c);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( 2*c*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[c*portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const int thisRank = row+k*r;
            const int thisRowShift = RawShift(thisRank,rowAlignment,p);
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, g.MRComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisColShift = RawShift(k,colAlignmentOfA,c);
            const int thisLocalHeight = RawLocalLength(height,thisColShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    thisLocalBuffer[(thisColShift+iLocal*c)+jLocal*thisLDim] =
                        data[iLocal+jLocal*thisLocalHeight];
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [* ,VC] <- [MR,MC]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+(rowAlignment%r)-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-(rowAlignment%r)) % r;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,c);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( 2*c*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[c*portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const int thisRank = sendRow+k*r;
            const int thisRowShift = RawShift(thisRank,rowAlignment,p);
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
            }
        }

        // AllToAll to gather all of the unaligned [*,VC] data into firstBuffer
        mpi::AllToAll
        ( secondBuffer, portionSize,
          firstBuffer,  portionSize, g.MRComm() );

        // SendRecv: properly align the [*,VC] via a trade in the column
        mpi::SendRecv
        ( firstBuffer,  portionSize, sendRow, 0,
          secondBuffer, portionSize, recvRow, mpi::ANY_TAG, g.MCComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisColShift = RawShift(k,colAlignmentOfA,c);
            const int thisLocalHeight = RawLocalLength(height,thisColShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    thisLocalBuffer[(thisColShift+iLocal*c)+jLocal*thisLDim] =
                        data[iLocal+jLocal*thisLocalHeight];
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MR,* ]");
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
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MC]");
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
            this->rowShift_ = 
                Shift( g.VCRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int rowShift = this->RowShift();
        const int rowShiftOfA = A.RowShift();
        const int rowOffset = (rowShift-rowShiftOfA) / r;

        const int height = this->Height();
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
            const T* ACol = &ALocalBuffer[(rowOffset+jLocal*c)*ALDim];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            memcpy( thisCol, ACol, height*sizeof(T) );
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [* ,VC] <- [* ,MC]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int col = g.MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        // We will SendRecv A[*,VC] within our process column to fix alignments.
        const int sendRow = (row+r+(rowAlignment%r)-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-(rowAlignment%r)) % r;
        const int sendRank = sendRow + r*col;

        const int sendRowShift = Shift( sendRank, rowAlignment, p );
        const int sendRowOffset = (sendRowShift-rowShiftOfA) / r;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localWidthOfSend = LocalLength(width,sendRowShift,p);

        const int sendSize = height * localWidthOfSend;
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
        for( int jLocal=0; jLocal<localWidthOfSend; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[(sendRowOffset+jLocal*c)*ALDim];
            T* sendBufferCol = &sendBuffer[jLocal*height];
            memcpy( sendBufferCol, ACol, height*sizeof(T) );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendRow, 0,
          recvBuffer, recvSize, recvRow, mpi::ANY_TAG, g.MCComm() );

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
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR = A;

    std::auto_ptr< DistMatrix<T,STAR,VR> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(g) );
    *A_STAR_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,VC]");
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
            std::cerr << "Unaligned [* ,VC] <- [* ,VC]." << std::endl;
#endif
        const int rank = g.VCRank();
        const int p = g.Size();

        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

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
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.VCComm() );

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
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [VR,* ]");
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
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int localWidthOfA = A.LocalWidth();
    
    const int sendSize = height * localWidthOfA;
    const int recvSize = height * localWidth;

    const elemental::Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();
    const int rankCM = g.VCRank();
    const int rankRM = g.VRRank(); 

    const int rowShift = this->RowShift();
    const int rowShiftOfA = A.RowShift();

    // Compute which colmajor rank has the rowShift equal to our rowShiftOfA
    const int sendRankCM = (rankCM+(p+rowShiftOfA-rowShift)) % p;

    // Compute which colmajor rank has the A rowShift that we need
    const int recvRankRM = (rankRM+(p+rowShift-rowShiftOfA)) % p;
    const int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

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
    ( sendBuffer, sendSize, sendRankCM, 0,
      recvBuffer, recvSize, recvRankCM, mpi::ANY_TAG, g.VCComm() );

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
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VC>&
DistMatrix<T,STAR,VC>::operator=( const DistMatrix<T,STAR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int p = this->Grid().Size();
    const int rowShift = this->RowShift();

    const int height = this->Height();
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
        const T* ACol = &ALocalBuffer[(rowShift+jLocal*p)*ALDim];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        memcpy( thisCol, ACol, height*sizeof(T) );
    }

#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline void
DistMatrix<T,STAR,VC>::SumScatterFrom( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SumScatterFrom( [* ,MC] )");
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
            this->rowShift_ = 
                Shift( g.VCRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int row = g.MCRank();
        const int rowAlignment = this->RowAlignment();
        const int rowShiftOfA = A.RowShift();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int maxLocalWidth = MaxLocalLength( width, p );

        const int recvSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);
        const int sendSize = c*recvSize;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        std::vector<int> recvSizes(c);
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisRank = row+k*r;
            const int thisRowShift = RawShift( thisRank, rowAlignment, p );
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const int thisLocalWidth = RawLocalLength( width, thisRowShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*height];
                memcpy( dataCol, ACol, height*sizeof(T) );
            }
        }

        // Reduce-scatter over each process row
        mpi::ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.MRComm() );

        // Unpack our received data
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
    else
    {
        throw std::logic_error
        ("Unaligned [* ,VC]::ReduceScatterFrom( [* ,MC] ) not yet implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,VC>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SumScatterUpdate( [* ,MC] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int row = g.MCRank();
        const int rowAlignment = this->RowAlignment();
        const int rowShiftOfA = A.RowShift();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int maxLocalWidth = MaxLocalLength( width, p );

        const int recvSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);
        const int sendSize = c*recvSize;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        std::vector<int> recvSizes(c);
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisRank = row+k*r;
            const int thisRowShift = RawShift( thisRank, rowAlignment, p );
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const int thisLocalWidth = RawLocalLength( width, thisRowShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*height];
                memcpy( dataCol, ACol, height*sizeof(T) );
            }
        }

        // Reduce-scatter over each process row
        mpi::ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.MRComm() );

        // Unpack our received data
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufferCol = &recvBuffer[jLocal*height];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            for( int i=0; i<height; ++i )
                thisCol[i] += alpha*recvBufferCol[i];
        }
        this->auxMemory_.Release();
    }
    else
    {
        throw std::logic_error
        ("Unaligned [* ,VC]::ReduceScatterUpdate( [* ,MC] ) not implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
