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
using namespace std;

template<typename T>
inline void
DistMatrix<T,STAR,VR>::PrintBase
( ostream& os, const string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VRRank() == 0 && msg != "" )
        os << msg << endl;

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

    vector<T> sendBuf(height*width,0);
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
    vector<T> recvBuf;
    if( g.VRRank() == 0 )
        recvBuf.resize(height*width);

    // Sum the contributions and send to the root
    mpi::Reduce
    ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, g.VRComm() );

    if( g.VRRank() == 0 )
    {
        // Print the data
        for( int i=0; i<height; ++i )
        {
            for( int j=0; j<width; ++j )
                os << WrapScalar(recvBuf[i+j*height]) << " ";
            os << "\n";
        }
        os << endl;
    }
    mpi::Barrier( g.VRComm() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,VR>::Align( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::Align");
    this->AssertFreeRowAlignment();
#endif
    this->AlignRows( rowAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,VR>::AlignRows( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Size() )
        throw runtime_error( "Invalid row alignment for [* ,VR]" );
#endif
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift( g.VRRank(), rowAlignment, g.Size() );
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
DistMatrix<T,STAR,VR>::View( DistMatrix<T,STAR,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View");
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
DistMatrix<T,STAR,VR>::View
( int height, int width, int rowAlignment,
  T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift(grid.VRRank(),rowAlignment,grid.Size());
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
DistMatrix<T,STAR,VR>::LockedView( const DistMatrix<T,STAR,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView(A)");
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
DistMatrix<T,STAR,VR>::LockedView
( int height, int width, int rowAlignment,
  const T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->rowShift_ = Shift(grid.VRRank(),rowAlignment,grid.Size());
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
DistMatrix<T,STAR,VR>::View
( DistMatrix<T,STAR,VR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    {
        const elemental::Grid& g = this->Grid();
        const int rowMajorRank = g.VRRank();
        const int size = g.Size();

        this->rowAlignment_ = (A.RowAlignment()+j) % size;
        this->rowShift_ = Shift( rowMajorRank, this->RowAlignment(), size );

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
DistMatrix<T,STAR,VR>::LockedView
( const DistMatrix<T,STAR,VR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    {
        const elemental::Grid& g = this->Grid();
        const int rowMajorRank = g.VRRank();
        const int size = g.Size();

        this->rowAlignment_ = (A.RowAlignment()+j) % size;
        this->rowShift_ = Shift( rowMajorRank, this->RowAlignment(), size );

        const int localWidthBefore = LocalLength( j, A.RowShift(), size );
        const int localWidth = LocalLength( width, this->RowShift(), size );

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
DistMatrix<T,STAR,VR>::View1x2
( DistMatrix<T,STAR,VR>& AL, DistMatrix<T,STAR,VR>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View1x2");
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
DistMatrix<T,STAR,VR>::LockedView1x2
( const DistMatrix<T,STAR,VR>& AL, const DistMatrix<T,STAR,VR>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView1x2");
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
DistMatrix<T,STAR,VR>::View2x1
( DistMatrix<T,STAR,VR>& AT,
  DistMatrix<T,STAR,VR>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View2x1");
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
DistMatrix<T,STAR,VR>::LockedView2x1
( const DistMatrix<T,STAR,VR>& AT,
  const DistMatrix<T,STAR,VR>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView2x1");
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
DistMatrix<T,STAR,VR>::View2x2
( DistMatrix<T,STAR,VR>& ATL, DistMatrix<T,STAR,VR>& ATR,
  DistMatrix<T,STAR,VR>& ABL, DistMatrix<T,STAR,VR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View2x2");
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
DistMatrix<T,STAR,VR>::LockedView2x2
( const DistMatrix<T,STAR,VR>& ATL, const DistMatrix<T,STAR,VR>& ATR,
  const DistMatrix<T,STAR,VR>& ABL, const DistMatrix<T,STAR,VR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView2x2");
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
DistMatrix<T,STAR,VR>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    const elemental::Grid& g = this->Grid();
    this->height_ = height;
    this->width_ = width;
    this->localMatrix_.ResizeTo
    ( height, LocalLength(width,this->RowShift(),g.Size()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
DistMatrix<T,STAR,VR>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    T u;
    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocalEntry(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,STAR,VR>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
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
DistMatrix<T,STAR,VR>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
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
DistMatrix<T,STAR,VR>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const elemental::Grid& g = this->Grid();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int p = g.Size();
    const int rowShift = this->RowShift();

    if( shape == LOWER )
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
                int boundary = min( lastZeroRow+1, height );
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
            int firstZeroRow = ( side==LEFT ? max(j-offset+1,0)
                                            : max(j-offset+height-width+1,0) );
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
DistMatrix<T,STAR,VR>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const elemental::Grid& g = this->Grid();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int p = g.Size();
    const int rowShift = this->RowShift();

    if( shape == UPPER )
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
            int boundary = min( lastRow+1, height );
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
            int firstRow = ( side==LEFT ? max(j-offset,0)
                                        : max(j-offset+height-width,0) );
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
DistMatrix<T,STAR,VR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const elemental::Grid& g = this->Grid();
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int p = g.Size();
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
DistMatrix<T,STAR,VR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetToRandom");
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
inline void
DistMatrix<T,STAR,VR>::AdjointFrom( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[*, VR]::AdjointFrom");
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
            this->rowAlignment_ = A.ColAlignment();
            this->rowShift_ = 
                Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->RowAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int rowShift = this->RowShift();
        const int colShiftOfA = A.ColShift();
        const int rowOffset = (rowShift-colShiftOfA) / c;

        const int height = this->Height();
        const int localWidth = this->LocalWidth();

        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            for( int i=0; i<height; ++i )
                thisLocalBuffer[i+jLocal*thisLDim] = 
                    Conj( ALocalBuffer[(rowOffset+jLocal*r)+i*ALDim] );
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,VR]::AdjointFrom" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int col = g.MRRank();
        const int colShiftOfA = A.ColShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const int sendCol = (col+c+(rowAlignment%c)-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-(rowAlignment%c)) % c;
        const int sendRank = sendCol + c*row;

        const int sendRowShift = Shift( sendRank, rowAlignment, p );
        const int sendRowOffset = (sendRowShift-colShiftOfA) / c;

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
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localWidthOfSend; ++jLocal )
            for( int i=0; i<height; ++i )
                sendBuffer[i+jLocal*height] = 
                    Conj( ALocalBuffer[(sendRowOffset+jLocal*r)+i*ALDim] );

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, mpi::ANY_TAG, g.MRComm() );

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
}

template<typename T>
inline void
DistMatrix<T,STAR,VR>::TransposeFrom
( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR]::TransposeFrom");
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
            this->rowAlignment_ = A.ColAlignment();
            this->rowShift_ = 
                Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->RowAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int rowShift = this->RowShift();
        const int colShiftOfA = A.ColShift();
        const int rowOffset = (rowShift-colShiftOfA) / c;

        const int height = this->Height();
        const int localWidth = this->LocalWidth();

        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            for( int i=0; i<height; ++i )
                thisLocalBuffer[i+jLocal*thisLDim] = 
                    ALocalBuffer[(rowOffset+jLocal*r)+i*ALDim];
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,VR]::TransposeFrom" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int col = g.MRRank();
        const int colShiftOfA = A.ColShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const int sendCol = (col+c+(rowAlignment%c)-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-(rowAlignment%c)) % c;
        const int sendRank = sendCol + c*row;

        const int sendRowShift = Shift( sendRank, rowAlignment, p );
        const int sendRowOffset = (sendRowShift-colShiftOfA) / c;

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
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localWidthOfSend; ++jLocal )
            for( int i=0; i<height; ++i )
                sendBuffer[i+jLocal*height] = 
                    ALocalBuffer[(sendRowOffset+jLocal*r)+i*ALDim];

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, mpi::ANY_TAG, g.MRComm() );

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
}

template<typename T>
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MC,MR]");
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
                Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,r);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( 2*r*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[r*portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const int thisRank = col+k*c;
            const int thisRowShift = RawShift(thisRank,rowAlignment,p);
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*r)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, g.MCComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisColShift = RawShift(k,colAlignmentOfA,r);
            const int thisLocalHeight = RawLocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    thisLocalBuffer[(thisColShift+iLocal*r)+jLocal*thisLDim] =
                        data[iLocal+jLocal*thisLocalHeight];
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,VR] <- [MC,MR]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendCol = (col+c+(rowAlignment%c)-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-(rowAlignment%c)) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,r);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( 2*r*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[r*portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const int thisRank = sendCol+k*c;
            const int thisRowShift = RawShift(thisRank,rowAlignment,p);
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*r)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
            }
        }

        // AllToAll to gather all of the unaligned [*,VR] data into firstBuffer
        mpi::AllToAll
        ( secondBuffer, portionSize,
          firstBuffer,  portionSize, g.MCComm() );

        // SendRecv: properly align the [*,VR] via a trade in the column
        mpi::SendRecv
        ( firstBuffer,  portionSize, sendCol, 0,
          secondBuffer, portionSize, recvCol, mpi::ANY_TAG, g.MRComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisColShift = RawShift(k,colAlignmentOfA,r);
            const int thisLocalHeight = RawLocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    thisLocalBuffer[(thisColShift+iLocal*r)+jLocal*thisLDim] =
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
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,MR]");
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
                Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int rowShift = this->RowShift();
        const int rowShiftOfA = A.RowShift();
        const int rowOffset = (rowShift-rowShiftOfA) / c;

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
            const T* ACol = &ALocalBuffer[(rowOffset+jLocal*r)*ALDim];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            memcpy( thisCol, ACol, height*sizeof(T) );
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,VR] <- [* ,MR]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int col = g.MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const int sendCol = (col+c+(rowAlignment%c)-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-(rowAlignment%c)) % c;
        const int sendRank = sendCol + c*row;

        const int sendRowShift = Shift( sendRank, rowAlignment, p );
        const int sendRowOffset = (sendRowShift-rowShiftOfA) / c;

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
            const T* ACol = &ALocalBuffer[(sendRowOffset+jLocal*r)*ALDim];
            T* sendBufferCol = &sendBuffer[jLocal*height];
            memcpy( sendBufferCol, ACol, height*sizeof(T) );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, mpi::ANY_TAG, g.MRComm() );

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
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,VR] = [MD,* ] is not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,VR] = [* ,MD] is not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,STAR,VC> A_STAR_VC(g);

    A_STAR_VC = A;
    *this = A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
    ( new DistMatrix<T,MR,MC>(g) );
    *A_MR_MC = A;

    auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(g) );
    *A_STAR_VC = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,STAR,VC> A_STAR_VC(g);

    A_STAR_VC = A;
    *this = A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int localWidthOfA = A.LocalWidth();

    const int sendSize = height * localWidthOfA;
    const int recvSize = height * localWidth;

    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();
    const int rankCM = g.VCRank();
    const int rankRM = g.VRRank();

    const int rowShift = this->RowShift();
    const int rowShiftOfA = A.RowShift();

    // Compute which rowmajor rank has the rowShift equal to our rowShiftOfA
    const int sendRankRM = (rankRM+(p+rowShiftOfA-rowShift)) % p;

    // Compute which rowmajor rank has the A rowShift that we need
    const int recvRankCM = (rankCM+(p+rowShift-rowShiftOfA)) % p;
    const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

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
    ( sendBuffer, sendSize, sendRankRM, 0,
      recvBuffer, recvSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

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
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
    ( new DistMatrix<T,MR,MC>(g) );
    *A_MR_MC = A;

    auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(g) );
    *A_STAR_VC = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
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
        const elemental::Grid& g = this->Grid();
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,VR] <- [* ,VR]." << endl;
#endif
        const int rank = g.VRRank();
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
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.VRComm() );

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
inline const DistMatrix<T,STAR,VR>&
DistMatrix<T,STAR,VR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int p = this->Grid().Size();
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
        const T* ACol = &ALocalBuffer[(rowShift+jLocal*p)*ALDim];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        memcpy( thisCol, ACol, localHeight*sizeof(T) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline void
DistMatrix<T,STAR,VR>::SumScatterFrom
( const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SumScatterFrom( [* ,MR] )");
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
                Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();
        const int rowAlignment = this->RowAlignment();
        const int rowShiftOfA = A.RowShift();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int maxLocalWidth = MaxLocalLength( width, p );

        const int recvSize = max(height*maxLocalWidth,mpi::MIN_COLL_MSG);
        const int sendSize = r*recvSize;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        vector<int> recvSizes(r);
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisRank = col+k*c;
            const int thisRowShift = RawShift( thisRank, rowAlignment, p );
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const int thisLocalWidth = RawLocalLength( width, thisRowShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*r)*ALDim];
                T* dataCol = &data[jLocal*height];
                memcpy( dataCol, ACol, height*sizeof(T) );
            }
        }

        // Reduce-scatter over each process column
        mpi::ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.MCComm() );

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
        throw logic_error
              ( "Unaligned [* ,VR]::ReduceScatterFrom( [* ,MR] ) is not "
                "yet implemented." );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,VR>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SumScatterUpdate( [* ,MR] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();
        const int rowAlignment = this->RowAlignment();
        const int rowShiftOfA = A.RowShift();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int maxLocalWidth = MaxLocalLength( width, p );

        const int recvSize = max(height*maxLocalWidth,mpi::MIN_COLL_MSG);
        const int sendSize = r*recvSize;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        vector<int> recvSizes(r);
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisRank = col+k*c;
            const int thisRowShift = RawShift( thisRank, rowAlignment, p );
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const int thisLocalWidth = RawLocalLength( width, thisRowShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*r)*ALDim];
                T* dataCol = &data[jLocal*height];
                memcpy( dataCol, ACol, height*sizeof(T) );
            }
        }

        // Reduce-scatter over each process column
        mpi::ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.MCComm() );

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
        throw logic_error
              ( "Unaligned [* ,VR]::ReduceScatterUpdate( [* ,MR] ) is not "
                "yet implemented." );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
