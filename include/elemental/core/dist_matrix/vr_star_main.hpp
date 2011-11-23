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
DistMatrix<T,VR,STAR>::PrintBase
( ostream& os, const string msg ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VRRank() == 0 && msg != "" )
        os << msg << endl;

    const int height      = this->Height();
    const int width       = this->Width();
    const int localHeight = this->LocalHeight();
    const int p           = g.Size();
    const int colShift    = this->ColShift();

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
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        for( int j=0; j<width; ++j )
            sendBuf[colShift+iLocal*p+j*height] = 
                thisLocalBuffer[iLocal+j*thisLDim];

    // If we are the root, allocate a receive buffer
    vector<T> recvBuf;
    if( g.VRRank() == 0 )
        recvBuf.resize( height*width );

    // Sum the contributions and send to the root
    mpi::Reduce
    ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, g.VCComm() );

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
    mpi::Barrier( g.VCComm() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::Align( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Align");
    this->AssertFreeColAlignment();
#endif
    this->AlignCols( colAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::AlignCols( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Size() )
        throw std::runtime_error( "Invalid column alignment for [VR,* ]" );
#endif
    this->colAlignment_ = colAlignment;
    this->colShift_ = Shift( g.VRRank(), colAlignment, g.Size() );
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
DistMatrix<T,VR,STAR>::View( DistMatrix<T,VR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View(A)");
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
DistMatrix<T,VR,STAR>::View
( int height, int width, int colAlignment,
  T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->colShift_ = Shift(grid.VRRank(),colAlignment,grid.Size());
    const int localHeight = LocalLength(height,this->colShift_,grid.Size());
    this->localMatrix_.View( localHeight, width, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::LockedView( const DistMatrix<T,VR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView(A)");
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
DistMatrix<T,VR,STAR>::LockedView
( int height, int width, int colAlignment,
  const T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->colShift_ = Shift(grid.VRRank(),colAlignment,grid.Size());
    const int localHeight = LocalLength(height,this->colShift_,grid.Size());
    this->localMatrix_.LockedView( localHeight, width, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::View
( DistMatrix<T,VR,STAR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View(A,i,j,height,width)");
    this->AssertFreeColAlignment();
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

        this->colAlignment_ = (A.ColAlignment()+i) % size;
        this->colShift_ = Shift( rowMajorRank, this->ColAlignment(), size );

        const int localHeightBefore = LocalLength( i, A.ColShift(), size );
        const int localHeight = LocalLength( height, this->ColShift(), size );

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
DistMatrix<T,VR,STAR>::LockedView
( const DistMatrix<T,VR,STAR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView");
    this->AssertFreeColAlignment();
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

        this->colAlignment_ = (A.ColAlignment()+i) % size;
        this->colShift_ = Shift( rowMajorRank, this->ColAlignment(), size );

        const int localHeightBefore = LocalLength( i, A.ColShift(), size );
        const int localHeight = LocalLength( height, this->ColShift(), size );

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
DistMatrix<T,VR,STAR>::View1x2
( DistMatrix<T,VR,STAR>& AL, DistMatrix<T,VR,STAR>& AR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View1x2");
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
DistMatrix<T,VR,STAR>::LockedView1x2
( const DistMatrix<T,VR,STAR>& AL, const DistMatrix<T,VR,STAR>& AR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView1x2");
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
DistMatrix<T,VR,STAR>::View2x1
( DistMatrix<T,VR,STAR>& AT,
  DistMatrix<T,VR,STAR>& AB )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View2x1");
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
    this->localMatrix_.View2x1( AT.LocalMatrix(), AB.LocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::LockedView2x1
( const DistMatrix<T,VR,STAR>& AT,
  const DistMatrix<T,VR,STAR>& AB )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView2x1");
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
DistMatrix<T,VR,STAR>::View2x2
( DistMatrix<T,VR,STAR>& ATL, DistMatrix<T,VR,STAR>& ATR,
  DistMatrix<T,VR,STAR>& ABL, DistMatrix<T,VR,STAR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ABL.Width() + ABR.Width();
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
DistMatrix<T,VR,STAR>::LockedView2x2
( const DistMatrix<T,VR,STAR>& ATL, const DistMatrix<T,VR,STAR>& ATR,
  const DistMatrix<T,VR,STAR>& ABL, const DistMatrix<T,VR,STAR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ABL.Width() + ABR.Width();
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
DistMatrix<T,VR,STAR>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    const elemental::Grid& g = this->Grid();
    this->height_ = height;
    this->width_  = width;
    this->localMatrix_.ResizeTo
    ( LocalLength(height,this->ColShift(),g.Size()) ,width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
DistMatrix<T,VR,STAR>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    T u;
    if( g.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        u = this->GetLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        this->SetLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
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
DistMatrix<T,VR,STAR>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->Grid().Size();
    const int colShift = this->ColShift();

    if( shape == LOWER )
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
                int boundary = min( lastZeroRow+1, height );
                int numZeroRows = RawLocalLength( boundary, colShift, p );
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
            int firstZeroRow = ( side==LEFT ? max(j-offset+1,0)
                                            : max(j-offset+height-width+1,0) );
            int numNonzeroRows = RawLocalLength(firstZeroRow,colShift,p);
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
DistMatrix<T,VR,STAR>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->Grid().Size();
    const int colShift = this->ColShift();

    if( shape == UPPER )
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            int boundary = min( lastRow+1, height );
            int numRows = RawLocalLength( boundary, colShift, p );
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
            int firstRow = ( side==LEFT ? max(j-offset,0)
                                        : max(j-offset+height-width,0) );
            int numZeroRows = RawLocalLength(firstRow,colShift,p);
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
DistMatrix<T,VR,STAR>::SetToIdentity()
{   
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->Grid().Size();
    const int colShift = this->ColShift();

    this->SetToZero();

    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*p;
        if( i < width )
            thisLocalBuffer[iLocal+i*thisLDim] = 1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    for( int j=0; j<width; ++j )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            this->SetLocalEntry(iLocal,j,SampleUnitBall<T>());
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR> A_VC_STAR(g);

    A_VC_STAR = A;
    *this = A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR> A_VC_STAR(g);

    A_VC_STAR = A;
    *this = A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR = A;

    auto_ptr< DistMatrix<T,VC,STAR> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[VR,* ] = [MD,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[VR,* ] = [* ,MD] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MR,MC]");
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
                Shift( g.VRRank(), this->ColAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,r);
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
            const int thisColShift = RawShift(thisRank,colAlignment,p);
            const int thisColOffset = (thisColShift-colShiftOfA) / c;
            const int thisLocalHeight = RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColOffset+iLocal*r)+jLocal*ALDim];
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

            const int thisRowShift = RawShift(k,rowAlignmentOfA,r);
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(thisRowShift+jLocal*r)*thisLDim];
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [VR,* ] <- [MR,MC]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendCol = (col+c+(colAlignment%c)-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-(colAlignment%c)) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,r);
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
            const int thisColShift = RawShift(thisRank,colAlignment,p);
            const int thisColOffset = (thisColShift-colShiftOfA) / c;
            const int thisLocalHeight = RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColOffset+iLocal*r)+jLocal*ALDim];
        }

        // AllToAll to gather all of the unaligned [VR,*] data into firstBuffer
        mpi::AllToAll
        ( secondBuffer, portionSize,
          firstBuffer,  portionSize, g.MCComm() );

        // SendRecv: properly align the [VR,*] via a trade in the row
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

            const int thisRowShift = RawShift(k,rowAlignmentOfA,r);
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(thisRowShift+jLocal*r)*thisLDim];
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
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MR,* ]");
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
                Shift( g.VRRank(), this->ColAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int colShift = this->ColShift();
        const int colShiftOfA = A.ColShift();
        const int colOffset = (colShift-colShiftOfA) / c;

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
                    ALocalBuffer[(colOffset+iLocal*r)+j*ALDim];
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [VR,* ] <- [MR,* ]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int col = g.MRRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[VR,*] within our process row to fix alignments.
        const int sendCol = (col+c+(colAlignment%c)-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-(colAlignment%c)) % c;
        const int sendRank = sendCol + c*row;

        const int sendColShift = Shift( sendRank, colAlignment, p );
        const int sendColOffset = (sendColShift-colShiftOfA) / c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfSend = LocalLength(height,sendColShift,p);
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = maxLocalHeight * width;

        this->auxMemory_.Require( 2*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<width; ++j )
            for( int iLocal=0; iLocal<localHeightOfSend; ++iLocal )
                sendBuffer[iLocal+j*localHeightOfSend] = 
                    ALocalBuffer[(sendColOffset+iLocal*r)+j*ALDim];

        // Communicate
        mpi::SendRecv
        ( sendBuffer, portionSize, sendCol, 0,
          recvBuffer, portionSize, recvCol, mpi::ANY_TAG, g.MRComm() );

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
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MC]");
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
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    
    const elemental::Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();
    const int rankCM = g.VCRank();
    const int rankRM = g.VRRank();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localHeightOfA = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength(height,p);

    const int portionSize = maxLocalHeight * width;

    const int colShift = this->ColShift();
    const int colShiftOfA = A.ColShift();

    // Compute which rowmajor rank has the colShift equal to our colShiftOfA
    const int sendRankRM = (rankRM+(p+colShiftOfA-colShift)) % p;

    // Compute which rowmajor rank has the A colShift that we need
    const int recvRankCM = (rankCM+(p+colShift-colShiftOfA)) % p;
    const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

    this->auxMemory_.Require( 2*portionSize );

    T* buffer = this->auxMemory_.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[portionSize];

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
    ( sendBuffer, portionSize, sendRankRM, 0,
      recvBuffer, portionSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

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
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,VC]");
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
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
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
        const elemental::Grid& g = this->Grid();
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [VR,* ] <- [VR,* ]." << endl;
#endif
        const int rank = g.VRRank();
        const int p = g.Size();

        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

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
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.VRComm() );

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
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR = A;

    auto_ptr< DistMatrix<T,VC,STAR> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VR,STAR>&
DistMatrix<T,VR,STAR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int p = this->Grid().Size();
    const int colShift = this->ColShift();

    const int localHeight = this->LocalHeight();
    const int width = this->Width();

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
                ALocalBuffer[(colShift+iLocal*p)+j*ALDim];
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::SumScatterFrom
( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SumScatterFrom( [MR,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.VCRank() == 0 )
    {
        cerr <<
          "[VR,* ]::SumScatterFrom([MR,* ]) potentially causes a large amount "
          "of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MR,* ] matrix instead." << endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            this->colShift_ = 
                Shift( g.VRRank(), this->ColAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();
        const int colAlignment = this->ColAlignment();
        const int colShiftOfA = A.ColShift();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int maxLocalHeight = MaxLocalLength( height, p );

        const int recvSize = max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
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
            const int thisColShift = RawShift( thisRank, colAlignment, p );
            const int thisColOffset = (thisColShift-colShiftOfA) / r;
            const int thisLocalHeight = 
                RawLocalLength( height, thisColShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+j*thisLocalHeight] = 
                        ALocalBuffer[(thisColOffset+iLocal*r)+j*ALDim];
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
        for( int j=0; j<width; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*localHeight];
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
        this->auxMemory_.Release();
    }
    else
    {
        throw logic_error
              ( "Unaligned [VR,* ]::ReduceScatterFrom( [MR,* ] ) is not "
                "yet implemented." );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::SumScatterFrom
( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SumScatterFrom( [* ,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int p = g.Size();
    const int VRRank = g.VRRank();
    const int colAlignment = this->ColAlignment();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int maxLocalHeight = MaxLocalLength( height, p );

    const int recvSize = max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
    const int sendSize = p*recvSize;

    this->auxMemory_.Require( sendSize + recvSize );

    T* buffer = this->auxMemory_.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack
    vector<int> recvSizes(p);
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<p; ++k )
    {
        T* data = &sendBuffer[k*recvSize];
        recvSizes[k] = recvSize;

        const int thisColShift = RawShift( k, colAlignment, p );
        const int thisLocalHeight = RawLocalLength( height, thisColShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<width; ++j )
            for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                data[iLocal+j*thisLocalHeight] =
                    ALocalBuffer[(thisColShift+iLocal*p)+j*ALDim];
    }

    // Reduce-scatter over each process row
    mpi::ReduceScatter
    ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.VRComm() );

    // Unpack our received data
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::SumScatterUpdate
( T alpha, const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SumScatterUpdate( [MR,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.VCRank() == 0 )
    {
        cerr <<
          "[VR,* ]::SumScatterUpdate([MR,* ]) potentially causes a large amount"
          " of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MR,* ] matrix instead." << endl;
    }
#endif
    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();
        const int colAlignment = this->ColAlignment();
        const int colShiftOfA = A.ColShift();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int maxLocalHeight = MaxLocalLength( height, p );

        const int recvSize = max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
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
            const int thisColShift = RawShift( thisRank, colAlignment, p );
            const int thisColOffset = (thisColShift-colShiftOfA) / c;
            const int thisLocalHeight = 
                RawLocalLength( height, thisColShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+j*thisLocalHeight] = 
                        ALocalBuffer[(thisColOffset+iLocal*r)+j*ALDim];
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
        for( int j=0; j<width; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*localHeight];
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                thisCol[iLocal] += alpha*recvBufferCol[iLocal];
        }
        this->auxMemory_.Release();
    }
    else
    {
        throw logic_error
              ( "Unaligned [VR,* ]::ReduceScatterUpdate( [MR,* ] ) is not "
                "yet implemented." );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VR,STAR>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SumScatterUpdate( [* ,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();

    const int p = g.Size();
    const int VRRank = g.VRRank();
    const int colAlignment = this->ColAlignment();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int maxLocalHeight = MaxLocalLength( height, p );

    const int recvSize = max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
    const int sendSize = p*recvSize;

    this->auxMemory_.Require( sendSize + recvSize );

    T* buffer = this->auxMemory_.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack
    vector<int> recvSizes(p);
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<p; ++k )
    {
        T* data = &sendBuffer[k*recvSize];
        recvSizes[k] = recvSize;

        const int thisColShift = RawShift( k, colAlignment, p );
        const int thisLocalHeight = RawLocalLength( height, thisColShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<width; ++j )
            for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                data[iLocal+j*thisLocalHeight] =
                    ALocalBuffer[(thisColShift+iLocal*p)+j*ALDim];
    }

    // Reduce-scatter over each process row
    mpi::ReduceScatter
    ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.VRComm() );

    // Unpack our received data
    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<width; ++j )
    {
        const T* recvBufferCol = &recvBuffer[j*localHeight];
        T* thisCol = &thisLocalBuffer[j*thisLDim];
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            thisCol[iLocal] += alpha*recvBufferCol[iLocal];
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
