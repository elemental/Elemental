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
DistMatrix<T,MR,MC>::PrintBase( ostream& os, const string msg ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VCRank() == 0 && msg != "" )
        os << msg << endl;

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Fill the send buffer: zero it then place our entries into their
    // appropriate locations
    vector<T> sendBuf(height*width,0);
    const T* thisLocalBuffer = this->LockedLocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            sendBuf[(colShift+iLocal*c) + (rowShift+jLocal*r)*height] = 
                thisLocalBuffer[iLocal+jLocal*thisLDim];

    // If we are the root, allocate a receive buffer
    vector<T> recvBuf;
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
        os << endl;
    }
    mpi::Barrier( g.VCComm() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::Align( int colAlignment, int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Align");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Width() )
        throw runtime_error("Invalid column alignment for [MR,MC]");
    if( rowAlignment < 0 || rowAlignment >= g.Height() )
        throw runtime_error("Invalid row alignment for [MR,MC]");
#endif
    this->colAlignment_ = colAlignment;
    this->rowAlignment_ = rowAlignment;
    this->colShift_ = Shift( g.MRRank(), colAlignment, g.Width() );
    this->rowShift_ = Shift( g.MCRank(), rowAlignment, g.Height() );
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MR,MC>::AlignCols( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Width() )
        throw runtime_error("Invalid column alignment for [MR,MC]");
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
DistMatrix<T,MR,MC>::AlignRows( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Height() )
        throw runtime_error("Invalid row alignment for [MR,MC]");
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
DistMatrix<T,MR,MC>::View( DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::View");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_  = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->rowAlignment_ = A.RowAlignment();
    this->colShift_     = A.ColShift();
    this->rowShift_     = A.RowShift();
    this->localMatrix_.View( A.LocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::View
( int height, int width, int colAlignment, int rowAlignment,
  T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::View");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->rowAlignment_ = rowAlignment;
    this->colShift_ = Shift(grid.MRRank(),colAlignment,grid.Width());
    this->rowShift_ = Shift(grid.MCRank(),rowAlignment,grid.Height());
    const int localHeight = LocalLength(height,this->colShift_,grid.Width());
    const int localWidth = LocalLength(width,this->rowShift_,grid.Height());
    this->localMatrix_.View( localHeight, localWidth, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::LockedView( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_  = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->rowAlignment_ = A.RowAlignment();
    this->colShift_     = A.ColShift();
    this->rowShift_     = A.RowShift();
    this->localMatrix_.LockedView( A.LockedLocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::LockedView
( int height, int width, int colAlignment, int rowAlignment,
  const T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->rowAlignment_ = rowAlignment;
    this->colShift_ = Shift(grid.MRRank(),colAlignment,grid.Width());
    this->rowShift_ = Shift(grid.MCRank(),rowAlignment,grid.Height());
    const int localHeight = LocalLength(height,this->colShift_,grid.Width());
    const int localWidth = LocalLength(width,this->rowShift_,grid.Height());
    this->localMatrix_.LockedView( localHeight, localWidth, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::View
( DistMatrix<T,MR,MC>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::View");
    this->AssertFreeColAlignment();
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
        const int c   = g.Width();
        const int row = g.MCRank();
        const int col = g.MRRank();

        this->colAlignment_ = (A.ColAlignment()+i) % c;
        this->rowAlignment_ = (A.RowAlignment()+j) % r;
        
        this->colShift_ = Shift( col, this->ColAlignment(), c );
        this->rowShift_ = Shift( row, this->RowAlignment(), r );

        const int localHeightBefore = LocalLength( i, A.ColShift(), c );
        const int localWidthBefore  = LocalLength( j, A.RowShift(), r );

        const int localHeight = LocalLength( height, this->ColShift(), c );
        const int localWidth  = LocalLength( width,  this->RowShift(), r );

        this->localMatrix_.View
        ( A.LocalMatrix(),
          localHeightBefore, localWidthBefore, localHeight, localWidth );
    }
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::LockedView
( const DistMatrix<T,MR,MC>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::LockedView");
    this->AssertFreeColAlignment();
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
        const int c   = g.Width();
        const int row = g.MCRank();
        const int col = g.MRRank();

        this->colAlignment_ = (A.ColAlignment()+i) % c;
        this->rowAlignment_ = (A.RowAlignment()+j) % r;
        
        this->colShift_ = Shift( col, this->ColAlignment(), c );
        this->rowShift_ = Shift( row, this->RowAlignment(), r );

        const int localHeightBefore = LocalLength( i, A.ColShift(), c );
        const int localWidthBefore  = LocalLength( j, A.RowShift(), r );

        const int localHeight = LocalLength( height, this->ColShift(), c );
        const int localWidth  = LocalLength( width,  this->RowShift(), r );

        this->localMatrix_.LockedView
        ( A.LockedLocalMatrix(),
          localHeightBefore, localWidthBefore, localHeight, localWidth );
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::View1x2
( DistMatrix<T,MR,MC>& AL, DistMatrix<T,MR,MC>& AR )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::View1x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_  = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->rowAlignment_ = AL.RowAlignment();
    this->colShift_     = AL.ColShift();
    this->rowShift_     = AL.RowShift();
    this->localMatrix_.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::LockedView1x2
( const DistMatrix<T,MR,MC>& AL, const DistMatrix<T,MR,MC>& AR )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::LockedView1x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_  = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->rowAlignment_ = AL.RowAlignment();
    this->colShift_     = AL.ColShift();
    this->rowShift_     = AL.RowShift();
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
DistMatrix<T,MR,MC>::View2x1
( DistMatrix<T,MR,MC>& AT,
  DistMatrix<T,MR,MC>& AB )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::View2x1");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_  = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->rowAlignment_ = AT.RowAlignment();
    this->colShift_     = AT.ColShift();
    this->rowShift_     = AT.RowShift();
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
DistMatrix<T,MR,MC>::LockedView2x1
( const DistMatrix<T,MR,MC>& AT,
  const DistMatrix<T,MR,MC>& AB )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::LockedView2x1");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_  = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->rowAlignment_ = AT.RowAlignment();
    this->colShift_     = AT.ColShift();
    this->rowShift_     = AT.RowShift();
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
DistMatrix<T,MR,MC>::View2x2
( DistMatrix<T,MR,MC>& ATL, DistMatrix<T,MR,MC>& ATR,
  DistMatrix<T,MR,MC>& ABL, DistMatrix<T,MR,MC>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::View2x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_  = ATL.Width() + ATR.Width();
    this->colAlignment_ = ATL.ColAlignment();
    this->rowAlignment_ = ATL.RowAlignment();
    this->colShift_     = ATL.ColShift();
    this->rowShift_     = ATL.RowShift();
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
DistMatrix<T,MR,MC>::LockedView2x2
( const DistMatrix<T,MR,MC>& ATL, const DistMatrix<T,MR,MC>& ATR,
  const DistMatrix<T,MR,MC>& ABL, const DistMatrix<T,MR,MC>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::LockedView2x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_  = ATL.Width() + ATR.Width();
    this->colAlignment_ = ATL.ColAlignment();
    this->rowAlignment_ = ATL.RowAlignment();
    this->colShift_     = ATL.ColShift();
    this->rowShift_     = ATL.RowShift();
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
DistMatrix<T,MR,MC>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error("Height and width must be non-negative");
#endif
    const elemental::Grid& g = this->Grid();
    this->height_ = height;
    this->width_ = width;
    this->localMatrix_.ResizeTo
    ( LocalLength( height, this->ColShift(), g.Width() ),
      LocalLength( width,  this->RowShift(), g.Height() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
DistMatrix<T,MR,MC>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        u = this->GetLocalEntry(iLoc,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        this->SetLocalEntry(iLoc,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateLocalEntry(iLoc,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::GetDiagonal
( DistMatrix<T,MD,STAR>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetDiagonal([MD,* ])");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Height() || d.Width() != 1) )
        throw logic_error("d is not of the correct dimensions");
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw logic_error("d must be aligned with the 'offset' diagonal");
#endif
    const elemental::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( diagLength, 1 );
    }

    if( d.InDiagonal() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
        const int diagShift = d.ColShift();

        int iStart,jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalHeight();

        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
        T* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            dLocalBuffer[k] = thisLocalBuffer[iLocal+jLocal*thisLDim]; 
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::GetDiagonal
( DistMatrix<T,STAR,MD>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetDiagonal([* ,MD])");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Width() || d.Height() != 1) )
        throw logic_error("d is not of the correct dimensions");
    if( ( d.Viewing() && d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw logic_error("d must be aligned with the 'offset' diagonal");
#endif
    const elemental::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( 1, diagLength );
    }

    if( d.InDiagonal() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
        const int diagShift = d.RowShift();

        int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalWidth();

        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
        T* dLocalBuffer = d.LocalBuffer();
        const int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            dLocalBuffer[k*dLDim] = thisLocalBuffer[iLocal+jLocal*thisLDim];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SetDiagonal
( const DistMatrix<T,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetDiagonal([MD,* ])");
    this->AssertSameGrid( d );
    if( d.Width() != 1 )
        throw logic_error("d must be a column vector");
    const int diagLength = this->DiagonalLength(offset);
    if( diagLength != d.Height() )
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw logic_error("d must be aligned with the 'offset' diagonal");
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
        const int diagShift = d.ColShift();

        int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalHeight();

        const T* dLocalBuffer = d.LockedLocalBuffer();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SetDiagonal
( const DistMatrix<T,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetDiagonal([* ,MD])");
    this->AssertSameGrid( d );
    if( d.Height() != 1 )
        throw logic_error("d must be a row vector");
    const int diagLength = this->DiagonalLength(offset);
    if( diagLength != d.Width() )
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw logic_error("d must be aligned with the 'offset' diagonal");
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
        const int diagShift = d.RowShift();

        int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalWidth();

        const T* dLocalBuffer = d.LockedLocalBuffer();
        const int dLDim = d.LocalLDim();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k*dLDim];
        }
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
DistMatrix<T,MR,MC>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::MakeTrapezoidal");
    this->AssertNotLockedView(); 
#endif
    const elemental::Grid& g = this->Grid();
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
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
            int j = rowShift + jLocal*r;
            int lastZeroRow = ( side==LEFT ? j-offset-1
                                           : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                int boundary = min( lastZeroRow+1, height );
                int numZeroRows = RawLocalLength( boundary, colShift, c );
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
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
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*r;
            int firstZeroRow = ( side==LEFT ? max(j-offset+1,0)
                                            : max(j-offset+height-width+1,0) );
            int numNonzeroRows = RawLocalLength(firstZeroRow,colShift,c);
            if( numNonzeroRows < localHeight )
            {
                T* thisCol = &thisLocalBuffer[numNonzeroRows+jLocal*thisLDim];
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
DistMatrix<T,MR,MC>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::ScaleTrapezoidal");
    this->AssertNotLockedView(); 
#endif
    const elemental::Grid& g = this->Grid();
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
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
            int j = rowShift + jLocal*r;
            int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            int boundary = min( lastRow+1, height );
            int numRows = RawLocalLength( boundary, colShift, c );
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
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
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*r;
            int firstRow = ( side==LEFT ? max(j-offset,0)
                                        : max(j-offset+height-width,0) );
            int numZeroRows = RawLocalLength( firstRow, colShift, c );
            T* thisCol = &thisLocalBuffer[numZeroRows+jLocal*thisLDim];
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
DistMatrix<T,MR,MC>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const elemental::Grid& g = this->Grid();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToZero();

    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*c;
        if( i % r == rowShift )
        {
            const int jLocal = (i-rowShift) / r;
            if( jLocal < localWidth )
                thisLocalBuffer[iLocal+jLocal*thisLDim] = 1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->SetLocalEntry(i,j,SampleUnitBall<T>());
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( A.Width() == 1 )
    {
        if( !this->Viewing() )
            ResizeTo( A.Height(), 1 );

        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int myRow = g.MCRank();
        const int myCol = g.MRRank();
        const int rankCM = g.VCRank();
        const int rankRM = g.VRRank();
        const int ownerRow = this->RowAlignment();
        const int ownerCol = A.RowAlignment();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = A.Height();
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = max(maxLocalHeight,mpi::MIN_COLL_MSG);

        const int colShiftVR = Shift(rankRM,colAlignment,p);
        const int colShiftVCOfA = Shift(rankCM,colAlignmentOfA,p);
        const int sendRankRM = (rankRM+(p+colShiftVCOfA-colShiftVR)) % p;
        const int recvRankCM = (rankCM+(p+colShiftVR-colShiftVCOfA)) % p;
        const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

        this->auxMemory_.Require( (r+c)*portionSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        if( myCol == ownerCol )
        {
            // Pack
            const int AColShift = A.ColShift();
            const T* ALocalBuffer = A.LockedLocalBuffer();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int k=0; k<c; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const int shift = RawShift(myRow+r*k,colAlignmentOfA,p);
                const int offset = (shift-AColShift) / r;
                const int thisLocalHeight = RawLocalLength(height,shift,p);

                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal] = ALocalBuffer[offset+iLocal*c];
            }
        }

        // A[VC,* ] <- A[MC,MR]
        mpi::Scatter
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerCol, g.MRComm() );

        // A[VR,* ] <- A[VC,* ]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankRM, 0,
          recvBuf, portionSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

        // A[MR,MC] <- A[VR,* ]
        mpi::Gather
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerRow, g.MCComm() );

        if( myRow == ownerRow )
        {
            // Unpack
            const int thisColShift = this->ColShift();
            T* thisLocalBuffer = this->LocalBuffer();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const int shift = RawShift(myCol+c*k,colAlignment,p);
                const int offset = (shift-thisColShift) / c;
                const int thisLocalHeight = RawLocalLength(height,shift,p);

                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    thisLocalBuffer[offset+iLocal*r] = data[iLocal];
            }
        }

        this->auxMemory_.Release();
    }
    else if( A.Height() == 1 )
    {
        if( !this->Viewing() )
            ResizeTo( 1, A.Width() );

        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int myRow = g.MCRank();
        const int myCol = g.MRRank();
        const int rankCM = g.VCRank();
        const int rankRM = g.VRRank();
        const int ownerCol = this->ColAlignment();
        const int ownerRow = A.ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int width = A.Width();
        const int maxLocalWidth = MaxLocalLength(width,p);

        const int portionSize = max(maxLocalWidth,mpi::MIN_COLL_MSG);

        const int rowShiftVC = Shift(rankCM,rowAlignment,p);
        const int rowShiftVROfA = Shift(rankRM,rowAlignmentOfA,p);
        const int sendRankCM = (rankCM+(p+rowShiftVROfA-rowShiftVC)) % p;
        const int recvRankRM = (rankRM+(p+rowShiftVC-rowShiftVROfA)) % p;
        const int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        this->auxMemory_.Require( (r+c)*portionSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        if( myRow == ownerRow )
        {
            // Pack
            const int ARowShift = A.RowShift();
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const int shift = RawShift(myCol+c*k,rowAlignmentOfA,p);
                const int offset = (shift-ARowShift) / c;
                const int thisLocalWidth = RawLocalLength(width,shift,p);

                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    data[jLocal] = ALocalBuffer[(offset+jLocal*r)*ALDim];
            }
        }

        // A[* ,VR] <- A[MC,MR]
        mpi::Scatter
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerRow, g.MCComm() );

        // A[* ,VC] <- A[* ,VR]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankCM, 0,
          recvBuf, portionSize, recvRankCM, mpi::ANY_TAG, g.VCComm() );

        // A[MR,MC] <- A[* ,VC]
        mpi::Gather
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerCol, g.MRComm() );

        if( myCol == ownerCol )
        {
            // Unpack
            const int thisRowShift = this->RowShift();
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int k=0; k<c; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const int shift = RawShift(myRow+r*k,rowAlignment,p);
                const int offset = (shift-thisRowShift) / r;
                const int thisLocalWidth = RawLocalLength(width,shift,p);

                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    thisLocalBuffer[(offset+jLocal*c)*thisLDim] = data[jLocal];
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
        if( A.Height() >= A.Width() )
        {
            auto_ptr< DistMatrix<T,VC,STAR> > A_VC_STAR
            ( new DistMatrix<T,VC,STAR>(g) );
            *A_VC_STAR = A;

            auto_ptr< DistMatrix<T,VR,STAR> > A_VR_STAR
            ( new DistMatrix<T,VR,STAR>(true,this->ColAlignment(),g) );
            *A_VR_STAR = *A_VC_STAR;
            delete A_VC_STAR.release(); // lowers memory highwater

            *this = *A_VR_STAR;
        }
        else
        {
            auto_ptr< DistMatrix<T,STAR,VR> > A_STAR_VR
            ( new DistMatrix<T,STAR,VR>(g) );
            *A_STAR_VR = A;

            auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
            ( new DistMatrix<T,STAR,VC>(true,this->RowAlignment(),g) );
            *A_STAR_VC = *A_STAR_VR;
            delete A_STAR_VR.release(); // lowers memory highwater

            *this = *A_STAR_VC;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,VC,STAR> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VC_STAR = A;

    auto_ptr< DistMatrix<T,VR,STAR> > A_VR_STAR
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
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,STAR,VR> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR>(g) );
    *A_STAR_VR = A;

    auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
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
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error("[MR,MC] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error("[MR,MC] = [* ,MD] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MR,MC]");
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
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            this->rowShift_ = A.RowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() &&
        this->RowAlignment() == A.RowAlignment() )
    {
        this->localMatrix_ = A.LockedLocalMatrix();
    }
    else
    {
        const elemental::Grid& g = this->Grid();
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,MC] <- [MR,MC]" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int row = g.MCRank();
        const int col = g.MRRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;
        const int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;
        const int sendRank = sendRow + sendCol*r;
        const int recvRank = recvRow + recvCol*r;

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = localHeightOfA * localWidthOfA;
        const int recvSize = localHeight * localWidth;

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
            T* sendBufferCol = &sendBuffer[jLocal*localHeightOfA];
            memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
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
            const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
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
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MR,* ]");
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
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,MC] <- [MR,* ]" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int rowShift = this->RowShift();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (col+c+colAlignment-colAlignmentOfA) % c;
        const int recvRank = (col+c+colAlignmentOfA-colAlignment) % c;

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int sendSize = localHeightOfA * localWidth;
        const int recvSize = localHeight * localWidth;

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
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[(rowShift+jLocal*r)*ALDim];
            T* sendBufferCol = &sendBuffer[jLocal*localHeightOfA];
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
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
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
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,MC]");
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
                Shift( g.MCRank(), this->RowAlignment(), g.Height() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const int c = g.Width();
        const int colShift = this->ColShift();

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();

        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                thisLocalBuffer[iLocal+jLocal*thisLDim] = 
                    ALocalBuffer[(colShift+iLocal*c)+jLocal*ALDim];
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,MC] <- [* ,MC]" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int row = g.MCRank(); 

        const int colShift = this->ColShift();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = localHeight * localWidthOfA;
        const int recvSize = localHeight * localWidth;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                sendBuffer[iLocal+jLocal*localHeight] = 
                    ALocalBuffer[(colShift+iLocal*c)+jLocal];

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
            const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
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
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,VR,STAR> A_VR_STAR(g);

    A_VR_STAR = A;
    *this     = A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,VC]");
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
        const int p = r * c;
        const int row = g.MCRank();

        const int rowShift = this->RowShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int maxHeight = MaxLocalLength(height,c);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

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

            const int thisColShift = RawShift(k,colAlignment,c);
            const int thisLocalHeight = RawLocalLength(height,thisColShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+jLocal*ALDim];
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

            const int thisRank = row+k*r;
            const int thisRowShift = RawShift(thisRank,rowAlignmentOfA,p);
            const int thisRowOffset = (thisRowShift-rowShift) / r;
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = 
                    &thisLocalBuffer[(thisRowOffset+jLocal*c)*thisLDim];
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,MC] <- [* ,VC]" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int row = g.MCRank();

        const int rowShift = this->RowShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+rowAlignment-(rowAlignmentOfA%r)) % r;
        const int recvRow = (row+r+(rowAlignmentOfA%r)-rowAlignment) % r;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int maxHeight = MaxLocalLength(height,c);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

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

            const int thisColShift = RawShift(k,colAlignment,c);
            const int thisLocalHeight = RawLocalLength(height,thisColShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+jLocal*ALDim];
        }

        // SendRecv to align A[* ,VC] via a trade in the column
        mpi::SendRecv
        ( secondBuffer, c*portionSize, sendRow, 0,
          firstBuffer,  c*portionSize, recvRow, mpi::ANY_TAG, g.MCComm() );

        // AllToAll to gather all of the aligned [* ,VC] into secondBuffer
        mpi::AllToAll
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

            const int thisRank = recvRow+k*r;
            const int thisRowShift = RawShift(thisRank,rowAlignmentOfA,p);
            const int thisRowOffset = (thisRowShift-rowShift) / r;
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = 
                    &thisLocalBuffer[(thisRowOffset+jLocal*c)*thisLDim];
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
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [VR,* ]");
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

        const int colShift = this->ColShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

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

            const int thisRowShift = RawShift(k,rowAlignment,r);
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
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

            const int thisRank = col+k*c;
            const int thisColShift = RawShift(thisRank,colAlignmentOfA,p);
            const int thisColOffset = (thisColShift-colShift) / c;
            const int thisLocalHeight = RawLocalLength(height,thisColShift,p);
            
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    thisLocalBuffer[(thisColOffset+iLocal*r)+jLocal*thisLDim] =
                        data[iLocal+jLocal*thisLocalHeight];
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,MC] <- [* ,VC]" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();

        const int colShift = this->ColShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendCol = (col+c+colAlignment-(colAlignmentOfA%c)) % c;
        const int recvCol = (col+c+(colAlignmentOfA%c)-colAlignment) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

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

            const int thisRowShift = RawShift(k,rowAlignment,r);
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
            }
        }

        // SendRecv to align A[VR,* ] via a trade in the row
        mpi::SendRecv
        ( secondBuffer, r*portionSize, sendCol, 0,
          firstBuffer,  r*portionSize, recvCol, mpi::ANY_TAG, g.MRComm() );

        // AllToAll to gather all of the aligned [VR,* ] data into secondBuffer
        mpi::AllToAll
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

            const int thisRank = recvCol+k*c;
            const int thisColShift = RawShift(thisRank,colAlignmentOfA,p);
            const int thisColOffset = (thisColShift-colShift) / c;
            const int thisLocalHeight = RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    thisLocalBuffer[(thisColOffset+iLocal*r)+jLocal*thisLDim] = 
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
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,VR]");
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
inline const DistMatrix<T,MR,MC>&
DistMatrix<T,MR,MC>::operator=( const DistMatrix<T,STAR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,* ]");
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
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();

    const T* ALocalBuffer = A.LockedLocalBuffer();
    const int ALDim = A.LocalLDim();
    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            thisLocalBuffer[iLocal+jLocal*thisLDim] = 
                ALocalBuffer[(colShift+iLocal*c)+(rowShift+jLocal*r)*ALDim];
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SumScatterFrom( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterFrom([MR,* ])");
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
        if( this->Width() == 1 )
        {
            const int rowAlignment = this->RowAlignment();
            const int myRow = g.MCRank();

            const int localHeight = this->LocalHeight();

            const int portionSize = max(localHeight,mpi::MIN_COLL_MSG);

            this->auxMemory_.Require( 2*portionSize );

            T* buffer = this->auxMemory_.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack 
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeight*sizeof(T) );

            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuffer, recvBuffer, portionSize,
              mpi::SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
            {
                T* thisCol = this->LocalBuffer(0,0);
                memcpy( thisCol, recvBuffer, localHeight*sizeof(T) );
            }

            this->auxMemory_.Release();
        }
        else
        {
            const int r = g.Height();
            const int rowAlignment = this->RowAlignment();

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int recvSize = 
                max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
            const int sendSize = r * recvSize;

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

                const int thisRowShift = RawShift( k, rowAlignment, r );
                const int thisLocalWidth = 
                      RawLocalLength( width, thisRowShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
                    T* dataCol = &data[jLocal*localHeight];
                    memcpy( dataCol, ACol, localHeight*sizeof(T) );
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
            for( int j=0; j<localWidth; ++j )
            {
                const T* recvBufferCol = &recvBuffer[j*localHeight];
                T* thisCol = &thisLocalBuffer[j*thisLDim];
                memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned SumScatterFrom [MR,MC] <- [MR,* ]" << endl;
#endif
        if( this->Width() == 1 )
        {
            const int c = g.Width();
            const int rowAlignment = this->RowAlignment();
            const int myRow = g.MCRank();
            const int myCol = g.MRRank();

            const int height = this->Height();
            const int localHeight = this->LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,c);

            const int portionSize = max(maxLocalHeight,mpi::MIN_COLL_MSG);

            const int colAlignment = this->ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (myCol+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (myCol+c+colAlignmentOfA-colAlignment) % c;

            this->auxMemory_.Require( 2*portionSize );

            T* buffer = this->auxMemory_.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeightOfA*sizeof(T) );

            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuffer, recvBuffer, portionSize,
              mpi::SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
            {
                // Perform the realignment
                mpi::SendRecv
                ( recvBuffer, portionSize, sendCol, 0,
                  sendBuffer, portionSize, recvCol, 0, g.MRComm() );

                T* thisCol = this->LocalBuffer(0,0);
                memcpy( thisCol, sendBuffer, localHeight*sizeof(T) );
            }

            this->auxMemory_.Release();
        }
        else
        {
            const int r = g.Height();
            const int c = g.Width();
            const int col = g.MRRank();

            const int colAlignment = this->ColAlignment();
            const int rowAlignment = this->RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int recvSize_RS = 
                max(localHeightOfA*maxLocalWidth,mpi::MIN_COLL_MSG);
            const int sendSize_RS = r * recvSize_RS;
            const int recvSize_SR = localHeight * localWidth;

            this->auxMemory_.Require
            ( recvSize_RS + max(sendSize_RS,recvSize_SR) );

            T* buffer = this->auxMemory_.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[recvSize_RS];

            // Pack
            vector<int> recvSizes(r);
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &secondBuffer[k*recvSize_RS];
                recvSizes[k] = recvSize_RS;

                const int thisRowShift = RawShift( k, rowAlignment, r );
                const int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
                    T* dataCol = &data[jLocal*localHeightOfA];
                    memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                }
            }

            // Reduce-scatter over each process col
            mpi::ReduceScatter
            ( secondBuffer, firstBuffer, &recvSizes[0], mpi::SUM, g.MCComm() );

            // Trade reduced data with the appropriate process col
            mpi::SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendCol, 0,
              secondBuffer, localHeight*localWidth,    recvCol, mpi::ANY_TAG,
              g.MRComm() );

            // Unpack the received data
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* secondBufferCol = &secondBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                memcpy( thisCol, secondBufferCol, localHeight*sizeof(T) );
            }
            this->auxMemory_.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SumScatterFrom( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterFrom([* ,MC])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Height() == 1 && g.VCRank() == 0 )
    {
        cerr <<    
          "The vector version of [MR,MC].SumScatterFrom([* ,MC]) is not yet"
          " written, but it only requires a modification of the vector "
          "version of [MR,MC].SumScatterFrom([MR,* ])" << endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.VCRank() == 0 )
    {
        cerr << 
          "[MR,MC]::SumScatterFrom([* ,MC]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MC] matrix instead." << endl;
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
        const int colAlignment = this->ColAlignment();

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);

        const int recvSize = max(maxLocalHeight*localWidth,mpi::MIN_COLL_MSG);
        const int sendSize = c * recvSize;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack 
        vector<int> recvSizes(c);
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisColShift = RawShift( k, colAlignment, c );
            const int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+jLocal*ALDim];
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
            const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned SumScatterFrom [MR,MC] <- [* ,MC]" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int row = g.MCRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);
        
        const int recvSize_RS = 
            max(maxLocalHeight*localWidthOfA,mpi::MIN_COLL_MSG);
        const int sendSize_RS = c* recvSize_RS;
        const int recvSize_SR = localHeight * localWidth;

        this->auxMemory_.Require( recvSize_RS + max(sendSize_RS,recvSize_SR) );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[recvSize_RS];

        // Pack 
        vector<int> recvSizes(c);
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &secondBuffer[k*recvSize_RS];
            recvSizes[k] = recvSize_RS;

            const int thisColShift = RawShift( k, colAlignment, c );
            const int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+jLocal*ALDim];
        }

        // Reduce-scatter over each process row
        mpi::ReduceScatter
        ( secondBuffer, firstBuffer, &recvSizes[0], mpi::SUM, g.MRComm() );

        // Trade reduced data with the appropriate process row
        mpi::SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendRow, 0,
          secondBuffer, localHeight*localWidth,    recvRow, mpi::ANY_TAG, 
          g.MCComm() );

        // Unpack the received data
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* secondBufferCol = &secondBuffer[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            memcpy( thisCol, secondBufferCol, localHeight*sizeof(T) );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SumScatterFrom
( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterFrom([* ,* ])");
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
    const int colAlignment = this->ColAlignment();
    const int rowAlignment = this->RowAlignment();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int maxLocalHeight = MaxLocalLength(height,c);
    const int maxLocalWidth = MaxLocalLength(width,r);

    const int recvSize = max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
    const int sendSize = r * c * recvSize;

    this->auxMemory_.Require( sendSize + recvSize );

    T* buffer = this->auxMemory_.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack 
    vector<int> recvSizes(r*c);
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int l=0; l<r; ++l )
    {
        const int thisRowShift = RawShift( l, rowAlignment, r );
        const int thisLocalWidth = RawLocalLength( width, thisRowShift, r );

        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[(k+l*c)*recvSize];
            recvSizes[k+l*c] = recvSize;

            const int thisColShift = RawShift( k, colAlignment, c );
            const int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+
                                     (thisRowShift+jLocal*r)*ALDim];
        }
    }

    // Reduce-scatter over each process col
    mpi::ReduceScatter
    ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.VRComm() );

    // Unpack our received data
    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SumScatterUpdate
( T alpha, const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterUpdate([MR,* ])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A ); 
#endif
    const elemental::Grid& g = this->Grid();
    if( this->ColAlignment() == A.ColAlignment() )
    {
        if( this->Width() == 1 )
        {
            const int rowAlignment = this->RowAlignment();
            const int myRow = g.MCRank();

            const int localHeight = this->LocalHeight();

            const int portionSize = max(localHeight,mpi::MIN_COLL_MSG);

            this->auxMemory_.Require( 2*portionSize );

            T* buffer = this->auxMemory_.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeight*sizeof(T) );

            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuffer, recvBuffer, portionSize,
              mpi::SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
            {
                T* thisCol = this->LocalBuffer(0,0);
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
                #pragma omp parallel for
#endif
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisCol[iLocal] += alpha*recvBuffer[iLocal];
            }
            this->auxMemory_.Release();
        }
        else
        {
            const int r = g.Height();
            const int rowAlignment = this->RowAlignment();

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int portionSize = 
                max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

            this->auxMemory_.Require( (r+1)*portionSize );

            T* buffer = this->auxMemory_.Buffer();
            T* recvBuffer = &buffer[0];
            T* sendBuffer = &buffer[portionSize];

            // Pack 
            vector<int> recvSizes(r);
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &sendBuffer[k*portionSize];
                recvSizes[k] = portionSize;

                const int thisRowShift = RawShift( k, rowAlignment, r );
                const int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
                    T* dataCol = &data[jLocal*localHeight];
                    memcpy( dataCol, ACol, localHeight*sizeof(T) );
                }
            }

            // Reduce-scatter over each process column
            mpi::ReduceScatter
            ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.MCComm() );

            // Update with our received data
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisCol[iLocal] += alpha*recvBufferCol[iLocal];
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned SumScatterUpdate [MR,MC] <- [MR,* ]" << endl;
#endif
        if( this->Width() == 1 )
        {
            const int c = g.Width();
            const int rowAlignment = this->RowAlignment();
            const int myRow = g.MCRank();
            const int myCol = g.MRRank();

            const int height = this->Height();
            const int localHeight = this->LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,c);

            const int portionSize = max(maxLocalHeight,mpi::MIN_COLL_MSG);

            const int colAlignment = this->ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (myCol+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (myCol+c+colAlignmentOfA-colAlignment) % c;

            this->auxMemory_.Require( 2*portionSize );

            T* buffer = this->auxMemory_.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeightOfA*sizeof(T) );

            // Reduce to rowAlignment
            mpi::Reduce
            ( sendBuffer, recvBuffer, portionSize,
              mpi::SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
            {
                // Perform the realignment
                mpi::SendRecv
                ( recvBuffer, portionSize, sendCol, 0,
                  sendBuffer, portionSize, recvCol, 0, g.MRComm() );

                T* thisCol = this->LocalBuffer(0,0);
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
                #pragma omp parallel for
#endif
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisCol[iLocal] += alpha*sendBuffer[iLocal];
            }
            this->auxMemory_.Release();
        }
        else
        {
            const int r = g.Height();
            const int c = g.Width();
            const int col = g.MRRank();

            const int colAlignment = this->ColAlignment();
            const int rowAlignment = this->RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int recvSize_RS = 
                max(localHeightOfA*maxLocalWidth,mpi::MIN_COLL_MSG);
            const int sendSize_RS = r * recvSize_RS;
            const int recvSize_SR = localHeight * localWidth;

            this->auxMemory_.Require
            ( recvSize_RS + max(sendSize_RS,recvSize_SR) );

            T* buffer = this->auxMemory_.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[recvSize_RS];

            // Pack 
            vector<int> recvSizes(r);
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &secondBuffer[k*recvSize_RS];
                recvSizes[k] = recvSize_RS;

                const int thisRowShift = RawShift( k, rowAlignment, r );
                const int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
                    T* dataCol = &data[jLocal*localHeightOfA];
                    memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                }
            }

            // Reduce-scatter over each process col
            mpi::ReduceScatter
            ( secondBuffer, firstBuffer, &recvSizes[0], mpi::SUM, g.MCComm() );

            // Trade reduced data with the appropriate process col
            mpi::SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendCol, 0, 
              secondBuffer, localHeight*localWidth,    recvCol, mpi::ANY_TAG,
              g.MRComm() );

            // Update with our received data
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* secondBufferCol = &secondBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisCol[iLocal] += alpha*secondBufferCol[iLocal];
            }
            this->auxMemory_.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterUpdate([* ,MC])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Height() == 1 && g.VCRank() == 0 )
    {
        cerr <<    
          "The vector version of [MR,MC].SumScatterUpdate([* ,MC]) is not "
          "yet written, but it only requires a modification of the vector "
          "version of [MR,MC].SumScatterUpdate([MR,* ])." << endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.VCRank() == 0 )
    {
        cerr <<
          "[MR,MC]::SumScatterUpdate([* ,MC]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MC] matrix instead." << endl;
    }
#endif
    if( this->RowAlignment() == A.RowAlignment() )
    {
        const int c = g.Width();
        const int colAlignment = this->ColAlignment();

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);

        const int recvSize = max(maxLocalHeight*localWidth,mpi::MIN_COLL_MSG);
        const int sendSize = c * recvSize;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        vector<int> recvSizes(c);
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisColShift = RawShift( k, colAlignment, c );
            const int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+jLocal*ALDim];
        }

        // Reduce-scatter over each process row
        mpi::ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.MRComm() );

        // Update with our received data
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                thisCol[iLocal] += alpha*recvBufferCol[iLocal];
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned SumScatterUpdate [MR,MC] <- [* ,MC]" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int row = g.MCRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);

        const int recvSize_RS = 
            max(maxLocalHeight*localWidthOfA,mpi::MIN_COLL_MSG);
        const int sendSize_RS = c * recvSize_RS;
        const int recvSize_SR = localHeight * localWidth;

        this->auxMemory_.Require( recvSize_RS + max(sendSize_RS,recvSize_SR) );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[recvSize_RS];

        // Pack 
        vector<int> recvSizes(c);
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &secondBuffer[k*recvSize_RS];
            recvSizes[k] = recvSize_RS;

            const int thisColShift = RawShift( k, colAlignment, c );
            const int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+jLocal*ALDim];
        }

        // Reduce-scatter over each process row
        mpi::ReduceScatter
        ( secondBuffer, firstBuffer, &recvSizes[0], mpi::SUM, g.MRComm() );

        // Trade reduced data with the appropriate process row
        mpi::SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendRow, 0,
          secondBuffer, localHeight*localWidth,    recvRow, mpi::ANY_TAG,
          g.MRComm() );

        // Update with our received data
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* secondBufferCol = &secondBuffer[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                thisCol[iLocal] += alpha*secondBufferCol[iLocal];
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterUpdate([* ,* ])");
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
    const int colAlignment = this->ColAlignment();
    const int rowAlignment = this->RowAlignment();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int maxLocalHeight = MaxLocalLength(height,c);
    const int maxLocalWidth = MaxLocalLength(width,r);

    const int recvSize = max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
    const int sendSize = r * c * recvSize;

    this->auxMemory_.Require( sendSize + recvSize );

    T* buffer = this->auxMemory_.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack 
    vector<int> recvSizes(r*c);
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int l=0; l<r; ++l )
    {
        const int thisRowShift = RawShift( l, rowAlignment, r );
        const int thisLocalWidth = RawLocalLength( width, thisRowShift, r );

        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[(k+l*c)*recvSize];
            recvSizes[k+l*c] = recvSize;

            const int thisColShift = RawShift( k, colAlignment, c );
            const int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c) + 
                                     (thisRowShift+jLocal*r)*ALDim];
        }
    }

    // Reduce-scatter over each process col
    mpi::ReduceScatter
    ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.VRComm() );

    // Unpack our received data
    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        blas::Axpy( localHeight, alpha, recvBufferCol, 1, thisCol, 1 );
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
