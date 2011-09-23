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
DistMatrix<T,MC,MR>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();

    const int r = g.Height();
    const int c = g.Width();

    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;

    const int height = this->Height();
    const int width  = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth  = this->LocalWidth();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( g.InGrid() )
    {
        // Fill the send buffer: zero it then place our entries into their 
        // appropriate locations
        std::vector<T> sendBuf(height*width,0);
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                sendBuf[(colShift+iLocal*r)+(rowShift+jLocal*c)*height] = 
                    thisLocalBuffer[iLocal+jLocal*thisLDim];

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
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::Align( int colAlignment, int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::Align");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Height() )
        throw std::runtime_error("Invalid column alignment for [MC,MR]");
    if( rowAlignment < 0 || rowAlignment >= g.Width() )
        throw std::runtime_error("Invalid row alignment for [MC,MR]");
#endif
    this->_colAlignment = colAlignment;
    this->_rowAlignment = rowAlignment;
    this->_constrainedColAlignment = true;
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    if( g.InGrid() )
    {
        this->_localMatrix.ResizeTo( 0, 0 );
        this->_colShift = Shift( g.MCRank(), colAlignment, g.Height() );
        this->_rowShift = Shift( g.MRRank(), rowAlignment, g.Width() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::AlignCols( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Height() )
        throw std::runtime_error( "Invalid column alignment for [MC,MR]" );
#endif
    this->_colAlignment = colAlignment;
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    if( g.InGrid() )
    {
        this->_localMatrix.ResizeTo( 0, 0 );
        this->_colShift = Shift( g.MCRank(), colAlignment, g.Height() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::AlignRows( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Width() )
        throw std::runtime_error( "Invalid row alignment for [MC,MR]" );
#endif
    this->_rowAlignment = rowAlignment;
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    if( g.InGrid() )
    {
        this->_localMatrix.ResizeTo( 0, 0 );
        this->_rowShift = Shift( g.MRRank(), rowAlignment, g.Width() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::View( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->_grid = A._grid;
    this->_height = A.Height();
    this->_width = A.Width();
    this->_colAlignment = A.ColAlignment();
    this->_rowAlignment = A.RowAlignment();
    this->_viewing = true;
    this->_lockedView = false;
    if( this->Grid().InGrid() )
    {
        this->_colShift = A.ColShift();
        this->_rowShift = A.RowShift();
        this->_localMatrix.View( A.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::View
( int height, int width, int colAlignment, int rowAlignment, 
  T* buffer, int ldim, const elemental::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->_grid = &g;
    this->_height = height;
    this->_width = width;
    this->_colAlignment = colAlignment;
    this->_rowAlignment = rowAlignment;
    this->_viewing = true;
    this->_lockedView = false;
    if( this->_grid->InGrid() )
    {
        this->_colShift = Shift(g.MCRank(),colAlignment,g.Height());
        this->_rowShift = Shift(g.MRRank(),rowAlignment,g.Width());
        const int localHeight = LocalLength(height,this->_colShift,g.Height());
        const int localWidth = LocalLength(width,this->_rowShift,g.Width());
        this->_localMatrix.View( localHeight, localWidth, buffer, ldim );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::LockedView( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->_grid = A._grid;
    this->_height = A.Height();
    this->_width  = A.Width();
    this->_colAlignment = A.ColAlignment();
    this->_rowAlignment = A.RowAlignment();
    this->_viewing = true;
    this->_lockedView = true;
    if( this->Grid().InGrid() )
    {
        this->_colShift = A.ColShift();
        this->_rowShift = A.RowShift();
        this->_localMatrix.LockedView( A.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::LockedView
( int height, int width, int colAlignment, int rowAlignment, 
  const T* buffer, int ldim, const elemental::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->_grid = &g;
    this->_height = height;
    this->_width = width;
    this->_colAlignment = colAlignment;
    this->_rowAlignment = rowAlignment;
    this->_viewing = true;
    this->_lockedView = true;
    if( this->_grid->InGrid() )
    {
        this->_colShift = Shift(g.MCRank(),colAlignment,g.Height());
        this->_rowShift = Shift(g.MRRank(),rowAlignment,g.Width());
        const int localHeight = LocalLength(height,this->_colShift,g.Height());
        const int localWidth = LocalLength(width,this->_rowShift,g.Width());
        this->_localMatrix.LockedView( localHeight, localWidth, buffer, ldim );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::View
( DistMatrix<T,MC,MR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_grid = A._grid;
    this->_height = height;
    this->_width  = width;

    const elemental::Grid& g = this->Grid();
    const int r   = g.Height();
    const int c   = g.Width();
    const int row = g.MCRank();
    const int col = g.MRRank();

    this->_colAlignment = (A.ColAlignment()+i) % r;
    this->_rowAlignment = (A.RowAlignment()+j) % c;
    this->_viewing = true;
    this->_lockedView = false;
  
    if( g.InGrid() )
    {
        this->_colShift = Shift( row, this->ColAlignment(), r );
        this->_rowShift = Shift( col, this->RowAlignment(), c );

        const int localHeightBehind = LocalLength(i,A.ColShift(),r);
        const int localWidthBehind  = LocalLength(j,A.RowShift(),c);

        const int localHeight = LocalLength( height, this->ColShift(), r );
        const int localWidth  = LocalLength( width,  this->RowShift(), c );

        this->_localMatrix.View
        ( A.LocalMatrix(), localHeightBehind, localWidthBehind,
                           localHeight,       localWidth );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::LockedView
( const DistMatrix<T,MC,MR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_grid = A._grid;
    this->_height = height;
    this->_width  = width;

    const elemental::Grid& g = this->Grid();
    const int r   = g.Height();
    const int c   = g.Width();
    const int row = g.MCRank();
    const int col = g.MRRank();

    this->_colAlignment = (A.ColAlignment()+i) % r;
    this->_rowAlignment = (A.RowAlignment()+j) % c;
    this->_viewing = true;
    this->_lockedView = true;
  
    if( g.InGrid() )
    {
        this->_colShift = Shift( row, this->ColAlignment(), r );
        this->_rowShift = Shift( col, this->RowAlignment(), c );

        const int localHeightBehind = LocalLength(i,A.ColShift(),r);
        const int localWidthBehind  = LocalLength(j,A.RowShift(),c);

        const int localHeight = LocalLength( height, this->ColShift(), r );
        const int localWidth  = LocalLength( width,  this->RowShift(), c );

        this->_localMatrix.LockedView
        ( A.LockedLocalMatrix(), localHeightBehind, localWidthBehind,
                                 localHeight,       localWidth );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::View1x2
( DistMatrix<T,MC,MR>& AL, DistMatrix<T,MC,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View1x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->_grid = AL._grid;
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_colAlignment = AL.ColAlignment();
    this->_rowAlignment = AL.RowAlignment();
    this->_viewing = true;
    this->_lockedView = false;
    if( this->Grid().InGrid() )
    {
        this->_colShift = AL.ColShift();
        this->_rowShift = AL.RowShift();
        this->_localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::LockedView1x2
( const DistMatrix<T,MC,MR>& AL, const DistMatrix<T,MC,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView1x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->_grid = AL._grid;
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_colAlignment = AL.ColAlignment();
    this->_rowAlignment = AL.RowAlignment();
    this->_viewing = true;
    this->_lockedView = true;
    if( this->Grid().InGrid() )
    {
        this->_colShift = AL.ColShift();
        this->_rowShift = AL.RowShift();
        this->_localMatrix.LockedView1x2
        ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::View2x1
( DistMatrix<T,MC,MR>& AT,
  DistMatrix<T,MC,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View2x1");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->_grid = AT._grid;
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_colAlignment = AT.ColAlignment();
    this->_rowAlignment = AT.RowAlignment();
    this->_viewing = true;
    this->_lockedView = false;
    if( this->Grid().InGrid() )
    {
        this->_colShift = AT.ColShift();
        this->_rowShift = AT.RowShift();
        this->_localMatrix.View2x1
        ( AT.LocalMatrix(), 
          AB.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::LockedView2x1
( const DistMatrix<T,MC,MR>& AT,
  const DistMatrix<T,MC,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView2x1");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->_grid = AT._grid;
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_colAlignment = AT.ColAlignment();
    this->_rowAlignment = AT.RowAlignment();
    this->_viewing = true;
    this->_lockedView = true;
    if( this->Grid().InGrid() )
    {
        this->_colShift = AT.ColShift();
        this->_rowShift = AT.RowShift();
        this->_localMatrix.LockedView2x1
        ( AT.LockedLocalMatrix(), 
          AB.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::View2x2
( DistMatrix<T,MC,MR>& ATL, DistMatrix<T,MC,MR>& ATR,
  DistMatrix<T,MC,MR>& ABL, DistMatrix<T,MC,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View2x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->_grid = ATL._grid;
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_colAlignment = ATL.ColAlignment();
    this->_rowAlignment = ATL.RowAlignment();
    this->_viewing = true;
    this->_lockedView = false;
    if( this->Grid().InGrid() )
    {
        this->_colShift = ATL.ColShift();
        this->_rowShift = ATL.RowShift();
        this->_localMatrix.View2x2
        ( ATL.LocalMatrix(), ATR.LocalMatrix(),
          ABL.LocalMatrix(), ABR.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::LockedView2x2
( const DistMatrix<T,MC,MR>& ATL, const DistMatrix<T,MC,MR>& ATR,
  const DistMatrix<T,MC,MR>& ABL, const DistMatrix<T,MC,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView2x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->_grid = ATL._grid;
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_colAlignment = ATL.ColAlignment();
    this->_rowAlignment = ATL.RowAlignment();
    this->_viewing = true;
    this->_lockedView = true;
    if( this->Grid().InGrid() )
    {
        this->_colShift = ATL.ColShift();
        this->_rowShift = ATL.RowShift();
        this->_localMatrix.LockedView2x2
        ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
          ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::ResizeTo");
    this->AssertNotLockedView(); // this should be relaxed...
#endif
    this->_height = height;
    this->_width = width;
    if( this->Grid().InGrid() )
    {
        this->_localMatrix.ResizeTo
        ( LocalLength(height,this->ColShift(),this->Grid().Height()),
          LocalLength(width, this->RowShift(),this->Grid().Width()) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
DistMatrix<T,MC,MR>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const int iLocal = (i-this->ColShift()) / g.Height();
        const int jLocal = (j-this->RowShift()) / g.Width();
        u = this->GetLocalEntry(iLocal,jLocal);
    }
    mpi::Broadcast
    ( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLocal = (i-this->ColShift()) / g.Height();
        const int jLocal = (j-this->RowShift()) / g.Width();
        this->SetLocalEntry(iLocal,jLocal,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLocal = (i-this->ColShift()) / g.Height();
        const int jLocal = (j-this->RowShift()) / g.Width();
        this->UpdateLocalEntry(iLocal,jLocal,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::GetDiagonal
( DistMatrix<T,MD,STAR>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && diagLength != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
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

        const int iLocalStart = (iStart-colShift) / r;
        const int jLocalStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
        T* dLocalBuffer = d.LocalBuffer();
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k*(lcm/r);
            const int jLocal = jLocalStart+k*(lcm/c);
            dLocalBuffer[k] = thisLocalBuffer[iLocal+jLocal*thisLDim];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::GetDiagonal
( DistMatrix<T,STAR,MD>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && diagLength != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
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

        const int iLocalStart = (iStart-colShift) / r;
        const int jLocalStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
        T* dLocalBuffer = d.LocalBuffer();
        const int dLDim = d.LocalLDim();
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k*(lcm/r);
            const int jLocal = jLocalStart+k*(lcm/c);
            dLocalBuffer[k*dLDim] = thisLocalBuffer[iLocal+jLocal*thisLDim];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::SetDiagonal
( const DistMatrix<T,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetDiagonal");
    this->AssertSameGrid( d );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const int diagLength = this->DiagonalLength(offset);
    if( diagLength != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
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

        const int iLocalStart = (iStart-colShift) / r;
        const int jLocalStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
        const T* dLocalBuffer = d.LockedLocalBuffer();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k*(lcm/r);
            const int jLocal = jLocalStart+k*(lcm/c);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::SetDiagonal
( const DistMatrix<T,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetDiagonal");
    this->AssertSameGrid( d );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const int diagLength = this->DiagonalLength(offset);
    if( diagLength != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
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

        const int iLocalStart = (iStart-colShift) / r;
        const int jLocalStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
        const T* dLocalBuffer = d.LockedLocalBuffer();
        T* thisLocalBuffer = this->LocalBuffer();
        const int dLDim = d.LocalLDim();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k*(lcm/r);
            const int jLocal = jLocalStart+k*(lcm/c);
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
DistMatrix<T,MC,MR>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = this->Grid().Height();
    const int c = this->Grid().Width();
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
            int j = rowShift + jLocal*c;
            int lastZeroRow = ( side==LEFT ? j-offset-1
                                           : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                int boundary = std::min( lastZeroRow+1, height );
                int numZeroRows = RawLocalLength( boundary, colShift, r );
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
            int j = rowShift + jLocal*c;
            int firstZeroRow = 
                ( side==LEFT ? std::max(j-offset+1,0)
                             : std::max(j-offset+height-width+1,0) );
            int numNonzeroRows = RawLocalLength(firstZeroRow,colShift,r);
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
DistMatrix<T,MC,MR>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = this->Grid().Height();
    const int c = this->Grid().Width();
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
            int j = rowShift + jLocal*c;
            int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            int boundary = std::min( lastRow+1, height );
            int numRows = RawLocalLength( boundary, colShift, r );
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
            int j = rowShift + jLocal*c;
            int firstRow = ( side==LEFT ? std::max(j-offset,0) 
                                        : std::max(j-offset+height-width,0) );
            int numZeroRows = RawLocalLength( firstRow, colShift, r );
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
DistMatrix<T,MC,MR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = this->Grid().Height();
    const int c = this->Grid().Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->_localMatrix.SetToZero();
    T* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*r;                
        if( i % c == rowShift )
        {
            const int jLocal = (i-rowShift) / c;
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
DistMatrix<T,MC,MR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            this->SetLocalEntry(iLocal,jLocal,SampleUnitBall<T>());
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::AdjointFrom( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AdjointFrom");
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
            this->_colAlignment = A.RowAlignment();
            if( g.InGrid() )
            {
                this->_colShift = 
                    Shift( g.MCRank(), this->ColAlignment(), g.Height() );
            }
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( g.InGrid() ) 
    { 
        if( this->ColAlignment() == A.RowAlignment() )
        {
            const int c = g.Width();
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
                        Conj( ALocalBuffer[(rowShift+jLocal*c)+iLocal*ALDim] );
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR]::AdjointFrom." << std::endl;
#endif
            const int r = g.Height();
            const int c = g.Width();
            const int rank = g.MCRank();
            const int rowShift = this->RowShift();
            const int colAlignment = this->ColAlignment();
            const int rowAlignmentOfA = A.RowAlignment();

            const int sendRank = (rank+r+colAlignment-rowAlignmentOfA) % r;
            const int recvRank = (rank+r+rowAlignmentOfA-colAlignment) % r;

            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localWidthOfA = A.LocalWidth();

            const int sendSize = localWidthOfA * localWidth;
            const int recvSize = localHeight * localWidth;

            this->_auxMemory.Require( sendSize + recvSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];

            // Pack
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<localWidthOfA; ++iLocal )
                    sendBuffer[iLocal+jLocal*localWidth] = 
                        Conj( ALocalBuffer[(rowShift+jLocal*c)+iLocal*ALDim] );

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
                const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->_auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::AdjointFrom( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AdjointFrom");
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
            this->_rowAlignment = A.ColAlignment();
            if( g.InGrid() )
                this->_rowShift = 
                    Shift( g.MRRank(), this->RowAlignment(), g.Height() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }
    if( this->_rowAlignment != A.ColAlignment() )
        throw std::logic_error("Unaligned AdjointFrom");

    if( g.InGrid() ) 
    { 
        const int r = g.Height();
        const int colShift = this->ColShift();

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
                    Conj(ALocalBuffer[jLocal+(colShift+iLocal*r)*ALDim]);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::TransposeFrom( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::TransposeFrom");
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
            this->_colAlignment = A.RowAlignment();
            if( g.InGrid() )
            {
                this->_colShift = 
                    Shift( g.MCRank(), this->ColAlignment(), g.Height() );
            }
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( g.InGrid() ) 
    { 
        if( this->ColAlignment() == A.RowAlignment() )
        {
            const int c = g.Width();
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
                        ALocalBuffer[(rowShift+jLocal*c)+iLocal*ALDim];
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR]::TransposeFrom." << std::endl;
#endif
            const int r = g.Height();
            const int c = g.Width();
            const int rank = g.MCRank();
            const int rowShift = this->RowShift();
            const int colAlignment = this->ColAlignment();
            const int rowAlignmentOfA = A.RowAlignment();

            const int sendRank = (rank+r+colAlignment-rowAlignmentOfA) % r;
            const int recvRank = (rank+r+rowAlignmentOfA-colAlignment) % r;

            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localWidthOfA = A.LocalWidth();

            const int sendSize = localWidthOfA * localWidth;
            const int recvSize = localHeight * localWidth;

            this->_auxMemory.Require( sendSize + recvSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];

            // Pack
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<localWidthOfA; ++iLocal )
                    sendBuffer[iLocal+jLocal*localWidth] = 
                        ALocalBuffer[(rowShift+jLocal*c)+iLocal*ALDim];

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
                const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->_auxMemory.Release();
        }
    } 
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::TransposeFrom( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::TransposeFrom");
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
            this->_rowAlignment = A.ColAlignment();
            if( g.InGrid() )
            {
                this->_rowShift = 
                    Shift( g.MRRank(), this->RowAlignment(), g.Height() );
            }
        }
        this->ResizeTo( A.Width(), A.Height() );
    }
    if( this->_rowAlignment != A.ColAlignment() )
        throw std::logic_error("Unaligned TransposeFrom");

    if( g.InGrid() ) 
    { 
        const int r = g.Height();
        const int colShift = this->ColShift();

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
                    ALocalBuffer[jLocal+(colShift+iLocal*r)*ALDim];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR] = [MC,MR]");
    this->AssertNotLockedView();
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( this->Grid() == A.Grid() )
    {
        if( !this->Viewing() )
        {
            if( !this->ConstrainedColAlignment() )
            {
                this->_colAlignment = A.ColAlignment();
                if( this->Grid().InGrid() )
                    this->_colShift = A.ColShift();
            }
            if( !this->ConstrainedRowAlignment() )
            {
                this->_rowAlignment = A.RowAlignment();
                if( this->Grid().InGrid() )
                    this->_rowShift = A.RowShift();
            }
            this->ResizeTo( A.Height(), A.Width() );
        }

        if( this->Grid().InGrid() ) 
        {
            if( this->ColAlignment() == A.ColAlignment() &&
                this->RowAlignment() == A.RowAlignment() )
            {
                this->_localMatrix = A.LockedLocalMatrix();
            }
            else
            {
                const elemental::Grid& g = this->Grid();
#ifdef UNALIGNED_WARNINGS
                if( g.VCRank() == 0 )
                    std::cerr << "Unaligned [MC,MR] <- [MC,MR]." << std::endl;
#endif
                const int r = g.Height();
                const int c = g.Width();
                const int row = g.MCRank();
                const int col = g.MRRank();

                const int colAlignment = this->ColAlignment();
                const int rowAlignment = this->RowAlignment();
                const int colAlignmentOfA = A.ColAlignment();
                const int rowAlignmentOfA = A.RowAlignment();

                const int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
                const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
                const int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;
                const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;
                const int sendRank = sendRow + sendCol*r;
                const int recvRank = recvRow + recvCol*r;

                const int localHeight = this->LocalHeight();
                const int localWidth = this->LocalWidth();
                const int localHeightOfA = A.LocalHeight();
                const int localWidthOfA = A.LocalWidth();

                const int sendSize = localHeightOfA * localWidthOfA;
                const int recvSize = localHeight * localWidth;

                this->_auxMemory.Require( sendSize + recvSize );

                T* buffer = this->_auxMemory.Buffer();
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
                this->_auxMemory.Release();
            }
        } 
    }
    else // the grids don't match
    {
        if( !mpi::CongruentComms( A.Grid().ViewingComm(), 
                                  this->Grid().ViewingComm() ) )
            throw std::logic_error
            ("Redistributing between nonmatching grids currently requires"
             " the viewing communicators to match.");
                     
        if( !this->Viewing() )
            this->ResizeTo( A.Height(), A.Width() );

        // Compute the number of process rows and columns that each process 
        // needs to send to.
        const int r0 = this->Grid().Height();
        const int c0 = this->Grid().Width();
        const int rA = A.Grid().Height();
        const int cA = A.Grid().Width();
        const int myRow0 = this->Grid().MCRank();
        const int myCol0 = this->Grid().MRRank();
        const int myRowA = A.Grid().MCRank();
        const int myColA = A.Grid().MRRank();
        const int rowGCD = GCD( r0, rA );
        const int colGCD = GCD( c0, cA );
        const int rowLCM = r0*rA / rowGCD;
        const int colLCM = c0*cA / colGCD;
        const int numRowSends = r0 / rowGCD;
        const int numColSends = c0 / colGCD;
        const int localColStride0 = rowLCM / r0;
        const int localRowStride0 = colLCM / c0;
        const int localColStrideA = numRowSends;
        const int localRowStrideA = numColSends;

        const int colAlign0 = this->ColAlignment();
        const int rowAlign0 = this->RowAlignment();
        const int colAlignA = A.ColAlignment();
        const int rowAlignA = A.RowAlignment();

        const bool inThisGrid = this->Grid().InGrid();
        const bool inAGrid = A.Grid().InGrid();

        int maxSendSize = 
            (A.Height()/(rA*localColStrideA)+1) * 
            (A.Width()/(cA*localRowStrideA)+1);

        // Have each member of A's grid individually send to all numRow x numCol
        // processes in order, while the members of this grid receive from all 
        // necessary processes at each step.
        int requiredMemory = 0;
        if( inAGrid )
            requiredMemory += maxSendSize;
        if( inThisGrid )
            requiredMemory += maxSendSize;
        this->_auxMemory.Require( requiredMemory );
        T* buffer = this->_auxMemory.Buffer();
        int offset = 0;
        T* sendBuffer = &buffer[offset];
        if( inAGrid )
            offset += maxSendSize;
        T* recvBuffer = &buffer[offset];

        int recvRow = 0; // avoid compiler warnings...
        if( inAGrid )
            recvRow = (((myRowA+rA-colAlignA)%rA)+colAlign0)%r0;
        for( int rowSendCycle=0; rowSendCycle<numRowSends; ++rowSendCycle )
        {
            int recvCol = 0; // avoid compiler warnings...
            if( inAGrid )
                recvCol = (((myColA+cA-rowAlignA)%cA)+rowAlign0)%c0;

            for( int colSendCycle=0; colSendCycle<numColSends; ++colSendCycle )
            {
                mpi::Request sendRequest;

                // Fire off this round of non-blocking sends
                if( inAGrid )
                {
                    // Pack the data
                    int sendHeight = LocalLength
                        ( A.LocalHeight(), rowSendCycle, numRowSends );
                    int sendWidth = LocalLength
                        ( A.LocalWidth(), colSendCycle, numColSends );
                    const T* ALocalBuffer = A.LockedLocalBuffer();
                    const int ALDim = A.LocalLDim();
#ifdef _OPENMP
                    #pragma omp parallel for
#endif
                    for( int jLocal=0; jLocal<sendWidth; ++jLocal )
                    {
                        const int j = colSendCycle+jLocal*localRowStrideA;
                        for( int iLocal=0; iLocal<sendHeight; ++iLocal )
                        {
                            const int i = rowSendCycle+iLocal*localColStrideA;
                            sendBuffer[iLocal+jLocal*sendHeight] = 
                                ALocalBuffer[i+j*ALDim];
                        }
                    }
                    // Send data
                    int recvVCRank = recvRow + recvCol*r0;
                    int recvViewingRank = 
                        this->Grid().VCToViewingMap( recvVCRank );
                    mpi::ISend
                    ( sendBuffer, sendHeight*sendWidth, recvViewingRank,
                      0, this->Grid().ViewingComm(), sendRequest );
                }
                // Perform this round of recv's
                if( inThisGrid )
                {
                    const int sendRowOffset = (rowSendCycle*rA+colAlignA) % rA;
                    const int sendColOffset = (colSendCycle*cA+rowAlignA) % cA;
                    const int recvRowOffset = (rowSendCycle*rA+colAlign0) % r0;
                    const int recvColOffset = (colSendCycle*cA+rowAlign0) % c0;
                    const int firstSendRow = 
                        (((myRow0+r0-recvRowOffset)%r0)+sendRowOffset)%rA;
                    const int firstSendCol = 
                        (((myCol0+c0-recvColOffset)%c0)+sendColOffset)%cA;

                    const int rowShift = (myRow0+r0-recvRowOffset)%r0;
                    const int colShift = (myCol0+c0-recvColOffset)%c0;
                    const int numRowRecvs = LocalLength( rA, rowShift, r0 ); 
                    const int numColRecvs = LocalLength( cA, colShift, c0 );

                    // Recv data
                    // For now, simply receive sequentially. Until we switch to 
                    // nonblocking recv's, we won't be using much of the 
                    // recvBuffer
                    int sendRow = firstSendRow;
                    for( int rowRecvCycle=0; rowRecvCycle<numRowRecvs; 
                         ++rowRecvCycle )
                    {
                        const int sendRowShift = 
                            Shift( sendRow, colAlignA, rA ) + rowSendCycle*rA;
                        const int sendHeight = 
                            LocalLength( A.Height(), sendRowShift, rowLCM );
                        const int localColOffset = 
                            (sendRowShift-this->ColShift()) / r0;

                        int sendCol = firstSendCol;
                        for( int colRecvCycle=0; 
                             colRecvCycle<numColRecvs;  ++colRecvCycle )
                        {
                            const int sendColShift = 
                                Shift( sendCol, rowAlignA, cA ) + 
                                colSendCycle*cA;
                            const int sendWidth = 
                                LocalLength( A.Width(), sendColShift, colLCM );
                            const int localRowOffset = 
                                (sendColShift-this->RowShift()) / c0;

                            const int sendVCRank = sendRow+sendCol*rA;
                            const int sendViewingRank = 
                                A.Grid().VCToViewingMap( sendVCRank );

                            mpi::Recv
                            ( recvBuffer, sendHeight*sendWidth, sendViewingRank,
                              0, this->Grid().ViewingComm() );
                            
                            // Unpack the data
                            T* thisLocalBuffer = this->LocalBuffer();
                            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
                            #pragma omp parallel for
#endif
                            for( int jLocal=0; jLocal<sendWidth; ++jLocal )
                            {
                                const int j =
                                    localRowOffset+jLocal*localRowStride0;
                                for( int iLocal=0; iLocal<sendHeight; ++iLocal )
                                {
                                    const int i = 
                                        localColOffset+iLocal*localColStride0;
                                    thisLocalBuffer[i+j*thisLDim] = 
                                        recvBuffer[iLocal+jLocal*sendHeight];
                                }
                            }
                            // Set up the next send col
                            sendCol = (sendCol + c0) % cA;
                        }
                        // Set up the next send row
                        sendRow = (sendRow + r0) % rA;
                    }
                }
                // Ensure that this round of non-blocking sends completes
                if( inAGrid )
                {
                    mpi::Wait( sendRequest );
                    recvCol = (recvCol + cA) % c0;
                }
            }
            if( inAGrid )
                recvRow = (recvRow + rA) % r0;
        }
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR] = [MC,* ]");
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
            this->_colAlignment = A.ColAlignment();
            if( g.InGrid() )
            {
                this->_colShift = 
                    Shift( g.MCRank(), this->ColAlignment(), g.Height() );
            }
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( g.InGrid() ) 
    {
        if( this->ColAlignment() == A.ColAlignment() )
        {
            const int c = g.Width();
            const int rowShift = this->RowShift();

            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();

            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(rowShift+jLocal*c)*ALDim];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                memcpy( thisCol, ACol, localHeight*sizeof(T) );
            }
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR] <- [MC,* ]." << std::endl;
#endif
            const int r = g.Height();
            const int c = g.Width();
            const int rank = g.MCRank();
            const int rowShift = this->RowShift();
            const int colAlignment = this->ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();

            const int sendRank = (rank+r+colAlignment-colAlignmentOfA) % r;
            const int recvRank = (rank+r+colAlignmentOfA-colAlignment) % r;

            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();

            const int sendSize = localHeightOfA * localWidth;
            const int recvSize = localHeight * localWidth;

            this->_auxMemory.Require( sendSize + recvSize );

            T* buffer = this->_auxMemory.Buffer();
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
                const T* ACol = &ALocalBuffer[(rowShift+jLocal*c)*ALDim];
                T* sendBufferCol = &sendBuffer[jLocal*localHeightOfA];
                memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
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
                const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->_auxMemory.Release();
        }
    } 
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,MR]");
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
            this->_rowAlignment = A.RowAlignment();
            if( g.InGrid() )
            {
                this->_rowShift = 
                    Shift( g.MRRank(), this->RowAlignment(), g.Width() );
            }
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( g.InGrid() ) 
    {
        if( this->RowAlignment() == A.RowAlignment() )
        {
            const int r = g.Height();
            const int colShift = this->ColShift();

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
                        ALocalBuffer[(colShift+iLocal*r)+jLocal*ALDim];
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR] <- [* ,MR]." << std::endl;
#endif
            const int r = g.Height(); 
            const int c = g.Width();
            const int col = g.MRRank();
            const int colShift = this->ColShift();
            const int rowAlignment = this->RowAlignment();
            const int rowAlignmentOfA = A.RowAlignment();

            const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
            const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localWidthOfA = A.LocalWidth();

            const int sendSize = localHeight * localWidthOfA;
            const int recvSize = localHeight * localWidth;

            this->_auxMemory.Require( sendSize + recvSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];

            // Pack
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    sendBuffer[iLocal+jLocal*localHeight] = 
                        ALocalBuffer[(colShift+iLocal*r)+jLocal*ALDim];

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
                const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->_auxMemory.Release();
        }
    } 
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MC,MR] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MC,MR] = [* ,MD] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( g.InGrid() ) 
    {
        if( A.Width() == 1 )
        {
            if( !this->Viewing() )
                this->ResizeTo( A.Height(), 1 );

            const int r = g.Height();
            const int c = g.Width();
            const int p = g.Size();
            const int myRow = g.MCRank();
            const int myCol = g.MRRank();
            const int rankCM = g.VCRank();
            const int rankRM = g.VRRank();
            const int ownerCol = this->RowAlignment();
            const int ownerRow = A.RowAlignment();
            const int colAlignment = this->ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int colShift = this->ColShift();
            const int colShiftOfA = A.ColShift();

            const int height = A.Height();
            const int maxLocalHeight = MaxLocalLength(height,p);

            const int portionSize = std::max(maxLocalHeight,mpi::MIN_COLL_MSG);

            const int colShiftVC = Shift(rankCM,colAlignment,p);
            const int colShiftVROfA = Shift(rankRM,colAlignmentOfA,p);
            const int sendRankCM = (rankCM+(p+colShiftVROfA-colShiftVC)) % p;
            const int recvRankRM = (rankRM+(p+colShiftVC-colShiftVROfA)) % p;
            const int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

            this->_auxMemory.Require( (r+c)*portionSize );
            T* buffer = this->_auxMemory.Buffer();
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[c*portionSize];

            if( myRow == ownerRow )
            {
                // Pack
                const T* ALocalBuffer = A.LockedLocalBuffer();
#ifdef _OPENMP
                #pragma omp parallel for  
#endif
                for( int k=0; k<r; ++k )
                {
                    T* data = &recvBuf[k*portionSize];

                    const int shift = RawShift(myCol+c*k,colAlignmentOfA,p);
                    const int offset = (shift-colShiftOfA) / c;
                    const int thisLocalHeight = RawLocalLength(height,shift,p);

                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        data[iLocal] = ALocalBuffer[offset+iLocal*r];
                }
            }

            // A[VR,* ] <- A[MR,MC]
            mpi::Scatter
            ( recvBuf, portionSize, 
              sendBuf, portionSize, ownerRow, g.MCComm() );

            // A[VC,* ] <- A[VR,* ]
            mpi::SendRecv
            ( sendBuf, portionSize, sendRankCM, 0,
              recvBuf, portionSize, recvRankCM, mpi::ANY_TAG, g.VCComm() );

            // A[MC,MR] <- A[VC,* ]
            mpi::Gather
            ( recvBuf, portionSize, 
              sendBuf, portionSize, ownerCol, g.MRComm() );

            if( myCol == ownerCol )
            {
                // Unpack
                T* thisLocalBuffer = this->LocalBuffer();
#ifdef _OPENMP
                #pragma omp parallel for  
#endif
                for( int k=0; k<c; ++k )
                {
                    const T* data = &sendBuf[k*portionSize];

                    const int shift = RawShift(myRow+r*k,colAlignment,p);
                    const int offset = (shift-colShift) / r;
                    const int thisLocalHeight = RawLocalLength(height,shift,p);

                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        thisLocalBuffer[offset+iLocal*c] = data[iLocal];
                }
            }
            this->_auxMemory.Release();
        }
        else if( A.Height() == 1 )
        {
            if( !this->Viewing() )
                this->ResizeTo( 1, A.Width() );

            const int r = g.Height();
            const int c = g.Width();
            const int p = g.Size();
            const int myRow = g.MCRank();
            const int myCol = g.MRRank();
            const int rankCM = g.VCRank();
            const int rankRM = g.VRRank();
            const int ownerRow = this->ColAlignment();
            const int ownerCol = A.ColAlignment();
            const int rowAlignment = this->RowAlignment();
            const int rowAlignmentOfA = A.RowAlignment();
            const int rowShift = this->RowShift();
            const int rowShiftOfA = A.RowShift();

            const int width = A.Width();
            const int maxLocalWidth = MaxLocalLength(width,p);

            const int portionSize = std::max(maxLocalWidth,mpi::MIN_COLL_MSG);

            const int rowShiftVR = Shift(rankRM,rowAlignment,p);
            const int rowShiftVCOfA = Shift(rankCM,rowAlignmentOfA,p);
            const int sendRankRM = (rankRM+(p+rowShiftVCOfA-rowShiftVR)) % p;
            const int recvRankCM = (rankCM+(p+rowShiftVR-rowShiftVCOfA)) % p;
            const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

            this->_auxMemory.Require( (r+c)*portionSize );
            T* buffer = this->_auxMemory.Buffer();
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[r*portionSize];

            if( myCol == ownerCol )
            {
                // Pack
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const int ALDim = A.LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for  
#endif
                for( int k=0; k<c; ++k )
                {
                    T* data = &recvBuf[k*portionSize];

                    const int shift = RawShift(myRow+r*k,rowAlignmentOfA,p);
                    const int offset = (shift-rowShiftOfA) / r;
                    const int thisLocalWidth = RawLocalLength(width,shift,p);

                    for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                        data[jLocal] = ALocalBuffer[(offset+jLocal*c)*ALDim];
                }
            }

            // A[* ,VC] <- A[MR,MC]
            mpi::Scatter
            ( recvBuf, portionSize, 
              sendBuf, portionSize, ownerCol, g.MRComm() );

            // A[* ,VR] <- A[* ,VC]
            mpi::SendRecv
            ( sendBuf, portionSize, sendRankRM, 0,
              recvBuf, portionSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

            // A[MC,MR] <- A[* ,VR]
            mpi::Gather
            ( recvBuf, portionSize, 
              sendBuf, portionSize, ownerRow, g.MCComm() );
    
            if( myRow == ownerRow )
            {
                // Unpack
                T* thisLocalBuffer = this->LocalBuffer();
                const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for  
#endif
                for( int k=0; k<r; ++k )
                {
                    const T* data = &sendBuf[k*portionSize];

                    const int shift = RawShift(myCol+c*k,rowAlignment,p);
                    const int offset = (shift-rowShift) / c;
                    const int thisLocalWidth = RawLocalLength(width,shift,p);

                    for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                        thisLocalBuffer[(offset+jLocal*r)*thisLDim] = 
                            data[jLocal];
                }
            }

            this->_auxMemory.Release();
        }
        else
        {
            if( A.Height() >= A.Width() )
            {
                auto_ptr< DistMatrix<T,VR,STAR> > A_VR_STAR
                ( new DistMatrix<T,VR,STAR>(g) );

                *A_VR_STAR = A;

                auto_ptr< DistMatrix<T,VC,STAR> > A_VC_STAR
                ( new DistMatrix<T,VC,STAR>(true,this->ColAlignment(),g) );
                *A_VC_STAR = *A_VR_STAR;
                delete A_VR_STAR.release(); // lowers memory highwater

                *this = *A_VC_STAR;
            }
            else
            {
                auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
                ( new DistMatrix<T,STAR,VC>(g) );
                *A_STAR_VC = A;

                auto_ptr< DistMatrix<T,STAR,VR> > A_STAR_VR
                ( new DistMatrix<T,STAR,VR>(true,this->RowAlignment(),g) );
                *A_STAR_VR = *A_STAR_VC;
                delete A_STAR_VC.release(); // lowers memory highwater
    
                *this = *A_STAR_VR;
                this->ResizeTo( A_STAR_VR->Height(), A_STAR_VR->Width() );
            }
        }
    } 
    else
    {
        this->ResizeTo( A.Height(), A.Width() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();

    auto_ptr< DistMatrix<T,VR,STAR> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(g) );
    *A_VR_STAR = A;

    auto_ptr< DistMatrix<T,VC,STAR> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();

    auto_ptr< DistMatrix<T,STAR,VC> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC>(g) );
    *A_STAR_VC = A;

    auto_ptr< DistMatrix<T,STAR,VR> > A_STAR_VR
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
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [VC,* ]");
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
            this->_colAlignment = A.ColAlignment() % g.Height();
            if( g.InGrid() )
            {
                this->_colShift = 
                    Shift( g.MCRank(), this->ColAlignment(), g.Height() );
            }
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( g.InGrid() )
    {
        if( this->ColAlignment() == A.ColAlignment() % g.Height() )
        {
            const int r = g.Height();
            const int c = g.Width();
            const int p = r * c;
            const int row = g.MCRank();
            const int colShift = this->ColShift();
            const int rowAlignment = this->RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();

            const int height = this->Height();
            const int width = this->Width();
            const int localWidth = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();

            const int maxHeight = MaxLocalLength(height,p);
            const int maxWidth = MaxLocalLength(width,c);
            const int portionSize = 
                std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

            this->_auxMemory.Require( 2*c*portionSize );

            T* buffer = this->_auxMemory.Buffer();
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

                const int thisRowShift = RawShift(k,rowAlignment,c);
                const int thisLocalWidth = RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for 
#endif
                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*c)*ALDim];
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

                const int thisRank = row+k*r;
                const int thisColShift = RawShift(thisRank,colAlignmentOfA,p);
                const int thisColOffset = (thisColShift-colShift) / r;
                const int thisLocalHeight = 
                    RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const int i = thisColOffset + iLocal*c;
                        thisLocalBuffer[i+jLocal*thisLDim] = 
                            data[iLocal+jLocal*thisLocalHeight];
                    }
                }
            }
            this->_auxMemory.Release();
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR] <- [VC,* ]." << std::endl;
#endif
            const int r = g.Height();
            const int c = g.Width();
            const int p = r * c;
            const int row = g.MCRank();
            const int colShift = this->ColShift();
            const int colAlignment = this->ColAlignment();
            const int rowAlignment = this->RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();

            const int sendRow = (row+r+colAlignment-(colAlignmentOfA%r)) % r;
            const int recvRow = (row+r+(colAlignmentOfA%r)-colAlignment) % r;

            const int height = this->Height();
            const int width = this->Width();
            const int localWidth = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();

            const int maxHeight = MaxLocalLength(height,p);
            const int maxWidth = MaxLocalLength(width,c);
            const int portionSize = 
                std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

            this->_auxMemory.Require( 2*c*portionSize );

            T* buffer = this->_auxMemory.Buffer();
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

                const int thisRowShift = RawShift(k,rowAlignment,c);
                const int thisLocalWidth = RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*c)*ALDim];
                    T* dataCol = &data[jLocal*localHeightOfA];
                    memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                }
            }

            // SendRecv: properly align A[VC,*] via a trade in the column
            mpi::SendRecv
            ( secondBuffer, c*portionSize, sendRow, 0,
              firstBuffer,  c*portionSize, recvRow, 0, g.MCComm() );

            // AllToAll to gather all of the aligned A[VC,*] data into 
            // secondBuff.
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
                const int thisColShift = RawShift(thisRank,colAlignmentOfA,p);
                const int thisColOffset = (thisColShift-colShift) / r;
                const int thisLocalHeight = 
                    RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const int i = thisColOffset + iLocal*c;
                        thisLocalBuffer[i+jLocal*thisLDim] = 
                            data[iLocal+jLocal*thisLocalHeight];
                    }
                }
            }
            this->_auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,VC]");
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
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR> A_VC_STAR(true,this->ColAlignment(),g);

    A_VC_STAR = A;
    *this = A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,VR]");
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
            this->_rowAlignment = A.RowAlignment() % g.Width();
            if( g.InGrid() )
            {
                this->_rowShift = 
                    Shift( g.MRRank(), this->RowAlignment(), g.Width() );
            }
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( g.InGrid() )
    {
        if( this->RowAlignment() == A.RowAlignment() % g.Width() )
        {
            const int r = g.Height();
            const int c = g.Width();
            const int p = r * c;
            const int col = g.MRRank();
            const int rowShift = this->RowShift();
            const int colAlignment = this->ColAlignment();
            const int rowAlignmentOfA = A.RowAlignment();
    
            const int height = this->Height();
            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidthOfA = A.LocalWidth();

            const int maxHeight = MaxLocalLength(height,r);
            const int maxWidth = MaxLocalLength(width,p);
            const int portionSize = 
                std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

            this->_auxMemory.Require( 2*r*portionSize );

            T* buffer = this->_auxMemory.Buffer();
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

                const int thisColShift = RawShift(k,colAlignment,r);
                const int thisLocalHeight = 
                    RawLocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                {
                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const int i = thisColShift + iLocal*r;
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[i+jLocal*ALDim];
                    }
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
                const int thisRowShift = RawShift(thisRank,rowAlignmentOfA,p);
                const int thisRowOffset = (thisRowShift-rowShift) / c;
                const int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* dataCol = &data[jLocal*localHeight];
                    T* thisCol = 
                        &thisLocalBuffer[(thisRowOffset+jLocal*r)*thisLDim];
                    memcpy( thisCol, dataCol, localHeight*sizeof(T) );
                }
            }
            this->_auxMemory.Release();
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR] <- [* ,VR]." << std::endl;
#endif
            const int r = g.Height();
            const int c = g.Width();
            const int p = r * c;
            const int col = g.MRRank();
            const int rowShift = this->RowShift();
            const int colAlignment = this->ColAlignment();
            const int rowAlignment = this->RowAlignment();
            const int rowAlignmentOfA = A.RowAlignment();

            const int sendCol = (col+c+rowAlignment-(rowAlignmentOfA%c)) % c;
            const int recvCol = (col+c+(rowAlignmentOfA%c)-rowAlignment) % c;

            const int height = this->Height();
            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidthOfA = A.LocalWidth();
            
            const int maxHeight = MaxLocalLength(height,r);
            const int maxWidth = MaxLocalLength(width,p);
            const int portionSize = 
                std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

            this->_auxMemory.Require( 2*r*portionSize );

            T* buffer = this->_auxMemory.Buffer();
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

                const int thisColShift = RawShift(k,colAlignment,r);
                const int thisLocalHeight = 
                    RawLocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                {
                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const int i = thisColShift + iLocal*r;
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[i+jLocal*ALDim];
                    }
                }
            }

            // SendRecv: properly align A[*,VR] via a trade in the column
            mpi::SendRecv
            ( secondBuffer, r*portionSize, sendCol, 0,
              firstBuffer,  r*portionSize, recvCol, 0, g.MRComm() );

            // AllToAll to gather all of the aligned [*,VR] data into 
            // secondBuffer
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
                const int thisRowShift = RawShift(thisRank,rowAlignmentOfA,p);
                const int thisRowOffset = (thisRowShift-rowShift) / c;
                const int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* dataCol = &data[jLocal*localHeight];
                    T* thisCol = 
                        &thisLocalBuffer[(thisRowOffset+jLocal*r)*thisLDim];
                    memcpy( thisCol, dataCol, localHeight*sizeof(T) );
                }
            }
            this->_auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,MC,MR>&
DistMatrix<T,MC,MR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int r = this->Grid().Height();
    const int c = this->Grid().Width();
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
                ALocalBuffer[(colShift+iLocal*r)+(rowShift+jLocal*c)*ALDim];
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::SumScatterFrom( const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterFrom([MC,* ])");
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
            this->_colAlignment = A.ColAlignment();
            if( g.InGrid() )
            {
                this->_colShift = 
                    Shift( g.MCRank(), this->ColAlignment(), g.Height() );
            }
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( g.InGrid() )
    {
        if( this->ColAlignment() == A.ColAlignment() )
        {
            if( this->Width() == 1 )
            {
                const int rowAlignment = this->RowAlignment();
                const int myCol = g.MRRank();

                const int localHeight = this->LocalHeight();

                const int recvSize = std::max(localHeight,mpi::MIN_COLL_MSG);
                const int sendSize = recvSize;

                this->_auxMemory.Require( sendSize + recvSize );

                T* buffer = this->_auxMemory.Buffer();
                T* sendBuffer = &buffer[0];
                T* recvBuffer = &buffer[sendSize];

                // Pack 
                const T* ACol = A.LockedLocalBuffer(0,0);
                memcpy( sendBuffer, ACol, localHeight*sizeof(T) );

                // Reduce to rowAlignment
                mpi::Reduce
                ( sendBuffer, recvBuffer, sendSize, 
                  mpi::SUM, rowAlignment, g.MRComm() );

                if( myCol == rowAlignment )
                {
                    T* thisCol = this->LocalBuffer(0,0);
                    memcpy( thisCol, recvBuffer, localHeight*sizeof(T) );
                }

                this->_auxMemory.Release();
            }
            else
            {
                const int c = g.Width();
                const int rowAlignment = this->RowAlignment();
            
                const int width = this->Width();
                const int localHeight = this->LocalHeight();
                const int localWidth = this->LocalWidth();
                const int maxLocalWidth = MaxLocalLength(width,c);

                const int recvSize = 
                    std::max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
                const int sendSize = c * recvSize;

                this->_auxMemory.Require( sendSize + recvSize );

                T* buffer = this->_auxMemory.Buffer();
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

                    const int thisRowShift = RawShift( k, rowAlignment, c );
                    const int thisLocalWidth = 
                        RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                    #pragma omp parallel for
#endif
                    for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    {
                        const T* ACol = 
                            &ALocalBuffer[(thisRowShift+jLocal*c)*ALDim];
                        T* dataCol = &data[jLocal*localHeight];
                        memcpy( dataCol, ACol, localHeight*sizeof(T) );
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
                    const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                    T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                    memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
                }
                this->_auxMemory.Release();
            }
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned SumScatterFrom [MC,MR] <- [MC,* ]." 
                          << std::endl;
#endif
            if( this->Width() == 1 )
            {
                const int r = g.Height();
                const int rowAlignment = this->RowAlignment();
                const int myRow = g.MCRank();
                const int myCol = g.MRRank();

                const int height = this->Height();
                const int localHeight = this->LocalHeight();
                const int localHeightOfA = A.LocalHeight();
                const int maxLocalHeight = MaxLocalLength(height,r);

                const int portionSize = 
                    std::max(maxLocalHeight,mpi::MIN_COLL_MSG);

                const int colAlignment = this->ColAlignment();
                const int colAlignmentOfA = A.ColAlignment();
                const int sendRow = (myRow+r+colAlignment-colAlignmentOfA) % r;
                const int recvRow = (myRow+r+colAlignmentOfA-colAlignment) % r;

                this->_auxMemory.Require( 2*portionSize );

                T* buffer = this->_auxMemory.Buffer();
                T* sendBuffer = &buffer[0];
                T* recvBuffer = &buffer[portionSize];

                // Pack 
                const T* ACol = A.LockedLocalBuffer(0,0);
                memcpy( sendBuffer, ACol, localHeightOfA*sizeof(T) );
            
                // Reduce to rowAlignment
                mpi::Reduce
                ( sendBuffer, recvBuffer, portionSize, 
                  mpi::SUM, rowAlignment, g.MRComm() );

                if( myCol == rowAlignment )
                {
                    // Perform the realignment
                    mpi::SendRecv
                    ( recvBuffer, portionSize, sendRow, 0,
                      sendBuffer, portionSize, recvRow, 0, g.MCComm() );

                    T* thisCol = this->LocalBuffer(0,0);
                    memcpy( thisCol, sendBuffer, localHeight*sizeof(T) );
                }

                this->_auxMemory.Release();
            }
            else
            {
                const int r = g.Height();
                const int c = g.Width();
                const int row = g.MCRank();

                const int colAlignment = this->ColAlignment();
                const int rowAlignment = this->RowAlignment();
                const int colAlignmentOfA = A.ColAlignment();
                const int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
                const int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

                const int width = this->Width();
                const int localHeight = this->LocalHeight();
                const int localWidth = this->LocalWidth();
                const int localHeightOfA = A.LocalHeight();
                const int maxLocalWidth = MaxLocalLength(width,c);

                const int recvSize_RS = 
                    std::max(localHeightOfA*maxLocalWidth,mpi::MIN_COLL_MSG);
                const int sendSize_RS = c * recvSize_RS;
                const int recvSize_SR = localHeight * localWidth;

                this->_auxMemory.Require
                ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );

                T* buffer = this->_auxMemory.Buffer();
                T* firstBuffer = &buffer[0];
                T* secondBuffer = &buffer[recvSize_RS];

                // Pack 
                std::vector<int> recvSizes(c);
                // TODO: Stick an optional outer parallelization here?
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const int ALDim = A.LocalLDim();
                for( int k=0; k<c; ++k )
                {
                    T* data = &secondBuffer[k*recvSize_RS];
                    recvSizes[k] = recvSize_RS;

                    const int thisRowShift = RawShift( k, rowAlignment, c );
                    const int thisLocalWidth = 
                        RawLocalLength(width,thisRowShift,c);

#ifdef _OPENMP
                    #pragma omp parallel for
#endif
                    for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    {
                        const int j = thisRowShift + jLocal*c;
                        const T* ACol = &ALocalBuffer[j*ALDim];
                        T* dataCol = &data[jLocal*localHeightOfA];
                        memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                    }
                }

                // Reduce-scatter over each process row
                mpi::ReduceScatter
                ( secondBuffer, firstBuffer, &recvSizes[0], mpi::SUM, 
                  g.MRComm() );

                // Trade reduced data with the appropriate process row
                mpi::SendRecv
                ( firstBuffer,  localHeightOfA*localWidth, sendRow, 0,
                  secondBuffer, localHeight*localWidth,    recvRow, 0, 
                  g.MCComm() );

                // Unpack the received data
                T* thisLocalBuffer = this->LocalBuffer();
                const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    const T* secondBufferCol = 
                        &secondBuffer[jLocal*localHeight];
                    T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                    memcpy( thisCol, secondBufferCol, localHeight*sizeof(T) );
                }
                this->_auxMemory.Release();
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::SumScatterFrom( const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterFrom([* ,MR])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Width() == 1 && g.VCRank() == 0 )
    {
        std::cerr <<
          "The vector version of [MC,MR].SumScatterFrom([* ,MR]) does not "
          "yet have a vector version implemented, but it would only require"
          " a modification of the vector version of "
          "[MC,MR].SumScatterFrom([MC,* ])" << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.VCRank() == 0 )
    {
        std::cerr << 
          "[MC,MR]::SumScatterFrom([* ,MR]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MR] matrix instead." << std::endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            if( g.InGrid() )
            {
                this->_rowShift = 
                    Shift( g.MRRank(), this->RowAlignment(), g.Width() );
            }
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( g.InGrid() )
    {
        if( this->RowAlignment() == A.RowAlignment() )
        {
            const int r = g.Height();
            const int colAlignment = this->ColAlignment();

            const int height = this->Height();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int recvSize = 
                std::max(maxLocalHeight*localWidth,mpi::MIN_COLL_MSG);
            const int sendSize = r * recvSize;

            this->_auxMemory.Require( sendSize + recvSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];

            // Pack 
            std::vector<int> recvSizes(r);
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &sendBuffer[k*recvSize];
                recvSizes[k] = recvSize;

                const int thisColShift = RawShift( k, colAlignment, r );
                const int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const int i = thisColShift + iLocal*r;
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[i+jLocal*ALDim];
                    }
                }
            }

            // Reduce-scatter over each process col
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
                const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->_auxMemory.Release();
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned SumScatterFrom [MC,MR] <- [* ,MR]." 
                          << std::endl;
#endif
            const int r = g.Height();
            const int c = g.Width();
            const int col = g.MRRank();

            const int colAlignment = this->ColAlignment();
            const int rowAlignment = this->RowAlignment();
            const int rowAlignmentOfA = A.RowAlignment();
            const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
            const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

            const int height = this->Height();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localWidthOfA = A.LocalWidth();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int recvSize_RS = 
                std::max(maxLocalHeight*localWidthOfA,mpi::MIN_COLL_MSG);
            const int sendSize_RS = r * recvSize_RS;
            const int recvSize_SR = localHeight * localWidth;

            this->_auxMemory.Require
                ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );

            T* buffer = this->_auxMemory.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[recvSize_RS];

            // Pack 
            std::vector<int> recvSizes(r);
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &secondBuffer[k*recvSize_RS];
                recvSizes[k] = recvSize_RS;

                const int thisColShift = RawShift( k, colAlignment, r );
                const int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                {
                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const int i = thisColShift + iLocal*r;
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[i+jLocal*ALDim];
                    }
                }
            }

            // Reduce-scatter over each process col
            mpi::ReduceScatter
            ( secondBuffer, firstBuffer, &recvSizes[0], mpi::SUM, g.MCComm() );

            // Trade reduced data with the appropriate process col
            mpi::SendRecv
            ( firstBuffer,  localHeight*localWidthOfA, sendCol, 0,
              secondBuffer, localHeight*localWidth,    recvCol, 0, g.MRComm() );

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
            this->_auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::SumScatterFrom
( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterFrom([* ,* ])");
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
        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);
        const int maxLocalWidth = MaxLocalLength(width,c);

        const int recvSize = 
            std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
        const int sendSize = r*c*recvSize;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack 
        std::vector<int> recvSizes(r*c);
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int l=0; l<c; ++l )
        {
            const int thisRowShift = RawShift( l, rowAlignment, c );
            const int thisLocalWidth = RawLocalLength( width, thisRowShift, c );

            for( int k=0; k<r; ++k )
            {
                T* data = &sendBuffer[(k+l*r)*recvSize];
                recvSizes[k+l*r] = recvSize;

                const int thisColShift = RawShift( k, colAlignment, r );
                const int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[(thisColShift+iLocal*r)+
                                         (thisRowShift+jLocal*c)*ALDim];
            }
        }

        // Reduce-scatter over each process col
        mpi::ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.VCComm() );

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
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::SumScatterUpdate
( T alpha, const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterUpdate([MC,* ])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        if( this->ColAlignment() == A.ColAlignment() )
        {
            if( this->Width() == 1 )
            {
                const int rowAlignment = this->RowAlignment();
                const int myCol = g.MRRank();

                const int localHeight = this->LocalHeight();

                const int portionSize = std::max(localHeight,mpi::MIN_COLL_MSG);

                this->_auxMemory.Require( 2*portionSize );

                T* buffer = this->_auxMemory.Buffer();
                T* sendBuffer = &buffer[0];
                T* recvBuffer = &buffer[portionSize];

                // Pack 
                const T* ACol = A.LockedLocalBuffer(0,0);
                memcpy( sendBuffer, ACol, localHeight*sizeof(T) );
            
                // Reduce to rowAlignment
                mpi::Reduce
                ( sendBuffer, recvBuffer, portionSize, 
                  mpi::SUM, rowAlignment, g.MRComm() );

                if( myCol == rowAlignment )
                {
                    T* thisCol = this->LocalBuffer(0,0);
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
                    #pragma omp parallel for
#endif
                    for( int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisCol[iLocal] += alpha*recvBuffer[iLocal];
                }

                this->_auxMemory.Release();
            }
            else
            {
                const int c = g.Width();
                const int rowAlignment = this->RowAlignment();

                const int width = this->Width();
                const int localHeight = this->LocalHeight();
                const int localWidth = this->LocalWidth();
                const int maxLocalWidth = MaxLocalLength(width,c);

                const int portionSize = 
                    std::max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

                this->_auxMemory.Require( (c+1)*portionSize );

                T* buffer = this->_auxMemory.Buffer();
                T* sendBuffer = &buffer[0];
                T* recvBuffer = &buffer[c*portionSize];

                // Pack 
                std::vector<int> recvSizes(c);
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int k=0; k<c; ++k )
                {
                    T* data = &sendBuffer[k*portionSize];
                    recvSizes[k] = portionSize;

                    const int thisRowShift = RawShift( k, rowAlignment, c );
                    const int thisLocalWidth = 
                        RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                    #pragma omp parallel for
#endif
                    for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    {
                        const T* ACol = 
                            &ALocalBuffer[(thisRowShift+jLocal*c)*ALDim];
                        T* dataCol = &data[jLocal*localHeight];
                        memcpy( dataCol, ACol, localHeight*sizeof(T) );
                    }
                }
            
                // Reduce-scatter over each process row
                mpi::ReduceScatter
                ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.MRComm() );

                // Update with our received data
                T* thisLocalBuffer = this->LocalBuffer();
                const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                    T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                    blas::Axpy
                    ( localHeight, alpha, recvBufferCol, 1, thisCol, 1 );
                }
                this->_auxMemory.Release();
            }
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
            {
                std::cerr << "Unaligned SumScatterUpdate [MC,MR] <- [MC,* ]." 
                          << std::endl;
            }
#endif
            if( this->Width() == 1 )
            {
                const int r = g.Height();
                const int rowAlignment = this->RowAlignment();
                const int myRow = g.MCRank();
                const int myCol = g.MRRank();

                const int height = this->Height();
                const int localHeight = this->LocalHeight();
                const int localHeightOfA = A.LocalHeight();
                const int maxLocalHeight = MaxLocalLength(height,r);

                const int portionSize = 
                    std::max(maxLocalHeight,mpi::MIN_COLL_MSG);

                const int colAlignment = this->ColAlignment();
                const int colAlignmentOfA = A.ColAlignment();
                const int sendRow = (myRow+r+colAlignment-colAlignmentOfA) % r;
                const int recvRow = (myRow+r+colAlignmentOfA-colAlignment) % r;

                this->_auxMemory.Require( 2*portionSize );

                T* buffer = this->_auxMemory.Buffer();
                T* sendBuffer = &buffer[0];
                T* recvBuffer = &buffer[portionSize];

                // Pack 
                const T* ACol = A.LockedLocalBuffer(0,0);
                memcpy( sendBuffer, ACol, localHeightOfA*sizeof(T) );
            
                // Reduce to rowAlignment
                mpi::Reduce
                ( sendBuffer, recvBuffer, portionSize, 
                  mpi::SUM, rowAlignment, g.MRComm() );

                if( myCol == rowAlignment )
                {
                    // Perform the realignment
                    mpi::SendRecv
                    ( recvBuffer, portionSize, sendRow, 0,
                      sendBuffer, portionSize, recvRow, 0, g.MCComm() );

                    T* thisCol = this->LocalBuffer(0,0);
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
                    #pragma omp parallel for
#endif
                    for( int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisCol[iLocal] += alpha*sendBuffer[iLocal];
                }
                this->_auxMemory.Release();
            }
            else
            {
                const int r = g.Height();
                const int c = g.Width();
                const int row = g.MCRank();

                const int colAlignment = this->ColAlignment();
                const int rowAlignment = this->RowAlignment();
                const int colAlignmentOfA = A.ColAlignment();
                const int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
                const int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

                const int width = this->Width();
                const int localHeight = this->LocalHeight();
                const int localWidth = this->LocalWidth();
                const int localHeightOfA = A.LocalHeight();
                const int maxLocalWidth = MaxLocalLength(width,c);

                const int recvSize_RS = 
                    std::max(localHeightOfA*maxLocalWidth,mpi::MIN_COLL_MSG);
                const int sendSize_RS = c * recvSize_RS;
                const int recvSize_SR = localHeight * localWidth;

                this->_auxMemory.Require
                ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );

                T* buffer = this->_auxMemory.Buffer();
                T* firstBuffer = &buffer[0];
                T* secondBuffer = &buffer[recvSize_RS];

                // Pack 
                std::vector<int> recvSizes(c);
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int k=0; k<c; ++k )
                {
                    T* data = &secondBuffer[k*recvSize_RS];
                    recvSizes[k] = recvSize_RS;

                    const int thisRowShift = RawShift( k, rowAlignment, c );
                    const int thisLocalWidth = 
                        RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                    #pragma omp parallel for
#endif
                    for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    {
                        const T* ACol = 
                            &ALocalBuffer[(thisRowShift+jLocal*c)*ALDim];
                        T* dataCol = &data[jLocal*localHeightOfA];
                        memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                    }
                }

                // Reduce-scatter over each process row
                mpi::ReduceScatter
                ( secondBuffer, firstBuffer, &recvSizes[0], mpi::SUM, 
                  g.MRComm() );

                // Trade reduced data with the appropriate process row
                mpi::SendRecv
                ( firstBuffer,  localHeightOfA*localWidth, sendRow, 0,
                  secondBuffer, localHeight*localWidth,    recvRow, 0, 
                  g.MCComm() );

                // Update with our received data
                T* thisLocalBuffer = this->LocalBuffer();
                const int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
                #pragma omp parallel for
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    const T* secondBufferCol = 
                        &secondBuffer[jLocal*localHeight];
                    T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                    for( int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisCol[iLocal] += alpha*secondBufferCol[iLocal];
                }
                this->_auxMemory.Release();
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterUpdate([* ,MR])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Width() == 1 && g.VCRank() == 0 )
    {
        std::cerr <<
          "The vector version of [MC,MR].SumScatterUpdate([* ,MR]) does not"
          " yet have a vector version implemented, but it would only "
          "require a modification of the vector version of "
          "[MC,MR].SumScatterUpdate([MC,* ])" << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.VCRank() == 0 )
    {
        std::cerr << 
          "[MC,MR]::SumScatterUpdate([* ,MR]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MR] matrix instead." << std::endl;
    }
#endif
    if( g.InGrid() )
    {
        if( this->RowAlignment() == A.RowAlignment() )
        {
            const int r = g.Height();
            const int colAlignment = this->ColAlignment();

            const int height = this->Height();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int recvSize = 
                std::max(maxLocalHeight*localWidth,mpi::MIN_COLL_MSG);
            const int sendSize = r*recvSize;

            this->_auxMemory.Require( sendSize + recvSize );

            T* buffer = this->_auxMemory.Buffer();
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

                const int thisColShift = RawShift( k, colAlignment, r );
                const int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[(thisColShift+iLocal*r)+jLocal*ALDim];
            }

            // Reduce-scatter over each process col
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
            this->_auxMemory.Release();
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
            {
                std::cerr << "Unaligned SumScatterUpdate [MC,MR] <- [* ,MR]." 
                          << std::endl;
            }
#endif
            const int r = g.Height();
            const int c = g.Width();
            const int col = g.MRRank();

            const int colAlignment = this->ColAlignment();
            const int rowAlignment = this->RowAlignment();
            const int rowAlignmentOfA = A.RowAlignment();
            const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
            const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

            const int height = this->Height();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localWidthOfA = A.LocalWidth();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int recvSize_RS = 
                std::max(maxLocalHeight*localWidthOfA,mpi::MIN_COLL_MSG);
            const int sendSize_RS = r * recvSize_RS;
            const int recvSize_SR = localHeight * localWidth;

            this->_auxMemory.Require
                ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );

            T* buffer = this->_auxMemory.Buffer();
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

                const int thisColShift = RawShift( k, colAlignment, r );
                const int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[(thisColShift+iLocal*r)+jLocal*ALDim];
            }

            // Reduce-scatter over each process col
            mpi::ReduceScatter
            ( secondBuffer, firstBuffer, &recvSizes[0], mpi::SUM, g.MCComm() );

            // Trade reduced data with the appropriate process col
            mpi::SendRecv
            ( firstBuffer,  localHeight*localWidthOfA, sendCol, 0,
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
            this->_auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MC,MR>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterUpdate([* ,* ])");
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
        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);
        const int maxLocalWidth = MaxLocalLength(width,c);

        const int recvSize = 
            std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
        const int sendSize = r * c * recvSize;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack 
        std::vector<int> recvSizes(r*c);
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int l=0; l<c; ++l )
        {
            const int thisRowShift = RawShift( l, rowAlignment, c );
            const int thisLocalWidth = RawLocalLength( width, thisRowShift, c );

            for( int k=0; k<r; ++k )
            {
                T* data = &sendBuffer[(k+l*r)*recvSize];
                recvSizes[k+l*r] = recvSize;

                const int thisColShift = RawShift( k, colAlignment, r );
                const int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[(thisColShift+iLocal*r)+
                                         (thisRowShift+jLocal*c)*ALDim];
            }
        }

        // Reduce-scatter over each process col
        mpi::ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.VCComm() );

        // Unpack our received data
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
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
