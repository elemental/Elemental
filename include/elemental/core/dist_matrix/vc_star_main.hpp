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
DistMatrix<T,VC,STAR>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;

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

    std::vector<T> sendBuf(height*width,0);
    const T* thisLocalBuffer = this->LockedLocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        for( int j=0; j<width; ++j )
            sendBuf[(colShift+iLocal*p)+j*height] = 
                thisLocalBuffer[iLocal+j*thisLDim];

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
DistMatrix<T,VC,STAR>::Align( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::Align");
    this->AssertFreeColAlignment();
#endif
    this->AlignCols( colAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VC,STAR>::AlignCols( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Size() )
        throw std::runtime_error("Invalid column alignment for [VC,* ]");
#endif
    this->colAlignment_ = colAlignment;
    this->colShift_ = Shift( g.VCRank(), colAlignment, g.Size() );
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
DistMatrix<T,VC,STAR>::View( DistMatrix<T,VC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::View");
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
DistMatrix<T,VC,STAR>::View
( int height, int width, int colAlignment,
  T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->colShift_ = Shift(grid.VCRank(),colAlignment,grid.Size());
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
DistMatrix<T,VC,STAR>::LockedView( const DistMatrix<T,VC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::LockedView");
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
DistMatrix<T,VC,STAR>::LockedView
( int height, int width, int colAlignment,
  const T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->colShift_ = Shift(grid.VCRank(),colAlignment,grid.Size());
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
DistMatrix<T,VC,STAR>::View
( DistMatrix<T,VC,STAR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::View");
    this->AssertFreeColAlignment();
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

        this->colAlignment_ = (A.ColAlignment()+i) % size;
        this->colShift_ = Shift( colMajorRank, this->ColAlignment(), size );

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
DistMatrix<T,VC,STAR>::LockedView
( const DistMatrix<T,VC,STAR>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::LockedView");
    this->AssertFreeColAlignment();
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

        this->colAlignment_ = (A.ColAlignment()+i) % size;
        this->colShift_ = Shift( colMajorRank, this->ColAlignment(), size );

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
DistMatrix<T,VC,STAR>::View1x2
( DistMatrix<T,VC,STAR>& AL, DistMatrix<T,VC,STAR>& AR )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::View1x2");
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
DistMatrix<T,VC,STAR>::LockedView1x2
( const DistMatrix<T,VC,STAR>& AL, const DistMatrix<T,VC,STAR>& AR )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::LockedView1x2");
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
DistMatrix<T,VC,STAR>::View2x1
( DistMatrix<T,VC,STAR>& AT,
  DistMatrix<T,VC,STAR>& AB )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::View2x1");
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
DistMatrix<T,VC,STAR>::LockedView2x1
( const DistMatrix<T,VC,STAR>& AT,
  const DistMatrix<T,VC,STAR>& AB )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::LockedView2x1");
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
DistMatrix<T,VC,STAR>::View2x2
( DistMatrix<T,VC,STAR>& ATL, DistMatrix<T,VC,STAR>& ATR,
  DistMatrix<T,VC,STAR>& ABL, DistMatrix<T,VC,STAR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::View2x2");
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
DistMatrix<T,VC,STAR>::LockedView2x2
( const DistMatrix<T,VC,STAR>& ATL, const DistMatrix<T,VC,STAR>& ATR,
  const DistMatrix<T,VC,STAR>& ABL, const DistMatrix<T,VC,STAR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::LockedView2x2");
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
DistMatrix<T,VC,STAR>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
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
DistMatrix<T,VC,STAR>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        u = this->GetLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,VC,STAR>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
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
DistMatrix<T,VC,STAR>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        this->UpdateLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VC,STAR>::GetDiagonal
( DistMatrix<T,VC,STAR>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Height() || d.Width() != 1) )
    {
        std::ostringstream msg;
        msg << "d is not a column vec of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( (d.Viewing() || d.ConstrainedColAlignment() ) &&
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

    if( g.InGrid() )
    {
        const int p = g.Size();
        const int colShift = this->ColShift();
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

        const int iLocalStart = (iStart-colShift) / p;
        const int localDiagLength = d.LocalHeight();
        T* dLocalBuffer = d.LocalBuffer();
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
            dLocalBuffer[k] = thisLocalBuffer[iLocal+jLocal*thisLDim];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VC,STAR>::GetDiagonal
( DistMatrix<T,STAR,VC>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Width() || d.Height() != 1) )
    {
        std::ostringstream msg;
        msg << "d is not a row vec of the same length as the diagonal:\n"
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

    if( g.InGrid() )
    {
        const int p = g.Size();
        const int colShift = this->ColShift();
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

        const int iLocalStart = (iStart-colShift) / p;
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
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
            dLocalBuffer[k*dLDim] = thisLocalBuffer[iLocal+jLocal*thisLDim];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VC,STAR>::SetDiagonal
( const DistMatrix<T,VC,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetDiagonal");
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
    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int colShift = this->ColShift();
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

        const int iLocalStart = (iStart-colShift)/p;
        const int localDiagLength = d.LocalHeight();
        const T* dLocalBuffer = d.LockedLocalBuffer();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
            thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VC,STAR>::SetDiagonal
( const DistMatrix<T,STAR,VC>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetDiagonal");
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
    const elemental::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int colShift = this->ColShift();
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

        const int iLocalStart = (iStart-colShift)/p;
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
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
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
DistMatrix<T,VC,STAR>::MakeTrapezoidal
( Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->Grid().Size();
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
            int firstZeroRow = 
                ( side==LEFT ? std::max(j-offset+1,0)
                             : std::max(j-offset+height-width+1,0) );
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
DistMatrix<T,VC,STAR>::ScaleTrapezoidal
( T alpha, Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->Grid().Size();
    const int colShift = this->ColShift();

    if( uplo == UPPER )
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            int boundary = std::min( lastRow+1, height );
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
            int firstRow = ( side==LEFT ? std::max(j-offset,0)
                                        : std::max(j-offset+height-width,0) );
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
DistMatrix<T,VC,STAR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int width       = this->Width();
    const int localHeight = this->LocalHeight();
    const int p           = this->Grid().Size();
    const int colShift    = this->ColShift();

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
DistMatrix<T,VC,STAR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetToRandom");
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
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [MC,MR]");
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
                Shift( g.VCRank(), this->ColAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() % g.Height() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,c);
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
            const int thisColShift = RawShift(thisRank,colAlignment,p); 
            const int thisColOffset = (thisColShift-colShiftOfA) / r;
            const int thisLocalHeight = RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColOffset+iLocal*c)+jLocal*ALDim];
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

            const int thisRowShift = RawShift(k,rowAlignmentOfA,c);
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(thisRowShift+jLocal*c)*thisLDim];
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [VC,* ] <- [MC,MR]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        
        const int sendRow = (row+r+(colAlignment%r)-colAlignmentOfA) % r;
        const int recvRow = (row+r+colAlignmentOfA-(colAlignment%r)) % r;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,c);
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
            const int thisColShift = RawShift(thisRank,colAlignment,p);
            const int thisColOffset = (thisColShift-colShiftOfA) / r;
            const int thisLocalHeight = RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColOffset+iLocal*c)+jLocal*ALDim];
        }

        // AllToAll to gather all of the unaligned [VC,*] data into firstBuffer
        mpi::AllToAll
        ( secondBuffer, portionSize, 
          firstBuffer,  portionSize, g.MRComm() );

        // SendRecv: properly align the [VC,*] via a trade in the column
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

            const int thisRowShift = RawShift(k,rowAlignmentOfA,c);
            const int thisLocalWidth = RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(thisRowShift+jLocal*c)*thisLDim];
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
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [MC,* ]");
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
                Shift( g.VCRank(), this->ColAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() % g.Height() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int colShift = this->ColShift();
        const int colShiftOfA = A.ColShift();
        const int colOffset = (colShift-colShiftOfA) / r;
        
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
                    ALocalBuffer[(colOffset+iLocal*c)+j*ALDim];
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [VC,* ] <- [MC,* ]." << std::endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int col = g.MRRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[VC,*] within our process column to fix alignments.
        const int sendRow = (row+r+(colAlignment%r)-colAlignmentOfA) % r;
        const int recvRow = (row+r+colAlignmentOfA-(colAlignment%r)) % r;
        const int sendRank = sendRow + r*col;

        const int sendColShift = Shift( sendRank, colAlignment, p );
        const int sendColOffset = (sendColShift-colShiftOfA) / r;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfSend = LocalLength(height,sendColShift,p);

        const int sendSize = localHeightOfSend * width;
        const int recvSize = localHeight * width;

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
        for( int j=0; j<width; ++j )
            for( int iLocal=0; iLocal<localHeightOfSend; ++iLocal )
                sendBuffer[iLocal+j*localHeightOfSend] = 
                    ALocalBuffer[(sendColOffset+iLocal*c)+j*ALDim];

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
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,MR]");
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
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[VC,* ] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,STAR,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[VC,* ] = [* ,MD] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,VR,STAR> A_VR_STAR(g);

    A_VR_STAR = A;
    *this = A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,VR,STAR> A_VR_STAR(g);

    A_VR_STAR = A;
    *this = A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
    ( new DistMatrix<T,MR,MC>(g) );
    *A_MR_MC = A;

    std::auto_ptr< DistMatrix<T,VR,STAR> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR>(g) );
    *A_VR_STAR = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [VC,* ]");
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
            std::cerr << "Unaligned [VC,* ] <- [VC,* ]." << std::endl;
#endif
        const int rank = g.VCRank();
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
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.VCComm() );

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
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
    ( new DistMatrix<T,MR,MC>(g) );
    *A_MR_MC = A;

    std::auto_ptr< DistMatrix<T,VC,STAR> > A_VR_STAR
    ( new DistMatrix<T,VC,STAR>(g) );
    *A_VR_STAR = *A_MR_MC; 
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localHeightOfA = A.LocalHeight();

    const int sendSize = localHeightOfA * width;
    const int recvSize = localHeight * width;

    const elemental::Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();
    const int rankCM = g.VCRank();
    const int rankRM = g.VRRank();

    const int colShift = this->ColShift();
    const int colShiftOfA = A.ColShift();

    // Compute which colmajor rank has the colShift equal to our colShiftOfA
    const int sendRankCM = (rankCM+(p+colShiftOfA-colShift)) % p;

    // Compute which colmajor rank has the A colShift that we need
    const int recvRankRM = (rankRM+(p+colShift-colShiftOfA)) % p;
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
    for( int j=0; j<width; ++j )
    {
        const T* ACol = &ALocalBuffer[j*ALDim];
        T* sendBufferCol = &sendBuffer[j*localHeightOfA];
        memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
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
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,VR]");
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
inline const DistMatrix<T,VC,STAR>&
DistMatrix<T,VC,STAR>::operator=( const DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,* ]");
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
DistMatrix<T,VC,STAR>::SumScatterFrom
( const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ]::SumScatterFrom( [MC,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && A.Grid().VCRank() == 0 )
    {
        std::cerr <<
          "[VC,* ]::SumScatterFrom([MC,* ]) potentially causes a large amount "
          "of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MC,* ] matrix instead." << std::endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            this->colShift_ = 
                Shift( g.VCRank(), this->ColAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() % g.Height() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int row = g.MCRank();
        const int colAlignment = this->ColAlignment();
        const int colShiftOfA = A.ColShift();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int maxLocalHeight = MaxLocalLength( height, p );

        const int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
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
                        ALocalBuffer[(thisColOffset+iLocal*c)+j*ALDim];
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
        throw std::logic_error
        ("Unaligned [VC,* ]::ReduceScatterFrom( [MC,* ] ) not implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VC,STAR>::SumScatterFrom
( const DistMatrix<T,STAR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ]::SumScatterFrom( [* ,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int p = g.Size();
    const int VCRank = g.VCRank();
    const int colAlignment = this->ColAlignment();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int maxLocalHeight = MaxLocalLength( height, p );

    const int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
    const int sendSize = p*recvSize;

    this->auxMemory_.Require( sendSize + recvSize );

    T* buffer = this->auxMemory_.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack
    std::vector<int> recvSizes(p);
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
    ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.VCComm() );

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
DistMatrix<T,VC,STAR>::SumScatterUpdate
( T alpha, const DistMatrix<T,MC,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ]::SumScatterUpdate( [MC,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && A.Grid().VCRank() == 0 )
    {
        std::cerr <<
          "[VC,* ]::SumScatterUpdate([MC,* ]) potentially causes a large amount"
          " of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MC,* ] matrix instead." << std::endl;
    }
#endif
    if( this->ColAlignment() % g.Height() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int row = g.MCRank();
        const int colAlignment = this->ColAlignment();
        const int colShiftOfA = A.ColShift();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int maxLocalHeight = MaxLocalLength( height, p );

        const int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
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
                        ALocalBuffer[(thisColOffset+iLocal*c)+j*ALDim];
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
        throw std::logic_error
        ("Unaligned [VC,* ]::ReduceScatterUpdate( [MC,* ] ) not implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,VC,STAR>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ]::SumScatterUpdate( [* ,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();

    const int p = g.Size();
    const int VCRank = g.VCRank();
    const int colAlignment = this->ColAlignment();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int maxLocalHeight = MaxLocalLength( height, p );

    const int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
    const int sendSize = p*recvSize;

    this->auxMemory_.Require( sendSize + recvSize );

    T* buffer = this->auxMemory_.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack
    std::vector<int> recvSizes(p);
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
    ( sendBuffer, recvBuffer, &recvSizes[0], mpi::SUM, g.VCComm() );

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
