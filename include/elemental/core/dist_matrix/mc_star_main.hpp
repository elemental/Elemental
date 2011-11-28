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

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;
        
    const Int height      = this->Height();
    const Int width       = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int r           = g.Height();
    const Int colShift    = this->ColShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( g.InGrid() )
    {
        // Only one process column needs to participate
        if( g.MRRank() == 0 )
        {
            std::vector<T> sendBuf(height*width,0);
            const T* thisLocalBuffer = this->LockedLocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                for( Int jLocal=0; jLocal<width; ++jLocal )
                    sendBuf[(colShift+iLocal*r)+jLocal*height] = 
                        thisLocalBuffer[iLocal+jLocal*thisLDim];

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
                for( Int i=0; i<height; ++i )
                {
                    for( Int j=0; j<width; ++j )
                        os << WrapScalar(recvBuf[i+j*height]) << " ";
                    os << "\n";
                }
                os << std::endl;
            }
        }
        mpi::Barrier( g.VCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::Align( Int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Align");
    this->AssertFreeColAlignment();
#endif
    this->AlignCols( colAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::AlignCols( Int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Height() )
        throw std::logic_error("Invalid column alignment for [MC,* ]");
#endif
    this->colAlignment_ = colAlignment;
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->colShift_ = Shift( g.MCRank(), colAlignment, g.Height() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::View( DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.ColShift();
        this->localMatrix_.View( A.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::View
( Int height, Int width, Int colAlignment,
  T* buffer, Int ldim, const elemental::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->grid_->InGrid() )
    {
        this->colShift_ = Shift(g.MCRank(),colAlignment,g.Height());
        const Int localHeight = LocalLength(height,this->colShift_,g.Height());
        this->localMatrix_.View( localHeight, width, buffer, ldim );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::LockedView( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.ColShift();
        this->localMatrix_.LockedView( A.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::LockedView
( Int height, Int width, Int colAlignment,
  const T* buffer, Int ldim, const elemental::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->viewing_ = true;
    this->lockedView_ = true; 
    if( this->grid_->InGrid() )
    {
        this->colShift_ = Shift(g.MCRank(),colAlignment,g.Height());
        const Int localHeight = LocalLength(height,this->colShift_,g.Height());
        this->localMatrix_.LockedView( localHeight, width, buffer, ldim );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::View
( DistMatrix<T,MC,STAR,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = false;

    const Int r = this->Grid().Height();
    const Int row = this->Grid().MCRank();

    this->colAlignment_ = (A.ColAlignment()+i) % r;

    if( this->Grid().InGrid() )
    {
        this->colShift_ = Shift( row, this->ColAlignment(), r );
        const Int localHeightBefore = LocalLength( i, A.ColShift(), r );
        const Int localHeight = LocalLength( height, this->colShift_, r );
        this->localMatrix_.View
        ( A.LocalMatrix(), localHeightBefore, j, localHeight, width );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::LockedView
( const DistMatrix<T,MC,STAR,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = true;

    const Int r   = this->Grid().Height();
    const Int row = this->Grid().MCRank();

    this->colAlignment_ = (A.ColAlignment()+i) % r;

    if( this->Grid().InGrid() )
    {
        this->colShift_ = Shift( row, this->ColAlignment(), r );
        const Int localHeightBefore = LocalLength( i, A.ColShift(), r );
        const Int localHeight = LocalLength( height, this->colShift_, r );

        this->localMatrix_.LockedView
        ( A.LockedLocalMatrix(), localHeightBefore, j, localHeight, width );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::View1x2
( DistMatrix<T,MC,STAR,Int>& AL, DistMatrix<T,MC,STAR,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::View1x2");    
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = AL.ColShift();
        this->localMatrix_.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::LockedView1x2
( const DistMatrix<T,MC,STAR,Int>& AL, const DistMatrix<T,MC,STAR,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedView1x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = AL.ColShift();
        this->localMatrix_.LockedView1x2
        ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::View2x1
( DistMatrix<T,MC,STAR,Int>& AT,
  DistMatrix<T,MC,STAR,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::View2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = AT.ColShift();
        this->localMatrix_.View2x1
        ( AT.LocalMatrix(),
          AB.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::LockedView2x1
( const DistMatrix<T,MC,STAR,Int>& AT,
  const DistMatrix<T,MC,STAR,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedView2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = AT.ColShift();
        this->localMatrix_.LockedView2x1
        ( AT.LockedLocalMatrix(),
          AB.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::View2x2
( DistMatrix<T,MC,STAR,Int>& ATL, DistMatrix<T,MC,STAR,Int>& ATR,
  DistMatrix<T,MC,STAR,Int>& ABL, DistMatrix<T,MC,STAR,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::View2x2");
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
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = ATL.ColShift();
        this->localMatrix_.View2x2
        ( ATL.LocalMatrix(), ATR.LocalMatrix(),
          ABL.LocalMatrix(), ABR.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::LockedView2x2
( const DistMatrix<T,MC,STAR,Int>& ATL, const DistMatrix<T,MC,STAR,Int>& ATR,
  const DistMatrix<T,MC,STAR,Int>& ABL, const DistMatrix<T,MC,STAR,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedView2x2");
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
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = ATL.ColShift();
        this->localMatrix_.LockedView2x2
        ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
          ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.ResizeTo
        ( LocalLength(height,this->ColShift(),this->Grid().Height()), width );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline T
DistMatrix<T,MC,STAR,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that 
    // row within each process column
    const elemental::Grid& g = this->Grid();
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();

    T u;
    if( g.VCRank() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        u = this->GetLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRow), g.ViewingComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        this->SetLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        this->UpdateLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetDiagonal
( DistMatrix<T,MC,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const Int diagLength = this->DiagonalLength(offset);
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
        const Int r = g.Height();
        const Int colShift = this->ColShift();
        const Int diagShift = d.ColShift();

        Int iStart, jStart;
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

        const Int iLocalStart = (iStart-colShift) / r;
        const Int localDiagLength = d.LocalHeight();
        T* dLocalBuffer = d.LocalBuffer();
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            dLocalBuffer[k] = thisLocalBuffer[iLocal+jLocal*thisLDim];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetDiagonal
( DistMatrix<T,STAR,MC,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const Int diagLength = this->DiagonalLength(offset);
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
        const Int r = g.Height();
        const Int colShift = this->ColShift();
        const Int diagShift = d.RowShift();

        Int iStart, jStart;
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

        const Int iLocalStart = (iStart-colShift) / r;
        const Int localDiagLength = d.LocalWidth();
        T* dLocalBuffer = d.LocalBuffer();
        const Int dLDim = d.LocalLDim();
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            dLocalBuffer[k*dLDim] = thisLocalBuffer[iLocal+jLocal*thisLDim];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetDiagonal
( const DistMatrix<T,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetDiagonal");
    this->AssertSameGrid( d );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const Int diagLength = this->DiagonalLength(offset);
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
        const Int r = g.Height();
        const Int colShift = this->ColShift();
        const Int diagShift = d.ColShift();

        Int iStart, jStart;
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

        const Int iLocalStart = (iStart-colShift)/r;
        const Int localDiagLength = d.LocalHeight();
        const T* dLocalBuffer = d.LockedLocalBuffer();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetDiagonal
( const DistMatrix<T,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetDiagonal");
    this->AssertSameGrid( d );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const Int diagLength = this->DiagonalLength(offset);
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
        const Int r = g.Height();
        const Int colShift = this->ColShift();
        const Int diagShift = d.RowShift();

        Int iStart, jStart;
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

        const Int iLocalStart = (iStart-colShift)/r;
        const Int localDiagLength = d.LocalWidth();
        const T* dLocalBuffer = d.LockedLocalBuffer();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int dLDim = d.LocalLDim();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
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

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::MakeTrapezoidal
( Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int r = this->Grid().Height();
    const Int colShift = this->ColShift();

    if( this->Grid().InGrid() )
    {
        if( uplo == LOWER )
        {
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                Int lastZeroRow = ( side==LEFT ? j-offset-1
                                               : j-offset+height-width-1 );
                if( lastZeroRow >= 0 )
                {
                    Int boundary = std::min( lastZeroRow+1, height );
                    Int numZeroRows = RawLocalLength( boundary, colShift, r );
                    T* thisCol = &thisLocalBuffer[j*thisLDim];
                    memset( thisCol, 0, numZeroRows*sizeof(T) );
                }
            }
        }
        else
        {
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                Int firstZeroRow = 
                    ( side==LEFT ? std::max(j-offset+1,0)
                      : std::max(j-offset+height-width+1,0) );
                Int numNonzeroRows = RawLocalLength(firstZeroRow,colShift,r);
                if( numNonzeroRows < localHeight )
                {
                    T* thisCol = &thisLocalBuffer[numNonzeroRows+j*thisLDim];
                    memset
                    ( thisCol, 0, (localHeight-numNonzeroRows)*sizeof(T) );
                }
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::ScaleTrapezoid
( T alpha, Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::ScaleTrapezoid");
    this->AssertNotLockedView();
#endif

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int r = this->Grid().Height();
    const Int colShift = this->ColShift();

    if( uplo == UPPER )
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            Int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            Int boundary = std::min( lastRow+1, height );
            Int numRows = RawLocalLength( boundary, colShift, r );
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            for( Int iLocal=0; iLocal<numRows; ++iLocal )
                thisCol[iLocal] *= alpha;
        }
    }
    else
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            Int firstRow = ( side==LEFT ? std::max(j-offset,0)
                                        : std::max(j+height-width-offset,0) );
            Int numZeroRows = RawLocalLength( firstRow, colShift, r );
            T* thisCol = &thisLocalBuffer[numZeroRows+j*thisLDim];
            for( Int iLocal=0; iLocal<(localHeight-numZeroRows); ++iLocal )
                thisCol[iLocal] *= alpha;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif

    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int r = this->Grid().Height();
    const Int colShift = this->ColShift();

    this->localMatrix_.SetToZero();

    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const Int i = colShift + iLocal*r;
        if( i < width )
            thisLocalBuffer[iLocal+i*thisLDim] = 1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif

    if( this->Grid().InGrid() )
    {
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int bufSize = localHeight*width;

        this->auxMemory_.Require( bufSize );

        // Create random matrix on process column 0, then broadcast
        T* buffer = this->auxMemory_.Buffer();
        if( this->grid_->MRRank() == 0 )
        {
            for( Int j=0; j<width; ++j )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    buffer[iLocal+j*localHeight] = SampleUnitBall<T>();
        }
        mpi::Broadcast( buffer, bufSize, 0, this->Grid().MRComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* bufferCol = &buffer[j*localHeight];
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            memcpy( thisCol, bufferCol, localHeight*sizeof(T) );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void 
DistMatrix<T,MC,STAR,Int>::SumOverRow()
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SumOverRow");
    this->AssertNotLockedView();
#endif

    if( this->Grid().InGrid() )
    {
        const Int localHeight = this->LocalHeight(); 
        const Int localWidth = this->LocalWidth();
        const Int localSize = 
            std::max( localHeight*localWidth, mpi::MIN_COLL_MSG );

        this->auxMemory_.Require( 2*localSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[localSize];

        // Pack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
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
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
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

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MC,MR]");
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
            if( g.InGrid() )
            {
                this->colShift_ = 
                    Shift( g.MCRank(), this->ColAlignment(), g.Height() );
            }
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( g.InGrid() )
    {
        if( this->ColAlignment() == A.ColAlignment() )
        {
            if( A.Width() == 1 )
            {
                if( g.MRRank() == A.RowAlignment() )
                    this->localMatrix_ = A.LockedLocalMatrix();

                // Communicate
                mpi::Broadcast
                ( this->localMatrix_.Buffer(), this->LocalHeight(),
                  A.RowAlignment(), g.MRComm() );
            }
            else
            {
                const Int c = g.Width();

                const Int width = this->Width();
                const Int localHeight = this->LocalHeight();
                const Int localWidthOfA = A.LocalWidth();
                const Int maxLocalWidth = MaxLocalLength(width,c);

                const Int portionSize = 
                    std::max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

                this->auxMemory_.Require( (c+1)*portionSize );

                T* buffer = this->auxMemory_.Buffer();
                T* originalData = &buffer[0];
                T* gatheredData = &buffer[portionSize];

                // Pack
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                {
                    const T* ACol = &ALocalBuffer[jLocal*ALDim];
                    T* originalDataCol = &originalData[jLocal*localHeight];
                    memcpy( originalDataCol, ACol, localHeight*sizeof(T) );
                }

                // Communicate
                mpi::AllGather
                ( originalData, portionSize,
                  gatheredData, portionSize, g.MRComm() );

                // Unpack
                const Int rowAlignmentOfA = A.RowAlignment();
                T* thisLocalBuffer = this->LocalBuffer();
                const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int k=0; k<c; ++k )
                {
                    const T* data = &gatheredData[k*portionSize];

                    const Int rowShift = RawShift( k, rowAlignmentOfA, c );
                    const Int localWidth = RawLocalLength( width, rowShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                    #pragma omp parallel for
#endif
                    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                    {
                        const T* dataCol = &data[jLocal*localHeight];
                        T* thisCol = &thisLocalBuffer[(rowShift+jLocal*c)*thisLDim];
                        memcpy( thisCol, dataCol, localHeight*sizeof(T) );
                    }
                }
                this->auxMemory_.Release();
            }
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,* ] <- [MC,MR]." << std::endl;
#endif
            const Int r = g.Height();
            const Int c = g.Width();
            const Int row = g.MCRank();

            const Int colAlignment = this->ColAlignment();
            const Int colAlignmentOfA = A.ColAlignment();
            const Int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
            const Int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

            if( A.Width()==1 )
            {
                const Int localHeight = this->LocalHeight();

                if( this->grid_->MRRank() == A.RowAlignment() )
                {
                    const Int localHeightOfA = A.LocalHeight();

                    this->auxMemory_.Require( localHeightOfA );
                    T* buffer = this->auxMemory_.Buffer();

                    // Pack
                    const T* ACol = A.LockedLocalBuffer(0,0);
                    memcpy( buffer, ACol, localHeightOfA*sizeof(T) );

                    // Communicate
                    mpi::SendRecv
                    ( buffer, localHeightOfA, sendRow, 0,
                      this->localMatrix_.Buffer(), localHeight, recvRow, 
                      mpi::ANY_TAG, g.MCComm() );

                    this->auxMemory_.Release();
                }

                // Communicate
                mpi::Broadcast
                ( this->localMatrix_.Buffer(), localHeight, A.RowAlignment(),
                  g.MRComm() );
            }
            else
            {
                const Int height = this->Height();
                const Int width = this->Width();
                const Int localHeight = this->LocalHeight();
                const Int localHeightOfA = A.LocalHeight();
                const Int localWidthOfA  = A.LocalWidth();
                const Int maxLocalHeight = MaxLocalLength(height,r);
                const Int maxLocalWidth  = MaxLocalLength(width,c);

                const Int portionSize = 
                    std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

                this->auxMemory_.Require( (c+1)*portionSize );

                T* buffer = this->auxMemory_.Buffer();
                T* firstBuffer = &buffer[0];
                T* secondBuffer = &buffer[portionSize];

                // Pack the currently owned local data of A into the second 
                // buffer
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                {
                    const T* ACol = &ALocalBuffer[jLocal*ALDim];
                    T* secondBufferCol = &secondBuffer[jLocal*localHeightOfA];
                    memcpy( secondBufferCol, ACol, localHeightOfA*sizeof(T) );
                }

                // Perform the SendRecv: puts the new data into the first buffer
                mpi::SendRecv
                ( secondBuffer, portionSize, sendRow, 0,
                  firstBuffer,  portionSize, recvRow, mpi::ANY_TAG, 
                  g.MCComm() );

                // Use the output of the SendRecv as the input to the AllGather
                mpi::AllGather
                ( firstBuffer,  portionSize, 
                  secondBuffer, portionSize, g.MRComm() );

                // Unpack the contents of each member of the process row
                const Int rowAlignmentOfA = A.RowAlignment();
                T* thisLocalBuffer = this->LocalBuffer();
                const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int k=0; k<c; ++k )
                {
                    const T* data = &secondBuffer[k*portionSize];

                    const Int rowShift = RawShift( k, rowAlignmentOfA, c ); 
                    const Int localWidth = RawLocalLength( width, rowShift, c );
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                    #pragma omp parallel for
#endif
                    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                    {
                        const T* dataCol = &data[jLocal*localHeight];
                        T* thisCol = &thisLocalBuffer[(rowShift+jLocal*c)*thisLDim];
                        memcpy( thisCol, dataCol, localHeight*sizeof(T) );
                    }
                }
                this->auxMemory_.Release();
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

// LEFT OFF HERE FOR IN/OUT OF GRID MODIFICATIONS

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MC,* ]");
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
            std::cerr << "Unaligned [MC,* ] <- [MC,* ]." << std::endl;
#endif
        const Int rank = g.MCRank();
        const Int r = g.Height();

        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int sendRank = (rank+r+colAlignment-colAlignmentOfA) % r;
        const Int recvRank = (rank+r+colAlignmentOfA-colAlignment) % r;

        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localHeightOfA = A.LocalHeight();

        const Int sendSize = localHeightOfA * width;
        const Int recvSize = localHeight * width;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ALocalBuffer[j*ALDim];
            T* sendBufferCol = &sendBuffer[j*localHeightOfA];
            memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.MCComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
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

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MC,MR,Int> A_MC_MR(true,false,this->ColAlignment(),0,g);

    A_MC_MR = A;
    *this   = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MC,* ] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MC,* ] = [* ,MD] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = A;

    std::auto_ptr< DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = A;

    std::auto_ptr< DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,MR,MC,Int> > 
        A_MR_MC( new DistMatrix<T,MR,MC,Int>(g) );
    *A_MR_MC = A;

    std::auto_ptr< DistMatrix<T,VR,STAR,Int> > 
        A_VR_STAR( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    std::auto_ptr< DistMatrix<T,VC,STAR,Int> > 
        A_VC_STAR( new DistMatrix<T,VC,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [VC,* ]");
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
          "The vector version of [MC,* ] <- [VC,* ] is not yet written, but"
          " it only requires a modification of the vector version of "
          "[* ,MR] <- [* ,VR]" << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.VCRank() == 0 )
    {
        std::cerr << 
          "[MC,* ] <- [VC,* ] potentially causes a large amount of cache-"
          "thrashing. If possible avoid it by performing the redistribution"
          " with a (conjugate)-transpose: \n"
          "  [* ,MC].(Conjugate)TransposeFrom([VC,* ])" << std::endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment() % g.Height();
            this->colShift_ = 
                Shift( g.MCRank(), this->ColAlignment(), g.Height() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() % g.Height() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int row = g.MCRank();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLocalLength(height,p);

        const Int portionSize = 
            std::max(maxLocalHeightOfA*width,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack 
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
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
        const Int colShift = this->ColShift();
        const Int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];    

            const Int colShiftOfA = RawShift( row+r*k, colAlignmentOfA, p );
            const Int colOffset = (colShiftOfA-colShift) / r;
            const Int localHeight = RawLocalLength( height, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int j=0; j<width; ++j )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colOffset+iLocal*c)+j*thisLDim] = 
                        data[iLocal+j*localHeight];
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [MC,* ] <- [VC,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int row = g.MCRank();
        const Int rank = g.VCRank();

        // Perform the SendRecv to make A have the same colAlignment
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int colShift = this->ColShift();

        const Int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const Int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLocalLength(height,p);

        const Int portionSize = 
            std::max(maxLocalHeightOfA*width,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ALocalBuffer[j*ALDim];
            T* secondBufferCol = &secondBuffer[j*localHeightOfA];
            memcpy( secondBufferCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, mpi::ANY_TAG, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MRComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const Int colShiftOfA = RawShift( row+r*k, colAlignment, p );
            const Int colOffset = (colShiftOfA-colShift) / r;
            const Int localHeight = RawLocalLength( height, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int j=0; j<width; ++j )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colOffset+iLocal*c)+j*thisLDim] = 
                        data[iLocal+j*localHeight];
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr< DistMatrix<T,STAR,VR,Int> > 
        A_STAR_VR( new DistMatrix<T,STAR,VR,Int>(g) );
    *A_STAR_VR = A;

    std::auto_ptr< DistMatrix<T,MC,MR,Int> > 
        A_MC_MR
        ( new DistMatrix<T,MC,MR,Int>(true,false,this->ColAlignment(),0,g) );
    *A_MC_MR = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR,Int> A_VC_STAR(true,this->ColAlignment(),g);

    A_VC_STAR = A;
    *this = A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MC,MR,Int> A_MC_MR(true,false,this->ColAlignment(),0,g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Int r = this->Grid().Height(); 
    const Int colShift = this->ColShift();

    const Int localHeight = this->LocalHeight();
    const Int width = this->Width();

    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( Int j=0; j<width; ++j )
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
            thisLocalBuffer[iLocal+j*thisLDim] = 
                ALocalBuffer[(colShift+iLocal*r)+j*ALDim];
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

} // namespace elemental
