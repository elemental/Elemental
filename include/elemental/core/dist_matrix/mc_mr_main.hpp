/*
   Copyright (c) 2009-2012, Jack Poulson
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
inline
DistMatrix<T,MC,MR,Int>::DistMatrix( const elemental::Grid& grid )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   (grid.InGrid() ? grid.MCRank() : 0),
   (grid.InGrid() ? grid.MRRank() : 0),
    0,0,grid)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width, const elemental::Grid& grid )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   (grid.InGrid() ? grid.MCRank() : 0),
   (grid.InGrid() ? grid.MRRank() : 0),
   (grid.InGrid() ? LocalLength(height,grid.MCRank(),0,grid.Height()) : 0),
   (grid.InGrid() ? LocalLength(width,grid.MRRank(),0,grid.Width()) : 0),
    grid)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( bool constrainedColAlignment, bool constrainedRowAlignment,
  Int colAlignment, Int rowAlignment, const elemental::Grid& grid )
: AbstractDistMatrix<T,Int>
  (0,0,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (grid.InGrid() ? Shift(grid.MCRank(),colAlignment,grid.Height()) : 0),
   (grid.InGrid() ? Shift(grid.MRRank(),rowAlignment,grid.Width()) : 0),
   0,0,grid)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  Int colAlignment, Int rowAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  Int colAlignment, Int rowAlignment, Int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment,
  const T* buffer, Int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment,
  T* buffer, Int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::DistMatrix");
#endif
    if( MC != U || MR != V ||
        reinterpret_cast<const DistMatrix<T,MC,MR,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [MC,MR] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetGrid( const elemental::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->colAlignment_ = 0;
    this->rowAlignment_ = 0;
    if( grid.InGrid() )
    {
        this->colShift_ = grid.MCRank();
        this->rowShift_ = grid.MRRank();
    }
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->rowAlignment_ = A.RowAlignment();
    this->constrainedColAlignment_ = true;
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.ColShift();
        this->rowShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,MC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.ColShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,STAR,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->rowAlignment_ = A.ColAlignment();
    this->constrainedColAlignment_ = true;
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.RowShift();
        this->rowShift_ = A.ColShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,MR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.ColShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,STAR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,VC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([VC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.ColAlignment() % g.Height();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->colShift_ =
            Shift<Int>( g.MCRank(), this->ColAlignment(), g.Height() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,STAR,VC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,VC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.RowAlignment() % g.Height();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->colShift_ =
            Shift( g.MCRank(), this->ColAlignment(), g.Height() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,VR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([VR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.ColAlignment() % g.Width();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->rowShift_ =
            Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,STAR,VR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,VR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.RowAlignment() % g.Width();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->rowShift_ =
            Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.ColShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,MC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([MC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.ColShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,STAR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([* ,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,VC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([VC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.ColAlignment() % g.Height();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->colShift_ =
            Shift( g.MCRank(), this->ColAlignment(), g.Height() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,STAR,VC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([* ,VC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.RowAlignment() % g.Height();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->colShift_ =
            Shift( g.MCRank(), this->ColAlignment(), g.Height() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,STAR,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([* ,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.ColShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,MR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([MR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.ColShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,VR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([VR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.ColAlignment() % g.Width();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->rowShift_ =
            Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,STAR,VR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([* ,VR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.RowAlignment() % g.Width();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->rowShift_ =
            Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();

    const Int r = g.Height();
    const Int c = g.Width();

    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;

    const Int height = this->Height();
    const Int width  = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localWidth  = this->LocalWidth();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

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
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
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
            for( Int i=0; i<height; ++i )
            {
                for( Int j=0; j<width; ++j )
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

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::Align( Int colAlignment, Int rowAlignment )
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
    this->colAlignment_ = colAlignment;
    this->rowAlignment_ = rowAlignment;
    this->constrainedColAlignment_ = true;
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->localMatrix_.ResizeTo( 0, 0 );
        this->colShift_ = Shift( g.MCRank(), colAlignment, g.Height() );
        this->rowShift_ = Shift( g.MRRank(), rowAlignment, g.Width() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::AlignCols( Int colAlignment )
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
    this->colAlignment_ = colAlignment;
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->localMatrix_.ResizeTo( 0, 0 );
        this->colShift_ = Shift( g.MCRank(), colAlignment, g.Height() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::AlignRows( Int rowAlignment )
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
    this->rowAlignment_ = rowAlignment;
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->localMatrix_.ResizeTo( 0, 0 );
        this->rowShift_ = Shift( g.MRRank(), rowAlignment, g.Width() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::View( DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->rowAlignment_ = A.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.ColShift();
        this->rowShift_ = A.RowShift();
        this->localMatrix_.View( A.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::View
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  T* buffer, Int ldim, const elemental::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->rowAlignment_ = rowAlignment;
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->grid_->InGrid() )
    {
        this->colShift_ = Shift(g.MCRank(),colAlignment,g.Height());
        this->rowShift_ = Shift(g.MRRank(),rowAlignment,g.Width());
        Int localHeight = LocalLength(height,this->colShift_,g.Height());
        Int localWidth = LocalLength(width,this->rowShift_,g.Width());
        this->localMatrix_.View( localHeight, localWidth, buffer, ldim );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::LockedView( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_  = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->rowAlignment_ = A.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.ColShift();
        this->rowShift_ = A.RowShift();
        this->localMatrix_.LockedView( A.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::LockedView
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  const T* buffer, Int ldim, const elemental::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->rowAlignment_ = rowAlignment;
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->grid_->InGrid() )
    {
        this->colShift_ = Shift(g.MCRank(),colAlignment,g.Height());
        this->rowShift_ = Shift(g.MRRank(),rowAlignment,g.Width());
        Int localHeight = LocalLength(height,this->colShift_,g.Height());
        Int localWidth = LocalLength(width,this->rowShift_,g.Width());
        this->localMatrix_.LockedView( localHeight, localWidth, buffer, ldim );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::View
( DistMatrix<T,MC,MR,Int>& A, 
  Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_  = width;

    const elemental::Grid& g = this->Grid();
    const Int r   = g.Height();
    const Int c   = g.Width();
    const Int row = g.MCRank();
    const Int col = g.MRRank();

    this->colAlignment_ = (A.ColAlignment()+i) % r;
    this->rowAlignment_ = (A.RowAlignment()+j) % c;
    this->viewing_ = true;
    this->lockedView_ = false;
  
    if( g.InGrid() )
    {
        this->colShift_ = Shift( row, this->ColAlignment(), r );
        this->rowShift_ = Shift( col, this->RowAlignment(), c );

        const Int localHeightBehind = LocalLength(i,A.ColShift(),r);
        const Int localWidthBehind  = LocalLength(j,A.RowShift(),c);

        const Int localHeight = LocalLength( height, this->ColShift(), r );
        const Int localWidth  = LocalLength( width,  this->RowShift(), c );

        this->localMatrix_.View
        ( A.LocalMatrix(), localHeightBehind, localWidthBehind,
                           localHeight,       localWidth );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::LockedView
( const DistMatrix<T,MC,MR,Int>& A, 
  Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_  = width;

    const elemental::Grid& g = this->Grid();
    const Int r   = g.Height();
    const Int c   = g.Width();
    const Int row = g.MCRank();
    const Int col = g.MRRank();

    this->colAlignment_ = (A.ColAlignment()+i) % r;
    this->rowAlignment_ = (A.RowAlignment()+j) % c;
    this->viewing_ = true;
    this->lockedView_ = true;
  
    if( g.InGrid() )
    {
        this->colShift_ = Shift( row, this->ColAlignment(), r );
        this->rowShift_ = Shift( col, this->RowAlignment(), c );

        const Int localHeightBehind = LocalLength(i,A.ColShift(),r);
        const Int localWidthBehind  = LocalLength(j,A.RowShift(),c);

        const Int localHeight = LocalLength( height, this->ColShift(), r );
        const Int localWidth  = LocalLength( width,  this->RowShift(), c );

        this->localMatrix_.LockedView
        ( A.LockedLocalMatrix(), localHeightBehind, localWidthBehind,
                                 localHeight,       localWidth );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::View1x2
( DistMatrix<T,MC,MR,Int>& AL, DistMatrix<T,MC,MR,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View1x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->rowAlignment_ = AL.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = AL.ColShift();
        this->rowShift_ = AL.RowShift();
        this->localMatrix_.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::LockedView1x2
( const DistMatrix<T,MC,MR,Int>& AL, const DistMatrix<T,MC,MR,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView1x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->rowAlignment_ = AL.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = AL.ColShift();
        this->rowShift_ = AL.RowShift();
        this->localMatrix_.LockedView1x2
        ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::View2x1
( DistMatrix<T,MC,MR,Int>& AT,
  DistMatrix<T,MC,MR,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View2x1");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->rowAlignment_ = AT.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = AT.ColShift();
        this->rowShift_ = AT.RowShift();
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
DistMatrix<T,MC,MR,Int>::LockedView2x1
( const DistMatrix<T,MC,MR,Int>& AT,
  const DistMatrix<T,MC,MR,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView2x1");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->rowAlignment_ = AT.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = AT.ColShift();
        this->rowShift_ = AT.RowShift();
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
DistMatrix<T,MC,MR,Int>::View2x2
( DistMatrix<T,MC,MR,Int>& ATL, DistMatrix<T,MC,MR,Int>& ATR,
  DistMatrix<T,MC,MR,Int>& ABL, DistMatrix<T,MC,MR,Int>& ABR )
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
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->colAlignment_ = ATL.ColAlignment();
    this->rowAlignment_ = ATL.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = ATL.ColShift();
        this->rowShift_ = ATL.RowShift();
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
DistMatrix<T,MC,MR,Int>::LockedView2x2
( const DistMatrix<T,MC,MR,Int>& ATL, 
  const DistMatrix<T,MC,MR,Int>& ATR,
  const DistMatrix<T,MC,MR,Int>& ABL, 
  const DistMatrix<T,MC,MR,Int>& ABR )
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
    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->colAlignment_ = ATL.ColAlignment();
    this->rowAlignment_ = ATL.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = ATL.ColShift();
        this->rowShift_ = ATL.RowShift();
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
DistMatrix<T,MC,MR,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::ResizeTo");
    this->AssertNotLockedView(); // this should be relaxed...
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.ResizeTo
        ( LocalLength(height,this->ColShift(),this->Grid().Height()),
          LocalLength(width, this->RowShift(),this->Grid().Width()) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline T
DistMatrix<T,MC,MR,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const elemental::Grid& g = this->Grid();
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol * g.Height();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLocal = (i-this->ColShift()) / g.Height();
        const Int jLocal = (j-this->RowShift()) / g.Width();
        u = this->GetLocalEntry(iLocal,jLocal);
    }
    mpi::Broadcast
    ( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const Int iLocal = (i-this->ColShift()) / g.Height();
        const Int jLocal = (j-this->RowShift()) / g.Width();
        this->SetLocalEntry(iLocal,jLocal,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();
    const Int ownerCol = (j + this->RowAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const Int iLocal = (i-this->ColShift()) / g.Height();
        const Int jLocal = (j-this->RowShift()) / g.Width();
        this->UpdateLocalEntry(iLocal,jLocal,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::GetDiagonal
( DistMatrix<T,MD,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const Int diagLength = this->DiagonalLength(offset);
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
        const Int r = g.Height();
        const Int c = g.Width();
        const Int lcm = g.LCM();
        const Int colShift = this->ColShift();
        const Int rowShift = this->RowShift();
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
        const Int jLocalStart = (jStart-rowShift) / c;

        const Int localDiagLength = d.LocalHeight();
        T* dLocalBuffer = d.LocalBuffer();
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k*(lcm/r);
            const Int jLocal = jLocalStart+k*(lcm/c);
            dLocalBuffer[k] = thisLocalBuffer[iLocal+jLocal*thisLDim];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::GetDiagonal
( DistMatrix<T,STAR,MD,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const Int diagLength = this->DiagonalLength(offset);
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
        const Int r = g.Height();
        const Int c = g.Width();
        const Int lcm = g.LCM();
        const Int colShift = this->ColShift();
        const Int rowShift = this->RowShift();
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
        const Int jLocalStart = (jStart-rowShift) / c;

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
            const Int iLocal = iLocalStart+k*(lcm/r);
            const Int jLocal = jLocalStart+k*(lcm/c);
            dLocalBuffer[k*dLDim] = thisLocalBuffer[iLocal+jLocal*thisLDim];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetDiagonal
( const DistMatrix<T,MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetDiagonal");
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
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const Int r = g.Height();
        const Int c = g.Width();
        const Int lcm = g.LCM();
        const Int colShift = this->ColShift();
        const Int rowShift = this->RowShift();
        const Int diagShift = d.ColShift();

        Int iStart,jStart;
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
        const Int jLocalStart = (jStart-rowShift) / c;

        const Int localDiagLength = d.LocalHeight();
        const T* dLocalBuffer = d.LockedLocalBuffer();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k*(lcm/r);
            const Int jLocal = jLocalStart+k*(lcm/c);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetDiagonal
( const DistMatrix<T,STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetDiagonal");
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
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const Int r = g.Height();
        const Int c = g.Width();
        const Int lcm = g.LCM();
        const Int colShift = this->ColShift();
        const Int rowShift = this->RowShift();
        const Int diagShift = d.RowShift();

        Int iStart,jStart;
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
        const Int jLocalStart = (jStart-rowShift) / c;

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
            const Int iLocal = iLocalStart+k*(lcm/r);
            const Int jLocal = jLocalStart+k*(lcm/c);
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
DistMatrix<T,MC,MR,Int>::MakeTrapezoidal
( Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int r = this->Grid().Height();
    const Int c = this->Grid().Width();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();

    if( uplo == LOWER )
    {

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const Int j = rowShift + jLocal*c;
            const Int lastZeroRow = 
                ( side==LEFT ? j-offset-1
                             : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                const Int boundary = std::min( lastZeroRow+1, height );
                const Int numZeroRows = RawLocalLength( boundary, colShift, r );
                std::memset
                ( &thisLocalBuffer[jLocal*thisLDim], 0, numZeroRows*sizeof(T) );
            }
        }
    }
    else
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const Int j = rowShift + jLocal*c;
            const Int firstZeroRow = 
                ( side==LEFT ? std::max(j-offset+1,0)
                             : std::max(j-offset+height-width+1,0) );
            const Int numNonzeroRows = RawLocalLength(firstZeroRow,colShift,r);
            if( numNonzeroRows < localHeight )
            {
                T* thisCol = &thisLocalBuffer[numNonzeroRows+jLocal*thisLDim];
                std::memset
                ( thisCol, 0, (localHeight-numNonzeroRows)*sizeof(T) );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::ScaleTrapezoid
( T alpha, Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::ScaleTrapezoid");
    this->AssertNotLockedView();
#endif
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int r = this->Grid().Height();
    const Int c = this->Grid().Width();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

    if( uplo == UPPER )
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            Int j = rowShift + jLocal*c;
            Int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            Int boundary = std::min( lastRow+1, height );
            Int numRows = RawLocalLength( boundary, colShift, r );
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
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
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            Int j = rowShift + jLocal*c;
            Int firstRow = 
                ( side==LEFT ? std::max(j-offset,0) 
                             : std::max(j-offset+height-width,0) );
            Int numZeroRows = RawLocalLength( firstRow, colShift, r );
            T* thisCol = &thisLocalBuffer[numZeroRows+jLocal*thisLDim];
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
DistMatrix<T,MC,MR,Int>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int r = this->Grid().Height();
    const Int c = this->Grid().Width();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

    this->localMatrix_.SetToZero();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const Int i = colShift + iLocal*r;                
        if( i % c == rowShift )
        {
            const Int jLocal = (i-rowShift) / c;
            if( jLocal < localWidth )
                thisLocalBuffer[iLocal+jLocal*thisLDim] = 1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            this->SetLocalEntry(iLocal,jLocal,SampleUnitBall<T>());
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::AdjointFrom
( const DistMatrix<T,STAR,MC,Int>& A )
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
            this->colAlignment_ = A.RowAlignment();
            if( g.InGrid() )
            {
                this->colShift_ = 
                    Shift( g.MCRank(), this->ColAlignment(), g.Height() );
            }
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( g.InGrid() ) 
    { 
        if( this->ColAlignment() == A.RowAlignment() )
        {
            const Int c = g.Width();
            const Int rowShift = this->RowShift();

            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();

            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[iLocal+jLocal*thisLDim] = 
                        Conj( ALocalBuffer[(rowShift+jLocal*c)+iLocal*ALDim] );
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR]::AdjointFrom." << std::endl;
#endif
            const Int r = g.Height();
            const Int c = g.Width();
            const Int rank = g.MCRank();
            const Int rowShift = this->RowShift();
            const Int colAlignment = this->ColAlignment();
            const Int rowAlignmentOfA = A.RowAlignment();

            const Int sendRank = (rank+r+colAlignment-rowAlignmentOfA) % r;
            const Int recvRank = (rank+r+rowAlignmentOfA-colAlignment) % r;

            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localWidthOfA = A.LocalWidth();

            const Int sendSize = localWidthOfA * localWidth;
            const Int recvSize = localHeight * localWidth;

            this->auxMemory_.Require( sendSize + recvSize );

            T* buffer = this->auxMemory_.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];

            // Pack
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                for( Int iLocal=0; iLocal<localWidthOfA; ++iLocal )
                    sendBuffer[iLocal+jLocal*localWidth] = 
                        Conj( ALocalBuffer[(rowShift+jLocal*c)+iLocal*ALDim] );

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
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                std::memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->auxMemory_.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::AdjointFrom
( const DistMatrix<T,MR,STAR,Int>& A )
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
            this->rowAlignment_ = A.ColAlignment();
            if( g.InGrid() )
                this->rowShift_ = 
                    Shift( g.MRRank(), this->RowAlignment(), g.Height() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }
    if( this->rowAlignment_ != A.ColAlignment() )
        throw std::logic_error("Unaligned AdjointFrom");

    if( g.InGrid() ) 
    { 
        const Int r = g.Height();
        const Int colShift = this->ColShift();

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();

        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                thisLocalBuffer[iLocal+jLocal*thisLDim] = 
                    Conj(ALocalBuffer[jLocal+(colShift+iLocal*r)*ALDim]);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::TransposeFrom
( const DistMatrix<T,STAR,MC,Int>& A )
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
            this->colAlignment_ = A.RowAlignment();
            if( g.InGrid() )
            {
                this->colShift_ = 
                    Shift( g.MCRank(), this->ColAlignment(), g.Height() );
            }
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( g.InGrid() ) 
    { 
        if( this->ColAlignment() == A.RowAlignment() )
        {
            const Int c = g.Width();
            const Int rowShift = this->RowShift();

            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();

            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[iLocal+jLocal*thisLDim] = 
                        ALocalBuffer[(rowShift+jLocal*c)+iLocal*ALDim];
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR]::TransposeFrom." << std::endl;
#endif
            const Int r = g.Height();
            const Int c = g.Width();
            const Int rank = g.MCRank();
            const Int rowShift = this->RowShift();
            const Int colAlignment = this->ColAlignment();
            const Int rowAlignmentOfA = A.RowAlignment();

            const Int sendRank = (rank+r+colAlignment-rowAlignmentOfA) % r;
            const Int recvRank = (rank+r+rowAlignmentOfA-colAlignment) % r;

            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localWidthOfA = A.LocalWidth();

            const Int sendSize = localWidthOfA * localWidth;
            const Int recvSize = localHeight * localWidth;

            this->auxMemory_.Require( sendSize + recvSize );

            T* buffer = this->auxMemory_.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];

            // Pack
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                for( Int iLocal=0; iLocal<localWidthOfA; ++iLocal )
                    sendBuffer[iLocal+jLocal*localWidth] = 
                        ALocalBuffer[(rowShift+jLocal*c)+iLocal*ALDim];

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
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                std::memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->auxMemory_.Release();
        }
    } 
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::TransposeFrom
( const DistMatrix<T,MR,STAR,Int>& A )
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
            this->rowAlignment_ = A.ColAlignment();
            if( g.InGrid() )
            {
                this->rowShift_ = 
                    Shift( g.MRRank(), this->RowAlignment(), g.Height() );
            }
        }
        this->ResizeTo( A.Width(), A.Height() );
    }
    if( this->rowAlignment_ != A.ColAlignment() )
        throw std::logic_error("Unaligned TransposeFrom");

    if( g.InGrid() ) 
    { 
        const Int r = g.Height();
        const Int colShift = this->ColShift();

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();

        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                thisLocalBuffer[iLocal+jLocal*thisLDim] = 
                    ALocalBuffer[jLocal+(colShift+iLocal*r)*ALDim];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
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
                this->colAlignment_ = A.ColAlignment();
                if( this->Grid().InGrid() )
                    this->colShift_ = A.ColShift();
            }
            if( !this->ConstrainedRowAlignment() )
            {
                this->rowAlignment_ = A.RowAlignment();
                if( this->Grid().InGrid() )
                    this->rowShift_ = A.RowShift();
            }
            this->ResizeTo( A.Height(), A.Width() );
        }

        if( this->Grid().InGrid() ) 
        {
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
                    std::cerr << "Unaligned [MC,MR] <- [MC,MR]." << std::endl;
#endif
                const Int r = g.Height();
                const Int c = g.Width();
                const Int row = g.MCRank();
                const Int col = g.MRRank();

                const Int colAlignment = this->ColAlignment();
                const Int rowAlignment = this->RowAlignment();
                const Int colAlignmentOfA = A.ColAlignment();
                const Int rowAlignmentOfA = A.RowAlignment();

                const Int sendRow = (row+r+colAlignment-colAlignmentOfA)%r;
                const Int sendCol = (col+c+rowAlignment-rowAlignmentOfA)%c;
                const Int recvRow = (row+r+colAlignmentOfA-colAlignment)%r;
                const Int recvCol = (col+c+rowAlignmentOfA-rowAlignment)%c;
                const Int sendRank = sendRow + sendCol*r;
                const Int recvRank = recvRow + recvCol*r;

                const Int localHeight = this->LocalHeight();
                const Int localWidth = this->LocalWidth();
                const Int localHeightOfA = A.LocalHeight();
                const Int localWidthOfA = A.LocalWidth();

                const Int sendSize = localHeightOfA * localWidthOfA;
                const Int recvSize = localHeight * localWidth;

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
                for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                {
                    const T* ACol = &ALocalBuffer[jLocal*ALDim];
                    T* sendBufferCol = &sendBuffer[jLocal*localHeightOfA];
                    std::memcpy
                    ( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
                }

                // Communicate
                mpi::SendRecv
                ( sendBuffer, sendSize, sendRank, 0,
                  recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.VCComm() );

                // Unpack
                T* thisLocalBuffer = this->LocalBuffer();
                const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                    T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                    std::memcpy
                    ( thisCol, recvBufferCol, localHeight*sizeof(T) );
                }
                this->auxMemory_.Release();
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
        const Int r0 = this->Grid().Height();
        const Int c0 = this->Grid().Width();
        const Int rA = A.Grid().Height();
        const Int cA = A.Grid().Width();
        const Int myRow0 = this->Grid().MCRank();
        const Int myCol0 = this->Grid().MRRank();
        const Int myRowA = A.Grid().MCRank();
        const Int myColA = A.Grid().MRRank();
        const Int rowGCD = GCD( r0, rA );
        const Int colGCD = GCD( c0, cA );
        const Int rowLCM = r0*rA / rowGCD;
        const Int colLCM = c0*cA / colGCD;
        const Int numRowSends = r0 / rowGCD;
        const Int numColSends = c0 / colGCD;
        const Int localColStride0 = rowLCM / r0;
        const Int localRowStride0 = colLCM / c0;
        const Int localColStrideA = numRowSends;
        const Int localRowStrideA = numColSends;

        const Int colAlign0 = this->ColAlignment();
        const Int rowAlign0 = this->RowAlignment();
        const Int colAlignA = A.ColAlignment();
        const Int rowAlignA = A.RowAlignment();

        const bool inThisGrid = this->Grid().InGrid();
        const bool inAGrid = A.Grid().InGrid();

        const Int maxSendSize = 
            (A.Height()/(rA*localColStrideA)+1) * 
            (A.Width()/(cA*localRowStrideA)+1);

        // Have each member of A's grid individually send to all numRow x numCol
        // processes in order, while the members of this grid receive from all 
        // necessary processes at each step.
        Int requiredMemory = 0;
        if( inAGrid )
            requiredMemory += maxSendSize;
        if( inThisGrid )
            requiredMemory += maxSendSize;
        this->auxMemory_.Require( requiredMemory );
        T* buffer = this->auxMemory_.Buffer();
        Int offset = 0;
        T* sendBuffer = &buffer[offset];
        if( inAGrid )
            offset += maxSendSize;
        T* recvBuffer = &buffer[offset];

        Int recvRow = 0; // avoid compiler warnings...
        if( inAGrid )
            recvRow = (((myRowA+rA-colAlignA)%rA)+colAlign0)%r0;
        for( Int rowSendCycle=0; rowSendCycle<numRowSends; ++rowSendCycle )
        {
            Int recvCol = 0; // avoid compiler warnings...
            if( inAGrid )
                recvCol = (((myColA+cA-rowAlignA)%cA)+rowAlign0)%c0;

            for( Int colSendCycle=0; colSendCycle<numColSends; 
                 ++colSendCycle )
            {
                mpi::Request sendRequest;

                // Fire off this round of non-blocking sends
                if( inAGrid )
                {
                    // Pack the data
                    Int sendHeight = LocalLength
                        ( A.LocalHeight(), rowSendCycle, numRowSends );
                    Int sendWidth = LocalLength
                        ( A.LocalWidth(), colSendCycle, numColSends );
                    const T* ALocalBuffer = A.LockedLocalBuffer();
                    const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
                    #pragma omp parallel for
#endif
                    for( Int jLocal=0; jLocal<sendWidth; ++jLocal )
                    {
                        const Int j = colSendCycle+jLocal*localRowStrideA;
                        for( Int iLocal=0; iLocal<sendHeight; ++iLocal )
                        {
                            const Int i = rowSendCycle+iLocal*localColStrideA;
                            sendBuffer[iLocal+jLocal*sendHeight] = 
                                ALocalBuffer[i+j*ALDim];
                        }
                    }
                    // Send data
                    const Int recvVCRank = recvRow + recvCol*r0;
                    const Int recvViewingRank = 
                        this->Grid().VCToViewingMap( recvVCRank );
                    mpi::ISend
                    ( sendBuffer, sendHeight*sendWidth, recvViewingRank,
                      0, this->Grid().ViewingComm(), sendRequest );
                }
                // Perform this round of recv's
                if( inThisGrid )
                {
                    const Int sendRowOffset = (rowSendCycle*rA+colAlignA) % rA;
                    const Int sendColOffset = (colSendCycle*cA+rowAlignA) % cA;
                    const Int recvRowOffset = (rowSendCycle*rA+colAlign0) % r0;
                    const Int recvColOffset = (colSendCycle*cA+rowAlign0) % c0;
                    const Int firstSendRow = 
                        (((myRow0+r0-recvRowOffset)%r0)+sendRowOffset)%rA;
                    const Int firstSendCol = 
                        (((myCol0+c0-recvColOffset)%c0)+sendColOffset)%cA;

                    const Int rowShift = (myRow0+r0-recvRowOffset)%r0;
                    const Int colShift = (myCol0+c0-recvColOffset)%c0;
                    const Int numRowRecvs = LocalLength( rA, rowShift, r0 ); 
                    const Int numColRecvs = LocalLength( cA, colShift, c0 );

                    // Recv data
                    // For now, simply receive sequentially. Until we switch to 
                    // nonblocking recv's, we won't be using much of the 
                    // recvBuffer
                    Int sendRow = firstSendRow;
                    for( Int rowRecvCycle=0; rowRecvCycle<numRowRecvs; 
                         ++rowRecvCycle )
                    {
                        const Int sendRowShift = 
                            Shift( sendRow, colAlignA, rA ) + rowSendCycle*rA;
                        const Int sendHeight = 
                            LocalLength( A.Height(), sendRowShift, rowLCM );
                        const Int localColOffset = 
                            (sendRowShift-this->ColShift()) / r0;

                        Int sendCol = firstSendCol;
                        for( Int colRecvCycle=0; 
                             colRecvCycle<numColRecvs;  ++colRecvCycle )
                        {
                            const Int sendColShift = 
                                Shift( sendCol, rowAlignA, cA ) + 
                                colSendCycle*cA;
                            const Int sendWidth = 
                                LocalLength( A.Width(), sendColShift, colLCM );
                            const Int localRowOffset = 
                                (sendColShift-this->RowShift()) / c0;

                            const Int sendVCRank = sendRow+sendCol*rA;
                            const Int sendViewingRank = 
                                A.Grid().VCToViewingMap( sendVCRank );

                            mpi::Recv
                            ( recvBuffer, sendHeight*sendWidth, sendViewingRank,
                              0, this->Grid().ViewingComm() );
                            
                            // Unpack the data
                            T* thisLocalBuffer = this->LocalBuffer();
                            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
                            #pragma omp parallel for
#endif
                            for( Int jLocal=0; jLocal<sendWidth; ++jLocal )
                            {
                                const Int j =
                                    localRowOffset+jLocal*localRowStride0;
                                for( Int iLocal=0; iLocal<sendHeight; ++iLocal )
                                {
                                    const Int i = 
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
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
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
            const Int c = g.Width();
            const Int rowShift = this->RowShift();

            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();

            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(rowShift+jLocal*c)*ALDim];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                std::memcpy( thisCol, ACol, localHeight*sizeof(T) );
            }
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR] <- [MC,* ]." << std::endl;
#endif
            const Int r = g.Height();
            const Int c = g.Width();
            const Int rank = g.MCRank();
            const Int rowShift = this->RowShift();
            const Int colAlignment = this->ColAlignment();
            const Int colAlignmentOfA = A.ColAlignment();

            const Int sendRank = (rank+r+colAlignment-colAlignmentOfA) % r;
            const Int recvRank = (rank+r+colAlignmentOfA-colAlignment) % r;

            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localHeightOfA = A.LocalHeight();

            const Int sendSize = localHeightOfA * localWidth;
            const Int recvSize = localHeight * localWidth;

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
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(rowShift+jLocal*c)*ALDim];
                T* sendBufferCol = &sendBuffer[jLocal*localHeightOfA];
                std::memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
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
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                std::memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->auxMemory_.Release();
        }
    } 
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
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
            this->rowAlignment_ = A.RowAlignment();
            if( g.InGrid() )
            {
                this->rowShift_ = 
                    Shift( g.MRRank(), this->RowAlignment(), g.Width() );
            }
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( g.InGrid() ) 
    {
        if( this->RowAlignment() == A.RowAlignment() )
        {
            const Int r = g.Height();
            const Int colShift = this->ColShift();

            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();

            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[iLocal+jLocal*thisLDim] = 
                        ALocalBuffer[(colShift+iLocal*r)+jLocal*ALDim];
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR] <- [* ,MR]." << std::endl;
#endif
            const Int r = g.Height(); 
            const Int c = g.Width();
            const Int col = g.MRRank();
            const Int colShift = this->ColShift();
            const Int rowAlignment = this->RowAlignment();
            const Int rowAlignmentOfA = A.RowAlignment();

            const Int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
            const Int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localWidthOfA = A.LocalWidth();

            const Int sendSize = localHeight * localWidthOfA;
            const Int recvSize = localHeight * localWidth;

            this->auxMemory_.Require( sendSize + recvSize );

            T* buffer = this->auxMemory_.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];

            // Pack
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    sendBuffer[iLocal+jLocal*localHeight] = 
                        ALocalBuffer[(colShift+iLocal*r)+jLocal*ALDim];

            // Communicate
            mpi::SendRecv
            ( sendBuffer, sendSize, sendCol, 0,
              recvBuffer, recvSize, recvCol, mpi::ANY_TAG, g.MRComm() );

            // Unpack
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for  
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                std::memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->auxMemory_.Release();
        }
    } 
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
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

template<typename T,typename Int>
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
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

template<typename T,typename Int>
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
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

            const Int r = g.Height();
            const Int c = g.Width();
            const Int p = g.Size();
            const Int myRow = g.MCRank();
            const Int myCol = g.MRRank();
            const Int rankCM = g.VCRank();
            const Int rankRM = g.VRRank();
            const Int ownerCol = this->RowAlignment();
            const Int ownerRow = A.RowAlignment();
            const Int colAlignment = this->ColAlignment();
            const Int colAlignmentOfA = A.ColAlignment();
            const Int colShift = this->ColShift();
            const Int colShiftOfA = A.ColShift();

            const Int height = A.Height();
            const Int maxLocalHeight = MaxLocalLength(height,p);

            const Int portionSize = std::max(maxLocalHeight,mpi::MIN_COLL_MSG);

            const Int colShiftVC = Shift(rankCM,colAlignment,p);
            const Int colShiftVROfA = Shift(rankRM,colAlignmentOfA,p);
            const Int sendRankCM = (rankCM+(p+colShiftVROfA-colShiftVC)) % p;
            const Int recvRankRM = (rankRM+(p+colShiftVC-colShiftVROfA)) % p;
            const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

            this->auxMemory_.Require( (r+c)*portionSize );
            T* buffer = this->auxMemory_.Buffer();
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[c*portionSize];

            if( myRow == ownerRow )
            {
                // Pack
                const T* ALocalBuffer = A.LockedLocalBuffer();
#ifdef _OPENMP
                #pragma omp parallel for  
#endif
                for( Int k=0; k<r; ++k )
                {
                    T* data = &recvBuf[k*portionSize];

                    const Int shift = RawShift(myCol+c*k,colAlignmentOfA,p);
                    const Int offset = (shift-colShiftOfA) / c;
                    const Int thisLocalHeight = RawLocalLength(height,shift,p);

                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
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
                for( Int k=0; k<c; ++k )
                {
                    const T* data = &sendBuf[k*portionSize];

                    const Int shift = RawShift(myRow+r*k,colAlignment,p);
                    const Int offset = (shift-colShift) / r;
                    const Int thisLocalHeight = RawLocalLength(height,shift,p);

                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        thisLocalBuffer[offset+iLocal*c] = data[iLocal];
                }
            }
            this->auxMemory_.Release();
        }
        else if( A.Height() == 1 )
        {
            if( !this->Viewing() )
                this->ResizeTo( 1, A.Width() );

            const Int r = g.Height();
            const Int c = g.Width();
            const Int p = g.Size();
            const Int myRow = g.MCRank();
            const Int myCol = g.MRRank();
            const Int rankCM = g.VCRank();
            const Int rankRM = g.VRRank();
            const Int ownerRow = this->ColAlignment();
            const Int ownerCol = A.ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int rowAlignmentOfA = A.RowAlignment();
            const Int rowShift = this->RowShift();
            const Int rowShiftOfA = A.RowShift();

            const Int width = A.Width();
            const Int maxLocalWidth = MaxLocalLength(width,p);

            const Int portionSize = std::max(maxLocalWidth,mpi::MIN_COLL_MSG);

            const Int rowShiftVR = Shift(rankRM,rowAlignment,p);
            const Int rowShiftVCOfA = Shift(rankCM,rowAlignmentOfA,p);
            const Int sendRankRM = (rankRM+(p+rowShiftVCOfA-rowShiftVR)) % p;
            const Int recvRankCM = (rankCM+(p+rowShiftVR-rowShiftVCOfA)) % p;
            const Int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

            this->auxMemory_.Require( (r+c)*portionSize );
            T* buffer = this->auxMemory_.Buffer();
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[r*portionSize];

            if( myCol == ownerCol )
            {
                // Pack
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for  
#endif
                for( Int k=0; k<c; ++k )
                {
                    T* data = &recvBuf[k*portionSize];

                    const Int shift = RawShift(myRow+r*k,rowAlignmentOfA,p);
                    const Int offset = (shift-rowShiftOfA) / r;
                    const Int thisLocalWidth = RawLocalLength(width,shift,p);

                    for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
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
                const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for  
#endif
                for( Int k=0; k<r; ++k )
                {
                    const T* data = &sendBuf[k*portionSize];

                    const Int shift = RawShift(myCol+c*k,rowAlignment,p);
                    const Int offset = (shift-rowShift) / c;
                    const Int thisLocalWidth = RawLocalLength(width,shift,p);

                    for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                        thisLocalBuffer[(offset+jLocal*r)*thisLDim] = 
                            data[jLocal];
                }
            }

            this->auxMemory_.Release();
        }
        else
        {
            if( A.Height() >= A.Width() )
            {
                std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
                ( new DistMatrix<T,VR,STAR,Int>(g) );

                *A_VR_STAR = A;

                std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
                ( new DistMatrix<T,VC,STAR,Int>(true,this->ColAlignment(),g) );
                *A_VC_STAR = *A_VR_STAR;
                delete A_VR_STAR.release(); // lowers memory highwater

                *this = *A_VC_STAR;
            }
            else
            {
                std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
                ( new DistMatrix<T,STAR,VC,Int>(g) );
                *A_STAR_VC = A;

                std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
                ( new DistMatrix<T,STAR,VR,Int>(true,this->RowAlignment(),g) );
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

template<typename T,typename Int>
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();

    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = A;

    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
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
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();

    std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC,Int>(g) );
    *A_STAR_VC = A;

    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR,Int>(true,this->RowAlignment(),g) );
    *A_STAR_VR = *A_STAR_VC;
    delete A_STAR_VC.release(); // lowers memory highwater

    *this = *A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
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
            this->colAlignment_ = A.ColAlignment() % g.Height();
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
        if( this->ColAlignment() == A.ColAlignment() % g.Height() )
        {
            const Int r = g.Height();
            const Int c = g.Width();
            const Int p = r * c;
            const Int row = g.MCRank();
            const Int colShift = this->ColShift();
            const Int rowAlignment = this->RowAlignment();
            const Int colAlignmentOfA = A.ColAlignment();

            const Int height = this->Height();
            const Int width = this->Width();
            const Int localWidth = this->LocalWidth();
            const Int localHeightOfA = A.LocalHeight();

            const Int maxHeight = MaxLocalLength(height,p);
            const Int maxWidth = MaxLocalLength(width,c);
            const Int portionSize = 
                std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

            this->auxMemory_.Require( 2*c*portionSize );

            T* buffer = this->auxMemory_.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[c*portionSize];

            // Pack
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for  
#endif
            for( Int k=0; k<c; ++k )
            {
                T* data = &sendBuffer[k*portionSize];

                const Int thisRowShift = RawShift(k,rowAlignment,c);
                const Int thisLocalWidth = RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for 
#endif
                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*c)*ALDim];
                    T* dataCol = &data[jLocal*localHeightOfA];
                    std::memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                }
            }

            // Communicate
            mpi::AllToAll
            ( sendBuffer, portionSize,
              recvBuffer, portionSize, g.MRComm() );

            // Unpack
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int k=0; k<c; ++k )
            {
                const T* data = &recvBuffer[k*portionSize];

                const Int thisRank = row+k*r;
                const Int thisColShift = RawShift(thisRank,colAlignmentOfA,p);
                const Int thisColOffset = (thisColShift-colShift) / r;
                const Int thisLocalHeight = 
                    RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const Int i = thisColOffset + iLocal*c;
                        thisLocalBuffer[i+jLocal*thisLDim] = 
                            data[iLocal+jLocal*thisLocalHeight];
                    }
                }
            }
            this->auxMemory_.Release();
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR] <- [VC,* ]." << std::endl;
#endif
            const Int r = g.Height();
            const Int c = g.Width();
            const Int p = r * c;
            const Int row = g.MCRank();
            const Int colShift = this->ColShift();
            const Int colAlignment = this->ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int colAlignmentOfA = A.ColAlignment();

            const Int sendRow = (row+r+colAlignment-(colAlignmentOfA%r)) % r;
            const Int recvRow = (row+r+(colAlignmentOfA%r)-colAlignment) % r;

            const Int height = this->Height();
            const Int width = this->Width();
            const Int localWidth = this->LocalWidth();
            const Int localHeightOfA = A.LocalHeight();

            const Int maxHeight = MaxLocalLength(height,p);
            const Int maxWidth = MaxLocalLength(width,c);
            const Int portionSize = 
                std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

            this->auxMemory_.Require( 2*c*portionSize );

            T* buffer = this->auxMemory_.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[c*portionSize];

            // Pack
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                T* data = &secondBuffer[k*portionSize];

                const Int thisRowShift = RawShift(k,rowAlignment,c);
                const Int thisLocalWidth = RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*c)*ALDim];
                    T* dataCol = &data[jLocal*localHeightOfA];
                    std::memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
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
            const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                const T* data = &secondBuffer[k*portionSize];

                const Int thisRank = recvRow+k*r;
                const Int thisColShift = RawShift(thisRank,colAlignmentOfA,p);
                const Int thisColOffset = (thisColShift-colShift) / r;
                const Int thisLocalHeight = 
                    RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const Int i = thisColOffset + iLocal*c;
                        thisLocalBuffer[i+jLocal*thisLDim] = 
                            data[iLocal+jLocal*thisLocalHeight];
                    }
                }
            }
            this->auxMemory_.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR,Int> A_STAR_VR(true,this->RowAlignment(),g);

    A_STAR_VR = A;
    *this = A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [VR,* ]");
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
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
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
            this->rowAlignment_ = A.RowAlignment() % g.Width();
            if( g.InGrid() )
            {
                this->rowShift_ = 
                    Shift( g.MRRank(), this->RowAlignment(), g.Width() );
            }
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( g.InGrid() )
    {
        if( this->RowAlignment() == A.RowAlignment() % g.Width() )
        {
            const Int r = g.Height();
            const Int c = g.Width();
            const Int p = r * c;
            const Int col = g.MRRank();
            const Int rowShift = this->RowShift();
            const Int colAlignment = this->ColAlignment();
            const Int rowAlignmentOfA = A.RowAlignment();
    
            const Int height = this->Height();
            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidthOfA = A.LocalWidth();

            const Int maxHeight = MaxLocalLength(height,r);
            const Int maxWidth = MaxLocalLength(width,p);
            const Int portionSize = 
                std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

            this->auxMemory_.Require( 2*r*portionSize );

            T* buffer = this->auxMemory_.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[r*portionSize];

            // Pack
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &sendBuffer[k*portionSize];

                const Int thisColShift = RawShift(k,colAlignment,r);
                const Int thisLocalHeight = 
                    RawLocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                {
                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const Int i = thisColShift + iLocal*r;
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
            const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                const T* data = &recvBuffer[k*portionSize];

                const Int thisRank = col+k*c;
                const Int thisRowShift = RawShift(thisRank,rowAlignmentOfA,p);
                const Int thisRowOffset = (thisRowShift-rowShift) / c;
                const Int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* dataCol = &data[jLocal*localHeight];
                    T* thisCol = 
                        &thisLocalBuffer[(thisRowOffset+jLocal*r)*thisLDim];
                    std::memcpy( thisCol, dataCol, localHeight*sizeof(T) );
                }
            }
            this->auxMemory_.Release();
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned [MC,MR] <- [* ,VR]." << std::endl;
#endif
            const Int r = g.Height();
            const Int c = g.Width();
            const Int p = r * c;
            const Int col = g.MRRank();
            const Int rowShift = this->RowShift();
            const Int colAlignment = this->ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int rowAlignmentOfA = A.RowAlignment();

            const Int sendCol = (col+c+rowAlignment-(rowAlignmentOfA%c)) % c;
            const Int recvCol = (col+c+(rowAlignmentOfA%c)-rowAlignment) % c;

            const Int height = this->Height();
            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidthOfA = A.LocalWidth();
            
            const Int maxHeight = MaxLocalLength(height,r);
            const Int maxWidth = MaxLocalLength(width,p);
            const Int portionSize = 
                std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

            this->auxMemory_.Require( 2*r*portionSize );

            T* buffer = this->auxMemory_.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[r*portionSize];

            // Pack
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &secondBuffer[k*portionSize];

                const Int thisColShift = RawShift(k,colAlignment,r);
                const Int thisLocalHeight = 
                    RawLocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                {
                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const Int i = thisColShift + iLocal*r;
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
            const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                const T* data = &secondBuffer[k*portionSize];

                const Int thisRank = recvCol+k*c;
                const Int thisRowShift = RawShift(thisRank,rowAlignmentOfA,p);
                const Int thisRowOffset = (thisRowShift-rowShift) / c;
                const Int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* dataCol = &data[jLocal*localHeight];
                    T* thisCol = 
                        &thisLocalBuffer[(thisRowOffset+jLocal*r)*thisLDim];
                    std::memcpy( thisCol, dataCol, localHeight*sizeof(T) );
                }
            }
            this->auxMemory_.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,MR,Int>&
DistMatrix<T,MC,MR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
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

    const Int r = this->Grid().Height();
    const Int c = this->Grid().Width();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();

    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
            thisLocalBuffer[iLocal+jLocal*thisLDim] = 
                ALocalBuffer[(colShift+iLocal*r)+(rowShift+jLocal*c)*ALDim];
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SumScatterFrom( const DistMatrix<T,MC,STAR,Int>& A )
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
            if( this->Width() == 1 )
            {
                const Int rowAlignment = this->RowAlignment();
                const Int myCol = g.MRRank();

                const Int localHeight = this->LocalHeight();

                const Int recvSize = std::max(localHeight,mpi::MIN_COLL_MSG);
                const Int sendSize = recvSize;

                this->auxMemory_.Require( sendSize + recvSize );

                T* buffer = this->auxMemory_.Buffer();
                T* sendBuffer = &buffer[0];
                T* recvBuffer = &buffer[sendSize];

                // Pack 
                const T* ACol = A.LockedLocalBuffer(0,0);
                std::memcpy( sendBuffer, ACol, localHeight*sizeof(T) );

                // Reduce to rowAlignment
                mpi::Reduce
                ( sendBuffer, recvBuffer, sendSize, 
                  mpi::SUM, rowAlignment, g.MRComm() );

                if( myCol == rowAlignment )
                {
                    T* thisCol = this->LocalBuffer(0,0);
                    std::memcpy( thisCol, recvBuffer, localHeight*sizeof(T) );
                }

                this->auxMemory_.Release();
            }
            else
            {
                const Int c = g.Width();
                const Int rowAlignment = this->RowAlignment();
            
                const Int width = this->Width();
                const Int localHeight = this->LocalHeight();
                const Int localWidth = this->LocalWidth();
                const Int maxLocalWidth = MaxLocalLength(width,c);

                const Int recvSize = 
                    std::max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
                const Int sendSize = c * recvSize;

                this->auxMemory_.Require( sendSize );
                T* buffer = this->auxMemory_.Buffer();
            
                // Pack 
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int k=0; k<c; ++k )
                {
                    T* data = &buffer[k*recvSize];

                    const Int thisRowShift = RawShift( k, rowAlignment, c );
                    const Int thisLocalWidth = 
                        RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                    #pragma omp parallel for
#endif
                    for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    {
                        const T* ACol = 
                            &ALocalBuffer[(thisRowShift+jLocal*c)*ALDim];
                        T* dataCol = &data[jLocal*localHeight];
                        std::memcpy( dataCol, ACol, localHeight*sizeof(T) );
                    }
                }

                // Communicate
                mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.MRComm() );

                // Unpack our received data
                T* thisLocalBuffer = this->LocalBuffer();
                const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    const T* bufferCol = &buffer[jLocal*localHeight];
                    T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                    std::memcpy( thisCol, bufferCol, localHeight*sizeof(T) );
                }
                this->auxMemory_.Release();
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
                const Int r = g.Height();
                const Int rowAlignment = this->RowAlignment();
                const Int myRow = g.MCRank();
                const Int myCol = g.MRRank();

                const Int height = this->Height();
                const Int localHeight = this->LocalHeight();
                const Int localHeightOfA = A.LocalHeight();
                const Int maxLocalHeight = MaxLocalLength(height,r);

                const Int portionSize = 
                    std::max(maxLocalHeight,mpi::MIN_COLL_MSG);

                const Int colAlignment = this->ColAlignment();
                const Int colAlignmentOfA = A.ColAlignment();
                const Int sendRow = (myRow+r+colAlignment-colAlignmentOfA) % r;
                const Int recvRow = (myRow+r+colAlignmentOfA-colAlignment) % r;

                this->auxMemory_.Require( 2*portionSize );

                T* buffer = this->auxMemory_.Buffer();
                T* sendBuffer = &buffer[0];
                T* recvBuffer = &buffer[portionSize];

                // Pack 
                const T* ACol = A.LockedLocalBuffer(0,0);
                std::memcpy( sendBuffer, ACol, localHeightOfA*sizeof(T) );
            
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
                    std::memcpy( thisCol, sendBuffer, localHeight*sizeof(T) );
                }

                this->auxMemory_.Release();
            }
            else
            {
                const Int r = g.Height();
                const Int c = g.Width();
                const Int row = g.MCRank();

                const Int colAlignment = this->ColAlignment();
                const Int rowAlignment = this->RowAlignment();
                const Int colAlignmentOfA = A.ColAlignment();
                const Int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
                const Int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

                const Int width = this->Width();
                const Int localHeight = this->LocalHeight();
                const Int localWidth = this->LocalWidth();
                const Int localHeightOfA = A.LocalHeight();
                const Int maxLocalWidth = MaxLocalLength(width,c);

                const Int recvSize_RS = 
                    std::max(localHeightOfA*maxLocalWidth,mpi::MIN_COLL_MSG);
                const Int sendSize_RS = c * recvSize_RS;
                const Int recvSize_SR = localHeight * localWidth;

                this->auxMemory_.Require
                ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );

                T* buffer = this->auxMemory_.Buffer();
                T* firstBuffer = &buffer[0];
                T* secondBuffer = &buffer[recvSize_RS];

                // Pack 
                // TODO: Stick an optional outer parallelization here?
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const Int ALDim = A.LocalLDim();
                for( Int k=0; k<c; ++k )
                {
                    T* data = &secondBuffer[k*recvSize_RS];

                    const Int thisRowShift = RawShift( k, rowAlignment, c );
                    const Int thisLocalWidth = 
                        RawLocalLength(width,thisRowShift,c);

#ifdef _OPENMP
                    #pragma omp parallel for
#endif
                    for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    {
                        const Int j = thisRowShift + jLocal*c;
                        const T* ACol = &ALocalBuffer[j*ALDim];
                        T* dataCol = &data[jLocal*localHeightOfA];
                        std::memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                    }
                }

                // Reduce-scatter over each process row
                mpi::ReduceScatter
                ( secondBuffer, firstBuffer, recvSize_RS, mpi::SUM, 
                  g.MRComm() );

                // Trade reduced data with the appropriate process row
                mpi::SendRecv
                ( firstBuffer,  localHeightOfA*localWidth, sendRow, 0,
                  secondBuffer, localHeight*localWidth,    recvRow, 0, 
                  g.MCComm() );

                // Unpack the received data
                T* thisLocalBuffer = this->LocalBuffer();
                const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    const T* secondBufferCol = 
                        &secondBuffer[jLocal*localHeight];
                    T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                    std::memcpy
                    ( thisCol, secondBufferCol, localHeight*sizeof(T) );
                }
                this->auxMemory_.Release();
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SumScatterFrom( const DistMatrix<T,STAR,MR,Int>& A )
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
            this->rowAlignment_ = A.RowAlignment();
            if( g.InGrid() )
            {
                this->rowShift_ = 
                    Shift( g.MRRank(), this->RowAlignment(), g.Width() );
            }
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( g.InGrid() )
    {
        if( this->RowAlignment() == A.RowAlignment() )
        {
            const Int r = g.Height();
            const Int colAlignment = this->ColAlignment();

            const Int height = this->Height();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int maxLocalHeight = MaxLocalLength(height,r);

            const Int recvSize = 
                std::max(maxLocalHeight*localWidth,mpi::MIN_COLL_MSG);
            const Int sendSize = r * recvSize;

            this->auxMemory_.Require( sendSize );
            T* buffer = this->auxMemory_.Buffer();

            // Pack 
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &buffer[k*recvSize];

                const Int thisColShift = RawShift( k, colAlignment, r );
                const Int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const Int i = thisColShift + iLocal*r;
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[i+jLocal*ALDim];
                    }
                }
            }

            // Communicate
            mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.MCComm() );

            // Unpack our received data
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* bufferCol = &buffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                std::memcpy( thisCol, bufferCol, localHeight*sizeof(T) );
            }
            this->auxMemory_.Release();
        }
        else
        {
#ifdef UNALIGNED_WARNINGS
            if( g.VCRank() == 0 )
                std::cerr << "Unaligned SumScatterFrom [MC,MR] <- [* ,MR]." 
                          << std::endl;
#endif
            const Int r = g.Height();
            const Int c = g.Width();
            const Int col = g.MRRank();

            const Int colAlignment = this->ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int rowAlignmentOfA = A.RowAlignment();
            const Int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
            const Int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

            const Int height = this->Height();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localWidthOfA = A.LocalWidth();
            const Int maxLocalHeight = MaxLocalLength(height,r);

            const Int recvSize_RS = 
                std::max(maxLocalHeight*localWidthOfA,mpi::MIN_COLL_MSG);
            const Int sendSize_RS = r * recvSize_RS;
            const Int recvSize_SR = localHeight * localWidth;

            this->auxMemory_.Require
                ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );

            T* buffer = this->auxMemory_.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[recvSize_RS];

            // Pack 
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &secondBuffer[k*recvSize_RS];

                const Int thisColShift = RawShift( k, colAlignment, r );
                const Int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                {
                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    {
                        const Int i = thisColShift + iLocal*r;
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[i+jLocal*ALDim];
                    }
                }
            }

            // Reduce-scatter over each process col
            mpi::ReduceScatter
            ( secondBuffer, firstBuffer, recvSize_RS, mpi::SUM, g.MCComm() );

            // Trade reduced data with the appropriate process col
            mpi::SendRecv
            ( firstBuffer,  localHeight*localWidthOfA, sendCol, 0,
              secondBuffer, localHeight*localWidth,    recvCol, 0, g.MRComm() );

            // Unpack the received data
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* secondBufferCol = &secondBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                std::memcpy( thisCol, secondBufferCol, localHeight*sizeof(T) );
            }
            this->auxMemory_.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SumScatterFrom
( const DistMatrix<T,STAR,STAR,Int>& A )
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
        const Int r = g.Height();
        const Int c = g.Width();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalHeight = MaxLocalLength(height,r);
        const Int maxLocalWidth = MaxLocalLength(width,c);

        const Int recvSize = 
            std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
        const Int sendSize = r*c*recvSize;

        this->auxMemory_.Require( sendSize );
        T* buffer = this->auxMemory_.Buffer();

        // Pack 
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int l=0; l<c; ++l )
        {
            const Int thisRowShift = RawShift( l, rowAlignment, c );
            const Int thisLocalWidth = RawLocalLength( width, thisRowShift, c );

            for( Int k=0; k<r; ++k )
            {
                T* data = &buffer[(k+l*r)*recvSize];

                const Int thisColShift = RawShift( k, colAlignment, r );
                const Int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[(thisColShift+iLocal*r)+
                                         (thisRowShift+jLocal*c)*ALDim];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.VCComm() );

        // Unpack our received data
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* bufferCol = &buffer[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            std::memcpy( thisCol, bufferCol, localHeight*sizeof(T) );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,MC,STAR,Int>& A )
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
                const Int rowAlignment = this->RowAlignment();
                const Int myCol = g.MRRank();

                const Int localHeight = this->LocalHeight();

                const Int portionSize = std::max(localHeight,mpi::MIN_COLL_MSG);

                this->auxMemory_.Require( 2*portionSize );

                T* buffer = this->auxMemory_.Buffer();
                T* sendBuffer = &buffer[0];
                T* recvBuffer = &buffer[portionSize];

                // Pack 
                const T* ACol = A.LockedLocalBuffer(0,0);
                std::memcpy( sendBuffer, ACol, localHeight*sizeof(T) );
            
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
                    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisCol[iLocal] += alpha*recvBuffer[iLocal];
                }

                this->auxMemory_.Release();
            }
            else
            {
                const Int c = g.Width();
                const Int rowAlignment = this->RowAlignment();

                const Int width = this->Width();
                const Int localHeight = this->LocalHeight();
                const Int localWidth = this->LocalWidth();
                const Int maxLocalWidth = MaxLocalLength(width,c);

                const Int portionSize = 
                    std::max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
                const Int sendSize = c*portionSize;

                this->auxMemory_.Require( sendSize );
                T* buffer = this->auxMemory_.Buffer();

                // Pack 
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int k=0; k<c; ++k )
                {
                    T* data = &buffer[k*portionSize];

                    const Int thisRowShift = RawShift( k, rowAlignment, c );
                    const Int thisLocalWidth = 
                        RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                    #pragma omp parallel for
#endif
                    for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    {
                        const T* ACol = 
                            &ALocalBuffer[(thisRowShift+jLocal*c)*ALDim];
                        T* dataCol = &data[jLocal*localHeight];
                        std::memcpy( dataCol, ACol, localHeight*sizeof(T) );
                    }
                }
            
                // Communicate
                mpi::ReduceScatter( buffer, portionSize, mpi::SUM, g.MRComm() );

                // Update with our received data
                T* thisLocalBuffer = this->LocalBuffer();
                const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    const T* bufferCol = &buffer[jLocal*localHeight];
                    T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                    blas::Axpy( localHeight, alpha, bufferCol, 1, thisCol, 1 );
                }
                this->auxMemory_.Release();
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
                const Int r = g.Height();
                const Int rowAlignment = this->RowAlignment();
                const Int myRow = g.MCRank();
                const Int myCol = g.MRRank();

                const Int height = this->Height();
                const Int localHeight = this->LocalHeight();
                const Int localHeightOfA = A.LocalHeight();
                const Int maxLocalHeight = MaxLocalLength(height,r);

                const Int portionSize = 
                    std::max(maxLocalHeight,mpi::MIN_COLL_MSG);

                const Int colAlignment = this->ColAlignment();
                const Int colAlignmentOfA = A.ColAlignment();
                const Int sendRow = (myRow+r+colAlignment-colAlignmentOfA) % r;
                const Int recvRow = (myRow+r+colAlignmentOfA-colAlignment) % r;

                this->auxMemory_.Require( 2*portionSize );

                T* buffer = this->auxMemory_.Buffer();
                T* sendBuffer = &buffer[0];
                T* recvBuffer = &buffer[portionSize];

                // Pack 
                const T* ACol = A.LockedLocalBuffer(0,0);
                std::memcpy( sendBuffer, ACol, localHeightOfA*sizeof(T) );
            
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
                    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisCol[iLocal] += alpha*sendBuffer[iLocal];
                }
                this->auxMemory_.Release();
            }
            else
            {
                const Int r = g.Height();
                const Int c = g.Width();
                const Int row = g.MCRank();

                const Int colAlignment = this->ColAlignment();
                const Int rowAlignment = this->RowAlignment();
                const Int colAlignmentOfA = A.ColAlignment();
                const Int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
                const Int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

                const Int width = this->Width();
                const Int localHeight = this->LocalHeight();
                const Int localWidth = this->LocalWidth();
                const Int localHeightOfA = A.LocalHeight();
                const Int maxLocalWidth = MaxLocalLength(width,c);

                const Int recvSize_RS = 
                    std::max(localHeightOfA*maxLocalWidth,mpi::MIN_COLL_MSG);
                const Int sendSize_RS = c * recvSize_RS;
                const Int recvSize_SR = localHeight * localWidth;

                this->auxMemory_.Require
                ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );

                T* buffer = this->auxMemory_.Buffer();
                T* firstBuffer = &buffer[0];
                T* secondBuffer = &buffer[recvSize_RS];

                // Pack 
                const T* ALocalBuffer = A.LockedLocalBuffer();
                const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int k=0; k<c; ++k )
                {
                    T* data = &secondBuffer[k*recvSize_RS];

                    const Int thisRowShift = RawShift( k, rowAlignment, c );
                    const Int thisLocalWidth = 
                        RawLocalLength(width,thisRowShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                    #pragma omp parallel for
#endif
                    for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    {
                        const T* ACol = 
                            &ALocalBuffer[(thisRowShift+jLocal*c)*ALDim];
                        T* dataCol = &data[jLocal*localHeightOfA];
                        std::memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                    }
                }

                // Reduce-scatter over each process row
                mpi::ReduceScatter
                ( secondBuffer, firstBuffer, recvSize_RS, mpi::SUM, 
                  g.MRComm() );

                // Trade reduced data with the appropriate process row
                mpi::SendRecv
                ( firstBuffer,  localHeightOfA*localWidth, sendRow, 0,
                  secondBuffer, localHeight*localWidth,    recvRow, 0, 
                  g.MCComm() );

                // Update with our received data
                T* thisLocalBuffer = this->LocalBuffer();
                const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    const T* secondBufferCol = 
                        &secondBuffer[jLocal*localHeight];
                    T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisCol[iLocal] += alpha*secondBufferCol[iLocal];
                }
                this->auxMemory_.Release();
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MR,Int>& A )
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
            const Int r = g.Height();
            const Int colAlignment = this->ColAlignment();

            const Int height = this->Height();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int maxLocalHeight = MaxLocalLength(height,r);

            const Int recvSize = 
                std::max(maxLocalHeight*localWidth,mpi::MIN_COLL_MSG);
            const Int sendSize = r*recvSize;

            this->auxMemory_.Require( sendSize );
            T* buffer = this->auxMemory_.Buffer();

            // Pack 
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &buffer[k*recvSize];

                const Int thisColShift = RawShift( k, colAlignment, r );
                const Int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[(thisColShift+iLocal*r)+jLocal*ALDim];
            }

            // Communicate
            mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.MCComm() );

            // Update with our received data
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* bufferCol = &buffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisCol[iLocal] += alpha*bufferCol[iLocal];
            }
            this->auxMemory_.Release();
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
            const Int r = g.Height();
            const Int c = g.Width();
            const Int col = g.MRRank();

            const Int colAlignment = this->ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int rowAlignmentOfA = A.RowAlignment();
            const Int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
            const Int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

            const Int height = this->Height();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localWidthOfA = A.LocalWidth();
            const Int maxLocalHeight = MaxLocalLength(height,r);

            const Int recvSize_RS = 
                std::max(maxLocalHeight*localWidthOfA,mpi::MIN_COLL_MSG);
            const Int sendSize_RS = r * recvSize_RS;
            const Int recvSize_SR = localHeight * localWidth;

            this->auxMemory_.Require
                ( recvSize_RS + std::max(sendSize_RS,recvSize_SR) );

            T* buffer = this->auxMemory_.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[recvSize_RS];

            // Pack
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &secondBuffer[k*recvSize_RS];

                const Int thisColShift = RawShift( k, colAlignment, r );
                const Int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[(thisColShift+iLocal*r)+jLocal*ALDim];
            }

            // Reduce-scatter over each process col
            mpi::ReduceScatter
            ( secondBuffer, firstBuffer, recvSize_RS, mpi::SUM, g.MCComm() );

            // Trade reduced data with the appropriate process col
            mpi::SendRecv
            ( firstBuffer,  localHeight*localWidthOfA, sendCol, 0,
              secondBuffer, localHeight*localWidth,    recvCol, mpi::ANY_TAG,
              g.MRComm() );

            // Update with our received data
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* secondBufferCol = &secondBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisCol[iLocal] += alpha*secondBufferCol[iLocal];
            }
            this->auxMemory_.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR,Int>& A )
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
        const Int r = g.Height();
        const Int c = g.Width();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalHeight = MaxLocalLength(height,r);
        const Int maxLocalWidth = MaxLocalLength(width,c);

        const Int recvSize = 
            std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
        const Int sendSize = r * c * recvSize;

        this->auxMemory_.Require( sendSize );
        T* buffer = this->auxMemory_.Buffer();

        // Pack 
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int l=0; l<c; ++l )
        {
            const Int thisRowShift = RawShift( l, rowAlignment, c );
            const Int thisLocalWidth = RawLocalLength( width, thisRowShift, c );

            for( Int k=0; k<r; ++k )
            {
                T* data = &buffer[(k+l*r)*recvSize];

                const Int thisColShift = RawShift( k, colAlignment, r );
                const Int thisLocalHeight = 
                    RawLocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        data[iLocal+jLocal*thisLocalHeight] = 
                            ALocalBuffer[(thisColShift+iLocal*r)+
                                         (thisRowShift+jLocal*c)*ALDim];
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.VCComm() );

        // Unpack our received data
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* bufferCol = &buffer[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                thisCol[iLocal] += alpha*bufferCol[iLocal];
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
