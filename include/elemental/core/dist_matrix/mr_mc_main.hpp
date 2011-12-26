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
inline
DistMatrix<T,MR,MC,Int>::DistMatrix( const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   (g.InGrid() ? g.MRRank() : 0),
   (g.InGrid() ? g.MCRank() : 0),
  0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MR,MC,Int>::DistMatrix
( Int height, Int width, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   (g.InGrid() ? g.MRRank() : 0),
   (g.InGrid() ? g.MCRank() : 0),
   (g.InGrid() ? LocalLength(height,g.MRRank(),0,g.Width()) : 0),
   (g.InGrid() ? LocalLength(width,g.MCRank(),0,g.Height()) : 0),
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MR,MC,Int>::DistMatrix
( bool constrainedColAlignment, bool constrainedRowAlignment,
  Int colAlignment, Int rowAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MR,MC,Int>::DistMatrix
( Int height, Int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  Int colAlignment, Int rowAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MR,MC,Int>::DistMatrix
( Int height, Int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  Int colAlignment, Int rowAlignment, Int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MR,MC,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  const T* buffer, Int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MR,MC,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  T* buffer, Int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,MR,MC,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix");
#endif
    if( MR != U || MC != V || 
        reinterpret_cast<const DistMatrix<T,MR,MC,Int>*>(&A) != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [MR,MC] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,MR,MC,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::SetGrid( const elemental::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->colAlignment_ = 0;
    this->rowAlignment_ = 0;
    this->colShift_ = grid.MRRank();
    this->rowShift_ = grid.MCRank();
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->rowAlignment_ = A.RowAlignment();
    this->colShift_     = A.ColShift();
    this->rowShift_     = A.RowShift();
    this->constrainedColAlignment_ = true;
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignWith( const DistMatrix<S,MR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ = A.ColShift();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignWith( const DistMatrix<S,STAR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([* ,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment();
    this->rowShift_ = A.RowShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->rowAlignment_ = A.ColAlignment();
    this->colShift_     = A.RowShift();
    this->rowShift_     = A.ColShift();
    this->constrainedColAlignment_ = true;
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignWith( const DistMatrix<S,MC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MC,*])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment();
    this->rowShift_ = A.ColShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignWith( const DistMatrix<S,STAR,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([* ,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ = A.RowShift();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int> 
template<typename S,typename N>
inline void 
DistMatrix<T,MR,MC,Int>::AlignWith( const DistMatrix<S,VC,STAR,N>& A ) 
{ 
#ifndef RELEASE 
    PushCallStack("[MR,MC]::AlignWith([VC,* ])"); 
    this->AssertFreeRowAlignment(); 
    this->AssertSameGrid( A ); 
#endif 
    const elemental::Grid& g = this->Grid(); 
    this->rowAlignment_ = A.ColAlignment(); 
    this->rowShift_ =  
        Shift( g.MCRank(), this->RowAlignment(), g.Height() ); 
    this->constrainedRowAlignment_ = true; 
    this->height_ = 0; 
    this->width_ = 0; 
    this->localMatrix_.ResizeTo( 0, 0 ); 
#ifndef RELEASE 
    PopCallStack(); 
#endif 
} 
 
template<typename T,typename Int> 
template<typename S,typename N>
inline void 
DistMatrix<T,MR,MC,Int>::AlignWith ( const DistMatrix<S,STAR,VC,N>& A ) 
{ 
#ifndef RELEASE 
    PushCallStack("[MR,MC]:AlignWith([* ,VC])"); 
    this->AssertFreeRowAlignment(); 
    this->AssertSameGrid( A ); 
#endif 
    const elemental::Grid& g = this->Grid(); 
    this->rowAlignment_ = A.RowAlignment(); 
    this->rowShift_ =  
        Shift( g.MCRank(), this->RowAlignment(), g.Height() ); 
    this->constrainedRowAlignment_ = true; 
    this->height_ = 0; 
    this->width_ = 0; 
    this->localMatrix_.ResizeTo( 0, 0 ); 
#ifndef RELEASE 
    PopCallStack(); 
#endif 
} 
 
template<typename T,typename Int> 
template<typename S,typename N>
inline void 
DistMatrix<T,MR,MC,Int>::AlignWith ( const DistMatrix<S,VR,STAR,N>& A ) 
{ 
#ifndef RELEASE 
    PushCallStack("[MR,MC]::AlignWith([VR,* ])"); 
    this->AssertFreeColAlignment(); 
    this->AssertSameGrid( A ); 
#endif 
    const elemental::Grid& g = this->Grid(); 
    this->colAlignment_ = A.ColAlignment(); 
    this->colShift_ =  
        Shift( g.MRRank(), this->ColAlignment(), g.Width() ); 
    this->constrainedColAlignment_ = true; 
    this->height_ = 0; 
    this->width_ = 0; 
    this->localMatrix_.ResizeTo( 0, 0 ); 
#ifndef RELEASE 
    PopCallStack(); 
#endif 
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignWith( const DistMatrix<S,STAR,VR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignWith([* ,VR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignColsWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ = A.ColShift();
    this->constrainedColAlignment_ = true; 
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignColsWith( const DistMatrix<S,MR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([MR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ = A.ColShift();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignColsWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ = A.RowShift();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignColsWith( const DistMatrix<S,STAR,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([* ,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ = A.RowShift();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignColsWith( const DistMatrix<S,VR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([VR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ =
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignColsWith( const DistMatrix<S,STAR,VR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignColsWith([* ,VR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ =
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignRowsWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment();
    this->rowShift_ = A.RowShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignRowsWith( const DistMatrix<S,STAR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([* ,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment();
    this->rowShift_ = A.RowShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignRowsWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment();
    this->rowShift_ = A.ColShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignRowsWith( const DistMatrix<S,MC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([MC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment();
    this->rowShift_ = A.ColShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignRowsWith( const DistMatrix<S,VC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([VC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.ColAlignment();
    this->rowShift_ =
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MR,MC,Int>::AlignRowsWith( const DistMatrix<S,STAR,VC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignRowsWith([* ,VC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.RowAlignment();
    this->rowShift_ =
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

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
            sendBuf[(colShift+iLocal*c) + (rowShift+jLocal*r)*height] = 
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::Align( Int colAlignment, Int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Align");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Width() )
        throw std::runtime_error("Invalid column alignment for [MR,MC]");
    if( rowAlignment < 0 || rowAlignment >= g.Height() )
        throw std::runtime_error("Invalid row alignment for [MR,MC]");
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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::AlignCols( Int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Width() )
        throw std::runtime_error("Invalid column alignment for [MR,MC]");
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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::AlignRows( Int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Height() )
        throw std::runtime_error("Invalid row alignment for [MR,MC]");
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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::View( DistMatrix<T,MR,MC,Int>& A )
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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::View
( Int height, Int width, Int colAlignment, Int rowAlignment,
  T* buffer, Int ldim, const elemental::Grid& grid )
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
    const Int localHeight = LocalLength(height,this->colShift_,grid.Width());
    const Int localWidth = LocalLength(width,this->rowShift_,grid.Height());
    this->localMatrix_.View( localHeight, localWidth, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::LockedView( const DistMatrix<T,MR,MC,Int>& A )
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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::LockedView
( Int height, Int width, Int colAlignment, Int rowAlignment,
  const T* buffer, Int ldim, const elemental::Grid& grid )
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
    const Int localHeight = LocalLength(height,this->colShift_,grid.Width());
    const Int localWidth = LocalLength(width,this->rowShift_,grid.Height());
    this->localMatrix_.LockedView( localHeight, localWidth, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::View
( DistMatrix<T,MR,MC,Int>& A, Int i, Int j, Int height, Int width )
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
        const Int r   = g.Height();
        const Int c   = g.Width();
        const Int row = g.MCRank();
        const Int col = g.MRRank();

        this->colAlignment_ = (A.ColAlignment()+i) % c;
        this->rowAlignment_ = (A.RowAlignment()+j) % r;
        
        this->colShift_ = Shift( col, this->ColAlignment(), c );
        this->rowShift_ = Shift( row, this->RowAlignment(), r );

        const Int localHeightBefore = LocalLength( i, A.ColShift(), c );
        const Int localWidthBefore  = LocalLength( j, A.RowShift(), r );

        const Int localHeight = LocalLength( height, this->ColShift(), c );
        const Int localWidth  = LocalLength( width,  this->RowShift(), r );

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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::LockedView
( const DistMatrix<T,MR,MC,Int>& A, Int i, Int j, Int height, Int width )
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
        const Int r   = g.Height();
        const Int c   = g.Width();
        const Int row = g.MCRank();
        const Int col = g.MRRank();

        this->colAlignment_ = (A.ColAlignment()+i) % c;
        this->rowAlignment_ = (A.RowAlignment()+j) % r;
        
        this->colShift_ = Shift( col, this->ColAlignment(), c );
        this->rowShift_ = Shift( row, this->RowAlignment(), r );

        const Int localHeightBefore = LocalLength( i, A.ColShift(), c );
        const Int localWidthBefore  = LocalLength( j, A.RowShift(), r );

        const Int localHeight = LocalLength( height, this->ColShift(), c );
        const Int localWidth  = LocalLength( width,  this->RowShift(), r );

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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::View1x2
( DistMatrix<T,MR,MC,Int>& AL, DistMatrix<T,MR,MC,Int>& AR )
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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::LockedView1x2
( const DistMatrix<T,MR,MC,Int>& AL, const DistMatrix<T,MR,MC,Int>& AR )
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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::View2x1
( DistMatrix<T,MR,MC,Int>& AT,
  DistMatrix<T,MR,MC,Int>& AB )
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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::LockedView2x1
( const DistMatrix<T,MR,MC,Int>& AT,
  const DistMatrix<T,MR,MC,Int>& AB )
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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::View2x2
( DistMatrix<T,MR,MC,Int>& ATL, DistMatrix<T,MR,MC,Int>& ATR,
  DistMatrix<T,MR,MC,Int>& ABL, DistMatrix<T,MR,MC,Int>& ABR )
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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::LockedView2x2
( const DistMatrix<T,MR,MC,Int>& ATL, const DistMatrix<T,MR,MC,Int>& ATR,
  const DistMatrix<T,MR,MC,Int>& ABL, const DistMatrix<T,MR,MC,Int>& ABR )
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

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
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

template<typename T,typename Int>
inline T
DistMatrix<T,MR,MC,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const elemental::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol * g.Height();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        const Int jLoc = (j-this->RowShift()) / g.Height();
        u = this->GetLocalEntry(iLoc,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->SetLocalEntry(iLoc,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();
    const Int ownerCol = (i + this->ColAlignment()) % g.Width();
    const Int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Width();
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateLocalEntry(iLoc,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::GetDiagonal
( DistMatrix<T,MD,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetDiagonal([MD,* ])");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const Int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Height() || d.Width() != 1) )
        throw std::logic_error("d is not of the correct dimensions");
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

        const Int iLocalStart = (iStart-colShift) / c;
        const Int jLocalStart = (jStart-rowShift) / r;

        const Int localDiagLength = d.LocalHeight();

        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const Int thisLDim = this->LocalLDim();
        T* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart + k*(lcm/c);
            const Int jLocal = jLocalStart + k*(lcm/r);
            dLocalBuffer[k] = thisLocalBuffer[iLocal+jLocal*thisLDim]; 
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::GetDiagonal
( DistMatrix<T,STAR,MD,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetDiagonal([* ,MD])");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const Int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Width() || d.Height() != 1) )
        throw std::logic_error("d is not of the correct dimensions");
    if( ( d.Viewing() && d.ConstrainedRowAlignment() ) &&
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

        const Int iLocalStart = (iStart-colShift) / c;
        const Int jLocalStart = (jStart-rowShift) / r;

        const Int localDiagLength = d.LocalWidth();

        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const Int thisLDim = this->LocalLDim();
        T* dLocalBuffer = d.LocalBuffer();
        const Int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart + k*(lcm/c);
            const Int jLocal = jLocalStart + k*(lcm/r);
            dLocalBuffer[k*dLDim] = thisLocalBuffer[iLocal+jLocal*thisLDim];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::SetDiagonal
( const DistMatrix<T,MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetDiagonal([MD,* ])");
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

        const Int iLocalStart = (iStart-colShift) / c;
        const Int jLocalStart = (jStart-rowShift) / r;

        const Int localDiagLength = d.LocalHeight();

        const T* dLocalBuffer = d.LockedLocalBuffer();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart + k*(lcm/c);
            const Int jLocal = jLocalStart + k*(lcm/r);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::SetDiagonal
( const DistMatrix<T,STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetDiagonal([* ,MD])");
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

        const Int iLocalStart = (iStart-colShift) / c;
        const Int jLocalStart = (jStart-rowShift) / r;

        const Int localDiagLength = d.LocalWidth();

        const T* dLocalBuffer = d.LockedLocalBuffer();
        const Int dLDim = d.LocalLDim();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart + k*(lcm/c);
            const Int jLocal = jLocalStart + k*(lcm/r);
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
DistMatrix<T,MR,MC,Int>::MakeTrapezoidal
( Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::MakeTrapezoidal");
    this->AssertNotLockedView(); 
#endif
    const elemental::Grid& g = this->Grid();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

    if( uplo == LOWER )
    {
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            Int j = rowShift + jLocal*r;
            Int lastZeroRow = ( side==LEFT ? j-offset-1
                                           : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                Int boundary = std::min( lastZeroRow+1, height );
                Int numZeroRows = RawLocalLength( boundary, colShift, c );
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                std::memset( thisCol, 0, numZeroRows*sizeof(T) );
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
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            Int j = rowShift + jLocal*r;
            Int firstZeroRow = 
                ( side==LEFT ? std::max(j-offset+1,0)
                             : std::max(j-offset+height-width+1,0) );
            Int numNonzeroRows = RawLocalLength(firstZeroRow,colShift,c);
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
DistMatrix<T,MR,MC,Int>::ScaleTrapezoid
( T alpha, Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::ScaleTrapezoid");
    this->AssertNotLockedView(); 
#endif
    const elemental::Grid& g = this->Grid();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int r = g.Height();
    const Int c = g.Width();
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
            Int j = rowShift + jLocal*r;
            Int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            Int boundary = std::min( lastRow+1, height );
            Int numRows = RawLocalLength( boundary, colShift, c );
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
            Int j = rowShift + jLocal*r;
            Int firstRow = ( side==LEFT ? std::max(j-offset,0)
                                        : std::max(j-offset+height-width,0) );
            Int numZeroRows = RawLocalLength( firstRow, colShift, c );
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
DistMatrix<T,MR,MC,Int>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const elemental::Grid& g = this->Grid();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int colShift = this->ColShift();
    const Int rowShift = this->RowShift();

    this->SetToZero();

    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const Int i = colShift + iLocal*c;
        if( i % r == rowShift )
        {
            const Int jLocal = (i-rowShift) / r;
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
DistMatrix<T,MR,MC,Int>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            this->SetLocalEntry(i,j,SampleUnitBall<T>());
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
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

        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int myRow = g.MCRank();
        const Int myCol = g.MRRank();
        const Int rankCM = g.VCRank();
        const Int rankRM = g.VRRank();
        const Int ownerRow = this->RowAlignment();
        const Int ownerCol = A.RowAlignment();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int height = A.Height();
        const Int maxLocalHeight = MaxLocalLength(height,p);

        const Int portionSize = std::max(maxLocalHeight,mpi::MIN_COLL_MSG);

        const Int colShiftVR = Shift(rankRM,colAlignment,p);
        const Int colShiftVCOfA = Shift(rankCM,colAlignmentOfA,p);
        const Int sendRankRM = (rankRM+(p+colShiftVCOfA-colShiftVR)) % p;
        const Int recvRankCM = (rankCM+(p+colShiftVR-colShiftVCOfA)) % p;
        const Int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

        this->auxMemory_.Require( (r+c)*portionSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        if( myCol == ownerCol )
        {
            // Pack
            const Int AColShift = A.ColShift();
            const T* ALocalBuffer = A.LockedLocalBuffer();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const Int shift = RawShift(myRow+r*k,colAlignmentOfA,p);
                const Int offset = (shift-AColShift) / r;
                const Int thisLocalHeight = RawLocalLength(height,shift,p);

                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
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
            const Int thisColShift = this->ColShift();
            T* thisLocalBuffer = this->LocalBuffer();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const Int shift = RawShift(myCol+c*k,colAlignment,p);
                const Int offset = (shift-thisColShift) / c;
                const Int thisLocalHeight = RawLocalLength(height,shift,p);

                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    thisLocalBuffer[offset+iLocal*r] = data[iLocal];
            }
        }

        this->auxMemory_.Release();
    }
    else if( A.Height() == 1 )
    {
        if( !this->Viewing() )
            ResizeTo( 1, A.Width() );

        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int myRow = g.MCRank();
        const Int myCol = g.MRRank();
        const Int rankCM = g.VCRank();
        const Int rankRM = g.VRRank();
        const Int ownerCol = this->ColAlignment();
        const Int ownerRow = A.ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int width = A.Width();
        const Int maxLocalWidth = MaxLocalLength(width,p);

        const Int portionSize = std::max(maxLocalWidth,mpi::MIN_COLL_MSG);

        const Int rowShiftVC = Shift(rankCM,rowAlignment,p);
        const Int rowShiftVROfA = Shift(rankRM,rowAlignmentOfA,p);
        const Int sendRankCM = (rankCM+(p+rowShiftVROfA-rowShiftVC)) % p;
        const Int recvRankRM = (rankRM+(p+rowShiftVC-rowShiftVROfA)) % p;
        const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        this->auxMemory_.Require( (r+c)*portionSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        if( myRow == ownerRow )
        {
            // Pack
            const Int ARowShift = A.RowShift();
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const Int shift = RawShift(myCol+c*k,rowAlignmentOfA,p);
                const Int offset = (shift-ARowShift) / c;
                const Int thisLocalWidth = RawLocalLength(width,shift,p);

                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
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
            const Int thisRowShift = this->RowShift();
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const Int shift = RawShift(myRow+r*k,rowAlignment,p);
                const Int offset = (shift-thisRowShift) / r;
                const Int thisLocalWidth = RawLocalLength(width,shift,p);

                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                    thisLocalBuffer[(offset+jLocal*c)*thisLDim] = data[jLocal];
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
        if( A.Height() >= A.Width() )
        {
            std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
            ( new DistMatrix<T,VC,STAR,Int>(g) );
            *A_VC_STAR = A;

            std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
            ( new DistMatrix<T,VR,STAR,Int>(true,this->ColAlignment(),g) );
            *A_VR_STAR = *A_VC_STAR;
            delete A_VC_STAR.release(); // lowers memory highwater

            *this = *A_VR_STAR;
        }
        else
        {
            std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
            ( new DistMatrix<T,STAR,VR,Int>(g) );
            *A_STAR_VR = A;

            std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
            ( new DistMatrix<T,STAR,VC,Int>(true,this->RowAlignment(),g) );
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

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(g) );
    *A_VC_STAR = A;

    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VR_STAR = *A_VC_STAR;
    delete A_VC_STAR.release(); // lowers memory highwater

    *this = *A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR,Int>(g) );
    *A_STAR_VR = A;

    std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC,Int>(true,this->RowAlignment(),g) );
    *A_STAR_VC = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MR,MC] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MR,MC] = [* ,MD] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
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
            std::cerr << "Unaligned [MR,MC] <- [MR,MC]" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.MCRank();
        const Int col = g.MRRank();

        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const Int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
        const Int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;
        const Int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;
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
            std::memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
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
            std::memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
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
        const Int r = g.Height();
        const Int rowShift = this->RowShift();

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();

        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[(rowShift+jLocal*r)*ALDim];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            std::memcpy( thisCol, ACol, localHeight*sizeof(T) );
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [MR,MC] <- [MR,* ]" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int col = g.MRRank();

        const Int rowShift = this->RowShift();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int sendRank = (col+c+colAlignment-colAlignmentOfA) % c;
        const Int recvRank = (col+c+colAlignmentOfA-colAlignment) % c;

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
            const T* ACol = &ALocalBuffer[(rowShift+jLocal*r)*ALDim];
            T* sendBufferCol = &sendBuffer[jLocal*localHeightOfA];
            std::memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.MRComm() );

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
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
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
        const Int c = g.Width();
        const Int colShift = this->ColShift();

        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();

        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                thisLocalBuffer[iLocal+jLocal*thisLDim] = 
                    ALocalBuffer[(colShift+iLocal*c)+jLocal*ALDim];
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [MR,MC] <- [* ,MC]" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.MCRank(); 

        const Int colShift = this->ColShift();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

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
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                sendBuffer[iLocal+jLocal*localHeight] = 
                    ALocalBuffer[(colShift+iLocal*c)+jLocal];

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendRow, 0,
          recvBuffer, recvSize, recvRow, mpi::ANY_TAG, g.MCComm() );

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
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,VR,STAR,Int> A_VR_STAR(g);

    A_VR_STAR = A;
    *this = A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
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
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int row = g.MCRank();

        const Int rowShift = this->RowShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLocalLength(height,c);
        const Int maxWidth = MaxLocalLength(width,p);
        const Int portionSize = std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

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

            const Int thisColShift = RawShift(k,colAlignment,c);
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+jLocal*ALDim];
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
            const Int thisRowShift = RawShift(thisRank,rowAlignmentOfA,p);
            const Int thisRowOffset = (thisRowShift-rowShift) / r;
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = 
                    &thisLocalBuffer[(thisRowOffset+jLocal*c)*thisLDim];
                std::memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [MR,MC] <- [* ,VC]" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int row = g.MCRank();

        const Int rowShift = this->RowShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRow = (row+r+rowAlignment-(rowAlignmentOfA%r)) % r;
        const Int recvRow = (row+r+(rowAlignmentOfA%r)-rowAlignment) % r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLocalLength(height,c);
        const Int maxWidth = MaxLocalLength(width,p);
        const Int portionSize = std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

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

            const Int thisColShift = RawShift(k,colAlignment,c);
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
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
        const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const Int thisRank = recvRow+k*r;
            const Int thisRowShift = RawShift(thisRank,rowAlignmentOfA,p);
            const Int thisRowOffset = (thisRowShift-rowShift) / r;
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = 
                    &thisLocalBuffer[(thisRowOffset+jLocal*c)*thisLDim];
                std::memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
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
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int col = g.MRRank();

        const Int colShift = this->ColShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int maxHeight = MaxLocalLength(height,p);
        const Int maxWidth = MaxLocalLength(width,r);
        const Int portionSize = std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

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

            const Int thisRowShift = RawShift(k,rowAlignment,r);
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                std::memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
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
            const Int thisColShift = RawShift(thisRank,colAlignmentOfA,p);
            const Int thisColOffset = (thisColShift-colShift) / c;
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,p);
            
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    thisLocalBuffer[(thisColOffset+iLocal*r)+jLocal*thisLDim] =
                        data[iLocal+jLocal*thisLocalHeight];
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [MR,MC] <- [* ,VC]" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int col = g.MRRank();

        const Int colShift = this->ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int sendCol = (col+c+colAlignment-(colAlignmentOfA%c)) % c;
        const Int recvCol = (col+c+(colAlignmentOfA%c)-colAlignment) % c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int maxHeight = MaxLocalLength(height,p);
        const Int maxWidth = MaxLocalLength(width,r);
        const Int portionSize = std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

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

            const Int thisRowShift = RawShift(k,rowAlignment,r);
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                std::memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
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
        const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const Int thisRank = recvCol+k*c;
            const Int thisColShift = RawShift(thisRank,colAlignmentOfA,p);
            const Int thisColOffset = (thisColShift-colShift) / c;
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
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

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,STAR,VC,Int> A_STAR_VC(g);

    A_STAR_VC = A;
    *this = A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MR,MC,Int>&
DistMatrix<T,MR,MC,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
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
    const Int r = g.Height();
    const Int c = g.Width();
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
                ALocalBuffer[(colShift+iLocal*c)+(rowShift+jLocal*r)*ALDim];
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::SumScatterFrom( const DistMatrix<T,MR,STAR,Int>& A )
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
            const Int rowAlignment = this->RowAlignment();
            const Int myRow = g.MCRank();

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
              mpi::SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
            {
                T* thisCol = this->LocalBuffer(0,0);
                std::memcpy( thisCol, recvBuffer, localHeight*sizeof(T) );
            }

            this->auxMemory_.Release();
        }
        else
        {
            const Int r = g.Height();
            const Int myRow = g.MCRank();
            const Int rowAlignment = this->RowAlignment();

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int maxLocalWidth = MaxLocalLength(width,r);

            const Int recvSize = 
                std::max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
            const Int sendSize = r * recvSize;

            this->auxMemory_.Require( sendSize );
            T* sendBuffer = this->auxMemory_.Buffer();

            // Pack 
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &sendBuffer[k*recvSize];

                const Int thisRowShift = RawShift( k, rowAlignment, r );
                const Int thisLocalWidth = 
                      RawLocalLength( width, thisRowShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
                    T* dataCol = &data[jLocal*localHeight];
                    std::memcpy( dataCol, ACol, localHeight*sizeof(T) );
                }
            }

            // AllReduce over each process column
            mpi::AllReduce( sendBuffer, sendSize, mpi::SUM, g.MCComm() );

            // Unpack our received data
            T* thisLocalBuffer = this->LocalBuffer();
            const T* recvBuffer = &sendBuffer[myRow*recvSize];
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int j=0; j<localWidth; ++j )
            {
                const T* recvBufferCol = &recvBuffer[j*localHeight];
                T* thisCol = &thisLocalBuffer[j*thisLDim];
                std::memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned SumScatterFrom [MR,MC] <- [MR,* ]" 
                      << std::endl;
#endif
        if( this->Width() == 1 )
        {
            const Int c = g.Width();
            const Int rowAlignment = this->RowAlignment();
            const Int myRow = g.MCRank();
            const Int myCol = g.MRRank();

            const Int height = this->Height();
            const Int localHeight = this->LocalHeight();
            const Int localHeightOfA = A.LocalHeight();
            const Int maxLocalHeight = MaxLocalLength(height,c);

            const Int portionSize = std::max(maxLocalHeight,mpi::MIN_COLL_MSG);

            const Int colAlignment = this->ColAlignment();
            const Int colAlignmentOfA = A.ColAlignment();
            const Int sendCol = (myCol+c+colAlignment-colAlignmentOfA) % c;
            const Int recvCol = (myCol+c+colAlignmentOfA-colAlignment) % c;

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
              mpi::SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
            {
                // Perform the realignment
                mpi::SendRecv
                ( recvBuffer, portionSize, sendCol, 0,
                  sendBuffer, portionSize, recvCol, 0, g.MRComm() );

                T* thisCol = this->LocalBuffer(0,0);
                std::memcpy( thisCol, sendBuffer, localHeight*sizeof(T) );
            }

            this->auxMemory_.Release();
        }
        else
        {
            const Int r = g.Height();
            const Int c = g.Width();
            const Int col = g.MRRank();

            const Int colAlignment = this->ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int colAlignmentOfA = A.ColAlignment();
            const Int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
            const Int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localHeightOfA = A.LocalHeight();
            const Int maxLocalWidth = MaxLocalLength(width,r);

            const Int recvSize_RS = 
                std::max(localHeightOfA*maxLocalWidth,mpi::MIN_COLL_MSG);
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

                const Int thisRowShift = RawShift( k, rowAlignment, r );
                const Int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
                    T* dataCol = &data[jLocal*localHeightOfA];
                    std::memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                }
            }

            // Reduce-scatter over each process col
            mpi::ReduceScatter
            ( secondBuffer, firstBuffer, recvSize_RS, mpi::SUM, g.MCComm() );

            // Trade reduced data with the appropriate process col
            mpi::SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendCol, 0,
              secondBuffer, localHeight*localWidth,    recvCol, mpi::ANY_TAG,
              g.MRComm() );

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
DistMatrix<T,MR,MC,Int>::SumScatterFrom( const DistMatrix<T,STAR,MC,Int>& A )
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
        std::cerr <<    
          "The vector version of [MR,MC].SumScatterFrom([* ,MC]) is not yet"
          " written, but it only requires a modification of the vector "
          "version of [MR,MC].SumScatterFrom([MR,* ])" << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.VCRank() == 0 )
    {
        std::cerr << 
          "[MR,MC]::SumScatterFrom([* ,MC]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MC] matrix instead." << std::endl;
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
        const Int c = g.Width();
        const Int myCol = g.MRRank();
        const Int colAlignment = this->ColAlignment();

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalHeight = MaxLocalLength(height,c);

        const Int recvSize = 
            std::max(maxLocalHeight*localWidth,mpi::MIN_COLL_MSG);
        const Int sendSize = c * recvSize;

        this->auxMemory_.Require( sendSize );
        T* sendBuffer = this->auxMemory_.Buffer();

        // Pack 
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*recvSize];

            const Int thisColShift = RawShift( k, colAlignment, c );
            const Int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+jLocal*ALDim];
        }

        // AllReduce over each process row
        mpi::AllReduce( sendBuffer, sendSize, mpi::SUM, g.MRComm() );

        // Unpack our received data
        T* thisLocalBuffer = this->LocalBuffer();
        const T* recvBuffer = &sendBuffer[myCol*recvSize];
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
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned SumScatterFrom [MR,MC] <- [* ,MC]" 
                      << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.MCRank();

        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalHeight = MaxLocalLength(height,c);
        
        const Int recvSize_RS = 
            std::max(maxLocalHeight*localWidthOfA,mpi::MIN_COLL_MSG);
        const Int sendSize_RS = c* recvSize_RS;
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

            const Int thisColShift = RawShift( k, colAlignment, c );
            const Int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+jLocal*ALDim];
        }

        // Reduce-scatter over each process row
        mpi::ReduceScatter
        ( secondBuffer, firstBuffer, recvSize_RS, mpi::SUM, g.MRComm() );

        // Trade reduced data with the appropriate process row
        mpi::SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendRow, 0,
          secondBuffer, localHeight*localWidth,    recvRow, mpi::ANY_TAG, 
          g.MCComm() );

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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::SumScatterFrom( const DistMatrix<T,STAR,STAR,Int>& A )
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
    const Int r = g.Height();
    const Int c = g.Width();
    const Int myRank = g.VRRank();
    const Int colAlignment = this->ColAlignment();
    const Int rowAlignment = this->RowAlignment();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int maxLocalHeight = MaxLocalLength(height,c);
    const Int maxLocalWidth = MaxLocalLength(width,r);

    const Int recvSize = 
        std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
    const Int sendSize = r * c * recvSize;

    this->auxMemory_.Require( sendSize );
    T* sendBuffer = this->auxMemory_.Buffer();

    // Pack 
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int l=0; l<r; ++l )
    {
        const Int thisRowShift = RawShift( l, rowAlignment, r );
        const Int thisLocalWidth = RawLocalLength( width, thisRowShift, r );

        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[(k+l*c)*recvSize];

            const Int thisColShift = RawShift( k, colAlignment, c );
            const Int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+
                                     (thisRowShift+jLocal*r)*ALDim];
        }
    }

    // AllReduce over the grid
    mpi::AllReduce( sendBuffer, sendSize, mpi::SUM, g.VRComm() );

    // Unpack our received data
    T* thisLocalBuffer = this->LocalBuffer();
    const T* recvBuffer = &sendBuffer[myRank*recvSize];
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,MR,STAR,Int>& A )
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
            const Int rowAlignment = this->RowAlignment();
            const Int myRow = g.MCRank();

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
              mpi::SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
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
            const Int r = g.Height();
            const Int myRow = g.MCRank();
            const Int rowAlignment = this->RowAlignment();

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int maxLocalWidth = MaxLocalLength(width,r);

            const Int portionSize = 
                std::max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
            const Int sendSize = r*portionSize;

            this->auxMemory_.Require( sendSize );
            T* sendBuffer = this->auxMemory_.Buffer();

            // Pack 
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<r; ++k )
            {
                T* data = &sendBuffer[k*portionSize];

                const Int thisRowShift = RawShift( k, rowAlignment, r );
                const Int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
                    T* dataCol = &data[jLocal*localHeight];
                    std::memcpy( dataCol, ACol, localHeight*sizeof(T) );
                }
            }

            // AllReduce over each process column
            mpi::AllReduce( sendBuffer, sendSize, mpi::SUM, g.MCComm() );

            // Update with our received data
            T* thisLocalBuffer = this->LocalBuffer();
            const T* recvBuffer = &sendBuffer[myRow*portionSize];
            const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisCol[iLocal] += alpha*recvBufferCol[iLocal];
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned SumScatterUpdate [MR,MC] <- [MR,* ]" 
                      << std::endl;
#endif
        if( this->Width() == 1 )
        {
            const Int c = g.Width();
            const Int rowAlignment = this->RowAlignment();
            const Int myRow = g.MCRank();
            const Int myCol = g.MRRank();

            const Int height = this->Height();
            const Int localHeight = this->LocalHeight();
            const Int localHeightOfA = A.LocalHeight();
            const Int maxLocalHeight = MaxLocalLength(height,c);

            const Int portionSize = std::max(maxLocalHeight,mpi::MIN_COLL_MSG);

            const Int colAlignment = this->ColAlignment();
            const Int colAlignmentOfA = A.ColAlignment();
            const Int sendCol = (myCol+c+colAlignment-colAlignmentOfA) % c;
            const Int recvCol = (myCol+c+colAlignmentOfA-colAlignment) % c;

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
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisCol[iLocal] += alpha*sendBuffer[iLocal];
            }
            this->auxMemory_.Release();
        }
        else
        {
            const Int r = g.Height();
            const Int c = g.Width();
            const Int col = g.MRRank();

            const Int colAlignment = this->ColAlignment();
            const Int rowAlignment = this->RowAlignment();
            const Int colAlignmentOfA = A.ColAlignment();
            const Int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
            const Int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidth = this->LocalWidth();
            const Int localHeightOfA = A.LocalHeight();
            const Int maxLocalWidth = MaxLocalLength(width,r);

            const Int recvSize_RS = 
                std::max(localHeightOfA*maxLocalWidth,mpi::MIN_COLL_MSG);
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

                const Int thisRowShift = RawShift( k, rowAlignment, r );
                const Int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                {
                    const T* ACol = 
                        &ALocalBuffer[(thisRowShift+jLocal*r)*ALDim];
                    T* dataCol = &data[jLocal*localHeightOfA];
                    std::memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                }
            }

            // Reduce-scatter over each process col
            mpi::ReduceScatter
            ( secondBuffer, firstBuffer, recvSize_RS, mpi::SUM, g.MCComm() );

            // Trade reduced data with the appropriate process col
            mpi::SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendCol, 0, 
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
DistMatrix<T,MR,MC,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MC,Int>& A )
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
        std::cerr <<    
          "The vector version of [MR,MC].SumScatterUpdate([* ,MC]) is not "
          "yet written, but it only requires a modification of the vector "
          "version of [MR,MC].SumScatterUpdate([MR,* ])." << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.VCRank() == 0 )
    {
        std::cerr <<
          "[MR,MC]::SumScatterUpdate([* ,MC]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MC] matrix instead." << std::endl;
    }
#endif
    if( this->RowAlignment() == A.RowAlignment() )
    {
        const Int c = g.Width();
        const Int myCol = g.MRRank();
        const Int colAlignment = this->ColAlignment();

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalHeight = MaxLocalLength(height,c);

        const Int recvSize = 
            std::max(maxLocalHeight*localWidth,mpi::MIN_COLL_MSG);
        const Int sendSize = c * recvSize;

        this->auxMemory_.Require( sendSize );
        T* sendBuffer = this->auxMemory_.Buffer();

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*recvSize];

            const Int thisColShift = RawShift( k, colAlignment, c );
            const Int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+jLocal*ALDim];
        }

        // AllReduce over each process row
        mpi::AllReduce( sendBuffer, sendSize, mpi::SUM, g.MRComm() );

        // Update with our received data
        T* thisLocalBuffer = this->LocalBuffer();
        const T* recvBuffer = &sendBuffer[myCol*recvSize];
        const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(AVOID_OMP_FMA)
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufferCol = &recvBuffer[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                thisCol[iLocal] += alpha*recvBufferCol[iLocal];
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned SumScatterUpdate [MR,MC] <- [* ,MC]" 
                      << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.MCRank();

        const Int colAlignment = this->ColAlignment();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const Int height = this->Height();
        const Int localHeight = this->LocalHeight();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalHeight = MaxLocalLength(height,c);

        const Int recvSize_RS = 
            std::max(maxLocalHeight*localWidthOfA,mpi::MIN_COLL_MSG);
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

            const Int thisColShift = RawShift( k, colAlignment, c );
            const Int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c)+jLocal*ALDim];
        }

        // Reduce-scatter over each process row
        mpi::ReduceScatter
        ( secondBuffer, firstBuffer, recvSize_RS, mpi::SUM, g.MRComm() );

        // Trade reduced data with the appropriate process row
        mpi::SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendRow, 0,
          secondBuffer, localHeight*localWidth,    recvRow, mpi::ANY_TAG,
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MR,MC,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR,Int>& A )
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
    const Int r = g.Height();
    const Int c = g.Width();
    const Int myRank = g.VRRank();
    const Int colAlignment = this->ColAlignment();
    const Int rowAlignment = this->RowAlignment();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int maxLocalHeight = MaxLocalLength(height,c);
    const Int maxLocalWidth = MaxLocalLength(width,r);

    const Int recvSize = 
        std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);
    const Int sendSize = r * c * recvSize;

    this->auxMemory_.Require( sendSize );
    T* sendBuffer = this->auxMemory_.Buffer();

    // Pack 
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int l=0; l<r; ++l )
    {
        const Int thisRowShift = RawShift( l, rowAlignment, r );
        const Int thisLocalWidth = RawLocalLength( width, thisRowShift, r );

        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[(k+l*c)*recvSize];

            const Int thisColShift = RawShift( k, colAlignment, c );
            const Int thisLocalHeight = 
                RawLocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColShift+iLocal*c) + 
                                     (thisRowShift+jLocal*r)*ALDim];
        }
    }

    // AllReduce over the grid
    mpi::AllReduce( sendBuffer, sendSize, mpi::SUM, g.VRComm() );

    // Unpack our received data
    T* thisLocalBuffer = this->LocalBuffer();
    const T* recvBuffer = &sendBuffer[myRank*recvSize];
    const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
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
