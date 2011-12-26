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
DistMatrix<T,VR,STAR,Int>::DistMatrix( const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   (g.InGrid() ? g.VRRank() : 0 ),0,
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VR,STAR,Int>::DistMatrix
( Int height, Int width, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   (g.InGrid() ? g.VRRank() : 0),0,
   (g.InGrid() ? LocalLength(height,g.VRRank(),0,g.Size()) : 0),width,
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VR,STAR,Int>::DistMatrix
( bool constrainedColAlignment, Int colAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.VRRank(),colAlignment,g.Size()) : 0),0,
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VR,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.VRRank(),colAlignment,g.Size()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.VRRank(),colAlignment,g.Size()) : 0),
   width,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VR,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
  Int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.VRRank(),colAlignment,g.Size()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.VRRank(),colAlignment,g.Size()) : 0),
   width,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VR,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, const T* buffer, Int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,colAlignment,0,
   (g.InGrid() ? Shift(g.VRRank(),colAlignment,g.Size()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.VRRank(),colAlignment,g.Size()) : 0),
   width,buffer,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VR,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, T* buffer, Int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,colAlignment,0,
   (g.InGrid() ? Shift(g.VRRank(),colAlignment,g.Size()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.VRRank(),colAlignment,g.Size()) : 0),
   width,buffer,ldim,g)
{ }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,VR,STAR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix");
#endif
    if( VR != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,VR,STAR,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [VR,* ] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,VR,STAR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::SetGrid( const elemental::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->colAlignment_ = 0;
    this->colShift_ = grid.VRRank();
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VR,STAR,Int>::AlignWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ = Shift( g.VRRank(), this->ColAlignment(), g.Size() );
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
DistMatrix<T,VR,STAR,Int>::AlignWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ = Shift( g.VRRank(), this->ColAlignment(), g.Size() );
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
DistMatrix<T,VR,STAR,Int>::AlignWith( const DistMatrix<S,MR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([MR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ = Shift( g.VRRank(), this->ColAlignment(), g.Size() );
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
DistMatrix<T,VR,STAR,Int>::AlignWith( const DistMatrix<S,STAR,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([* ,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ = Shift( g.VRRank(), this->ColAlignment(), g.Size() );
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
DistMatrix<T,VR,STAR,Int>::AlignWith( const DistMatrix<S,VR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith(DistMatrix[VR,* ])");
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
DistMatrix<T,VR,STAR,Int>::AlignWith( const DistMatrix<S,STAR,VR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith(DistMatrix[* ,VR])");
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
DistMatrix<T,VR,STAR,Int>::AlignColsWith( const DistMatrix<S,MC,MR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VR,STAR,Int>::AlignColsWith( const DistMatrix<S,MR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VR,STAR,Int>::AlignColsWith( const DistMatrix<S,MR,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VR,STAR,Int>::AlignColsWith( const DistMatrix<S,STAR,MR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VR,STAR,Int>::AlignColsWith( const DistMatrix<S,VR,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VR,STAR,Int>::AlignColsWith( const DistMatrix<S,STAR,VR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VRRank() == 0 && msg != "" )
        os << msg << std::endl;

    const Int height      = this->Height();
    const Int width       = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int p           = g.Size();
    const Int colShift    = this->ColShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    std::vector<T> sendBuf(height*width,0);
    const T* thisLocalBuffer = this->LockedLocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
        for( Int j=0; j<width; ++j )
            sendBuf[colShift+iLocal*p+j*height] = 
                thisLocalBuffer[iLocal+j*thisLDim];

    // If we are the root, allocate a receive buffer
    std::vector<T> recvBuf;
    if( g.VRRank() == 0 )
        recvBuf.resize( height*width );

    // Sum the contributions and send to the root
    mpi::Reduce
    ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, g.VCComm() );

    if( g.VRRank() == 0 )
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
DistMatrix<T,VR,STAR,Int>::Align( Int colAlignment )
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

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::AlignCols( Int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Size() )
        throw std::runtime_error("Invalid column alignment for [VR,* ]");
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

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::View( DistMatrix<T,VR,STAR,Int>& A )
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

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::View
( Int height, Int width, Int colAlignment,
  T* buffer, Int ldim, const elemental::Grid& grid )
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
    const Int localHeight = LocalLength(height,this->colShift_,grid.Size());
    this->localMatrix_.View( localHeight, width, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::LockedView( const DistMatrix<T,VR,STAR,Int>& A )
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

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::LockedView
( Int height, Int width, Int colAlignment,
  const T* buffer, Int ldim, const elemental::Grid& grid )
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
    const Int localHeight = LocalLength(height,this->colShift_,grid.Size());
    this->localMatrix_.LockedView( localHeight, width, buffer, ldim );
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::View
( DistMatrix<T,VR,STAR,Int>& A, Int i, Int j, Int height, Int width )
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
        const Int rowMajorRank = g.VRRank();
        const Int size = g.Size();

        this->colAlignment_ = (A.ColAlignment()+i) % size;
        this->colShift_ = Shift( rowMajorRank, this->ColAlignment(), size );

        const Int localHeightBefore = LocalLength( i, A.ColShift(), size );
        const Int localHeight = LocalLength( height, this->ColShift(), size );

        this->localMatrix_.View
        ( A.LocalMatrix(), localHeightBefore, j, localHeight, width );
    }
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::LockedView
( const DistMatrix<T,VR,STAR,Int>& A, Int i, Int j, Int height, Int width )
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
        const Int rowMajorRank = g.VRRank();
        const Int size = g.Size();

        this->colAlignment_ = (A.ColAlignment()+i) % size;
        this->colShift_ = Shift( rowMajorRank, this->ColAlignment(), size );

        const Int localHeightBefore = LocalLength( i, A.ColShift(), size );
        const Int localHeight = LocalLength( height, this->ColShift(), size );

        this->localMatrix_.LockedView
        ( A.LockedLocalMatrix(), localHeightBefore, j, localHeight, width );
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::View1x2
( DistMatrix<T,VR,STAR,Int>& AL, DistMatrix<T,VR,STAR,Int>& AR )
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

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::LockedView1x2
( const DistMatrix<T,VR,STAR,Int>& AL, const DistMatrix<T,VR,STAR,Int>& AR )
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

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::View2x1
( DistMatrix<T,VR,STAR,Int>& AT,
  DistMatrix<T,VR,STAR,Int>& AB )
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

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::LockedView2x1
( const DistMatrix<T,VR,STAR,Int>& AT,
  const DistMatrix<T,VR,STAR,Int>& AB )
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

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::View2x2
( DistMatrix<T,VR,STAR,Int>& ATL, DistMatrix<T,VR,STAR,Int>& ATR,
  DistMatrix<T,VR,STAR,Int>& ABL, DistMatrix<T,VR,STAR,Int>& ABR )
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

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::LockedView2x2
( const DistMatrix<T,VR,STAR,Int>& ATL, const DistMatrix<T,VR,STAR,Int>& ATR,
  const DistMatrix<T,VR,STAR,Int>& ABL, const DistMatrix<T,VR,STAR,Int>& ABR )
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

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::ResizeTo");
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

template<typename T,typename Int>
inline T
DistMatrix<T,VR,STAR,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    T u;
    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        u = this->GetLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->SetLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->UpdateLocalEntry(iLoc,j,u);
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
DistMatrix<T,VR,STAR,Int>::MakeTrapezoidal
( Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int p = this->Grid().Size();
    const Int colShift = this->ColShift();

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
                Int numZeroRows = RawLocalLength( boundary, colShift, p );
                T* thisCol = &thisLocalBuffer[j*thisLDim];
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
        for( Int j=0; j<width; ++j )
        {
            Int firstZeroRow = 
                ( side==LEFT ? std::max(j-offset+1,0)
                             : std::max(j-offset+height-width+1,0) );
            Int numNonzeroRows = RawLocalLength(firstZeroRow,colShift,p);
            if( numNonzeroRows < localHeight )
            {
                T* thisCol = &thisLocalBuffer[numNonzeroRows+j*thisLDim];
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
DistMatrix<T,VR,STAR,Int>::ScaleTrapezoid
( T alpha, Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::ScaleTrapezoid");
    this->AssertNotLockedView();
#endif
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int p = this->Grid().Size();
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
            Int numRows = RawLocalLength( boundary, colShift, p );
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
                                        : std::max(j-offset+height-width,0) );
            Int numZeroRows = RawLocalLength(firstRow,colShift,p);
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
DistMatrix<T,VR,STAR,Int>::SetToIdentity()
{   
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int p = this->Grid().Size();
    const Int colShift = this->ColShift();

    this->SetToZero();

    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const Int i = colShift + iLocal*p;
        if( i < width )
            thisLocalBuffer[iLocal+i*thisLDim] = 1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    for( Int j=0; j<width; ++j )
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
            this->SetLocalEntry(iLocal,j,SampleUnitBall<T>());
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR,Int> A_VC_STAR(g);

    A_VC_STAR = A;
    *this = A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR,Int> A_VC_STAR(g);

    A_VC_STAR = A;
    *this = A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MC,MR,Int> > A_MC_MR
    ( new DistMatrix<T,MC,MR,Int>(g) );
    *A_MC_MR = A;

    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(g) );
    *A_VC_STAR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[VR,* ] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[VR,* ] = [* ,MD] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
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
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.MRRank();
        const Int colShiftOfA = A.ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

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

            const Int thisRank = col+k*c;
            const Int thisColShift = RawShift(thisRank,colAlignment,p);
            const Int thisColOffset = (thisColShift-colShiftOfA) / c;
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+jLocal*thisLocalHeight] = 
                        ALocalBuffer[(thisColOffset+iLocal*r)+jLocal*ALDim];
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

            const Int thisRowShift = RawShift(k,rowAlignmentOfA,r);
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(thisRowShift+jLocal*r)*thisLDim];
                std::memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [VR,* ] <- [MR,MC]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.MRRank();
        const Int colShiftOfA = A.ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendCol = (col+c+(colAlignment%c)-colAlignmentOfA) % c;
        const Int recvCol = (col+c+colAlignmentOfA-(colAlignment%c)) % c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

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

            const Int thisRank = sendCol+k*c;
            const Int thisColShift = RawShift(thisRank,colAlignment,p);
            const Int thisColOffset = (thisColShift-colShiftOfA) / c;
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
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
        const Int thisLDim = this->LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const Int thisRowShift = RawShift(k,rowAlignmentOfA,r);
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(thisRowShift+jLocal*r)*thisLDim];
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
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
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
        const Int r = g.Height();
        const Int c = g.Width();
        const Int colShift = this->ColShift();
        const Int colShiftOfA = A.ColShift();
        const Int colOffset = (colShift-colShiftOfA) / c;

        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();

        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( Int j=0; j<width; ++j )
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                thisLocalBuffer[iLocal+j*thisLDim] = 
                    ALocalBuffer[(colOffset+iLocal*r)+j*ALDim];
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            std::cerr << "Unaligned [VR,* ] <- [MR,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.MCRank();
        const Int col = g.MRRank();
        const Int colShiftOfA = A.ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[VR,*] within our process row to fix alignments.
        const Int sendCol = (col+c+(colAlignment%c)-colAlignmentOfA) % c;
        const Int recvCol = (col+c+colAlignmentOfA-(colAlignment%c)) % c;
        const Int sendRank = sendCol + c*row;

        const Int sendColShift = Shift( sendRank, colAlignment, p );
        const Int sendColOffset = (sendColShift-colShiftOfA) / c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localHeightOfSend = LocalLength(height,sendColShift,p);
        const Int maxLocalHeight = MaxLocalLength(height,p);

        const Int portionSize = maxLocalHeight * width;

        this->auxMemory_.Require( 2*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( Int j=0; j<width; ++j )
            for( Int iLocal=0; iLocal<localHeightOfSend; ++iLocal )
                sendBuffer[iLocal+j*localHeightOfSend] = 
                    ALocalBuffer[(sendColOffset+iLocal*r)+j*ALDim];

        // Communicate
        mpi::SendRecv
        ( sendBuffer, portionSize, sendCol, 0,
          recvBuffer, portionSize, recvCol, mpi::ANY_TAG, g.MRComm() );

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
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MR,MC,Int> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
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
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int rankCM = g.VCRank();
    const Int rankRM = g.VRRank();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localHeightOfA = A.LocalHeight();
    const Int maxLocalHeight = MaxLocalLength(height,p);

    const Int portionSize = maxLocalHeight * width;

    const Int colShift = this->ColShift();
    const Int colShiftOfA = A.ColShift();

    // Compute which rowmajor rank has the colShift equal to our colShiftOfA
    const Int sendRankRM = (rankRM+(p+colShiftOfA-colShift)) % p;

    // Compute which rowmajor rank has the A colShift that we need
    const Int recvRankCM = (rankCM+(p+colShift-colShiftOfA)) % p;
    const Int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

    this->auxMemory_.Require( 2*portionSize );

    T* buffer = this->auxMemory_.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[portionSize];

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
        std::memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
    }

    // Communicate
    mpi::SendRecv
    ( sendBuffer, portionSize, sendRankRM, 0,
      recvBuffer, portionSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

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
        std::memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    DistMatrix<T,MR,MC,Int> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
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
            std::cerr << "Unaligned [VR,* ] <- [VR,* ]." << std::endl;
#endif
        const Int rank = g.VRRank();
        const Int p = g.Size();

        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const Int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

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
            std::memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.VRComm() );

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
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MC,MR,Int> > A_MC_MR
    ( new DistMatrix<T,MC,MR,Int>(g) );
    *A_MC_MR = A;

    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(g) );
    *A_VC_STAR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VR,STAR,Int>&
DistMatrix<T,VR,STAR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
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

    const Int p = this->Grid().Size();
    const Int colShift = this->ColShift();

    const Int localHeight = this->LocalHeight();
    const Int width = this->Width();

    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( Int j=0; j<width; ++j )
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
            thisLocalBuffer[iLocal+j*thisLDim] = 
                ALocalBuffer[(colShift+iLocal*p)+j*ALDim];
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::SumScatterFrom
( const DistMatrix<T,MR,STAR,Int>& A )
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
        std::cerr <<
          "[VR,* ]::SumScatterFrom([MR,* ]) potentially causes a large amount "
          "of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MR,* ] matrix instead." << std::endl;
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
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.MCRank();
        const Int myCol = g.MRRank();
        const Int colAlignment = this->ColAlignment();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int maxLocalHeight = MaxLocalLength( height, p );

        const Int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
        const Int sendSize = r*recvSize;

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

            const Int thisRank = myCol+k*c;
            const Int thisColShift = RawShift( thisRank, colAlignment, p );
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = 
                RawLocalLength( height, thisColShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int j=0; j<width; ++j )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+j*thisLocalHeight] = 
                        ALocalBuffer[(thisColOffset+iLocal*r)+j*ALDim];
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
        for( Int j=0; j<width; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*localHeight];
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            std::memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
        this->auxMemory_.Release();
    }
    else
    {
        throw std::logic_error
        ("Unaligned [VR,* ]::ReduceScatterFrom( [MR,* ] ) not implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::SumScatterFrom
( const DistMatrix<T,STAR,STAR,Int>& A )
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

    const Int p = g.Size();
    const Int VRRank = g.VRRank();
    const Int colAlignment = this->ColAlignment();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int maxLocalHeight = MaxLocalLength( height, p );

    const Int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
    const Int sendSize = p*recvSize;

    this->auxMemory_.Require( sendSize );
    T* sendBuffer = this->auxMemory_.Buffer();

    // Pack
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        T* data = &sendBuffer[k*recvSize];

        const Int thisColShift = RawShift( k, colAlignment, p );
        const Int thisLocalHeight = RawLocalLength( height, thisColShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( Int j=0; j<width; ++j )
            for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                data[iLocal+j*thisLocalHeight] =
                    ALocalBuffer[(thisColShift+iLocal*p)+j*ALDim];
    }

    // AllReduce over the grid
    mpi::AllReduce( sendBuffer, sendSize, mpi::SUM, g.VRComm() );

    // Unpack our received data
    T* thisLocalBuffer = this->LocalBuffer();
    const T* recvBuffer = &sendBuffer[VRRank*recvSize];
    const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
    {
        const T* recvBufferCol = &recvBuffer[j*localHeight];
        T* thisCol = &thisLocalBuffer[j*thisLDim];
        std::memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,MR,STAR,Int>& A )
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
        std::cerr <<
          "[VR,* ]::SumScatterUpdate([MR,* ]) potentially causes a large amount"
          " of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MR,* ] matrix instead." << std::endl;
    }
#endif
    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.MCRank();
        const Int myCol = g.MRRank();
        const Int colAlignment = this->ColAlignment();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int maxLocalHeight = MaxLocalLength( height, p );

        const Int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
        const Int sendSize = r*recvSize;

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

            const Int thisRank = myCol+k*c;
            const Int thisColShift = RawShift( thisRank, colAlignment, p );
            const Int thisColOffset = (thisColShift-colShiftOfA) / c;
            const Int thisLocalHeight = 
                RawLocalLength( height, thisColShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int j=0; j<width; ++j )
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    data[iLocal+j*thisLocalHeight] = 
                        ALocalBuffer[(thisColOffset+iLocal*r)+j*ALDim];
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
        for( Int j=0; j<width; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*localHeight];
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                thisCol[iLocal] += alpha*recvBufferCol[iLocal];
        }
        this->auxMemory_.Release();
    }
    else
    {
        throw std::logic_error
        ("Unaligned [VR,* ]::ReduceScatterUpdate( [MR,* ] ) not implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SumScatterUpdate( [* ,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const elemental::Grid& g = this->Grid();

    const Int p = g.Size();
    const Int VRRank = g.VRRank();
    const Int colAlignment = this->ColAlignment();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int maxLocalHeight = MaxLocalLength( height, p );

    const Int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
    const Int sendSize = p*recvSize;

    this->auxMemory_.Require( sendSize );
    T* sendBuffer = this->auxMemory_.Buffer();

    // Pack
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        T* data = &sendBuffer[k*recvSize];

        const Int thisColShift = RawShift( k, colAlignment, p );
        const Int thisLocalHeight = RawLocalLength( height, thisColShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( Int j=0; j<width; ++j )
            for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                data[iLocal+j*thisLocalHeight] =
                    ALocalBuffer[(thisColShift+iLocal*p)+j*ALDim];
    }

    // AllReduce over the grid
    mpi::AllReduce( sendBuffer, sendSize, mpi::SUM, g.VRComm() );

    // Unpack our received data
    T* thisLocalBuffer = this->LocalBuffer();
    const T* recvBuffer = &sendBuffer[VRRank*recvSize];
    const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
    {
        const T* recvBufferCol = &recvBuffer[j*localHeight];
        T* thisCol = &thisLocalBuffer[j*thisLDim];
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
            thisCol[iLocal] += alpha*recvBufferCol[iLocal];
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
