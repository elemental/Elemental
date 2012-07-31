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

namespace elem {

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VR,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   0,(g.InGrid() ? g.VRRank() : 0 ),
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VR,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   0,(g.InGrid() ? g.VRRank() : 0),
   height,(g.InGrid() ? LocalLength(width,g.VRRank(),0,g.Size()) : 0),
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VR,Int>::DistMatrix
( bool constrainedRowAlignment, Int rowAlignment, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VRRank(),rowAlignment,g.Size()) : 0),
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VR,Int>::DistMatrix
( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VRRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.VRRank(),rowAlignment,g.Size()) : 0),
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VR,Int>::DistMatrix
( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
  Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VRRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.VRRank(),rowAlignment,g.Size()) : 0),
   ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VR,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VRRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.VRRank(),rowAlignment,g.Size()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VR,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VRRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.VRRank(),rowAlignment,g.Size()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,STAR,VR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  0,(A.Grid().InGrid() ? A.VRRank() : 0),
  0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::DistMatrix");
#endif
    if( STAR != U || VR != V || 
        reinterpret_cast<const DistMatrix<T,STAR,VR,Int>*>(&A) != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,VR] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::SetGrid( const elem::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->rowAlignment_ = 0;
    if( grid.InGrid() )
        this->rowShift_ = grid.VRRank();
    else
        this->rowShift_ = 0;
}

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,VR,Int>::ColStride() const
{ return 1; }

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,VR,Int>::RowStride() const
{ return this->grid_->Size(); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->rowAlignment_ = A.RowAlignment();
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift( g.VRRank(), this->RowAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->rowAlignment_ = A.ColAlignment();
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift( g.VRRank(), this->RowAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignWith( const DistMatrix<S,MR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([MR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->rowAlignment_ = A.ColAlignment();
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift( g.VRRank(), this->RowAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignWith( const DistMatrix<S,STAR,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([* ,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->rowAlignment_ = A.RowAlignment();
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift( g.VRRank(), this->RowAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignWith( const DistMatrix<S,STAR,VR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([* ,VR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->Empty();
    this->rowAlignment_ = A.RowAlignment();
    this->constrainedRowAlignment_ = true;
    if( this->Grid().InGrid() )
        this->rowShift_ = A.RowShift();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignWith( const DistMatrix<S,VR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([VR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->Empty();
    this->rowAlignment_ = A.ColAlignment();
    this->constrainedRowAlignment_ = true;
    if( this->Grid().InGrid() )
        this->rowShift_ = A.ColShift();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignRowsWith( const DistMatrix<S,MC,MR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignRowsWith( const DistMatrix<S,MR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignRowsWith( const DistMatrix<S,MR,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignRowsWith( const DistMatrix<S,STAR,MR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignRowsWith( const DistMatrix<S,STAR,VR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VR,Int>::AlignRowsWith( const DistMatrix<S,VR,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::PrintBase");
#endif
    const elem::Grid& g = this->Grid();
    if( g.VRRank() == 0 && msg != "" )
        os << msg << std::endl;

    const Int height     = this->Height();
    const Int width      = this->Width();
    const Int localWidth = this->LocalWidth();
    const Int p          = g.Size();
    const Int rowShift   = this->RowShift();

    if( height == 0 || width == 0 || !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    std::vector<T> sendBuf(height*width,0);
    const T* thisLocalBuffer = this->LockedLocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        T* destCol = &sendBuf[(rowShift+jLocal*p)*height];
        const T* sourceCol = &thisLocalBuffer[jLocal*thisLDim];
        MemCopy( destCol, sourceCol, height );
    }

    // If we are the root, allocate a receive buffer
    std::vector<T> recvBuf;
    if( g.VRRank() == 0 )
        recvBuf.resize(height*width);

    // Sum the contributions and send to the root
    mpi::Reduce
    ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, g.VRComm() );

    if( g.VRRank() == 0 )
    {
        // Print the data
        for( Int i=0; i<height; ++i )
        {
            for( Int j=0; j<width; ++j )
                os << recvBuf[i+j*height] << " ";
            os << "\n";
        }
        os << std::endl;
    }
    mpi::Barrier( g.VRComm() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::Align( Int rowAlignment )
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

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::AlignRows( Int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elem::Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Size() )
        throw std::runtime_error("Invalid row alignment for [* ,VR]");
#endif
    this->Empty();
    this->rowAlignment_ = rowAlignment;
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift( g.VRRank(), rowAlignment, g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::View( DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View");
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->rowAlignment_ = A.RowAlignment();
    this->viewing_ = true;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.RowShift();
        this->localMatrix_.View( A.LocalMatrix() );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::View
( Int height, Int width, Int rowAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->viewing_ = true;
    if( g.InGrid() )
    {
        this->rowShift_ = Shift(g.VRRank(),rowAlignment,g.Size());
        const Int localWidth = LocalLength(width,this->rowShift_,g.Size());
        this->localMatrix_.View( height, localWidth, buffer, ldim );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::LockedView( const DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView(A)");
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->rowAlignment_ = A.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.RowShift();
        this->localMatrix_.LockedView( A.LockedLocalMatrix() );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::LockedView
( Int height, Int width, Int rowAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->viewing_ = true;
    this->lockedView_ = true;
    if( g.InGrid() )
    {
        this->rowShift_ = Shift(g.VRRank(),rowAlignment,g.Size());
        const Int localWidth = LocalLength(width,this->rowShift_,g.Size());
        this->localMatrix_.LockedView( height, localWidth, buffer, ldim );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::View
( DistMatrix<T,STAR,VR,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View");
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;

    const elem::Grid& g = this->Grid();
    const Int rowMajorRank = g.VRRank();
    const Int size = g.Size();

    this->rowAlignment_ = (A.RowAlignment()+j) % size;
    this->viewing_ = true;
    
    if( g.InGrid() )
    {
        this->rowShift_ = Shift( rowMajorRank, this->RowAlignment(), size );
        const Int localWidthBefore = LocalLength( j, A.RowShift(), size );
        const Int localWidth = LocalLength( width, this->RowShift(), size );
        this->localMatrix_.View
        ( A.LocalMatrix(), i, localWidthBefore, height, localWidth );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::LockedView
( const DistMatrix<T,STAR,VR,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView");
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = true;

    const elem::Grid& g = this->Grid();
    const Int rowMajorRank = g.VRRank();
    const Int size = g.Size();

    this->rowAlignment_ = (A.RowAlignment()+j) % size;

    if( g.InGrid() )
    {
        this->rowShift_ = Shift( rowMajorRank, this->RowAlignment(), size );
        const Int localWidthBefore = LocalLength( j, A.RowShift(), size );
        const Int localWidth = LocalLength( width, this->RowShift(), size );
        this->localMatrix_.LockedView
        ( A.LockedLocalMatrix(), i, localWidthBefore, height, localWidth );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::View1x2
( DistMatrix<T,STAR,VR,Int>& AL, DistMatrix<T,STAR,VR,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View1x2");
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->Empty();

    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->rowAlignment_ = AL.RowAlignment();
    this->viewing_ = true;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = AL.RowShift();
        this->localMatrix_.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::LockedView1x2
( const DistMatrix<T,STAR,VR,Int>& AL, const DistMatrix<T,STAR,VR,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView1x2");
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->Empty();

    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->rowAlignment_ = AL.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = AL.RowShift();
        this->localMatrix_.LockedView1x2
        ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::View2x1
( DistMatrix<T,STAR,VR,Int>& AT,
  DistMatrix<T,STAR,VR,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View2x1");
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->Empty();

    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->rowAlignment_ = AT.RowAlignment();
    this->viewing_ = true;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = AT.RowShift();
        this->localMatrix_.View2x1( AT.LocalMatrix(), AB.LocalMatrix() );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::LockedView2x1
( const DistMatrix<T,STAR,VR,Int>& AT,
  const DistMatrix<T,STAR,VR,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView2x1");
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->Empty();

    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->rowAlignment_ = AT.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = AT.RowShift();
        this->localMatrix_.LockedView2x1
        ( AT.LockedLocalMatrix(), 
          AB.LockedLocalMatrix() );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::View2x2
( DistMatrix<T,STAR,VR,Int>& ATL, DistMatrix<T,STAR,VR,Int>& ATR,
  DistMatrix<T,STAR,VR,Int>& ABL, DistMatrix<T,STAR,VR,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View2x2");
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->Empty();

    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->rowAlignment_ = ATL.RowAlignment();
    this->viewing_ = true;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = ATL.RowShift();
        this->localMatrix_.View2x2
        ( ATL.LocalMatrix(), ATR.LocalMatrix(),
          ABL.LocalMatrix(), ABR.LocalMatrix() );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::LockedView2x2
( const DistMatrix<T,STAR,VR,Int>& ATL, const DistMatrix<T,STAR,VR,Int>& ATR,
  const DistMatrix<T,STAR,VR,Int>& ABL, const DistMatrix<T,STAR,VR,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView2x2");
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->Empty();

    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->rowAlignment_ = ATL.RowAlignment();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = ATL.RowShift();
        this->localMatrix_.LockedView2x2
        ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
          ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    const elem::Grid& g = this->Grid();
    this->height_ = height;
    this->width_ = width;
    if( g.InGrid() )
        this->localMatrix_.ResizeTo
        ( height, LocalLength(width,this->RowShift(),g.Size()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline T
DistMatrix<T,STAR,VR,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::Get");
    this->AssertValidEntry( i, j );
    // TODO: Generalize this function to always work...
    if( !this->Grid().InGrid() )
        throw std::logic_error("Should only call with processes in grid");
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    T u;
    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocal(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VRComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->SetLocal(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->UpdateLocal(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., AdjointFrom
//

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::AdjointFrom( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[*, VR]::AdjointFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.ColAlignment();
            if( g.InGrid() )
                this->rowShift_ = 
                    Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( this->RowAlignment() % g.Width() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int rowShift = this->RowShift();
        const Int colShiftOfA = A.ColShift();
        const Int rowOffset = (rowShift-colShiftOfA) / c;

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();

        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            T* destCol = &thisLocalBuffer[jLocal*thisLDim];
            const T* sourceCol = &ALocalBuffer[rowOffset+jLocal*r];
            for( Int i=0; i<height; ++i )
                destCol[i] = Conj( sourceCol[i*ALDim] );
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VR]::AdjointFrom" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int col = g.Col();
        const Int colShiftOfA = A.ColShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const Int sendCol = (col+c+(rowAlignment%c)-colAlignmentOfA) % c;
        const Int recvCol = (col+c+colAlignmentOfA-(rowAlignment%c)) % c;
        const Int sendRank = sendCol + c*row;

        const Int sendRowShift = Shift( sendRank, rowAlignment, p );
        const Int sendRowOffset = (sendRowShift-colShiftOfA) / c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfSend = LocalLength(width,sendRowShift,p);

        const Int sendSize = height * localWidthOfSend;
        const Int recvSize = height * localWidth;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidthOfSend; ++jLocal )
        {
            T* destCol = &sendBuffer[jLocal*height];
            const T* sourceCol = &ALocalBuffer[sendRowOffset+jLocal*r];
            for( Int i=0; i<height; ++i )
                destCol[i] = Conj( sourceCol[i*ALDim] );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, mpi::ANY_TAG, g.RowComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufferCol = &recvBuffer[jLocal*height];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( thisCol, recvBufferCol, height );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::TransposeFrom( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR]::TransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.ColAlignment();
            if( g.InGrid() )
                this->rowShift_ = 
                    Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( this->RowAlignment() % g.Width() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int rowShift = this->RowShift();
        const Int colShiftOfA = A.ColShift();
        const Int rowOffset = (rowShift-colShiftOfA) / c;

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();

        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for 
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            T* destCol = &thisLocalBuffer[jLocal*thisLDim];
            const T* sourceCol = &ALocalBuffer[rowOffset+jLocal*r];
            for( Int i=0; i<height; ++i )
                destCol[i] = sourceCol[i*ALDim];
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VR]::TransposeFrom" << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int col = g.Col();
        const Int colShiftOfA = A.ColShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const Int sendCol = (col+c+(rowAlignment%c)-colAlignmentOfA) % c;
        const Int recvCol = (col+c+colAlignmentOfA-(rowAlignment%c)) % c;
        const Int sendRank = sendCol + c*row;

        const Int sendRowShift = Shift( sendRank, rowAlignment, p );
        const Int sendRowOffset = (sendRowShift-colShiftOfA) / c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfSend = LocalLength(width,sendRowShift,p);

        const Int sendSize = height * localWidthOfSend;
        const Int recvSize = height * localWidth;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for 
#endif
        for( Int jLocal=0; jLocal<localWidthOfSend; ++jLocal )
        {
            T* destCol = &sendBuffer[jLocal*height];
            const T* sourceCol = &ALocalBuffer[sendRowOffset+jLocal*r];
            for( Int i=0; i<height; ++i )
                destCol[i] = sourceCol[i*ALDim];
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, mpi::ANY_TAG, g.RowComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufferCol = &recvBuffer[jLocal*height];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( thisCol, recvBufferCol, height );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            if( g.InGrid() )
                this->rowShift_ = 
                    Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return *this;
    }

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int maxHeight = MaxLocalLength(height,r);
        const Int maxWidth = MaxLocalLength(width,p);
        const Int portionSize = std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( 2*r*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[r*portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const Int thisRank = col+k*c;
            const Int thisRowShift = RawShift(thisRank,rowAlignment,p);
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*r)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, g.ColComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const Int thisColShift = RawShift(k,colAlignmentOfA,r);
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,r);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                T* destCol = &thisLocalBuffer[thisColShift+jLocal*thisLDim];
                const T* sourceCol = &data[jLocal*thisLocalHeight];
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    destCol[iLocal*r] = sourceCol[iLocal];
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VR] <- [MC,MR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int col = g.Col();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendCol = (col+c+(rowAlignment%c)-rowAlignmentOfA) % c;
        const Int recvCol = (col+c+rowAlignmentOfA-(rowAlignment%c)) % c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

        const Int maxHeight = MaxLocalLength(height,r);
        const Int maxWidth = MaxLocalLength(width,p);
        const Int portionSize = std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( 2*r*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[r*portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const Int thisRank = sendCol+k*c;
            const Int thisRowShift = RawShift(thisRank,rowAlignment,p);
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*r)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // AllToAll to gather all of the unaligned [*,VR] data into firstBuffer
        mpi::AllToAll
        ( secondBuffer, portionSize,
          firstBuffer,  portionSize, g.ColComm() );

        // SendRecv: properly align the [*,VR] via a trade in the column
        mpi::SendRecv
        ( firstBuffer,  portionSize, sendCol, 0,
          secondBuffer, portionSize, recvCol, mpi::ANY_TAG, g.RowComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const Int thisColShift = RawShift(k,colAlignmentOfA,r);
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,r);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                T* destCol = &thisLocalBuffer[thisColShift+jLocal*thisLDim];
                const T* sourceCol = &data[jLocal*thisLocalHeight];
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    destCol[iLocal*r] = sourceCol[iLocal];
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
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR,Int> A_MC_MR(g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            if( g.InGrid() )
                this->rowShift_ = 
                    Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return *this;
    }

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int rowShift = this->RowShift();
        const Int rowShiftOfA = A.RowShift();
        const Int rowOffset = (rowShift-rowShiftOfA) / c;

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();

        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[(rowOffset+jLocal*r)*ALDim];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( thisCol, ACol, height );
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VR] <- [* ,MR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int col = g.Col();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const Int sendCol = (col+c+(rowAlignment%c)-rowAlignmentOfA) % c;
        const Int recvCol = (col+c+rowAlignmentOfA-(rowAlignment%c)) % c;
        const Int sendRank = sendCol + c*row;

        const Int sendRowShift = Shift( sendRank, rowAlignment, p );
        const Int sendRowOffset = (sendRowShift-rowShiftOfA) / c;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfSend = LocalLength(width,sendRowShift,p);

        const Int sendSize = height * localWidthOfSend;
        const Int recvSize = height * localWidth;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidthOfSend; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[(sendRowOffset+jLocal*r)*ALDim];
            T* sendBufferCol = &sendBuffer[jLocal*height];
            MemCopy( sendBufferCol, ACol, height );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, mpi::ANY_TAG, g.RowComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufferCol = &recvBuffer[jLocal*height];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( thisCol, recvBufferCol, height );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    // TODO: Optimize this later if important
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VC,Int> A_STAR_VC(g);

    A_STAR_VC = A;
    *this = A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MR,MC,Int> > A_MR_MC
    ( new DistMatrix<T,MR,MC,Int>(g) );
    *A_MR_MC = A;

    std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC,Int>(g) );
    *A_STAR_VC = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VC,Int> A_STAR_VC(g);

    A_STAR_VC = A;
    *this = A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR,Int> A_MC_MR(g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return *this;
    }
    
    const Int height = this->Height();
    const Int localWidth = this->LocalWidth();
    const Int localWidthOfA = A.LocalWidth();

    const Int sendSize = height * localWidthOfA;
    const Int recvSize = height * localWidth;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int rankCM = g.VCRank();
    const Int rankRM = g.VRRank();

    const Int rowShift = this->RowShift();
    const Int rowShiftOfA = A.RowShift();

    // Compute which rowmajor rank has the rowShift equal to our rowShiftOfA
    const Int sendRankRM = (rankRM+(p+rowShiftOfA-rowShift)) % p;

    // Compute which rowmajor rank has the A rowShift that we need
    const Int recvRankCM = (rankCM+(p+rowShift-rowShiftOfA)) % p;
    const Int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

    this->auxMemory_.Require( sendSize + recvSize );

    T* buffer = this->auxMemory_.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
    {
        const T* ACol = &ALocalBuffer[jLocal*ALDim];
        T* sendBufferCol = &sendBuffer[jLocal*height];
        MemCopy( sendBufferCol, ACol, height );
    }

    // Communicate
    mpi::SendRecv
    ( sendBuffer, sendSize, sendRankRM, 0,
      recvBuffer, recvSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

    // Unpack
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* recvBufferCol = &recvBuffer[jLocal*height];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        MemCopy( thisCol, recvBufferCol, height );
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MR,MC,Int> > A_MR_MC
    ( new DistMatrix<T,MR,MC,Int>(g) );
    *A_MR_MC = A;

    std::auto_ptr<DistMatrix<T,STAR,VC,Int> > A_STAR_VC
    ( new DistMatrix<T,STAR,VC,Int>(g) );
    *A_STAR_VC = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_STAR_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            if( g.InGrid() )
                this->rowShift_ = A.RowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return *this;
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        this->localMatrix_ = A.LockedLocalMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VR] <- [* ,VR]." << std::endl;
#endif
        const Int rank = g.VRRank();
        const Int p = g.Size();

        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const Int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();
        const Int localWidthOfA = A.LocalWidth();

        const Int sendSize = height * localWidthOfA;
        const Int recvSize = height * localWidth;

        this->auxMemory_.Require( sendSize + recvSize );

        T* buffer = this->auxMemory_.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[jLocal*ALDim];
            T* sendBufferCol = &sendBuffer[jLocal*height];
            MemCopy( sendBufferCol, ACol, height );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.VRComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufferCol = &recvBuffer[jLocal*height];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( thisCol, recvBufferCol, height );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VR,Int>&
DistMatrix<T,STAR,VR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return *this;
    }

    const Int p = g.Size();
    const Int rowShift = this->RowShift();

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();

    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* ACol = &ALocalBuffer[(rowShift+jLocal*p)*ALDim];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        MemCopy( thisCol, ACol, localHeight );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::SumScatterFrom( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SumScatterFrom( [* ,MR] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            if( g.InGrid() )
                this->rowShift_ = 
                    Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myCol = g.Col();
        const Int rowAlignment = this->RowAlignment();
        const Int rowShiftOfA = A.RowShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalWidth = MaxLocalLength( width, p );

        const Int recvSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);
        const Int sendSize = r*recvSize;

        this->auxMemory_.Require( sendSize );
        T* buffer = this->auxMemory_.Buffer();

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[k*recvSize];

            const Int thisRank = myCol+k*c;
            const Int thisRowShift = RawShift( thisRank, rowAlignment, p );
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const Int thisLocalWidth = RawLocalLength( width, thisRowShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*r)*ALDim];
                T* dataCol = &data[jLocal*height];
                MemCopy( dataCol, ACol, height );
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.ColComm() );

        // Unpack our received data
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* bufferCol = &buffer[jLocal*height];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( thisCol, bufferCol, height );
        }
        this->auxMemory_.Release();
    }
    else
    {
        throw std::logic_error
        ("Unaligned [* ,VR]::ReduceScatterFrom( [* ,MR] ) not implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SumScatterUpdate( [* ,MR] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myCol = g.Col();
        const Int rowAlignment = this->RowAlignment();
        const Int rowShiftOfA = A.RowShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalWidth = MaxLocalLength( width, p );

        const Int recvSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);
        const Int sendSize = r*recvSize;

        this->auxMemory_.Require( sendSize );
        T* buffer = this->auxMemory_.Buffer();

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            T* data = &buffer[k*recvSize];

            const Int thisRank = myCol+k*c;
            const Int thisRowShift = RawShift( thisRank, rowAlignment, p );
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const Int thisLocalWidth = RawLocalLength( width, thisRowShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*r)*ALDim];
                T* dataCol = &data[jLocal*height];
                MemCopy( dataCol, ACol, height );
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.ColComm() );

        // Unpack our received data
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* bufferCol = &buffer[jLocal*height];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            for( Int i=0; i<height; ++i )
                thisCol[i] += alpha*bufferCol[i];
        }
        this->auxMemory_.Release();
    }
    else
    {
        throw std::logic_error
        ("Unaligned [* ,VR]::ReduceScatterUpdate( [* ,MR] ) not implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Routines which explicitly work in the complex plane
//

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,VR,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::GetRealPart");
    this->AssertValidEntry( i, j );
    // TODO: Generalize this function to always work...
    if( !this->Grid().InGrid() )
        throw std::logic_error("Should only call with processes in grid");
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    R u;
    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocalRealPart(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VRComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,VR,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::GetImagPart");
    this->AssertValidEntry( i, j );
    // TODO: Generalize this function to always work...
    if( !this->Grid().InGrid() )
        throw std::logic_error("Should only call with processes in grid");
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    R u;
    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocalImagPart(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VRComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::SetRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->SetLocalRealPart(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::SetImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->SetLocalImagPart(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::UpdateRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->UpdateLocalRealPart(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VR,Int>::UpdateImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->UpdateLocalImagPart(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
