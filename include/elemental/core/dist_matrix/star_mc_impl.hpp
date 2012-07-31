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
DistMatrix<T,STAR,MC,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   0,(g.InGrid() ? g.Row() : 0 ),
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MC,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   0,(g.InGrid() ? g.Row() : 0),
   height,(g.InGrid() ? LocalLength(width,g.Row(),0,g.Height()) : 0),
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MC,Int>::DistMatrix
( bool constrainedRowAlignment, Int rowAlignment, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.Row(),rowAlignment,g.Height()) : 0),
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MC,Int>::DistMatrix
( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.Row(),rowAlignment,g.Height()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.Row(),rowAlignment,g.Height()) : 0),
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MC,Int>::DistMatrix
( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
  Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.Row(),rowAlignment,g.Height()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.Row(),rowAlignment,g.Height()) : 0),
   ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, const T* buffer, Int ldim, 
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.Row(),rowAlignment,g.Height()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.Row(),rowAlignment,g.Height()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, T* buffer, Int ldim, 
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.Row(),rowAlignment,g.Height()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.Row(),rowAlignment,g.Height()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,STAR,MC,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  0,(A.Grid().InGrid() ? A.Grid().Row() : 0),
  0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::DistMatrix");
#endif
    if( STAR != U || MC != V || 
        reinterpret_cast<const DistMatrix<T,STAR,MC,Int>*>(&A) != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MC] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MC,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::SetGrid( const elem::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->rowAlignment_ = 0;
    if( grid.InGrid() )
        this->rowShift_ = grid.Row();
    else
        this->rowShift_ = 0;
}

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,MC,Int>::ColStride() const
{ return 1; }

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,MC,Int>::RowStride() const
{ return this->grid_->Height(); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MC,Int>::AlignWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWith([MR,MC])");
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
DistMatrix<T,STAR,MC,Int>::AlignWith( const DistMatrix<S,STAR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWith([* ,MC])");
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
DistMatrix<T,STAR,MC,Int>::AlignWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWith([MC,MR])");
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
DistMatrix<T,STAR,MC,Int>::AlignWith( const DistMatrix<S,MC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWith([MC,* ])");
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
DistMatrix<T,STAR,MC,Int>::AlignRowsWith( const DistMatrix<S,MC,MR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MC,Int>::AlignRowsWith( const DistMatrix<S,MC,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MC,Int>::AlignRowsWith( const DistMatrix<S,STAR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MC,Int>::AlignRowsWith( const DistMatrix<S,MR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,STAR,MC,Int>::AlignedWithDiagonal
( const DistMatrix<S,MC,STAR,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignedWithDiagonal([* ,MC])");
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int colAlignment = A.ColAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        aligned = ( this->RowAlignment() == ownerRow );
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        aligned = ( this->RowAlignment() == ownerRow );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,STAR,MC,Int>::AlignedWithDiagonal
( const DistMatrix<S,STAR,MC,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignedWithDiagonal([* ,MC])");
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRow = (rowAlignment + offset) % r;
        aligned = ( this->RowAlignment() == ownerRow );
    }
    else
    {
        const Int ownerRow = rowAlignment;
        aligned = ( this->RowAlignment() == ownerRow );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MC,Int>::AlignWithDiagonal
( const DistMatrix<S,MC,STAR,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWithDiagonal([MC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int colAlignment = A.ColAlignment();

    this->Empty();
    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        this->rowAlignment_ = ownerRow;
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        this->rowAlignment_ = ownerRow;
    }
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift(g.Row(),this->rowAlignment_,r);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MC,Int>::AlignWithDiagonal
( const DistMatrix<S,STAR,MC,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWithDiagonal([* ,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int rowAlignment = A.RowAlignment();

    this->Empty();
    if( offset >= 0 )
    {
        const Int ownerRow = (rowAlignment+offset) % r;
        this->rowAlignment_ = ownerRow;
    }
    else
    {
        const Int ownerRow = rowAlignment;
        this->rowAlignment_ = ownerRow;
    }
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift(g.Row(),this->rowAlignment_,r);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::PrintBase");
#endif
    const elem::Grid& g = this->Grid();
    if( g.Rank() == 0 && msg != "" )
        os << msg << std::endl;

    const Int height     = this->Height();
    const Int width      = this->Width();
    const Int localWidth = this->LocalWidth();
    const Int r          = g.Height();
    const Int rowShift   = this->RowShift();

    if( height == 0 || width == 0 || !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Only one process col needs to participate
    if( g.Col() == 0 )
    {
        std::vector<T> sendBuf(height*width,0);
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            T* destCol = &sendBuf[(rowShift+jLocal*r)*height];
            const T* sourceCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( destCol, sourceCol, height );
        }

        // If we are the root, allocate a receive buffer
        std::vector<T> recvBuf;
        if( g.Row() == 0 )
            recvBuf.resize( height*width );

        // Sum the contributions and send to the root
        mpi::Reduce
        ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, g.ColComm() );

        if( g.Row() == 0 )
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
    }
    mpi::Barrier( g.VCComm() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::Align( Int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::Align");
    this->AssertFreeRowAlignment();
#endif
    this->AlignRows( rowAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::AlignRows( Int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elem::Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Height() )
        throw std::runtime_error("Invalid row alignment for [* ,MC]");
#endif
    this->Empty();
    this->rowAlignment_ = rowAlignment;
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift( g.Row(), rowAlignment, g.Height() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::View( DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View");
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
DistMatrix<T,STAR,MC,Int>::View
( Int height, Int width, Int rowAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->viewing_ = true;
    if( g.InGrid() )
    {
        this->rowShift_ = Shift(g.Row(),rowAlignment,g.Height());
        const Int localWidth = LocalLength(width,this->rowShift_,g.Height());
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
DistMatrix<T,STAR,MC,Int>::LockedView( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE 
    PushCallStack("[* ,MC]::LockedView");
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
DistMatrix<T,STAR,MC,Int>::LockedView
( Int height, Int width, Int rowAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::LockedView");
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
        this->rowShift_ = Shift(g.Row(),rowAlignment,g.Height());
        const Int localWidth = LocalLength(width,this->rowShift_,g.Height());
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
DistMatrix<T,STAR,MC,Int>::View
( DistMatrix<T,STAR,MC,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View");
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;

    const elem::Grid& g = this->Grid();
    const Int r   = g.Height();
    const Int row = g.Row();

    this->rowAlignment_ = (A.RowAlignment()+j) % r;
    this->viewing_ = true;

    if( g.InGrid() )
    {
        this->rowShift_ = Shift( row, this->RowAlignment(), r ); 

        const Int localWidthBefore = LocalLength( j, A.RowShift(), r );
        const Int localWidth = LocalLength( width, this->RowShift(), r );

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
DistMatrix<T,STAR,MC,Int>::LockedView
( const DistMatrix<T,STAR,MC,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::LockedView");
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;

    const elem::Grid& g = this->Grid();
    const Int r   = g.Height();
    const Int row = g.Row();

    this->rowAlignment_ = (A.RowAlignment()+j) % r;
    this->viewing_ = true;
    this->lockedView_ = true;

    if( g.InGrid() )
    {
        this->rowShift_ = Shift( row, this->RowAlignment(), r );

        const Int localWidthBefore = LocalLength( j, A.RowShift(), r );
        const Int localWidth = LocalLength( width, this->RowShift(), r );

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
DistMatrix<T,STAR,MC,Int>::View1x2
( DistMatrix<T,STAR,MC,Int>& AL, DistMatrix<T,STAR,MC,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View1x2");
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
DistMatrix<T,STAR,MC,Int>::LockedView1x2
( const DistMatrix<T,STAR,MC,Int>& AL, const DistMatrix<T,STAR,MC,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::LockedView1x2");
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
DistMatrix<T,STAR,MC,Int>::View2x1
( DistMatrix<T,STAR,MC,Int>& AT,
  DistMatrix<T,STAR,MC,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View2x1");
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
        this->localMatrix_.View2x1
        ( AT.LocalMatrix(),
          AB.LocalMatrix() );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::LockedView2x1
( const DistMatrix<T,STAR,MC,Int>& AT,
  const DistMatrix<T,STAR,MC,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::LockedView2x1");
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
DistMatrix<T,STAR,MC,Int>::View2x2
( DistMatrix<T,STAR,MC,Int>& ATL, DistMatrix<T,STAR,MC,Int>& ATR,
  DistMatrix<T,STAR,MC,Int>& ABL, DistMatrix<T,STAR,MC,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::View2x2");
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
DistMatrix<T,STAR,MC,Int>::LockedView2x2
( const DistMatrix<T,STAR,MC,Int>& ATL, const DistMatrix<T,STAR,MC,Int>& ATR,
  const DistMatrix<T,STAR,MC,Int>& ABL, const DistMatrix<T,STAR,MC,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::LockedView2x2");
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
DistMatrix<T,STAR,MC,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Grid().InGrid() )
        this->localMatrix_.ResizeTo
        ( height, LocalLength(width,this->RowShift(),this->Grid().Height()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline T
DistMatrix<T,STAR,MC,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Get");
    this->AssertValidEntry( i, j );
    if( !this->Grid().InGrid() )
        throw std::logic_error("Should only be called by processes in grid");
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that
    // row within each process column
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    T u;
    if( g.Row() == ownerRow )
    {
        const Int jLoc = (j-this->RowShift()) / g.Height();
        u = this->GetLocal(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRow, g.ColComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->SetLocal(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateLocal(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., SumOverRow
//

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::SumOverRow()
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::SumOverRow");
    this->AssertNotLockedView();
#endif
    if( !this->Grid().InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int localSize = std::max( localHeight*localWidth, mpi::MIN_COLL_MSG );

    this->auxMemory_.Require( 2*localSize );
    T* buffer = this->auxMemory_.Buffer();
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[localSize];

    // Pack
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        T* sendBufCol = &sendBuf[jLocal*localHeight];
        MemCopy( sendBufCol, thisCol, localHeight );
    }

    // AllReduce sum
    mpi::AllReduce
    ( sendBuf, recvBuf, localSize, mpi::SUM, this->Grid().RowComm() );

    // Unpack
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const T* recvBufCol = &recvBuf[jLocal*localHeight];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        MemCopy( thisCol, recvBufCol, localHeight );
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::AdjointFrom( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC]::AdjointFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.ColAlignment() % g.Height();
            if( g.InGrid() )
                this->rowShift_ = 
                    Shift( g.Row(), this->RowAlignment(), g.Height() );
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

    if( this->RowAlignment() == A.ColAlignment() % g.Height() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLocalLength(width,p);

        const Int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for 
#endif
        for( Int jLocal=0; jLocal<localHeightOfA; ++jLocal )
        {
            T* destCol = &originalData[jLocal*height];
            const T* sourceCol = &ALocalBuffer[jLocal];
            for( Int i=0; i<height; ++i )
                destCol[i] = Conj( sourceCol[i*ALDim] );
        }

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.RowComm() );

        // Unpack
        const Int rowShift = this->RowShift();
        const Int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int colShiftOfA = RawShift( row+k*r, colAlignmentOfA, p );
            const Int rowOffset = (colShiftOfA-rowShift) / r;
            const Int localWidth = RawLocalLength( width, colShiftOfA, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MC]::AdjointFrom." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int rank = g.VCRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowShift = this->RowShift();

        const Int sendRank = (rank+p+rowAlignment-colAlignmentOfA) % p;
        const Int recvRank = (rank+p+colAlignmentOfA-rowAlignment) % p;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLocalLength(width,p);

        const Int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for 
#endif
        for( Int jLocal=0; jLocal<localHeightOfA; ++jLocal )
        {
            T* destCol = &secondBuffer[jLocal*height];
            const T* sourceCol = &ALocalBuffer[jLocal];
            for( Int i=0; i<height; ++i )
                destCol[i] = Conj( sourceCol[i*ALDim] );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, 0, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.RowComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const Int colShiftOfA = RawShift(row+r*k,rowAlignment,p);
            const Int rowOffset = (colShiftOfA-rowShift) / r;
            const Int localWidth = RawLocalLength( width, colShiftOfA, p );
            
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::TransposeFrom
( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC]::TransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.ColAlignment() % g.Height();
            if( g.InGrid() )
                this->rowShift_ = 
                    Shift( g.Row(), this->RowAlignment(), g.Height() );
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

    if( this->RowAlignment() == A.ColAlignment() % g.Height() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLocalLength(width,p);

        const Int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localHeightOfA; ++jLocal )
        {
            T* destCol = &originalData[jLocal*height];
            const T* sourceCol = &ALocalBuffer[jLocal];
            for( Int i=0; i<height; ++i )
                destCol[i] = sourceCol[i*ALDim];
        }

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.RowComm() );

        // Unpack
        const Int rowShift = this->RowShift();
        const Int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int colShiftOfA = RawShift( row+k*r, colAlignmentOfA, p );
            const Int rowOffset = (colShiftOfA-rowShift) / r;
            const Int localWidth = RawLocalLength( width, colShiftOfA, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MC]::TransposeFrom." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int rank = g.VCRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowShift = this->RowShift();

        const Int sendRank = (rank+p+rowAlignment-colAlignmentOfA) % p;
        const Int recvRank = (rank+p+colAlignmentOfA-rowAlignment) % p;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLocalLength(width,p);

        const Int portionSize = 
            std::max(height*maxLocalHeightOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localHeightOfA; ++jLocal )
        {
            T* destCol = &secondBuffer[jLocal*height];
            const T* sourceCol = &ALocalBuffer[jLocal];
            for( Int i=0; i<height; ++i )
                destCol[i] = sourceCol[i*ALDim];
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, 0, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.RowComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const Int colShiftOfA = RawShift(row+r*k,rowAlignment,p);
            const Int rowOffset = (colShiftOfA-rowShift) / r;
            const Int localWidth = RawLocalLength( width, colShiftOfA, p );
            
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
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
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MC,MR,Int> > A_MC_MR
    ( new DistMatrix<T,MC,MR,Int>(g) );
    *A_MC_MR   = A;

    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR,Int>(g) );
    *A_STAR_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

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
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    if( A.Height() == 1 )
    {
        if( !this->Viewing() )
            this->ResizeTo( 1, A.Width() );
        if( !g.InGrid() )
        {
#ifndef RELEASE
            PopCallStack();
#endif
            return *this;
        }

        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int myRow = g.Row();
        const Int rankCM = g.VCRank();
        const Int rankRM = g.VRRank();
        const Int rowAlignment = this->RowAlignment();
        const Int rowShift = this->RowShift();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int rowShiftOfA = A.RowShift();

        const Int width = this->Width();
        const Int maxLocalVectorWidth = MaxLocalLength(width,p);
        const Int portionSize = std::max(maxLocalVectorWidth,mpi::MIN_COLL_MSG);

        const Int rowShiftVC = Shift(rankCM,rowAlignment,p);
        const Int rowShiftVROfA = Shift(rankRM,rowAlignmentOfA,p);
        const Int sendRankCM = (rankCM+(p+rowShiftVROfA-rowShiftVC)) % p;
        const Int recvRankRM = (rankRM+(p+rowShiftVC-rowShiftVROfA)) % p;
        const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        this->auxMemory_.Require( (c+1)*portionSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        // A[* ,VR] <- A[* ,MR]
        {
            const Int shift = Shift(rankRM,rowAlignmentOfA,p);
            const Int offset = (shift-rowShiftOfA) / c;
            const Int thisLocalWidth = LocalLength(width,shift,p);

            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                sendBuf[jLocal] = ALocalBuffer[(offset+jLocal*r)*ALDim];
        }

        // A[* ,VC] <- A[* ,VR]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankCM, 0,
          recvBuf, portionSize, recvRankCM, mpi::ANY_TAG, g.VCComm() );

        // A[* ,MC] <- A[* ,VC]
        mpi::AllGather
        ( recvBuf, portionSize,
          sendBuf, portionSize, g.RowComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &sendBuf[k*portionSize];

            const Int shift = RawShift(myRow+r*k,rowAlignment,p);
            const Int offset = (shift-rowShift) / r;
            const Int thisLocalWidth = RawLocalLength(width,shift,p);

            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
                thisLocalBuffer[(offset+jLocal*c)*thisLDim] = data[jLocal];
        }
        this->auxMemory_.Release();
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
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MC] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MC] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Height() == 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "The vector version of [* ,MC] <- [MR,MC] is not yet written, but"
          " it would only require a modification of the vector version of "
          "[* ,MR] <- [MC,MR]." << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "The redistribution [* ,MC] <- [MR,MC] potentially causes a large"
          " amount of cache-thrashing. If possible, avoid it. "
          "Unfortunately, the following routines are not yet implemented: \n"
          << "  [MC,* ].(Conjugate)TransposeFrom([MR,MC])" << std::endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            if( g.InGrid() )
                this->rowShift_ = 
                    Shift( g.Row(), this->RowAlignment(), g.Height() );
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
        const Int c = g.Width();
        const Int height = this->Height();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLocalLength(height,c);

        const Int portionSize = 
            std::max(maxLocalHeightOfA*localWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[jLocal*ALDim];
            T* originalDataCol = &originalData[jLocal*localHeightOfA];
            MemCopy( originalDataCol, ACol, localHeightOfA );
        }

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.RowComm() );

        // Unpack
        const Int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int colShift = RawShift( k, colAlignmentOfA, c );
            const Int localHeight = RawLocalLength( height, colShift, c );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                T* destCol = &thisLocalBuffer[colShift+jLocal*thisLDim];
                const T* sourceCol = &data[jLocal*localHeight];
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    destCol[iLocal*c] = sourceCol[iLocal];
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MC] <- [MR,MC]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.Row();

        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalHeightOfA = MaxLocalLength(height,c);
        const Int maxLocalWidth = MaxLocalLength(width,r);

        const Int portionSize = 
            std::max(maxLocalHeightOfA*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[jLocal*ALDim];
            T* secondBufferCol = &secondBuffer[jLocal*localHeightOfA];
            MemCopy( secondBufferCol, ACol, localHeightOfA );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuffer, portionSize, sendRow, 0,
          firstBuffer,  portionSize, recvRow, mpi::ANY_TAG, g.ColComm() );

        // Use the output of the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.RowComm() );

        // Unpack the contents of each member of the process row
        const Int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const Int colShift = RawShift( k, colAlignmentOfA, c );
            const Int localHeight = RawLocalLength( height, colShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                T* destCol = &thisLocalBuffer[colShift+jLocal*thisLDim];
                const T* sourceCol = &data[jLocal*localHeight];
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    destCol[iLocal*c] = sourceCol[iLocal];
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
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MR,MC,Int> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,MC]");
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
            std::cerr << "Unaligned [* ,MC] <- [* ,MC]." << std::endl;
#endif
        const Int rank = g.Row();
        const Int r = g.Height();

        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRank = (rank+r+rowAlignment-rowAlignmentOfA) % r;
        const Int recvRank = (rank+r+rowAlignmentOfA-rowAlignment) % r;

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
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.ColComm() );

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
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = A;

    std::auto_ptr<DistMatrix<T,MR,MC,Int> > A_MR_MC
    ( new DistMatrix<T,MR,MC,Int>(false,true,0,this->RowAlignment(),g) );
    *A_MR_MC = *A_VR_STAR;
    delete A_VR_STAR.release();

    *this = *A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,VC]");
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
            this->rowAlignment_ = A.RowAlignment() % g.Height();
            if( g.InGrid() )
                this->rowShift_ = 
                    Shift( g.Row(), this->RowAlignment(), g.Height() );
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

    if( this->RowAlignment() == A.RowAlignment() % g.Height() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidthOfA = MaxLocalLength(width,p);

        const Int portionSize = 
            std::max(height*maxLocalWidthOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[jLocal*ALDim];
            T* originalDataCol = &originalData[jLocal*height];
            MemCopy( originalDataCol, ACol, height );
        }

        // Communicate
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.RowComm() );

        // Unpack
        const Int rowShift = this->RowShift();
        const Int rowAlignmentOfA = A.RowAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int rowShiftOfA = RawShift( row+k*r, rowAlignmentOfA, p );
            const Int rowOffset = (rowShiftOfA-rowShift) / r;
            const Int localWidth = RawLocalLength( width, rowShiftOfA, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,MC] <- [* ,VC]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int rank = g.VCRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        const Int rowShift = this->RowShift();

        const Int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const Int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidthOfA = MaxLocalLength(width,p);

        const Int portionSize = 
            std::max(height*maxLocalWidthOfA,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[jLocal*ALDim];
            T* secondBufferCol = &secondBuffer[jLocal*height];
            MemCopy( secondBufferCol, ACol, height );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, 0, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        mpi::AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.RowComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const Int rowShiftOfA = RawShift(row+r*k,rowAlignment,p);
            const Int rowOffset = (rowShiftOfA-rowShift) / r;
            const Int localWidth = RawLocalLength( width, rowShiftOfA, p );
            
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = &thisLocalBuffer[(rowOffset+jLocal*c)*thisLDim];
                MemCopy( thisCol, dataCol, height );
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
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MR,MC,Int> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VC,Int> A_STAR_VC(true,this->RowAlignment(),g);
    *this = A_STAR_VC = A;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MC,Int>&
DistMatrix<T,STAR,MC,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MC] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( !this->Grid().InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return *this;
    }

    const Int r = this->Grid().Height();
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
        const T* ACol = &ALocalBuffer[(rowShift+jLocal*r)*ALDim];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        MemCopy( thisCol, ACol, localHeight );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

//
// Routines which explicitly work in the complex plane
//

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,MC,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetRealPart");
    this->AssertValidEntry( i, j );
    if( !this->Grid().InGrid() )
        throw std::logic_error("Should only be called by processes in grid");
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner row of entry (i,j) and broadcast from that
    // row within each process column
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    R u;
    if( g.Row() == ownerRow )
    {
        const Int jLocal = (j-this->RowShift()) / g.Height();
        u = this->GetLocalRealPart( i, jLocal );
    }
    mpi::Broadcast( &u, 1, ownerRow, g.ColComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,MC,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetImagPart");
    this->AssertValidEntry( i, j );
    if( !this->Grid().InGrid() )
        throw std::logic_error("Should only be called by processes in grid");
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner row of entry (i,j) and broadcast from that
    // row within each process column
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    R u;
    if( g.Row() == ownerRow )
    {
        const Int jLocal = (j-this->RowShift()) / g.Height();
        u = this->GetLocalImagPart( i, jLocal );
    }
    mpi::Broadcast( &u, 1, ownerRow, g.ColComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::SetRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int jLocal = (j-this->RowShift()) / g.Height();
        this->SetLocalRealPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::SetImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int jLocal = (j-this->RowShift()) / g.Height();
        this->SetLocalImagPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::UpdateRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int jLocal = (j-this->RowShift()) / g.Height();
        this->UpdateLocalRealPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MC,Int>::UpdateImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRow = (j + this->RowAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int jLocal = (j-this->RowShift()) / g.Height();
        this->UpdateLocalImagPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
