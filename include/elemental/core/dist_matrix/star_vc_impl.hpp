/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VC,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   0,(g.InGrid() ? g.VCRank() : 0 ),
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VC,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   0,(g.InGrid() ? g.VCRank() : 0),
   height,(g.InGrid() ? LocalLength(width,g.VCRank(),0,g.Size()) : 0),
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VC,Int>::DistMatrix
( bool constrainedRowAlignment, Int rowAlignment, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VCRank(),rowAlignment,g.Size()) : 0),
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VC,Int>::DistMatrix
( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VCRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.VCRank(),rowAlignment,g.Size()) : 0),
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VC,Int>::DistMatrix
( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
  Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VCRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.VCRank(),rowAlignment,g.Size()) : 0),
   ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VCRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.VCRank(),rowAlignment,g.Size()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VC,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.VCRank(),rowAlignment,g.Size()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.VCRank(),rowAlignment,g.Size()) : 0),
   buffer,ldim,g)
{ }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,STAR,VC,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  0,(A.Grid().InGrid() ? A.VCRank() : 0),
  0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::DistMatrix");
#endif
    if( STAR != U || VC != V || 
        reinterpret_cast<const DistMatrix<T,STAR,VC,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,VC] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,STAR,VC,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::SetGrid( const elem::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->rowAlignment_ = 0;
    if( grid.InGrid() )
        this->rowShift_ = grid.VCRank();
    else
        this->rowShift_ = 0;
}

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,VC,Int>::ColStride() const
{ return 1; }

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,VC,Int>::RowStride() const
{ return this->grid_->Size(); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->rowAlignment_ = A.ColAlignment();
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift( g.VCRank(), this->RowAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->rowAlignment_ = A.RowAlignment();
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift( g.VCRank(), this->RowAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignWith( const DistMatrix<S,MC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([MC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->rowAlignment_ = A.ColAlignment();
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift( g.VCRank(), this->RowAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignWith( const DistMatrix<S,STAR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([* ,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->rowAlignment_ = A.RowAlignment();
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift( g.VCRank(), this->RowAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignWith( const DistMatrix<S,STAR,VC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([* ,VC])");
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
DistMatrix<T,STAR,VC,Int>::AlignWith( const DistMatrix<S,VC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([VC,* ])");
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
DistMatrix<T,STAR,VC,Int>::AlignRowsWith( const DistMatrix<S,MC,MR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignRowsWith( const DistMatrix<S,MR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignRowsWith( const DistMatrix<S,MC,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignRowsWith( const DistMatrix<S,STAR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignRowsWith( const DistMatrix<S,STAR,VC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignRowsWith( const DistMatrix<S,VC,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,STAR,VC,Int>::AlignedWithDiagonal
( const DistMatrix<S,VC,STAR,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignedWithDiagonal([* ,MC])");
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int p = g.Size();
    const Int colAlignment = A.ColAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRank = colAlignment;
        aligned = ( this->RowAlignment() == ownerRank );
    }
    else
    {
        const Int ownerRank = (colAlignment-offset) % p;
        aligned = ( this->RowAlignment() == ownerRank );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,STAR,VC,Int>::AlignedWithDiagonal
( const DistMatrix<S,STAR,VC,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignedWithDiagonal([* ,VC])");
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int p = g.Size();
    const Int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRank = (rowAlignment + offset) % p;
        aligned = ( this->RowAlignment() == ownerRank );
    }
    else
    {
        const Int ownerRank = rowAlignment;
        aligned = ( this->RowAlignment() == ownerRank );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignWithDiagonal
( const DistMatrix<S,VC,STAR,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWithDiagonal([VC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int p = g.Size();
    const Int colAlignment = A.ColAlignment();

    this->Empty();
    if( offset >= 0 )
    {
        const Int ownerRank = colAlignment;
        this->rowAlignment_ = ownerRank;
    }
    else
    {
        const Int ownerRank = (colAlignment-offset) % p;
        this->rowAlignment_ = ownerRank;
    }
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift(g.VCRank(),this->rowAlignment_,p);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,VC,Int>::AlignWithDiagonal
( const DistMatrix<S,STAR,VC,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWithDiagonal([* ,VC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int p = g.Size();
    const Int rowAlignment = A.RowAlignment();

    this->Empty();
    if( offset >= 0 )
    {
        const Int ownerRank = (rowAlignment+offset) % p;
        this->rowAlignment_ = ownerRank;
    }
    else
    {
        const Int ownerRank = rowAlignment;
        this->rowAlignment_ = ownerRank;
    }
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift(g.VCRank(),this->rowAlignment_,p);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::PrintBase");
#endif
    const elem::Grid& g = this->Grid();
    if( g.Rank() == 0 && msg != "" )
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
    if( g.Rank() == 0 )
        recvBuf.resize( height*width );

    // Sum the contributions and send to the root
    mpi::Reduce
    ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, g.Comm() );

    if( g.Rank() == 0 )
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
    mpi::Barrier( g.Comm() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::Align( Int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Align");
    this->AssertFreeRowAlignment();
#endif
    this->AlignRows( rowAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::AlignRows( Int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elem::Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Size() )
        throw std::runtime_error("Invalid row alignment for [* ,VC]");
#endif
    this->Empty();
    this->rowAlignment_ = rowAlignment;
    this->constrainedRowAlignment_ = true;
    if( g.InGrid() )
        this->rowShift_ = Shift( g.VCRank(), rowAlignment, g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::View( DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View");
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
DistMatrix<T,STAR,VC,Int>::View
( Int height, Int width, Int rowAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->viewing_ = true;
    if( g.InGrid() )
    {
        this->rowShift_ = Shift(g.VCRank(),rowAlignment,g.Size());
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
DistMatrix<T,STAR,VC,Int>::LockedView( const DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView");
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
DistMatrix<T,STAR,VC,Int>::LockedView
( Int height, Int width, Int rowAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView");
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
        this->rowShift_ = Shift(g.VCRank(),rowAlignment,g.Size());
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
DistMatrix<T,STAR,VC,Int>::View
( DistMatrix<T,STAR,VC,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View");
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;

    const elem::Grid& g = this->Grid();
    const Int colMajorRank = g.VCRank();
    const Int size = g.Size();

    this->rowAlignment_ = (A.RowAlignment()+j) % size;
    this->viewing_ = true;

    if( g.InGrid() )
    {
        this->rowShift_ = Shift( colMajorRank, this->RowAlignment(), size );
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
DistMatrix<T,STAR,VC,Int>::LockedView
( const DistMatrix<T,STAR,VC,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView");
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;

    const elem::Grid& g = this->Grid();
    const Int colMajorRank = g.VCRank();
    const Int size = g.Size();

    this->rowAlignment_ = (A.RowAlignment()+j) % size;
    this->viewing_ = true;
    this->lockedView_ = true;

    if( g.InGrid() )
    {
        this->rowShift_ = Shift( colMajorRank, this->RowAlignment(), size );
        const Int localWidth = LocalLength( width, this->RowShift(), size );
        this->localMatrix_.LockedView
        ( A.LockedLocalMatrix(), i, A.RowShift(), height, localWidth );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::View1x2
( DistMatrix<T,STAR,VC,Int>& AL, DistMatrix<T,STAR,VC,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View1x2");
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
DistMatrix<T,STAR,VC,Int>::LockedView1x2
( const DistMatrix<T,STAR,VC,Int>& AL, const DistMatrix<T,STAR,VC,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView1x2");
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
DistMatrix<T,STAR,VC,Int>::View2x1
( DistMatrix<T,STAR,VC,Int>& AT,
  DistMatrix<T,STAR,VC,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View2x1");
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
DistMatrix<T,STAR,VC,Int>::LockedView2x1
( const DistMatrix<T,STAR,VC,Int>& AT,
  const DistMatrix<T,STAR,VC,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView2x1");
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
DistMatrix<T,STAR,VC,Int>::View2x2
( DistMatrix<T,STAR,VC,Int>& ATL, DistMatrix<T,STAR,VC,Int>& ATR,
  DistMatrix<T,STAR,VC,Int>& ABL, DistMatrix<T,STAR,VC,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View2x2");
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
DistMatrix<T,STAR,VC,Int>::LockedView2x2
( const DistMatrix<T,STAR,VC,Int>& ATL, const DistMatrix<T,STAR,VC,Int>& ATR,
  const DistMatrix<T,STAR,VC,Int>& ABL, const DistMatrix<T,STAR,VC,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView2x2");
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
DistMatrix<T,STAR,VC,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    const elem::Grid& g = this->Grid();
    this->height_ = height;
    this->width_  = width;
    if( g.InGrid() )
        this->localMatrix_.ResizeTo
        ( height, LocalLength( width, this->RowShift(), g.Size() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline T
DistMatrix<T,STAR,VC,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocal(i,jLoc);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
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
DistMatrix<T,STAR,VC,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        this->UpdateLocal(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., operator=
//

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR,Int> A_STAR_VR(g);

    A_STAR_VR = A;
    *this = A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MC,MR,Int> > A_MC_MR
    ( new DistMatrix<T,MC,MR,Int>(g) );
    *A_MC_MR = A;

    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR,Int>(g) );
    *A_STAR_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,STAR,VR,Int> A_STAR_VR(g);

    A_STAR_VR = A;
    *this = A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MD,* ]");
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
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MD]");
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
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MR,MC]");
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
                    Shift( g.VCRank(), this->RowAlignment(), g.Size() );
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

    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

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
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const Int thisRank = row+k*r;
            const Int thisRowShift = RawShift(thisRank,rowAlignment,p);
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // Communicate
        mpi::AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, g.RowComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const Int thisColShift = RawShift(k,colAlignmentOfA,c);
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,c);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                T* destCol = &thisLocalBuffer[thisColShift+jLocal*thisLDim];
                const T* sourceCol = &data[jLocal*thisLocalHeight];
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    destCol[iLocal*c] = sourceCol[iLocal];
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VC] <- [MR,MC]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlignment = this->RowAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int sendRow = (row+r+(rowAlignment%r)-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-(rowAlignment%r)) % r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int localHeightOfA = A.LocalHeight();

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
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const Int thisRank = sendRow+k*r;
            const Int thisRowShift = RawShift(thisRank,rowAlignment,p);
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*localHeightOfA];
                MemCopy( dataCol, ACol, localHeightOfA );
            }
        }

        // AllToAll to gather all of the unaligned [*,VC] data into firstBuffer
        mpi::AllToAll
        ( secondBuffer, portionSize,
          firstBuffer,  portionSize, g.RowComm() );

        // SendRecv: properly align the [*,VC] via a trade in the column
        mpi::SendRecv
        ( firstBuffer,  portionSize, sendRow, 0,
          secondBuffer, portionSize, recvRow, mpi::ANY_TAG, g.ColComm() );

        // Unpack
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const Int thisColShift = RawShift(k,colAlignmentOfA,c);
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,c);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                T* destCol = &thisLocalBuffer[thisColShift+jLocal*thisLDim];
                const T* sourceCol = &data[jLocal*thisLocalHeight];
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
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
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MR,* ]");
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
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MC]");
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
                    Shift( g.VCRank(), this->RowAlignment(), g.Size() );
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

    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int rowShift = this->RowShift();
        const Int rowShiftOfA = A.RowShift();
        const Int rowOffset = (rowShift-rowShiftOfA) / r;

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
            const T* ACol = &ALocalBuffer[(rowOffset+jLocal*c)*ALDim];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( thisCol, ACol, height );
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [* ,VC] <- [* ,MC]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int col = g.Col();
        const Int rowShiftOfA = A.RowShift();
        const Int rowAlignment = this->RowAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        // We will SendRecv A[*,VC] within our process column to fix alignments.
        const Int sendRow = (row+r+(rowAlignment%r)-rowAlignmentOfA) % r;
        const Int recvRow = (row+r+rowAlignmentOfA-(rowAlignment%r)) % r;
        const Int sendRank = sendRow + r*col;

        const Int sendRowShift = Shift( sendRank, rowAlignment, p );
        const Int sendRowOffset = (sendRowShift-rowShiftOfA) / r;

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
            const T* ACol = &ALocalBuffer[(sendRowOffset+jLocal*c)*ALDim];
            T* sendBufferCol = &sendBuffer[jLocal*height];
            MemCopy( sendBufferCol, ACol, height );
        }

        // Communicate
        mpi::SendRecv
        ( sendBuffer, sendSize, sendRow, 0,
          recvBuffer, recvSize, recvRow, mpi::ANY_TAG, g.ColComm() );

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
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MC,MR,Int> > A_MC_MR
    ( new DistMatrix<T,MC,MR,Int>(g) );
    *A_MC_MR = A;

    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > A_STAR_VR
    ( new DistMatrix<T,STAR,VR,Int>(g) );
    *A_STAR_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_STAR_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,VC]");
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
            std::cerr << "Unaligned [* ,VC] <- [* ,VC]." << std::endl;
#endif
        const Int rank = g.VCRank();
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
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.VCComm() );

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
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [VR,* ]");
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
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,VR]");
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

    // Compute which colmajor rank has the rowShift equal to our rowShiftOfA
    const Int sendRankCM = (rankCM+(p+rowShiftOfA-rowShift)) % p;

    // Compute which colmajor rank has the A rowShift that we need
    const Int recvRankRM = (rankRM+(p+rowShift-rowShiftOfA)) % p;
    const Int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

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
    ( sendBuffer, sendSize, sendRankCM, 0,
      recvBuffer, recvSize, recvRankCM, mpi::ANY_TAG, g.VCComm() );

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
inline const DistMatrix<T,STAR,VC,Int>&
DistMatrix<T,STAR,VC,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,* ]");
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
        const T* ACol = &ALocalBuffer[(rowShift+jLocal*p)*ALDim];
        T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
        MemCopy( thisCol, ACol, height );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::SumScatterFrom( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SumScatterFrom( [* ,MC] )");
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
                    Shift( g.VCRank(), this->RowAlignment(), g.Size() );
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

    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int rowAlignment = this->RowAlignment();
        const Int rowShiftOfA = A.RowShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalWidth = MaxLocalLength( width, p );

        const Int recvSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);
        const Int sendSize = c*recvSize;

        this->auxMemory_.Require( sendSize );
        T* buffer = this->auxMemory_.Buffer();

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[k*recvSize];

            const Int thisRank = myRow+k*r;
            const Int thisRowShift = RawShift( thisRank, rowAlignment, p );
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = RawLocalLength( width, thisRowShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*height];
                MemCopy( dataCol, ACol, height );
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.RowComm() );

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
        ("Unaligned [* ,VC]::ReduceScatterFrom( [* ,MC] ) not yet implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SumScatterUpdate( [* ,MC] )");
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

    if( this->RowAlignment() % g.Height() == A.RowAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int rowAlignment = this->RowAlignment();
        const Int rowShiftOfA = A.RowShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = this->LocalWidth();
        const Int maxLocalWidth = MaxLocalLength( width, p );

        const Int recvSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);
        const Int sendSize = c*recvSize;

        this->auxMemory_.Require( sendSize );
        T* buffer = this->auxMemory_.Buffer();

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            T* data = &buffer[k*recvSize];

            const Int thisRank = myRow+k*r;
            const Int thisRowShift = RawShift( thisRank, rowAlignment, p );
            const Int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const Int thisLocalWidth = RawLocalLength( width, thisRowShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* ACol = &ALocalBuffer[(thisRowOffset+jLocal*c)*ALDim];
                T* dataCol = &data[jLocal*height];
                MemCopy( dataCol, ACol, height );
            }
        }

        // Communicate
        mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.RowComm() );

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
        ("Unaligned [* ,VC]::ReduceScatterUpdate( [* ,MC] ) not implemented");
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
DistMatrix<T,STAR,VC,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::GetRealPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocalRealPart(i,jLoc);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,VC,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::GetImagPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocalImagPart(i,jLoc);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::SetRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.Size();
        this->SetLocalRealPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::SetImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.Size();
        this->SetLocalImagPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::UpdateRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.Size();
        this->UpdateLocalRealPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,VC,Int>::UpdateImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.Size();
        this->UpdateLocalImagPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
