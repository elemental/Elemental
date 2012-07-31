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
DistMatrix<T,STAR,MD,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),
   0,0,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(0) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),
   g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(0) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( bool constrainedRowAlignment, Int rowAlignment, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,constrainedRowAlignment,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   0,0,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
  Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   ldim,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   buffer,ldim,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   buffer,ldim,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  0,(A.Grid.InGrid() && A.Grid().DiagPath()==A.Grid().DiagPath(0) ?
     A.Grid().DiagPathRank() : 0),
  0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
    if( STAR != U || MD != V || 
        reinterpret_cast<const DistMatrix<T,STAR,MD,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MD] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::SetGrid( const elem::Grid& g )
{
    this->Empty();
    this->grid_ = &g;
    this->rowAlignment_ = 0;
    if( g.InGrid() && g.DiagPath()==g.DiagPath(0) )
    {
        inDiagonal_ = true;
        this->rowShift_ = Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM());
    }
    else
    {
        inDiagonal_ = false;
        this->rowShift_ = 0;
    }
}

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,MD,Int>::ColStride() const
{ return 1; }

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,MD,Int>::RowStride() const
{ return this->grid_->LCM(); }

template<typename T,typename Int>
inline bool
DistMatrix<T,STAR,MD,Int>::InDiagonal() const
{ return inDiagonal_; }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignWith( const DistMatrix<S,STAR,MD,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWith([* ,MD])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->Empty();
    this->rowAlignment_ = A.RowAlignment();
    this->inDiagonal_   = A.InDiagonal();
    this->constrainedRowAlignment_ = true;
    if( this->InDiagonal() )
        this->rowShift_ = A.RowShift();
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignWith( const DistMatrix<S,MD,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWith([MD,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->Empty();
    this->rowAlignment_ = A.ColAlignment();
    this->inDiagonal_   = A.InDiagonal();
    this->constrainedRowAlignment_ = true;
    if( this->InDiagonal() )
        this->rowShift_ = A.ColShift();
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignRowsWith( const DistMatrix<S,STAR,MD,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignRowsWith( const DistMatrix<S,MD,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,STAR,MD,Int>::AlignedWithDiagonal
( const DistMatrix<S,MC,MR,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignedWithDiagonal([MC,MR])");
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        const Int ownerCol = (rowAlignment + offset) % c;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        const Int ownerCol = rowAlignment;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,STAR,MD,Int>::AlignedWithDiagonal
( const DistMatrix<S,MR,MC,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignedWithDiagonal([MR,MC])");
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRow = rowAlignment;
        const Int ownerCol = (colAlignment + offset) % c;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const Int ownerRow = (rowAlignment-offset) % r;
        const Int ownerCol = colAlignment;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignWithDiagonal
( const DistMatrix<S,MC,MR,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiagonal([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();

    this->Empty();
    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        const Int ownerCol = (rowAlignment + offset) % c;
        this->rowAlignment_ = ownerRow + r*ownerCol;
        this->inDiagonal_ =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        const Int ownerCol = rowAlignment;
        this->rowAlignment_ = ownerRow + r*ownerCol;
        this->inDiagonal_ =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    this->constrainedRowAlignment_ = true;
    if( this->InDiagonal() )
        this->rowShift_ =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignWithDiagonal
( const DistMatrix<S,MR,MC,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiagonal([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();

    this->Empty();
    if( offset >= 0 )
    {
        const Int ownerRow = rowAlignment;
        const Int ownerCol = (colAlignment + offset) % c;
        this->rowAlignment_ = ownerRow + r*ownerCol;
        this->inDiagonal_ =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const Int ownerRow = (rowAlignment-offset) % r;
        const Int ownerCol = colAlignment;
        this->rowAlignment_ = ownerRow + r*ownerCol;
        this->inDiagonal_ =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    this->constrainedRowAlignment_ = true;
    if( this->InDiagonal() )
        this->rowShift_ =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::PrintBase");
#endif
    if( this->Grid().Rank() == 0 && msg != "" )
        os << msg << std::endl;
        
    const Int height     = this->Height();
    const Int width      = this->Width();
    const Int localWidth = this->LocalWidth();
    const Int lcm        = this->Grid().LCM();
    const Int inDiagonal = this->InDiagonal();

    if( height == 0 || width == 0 || !this->Grid().InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    std::vector<T> sendBuf(height*width,0);
    if( inDiagonal )
    {
        const Int colShift = this->ColShift();
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            T* destCol = &sendBuf[colShift+jLocal*lcm*height];
            const T* sourceCol = &thisLocalBuffer[jLocal*thisLDim];
            for( Int i=0; i<height; ++i )
                destCol[i] = sourceCol[i];
        }
    }

    // If we are the root, allocate a receive buffer
    std::vector<T> recvBuf;
    if( this->Grid().Rank() == 0 )
        recvBuf.resize( height*width );

    // Sum the contributions and send to the root
    mpi::Reduce
    ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, 
      this->Grid().Comm() );

    if( this->Grid().Rank() == 0 )
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::Align( Int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[STAR,MD]::Align");
    this->AssertFreeRowAlignment();
#endif
    this->AlignRows( rowAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::AlignRows( Int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[STAR,MD]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Size() )
        throw std::runtime_error("Invalid row alignment for [STAR,MD]");
#endif
    this->rowAlignment_ = rowAlignment;
    this->inDiagonal_ = ( g.DiagPath() == g.DiagPath(rowAlignment) );
    this->constrainedRowAlignment_ = true;
    if( this->inDiagonal_ )
        this->rowShift_ = Shift( g.DiagPathRank(), rowAlignment, g.Size() );
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::View( DistMatrix<T,STAR,MD,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View");
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_  = A.Width();
    this->rowAlignment_ = A.RowAlignment();
    this->inDiagonal_   = A.InDiagonal();
    this->viewing_ = true;
    if( this->InDiagonal() )
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
DistMatrix<T,STAR,MD,Int>::View
( Int height, Int width, Int rowAlignment,
  T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->inDiagonal_ = 
        grid.InGrid() && grid.DiagPath()==grid.DiagPath(rowAlignment);
    this->viewing_ = true;
    if( this->inDiagonal_ )
    {
        this->rowShift_ =
            Shift(grid.DiagPathRank(),
                  grid.DiagPathRank(rowAlignment),
                  grid.LCM());
        const Int localWidth = LocalLength(width,this->rowShift_,grid.LCM());
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
DistMatrix<T,STAR,MD,Int>::LockedView( const DistMatrix<T,STAR,MD,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView");
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_  = A.Width();
    this->rowAlignment_ = A.RowAlignment();
    this->inDiagonal_   = A.InDiagonal();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->InDiagonal() )
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
DistMatrix<T,STAR,MD,Int>::LockedView
( Int height, Int width, Int rowAlignment,
  const T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->rowAlignment_ = rowAlignment;
    this->inDiagonal_ = 
        grid.InGrid() && grid.DiagPath()==grid.DiagPath(rowAlignment);
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->inDiagonal_ )
    {
        this->rowShift_ =
            Shift(grid.DiagPathRank(),
                  grid.DiagPathRank(rowAlignment),
                  grid.LCM());
        const Int localWidth = LocalLength(width,this->rowShift_,grid.LCM());
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
DistMatrix<T,STAR,MD,Int>::View
( DistMatrix<T,STAR,MD,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View");
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_  = width;
    this->viewing_ = true;

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int diagPathRank = g.DiagPathRank(); 
    const Int alignmentRank = A.RowAlignment();
    const Int alignmentRow = alignmentRank % r;
    const Int alignmentCol = alignmentRank / r;
    const Int newAlignmentRow = (alignmentRow + i) % r;
    const Int newAlignmentCol = (alignmentCol + i) % c;
    const Int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

    this->rowAlignment_ = newAlignmentRank;
    this->inDiagonal_ = A.InDiagonal();

    if( this->InDiagonal() )
    {
        this->rowShift_ = 
            Shift( diagPathRank,
                   g.DiagPathRank(this->RowAlignment()),
                   lcm );
        Int localWidthBefore = LocalLength( j, A.RowShift(), lcm );
        Int localWidth = LocalLength( width, this->RowShift(), lcm );
    
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
DistMatrix<T,STAR,MD,Int>::LockedView
( const DistMatrix<T,STAR,MD,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView");
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_  = width;
    this->viewing_ = true;
    this->lockedView_ = true;

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int diagPathRank = g.DiagPathRank();
    const Int alignmentRank = A.RowAlignment();
    const Int alignmentRow = alignmentRank % r;
    const Int alignmentCol = alignmentRank / r;
    const Int newAlignmentRow = (alignmentRow + i) % r;
    const Int newAlignmentCol = (alignmentCol + i) % c;
    const Int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

    this->rowAlignment_ = newAlignmentRank;
    this->inDiagonal_ = A.InDiagonal();

    if( this->InDiagonal() )
    {
        this->rowShift_ = 
            Shift( diagPathRank,
                   g.DiagPathRank( this->RowAlignment() ),
                   lcm );
        Int localWidthBefore = LocalLength( j, A.RowShift(), lcm);
        Int localWidth = LocalLength( width, this->RowShift(), lcm );
    
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
DistMatrix<T,STAR,MD,Int>::View1x2
( DistMatrix<T,STAR,MD,Int>& AL, DistMatrix<T,STAR,MD,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View1x2");    
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->Empty();

    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_  = AL.Width() + AR.Width();
    this->rowAlignment_ = AL.RowAlignment();
    this->inDiagonal_ = AL.InDiagonal();
    this->viewing_ = true;
    if( this->InDiagonal() )
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
DistMatrix<T,STAR,MD,Int>::LockedView1x2
( const DistMatrix<T,STAR,MD,Int>& AL, const DistMatrix<T,STAR,MD,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView1x2");
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->Empty();

    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->rowAlignment_ = AL.RowAlignment();
    this->inDiagonal_ = AL.InDiagonal();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->InDiagonal() )
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
DistMatrix<T,STAR,MD,Int>::View2x1
( DistMatrix<T,STAR,MD,Int>& AT,
  DistMatrix<T,STAR,MD,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View2x1");
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->Empty();

    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->rowAlignment_ = AT.RowAlignment();
    this->inDiagonal_ = AT.InDiagonal();
    this->viewing_ = true;
    if( this->InDiagonal() )
    {
        this->rowShift_ = AT.RowShift();
        this->localMatrix_.View2x1
        ( AT.LocalMatrix(), AB.LocalMatrix() );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::LockedView2x1
( const DistMatrix<T,STAR,MD,Int>& AT,
  const DistMatrix<T,STAR,MD,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView2x1");
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->Empty();

    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->rowAlignment_ = AT.RowAlignment();
    this->inDiagonal_ = AT.InDiagonal();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->InDiagonal() )
    {
        this->rowShift_ = AT.RowShift();
        this->localMatrix_.LockedView2x1
        ( AT.LockedLocalMatrix(), AB.LockedLocalMatrix() );
    }
    else
        this->rowShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::View2x2
( DistMatrix<T,STAR,MD,Int>& ATL, DistMatrix<T,STAR,MD,Int>& ATR,
  DistMatrix<T,STAR,MD,Int>& ABL, DistMatrix<T,STAR,MD,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View2x2");
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
    this->inDiagonal_ = ATL.InDiagonal();
    this->viewing_ = true;
    if( this->InDiagonal() )
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
DistMatrix<T,STAR,MD,Int>::LockedView2x2
( const DistMatrix<T,STAR,MD,Int>& ATL, const DistMatrix<T,STAR,MD,Int>& ATR,
  const DistMatrix<T,STAR,MD,Int>& ABL, const DistMatrix<T,STAR,MD,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView2x2");
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
    this->inDiagonal_ = ATL.InDiagonal();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->InDiagonal() )
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
DistMatrix<T,STAR,MD,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->InDiagonal() )
    {
        const Int lcm = this->Grid().LCM();
        this->localMatrix_.ResizeTo
        ( height, LocalLength(width,this->RowShift(),lcm) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline T
DistMatrix<T,STAR,MD,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    Int ownerRank;
    const elem::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    T u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
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
DistMatrix<T,STAR,MD,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::Set");
    this->AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elem::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
        this->SetLocal(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::Update");
    this->AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elem::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-this->RowShift()) / g.LCM();
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
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MD] = [MC,MR] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MD] = [MC,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MD] = [* ,MR] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MD] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->rowAlignment_ = A.RowAlignment();
            this->inDiagonal_ = A.InDiagonal();
            if( this->InDiagonal() )
                this->rowShift_ = A.RowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        this->localMatrix_ = A.LockedLocalMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( this->Grid().Rank() == 0 )
            std::cerr << "Unaligned [* ,MD] <- [* ,MD]." << std::endl;
#endif
        throw std::logic_error
        ("Unaligned [* ,MD] = [* ,MD] not yet implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MD] = [MR,MC] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MD] = [MR,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MD] = [* ,MC] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MD] = [VC,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MD] = [* ,VC] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MD] = [VR,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[* ,MD] = [* ,VR] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,MD,Int>&
DistMatrix<T,STAR,MD,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    if( this->InDiagonal() )
    {
        const Int lcm = this->Grid().LCM();
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
            const T* ACol = &ALocalBuffer[(rowShift+jLocal*lcm)*ALDim];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( thisCol, ACol, height );
        }
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
DistMatrix<T,STAR,MD,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::GetRealPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner of entry (i,j) and broadcast from it
    Int ownerRank;
    const elem::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.LCM();
        u = this->GetLocalRealPart( i, jLocal );
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,MD,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::GetImagPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner of entry (i,j) and broadcast from it
    Int ownerRank;
    const elem::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.LCM();
        u = this->GetLocalImagPart( i, jLocal );
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::SetRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elem::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.LCM();
        this->SetLocalRealPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::SetImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    Int ownerRank;
    const elem::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.LCM();
        this->SetLocalImagPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::UpdateRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elem::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.LCM();
        this->UpdateLocalRealPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::UpdateImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    Int ownerRank;
    const elem::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int jLocal = (j-this->RowShift()) / g.LCM();
        this->UpdateLocalImagPart( i, jLocal, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
