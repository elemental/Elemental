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
DistMatrix<T,MD,STAR,Int>::DistMatrix( const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ? 
    Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),0,
   0,0,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(0) ); }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ? 
    Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),width,
   g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(0) ); }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( bool constrainedColAlignment, Int colAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) : 0),0,
   0,0,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) :
    0),
   width,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
  Int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) :
    0),
   width,ldim,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, const T* buffer, Int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,colAlignment,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) :
    0),
   width,buffer,ldim,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, T* buffer, Int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,colAlignment,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) :
    0),
   width,buffer,ldim,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ); }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::DistMatrix");
#endif
    if( MD != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,MD,STAR,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [MD,* ] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::SetGrid( const elemental::Grid& g )
{
    this->Empty();
    this->grid_ = &g;
    this->colAlignment_ = 0;
    if( g.InGrid() && g.DiagPath()==g.DiagPath(0) )
    {
        inDiagonal_ = true;
        this->colShift_ = Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM());
    }
    else
        inDiagonal_ = false;
}

template<typename T,typename Int>
inline bool
DistMatrix<T,MD,STAR,Int>::InDiagonal() const
{ return inDiagonal_; }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MD,STAR,Int>::AlignWith( const DistMatrix<S,MD,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWith([MD,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->inDiagonal_ = A.InDiagonal();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        if( this->InDiagonal() )
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
DistMatrix<T,MD,STAR,Int>::AlignWith( const DistMatrix<S,STAR,MD,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWith([* ,MD])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->inDiagonal_ = A.InDiagonal();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        if( this->InDiagonal() )
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
DistMatrix<T,MD,STAR,Int>::AlignColsWith( const DistMatrix<S,MD,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MD,STAR,Int>::AlignColsWith( const DistMatrix<S,STAR,MD,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,MD,STAR,Int>::AlignedWithDiagonal
( const DistMatrix<S,MC,MR,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignedWithDiagonal([MC,MR])");
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        const Int ownerCol = (rowAlignment + offset) % c;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        const Int ownerCol = rowAlignment;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,MD,STAR,Int>::AlignedWithDiagonal
( const DistMatrix<S,MR,MC,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignedWithDiagonal([MR,MC])");
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerCol = colAlignment;
        const Int ownerRow = (rowAlignment + offset) % r;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const Int ownerCol = (colAlignment-offset) % c;
        const Int ownerRow = rowAlignment;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MD,STAR,Int>::AlignWithDiagonal
( const DistMatrix<S,MC,MR,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWithDiagonal([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        const Int ownerCol = (rowAlignment + offset) % c;
        this->colAlignment_ = ownerRow + r*ownerCol;
        if( g.InGrid() )
        {
            this->inDiagonal_ =
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->inDiagonal_ = false;
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        const Int ownerCol = rowAlignment;
        this->colAlignment_ = ownerRow + r*ownerCol;
        if( g.InGrid() )
        {
            this->inDiagonal_ =
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->inDiagonal_ = false;
    }
    if( this->InDiagonal() )
    {
        this->colShift_ =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
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
DistMatrix<T,MD,STAR,Int>::AlignWithDiagonal
( const DistMatrix<S,MR,MC,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWithDiagonal([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const Int ownerRow = rowAlignment;
        const Int ownerCol = (colAlignment + offset) % c;
        this->colAlignment_ = ownerRow + r*ownerCol;
        if( g.InGrid() )
        {
            this->inDiagonal_ =
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->inDiagonal_ = false;
    }
    else
    {
        const Int ownerRow = (rowAlignment-offset) % r;
        const Int ownerCol = colAlignment;
        this->colAlignment_ = ownerRow + r*ownerCol;
        if( g.InGrid() )
        {
            this->inDiagonal_ =
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->inDiagonal_ = false;
    }
    if( this->InDiagonal() )
    {
        this->colShift_ =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
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
DistMatrix<T,MD,STAR,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;
        
    const Int height      = this->Height();
    const Int width       = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int inDiagonal  = this->InDiagonal();
    const Int lcm         = g.LCM();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( g.InGrid() )
    {
        std::vector<T> sendBuf(height*width,0);
        if( inDiagonal )
        {
            const Int colShift = this->ColShift();
            const T* thisLocalBuffer = this->LockedLocalBuffer();
            const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                for( Int j=0; j<width; ++j )
                    sendBuf[(colShift+iLocal*lcm)+j*height] = 
                        thisLocalBuffer[iLocal+j*thisLDim];
        }

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
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::Align( Int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MD,STAR]::Align");
    this->AssertFreeColAlignment();
#endif
    this->AlignCols( colAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::AlignCols( Int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MD,STAR]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Size() )
        throw std::runtime_error("Invalid column alignment for [MD,STAR]");
#endif
    this->colAlignment_ = colAlignment;
    this->inDiagonal_ = ( g.DiagPath() == g.DiagPath(colAlignment) );
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        if( this->inDiagonal_ )
            this->colShift_ = Shift( g.DiagPathRank(), colAlignment, g.Size() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::View( DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->inDiagonal_ = A.InDiagonal();
    if( this->InDiagonal() )
    {
        this->colShift_ = A.ColShift();
        this->localMatrix_.View( A.LocalMatrix() );
    }
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::View
( Int height, Int width, Int colAlignment,
  T* buffer, Int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->inDiagonal_ = 
        grid.InGrid() && grid.DiagPath()==grid.DiagPath(colAlignment);
    if( this->inDiagonal_ )
    {
        this->colShift_ = 
            Shift(grid.DiagPathRank(),
                  grid.DiagPathRank(colAlignment),
                  grid.LCM());
        const Int localHeight = LocalLength(height,this->colShift_,grid.LCM());
        this->localMatrix_.View( localHeight, width, buffer, ldim );
    }
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::LockedView( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->inDiagonal_ = A.InDiagonal();
    if( this->InDiagonal() )
    {
        this->colShift_ = A.ColShift();
        this->localMatrix_.LockedView( A.LockedLocalMatrix() );
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::LockedView
( Int height, Int width, Int colAlignment,
  const T* buffer, Int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
#endif
    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->inDiagonal_ = 
        grid.InGrid() && grid.DiagPath()==grid.DiagPath(colAlignment);
    if( this->inDiagonal_ )
    {
        this->colShift_ = 
            Shift(grid.DiagPathRank(),
                  grid.DiagPathRank(colAlignment),
                  grid.LCM());
        const Int localHeight = LocalLength(height,this->colShift_,grid.LCM());
        this->localMatrix_.LockedView( localHeight, width, buffer, ldim );
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::View
( DistMatrix<T,MD,STAR,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    {
        const elemental::Grid& g = this->Grid();
        const Int r = g.Height();
        const Int c = g.Width();
        const Int lcm = g.LCM();
        const Int diagPathRank = g.DiagPathRank();
        const Int alignmentRank = A.ColAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int newAlignmentRow = (alignmentRow + i) % r;
        const Int newAlignmentCol = (alignmentCol + i) % c;
        const Int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        this->colAlignment_ = newAlignmentRank;
        this->inDiagonal_ = A.InDiagonal();

        if( this->inDiagonal_ )
        {
            this->colShift_ = 
                Shift( diagPathRank,
                       g.DiagPathRank( this->ColAlignment() ),
                       lcm );
            Int localHeightBefore = LocalLength( i, A.ColShift(), lcm);
            Int localHeight = LocalLength( height, this->ColShift(), lcm );

            this->localMatrix_.View
            ( A.LocalMatrix(), localHeightBefore, j, localHeight, width );
        }
    }
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::LockedView
( const DistMatrix<T,MD,STAR,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_  = width;
    {
        const elemental::Grid& g = this->Grid();
        const Int r = g.Height();
        const Int c = g.Width();
        const Int lcm = g.LCM();
        const Int diagPathRank = g.DiagPathRank();
        const Int alignmentRank = A.ColAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int newAlignmentRow = (alignmentRow + i) % r;
        const Int newAlignmentCol = (alignmentCol + i) % c;
        const Int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        this->colAlignment_ = newAlignmentRank;
        this->inDiagonal_ = A.InDiagonal();

        if( this->InDiagonal() )
        {
            this->colShift_ = 
                Shift( diagPathRank,
                       g.DiagPathRank( this->ColAlignment() ),
                       lcm );
            Int localHeightBefore = LocalLength( i, A.ColShift(), lcm);
            Int localHeight = LocalLength( height, this->ColShift(), lcm );
        
            this->localMatrix_.LockedView
            ( A.LockedLocalMatrix(), localHeightBefore, j, localHeight, width );
        }
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::View1x2
( DistMatrix<T,MD,STAR,Int>& AL, DistMatrix<T,MD,STAR,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View1x2");    
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->inDiagonal_ = AL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->colShift_ = AL.ColShift();
        this->localMatrix_.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    }
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::LockedView1x2
( const DistMatrix<T,MD,STAR,Int>& AL, const DistMatrix<T,MD,STAR,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView1x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->inDiagonal_ = AL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->colShift_ = AL.ColShift();
        this->localMatrix_.LockedView1x2
        ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::View2x1
( DistMatrix<T,MD,STAR,Int>& AT,
  DistMatrix<T,MD,STAR,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->inDiagonal_ = AT.InDiagonal();
    if( this->InDiagonal() )
    {
        this->colShift_ = AT.ColShift();
        this->localMatrix_.View2x1
        ( AT.LocalMatrix(),
          AB.LocalMatrix() );
    }
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::LockedView2x1
( const DistMatrix<T,MD,STAR,Int>& AT,
  const DistMatrix<T,MD,STAR,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->inDiagonal_ = AT.InDiagonal();
    if( this->InDiagonal() )
    {
        this->colShift_ = AT.ColShift();
        this->localMatrix_.LockedView2x1
        ( AT.LockedLocalMatrix(),
          AB.LockedLocalMatrix() );
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::View2x2
( DistMatrix<T,MD,STAR,Int>& ATL, DistMatrix<T,MD,STAR,Int>& ATR,
  DistMatrix<T,MD,STAR,Int>& ABL, DistMatrix<T,MD,STAR,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View2x2");
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
    this->inDiagonal_ = ATL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->colShift_ = ATL.ColShift();
        this->localMatrix_.View2x2
        ( ATL.LocalMatrix(), ATR.LocalMatrix(),
          ABL.LocalMatrix(), ABR.LocalMatrix() );
    }
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::LockedView2x2
( const DistMatrix<T,MD,STAR,Int>& ATL, const DistMatrix<T,MD,STAR,Int>& ATR,
  const DistMatrix<T,MD,STAR,Int>& ABL, const DistMatrix<T,MD,STAR,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView2x2");
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
    this->inDiagonal_ = ATL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->colShift_ = ATL.ColShift();
        this->localMatrix_.LockedView2x2
        ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
          ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::ResizeTo");
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
        ( LocalLength(height,this->ColShift(),lcm), width );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline T
DistMatrix<T,MD,STAR,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    Int ownerRank;
    const elemental::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->ColAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + i) % r;
        const Int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    T u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.LCM();
        u = this->GetLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elemental::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->ColAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + i) % r;
        const Int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.LCM();
        this->SetLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elemental::Grid& g = this->Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = this->ColAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + i) % r;
        const Int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.LCM();
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
DistMatrix<T,MD,STAR,Int>::MakeTrapezoidal
( Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int lcm = this->Grid().LCM();
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
                    Int numZeroRows = RawLocalLength( boundary, colShift, lcm );
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
                Int numNonzeroRows = RawLocalLength(firstZeroRow,colShift,lcm);
                if( numNonzeroRows < localHeight )
                {
                    T* thisCol = &thisLocalBuffer[numNonzeroRows+j*thisLDim];
                    std::memset
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
DistMatrix<T,MD,STAR,Int>::ScaleTrapezoid
( T alpha, Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::ScaleTrapezoid");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int lcm = this->Grid().LCM();
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
                Int numRows = RawLocalLength( boundary, colShift, lcm );
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
                Int firstRow = 
                    ( side==LEFT ? std::max(j-offset,0)
                                 : std::max(j+height-width-offset,0) );
                Int numZeroRows = RawLocalLength( firstRow, colShift, lcm );
                T* thisCol = &thisLocalBuffer[numZeroRows+j*thisLDim];
                for( Int iLocal=0; iLocal<(localHeight-numZeroRows); ++iLocal )
                    thisCol[iLocal] *= alpha;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const Int lcm = this->Grid().LCM();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int colShift = this->ColShift();

        this->localMatrix_.SetToZero();

        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const Int i = colShift + iLocal*lcm;
            if( i < width )
                thisLocalBuffer[iLocal+i*thisLDim] = 1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        for( Int j=0; j<width; ++j )
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                this->SetLocalEntry(iLocal,j,SampleUnitBall<T>());
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MD,* ] = [MC,MR] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MD,* ] = [MC,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MD,* ] = [* ,MR] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [MD,* ]");
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
            this->inDiagonal_ = A.InDiagonal();
            if( this->InDiagonal() )
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
        if( this->Grid().VCRank() == 0 )
            std::cerr << "Unaligned [MD,* ] <- [MD,* ]." << std::endl;
#endif
        throw std::logic_error
        ("Unaligned [MD,* ] = [MD,* ] not yet implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MD,* ] = [* ,MD] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A ); 
#endif
    throw std::logic_error("[MD,* ] = [MR,MC] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MD,* ] = [MR,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MD,* ] = [* ,MC] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MD,* ] = [VC,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MD,* ] = [* ,VC] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MD,* ] = [VR,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MD,* ] = [* ,VR] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MD,STAR,Int>&
DistMatrix<T,MD,STAR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    if( this->InDiagonal() )
    {
        const Int lcm = this->grid_->LCM();
        const Int colShift = this->ColShift();

        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();

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
                    ALocalBuffer[(colShift+iLocal*lcm)+j*ALDim];
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

} // namespace elemental
