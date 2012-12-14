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
DistMatrix<T,MD,STAR,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ? 
    Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),0,
   0,0,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(0) ); }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
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
( bool constrainedColAlignment, Int colAlignment, const elem::Grid& g )
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
  const elem::Grid& g )
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
  Int ldim, const elem::Grid& g )
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
  const elem::Grid& g )
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
  const elem::Grid& g )
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
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  (A.Grid().InGrid() && A.Grid().DiagPath()==A.Grid().DiagPath(0) ?
   A.Grid().DiagPathRank() : 0),0,
  0,0,A.Grid())
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
DistMatrix<T,MD,STAR,Int>::SetGrid( const elem::Grid& g )
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
    {
        inDiagonal_ = false;
        this->colShift_ = 0;
    }
}

template<typename T,typename Int>
inline Int
DistMatrix<T,MD,STAR,Int>::ColStride() const
{ return this->grid_->LCM(); }

template<typename T,typename Int>
inline Int
DistMatrix<T,MD,STAR,Int>::RowStride() const
{ return 1; }

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
    this->Empty();
    this->colAlignment_ = A.ColAlignment();
    this->inDiagonal_ = A.InDiagonal();
    this->constrainedColAlignment_ = true;
    if( this->inDiagonal_ )
        this->colShift_ = A.ColShift();
    else
        this->colShift_ = 0;
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
    this->Empty();
    this->colAlignment_ = A.RowAlignment();
    this->inDiagonal_ = A.InDiagonal();
    this->constrainedColAlignment_ = true;
    if( this->inDiagonal_ )
        this->colShift_ = A.RowShift();
    else
        this->colShift_ = 0;
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
    const elem::Grid& g = this->Grid();
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
    this->constrainedColAlignment_ = true;
    if( this->InDiagonal() )
        this->colShift_ =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    else
        this->colShift_ = 0;
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
    this->constrainedColAlignment_ = true;
    if( this->InDiagonal() )
        this->colShift_ =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    else
        this->colShift_ = 0;
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
    const elem::Grid& g = this->Grid();
    if( g.Rank() == 0 && msg != "" )
        os << msg << std::endl;
        
    const Int height      = this->Height();
    const Int width       = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int inDiagonal  = this->InDiagonal();
    const Int lcm         = g.LCM();

    if( height == 0 || width == 0 || !g.InGrid() )
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
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &sendBuf[colShift+j*height];
            const T* sourceCol = &thisLocalBuffer[j*thisLDim];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                destCol[iLocal*lcm] = sourceCol[iLocal];
        }
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
    const elem::Grid& g = this->Grid();
    this->Empty();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Size() )
        throw std::runtime_error("Invalid column alignment for [MD,STAR]");
#endif
    this->colAlignment_ = colAlignment;
    this->inDiagonal_ = ( g.DiagPath() == g.DiagPath(colAlignment) );
    this->constrainedColAlignment_ = true;
    if( this->inDiagonal_ )
        this->colShift_ = Shift( g.DiagPathRank(), colAlignment, g.Size() );
    else
        this->colShift_ = 0;
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
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->inDiagonal_ = A.InDiagonal();
    this->viewing_ = true;
    if( this->InDiagonal() )
    {
        this->colShift_ = A.ColShift();
        this->localMatrix_.View( A.LocalMatrix() );
    }
    else
        this->colShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::View
( Int height, Int width, Int colAlignment,
  T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->inDiagonal_ = 
        grid.InGrid() && grid.DiagPath()==grid.DiagPath(colAlignment);
    this->viewing_ = true;
    if( this->inDiagonal_ )
    {
        this->colShift_ = 
            Shift(grid.DiagPathRank(),
                  grid.DiagPathRank(colAlignment),
                  grid.LCM());
        const Int localHeight = LocalLength(height,this->colShift_,grid.LCM());
        this->localMatrix_.View( localHeight, width, buffer, ldim );
    }
    else
        this->colShift_ = 0;
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
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->colAlignment_ = A.ColAlignment();
    this->inDiagonal_ = A.InDiagonal();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->InDiagonal() )
    {
        this->colShift_ = A.ColShift();
        this->localMatrix_.LockedView( A.LockedLocalMatrix() );
    }
    else
        this->colShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::LockedView
( Int height, Int width, Int colAlignment,
  const T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->inDiagonal_ = 
        grid.InGrid() && grid.DiagPath()==grid.DiagPath(colAlignment);
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->inDiagonal_ )
    {
        this->colShift_ = 
            Shift(grid.DiagPathRank(),
                  grid.DiagPathRank(colAlignment),
                  grid.LCM());
        const Int localHeight = LocalLength(height,this->colShift_,grid.LCM());
        this->localMatrix_.LockedView( localHeight, width, buffer, ldim );
    }
    else
        this->colShift_ = 0;
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
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;

    const elem::Grid& g = this->Grid();
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
    else
        this->colShift_ = 0;
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
    else
        this->colShift_ = 0;
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
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->Empty();

    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->inDiagonal_ = AL.InDiagonal();
    this->viewing_ = true;
    if( this->InDiagonal() )
    {
        this->colShift_ = AL.ColShift();
        this->localMatrix_.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    }
    else
        this->colShift_ = 0;
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
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->Empty();

    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->colAlignment_ = AL.ColAlignment();
    this->inDiagonal_ = AL.InDiagonal();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->InDiagonal() )
    {
        this->colShift_ = AL.ColShift();
        this->localMatrix_.LockedView1x2
        ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    }
    else
        this->colShift_ = 0;
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
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->Empty();

    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->inDiagonal_ = AT.InDiagonal();
    this->viewing_ = true;
    if( this->InDiagonal() )
    {
        this->colShift_ = AT.ColShift();
        this->localMatrix_.View2x1
        ( AT.LocalMatrix(),
          AB.LocalMatrix() );
    }
    else
        this->colShift_ = 0;
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
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->Empty();

    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->colAlignment_ = AT.ColAlignment();
    this->inDiagonal_ = AT.InDiagonal();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->InDiagonal() )
    {
        this->colShift_ = AT.ColShift();
        this->localMatrix_.LockedView2x1
        ( AT.LockedLocalMatrix(),
          AB.LockedLocalMatrix() );
    }
    else
        this->colShift_ = 0;
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
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->Empty();

    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->colAlignment_ = ATL.ColAlignment();
    this->inDiagonal_ = ATL.InDiagonal();
    this->viewing_ = true;
    if( this->InDiagonal() )
    {
        this->colShift_ = ATL.ColShift();
        this->localMatrix_.View2x2
        ( ATL.LocalMatrix(), ATR.LocalMatrix(),
          ABL.LocalMatrix(), ABR.LocalMatrix() );
    }
    else
        this->colShift_ = 0;
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
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->Empty();

    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->colAlignment_ = ATL.ColAlignment();
    this->inDiagonal_ = ATL.InDiagonal();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->InDiagonal() )
    {
        this->colShift_ = ATL.ColShift();
        this->localMatrix_.LockedView2x2
        ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
          ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    }
    else
        this->colShift_ = 0;
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
        this->localMatrix_.ResizeTo
        ( LocalLength(height,this->ColShift(),this->Grid().LCM()), width );
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
    const elem::Grid& g = this->Grid();
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
        u = this->GetLocal(iLoc,j);
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
    const elem::Grid& g = this->Grid();
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
        this->SetLocal(iLoc,j,u);
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
    const elem::Grid& g = this->Grid();
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
        this->UpdateLocal(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., operator=
//

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
        if( this->Grid().Rank() == 0 )
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
#ifdef HAVE_OPENMP
        #pragma omp parallel for 
#endif
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisLocalBuffer[j*thisLDim];
            const T* sourceCol = &ALocalBuffer[colShift+j*ALDim];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                destCol[iLocal] = sourceCol[iLocal*lcm];
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
DistMatrix<T,MD,STAR,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetRealPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner of entry (i,j) and broadcast from it
    const elem::Grid& g = this->Grid();
    Int ownerRank;
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

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLocal = (i-this->ColShift()) / g.LCM();
        u = this->GetLocalRealPart( iLocal, j );
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,MD,STAR,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetImagPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner of entry (i,j) and broadcast from it
    const elem::Grid& g = this->Grid();
    Int ownerRank;
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

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLocal = (i-this->ColShift()) / g.LCM();
        u = this->GetLocalImagPart( iLocal, j );
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::SetRealPart( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    Int ownerRank;
    const elem::Grid& g = this->Grid();
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
        const Int iLocal = (i-this->ColShift()) / g.LCM();
        this->SetLocalRealPart( iLocal, j, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::SetImagPart( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    Int ownerRank;
    const elem::Grid& g = this->Grid();
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
        const Int iLocal = (i-this->ColShift()) / g.LCM();
        this->SetLocalImagPart( iLocal, j, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::UpdateRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    Int ownerRank;
    const elem::Grid& g = this->Grid();
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
        const Int iLocal = (i-this->ColShift()) / g.LCM();
        this->UpdateLocalRealPart( iLocal, j, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::UpdateImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    Int ownerRank;
    const elem::Grid& g = this->Grid();
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
        const Int iLocal = (i-this->ColShift()) / g.LCM();
        this->UpdateLocalImagPart( iLocal, j, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
