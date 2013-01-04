/*
   Copyright (c) 2009-2013, Jack Poulson
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
   (g.InGrid() && g.DiagPath()==0 ? g.DiagPathRank() : 0),0,
   0,0,g), 
  diagPath_(0)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   (g.InGrid() && g.DiagPath()==0 ? g.DiagPathRank() : 0),0,
   (g.InGrid() && g.DiagPath()==0 ?
    LocalLength(height,g.DiagPathRank(),0,g.LCM()) : 0),width,g),
  diagPath_(0)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( bool constrainedColAlignment, Int colAlignmentVC, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,constrainedColAlignment,false,g.DiagPathRank(colAlignmentVC),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignmentVC) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignmentVC),g.LCM()) : 0),0,
   0,0,g),
  diagPath_(g.DiagPath(colAlignmentVC))
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignmentVC,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,g.DiagPathRank(colAlignmentVC),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignmentVC) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignmentVC),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignmentVC) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignmentVC),g.LCM())
    : 0),width,g),
  diagPath_(g.DiagPath(colAlignmentVC))
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignmentVC,
  Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,g.DiagPathRank(colAlignmentVC),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignmentVC) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignmentVC),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignmentVC) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignmentVC),g.LCM())
    : 0),width,ldim,g), 
  diagPath_(g.DiagPath(colAlignmentVC))
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignmentVC, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,g.DiagPathRank(colAlignmentVC),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignmentVC) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignmentVC),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignmentVC) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignmentVC),g.LCM())
    : 0),width,buffer,ldim,g),
  diagPath_(g.DiagPath(colAlignmentVC))
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignmentVC, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,g.DiagPathRank(colAlignmentVC),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignmentVC) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignmentVC),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignmentVC) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignmentVC),g.LCM())
    : 0),width,buffer,ldim,g),
  diagPath_(g.DiagPath(colAlignmentVC))
{ }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,MD,STAR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  (A.Participating() ? A.ColRank() : 0),0,
  0,0,A.Grid()),
  diagPath_(A.diagPath_)
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
    this->diagPath_ = 0;
    this->colAlignment_ = 0;
    if( g.InGrid() && g.DiagPath()==0 )
        this->colShift_ = g.DiagPathRank();
    else
        this->colShift_ = 0;
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
inline Int
DistMatrix<T,MD,STAR,Int>::ColRank() const
{ return this->grid_->DiagPathRank(); }

template<typename T,typename Int>
inline Int
DistMatrix<T,MD,STAR,Int>::RowRank() const
{ return 0; }

template<typename T,typename Int>
inline bool
DistMatrix<T,MD,STAR,Int>::Participating() const
{
    const Grid& g = this->Grid();
    return ( g.InGrid() && g.DiagPath()==this->diagPath_ );
}

template<typename T,typename Int>
inline Int
DistMatrix<T,MD,STAR,Int>::DiagPath() const
{ return this->diagPath_; }

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
    this->diagPath_ = A.diagPath_;
    this->colAlignment_ = A.colAlignment_;
    this->constrainedColAlignment_ = true;
    this->colShift_ = A.ColShift();
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
    this->diagPath_ = A.diagPath_;
    this->colAlignment_ = A.rowAlignment_;
    this->constrainedColAlignment_ = true;
    this->colShift_ = A.RowShift();
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

    const Int firstDiagRow = 0;
    const Int firstDiagCol = this->diagPath_;
    const Int diagRow = (firstDiagRow+this->ColAlignment()) % r;
    const Int diagCol = (firstDiagCol+this->ColAlignment()) % c;

    bool aligned;
    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        const Int ownerCol = (rowAlignment + offset) % c;
        aligned = ( ownerRow==diagRow && ownerCol==diagCol );
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        const Int ownerCol = rowAlignment;
        aligned = ( ownerRow==diagRow && ownerCol==diagCol );
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

    const Int firstDiagRow = 0;
    const Int firstDiagCol = this->diagPath_;
    const Int diagRow = (firstDiagRow+this->ColAlignment()) % r;
    const Int diagCol = (firstDiagCol+this->ColAlignment()) % c;

    bool aligned;
    if( offset >= 0 )
    {
        const Int ownerCol = colAlignment;
        const Int ownerRow = (rowAlignment + offset) % r;
        aligned = ( ownerRow==diagRow && ownerCol==diagCol );
    }
    else
    {
        const Int ownerCol = (colAlignment-offset) % c;
        const Int ownerRow = rowAlignment;
        aligned = ( ownerRow==diagRow && ownerCol==diagCol );
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
    Int owner;
    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        const Int ownerCol = (rowAlignment + offset) % c;
        owner = ownerRow + r*ownerCol;
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        const Int ownerCol = rowAlignment;
        owner = ownerRow + r*ownerCol;
    }
    this->diagPath_ = g.DiagPath(owner);
    this->colAlignment_ = g.DiagPathRank(owner);
    this->constrainedColAlignment_ = true;
    if( this->Participating() )
        this->colShift_ = (g.DiagPathRank()+lcm-this->colAlignment_) % lcm;
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
    Int owner;
    if( offset >= 0 )
    {
        const Int ownerRow = rowAlignment;
        const Int ownerCol = (colAlignment + offset) % c;
        owner = ownerRow + r*ownerCol;
    }
    else
    {
        const Int ownerRow = (rowAlignment-offset) % r;
        const Int ownerCol = colAlignment;
        owner = ownerRow + r*ownerCol;
    }
    this->diagPath_ = g.DiagPath(owner);
    this->colAlignment_ = g.DiagPathRank(owner);
    this->constrainedColAlignment_ = true;
    if( this->Participating() )
        this->colShift_ = (g.DiagPathRank()+lcm-this->colAlignment_) % lcm;
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
    const Int lcm         = g.LCM();

    if( height == 0 || width == 0 || !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    std::vector<T> sendBuf(height*width,0);
    if( this->Participating() )
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
DistMatrix<T,MD,STAR,Int>::Align( Int colAlignmentVC )
{
#ifndef RELEASE
    PushCallStack("[MD,STAR]::Align");
    this->AssertFreeColAlignment();
#endif
    this->AlignCols( colAlignmentVC );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::AlignCols( Int colAlignmentVC )
{
#ifndef RELEASE
    PushCallStack("[MD,STAR]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
#ifndef RELEASE
    if( colAlignmentVC < 0 || colAlignmentVC >= g.Size() )
        throw std::runtime_error("Invalid column alignment for [MD,STAR]");
#endif
    this->diagPath_ = g.DiagPath(colAlignmentVC);
    this->colAlignment_ = g.DiagPathRank(colAlignmentVC);
    this->constrainedColAlignment_ = true;
    if( this->Participating() )
        this->colShift_ = 
            Shift( g.DiagPathRank(), this->colAlignment_, g.LCM() );
    else
        this->colShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::Attach
( Int height, Int width, Int colAlignmentVC,
  T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Attach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->diagPath_ = grid.DiagPath(colAlignmentVC);
    this->colAlignment_ = grid.DiagPathRank(colAlignmentVC);
    this->viewing_ = true;
    if( this->Participating() )
    {
        this->colShift_ = 
            Shift(grid.DiagPathRank(),this->colAlignment_,grid.LCM());
        const Int localHeight = LocalLength(height,this->colShift_,grid.LCM());
        this->localMatrix_.Attach( localHeight, width, buffer, ldim );
    }
    else
        this->colShift_ = 0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::LockedAttach
( Int height, Int width, Int colAlignmentVC,
  const T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->diagPath_ = grid.DiagPath(colAlignmentVC);
    this->colAlignment_ = grid.DiagPathRank(colAlignmentVC);
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Participating() )
    {
        this->colShift_ = 
            Shift( grid.DiagPathRank(), this->colAlignment_, grid.LCM() );
        const Int localHeight = LocalLength(height,this->colShift_,grid.LCM());
        this->localMatrix_.LockedAttach( localHeight, width, buffer, ldim );
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
    if( this->Participating() )
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
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

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
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;
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
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;
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
            this->diagPath_ = A.diagPath_;
            this->colAlignment_ = A.colAlignment_;
            if( this->Participating() )
                this->colShift_ = A.ColShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->diagPath_ == A.diagPath_ && 
        this->colAlignment_ == A.colAlignment_ )
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

    if( this->Participating() )
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
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

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
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

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

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

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

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

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

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

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

    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int ownerRow = (i + this->colAlignment_) % r;
    const Int ownerCol = (i + this->colAlignment_ + this->diagPath_) % c;
    const Int ownerRank = ownerRow + r*ownerCol;

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
