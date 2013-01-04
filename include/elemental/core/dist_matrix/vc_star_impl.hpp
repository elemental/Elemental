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
DistMatrix<T,VC,STAR,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   (g.InGrid() ? g.VCRank() : 0 ),0,
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VC,STAR,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   (g.InGrid() ? g.VCRank() : 0),0,
   (g.InGrid() ? LocalLength(height,g.VCRank(),0,g.Size()) : 0),width,
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VC,STAR,Int>::DistMatrix
( bool constrainedColAlignment, Int colAlignment, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.VCRank(),colAlignment,g.Size()) : 0),0,
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VC,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.VCRank(),colAlignment,g.Size()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.VCRank(),colAlignment,g.Size()) : 0),
   width,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VC,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
  Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.VCRank(),colAlignment,g.Size()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.VCRank(),colAlignment,g.Size()) : 0),
   width,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VC,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,colAlignment,0,
   (g.InGrid() ? Shift(g.VCRank(),colAlignment,g.Size()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.VCRank(),colAlignment,g.Size()) : 0),
   width,buffer,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,VC,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,colAlignment,0,
   (g.InGrid() ? Shift(g.VCRank(),colAlignment,g.Size()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.VCRank(),colAlignment,g.Size()) : 0),
   width,buffer,ldim,g)
{ }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,VC,STAR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  (A.Participating() ? A.ColRank() : 0),0,
  0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::DistMatrix");
#endif
    if( VC != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,VC,STAR,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [VC,* ] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,VC,STAR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetGrid( const elem::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->colAlignment_ = 0;
    if( grid.InGrid() )
        this->colShift_ = grid.VCRank();
    else
        this->colShift_ = 0;
}

template<typename T,typename Int>
inline Int
DistMatrix<T,VC,STAR,Int>::ColStride() const
{ return this->grid_->Size(); }

template<typename T,typename Int>
inline Int
DistMatrix<T,VC,STAR,Int>::RowStride() const
{ return 1; }

template<typename T,typename Int>
inline Int
DistMatrix<T,VC,STAR,Int>::ColRank() const
{ return this->grid_->VCRank(); }

template<typename T,typename Int>
inline Int
DistMatrix<T,VC,STAR,Int>::RowRank() const
{ return 0; }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->colAlignment_ = A.ColAlignment();
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift( g.VCRank(), this->ColAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->colAlignment_ = A.RowAlignment();
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift( g.VCRank(), this->ColAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignWith( const DistMatrix<S,MC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignWith([MC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->colAlignment_ = A.ColAlignment();
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift( g.VCRank(), this->ColAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignWith( const DistMatrix<S,STAR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignWith([* ,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->colAlignment_ = A.RowAlignment();
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift( g.VCRank(), this->ColAlignment(), g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignWith( const DistMatrix<S,VC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignWith([VC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->Empty();
    this->colAlignment_ = A.ColAlignment();
    this->constrainedColAlignment_ = true;
    if( this->Participating() )
        this->colShift_ = A.ColShift();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignWith( const DistMatrix<S,STAR,VC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignWith([* ,VC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->Empty();
    this->colAlignment_ = A.RowAlignment();
    this->constrainedColAlignment_ = true;
    if( this->Participating() )
        this->colShift_ = A.RowShift();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignColsWith( const DistMatrix<S,MC,MR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignColsWith( const DistMatrix<S,MR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignColsWith( const DistMatrix<S,MC,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignColsWith( const DistMatrix<S,STAR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignColsWith( const DistMatrix<S,VC,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignColsWith( const DistMatrix<S,STAR,VC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,VC,STAR,Int>::AlignedWithDiagonal
( const DistMatrix<S,VC,STAR,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignedWithDiagonal([VC,* ])");
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int p = g.Size();
    const Int colAlignment = A.ColAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRank = colAlignment;
        aligned = ( this->ColAlignment() == ownerRank );
    }
    else
    {
        const Int ownerRank = (colAlignment-offset) % p;
        aligned = ( this->ColAlignment() == ownerRank );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,VC,STAR,Int>::AlignedWithDiagonal
( const DistMatrix<S,STAR,VC,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignedWithDiagonal([* ,VC])");
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int p = g.Size();
    const Int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRank = (rowAlignment + offset) % p;
        aligned = ( this->ColAlignment() == ownerRank );
    }
    else
    {
        const Int ownerRank = rowAlignment;
        aligned = ( this->ColAlignment() == ownerRank );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignWithDiagonal
( const DistMatrix<S,VC,STAR,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignWithDiagonal([VC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int p = g.Size();
    const Int colAlignment = A.ColAlignment();

    this->Empty();
    if( offset >= 0 )
    {
        const Int ownerRank = colAlignment;
        this->colAlignment_ = ownerRank;
    }
    else
    {
        const Int ownerRank = (colAlignment-offset) % p;
        this->colAlignment_ = ownerRank;
    }
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift(g.VCRank(),this->colAlignment_,p);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,VC,STAR,Int>::AlignWithDiagonal
( const DistMatrix<S,STAR,VC,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignWithDiagonal([* ,VC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int p = g.Size();
    const Int rowAlignment = A.RowAlignment();

    this->Empty();
    if( offset >= 0 )
    {
        const Int ownerRank = (rowAlignment+offset) % p;
        this->colAlignment_ = ownerRank;
    }
    else
    {
        const Int ownerRank = rowAlignment;
        this->colAlignment_ = ownerRank;
    }
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift(g.VCRank(),this->colAlignment_,p);
#ifndef RELEASE
    PopCallStack();
#endif
}


template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::PrintBase");
#endif
    const elem::Grid& g = this->Grid();
    if( g.Rank() == 0 && msg != "" )
        os << msg << std::endl;

    const Int height      = this->Height();
    const Int width       = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int p           = g.Size();
    const Int colShift    = this->ColShift();

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
    for( Int j=0; j<width; ++j )
    {
        T* destCol = &sendBuf[colShift+j*height];
        const T* sourceCol = &thisLocalBuffer[j*thisLDim];
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
            destCol[iLocal*p] = sourceCol[iLocal];
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
                os << recvBuf[i+j*height] << " ";
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
DistMatrix<T,VC,STAR,Int>::Align( Int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::Align");
    this->AssertFreeColAlignment();
#endif
    this->AlignCols( colAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::AlignCols( Int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Size() )
        throw std::runtime_error("Invalid column alignment for [VC,* ]");
#endif
    this->colAlignment_ = colAlignment;
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift( g.VCRank(), colAlignment, g.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::Attach
( Int height, Int width, Int colAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::Attach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->viewing_ = true;
    if( g.InGrid() )
    {
        this->colShift_ = Shift(g.VCRank(),colAlignment,g.Size());
        const Int localHeight = LocalLength(height,this->colShift_,g.Size());
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
DistMatrix<T,VC,STAR,Int>::LockedAttach
( Int height, Int width, Int colAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->viewing_ = true;
    this->lockedView_ = true;
    if( g.InGrid() )
    {
        this->colShift_ = Shift(g.VCRank(),colAlignment,g.Size());
        const Int localHeight = LocalLength(height,this->colShift_,g.Size());
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
DistMatrix<T,VC,STAR,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    const elem::Grid& g = this->Grid();
    this->height_ = height;
    this->width_  = width;
    if( g.InGrid() )
        this->localMatrix_.ResizeTo
        ( LocalLength(height,this->ColShift(),g.Size()) ,width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline T
DistMatrix<T,VC,STAR,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
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
DistMatrix<T,VC,STAR,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->SetLocal(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->UpdateLocal(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::GetDiagonal
( DistMatrix<T,VC,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const Int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Height() || d.Width() != 1) )
    {
        std::ostringstream msg;
        msg << "d is not a column vec of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( (d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( diagLength, 1 );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift) / p;
    const Int localDiagLength = d.LocalHeight();
    T* dLocalBuffer = d.LocalBuffer();
    const T* thisLocalBuffer = this->LockedLocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        dLocalBuffer[k] = thisLocalBuffer[iLocal+jLocal*thisLDim];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::GetDiagonal
( DistMatrix<T,STAR,VC>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const Int diagLength = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (diagLength != d.Width() || d.Height() != 1) )
    {
        std::ostringstream msg;
        msg << "d is not a row vec of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << diagLength << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( 1, diagLength );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift) / p;
    const Int localDiagLength = d.LocalWidth();
    T* dLocalBuffer = d.LocalBuffer();
    const Int dLDim = d.LocalLDim();
    const T* thisLocalBuffer = this->LockedLocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        dLocalBuffer[k*dLDim] = thisLocalBuffer[iLocal+jLocal*thisLDim];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetDiagonal
( const DistMatrix<T,VC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetDiagonal");
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
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift)/p;
    const Int localDiagLength = d.LocalHeight();
    const T* dLocalBuffer = d.LockedLocalBuffer();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetDiagonal
( const DistMatrix<T,STAR,VC>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetDiagonal");
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
    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift)/p;
    const Int localDiagLength = d.LocalWidth();
    const T* dLocalBuffer = d.LockedLocalBuffer();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int dLDim = d.LocalLDim();
    const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k*dLDim];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., operator=
//

template<typename T,typename Int>
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            if( g.InGrid() )
                this->colShift_ = 
                    Shift( g.VCRank(), this->ColAlignment(), g.Size() );
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

    if( this->ColAlignment() % g.Height() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int colShiftOfA = A.ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLocalLength(height,p);
        const Int maxWidth = MaxLocalLength(width,c);
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
            const Int thisColShift = RawShift(thisRank,colAlignment,p); 
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
            {
                T* destCol = &data[jLocal*thisLocalHeight];
                const T* sourceCol = &ALocalBuffer[thisColOffset+jLocal*ALDim];
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    destCol[iLocal] = sourceCol[iLocal*c];
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

            const Int thisRowShift = RawShift(k,rowAlignmentOfA,c);
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,c);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(thisRowShift+jLocal*c)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
        this->auxMemory_.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [VC,* ] <- [MC,MR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int colShiftOfA = A.ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int rowAlignmentOfA = A.RowAlignment();
        
        const Int sendRow = (row+r+(colAlignment%r)-colAlignmentOfA) % r;
        const Int recvRow = (row+r+colAlignmentOfA-(colAlignment%r)) % r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localWidthOfA = A.LocalWidth();

        const Int maxHeight = MaxLocalLength(height,p);
        const Int maxWidth = MaxLocalLength(width,c);
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
            const Int thisColShift = RawShift(thisRank,colAlignment,p);
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = RawLocalLength(height,thisColShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int jLocal=0; jLocal<localWidthOfA; ++jLocal )
            {
                T* destCol = &data[jLocal*thisLocalHeight];
                const T* sourceCol = &ALocalBuffer[thisColOffset+jLocal*ALDim];
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    destCol[iLocal] = sourceCol[iLocal*c];
            }
        }

        // AllToAll to gather all of the unaligned [VC,*] data into firstBuffer
        mpi::AllToAll
        ( secondBuffer, portionSize, 
          firstBuffer,  portionSize, g.RowComm() );

        // SendRecv: properly align the [VC,*] via a trade in the column
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

            const Int thisRowShift = RawShift(k,rowAlignmentOfA,c);
            const Int thisLocalWidth = RawLocalLength(width,thisRowShift,c);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*localHeight];
                T* thisCol = &thisLocalBuffer[(thisRowShift+jLocal*c)*thisLDim];
                MemCopy( thisCol, dataCol, localHeight );
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
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            if( g.InGrid() )
                this->colShift_ = 
                    Shift( g.VCRank(), this->ColAlignment(), g.Size() );
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

    if( this->ColAlignment() % g.Height() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int colShift = this->ColShift();
        const Int colShiftOfA = A.ColShift();
        const Int colOffset = (colShift-colShiftOfA) / r;
        
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();

        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for 
#endif
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisLocalBuffer[j*thisLDim];
            const T* sourceCol = &ALocalBuffer[colOffset+j*ALDim];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                destCol[iLocal] = sourceCol[iLocal*c];
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [VC,* ] <- [MC,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = g.Size();
        const Int row = g.Row();
        const Int col = g.Col();
        const Int colShiftOfA = A.ColShift();
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[VC,*] within our process column to fix alignments.
        const Int sendRow = (row+r+(colAlignment%r)-colAlignmentOfA) % r;
        const Int recvRow = (row+r+colAlignmentOfA-(colAlignment%r)) % r;
        const Int sendRank = sendRow + r*col;

        const Int sendColShift = Shift( sendRank, colAlignment, p );
        const Int sendColOffset = (sendColShift-colShiftOfA) / r;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int localHeightOfSend = LocalLength(height,sendColShift,p);

        const Int sendSize = localHeightOfSend * width;
        const Int recvSize = localHeight * width;

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
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &sendBuffer[j*localHeightOfSend];
            const T* sourceCol = &ALocalBuffer[sendColOffset+j*ALDim];
            for( Int iLocal=0; iLocal<localHeightOfSend; ++iLocal )
                destCol[iLocal] = sourceCol[iLocal*c];
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
        for( Int j=0; j<width; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*localHeight];
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            MemCopy( thisCol, recvBufferCol, localHeight );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,MR]");
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
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ] = [MD,* ]");
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
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,MD]");
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
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VR,STAR,Int> A_VR_STAR(g);

    A_VR_STAR = A;
    *this = A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VR,STAR,Int> A_VR_STAR(g);

    A_VR_STAR = A;
    *this = A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MR,MC,Int> > A_MR_MC
    ( new DistMatrix<T,MR,MC,Int>(g) );
    *A_MR_MC = A;

    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            if( g.InGrid() )
                this->colShift_ = A.ColShift();
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

    if( this->ColAlignment() == A.ColAlignment() )
    {
        this->localMatrix_ = A.LockedLocalMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [VC,* ] <- [VC,* ]." << std::endl;
#endif
        const Int rank = g.VCRank();
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
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ALocalBuffer[j*ALDim];
            T* sendBufferCol = &sendBuffer[j*localHeightOfA];
            MemCopy( sendBufferCol, ACol, localHeightOfA );
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
        for( Int j=0; j<width; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*localHeight];
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            MemCopy( thisCol, recvBufferCol, localHeight );
        }
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MR,MC,Int> > A_MR_MC
    ( new DistMatrix<T,MR,MC,Int>(g) );
    *A_MR_MC = A;

    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VC,STAR,Int>(g) );
    *A_VR_STAR = *A_MR_MC; 
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_VR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [VR,* ]");
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
    
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int localHeightOfA = A.LocalHeight();

    const Int sendSize = localHeightOfA * width;
    const Int recvSize = localHeight * width;

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int rankCM = g.VCRank();
    const Int rankRM = g.VRRank();

    const Int colShift = this->ColShift();
    const Int colShiftOfA = A.ColShift();

    // Compute which colmajor rank has the colShift equal to our colShiftOfA
    const Int sendRankCM = (rankCM+(p+colShiftOfA-colShift)) % p;

    // Compute which colmajor rank has the A colShift that we need
    const Int recvRankRM = (rankRM+(p+colShift-colShiftOfA)) % p;
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
    for( Int j=0; j<width; ++j )
    {
        const T* ACol = &ALocalBuffer[j*ALDim];
        T* sendBufferCol = &sendBuffer[j*localHeightOfA];
        MemCopy( sendBufferCol, ACol, localHeightOfA );
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
    for( Int j=0; j<width; ++j )
    {
        const T* recvBufferCol = &recvBuffer[j*localHeight];
        T* thisCol = &thisLocalBuffer[j*thisLDim];
        MemCopy( thisCol, recvBufferCol, localHeight );
    }

    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,VR]");
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
inline const DistMatrix<T,VC,STAR,Int>&
DistMatrix<T,VC,STAR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[VC,* ] = [* ,* ]");
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
    const Int colShift = this->ColShift();

    const Int localHeight = this->LocalHeight();
    const Int width = this->Width();

    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for 
#endif
    for( Int j=0; j<width; ++j )
    {
        T* destCol = &thisLocalBuffer[j*thisLDim];
        const T* sourceCol = &ALocalBuffer[colShift+j*ALDim];
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
            destCol[iLocal] = sourceCol[iLocal*p];
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SumScatterFrom
( const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ]::SumScatterFrom( [MC,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && A.Grid().Rank() == 0 )
    {
        std::cerr <<
          "[VC,* ]::SumScatterFrom([MC,* ]) potentially causes a large amount "
          "of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MC,* ] matrix instead." << std::endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment();
            if( g.InGrid() )
                this->colShift_ = 
                    Shift( g.VCRank(), this->ColAlignment(), g.Size() );
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

    if( this->ColAlignment() % g.Height() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int colAlignment = this->ColAlignment();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int maxLocalHeight = MaxLocalLength( height, p );

        const Int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
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
            const Int thisColShift = RawShift( thisRank, colAlignment, p );
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = 
                RawLocalLength( height, thisColShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &data[j*thisLocalHeight];
                const T* sourceCol = &ALocalBuffer[thisColOffset+j*ALDim];
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    destCol[iLocal] = sourceCol[iLocal*c];
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
        for( Int j=0; j<width; ++j )
        {
            const T* bufferCol = &buffer[j*localHeight];
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            MemCopy( thisCol, bufferCol, localHeight );
        }
        this->auxMemory_.Release();
    }
    else
    {
        throw std::logic_error
        ("Unaligned [VC,* ]::ReduceScatterFrom( [MC,* ] ) not implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SumScatterFrom
( const DistMatrix<T,STAR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ]::SumScatterFrom( [* ,* ] )");
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
        return;
    }

    const Int p = g.Size();
    const Int colAlignment = this->ColAlignment();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int maxLocalHeight = MaxLocalLength( height, p );

    const Int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
    const Int sendSize = p*recvSize;

    this->auxMemory_.Require( sendSize );
    T* buffer = this->auxMemory_.Buffer();

    // Pack
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        T* data = &buffer[k*recvSize];

        const Int thisColShift = RawShift( k, colAlignment, p );
        const Int thisLocalHeight = RawLocalLength( height, thisColShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &data[j*thisLocalHeight];
            const T* sourceCol = &ALocalBuffer[thisColShift+j*ALDim];
            for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                destCol[iLocal] = sourceCol[iLocal*p];
        }
    }

    // Communicate
    mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.VCComm() );

    // Unpack our received data
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
    {
        const T* bufferCol = &buffer[j*localHeight];
        T* thisCol = &thisLocalBuffer[j*thisLDim];
        MemCopy( thisCol, bufferCol, localHeight );
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,MC,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ]::SumScatterUpdate( [MC,* ] )");
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

#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && A.Grid().Rank() == 0 )
    {
        std::cerr <<
          "[VC,* ]::SumScatterUpdate([MC,* ]) potentially causes a large amount"
          " of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MC,* ] matrix instead." << std::endl;
    }
#endif
    if( this->ColAlignment() % g.Height() == A.ColAlignment() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int myRow = g.Row();
        const Int colAlignment = this->ColAlignment();
        const Int colShiftOfA = A.ColShift();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = this->LocalHeight();
        const Int maxLocalHeight = MaxLocalLength( height, p );

        const Int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
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
            const Int thisColShift = RawShift( thisRank, colAlignment, p );
            const Int thisColOffset = (thisColShift-colShiftOfA) / r;
            const Int thisLocalHeight = 
                RawLocalLength( height, thisColShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
#endif
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &data[j*thisLocalHeight];
                const T* sourceCol = &ALocalBuffer[thisColOffset+j*ALDim];
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    destCol[iLocal] = sourceCol[iLocal*c];
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
        for( Int j=0; j<width; ++j )
        {
            const T* bufferCol = &buffer[j*localHeight];
            T* thisCol = &thisLocalBuffer[j*thisLDim];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                thisCol[iLocal] += alpha*bufferCol[iLocal];
        }
        this->auxMemory_.Release();
    }
    else
    {
        throw std::logic_error
        ("Unaligned [VC,* ]::ReduceScatterUpdate( [MC,* ] ) not implemented");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SumScatterUpdate
( T alpha, const DistMatrix<T,STAR,STAR,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VC,* ]::SumScatterUpdate( [* ,* ] )");
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

    const Int p = g.Size();
    const Int colAlignment = this->ColAlignment();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int maxLocalHeight = MaxLocalLength( height, p );

    const Int recvSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);
    const Int sendSize = p*recvSize;

    this->auxMemory_.Require( sendSize );
    T* buffer = this->auxMemory_.Buffer();

    // Pack
    const T* ALocalBuffer = A.LockedLocalBuffer();
    const Int ALDim = A.LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        T* data = &buffer[k*recvSize];

        const Int thisColShift = RawShift( k, colAlignment, p );
        const Int thisLocalHeight = RawLocalLength( height, thisColShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &data[j*thisLocalHeight];
            const T* sourceCol = &ALocalBuffer[thisColShift+j*ALDim];
            for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                destCol[iLocal] = sourceCol[iLocal*p];
        }
    }

    // Communicate
    mpi::ReduceScatter( buffer, recvSize, mpi::SUM, g.VCComm() );

    // Unpack our received data
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
    {
        const T* bufferCol = &buffer[j*localHeight];
        T* thisCol = &thisLocalBuffer[j*thisLDim];
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
            thisCol[iLocal] += alpha*bufferCol[iLocal];
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Routines which explicitly work in the complex plane
//

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,VC,STAR,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetRealPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        u = this->GetLocalRealPart(iLoc,j);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,VC,STAR,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetImagPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        u = this->GetLocalImagPart(iLoc,j);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetRealPart( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->SetLocalRealPart(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetImagPart( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->SetLocalImagPart(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::UpdateRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->UpdateLocalRealPart(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::UpdateImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-this->ColShift()) / g.Size();
        this->UpdateLocalImagPart(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::GetRealPartOfDiagonal
( DistMatrix<typename Base<T>::type,VC,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetRealPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const Int length = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (length != d.Height() || d.Width() != 1) )
    {
        std::ostringstream msg;
        msg << "d is not a column vec of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef typename Base<T>::type R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( length, 1 );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift) / p;
    const Int localDiagLength = d.LocalHeight();
    const T* thisLocalBuffer = this->LockedLocalBuffer();
    const Int thisLDim = this->LocalLDim();
    R* dLocalBuffer = d.LocalBuffer();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        dLocalBuffer[k] = RealPart( thisLocalBuffer[iLocal+jLocal*thisLDim] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::GetImagPartOfDiagonal
( DistMatrix<typename Base<T>::type,VC,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetImagPartOfDiagonal");
    if( d.Viewing() )
        this->AssertSameGrid( d );
#endif
    const Int length = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (length != d.Height() || d.Width() != 1) )
    {
        std::ostringstream msg;
        msg << "d is not a column vec of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef typename Base<T>::type R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( length, 1 );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift) / p;
    const Int localDiagLength = d.LocalHeight();
    const T* thisLocalBuffer = this->LockedLocalBuffer();
    const Int thisLDim = this->LocalLDim();
    R* dLocalBuffer = d.LocalBuffer();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        dLocalBuffer[k] = ImagPart( thisLocalBuffer[iLocal+jLocal*thisLDim] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::GetRealPartOfDiagonal
( DistMatrix<typename Base<T>::type,STAR,VC,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetRealPartOfDiagonal");
    if( d.Viewing() )
        AssertSameGrid( d );
#endif
    const Int length = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef typename Base<T>::type R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( 1, length );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift) / p;
    const Int localDiagLength = d.LocalWidth();
    const T* thisLocalBuffer = this->LockedLocalBuffer();
    const Int thisLDim = this->LocalLDim();
    R* dLocalBuffer = d.LocalBuffer();
    const Int dLDim = d.LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        dLocalBuffer[k*dLDim] = 
            RealPart( thisLocalBuffer[iLocal+jLocal*thisLDim] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::GetImagPartOfDiagonal
( DistMatrix<typename Base<T>::type,STAR,VC,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetImagPartOfDiagonal");
    if( d.Viewing() )
        AssertSameGrid( d );
#endif
    const Int length = this->DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    typedef typename Base<T>::type R;

    const elem::Grid& g = this->Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( *this, offset );
        d.ResizeTo( 1, length );
    }
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift) / p;
    const Int localDiagLength = d.LocalWidth();
    const T* thisLocalBuffer = this->LockedLocalBuffer();
    const Int thisLDim = this->LocalLDim();
    R* dLocalBuffer = d.LocalBuffer();
    const Int dLDim = d.LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        dLocalBuffer[k*dLDim] = 
            ImagPart( thisLocalBuffer[iLocal+jLocal*thisLDim] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetRealPartOfDiagonal
( const DistMatrix<typename Base<T>::type,VC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetRealPartOfDiagonal");
    AssertSameGrid( d );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    typedef typename Base<T>::type R;

    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift)/p;
    const Int localDiagLength = d.LocalHeight();
    const R* dLocalBuffer = d.LockedLocalBuffer();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        this->SetLocalRealPart( iLocal, jLocal, dLocalBuffer[k] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetImagPartOfDiagonal
( const DistMatrix<typename Base<T>::type,VC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetImagPartOfDiagonal");
    AssertSameGrid( d );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    typedef typename Base<T>::type R;
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift)/p;
    const Int localDiagLength = d.LocalHeight();
    const R* dLocalBuffer = d.LockedLocalBuffer();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        this->SetLocalImagPart( iLocal, jLocal, dLocalBuffer[k] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetRealPartOfDiagonal
( const DistMatrix<typename Base<T>::type,STAR,VC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetRealPartOfDiagonal");
    AssertSameGrid( d );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    typedef typename Base<T>::type R;

    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift)/p;
    const Int localDiagLength = d.LocalWidth();
    const R* dLocalBuffer = d.LockedLocalBuffer();
    const Int dLDim = d.LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        this->SetLocalRealPart( iLocal, jLocal, dLocalBuffer[k*dLDim] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetImagPartOfDiagonal
( const DistMatrix<typename Base<T>::type,STAR,VC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetImagPartOfDiagonal");
    AssertSameGrid( d );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const Int length = this->DiagonalLength(offset);
    if( length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << this->Height() << " x " << this->Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( *this, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    typedef typename Base<T>::type R;
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int p = g.Size();
    const Int colShift = this->ColShift();
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

    const Int iLocalStart = (iStart-colShift)/p;
    const Int localDiagLength = d.LocalWidth();
    const R* dLocalBuffer = d.LockedLocalBuffer();
    const Int dLDim = d.LocalLDim();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLocal = iLocalStart+k;
        const Int jLocal = jStart+k*p;
        this->SetLocalImagPart( iLocal, jLocal, dLocalBuffer[k*dLDim] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
