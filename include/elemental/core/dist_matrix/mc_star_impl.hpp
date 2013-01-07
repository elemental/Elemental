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
DistMatrix<T,MC,STAR,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   (g.InGrid() ? g.Row() : 0),0,
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   (g.InGrid() ? g.Row() : 0),0,
   (g.InGrid() ? LocalLength(height,g.Row(),0,g.Height()) : 0),width,
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( bool constrainedColAlignment, Int colAlignment, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.Row(),colAlignment,g.Height()) : 0),0,
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.Row(),colAlignment,g.Height()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.Row(),colAlignment,g.Height()) : 0),
   width,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
  Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.Row(),colAlignment,g.Height()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.Row(),colAlignment,g.Height()) : 0),
   width,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, const T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,colAlignment,0,
   (g.InGrid() ? Shift(g.Row(),colAlignment,g.Height()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.Row(),colAlignment,g.Height()) : 0),
   width,buffer,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, T* buffer, Int ldim,
  const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,colAlignment,0,
   (g.InGrid() ? Shift(g.Row(),colAlignment,g.Height()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.Row(),colAlignment,g.Height()) : 0),
   width,buffer,ldim,g)
{ }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,
  (A.Participating() ? A.ColRank() : 0),0,
  0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::DistMatrix");
#endif
    if( MC != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,MC,STAR,Int>*>(&A) != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [MC,* ] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetGrid( const elem::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->colAlignment_ = 0;
    if( grid.InGrid() )
        this->colShift_ = grid.Row();
    else
        this->colShift_ = 0;
}

template<typename T,typename Int>
inline Int
DistMatrix<T,MC,STAR,Int>::ColStride() const
{ return this->grid_->Height(); }

template<typename T,typename Int>
inline Int
DistMatrix<T,MC,STAR,Int>::RowStride() const
{ return 1; }

template<typename T,typename Int>
inline Int
DistMatrix<T,MC,STAR,Int>::ColRank() const
{ return this->grid_->Row(); }

template<typename T,typename Int>
inline Int
DistMatrix<T,MC,STAR,Int>::RowRank() const
{ return 0; }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([MC,MR])");
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
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,MC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([MC,* ])");
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
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([MR,MC])");
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
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,STAR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([* ,MC])");
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
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,VC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([VC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->colAlignment_ = A.ColAlignment() % g.Height();
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift( g.Row(), this->ColAlignment(), g.Height() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,STAR,VC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([* ,VC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    this->Empty();
    this->colAlignment_ = A.RowAlignment() % g.Height();
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift( g.Row(), this->ColAlignment(), g.Height() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,MC,MR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,MC,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,MR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,STAR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,VC,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,STAR,VC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,MC,STAR,Int>::AlignedWithDiagonal
( const DistMatrix<S,MC,STAR,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignedWithDiagonal([MC,* ])");
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int colAlignment = A.ColAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        aligned = ( this->ColAlignment() == ownerRow );
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        aligned = ( this->ColAlignment() == ownerRow );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,MC,STAR,Int>::AlignedWithDiagonal
( const DistMatrix<S,STAR,MC,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignedWithDiagonal([* ,MC])");
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRow = (rowAlignment + offset) % r;
        aligned = ( this->ColAlignment() == ownerRow );
    }
    else
    {
        const Int ownerRow = rowAlignment;
        aligned = ( this->ColAlignment() == ownerRow );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignWithDiagonal
( const DistMatrix<S,MC,STAR,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWithDiagonal([MC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int colAlignment = A.ColAlignment();

    this->Empty();
    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        this->colAlignment_ = ownerRow;
    }
    else 
    {
        const Int ownerRow = (colAlignment-offset) % r;
        this->colAlignment_ = ownerRow;
    }
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift(g.Row(),this->colAlignment_,r);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignWithDiagonal
( const DistMatrix<S,STAR,MC,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWithDiagonal([* ,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elem::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int rowAlignment = A.RowAlignment();

    this->Empty();
    if( offset >= 0 )
    {
        const Int ownerRow = (rowAlignment+offset) % r;
        this->colAlignment_ = ownerRow;
    }
    else
    {
        const Int ownerRow = rowAlignment;
        this->colAlignment_ = ownerRow;
    }
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift(g.Row(),this->colAlignment_,r);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::PrintBase");
#endif
    const elem::Grid& g = this->Grid();
    if( g.Rank() == 0 && msg != "" )
        os << msg << std::endl;
        
    const Int height      = this->Height();
    const Int width       = this->Width();
    const Int localHeight = this->LocalHeight();
    const Int r           = g.Height();
    const Int colShift    = this->ColShift();

    if( height == 0 || width == 0 || !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Only one process column needs to participate
    if( g.Col() == 0 )
    {
        std::vector<T> sendBuf(height*width,0);
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const Int thisLDim = this->LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for 
#endif

        for( Int jLocal=0; jLocal<width; ++jLocal )
        {
            T* destCol = &sendBuf[colShift+jLocal*height];
            const T* sourceCol = &thisLocalBuffer[jLocal*thisLDim];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                destCol[iLocal*r] = sourceCol[iLocal];
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
DistMatrix<T,MC,STAR,Int>::Align( Int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Align");
    this->AssertFreeColAlignment();
#endif
    this->AlignCols( colAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::AlignCols( Int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const elem::Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Height() )
        throw std::logic_error("Invalid column alignment for [MC,* ]");
#endif
    this->Empty();
    this->colAlignment_ = colAlignment;
    this->constrainedColAlignment_ = true;
    if( g.InGrid() )
        this->colShift_ = Shift( g.Row(), colAlignment, g.Height() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::Attach
( Int height, Int width, Int colAlignment,
  T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Attach");
#endif
    this->Empty();

    this->grid_ = &g;
    this->height_ = height;
    this->width_ = width;
    this->colAlignment_ = colAlignment;
    this->viewing_ = true;
    if( g.InGrid() )
    {
        this->colShift_ = Shift(g.Row(),colAlignment,g.Height());
        const Int localHeight = LocalLength(height,this->colShift_,g.Height());
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
DistMatrix<T,MC,STAR,Int>::LockedAttach
( Int height, Int width, Int colAlignment,
  const T* buffer, Int ldim, const elem::Grid& g )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedAttach");
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
        this->colShift_ = Shift(g.Row(),colAlignment,g.Height());
        const Int localHeight = LocalLength(height,this->colShift_,g.Height());
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
DistMatrix<T,MC,STAR,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->localMatrix_.ResizeTo
        ( LocalLength(height,this->ColShift(),this->Grid().Height()), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline T DistMatrix<T,MC,STAR,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Get");
    this->AssertValidEntry( i, j );
    if( !this->Participating() )
        throw std::logic_error("Should only be called by grid members");
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that
    // row within each process column
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();

    T u;
    if( g.Row() == ownerRow )
    {
        const Int iLocal = (i-this->ColShift()) / g.Height();
        u = this->GetLocal( iLocal, j );
    }
    mpi::Broadcast( &u, 1, ownerRow, g.ColComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void DistMatrix<T,MC,STAR,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        this->SetLocal(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void DistMatrix<T,MC,STAR,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i + this->ColAlignment()) % g.Height();

    if( g.Row() == ownerRow )
    {
        const Int iLoc = (i-this->ColShift()) / g.Height();
        this->UpdateLocal(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetDiagonal
( DistMatrix<T,MC,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetDiagonal");
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

    const Int r = g.Height();
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

    const Int iLocalStart = (iStart-colShift) / r;
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
        const Int jLocal = jStart+k*r;
        dLocalBuffer[k] = thisLocalBuffer[iLocal+jLocal*thisLDim];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetDiagonal
( DistMatrix<T,STAR,MC,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetDiagonal");
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

    const Int r = g.Height();
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

    const Int iLocalStart = (iStart-colShift) / r;
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
        const Int jLocal = jStart+k*r;
        dLocalBuffer[k*dLDim] = thisLocalBuffer[iLocal+jLocal*thisLDim];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetDiagonal
( const DistMatrix<T,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetDiagonal");
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

    const Int r = g.Height();
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

    const Int iLocalStart = (iStart-colShift)/r;
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
        const Int jLocal = jStart+k*r;
        thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetDiagonal
( const DistMatrix<T,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetDiagonal");
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

    const Int r = g.Height();
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

    const Int iLocalStart = (iStart-colShift)/r;
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
        const Int jLocal = jStart+k*r;
        thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k*dLDim];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., SumOverRow
//

template<typename T,typename Int>
inline void DistMatrix<T,MC,STAR,Int>::SumOverRow()
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SumOverRow");
    this->AssertNotLockedView();
#endif
    if( !this->Participating() )
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
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MC,MR]");
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
                    Shift( g.Row(), this->ColAlignment(), this->ColStride() );
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
        if( A.Width() == 1 )
        {
            if( g.Col() == A.RowAlignment() )
                this->localMatrix_ = A.LockedLocalMatrix();

            // Communicate
            mpi::Broadcast
            ( this->localMatrix_.Buffer(), this->LocalHeight(),
              A.RowAlignment(), g.RowComm() );
        }
        else
        {
            const Int c = g.Width();

            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localWidthOfA = A.LocalWidth();
            const Int maxLocalWidth = MaxLocalLength(width,c);

            const Int portionSize = 
                std::max(localHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

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
                T* originalDataCol = &originalData[jLocal*localHeight];
                MemCopy( originalDataCol, ACol, localHeight );
            }

            // Communicate
            mpi::AllGather
            ( originalData, portionSize,
              gatheredData, portionSize, g.RowComm() );

            // Unpack
            const Int rowAlignmentOfA = A.RowAlignment();
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                const T* data = &gatheredData[k*portionSize];

                const Int rowShift = RawShift( k, rowAlignmentOfA, c );
                const Int localWidth = RawLocalLength( width, rowShift, c );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    const T* dataCol = &data[jLocal*localHeight];
                    T* thisCol = 
                        &thisLocalBuffer[(rowShift+jLocal*c)*thisLDim];
                    MemCopy( thisCol, dataCol, localHeight );
                }
            }
            this->auxMemory_.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [MC,* ] <- [MC,MR]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int row = g.Row();

        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
        const Int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

        if( A.Width()==1 )
        {
            const Int localHeight = this->LocalHeight();

            if( this->grid_->Col() == A.RowAlignment() )
            {
                const Int localHeightOfA = A.LocalHeight();

                this->auxMemory_.Require( localHeightOfA );
                T* buffer = this->auxMemory_.Buffer();

                // Pack
                const T* ACol = A.LockedLocalBuffer(0,0);
                MemCopy( buffer, ACol, localHeightOfA );

                // Communicate
                mpi::SendRecv
                ( buffer, localHeightOfA, sendRow, 0,
                  this->localMatrix_.Buffer(), localHeight, recvRow, 
                  mpi::ANY_TAG, g.ColComm() );

                this->auxMemory_.Release();
            }

            // Communicate
            mpi::Broadcast
            ( this->localMatrix_.Buffer(), localHeight, A.RowAlignment(),
              g.RowComm() );
        }
        else
        {
            const Int height = this->Height();
            const Int width = this->Width();
            const Int localHeight = this->LocalHeight();
            const Int localHeightOfA = A.LocalHeight();
            const Int localWidthOfA  = A.LocalWidth();
            const Int maxLocalHeight = MaxLocalLength(height,r);
            const Int maxLocalWidth  = MaxLocalLength(width,c);

            const Int portionSize = 
                std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

            this->auxMemory_.Require( (c+1)*portionSize );

            T* buffer = this->auxMemory_.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[portionSize];

            // Pack the currently owned local data of A into the second 
            // buffer
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
              firstBuffer,  portionSize, recvRow, mpi::ANY_TAG, 
              g.ColComm() );

            // Use the output of the SendRecv as the input to the AllGather
            mpi::AllGather
            ( firstBuffer,  portionSize, 
              secondBuffer, portionSize, g.RowComm() );

            // Unpack the contents of each member of the process row
            const Int rowAlignmentOfA = A.RowAlignment();
            T* thisLocalBuffer = this->LocalBuffer();
            const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int k=0; k<c; ++k )
            {
                const T* data = &secondBuffer[k*portionSize];

                const Int rowShift = RawShift( k, rowAlignmentOfA, c ); 
                const Int localWidth = RawLocalLength( width, rowShift, c );
#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    const T* dataCol = &data[jLocal*localHeight];
                    T* thisCol = 
                        &thisLocalBuffer[(rowShift+jLocal*c)*thisLDim];
                    MemCopy( thisCol, dataCol, localHeight );
                }
            }
            this->auxMemory_.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MC,* ]");
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
            std::cerr << "Unaligned [MC,* ] <- [MC,* ]." << std::endl;
#endif
        const Int rank = g.Row();
        const Int r = g.Height();

        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();

        const Int sendRank = (rank+r+colAlignment-colAlignmentOfA) % r;
        const Int recvRank = (rank+r+colAlignmentOfA-colAlignment) % r;

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
          recvBuffer, recvSize, recvRank, mpi::ANY_TAG, g.ColComm() );

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
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR,Int> A_MC_MR(true,false,this->ColAlignment(),0,g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MC,* ] = [MD,* ] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw std::logic_error("[MC,* ] = [* ,MD] not yet implemented");
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = A;

    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > A_VR_STAR
    ( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = A;

    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > A_VC_STAR
    ( new DistMatrix<T,VC,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,MR,MC,Int> > 
        A_MR_MC( new DistMatrix<T,MR,MC,Int>(g) );
    *A_MR_MC = A;

    std::auto_ptr<DistMatrix<T,VR,STAR,Int> > 
        A_VR_STAR( new DistMatrix<T,VR,STAR,Int>(g) );
    *A_VR_STAR = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    std::auto_ptr<DistMatrix<T,VC,STAR,Int> > 
        A_VC_STAR( new DistMatrix<T,VC,STAR,Int>(true,this->ColAlignment(),g) );
    *A_VC_STAR = *A_VR_STAR;
    delete A_VR_STAR.release(); // lowers memory highwater

    *this = *A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Width() == 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "The vector version of [MC,* ] <- [VC,* ] is not yet written, but"
          " it only requires a modification of the vector version of "
          "[* ,MR] <- [* ,VR]" << std::endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.Rank() == 0 )
    {
        std::cerr << 
          "[MC,* ] <- [VC,* ] potentially causes a large amount of cache-"
          "thrashing. If possible avoid it by performing the redistribution"
          " with a (conjugate)-transpose: \n"
          "  [* ,MC].TransposeFrom([VC,* ])" << std::endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->colAlignment_ = A.ColAlignment() % g.Height();
            if( g.InGrid() )
                this->colShift_ = 
                    Shift( g.Row(), this->ColAlignment(), g.Height() );
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

    if( this->ColAlignment() == A.ColAlignment() % g.Height() )
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int row = g.Row();

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLocalLength(height,p);

        const Int portionSize = 
            std::max(maxLocalHeightOfA*width,mpi::MIN_COLL_MSG);

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
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ALocalBuffer[j*ALDim];
            T* originalDataCol = &originalData[j*localHeightOfA];
            MemCopy( originalDataCol, ACol, localHeightOfA );
        }

        // Communicate 
        mpi::AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.RowComm() );

        // Unpack
        const Int colShift = this->ColShift();
        const Int colAlignmentOfA = A.ColAlignment();
        T* thisLocalBuffer = this->LocalBuffer();
        const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];    

            const Int colShiftOfA = RawShift( row+r*k, colAlignmentOfA, p );
            const Int colOffset = (colShiftOfA-colShift) / r;
            const Int localHeight = RawLocalLength( height, colShiftOfA, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &thisLocalBuffer[colOffset+j*thisLDim];
                const T* sourceCol = &data[j*localHeight];
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
            std::cerr << "Unaligned [MC,* ] <- [VC,* ]." << std::endl;
#endif
        const Int r = g.Height();
        const Int c = g.Width();
        const Int p = r * c;
        const Int row = g.Row();
        const Int rank = g.VCRank();

        // Perform the SendRecv to make A have the same colAlignment
        const Int colAlignment = this->ColAlignment();
        const Int colAlignmentOfA = A.ColAlignment();
        const Int colShift = this->ColShift();

        const Int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const Int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeightOfA = MaxLocalLength(height,p);

        const Int portionSize = 
            std::max(maxLocalHeightOfA*width,mpi::MIN_COLL_MSG);

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
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ALocalBuffer[j*ALDim];
            T* secondBufferCol = &secondBuffer[j*localHeightOfA];
            MemCopy( secondBufferCol, ACol, localHeightOfA );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        mpi::SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, mpi::ANY_TAG, g.VCComm() );

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

            const Int colShiftOfA = RawShift( row+r*k, colAlignment, p );
            const Int colOffset = (colShiftOfA-colShift) / r;
            const Int localHeight = RawLocalLength( height, colShiftOfA, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &thisLocalBuffer[colOffset+j*thisLDim];
                const T* sourceCol = &data[j*localHeight];
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
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    std::auto_ptr<DistMatrix<T,STAR,VR,Int> > 
        A_STAR_VR( new DistMatrix<T,STAR,VR,Int>(g) );
    *A_STAR_VR = A;

    std::auto_ptr<DistMatrix<T,MC,MR,Int> > 
        A_MC_MR
        ( new DistMatrix<T,MC,MR,Int>(true,false,this->ColAlignment(),0,g) );
    *A_MC_MR = *A_STAR_VR;
    delete A_STAR_VR.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,VC,STAR,Int> A_VC_STAR(true,this->ColAlignment(),g);

    A_VC_STAR = A;
    *this = A_VC_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const elem::Grid& g = this->Grid();
    DistMatrix<T,MC,MR,Int> A_MC_MR(true,false,this->ColAlignment(),0,g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,MC,STAR,Int>&
DistMatrix<T,MC,STAR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Int r = this->Grid().Height(); 
    const Int colShift = this->ColShift();

    const Int localHeight = this->LocalHeight();
    const Int width = this->Width();

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
            destCol[iLocal] = sourceCol[iLocal*r];
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
DistMatrix<T,MC,STAR,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealPart");
    this->AssertValidEntry( i, j );
    if( !this->Participating() )
        throw std::logic_error("Should only be called by grid members");
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner row of entry (i,j) and broadcast from that
    // row within each process column
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();

    R u;
    if( g.Row() == ownerRow )
    {
        const Int iLocal = (i-this->ColShift()) / g.Height();
        u = this->GetLocalRealPart( iLocal, j );
    }
    mpi::Broadcast( &u, 1, ownerRow, g.ColComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,MC,STAR,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImagPart");
    this->AssertValidEntry( i, j );
    if( !this->Participating() )
        throw std::logic_error("Should only be called by grid members");
#endif
    typedef typename Base<T>::type R;

    // We will determine the owner row of entry (i,j) and broadcast from that
    // row within each process column
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();

    R u;
    if( g.Row() == ownerRow )
    {
        const Int iLocal = (i-this->ColShift()) / g.Height();
        u = this->GetLocalImagPart( iLocal, j );
    }
    mpi::Broadcast( &u, 1, ownerRow, g.ColComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealPart( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetRealPart");
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int iLocal = (i-this->ColShift()) / g.Height();
        this->SetLocalRealPart( iLocal, j, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagPart( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImagPart");
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int iLocal = (i-this->ColShift()) / g.Height();
        this->SetLocalImagPart( iLocal, j, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::UpdateRealPart");
#endif
    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int iLocal = (i-this->ColShift()) / g.Height();
        this->UpdateLocalRealPart( iLocal, j, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::UpdateImagPart");
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");

    const elem::Grid& g = this->Grid();
    const Int ownerRow = (i+this->ColAlignment()) % g.Height();
    if( g.Row() == ownerRow )
    {
        const Int iLocal = (i-this->ColShift()) / g.Height();
        this->UpdateLocalImagPart( iLocal, j, u );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealPartOfDiagonal
( DistMatrix<typename Base<T>::type,MC,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealPartOfDiagonal");
    if( d.Viewing() )
        AssertSameGrid( d );
#endif
    const Int length = this->DiagonalLength( offset );
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

    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift) / r;
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
            const Int jLocal = jStart+k*r;
            dLocalBuffer[k] = RealPart(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagPartOfDiagonal
( DistMatrix<typename Base<T>::type,MC,STAR,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImagPartOfDiagonal");
    if( d.Viewing() )
        AssertSameGrid( d );
#endif
    const Int length = this->DiagonalLength( offset );
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

    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift) / r;
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
            const Int jLocal = jStart+k*r;
            dLocalBuffer[k] = ImagPart(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealPartOfDiagonal
( DistMatrix<typename Base<T>::type,STAR,MC,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealPartOfDiagonal");
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

    if( g.InGrid() )
    {
        const Int r = g.Height();
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
        const Int iLocalStart = (iStart-colShift) / r;
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
            const Int jLocal = jStart+k*r;
            dLocalBuffer[k*dLDim] =
                RealPart(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagPartOfDiagonal
( DistMatrix<typename Base<T>::type,STAR,MC,Int>& d, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImagPartOfDiagonal");
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

    if( g.InGrid() )
    {
        const Int r = g.Height();
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
        const Int iLocalStart = (iStart-colShift) / r;
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
            const Int jLocal = jStart+k*r;
            dLocalBuffer[k*dLDim] =
                ImagPart(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealPartOfDiagonal
( const DistMatrix<typename Base<T>::type,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetRealPartOfDiagonal");
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
    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift)/r;
        const Int localDiagLength = d.LocalHeight();
        const R* dLocalBuffer = d.LockedLocalBuffer();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            this->SetLocalRealPart( iLocal, jLocal, dLocalBuffer[k] );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagPartOfDiagonal
( const DistMatrix<typename Base<T>::type,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImagPartOfDiagonal");
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
    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift)/r;
        const Int localDiagLength = d.LocalHeight();
        const R* dLocalBuffer = d.LockedLocalBuffer();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            this->SetLocalImagPart( iLocal, jLocal, dLocalBuffer[k] );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealPartOfDiagonal
( const DistMatrix<typename Base<T>::type,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetRealPartOfDiagonal");
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
    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift)/r;
        const Int localDiagLength = d.LocalWidth();

        const R* dLocalBuffer = d.LockedLocalBuffer();
        const Int dLDim = d.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            this->SetLocalRealPart( iLocal, jLocal, dLocalBuffer[k*dLDim] );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagPartOfDiagonal
( const DistMatrix<typename Base<T>::type,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImagPartOfDiagonal");
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
    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift)/r;
        const Int localDiagLength = d.LocalWidth();

        const R* dLocalBuffer = d.LockedLocalBuffer();
        const Int dLDim = d.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            this->SetLocalImagPart( iLocal, jLocal, dLocalBuffer[k*dLDim] );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
