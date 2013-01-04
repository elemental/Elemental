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
DistMatrix<T,STAR,STAR,Int>::DistMatrix( const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,0,0,0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,STAR,Int>::DistMatrix
( Int height, Int width, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,0,0,height,width,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,STAR,Int>::DistMatrix
( Int height, Int width, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,0,0,height,width,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,STAR,Int>::DistMatrix
( Int height, Int width, const T* buffer, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,0,0,0,height,width,buffer,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,STAR,Int>::DistMatrix
( Int height, Int width, T* buffer, Int ldim, const elem::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,0,0,0,height,width,buffer,ldim,g)
{ }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,STAR,STAR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::DistMatrix");
#endif
    if( STAR != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,STAR,STAR,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,* ] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,STAR,STAR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::SetGrid( const elem::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
}

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,STAR,Int>::ColStride() const
{ return 1; }

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,STAR,Int>::RowStride() const
{ return 1; }

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,STAR,Int>::ColRank() const
{ return 0; }

template<typename T,typename Int>
inline Int
DistMatrix<T,STAR,STAR,Int>::RowRank() const
{ return 0; }

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::PrintBase");
#endif
    const elem::Grid& g = this->Grid();
    if( g.Rank() == 0 && msg != "" )
        os << msg << std::endl;

    const Int height = this->Height();
    const Int width  = this->Width();

    if( height == 0 || width == 0 || !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( g.Rank() == 0 )
    {
        for( Int i=0; i<height; ++i )
        {
            for( Int j=0; j<width; ++j )
                os << this->GetLocal(i,j) << " ";
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
DistMatrix<T,STAR,STAR,Int>::Attach
( Int height, Int width, 
  T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Attach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    if( this->Participating() )
        this->localMatrix_.Attach( height, width, buffer, ldim );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::LockedAttach
( Int height, Int width, 
  const T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedAttach");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Participating() )
        this->localMatrix_.LockedAttach( height, width, buffer, ldim );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    this->height_ = height;
    this->width_ = width;
    if( this->Participating() )
        this->localMatrix_.ResizeTo( height, width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline T
DistMatrix<T,STAR,STAR,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    const Int viewingSize = mpi::CommSize( this->Grid().ViewingComm() );
    const Int owningSize = mpi::GroupSize( this->Grid().OwningGroup() );
    T u;
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        u = this->GetLocal(i,j);
    }
    else
    {
        // Have the root broadcast its data
        if( this->Grid().VCRank() == 0 )
            u = this->GetLocal(i,j);
        mpi::Broadcast
        ( &u, 1, this->Grid().VCToViewingMap(0), 
          this->Grid().ViewingComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::Set( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->SetLocal(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::Update( Int i, Int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->UpdateLocal(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., operator=
//

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,MC,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MC,MR]");
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

    const Int r = g.Height();
    const Int c = g.Width(); 
    const Int p = g.Size();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalHeight = MaxLocalLength(height,r);
    const Int maxLocalWidth = MaxLocalLength(width,c);

    const Int portionSize = 
        std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

    this->auxMemory_.Require( (p+1)*portionSize );

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
        T* originalDataCol = &originalData[jLocal*localHeightOfA];
        MemCopy( originalDataCol, ACol, localHeightOfA );
    }

    // Communicate
    mpi::AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, g.VCComm() );

    // Unpack
    const Int colAlignmentOfA = A.ColAlignment();
    const Int rowAlignmentOfA = A.RowAlignment();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int l=0; l<c; ++l )
    {
        const Int rowShift = RawShift( l, rowAlignmentOfA, c );
        const Int localWidth = RawLocalLength( width, rowShift, c );

        for( Int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[(k+l*r)*portionSize];

            const Int colShift = RawShift( k, colAlignmentOfA, r );
            const Int localHeight = RawLocalLength( height, colShift, r );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                T* destCol = 
                    &thisLocalBuffer[colShift+(rowShift+jLocal*c)*thisLDim];
                const T* sourceCol = &data[jLocal*localHeight];
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    destCol[iLocal*r] = sourceCol[iLocal];
            }
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,MC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MC,* ]");
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

    const Int r = g.Height();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int maxLocalHeight = MaxLocalLength(height,r);

    const Int portionSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

    this->auxMemory_.Require( (r+1)*portionSize );

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
      gatheredData, portionSize, g.ColComm() );

    // Unpack
    const Int colAlignmentOfA = A.ColAlignment();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<r; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const Int colShift = RawShift( k, colAlignmentOfA, r );
        const Int localHeight = RawLocalLength( height, colShift, r );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisLocalBuffer[colShift+j*thisLDim];
            const T* sourceCol = &data[j*localHeight];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                destCol[iLocal*r] = sourceCol[iLocal];
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,MR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,MR]");
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

    const Int c = g.Width();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalWidth = MaxLocalLength(width,c);

    const Int portionSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

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
            const T* dataCol = &data[jLocal*height];
            T* thisCol = &thisLocalBuffer[(rowShift+jLocal*c)*thisLDim];
            MemCopy( thisCol, dataCol, height );
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,MD,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MD,* ]");
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

    const Int p = g.Size();
    const Int lcm = g.LCM();
    const Int ownerPath = A.diagPath_;
    const Int ownerPathRank = A.colAlignment_;

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeight = A.LocalHeight();
    const Int maxLocalHeight = MaxLocalLength( height, lcm );
    const Int portionSize = std::max( maxLocalHeight*width, mpi::MIN_COLL_MSG );

    // Since a MD communicator has not been implemented, we will take
    // the suboptimal route of 'rounding up' everyone's contribution over 
    // the VC communicator.
    this->auxMemory_.Require( (p+1)*portionSize );
    T* buffer = this->auxMemory_.Buffer();
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    if( A.Participating() )
    {
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<width; ++j )
        {
            const T* ACol = &ALocalBuffer[j*ALDim];
            T* sendBufCol = &sendBuf[j*localHeight];
            MemCopy( sendBufCol, ACol, localHeight );
        }
    }

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.VCComm() );

    // Unpack
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        if( g.DiagPath( k ) == ownerPath )
        {
            const T* data = &recvBuf[k*portionSize];

            const Int thisPathRank = g.DiagPathRank( k );
            const Int thisColShift = 
                RawShift( thisPathRank, ownerPathRank, lcm );
            const Int thisLocalHeight = 
                RawLocalLength( height, thisColShift, lcm );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                T* destCol = &thisLocalBuffer[thisColShift+j*thisLDim];
                const T* sourceCol = &data[j*thisLocalHeight];
                for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                    destCol[iLocal*lcm] = sourceCol[iLocal];
            }
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,MD,Int>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,MD]");
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

    const Int p = g.Size();
    const Int lcm = g.LCM();
    const Int ownerPath = A.diagPath_;
    const Int ownerPathRank = A.rowAlignment_;

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localWidth = A.LocalWidth();
    const Int maxLocalWidth = MaxLocalLength( width, lcm );
    const Int portionSize = std::max( height*maxLocalWidth, mpi::MIN_COLL_MSG );

    // Since a MD communicator has not been implemented, we will take
    // the suboptimal route of 'rounding up' everyone's contribution over 
    // the VC communicator.
    this->auxMemory_.Require( (p+1)*portionSize );
    T* buffer = this->auxMemory_.Buffer();
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    if( A.Participating() )
    {
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[jLocal*ALDim];
            T* sendBufCol = &sendBuf[jLocal*height];
            MemCopy( sendBufCol, ACol, height );
        }
    }

    // Communicate
    mpi::AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, g.VCComm() );

    // Unpack
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        if( g.DiagPath( k ) == ownerPath )
        {
            const T* data = &recvBuf[k*portionSize];

            const Int thisPathRank = g.DiagPathRank( k );
            const Int thisRowShift = 
                RawShift( thisPathRank, ownerPathRank, lcm );
            const Int thisLocalWidth = 
                RawLocalLength( width, thisRowShift, lcm );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<thisLocalWidth; ++jLocal )
            {
                const T* dataCol = &data[jLocal*height];
                T* thisCol = 
                    &thisLocalBuffer[(thisRowShift+jLocal*lcm)*thisLDim];
                MemCopy( thisCol, dataCol, height );
            }
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,MR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MR,MC]");
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

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();

    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalHeight = MaxLocalLength(height,c);
    const Int maxLocalWidth = MaxLocalLength(width,r);

    const Int portionSize = 
        std::max(maxLocalHeight*maxLocalWidth,mpi::MIN_COLL_MSG);

    this->auxMemory_.Require( (p+1)*portionSize );

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
        T* originalDataCol = &originalData[jLocal*localHeightOfA];
        MemCopy( originalDataCol, ACol, localHeightOfA );
    }

    // Communicate
    mpi::AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, g.VRComm() );

    // Unpack
    const Int colAlignmentOfA = A.ColAlignment();
    const Int rowAlignmentOfA = A.RowAlignment();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int l=0; l<r; ++l )
    {
        const Int rowShift = RawShift( l, rowAlignmentOfA, r );
        const Int localWidth = RawLocalLength( width, rowShift, r );

        for( Int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[(k+l*c)*portionSize];

            const Int colShift = RawShift( k, colAlignmentOfA, c );
            const Int localHeight = RawLocalLength( height, colShift, c );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( Int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                T* destCol = 
                    &thisLocalBuffer[colShift+(rowShift+jLocal*r)*thisLDim];
                const T* sourceCol = &data[jLocal*localHeight];
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    destCol[iLocal*c] = sourceCol[iLocal];
            }
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,MR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MR,* ]");
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

    const Int c = g.Width();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int maxLocalHeight = MaxLocalLength(height,c);

    const Int portionSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

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
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisLocalBuffer[colShift+j*thisLDim];
            const T* sourceCol = &data[j*localHeight];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                destCol[iLocal*c] = sourceCol[iLocal];
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,MC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,MC]");
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

    const Int r = g.Height();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalWidth = MaxLocalLength(width,r);

    const Int portionSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

    this->auxMemory_.Require( (r+1)*portionSize );

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
      gatheredData, portionSize, g.ColComm() );

    // Unpack
    const Int rowAlignmentOfA = A.RowAlignment();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<r; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const Int rowShift = RawShift( k, rowAlignmentOfA, r );
        const Int localWidth = RawLocalLength( width, rowShift, r );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* dataCol = &data[jLocal*height];
            T* thisCol = &thisLocalBuffer[(rowShift+jLocal*r)*thisLDim];
            MemCopy( thisCol, dataCol, height );
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,VC,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [VC,* ]");
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

    const Int p = g.Size();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int maxLocalHeight = MaxLocalLength(height,p);

    const Int portionSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

    this->auxMemory_.Require( (p+1)*portionSize );

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
      gatheredData, portionSize, g.VCComm() );

    // Unpack
    const Int colAlignmentOfA = A.ColAlignment();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const Int colShift = RawShift( k, colAlignmentOfA, p );
        const Int localHeight = RawLocalLength( height, colShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for 
#endif
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisLocalBuffer[colShift+j*thisLDim];
            const T* sourceCol = &data[j*localHeight];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                destCol[iLocal*p] = sourceCol[iLocal];
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,VC,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,* ]");
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

    const Int p = g.Size();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalWidth = MaxLocalLength(width,p);

    const Int portionSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

    this->auxMemory_.Require( (p+1)*portionSize );

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
      gatheredData, portionSize, g.VCComm() );

    // Unpack
    const Int rowAlignmentOfA = A.RowAlignment();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const Int rowShift = RawShift( k, rowAlignmentOfA, p );
        const Int localWidth = RawLocalLength( width, rowShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* dataCol = &data[jLocal*height];
            T* thisCol = &thisLocalBuffer[(rowShift+jLocal*p)*thisLDim];
            MemCopy( thisCol, dataCol, height );
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,VR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [VR,* ]");
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

    const Int p = g.Size();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localHeightOfA = A.LocalHeight();
    const Int maxLocalHeight = MaxLocalLength(height,p);

    const Int portionSize = std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

    this->auxMemory_.Require( (p+1)*portionSize );

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
      gatheredData, portionSize, g.VRComm() );

    // Unpack
    const Int colAlignmentOfA = A.ColAlignment();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const Int colShift = RawShift( k, colAlignmentOfA, p );
        const Int localHeight = RawLocalLength( height, colShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for 
#endif
        for( Int j=0; j<width; ++j )
        {
            T* destCol = &thisLocalBuffer[colShift+j*thisLDim];
            const T* sourceCol = &data[j*localHeight];
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                destCol[iLocal*p] = sourceCol[iLocal];
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,VR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,VR]");
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

    const Int p = g.Size();
    const Int height = this->Height();
    const Int width = this->Width();
    const Int localWidthOfA = A.LocalWidth();
    const Int maxLocalWidth = MaxLocalLength(width,p);

    const Int portionSize = std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

    this->auxMemory_.Require( (p+1)*portionSize );

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
      gatheredData, portionSize, g.VRComm() );

    // Unpack
    const Int rowAlignmentOfA = A.RowAlignment();
    T* thisLocalBuffer = this->LocalBuffer();
    const Int thisLDim = this->LocalLDim();
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( Int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const Int rowShift = RawShift( k, rowAlignmentOfA, p );
        const Int localWidth = RawLocalLength( width, rowShift, p );

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* dataCol = &data[jLocal*height];
            T* thisCol = &thisLocalBuffer[(rowShift+jLocal*p)*thisLDim];
            MemCopy( thisCol, dataCol, height );
        }
    }
    this->auxMemory_.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline const DistMatrix<T,STAR,STAR,Int>&
DistMatrix<T,STAR,STAR,Int>::operator=( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,* ]");
    this->AssertNotLockedView();
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    if( this->Grid() == A.Grid() )
    {
        this->localMatrix_ = A.LockedLocalMatrix();
    }
    else
    {
        // TODO: Remember why I wrote this...
        if( !mpi::CongruentComms( A.Grid().ViewingComm(),
                                  this->Grid().ViewingComm() ) )
            throw std::logic_error
            ("Redistributing between nonmatching grids currently requires"
             " the viewing communicators to match.");

        // Compute and allocate the amount of required memory
        Int requiredMemory = 0;
        if( A.Grid().VCRank() == 0 )
            requiredMemory += A.Height()*A.Width();
        if( this->Participating() )
            requiredMemory += A.Height()*A.Width();
        this->auxMemory_.Require( requiredMemory );
        T* buffer = this->auxMemory_.Buffer();
        Int offset = 0;
        T* sendBuffer = &buffer[offset];
        if( A.Grid().VCRank() == 0 )
            offset += A.Height()*A.Width();
        T* bcastBuffer = &buffer[offset];

        // Send from the root of A to the root of this matrix's grid
        mpi::Request sendRequest;
        if( A.Grid().VCRank() == 0 )
        {
            for( Int j=0; j<A.Width(); ++j ) 
                for( Int i=0; i<A.Height(); ++i )
                    sendBuffer[i+j*A.Height()] = A.GetLocal(i,j);
            const Int recvViewingRank = this->Grid().VCToViewingMap(0);
            mpi::ISend
            ( sendBuffer, A.Height()*A.Width(), recvViewingRank, 0,
              this->Grid().ViewingComm(), sendRequest );
        }

        // Receive on the root of this matrix's grid and then broadcast
        // over this matrix's owning communicator
        if( this->Participating() )
        {
            if( this->Grid().VCRank() == 0 )
            {
                const Int sendViewingRank = A.Grid().VCToViewingMap(0);
                mpi::Recv
                ( bcastBuffer, A.Height()*A.Width(), sendViewingRank, 0,
                  this->Grid().ViewingComm() );
            }

            mpi::Broadcast
            ( bcastBuffer, A.Height()*A.Width(), 0, this->Grid().VCComm() );

            for( Int j=0; j<A.Width(); ++j )
                for( Int i=0; i<A.Height(); ++i )
                    this->SetLocal(i,j,bcastBuffer[i+j*A.Height()]);
        }

        if( A.Grid().VCRank() == 0 )
            mpi::Wait( sendRequest );
        this->auxMemory_.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::SumOverCol()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SumOverCol");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->Grid();
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int localSize = localHeight*localWidth;
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

    // Sum
    mpi::AllReduce( sendBuf, recvBuf, localSize, mpi::SUM, g.ColComm() );

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
DistMatrix<T,STAR,STAR,Int>::SumOverRow()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SumOverRow");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->Grid();
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int localSize = localHeight*localWidth;
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

    // Sum
    mpi::AllReduce( sendBuf, recvBuf, localSize, mpi::SUM, g.RowComm() );

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
DistMatrix<T,STAR,STAR,Int>::SumOverGrid()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SumOverGrid");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->Grid();
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const Int localHeight = this->LocalHeight();
    const Int localWidth = this->LocalWidth();
    const Int localSize = localHeight*localWidth;
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

    // Sum
    mpi::AllReduce( sendBuf, recvBuf, localSize, mpi::SUM, g.VCComm() );

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

//
// Routines which explicitly work in the complex plane
//

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,STAR,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetRealPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    const Grid& g = this->Grid();
    const Int viewingSize = mpi::CommSize( g.ViewingComm() );
    const Int owningSize = mpi::GroupSize( g.OwningGroup() );

    R u;
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        u = this->GetLocalRealPart(i,j);
    }
    else
    {
        // Have the root broadcast its data
        if( g.VCRank() == 0 )
            u = this->GetLocalRealPart(i,j);
        mpi::Broadcast( &u, 1, g.VCToViewingMap(0), g.ViewingComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,STAR,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetImagPart");
    this->AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    const Grid& g = this->Grid();
    const Int viewingSize = mpi::CommSize( g.ViewingComm() );
    const Int owningSize = mpi::GroupSize( g.OwningGroup() );

    R u;
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        u = this->GetLocalImagPart(i,j);
    }
    else
    {
        // Have the root broadcast its data
        if( g.VCRank() == 0 )
            u = this->GetLocalImagPart(i,j);
        mpi::Broadcast( &u, 1, g.VCToViewingMap(0), g.ViewingComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::SetRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetRealPart");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->SetLocalRealPart(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::SetImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");
    if( this->Participating() )
        this->SetLocalImagPart(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::UpdateRealPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::UpdateRealPart");
    this->AssertValidEntry( i, j );
#endif
    if( this->Participating() )
        this->UpdateLocalRealPart(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::UpdateImagPart
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::UpdateImagPart");
    this->AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");
    if( this->Participating() )
        this->UpdateLocalImagPart(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
