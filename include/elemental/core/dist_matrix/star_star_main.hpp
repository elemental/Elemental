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

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( g.InGrid() )
    {
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
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::View( DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View");
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
        this->localMatrix_.View( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::View
( Int height, Int width, 
  T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
        this->localMatrix_.View( height, width, buffer, ldim );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::LockedView
( const DistMatrix<T,STAR,STAR,Int>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView");
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = A.Height();
    this->width_ = A.Width();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
        this->localMatrix_.LockedView( A.LockedLocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::LockedView
( Int height, Int width, 
  const T* buffer, Int ldim, const elem::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView");
#endif
    this->Empty();

    this->grid_ = &grid;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
        this->localMatrix_.LockedView( height, width, buffer, ldim );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::View
( DistMatrix<T,STAR,STAR,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View");
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
        this->localMatrix_.View( A.LocalMatrix(), i, j, height, width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::LockedView
( const DistMatrix<T,STAR,STAR,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView");
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->Empty();

    this->grid_ = A.grid_;
    this->height_ = height;
    this->width_ = width;
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.LockedView
        ( A.LockedLocalMatrix(), i, j, height, width );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::View1x2
( DistMatrix<T,STAR,STAR,Int>& AL, DistMatrix<T,STAR,STAR,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View1x2");
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->Empty();

    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
        this->localMatrix_.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::LockedView1x2
( const DistMatrix<T,STAR,STAR,Int>& AL, const DistMatrix<T,STAR,STAR,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView1x2");
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->Empty();

    this->grid_ = AL.grid_;
    this->height_ = AL.Height();
    this->width_ = AL.Width() + AR.Width();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.LockedView1x2
        ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::View2x1
( DistMatrix<T,STAR,STAR,Int>& AT,
  DistMatrix<T,STAR,STAR,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View2x1");
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->Empty();

    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.View2x1
        ( AT.LocalMatrix(),
          AB.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::LockedView2x1
( const DistMatrix<T,STAR,STAR,Int>& AT,
  const DistMatrix<T,STAR,STAR,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView2x1");
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->Empty();

    this->grid_ = AT.grid_;
    this->height_ = AT.Height() + AB.Height();
    this->width_ = AT.Width();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.LockedView2x1
        ( AT.LockedLocalMatrix(),
          AB.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::View2x2
( DistMatrix<T,STAR,STAR,Int>& ATL, DistMatrix<T,STAR,STAR,Int>& ATR,
  DistMatrix<T,STAR,STAR,Int>& ABL, DistMatrix<T,STAR,STAR,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View2x2");
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->Empty();

    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->viewing_ = true;
    this->lockedView_ = false;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.View2x2
        ( ATL.LocalMatrix(), ATR.LocalMatrix(),
          ABL.LocalMatrix(), ABR.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::LockedView2x2
( const DistMatrix<T,STAR,STAR,Int>& ATL, 
  const DistMatrix<T,STAR,STAR,Int>& ATR,
  const DistMatrix<T,STAR,STAR,Int>& ABL, 
  const DistMatrix<T,STAR,STAR,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView2x2");
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->Empty();

    this->grid_ = ATL.grid_;
    this->height_ = ATL.Height() + ABL.Height();
    this->width_ = ATL.Width() + ATR.Width();
    this->viewing_ = true;
    this->lockedView_ = true;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.LockedView2x2
        ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
          ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    }
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
    if( this->Grid().InGrid() )
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
    if( this->Grid().InGrid() )
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
    if( this->Grid().InGrid() )
    {
        this->UpdateLocal(i,j,u);
    }
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
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
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
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

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisLocalBuffer[(colShift+iLocal*r)+
                                        (rowShift+jLocal*c)*thisLDim] = 
                            data[iLocal+jLocal*localHeight];
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const Int r = g.Height();
        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeight = MaxLocalLength(height,r);

        const Int portionSize = 
            std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int colShift = RawShift( k, colAlignmentOfA, r );
            const Int localHeight = RawLocalLength( height, colShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int j=0; j<width; ++j )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colShift+iLocal*r)+j*thisLDim] =
                        data[iLocal+j*localHeight];
        }
        this->auxMemory_.Release();
    }
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const Int c = g.Width();
        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidth = MaxLocalLength(width,c);

        const Int portionSize = 
            std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int rowShift = RawShift( k, rowAlignmentOfA, c );
            const Int localWidth = RawLocalLength( width, rowShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
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
    }
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int lcm = g.LCM();
        const Int ownerPath = g.DiagPath( A.ColAlignment() );
        const Int ownerPathRank = g.DiagPathRank( A.ColAlignment() );

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeight = A.LocalHeight();
        const Int maxLocalHeight = MaxLocalLength( height, lcm );
        const Int portionSize = 
            std::max( maxLocalHeight*width, mpi::MIN_COLL_MSG );

        // Since a MD communicator has not been implemented, we will take
        // the suboptimal route of 'rounding up' everyone's contribution over 
        // the VC communicator.
        this->auxMemory_.Require( (p+1)*portionSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        if( A.InDiagonal() )
        {
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
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

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int j=0; j<width; ++j )
                    for( Int iLocal=0; iLocal<thisLocalHeight; ++iLocal )
                        thisLocalBuffer[(thisColShift+iLocal*lcm)+j*thisLDim] =
                            data[iLocal+j*thisLocalHeight];
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int lcm = g.LCM();
        const Int ownerPath = g.DiagPath( A.RowAlignment() );
        const Int ownerPathRank = g.DiagPathRank( A.RowAlignment() );

        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidth = A.LocalWidth();
        const Int maxLocalWidth = MaxLocalLength( width, lcm );
        const Int portionSize = 
            std::max( height*maxLocalWidth, mpi::MIN_COLL_MSG );

        // Since a MD communicator has not been implemented, we will take
        // the suboptimal route of 'rounding up' everyone's contribution over 
        // the VC communicator.
        this->auxMemory_.Require( (p+1)*portionSize );
        T* buffer = this->auxMemory_.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        if( A.InDiagonal() )
        {
            const T* ALocalBuffer = A.LockedLocalBuffer();
            const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
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

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
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
    }
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
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
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
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

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( Int jLocal=0; jLocal<localWidth; ++jLocal )
                    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                        thisLocalBuffer[(colShift+iLocal*c)+
                                        (rowShift+jLocal*r)*thisLDim] = 
                            data[iLocal+jLocal*localHeight];
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const Int c = g.Width();
        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeight = MaxLocalLength(height,c);

        const Int portionSize = 
            std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (c+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int colShift = RawShift( k, colAlignmentOfA, c );
            const Int localHeight = RawLocalLength( height, colShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int j=0; j<width; ++j )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colShift+iLocal*c)+j*thisLDim] = 
                        data[iLocal+j*localHeight];
        }
        this->auxMemory_.Release();
    }
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const Int r = g.Height();
        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidth = MaxLocalLength(width,r);

        const Int portionSize = 
            std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (r+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int rowShift = RawShift( k, rowAlignmentOfA, r );
            const Int localWidth = RawLocalLength( width, rowShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
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
    }
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeight = MaxLocalLength(height,p);

        const Int portionSize = 
            std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (p+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int colShift = RawShift( k, colAlignmentOfA, p );
            const Int localHeight = RawLocalLength( height, colShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int j=0; j<width; ++j )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colShift+iLocal*p)+j*thisLDim] = 
                        data[iLocal+j*localHeight];
        }
        this->auxMemory_.Release();
    }
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidth = MaxLocalLength(width,p);

        const Int portionSize = 
            std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (p+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int rowShift = RawShift( k, rowAlignmentOfA, p );
            const Int localWidth = RawLocalLength( width, rowShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
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
    }
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int height = this->Height();
        const Int width = this->Width();
        const Int localHeightOfA = A.LocalHeight();
        const Int maxLocalHeight = MaxLocalLength(height,p);

        const Int portionSize = 
            std::max(maxLocalHeight*width,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (p+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int colShift = RawShift( k, colAlignmentOfA, p );
            const Int localHeight = RawLocalLength( height, colShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( Int j=0; j<width; ++j )
                for( Int iLocal=0; iLocal<localHeight; ++iLocal )
                    thisLocalBuffer[(colShift+iLocal*p)+j*thisLDim] = 
                        data[iLocal+j*localHeight];
        }
        this->auxMemory_.Release();
    }
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
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const elem::Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int height = this->Height();
        const Int width = this->Width();
        const Int localWidthOfA = A.LocalWidth();
        const Int maxLocalWidth = MaxLocalLength(width,p);

        const Int portionSize = 
            std::max(height*maxLocalWidth,mpi::MIN_COLL_MSG);

        this->auxMemory_.Require( (p+1)*portionSize );

        T* buffer = this->auxMemory_.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const Int ALDim = A.LocalLDim();
#ifdef _OPENMP
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( Int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const Int rowShift = RawShift( k, rowAlignmentOfA, p );
            const Int localWidth = RawLocalLength( width, rowShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
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
    }
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
        if( !mpi::CongruentComms( A.Grid().ViewingComm(),
                                  this->Grid().ViewingComm() ) )
            throw std::logic_error
                ("Redistributing between nonmatching grids currently requires"
                 " the viewing communicators to match.");

        // Compute and allocate the amount of required memory
        Int requiredMemory = 0;
        if( A.Grid().VCRank() == 0 )
            requiredMemory += A.Height()*A.Width();
        if( this->Grid().InGrid() )
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
        if( this->Grid().InGrid() )
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
    if( this->Grid().InGrid() )
    {
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            T* sendBufCol = &sendBuf[jLocal*localHeight];
            MemCopy( sendBufCol, thisCol, localHeight );
        }

        // Sum
        mpi::AllReduce
        ( sendBuf, recvBuf, localSize, mpi::SUM, this->Grid().ColComm() );

        // Unpack
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufCol = &recvBuf[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( thisCol, recvBufCol, localHeight );
        }
        this->auxMemory_.Release();
    }
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
    if( this->Grid().InGrid() )
    {
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            T* sendBufCol = &sendBuf[jLocal*localHeight];
            MemCopy( sendBufCol, thisCol, localHeight );
        }

        // Sum
        mpi::AllReduce
        ( sendBuf, recvBuf, localSize, mpi::SUM, this->Grid().RowComm() );

        // Unpack
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufCol = &recvBuf[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( thisCol, recvBufCol, localHeight );
        }
        this->auxMemory_.Release();
    }
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
    if( this->Grid().InGrid() )
    {
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            T* sendBufCol = &sendBuf[jLocal*localHeight];
            MemCopy( sendBufCol, thisCol, localHeight );
        }

        // Sum
        mpi::AllReduce
        ( sendBuf, recvBuf, localSize, mpi::SUM, this->Grid().VCComm() );

        // Unpack
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* recvBufCol = &recvBuf[jLocal*localHeight];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            MemCopy( thisCol, recvBufCol, localHeight );
        }
        this->auxMemory_.Release();
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
DistMatrix<T,STAR,STAR,Int>::GetReal( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetReal");
    AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    const Int viewingSize = mpi::CommSize( this->Grid().ViewingComm() );
    const Int owningSize = mpi::GroupSize( this->Grid().OwningGroup() );

    R u;
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        u = this->GetRealLocal(i,j);
    }
    else
    {
        // Have the root broadcast its data
        if( this->Grid().VCRank() == 0 )
            u = this->GetRealLocal(i,j);
        mpi::Broadcast
        ( &u, 1, this->Grid().VCToViewingMap(0), this->Grid().ViewingComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,STAR,Int>::GetImag( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetImag");
    AssertValidEntry( i, j );
#endif
    typedef typename Base<T>::type R;

    const Int viewingSize = mpi::CommSize( this->Grid().ViewingComm() );
    const Int owningSize = mpi::GroupSize( this->Grid().OwningGroup() );

    R u;
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        u = this->GetImagLocal(i,j);
    }
    else
    {
        // Have the root broadcast its data
        if( this->Grid().VCRank() == 0 )
            u = this->GetImagLocal(i,j);
        mpi::Broadcast
        ( &u, 1, this->Grid().VCToViewingMap(0), this->Grid().ViewingComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::SetReal( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetReal");
    AssertValidEntry( i, j );
#endif
    if( this->Grid().InGrid() )
        this->SetRealLocal(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::SetImag( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetImag");
    AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");
    if( this->Grid().InGrid() )
        this->SetImagLocal(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::UpdateReal
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::UpdateReal");
    AssertValidEntry( i, j );
#endif
    if( this->Grid().InGrid() )
        this->UpdateRealLocal(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::UpdateImag
( Int i, Int j, typename Base<T>::type u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::UpdateImag");
    AssertValidEntry( i, j );
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");
    if( this->Grid().InGrid() )
        this->UpdateImagLocal(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
