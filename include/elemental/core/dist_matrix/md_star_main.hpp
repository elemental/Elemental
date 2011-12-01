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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::PrintBase
( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::PrintBase");
#endif
    const elemental::Grid& g = this->Grid();
    if( g.VCRank() == 0 && msg != "" )
        os << msg << std::endl;
        
    const int height      = this->Height();
    const int width       = this->Width();
    const int localHeight = this->LocalHeight();
    const int inDiagonal  = this->InDiagonal();
    const int lcm         = g.LCM();

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
            const int colShift = this->ColShift();
            const T* thisLocalBuffer = this->LockedLocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                for( int j=0; j<width; ++j )
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
            for( int i=0; i<height; ++i )
            {
                for( int j=0; j<width; ++j )
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::Align( int colAlignment )
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::AlignCols( int colAlignment )
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::View( DistMatrix<T,MD,STAR>& A )
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::View
( int height, int width, int colAlignment,
  T* buffer, int ldim, const elemental::Grid& grid )
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
        const int localHeight = LocalLength(height,this->colShift_,grid.LCM());
        this->localMatrix_.View( localHeight, width, buffer, ldim );
    }
    this->viewing_ = true;
    this->lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MD,STAR>::LockedView( const DistMatrix<T,MD,STAR>& A )
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::LockedView
( int height, int width, int colAlignment,
  const T* buffer, int ldim, const elemental::Grid& grid )
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
        const int localHeight = LocalLength(height,this->colShift_,grid.LCM());
        this->localMatrix_.LockedView( localHeight, width, buffer, ldim );
    }
    this->viewing_ = true;
    this->lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MD,STAR>::View
( DistMatrix<T,MD,STAR>& A, int i, int j, int height, int width )
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
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int diagPathRank = g.DiagPathRank();
        const int alignmentRank = A.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        this->colAlignment_ = newAlignmentRank;
        this->inDiagonal_ = A.InDiagonal();

        if( this->inDiagonal_ )
        {
            this->colShift_ = 
                Shift( diagPathRank,
                       g.DiagPathRank( this->ColAlignment() ),
                       lcm );
            int localHeightBefore = LocalLength( i, A.ColShift(), lcm);
            int localHeight = LocalLength( height, this->ColShift(), lcm );

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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::LockedView
( const DistMatrix<T,MD,STAR>& A, int i, int j, int height, int width )
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
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int diagPathRank = g.DiagPathRank();
        const int alignmentRank = A.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        this->colAlignment_ = newAlignmentRank;
        this->inDiagonal_ = A.InDiagonal();

        if( this->InDiagonal() )
        {
            this->colShift_ = 
                Shift( diagPathRank,
                       g.DiagPathRank( this->ColAlignment() ),
                       lcm );
            int localHeightBefore = LocalLength( i, A.ColShift(), lcm);
            int localHeight = LocalLength( height, this->ColShift(), lcm );
        
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::View1x2
( DistMatrix<T,MD,STAR>& AL, DistMatrix<T,MD,STAR>& AR )
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::LockedView1x2
( const DistMatrix<T,MD,STAR>& AL, const DistMatrix<T,MD,STAR>& AR )
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::View2x1
( DistMatrix<T,MD,STAR>& AT,
  DistMatrix<T,MD,STAR>& AB )
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::LockedView2x1
( const DistMatrix<T,MD,STAR>& AT,
  const DistMatrix<T,MD,STAR>& AB )
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::View2x2
( DistMatrix<T,MD,STAR>& ATL, DistMatrix<T,MD,STAR>& ATR,
  DistMatrix<T,MD,STAR>& ABL, DistMatrix<T,MD,STAR>& ABR )
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::LockedView2x2
( const DistMatrix<T,MD,STAR>& ATL, const DistMatrix<T,MD,STAR>& ATR,
  const DistMatrix<T,MD,STAR>& ABL, const DistMatrix<T,MD,STAR>& ABR )
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::ResizeTo( int height, int width )
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
        const int lcm = this->Grid().LCM();
        this->localMatrix_.ResizeTo
        ( LocalLength(height,this->ColShift(),lcm), width );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
DistMatrix<T,MD,STAR>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const elemental::Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    T u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        u = this->GetLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,MD,STAR>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const elemental::Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        this->SetLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MD,STAR>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const elemental::Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        this->UpdateLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., SetToIdentity and MakeTrapezoidal
//

template<typename T>
inline void
DistMatrix<T,MD,STAR>::MakeTrapezoidal
( Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int lcm = this->Grid().LCM();
        const int colShift = this->ColShift();

        if( uplo == LOWER )
        {
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int lastZeroRow = ( side==LEFT ? j-offset-1
                                               : j-offset+height-width-1 );
                if( lastZeroRow >= 0 )
                {
                    int boundary = std::min( lastZeroRow+1, height );
                    int numZeroRows = RawLocalLength( boundary, colShift, lcm );
                    T* thisCol = &thisLocalBuffer[j*thisLDim];
                    std::memset( thisCol, 0, numZeroRows*sizeof(T) );
                }
            }
        }
        else
        {
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int firstZeroRow = 
                    ( side==LEFT ? std::max(j-offset+1,0)
                                 : std::max(j-offset+height-width+1,0) );
                int numNonzeroRows = RawLocalLength(firstZeroRow,colShift,lcm);
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

template<typename T>
inline void
DistMatrix<T,MD,STAR>::ScaleTrapezoid
( T alpha, Side side, UpperOrLower uplo, int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::ScaleTrapezoid");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int lcm = this->Grid().LCM();
        const int colShift = this->ColShift();

        if( uplo == UPPER )
        {
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
                int boundary = std::min( lastRow+1, height );
                int numRows = RawLocalLength( boundary, colShift, lcm );
                T* thisCol = &thisLocalBuffer[j*thisLDim];
                for( int iLocal=0; iLocal<numRows; ++iLocal )
                    thisCol[iLocal] *= alpha;
            }
        }
        else
        {
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int firstRow = 
                    ( side==LEFT ? std::max(j-offset,0)
                                 : std::max(j+height-width-offset,0) );
                int numZeroRows = RawLocalLength( firstRow, colShift, lcm );
                T* thisCol = &thisLocalBuffer[numZeroRows+j*thisLDim];
                for( int iLocal=0; iLocal<(localHeight-numZeroRows); ++iLocal )
                    thisCol[iLocal] *= alpha;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MD,STAR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int lcm = this->Grid().LCM();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int colShift = this->ColShift();

        this->localMatrix_.SetToZero();

        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*lcm;
            if( i < width )
                thisLocalBuffer[iLocal+i*thisLDim] = 1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,MD,STAR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        for( int j=0; j<width; ++j )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                this->SetLocalEntry(iLocal,j,SampleUnitBall<T>());
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,MC,MR>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,MC,STAR>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,MR>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,MD,STAR>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,MD>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,MR,MC>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,MR,STAR>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,MC>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,VC,STAR>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,VC>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,VR,STAR>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,VR>& A )
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

template<typename T>
inline const DistMatrix<T,MD,STAR>&
DistMatrix<T,MD,STAR>::operator=( const DistMatrix<T,STAR,STAR>& A )
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
        const int lcm = this->grid_->LCM();
        const int colShift = this->ColShift();

        const int width = this->Width();
        const int localHeight = this->LocalHeight();

        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<width; ++j )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                thisLocalBuffer[iLocal+j*thisLDim] = 
                    ALocalBuffer[(colShift+iLocal*lcm)+j*ALDim];
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

} // namespace elemental
