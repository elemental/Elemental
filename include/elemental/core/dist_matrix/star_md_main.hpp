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
using namespace std;

template<typename T>
inline void
DistMatrix<T,STAR,MD>::PrintBase
( ostream& os, const string msg ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::PrintBase");
#endif
    if( this->Grid().VCRank() == 0 && msg != "" )
        os << msg << endl;
        
    const int height     = this->Height();
    const int width      = this->Width();
    const int localWidth = this->LocalWidth();
    const int lcm        = this->Grid().LCM();
    const int inDiagonal = this->InDiagonal();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    vector<T> sendBuf(height*width,0);
    if( inDiagonal )
    {
        const int colShift = this->ColShift();
        const T* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int i=0; i<height; ++i )
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                sendBuf[(colShift+i)+jLocal*lcm*height] = 
                    thisLocalBuffer[i+jLocal*thisLDim];
    }

    // If we are the root, allocate a receive buffer
    vector<T> recvBuf;
    if( this->Grid().VCRank() == 0 )
        recvBuf.resize( height*width );

    // Sum the contributions and send to the root
    mpi::Reduce
    ( &sendBuf[0], &recvBuf[0], height*width, mpi::SUM, 0, 
      this->Grid().VCComm() );

    if( this->Grid().VCRank() == 0 )
    {
        // Print the data
        for( int i=0; i<height; ++i )
        {
            for( int j=0; j<width; ++j )
                os << WrapScalar(recvBuf[i+j*height]) << " ";
            os << "\n";
        }
        os << endl;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::Align( int rowAlignment )
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

template<typename T>
inline void
DistMatrix<T,STAR,MD>::AlignRows( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[STAR,MD]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const elemental::Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Size() )
        throw runtime_error( "Invalid row alignment for [STAR,MD]" );
#endif
    this->_rowAlignment = rowAlignment;
    this->_inDiagonal = ( g.DiagPath() == g.DiagPath(rowAlignment) );
    if( this->_inDiagonal )
        this->_rowShift = Shift( g.DiagPathRank(), rowAlignment, g.Size() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::View( DistMatrix<T,STAR,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->_grid = A._grid;
    this->_height = A.Height();
    this->_width  = A.Width();
    this->_rowAlignment = A.RowAlignment();
    this->_inDiagonal   = A.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = A.RowShift();
        this->_localMatrix.View( A.LocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::View
( int height, int width, int rowAlignment,
  T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->_grid = &grid;
    this->_height = height;
    this->_width = width;
    this->_rowAlignment = rowAlignment;
    this->_inDiagonal = 
        grid.InGrid() && grid.DiagPath()==grid.DiagPath(rowAlignment);
    if( this->_inDiagonal )
    {
        this->_rowShift =
            Shift(grid.DiagPathRank(),
                  grid.DiagPathRank(rowAlignment),
                  grid.LCM());
        const int localWidth = LocalLength(width,this->_rowShift,grid.LCM());
        this->_localMatrix.View( height, localWidth, buffer, ldim );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::LockedView( const DistMatrix<T,STAR,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->_grid = A._grid;
    this->_height = A.Height();
    this->_width  = A.Width();
    this->_rowAlignment = A.RowAlignment();
    this->_inDiagonal   = A.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = A.RowShift();
        this->_localMatrix.LockedView( A.LockedLocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::LockedView
( int height, int width, int rowAlignment,
  const T* buffer, int ldim, const elemental::Grid& grid )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
#endif
    this->_grid = &grid;
    this->_height = height;
    this->_width = width;
    this->_rowAlignment = rowAlignment;
    this->_inDiagonal = 
        grid.InGrid() && grid.DiagPath()==grid.DiagPath(rowAlignment);
    if( this->_inDiagonal )
    {
        this->_rowShift =
            Shift(grid.DiagPathRank(),
                  grid.DiagPathRank(rowAlignment),
                  grid.LCM());
        const int localWidth = LocalLength(width,this->_rowShift,grid.LCM());
        this->_localMatrix.LockedView( height, localWidth, buffer, ldim );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::View
( DistMatrix<T,STAR,MD>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_grid = A._grid;
    this->_height = height;
    this->_width  = width;
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int diagPathRank = g.DiagPathRank(); 
        const int alignmentRank = A.RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        this->_rowAlignment = newAlignmentRank;
        this->_inDiagonal = A.InDiagonal();

        if( this->InDiagonal() )
        {
            this->_rowShift = 
                Shift( diagPathRank,
                       g.DiagPathRank(this->RowAlignment()),
                       lcm );
            int localWidthBefore = LocalLength( j, A.RowShift(), lcm );
            int localWidth = LocalLength( width, this->RowShift(), lcm );
        
            this->_localMatrix.View
            ( A.LocalMatrix(), i, localWidthBefore, height, localWidth );
        }

    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::LockedView
( const DistMatrix<T,STAR,MD>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_grid = A._grid;
    this->_height = height;
    this->_width  = width;
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int diagPathRank = g.DiagPathRank();
        const int alignmentRank = A.RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        this->_rowAlignment = newAlignmentRank;
        this->_inDiagonal = A.InDiagonal();

        if( this->InDiagonal() )
        {
            this->_rowShift = 
                Shift( diagPathRank,
                       g.DiagPathRank( this->RowAlignment() ),
                       lcm );
            int localWidthBefore = LocalLength( j, A.RowShift(), lcm);
            int localWidth = LocalLength( width, this->RowShift(), lcm );
        
            this->_localMatrix.LockedView
            ( A.LockedLocalMatrix(), i, localWidthBefore, height, localWidth );
        }
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::View1x2
( DistMatrix<T,STAR,MD>& AL, DistMatrix<T,STAR,MD>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View1x2");    
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->_grid = AL._grid;
    this->_height = AL.Height();
    this->_width  = AL.Width() + AR.Width();
    this->_rowAlignment = AL.RowAlignment();
    this->_inDiagonal = AL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = AL.RowShift();
        this->_localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::LockedView1x2
( const DistMatrix<T,STAR,MD>& AL, const DistMatrix<T,STAR,MD>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView1x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming1x2( AL, AR );
    AL.AssertSameGrid( AR );
#endif
    this->_grid = AL._grid;
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_rowAlignment = AL.RowAlignment();
    this->_inDiagonal = AL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = AL.RowShift();
        this->_localMatrix.LockedView1x2
        ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::View2x1
( DistMatrix<T,STAR,MD>& AT,
  DistMatrix<T,STAR,MD>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View2x1");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->_grid = AT._grid;
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_rowAlignment = AT.RowAlignment();
    this->_inDiagonal = AT.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = AT.RowShift();
        this->_localMatrix.View2x1
        ( AT.LocalMatrix(), AB.LocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::LockedView2x1
( const DistMatrix<T,STAR,MD>& AT,
  const DistMatrix<T,STAR,MD>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView2x1");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x1( AT, AB );
    AT.AssertSameGrid( AB );
#endif
    this->_grid = AT._grid;
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_rowAlignment = AT.RowAlignment();
    this->_inDiagonal = AT.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = AT.RowShift();
        this->_localMatrix.LockedView2x1
        ( AT.LockedLocalMatrix(), AB.LockedLocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::View2x2
( DistMatrix<T,STAR,MD>& ATL, DistMatrix<T,STAR,MD>& ATR,
  DistMatrix<T,STAR,MD>& ABL, DistMatrix<T,STAR,MD>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View2x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->_grid = ATL._grid;
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_rowAlignment = ATL.RowAlignment();
    this->_inDiagonal = ATL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = ATL.RowShift();
        this->_localMatrix.View2x2
        ( ATL.LocalMatrix(), ATR.LocalMatrix(),
          ABL.LocalMatrix(), ABR.LocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::LockedView2x2
( const DistMatrix<T,STAR,MD>& ATL, const DistMatrix<T,STAR,MD>& ATR,
  const DistMatrix<T,STAR,MD>& ABL, const DistMatrix<T,STAR,MD>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView2x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
    ATL.AssertSameGrid( ATR );
    ATL.AssertSameGrid( ABL );
    ATL.AssertSameGrid( ABR );
#endif
    this->_grid = ATL._grid;
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_rowAlignment = ATL.RowAlignment();
    this->_inDiagonal = ATL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = ATL.RowShift();
        this->_localMatrix.LockedView2x2
        ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
          ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    this->_height = height;
    this->_width = width;
    if( this->InDiagonal() )
    {
        const int lcm = this->Grid().LCM();
        this->_localMatrix.ResizeTo
        ( height, LocalLength(width,this->RowShift(),lcm) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
DistMatrix<T,STAR,MD>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const elemental::Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    T u;
    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        u = this->GetLocalEntry(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::Set( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::Set");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const elemental::Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        this->SetLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::Update( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::Update");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const elemental::Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        this->UpdateLocalEntry(i,jLoc,u);
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
DistMatrix<T,STAR,MD>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int lcm = this->Grid().LCM();
        const int rowShift = this->RowShift();

        if( shape == LOWER )
        {
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                int j = rowShift + jLocal*lcm;
                int lastZeroRow = ( side==LEFT ? j-offset-1
                                               : j-offset+height-width-1 );
                if( lastZeroRow >= 0 )
                {
                    int boundary = min( lastZeroRow+1, height );
                    T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                    memset( thisCol, 0, boundary*sizeof(T) );
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
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                int j = rowShift + jLocal*lcm;
                int firstZeroRow = 
                    ( side==LEFT ? max(j-offset+1,0)
                                 : max(j-offset+height-width+1,0) );
                if( firstZeroRow < height )
                {
                    T* thisCol = &thisLocalBuffer[firstZeroRow+jLocal*thisLDim];
                    memset( thisCol, 0, (height-firstZeroRow)*sizeof(T) );
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
DistMatrix<T,STAR,MD>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int lcm = this->Grid().LCM();
        const int rowShift = this->RowShift();

        if( shape == UPPER )
        {
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                int j = rowShift + jLocal*lcm;
                int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
                int boundary = min( lastRow+1, height );
                T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
                for( int i=0; i<boundary; ++i )
                    thisCol[i] *= alpha;
            }
        }
        else
        {
            T* thisLocalBuffer = this->LocalBuffer();
            const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                int j = rowShift + jLocal*lcm;
                int firstRow = ( side==LEFT ? max(j-offset,0)
                                            : max(j-offset+height-width,0) );
                T* thisCol = &thisLocalBuffer[firstRow+jLocal*thisLDim];
                for( int i=0; i<(height-firstRow); ++i )
                    thisCol[i] *= alpha;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR, MD>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int lcm = this->Grid().LCM();
        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        const int rowShift = this->RowShift();

        this->_localMatrix.SetToZero();

        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const int j = rowShift + jLocal*lcm;
            if( j < height )
                thisLocalBuffer[j+jLocal*thisLDim] = 1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistMatrix<T,STAR,MD>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToRandom");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            for( int i=0; i<height; ++i )
                this->SetLocalEntry(i,jLocal,SampleUnitBall<T>());
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [MC,MR] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,MC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [MC,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,STAR,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [* ,MR] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [MD,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,STAR,MD>& A )
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
            this->_rowAlignment = A.RowAlignment();
            this->_inDiagonal = A.InDiagonal();
            if( this->InDiagonal() )
                this->_rowShift = A.RowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        this->_localMatrix = A.LockedLocalMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( this->Grid().VCRank() == 0 )
            cerr << "Unaligned [* ,MD] <- [* ,MD]." << endl;
#endif
        throw logic_error( "Unaligned [* ,MD] = [* ,MD] not yet implemented." );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [MR,MC] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [MR,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [* ,MC] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,VC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [VC,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,STAR,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [* ,VC] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,VR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [VR,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,STAR,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [* ,VR] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const DistMatrix<T,STAR,MD>&
DistMatrix<T,STAR,MD>::operator=( const DistMatrix<T,STAR,STAR>& A )
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
        const int lcm = this->Grid().LCM();
        const int rowShift = this->RowShift();

        const int height = this->Height();
        const int localWidth = this->LocalWidth();

        T* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
        const T* ALocalBuffer = A.LockedLocalBuffer();
        const int ALDim = A.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* ACol = &ALocalBuffer[(rowShift+jLocal*lcm)*ALDim];
            T* thisCol = &thisLocalBuffer[jLocal*thisLDim];
            memcpy( thisCol, ACol, height*sizeof(T) );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

} // namespace elemental
