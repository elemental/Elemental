/*
   Copyright (c) 2009-2010, Jack Poulson
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
#include "elemental/dist_matrix.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::utilities;
using namespace elemental::wrappers::mpi;

//----------------------------------------------------------------------------//
// DistMatrixBase                                                             //
//----------------------------------------------------------------------------//

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::Print");
#endif
    const Grid& g = this->GetGrid();
    if( g.VRRank() == 0 && s != "" )
        cout << s << endl;

    const int height     = this->Height();
    const int width      = this->Width();
    const int localWidth = this->LocalWidth();
    const int p          = g.Size();
    const int rowShift   = this->RowShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    vector<T> sendBuf(height*width,0);
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int i=0; i<height; ++i )
        for( int j=0; j<localWidth; ++j )
            sendBuf[i+(rowShift+j*p)*height] = this->GetLocalEntry(i,j);

    // If we are the root, allocate a receive buffer
    vector<T> recvBuf;
    if( g.VRRank() == 0 )
        recvBuf.resize(height*width);

    // Sum the contributions and send to the root
    Reduce( &sendBuf[0], &recvBuf[0], height*width, MPI_SUM, 0, g.VRComm() );

    if( g.VRRank() == 0 )
    {
        // Print the data
        for( int i=0; i<height; ++i )
        {
            for( int j=0; j<width; ++j )
                cout << recvBuf[i+j*height] << " ";
            cout << "\n";
        }
        cout << endl;
    }
    Barrier( g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::Align
( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::Align");
    this->AssertFreeRowAlignment();
#endif
    this->AlignRows( rowAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignRows
( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const Grid& g = this->GetGrid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Size() )
        throw std::runtime_error( "Invalid row alignment for [* ,VR]" );
#endif
    this->_rowAlignment = rowAlignment;
    this->_rowShift = Shift( g.VRRank(), rowAlignment, g.Size() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = Shift( g.VRRank(), this->RowAlignment(), g.Size() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = Shift( g.VRRank(), this->RowAlignment(), g.Size() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignWith
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([MR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = Shift( g.VRRank(), this->RowAlignment(), g.Size() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignWith
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([* ,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = Shift( g.VRRank(), this->RowAlignment(), g.Size() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignWith
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([* ,VR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = A.RowShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignWith
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::AlignWith([VR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = A.ColShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignRowsWith
( const DistMatrixBase<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignRowsWith
( const DistMatrixBase<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignRowsWith
( const DistMatrixBase<T,MR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignRowsWith
( const DistMatrixBase<T,Star,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignRowsWith
( const DistMatrixBase<T,Star,VR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::AlignRowsWith
( const DistMatrixBase<T,VR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::View
( DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = A.RowShift();
    this->_localMatrix.View( A.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::LockedView
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView(A)");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = A.RowShift();
    this->_localMatrix.LockedView( A.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::View
( DistMatrixBase<T,Star,VR>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const Grid& g = this->GetGrid();
        const int rowMajorRank = g.VRRank();
        const int size = g.Size();

        this->_rowAlignment = (A.RowAlignment()+j) % size;
        this->_rowShift = Shift( rowMajorRank, this->RowAlignment(), size );

        const int localWidthBefore = LocalLength( j, A.RowShift(), size );
        const int localWidth = LocalLength( width, this->RowShift(), size );

        this->_localMatrix.View
        ( A.LocalMatrix(), i, localWidthBefore, height, localWidth );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::LockedView
( const DistMatrixBase<T,Star,VR>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const Grid& g = this->GetGrid();
        const int rowMajorRank = g.VRRank();
        const int size = g.Size();

        this->_rowAlignment = (A.RowAlignment()+j) % size;
        this->_rowShift = Shift( rowMajorRank, this->RowAlignment(), size );

        const int localWidthBefore = LocalLength( j, A.RowShift(), size );
        const int localWidth = LocalLength( width, this->RowShift(), size );

        this->_localMatrix.LockedView
        ( A.LockedLocalMatrix(), i, localWidthBefore, height, localWidth );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::View1x2
( DistMatrixBase<T,Star,VR>& AL,
  DistMatrixBase<T,Star,VR>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View1x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_rowAlignment = AL.RowAlignment();
    this->_rowShift = AL.RowShift();
    this->_localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::LockedView1x2
( const DistMatrixBase<T,Star,VR>& AL,
  const DistMatrixBase<T,Star,VR>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView1x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_rowAlignment = AL.RowAlignment();
    this->_rowShift = AL.RowShift();
    this->_localMatrix.LockedView1x2
    ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::View2x1
( DistMatrixBase<T,Star,VR>& AT,
  DistMatrixBase<T,Star,VR>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View2x1");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_rowAlignment = AT.RowAlignment();
    this->_rowShift = AT.RowShift();
    this->_localMatrix.View2x1( AT.LocalMatrix(), AB.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::LockedView2x1
( const DistMatrixBase<T,Star,VR>& AT,
  const DistMatrixBase<T,Star,VR>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView2x1");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_rowAlignment = AT.RowAlignment();
    this->_rowShift = AT.RowShift();
    this->_localMatrix.LockedView2x1
    ( AT.LockedLocalMatrix(), 
      AB.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::View2x2
( DistMatrixBase<T,Star,VR>& ATL,
  DistMatrixBase<T,Star,VR>& ATR,
  DistMatrixBase<T,Star,VR>& ABL,
  DistMatrixBase<T,Star,VR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::View2x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_rowAlignment = ATL.RowAlignment();
    this->_rowShift = ATL.RowShift();
    this->_localMatrix.View2x2
    ( ATL.LocalMatrix(), ATR.LocalMatrix(),
      ABL.LocalMatrix(), ABR.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::LockedView2x2
( const DistMatrixBase<T,Star,VR>& ATL,
  const DistMatrixBase<T,Star,VR>& ATR,
  const DistMatrixBase<T,Star,VR>& ABL,
  const DistMatrixBase<T,Star,VR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::LockedView2x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_rowAlignment = ATL.RowAlignment();
    this->_rowShift = ATL.RowShift();
    this->_localMatrix.LockedView2x2
    ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
      ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    const Grid& g = this->GetGrid();
    this->_height = height;
    this->_width = width;
    this->_localMatrix.ResizeTo
    ( height, LocalLength(width,this->RowShift(),g.Size()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,Star,VR>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const Grid& g = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    T u;
    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.Size();
        u = this->GetLocalEntry(i,jLoc);
    }
    Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::Set");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.Size();
        this->SetLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., SetToIdentity and MakeTrapezoidal
//

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->GetGrid();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int p = g.Size();
    const int rowShift = this->RowShift();

    if( shape == Lower )
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*p;
            int lastZeroRow = ( side==Left ? j-offset-1
                                           : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                int boundary = min( lastZeroRow+1, height );
#ifdef RELEASE
                T* thisCol = this->LocalBuffer(0,jLoc);
                memset( thisCol, 0, boundary*sizeof(T) );
#else
                for( int i=0; i<boundary; ++i )
                    this->SetLocalEntry(i,jLoc,0);
#endif
            }
        }
    }
    else
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*p;
            int firstZeroRow = ( side==Left ? max(j-offset+1,0)
                                            : max(j-offset+height-width+1,0) );
#ifdef RELEASE
            if( firstZeroRow < height )
            {
                T* thisCol = this->LocalBuffer(firstZeroRow,jLoc);
                memset( thisCol, 0, (height-firstZeroRow)*sizeof(T) );
            }
#else
            for( int i=firstZeroRow; i<height; ++i )
                this->SetLocalEntry(i,jLoc,0);
#endif
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->GetGrid();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int p = g.Size();
    const int rowShift = this->RowShift();

    if( shape == Upper )
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*p;
            int lastRow = ( side==Left ? j-offset : j-offset+height-width );
            int boundary = min( lastRow+1, height );
#ifdef RELEASE
            T* thisCol = this->LocalBuffer(0,jLoc);
            for( int i=0; i<boundary; ++i )
                thisCol[i] *= alpha;
#else
            for( int i=0; i<boundary; ++i )
            {
                const T value = this->GetLocalEntry(i,jLoc);
                this->SetLocalEntry(i,jLoc,alpha*value);
            }
#endif
        }
    }
    else
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*p;
            int firstRow = ( side==Left ? max(j-offset,0)
                                        : max(j-offset+height-width,0) );
#ifdef RELEASE
            T* thisCol = this->LocalBuffer(firstRow,jLoc);
            for( int i=0; i<(height-firstRow); ++i )
                thisCol[i] *= alpha;
#else
            for( int i=firstRow; i<height; ++i )
            {
                const T value = this->GetLocalEntry(i,jLoc);
                this->SetLocalEntry(i,jLoc,alpha*value);
            }
#endif
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->GetGrid();
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int p = g.Size();
    const int rowShift = this->RowShift();

    this->SetToZero();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*p;
        if( j < height )
            this->SetLocalEntry(j,jLoc,1);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<height; ++i )
            this->SetLocalEntry(i,j,Random<T>());
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::ConjugateTransposeFrom
( const DistMatrixBase<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[*, VR]::ConjugateTransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.ColAlignment();
            this->_rowShift = 
                Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->RowAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int rowShift = this->RowShift();
        const int colShiftOfA = A.ColShift();
        const int rowOffset = (rowShift-colShiftOfA) / c;

        const int height = this->Height();
        const int localWidth = this->LocalWidth();

#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->SetLocalEntry
                    (i,j,Conj( A.GetLocalEntry(rowOffset+j*r,i)));
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,VR]::ConjugateTransposeFrom" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int col = g.MRRank();
        const int colShiftOfA = A.ColShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const int sendCol = (col+c+(rowAlignment%c)-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-(rowAlignment%c)) % c;
        const int sendRank = sendCol + c*row;

        const int sendRowShift = Shift( sendRank, rowAlignment, p );
        const int sendRowOffset = (sendRowShift-colShiftOfA) / c;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localWidthOfSend = LocalLength(width,sendRowShift,p);

        const int sendSize = height * localWidthOfSend;
        const int recvSize = height * localWidth;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidthOfSend; ++j )
            for( int i=0; i<height; ++i )
                sendBuffer[i+j*height] = 
                    Conj( A.GetLocalEntry(sendRowOffset+j*r,i) );

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, MPI_ANY_TAG, g.MRComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*height]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*height]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::TransposeFrom
( const DistMatrixBase<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR]::TransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.ColAlignment();
            this->_rowShift = 
                Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->RowAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int rowShift = this->RowShift();
        const int colShiftOfA = A.ColShift();
        const int rowOffset = (rowShift-colShiftOfA) / c;

        const int height = this->Height();
        const int localWidth = this->LocalWidth();

#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->SetLocalEntry(i,j,A.GetLocalEntry(rowOffset+j*r,i));
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,VR]::TransposeFrom" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int col = g.MRRank();
        const int colShiftOfA = A.ColShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const int sendCol = (col+c+(rowAlignment%c)-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-(rowAlignment%c)) % c;
        const int sendRank = sendCol + c*row;

        const int sendRowShift = Shift( sendRank, rowAlignment, p );
        const int sendRowOffset = (sendRowShift-colShiftOfA) / c;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localWidthOfSend = LocalLength(width,sendRowShift,p);

        const int sendSize = height * localWidthOfSend;
        const int recvSize = height * localWidth;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidthOfSend; ++j )
            for( int i=0; i<height; ++i )
                sendBuffer[i+j*height] = A.GetLocalEntry(sendRowOffset+j*r,i);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, MPI_ANY_TAG, g.MRComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*height]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*height]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_rowShift = 
                Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,r);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        this->_auxMemory.Require( 2*r*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[r*portionSize];

        // Pack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const int thisRank = col+k*c;
            const int thisRowShift = Shift(thisRank,rowAlignment,p);
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,thisRowOffset+j*r);
                T* dataCol = &(data[j*localHeightOfA]);
                memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = 
                          A.GetLocalEntry(i,thisRowOffset+j*r);
#endif
        }

        // Communicate
        AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, g.MCComm() );

        // Unpack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignmentOfA,r);
            const int thisLocalHeight = LocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    this->SetLocalEntry
                        (thisColShift+i*r,j,data[i+j*thisLocalHeight]);
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,VR] <- [MC,MR]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendCol = (col+c+(rowAlignment%c)-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-(rowAlignment%c)) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,r);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        this->_auxMemory.Require( 2*r*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[r*portionSize];

        // Pack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const int thisRank = sendCol+k*c;
            const int thisRowShift = Shift(thisRank,rowAlignment,p);
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,thisRowOffset+j*r);
                T* dataCol = &(data[j*localHeightOfA]);
                memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = 
                          A.GetLocalEntry(i,thisRowOffset+j*r);
#endif
        }

        // AllToAll to gather all of the unaligned [*,VR] data into firstBuffer
        AllToAll
        ( secondBuffer, portionSize,
          firstBuffer,  portionSize, g.MCComm() );

        // SendRecv: properly align the [*,VR] via a trade in the column
        SendRecv
        ( firstBuffer,  portionSize, sendCol, 0,
          secondBuffer, portionSize, recvCol, MPI_ANY_TAG, g.MRComm() );

        // Unpack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignmentOfA,r);
            const int thisLocalHeight = LocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    this->SetLocalEntry
                        (thisColShift+i*r,j,data[i+j*thisLocalHeight]);
        }

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,MC,MR> A_MC_MR(g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_rowShift = 
                Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int rowShift = this->RowShift();
        const int rowShiftOfA = A.RowShift();
        const int rowOffset = (rowShift-rowShiftOfA) / c;

        const int height = this->Height();
        const int localWidth = this->LocalWidth();

#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,rowOffset+j*r);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, ACol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->SetLocalEntry(i,j,A.GetLocalEntry(i,rowOffset+j*r));
#endif
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,VR] <- [* ,MR]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int col = g.MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const int sendCol = (col+c+(rowAlignment%c)-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-(rowAlignment%c)) % c;
        const int sendRank = sendCol + c*row;

        const int sendRowShift = Shift( sendRank, rowAlignment, p );
        const int sendRowOffset = (sendRowShift-rowShiftOfA) / c;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localWidthOfSend = LocalLength(width,sendRowShift,p);

        const int sendSize = height * localWidthOfSend;
        const int recvSize = height * localWidth;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidthOfSend; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,sendRowOffset+j*r);
            T* sendBufferCol = &(sendBuffer[j*height]);
            memcpy( sendBufferCol, ACol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidthOfSend; ++j )
            for( int i=0; i<height; ++i )
                sendBuffer[i+j*height] = A.GetLocalEntry(i,sendRowOffset+j*r);
#endif

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, MPI_ANY_TAG, g.MRComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*height]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*height]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,VR] = [MD,* ] is not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,VR] = [* ,MD] is not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,Star,VC> A_Star_VC(g);

    A_Star_VC = A;
    *this = A_Star_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
    ( new DistMatrix<T,MR,MC>(g) );
    *A_MR_MC = A;

    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
    ( new DistMatrix<T,Star,VC>(g) );
    *A_Star_VC = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_Star_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,Star,VC> A_Star_VC(g);

    A_Star_VC = A;
    *this = A_Star_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,MC,MR> A_MC_MR(g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int localWidthOfA = A.LocalWidth();

    const int sendSize = height * localWidthOfA;
    const int recvSize = height * localWidth;

    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();
    const int rankCM = g.VCRank();
    const int rankRM = g.VRRank();

    const int rowShift = this->RowShift();
    const int rowShiftOfA = A.RowShift();

    // Compute which rowmajor rank has the rowShift equal to our rowShiftOfA
    const int sendRankRM = (rankRM+(p+rowShiftOfA-rowShift)) % p;

    // Compute which rowmajor rank has the A rowShift that we need
    const int recvRankCM = (rankCM+(p+rowShift-rowShiftOfA)) % p;
    const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

    this->_auxMemory.Require( sendSize + recvSize );

    T* buffer = this->_auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack
#ifdef RELEASE
# ifdef _OPENMP
    #pragma omp parallel for
# endif
    for( int j=0; j<localWidthOfA; ++j )
    {
        const T* ACol = A.LockedLocalBuffer(0,j);
        T* sendBufferCol = &(sendBuffer[j*height]);
        memcpy( sendBufferCol, ACol, height*sizeof(T) );
    }
#else
# ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
# endif
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<height; ++i )
            sendBuffer[i+j*height] = A.GetLocalEntry(i,j);
#endif

    // Communicate
    SendRecv
    ( sendBuffer, sendSize, sendRankRM, 0,
      recvBuffer, recvSize, recvRankRM, MPI_ANY_TAG, g.VRComm() );

    // Unpack
#ifdef RELEASE
# ifdef _OPENMP
    #pragma omp parallel for
# endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* recvBufferCol = &(recvBuffer[j*height]);
        T* thisCol = this->LocalBuffer(0,j);
        memcpy( thisCol, recvBufferCol, height*sizeof(T) );
    }
#else
# ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
# endif
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<height; ++i )
            this->SetLocalEntry(i,j,recvBuffer[i+j*height]);
#endif

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
    ( new DistMatrix<T,MR,MC>(g) );
    *A_MR_MC = A;

    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
    ( new DistMatrix<T,Star,VC>(g) );
    *A_Star_VC = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_Star_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,VR]");
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
        const Grid& g = this->GetGrid();
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,VR] <- [* ,VR]." << endl;
#endif
        const int rank = g.VRRank();
        const int p = g.Size();

        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = height * localWidthOfA;
        const int recvSize = height * localWidth;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* sendBufferCol = &(sendBuffer[j*height]);
            memcpy( sendBufferCol, ACol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<height; ++i )
                sendBuffer[i+j*height] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.VRComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*height]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*height]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VR>&
elemental::DistMatrixBase<T,Star,VR>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VR] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int p = this->GetGrid().Size();
    const int rowShift = this->RowShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
#ifdef RELEASE
# ifdef _OPENMP
    #pragma omp parallel for
# endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* ACol = A.LockedLocalBuffer(0,rowShift+j*p);
        T* thisCol = this->LocalBuffer(0,j);
        memcpy( thisCol, ACol, localHeight*sizeof(T) );
    }
#else
# ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
# endif
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->SetLocalEntry(i,j,A.GetLocalEntry(i,rowShift+j*p));
#endif

#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::SumScatterFrom
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SumScatterFrom( [* ,MR] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_rowShift = 
                Shift( g.VRRank(), this->RowAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();
        const int rowAlignment = this->RowAlignment();
        const int rowShiftOfA = A.RowShift();

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalWidth = MaxLocalLength( width, p );

        const int recvSize = max(localHeight*maxLocalWidth,MinCollectContrib);
        const int sendSize = r*recvSize;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        vector<int> recvSizes(r);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisRank = col+k*c;
            const int thisRowShift = Shift( thisRank, rowAlignment, p );
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const int thisLocalWidth = LocalLength( width, thisRowShift, p );

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,thisRowOffset+j*r);
                T* dataCol = &(data[j*localHeight]);
                memcpy( dataCol, ACol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    data[i+j*localHeight] = 
                        A.GetLocalEntry(i,thisRowOffset+j*r);
#endif
        }

        // Reduce-scatter over each process column
        ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

        // Unpack our received data
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*localHeight]);
#endif

        this->_auxMemory.Release();
    }
    else
    {
        throw logic_error
              ( "Unaligned [* ,VR]::ReduceScatterFrom( [* ,MR] ) is not "
                "yet implemented." );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VR>::SumScatterUpdate
( T alpha, const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SumScatterUpdate( [* ,MR] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( this->RowAlignment() % g.Width() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();
        const int rowAlignment = this->RowAlignment();
        const int rowShiftOfA = A.RowShift();

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalWidth = MaxLocalLength( width, p );

        const int recvSize = max(localHeight*maxLocalWidth,MinCollectContrib);
        const int sendSize = r*recvSize;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        vector<int> recvSizes(r);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisRank = col+k*c;
            const int thisRowShift = Shift( thisRank, rowAlignment, p );
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const int thisLocalWidth = LocalLength( width, thisRowShift, p );

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,thisRowOffset+j*r);
                T* dataCol = &(data[j*localHeight]);
                memcpy( dataCol, ACol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    data[i+j*localHeight] = 
                        A.GetLocalEntry(i,thisRowOffset+j*r);
#endif
        }

        // Reduce-scatter over each process column
        ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

        // Unpack our received data
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            for( int i=0; i<localHeight; ++i )
                thisCol[i] += alpha*recvBufferCol[i];
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
        {
            for( int i=0; i<localHeight; ++i )
            {
                const T value = this->GetLocalEntry(i,j);
                this->SetLocalEntry
                    (i,j,value+alpha*recvBuffer[i+j*localHeight]);
            }
        }
#endif

        this->_auxMemory.Release();
    }
    else
    {
        throw logic_error
              ( "Unaligned [* ,VR]::ReduceScatterUpdate( [* ,MR] ) is not "
                "yet implemented." );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// DistMatrix                                                                 //
//----------------------------------------------------------------------------//

template<typename R>
void
elemental::DistMatrix<R,Star,VR>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const Grid& g = this->GetGrid();
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int p = g.Size();
    const int rowShift = this->RowShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*p;
        if( j < height )
        {
            const R value = this->GetLocalEntry(j,jLoc);
            this->SetLocalEntry(j,jLoc,value+this->Width());
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,Star,VR>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const Grid& g = this->GetGrid();
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int p = g.Size();
    const int rowShift = this->RowShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*p;
        if( j < height )
        {
            const R value = real(this->GetLocalEntry(j,jLoc));
            this->SetLocalEntry(j,jLoc,value+this->Width());
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,VR>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const Grid& g = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    R u;
    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.Size();
        u = real(this->GetLocalEntry(i,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,VR>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const Grid& g = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    R u;
    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.Size();
        u = imag(this->GetLocalEntry(i,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,VR>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.Size();
        const R v = imag(this->GetLocalEntry(i,jLoc));
        this->SetLocalEntry(i,jLoc,complex<R>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,VR>::SetImag
( int i, int j, R v )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.Size();
        const R u = real(this->GetLocalEntry(i,jLoc));
        this->SetLocalEntry(i,jLoc,complex<R>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrixBase<int,   Star,VR>;
template class elemental::DistMatrixBase<float, Star,VR>;
template class elemental::DistMatrixBase<double,Star,VR>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,Star,VR>;
template class elemental::DistMatrixBase<dcomplex,Star,VR>;
#endif

template class elemental::DistMatrix<int,   Star,VR>;
template class elemental::DistMatrix<float, Star,VR>;
template class elemental::DistMatrix<double,Star,VR>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,Star,VR>;
template class elemental::DistMatrix<dcomplex,Star,VR>;
#endif

