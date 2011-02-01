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
#include "elemental/dist_matrix.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::utilities;
using namespace elemental::imports::mpi;

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Print");
#endif
    const Grid& g = this->Grid();
    if( g.VCRank() == 0 && s != "" )
        cout << s << endl;

    const int height     = this->Height();
    const int width      = this->Width();
    const int localWidth = this->LocalWidth();
    const int c          = g.Width();
    const int rowShift   = this->RowShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Only one process row needs to participate
    if( g.MCRank() == 0 )
    {
        vector<T> sendBuf(height*width,0);
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int i=0; i<height; ++i )
            for( int j=0; j<localWidth; ++j )
                sendBuf[i+(rowShift+j*c)*height] = this->GetLocalEntry(i,j);

        // If we are the root, allocate the receive buffer
        vector<T> recvBuf;
        if( g.MRRank() == 0 )
            recvBuf.resize( height*width );

        // Sum the contributions and send to the root
        Reduce
        ( &sendBuf[0], &recvBuf[0], height*width, MPI_SUM, 0, g.MRComm() );

        if( g.MRRank() == 0 )
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
    }
    Barrier( g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::Align
( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Align");
    this->AssertFreeRowAlignment();
#endif
    this->AlignRows( rowAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::AlignRows
( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Width() )
        throw std::runtime_error( "Invalid row alignment for [* ,MR]" );
#endif
    this->_rowAlignment = rowAlignment;
    this->_rowShift = Shift( g.MRRank(), rowAlignment, g.Width() );
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
elemental::DistMatrixBase<T,Star,MR>::AlignWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([MC,MR])");
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
elemental::DistMatrixBase<T,Star,MR>::AlignWith
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([* ,MR])");
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
elemental::DistMatrixBase<T,Star,MR>::AlignWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([MR,MC])");
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
elemental::DistMatrixBase<T,Star,MR>::AlignWith
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([MR,* ])");
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
elemental::DistMatrixBase<T,Star,MR>::AlignWith
( const DistMatrixBase<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([VR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    this->_rowAlignment = A.ColAlignment() % g.Width();
    this->_rowShift = 
        Shift( g.MRRank(), this->RowAlignment(), g.Width() );
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
elemental::DistMatrixBase<T,Star,MR>::AlignWith
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([* ,VR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    this->_rowAlignment = A.RowAlignment() % g.Width();
    this->_rowShift = 
        Shift( g.MRRank(), this->RowAlignment(), g.Width() );
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
elemental::DistMatrixBase<T,Star,MR>::AlignRowsWith
( const DistMatrixBase<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::AlignRowsWith
( const DistMatrixBase<T,Star,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::AlignRowsWith
( const DistMatrixBase<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::AlignRowsWith
( const DistMatrixBase<T,MR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::AlignRowsWith
( const DistMatrixBase<T,VR,Star>& A )
{ AlignWith( A ); } 

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::AlignRowsWith
( const DistMatrixBase<T,Star,VR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::View
( DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::View");
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
elemental::DistMatrixBase<T,Star,MR>::LockedView
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[*, MR]::LockedView");
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
elemental::DistMatrixBase<T,Star,MR>::View
( DistMatrixBase<T,Star,MR>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const Grid& g = this->Grid();
        const int c   = g.Width();
        const int col = g.MRRank();

        this->_rowAlignment = (A.RowAlignment()+j) % c;
        this->_rowShift = Shift( col, this->RowAlignment(), c );

        const int localWidthBefore = LocalLength( j, A.RowShift(), c );
        const int localWidth = LocalLength( width, this->RowShift(), c );

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
elemental::DistMatrixBase<T,Star,MR>::LockedView
( const DistMatrixBase<T,Star,MR>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const Grid& g = this->Grid();
        const int c = g.Width();
        const int col = g.MRRank();

        this->_rowAlignment = (A.RowAlignment()+j) % c;
        this->_rowShift = Shift( col, this->RowAlignment(), c );

        const int localWidthBefore = LocalLength( j, A.RowShift(), c );
        const int localWidth = LocalLength( width, this->RowShift(), c );

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
elemental::DistMatrixBase<T,Star,MR>::View1x2
( DistMatrixBase<T,Star,MR>& AL,
  DistMatrixBase<T,Star,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::View1x2");
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
elemental::DistMatrixBase<T,Star,MR>::LockedView1x2
( const DistMatrixBase<T,Star,MR>& AL,
  const DistMatrixBase<T,Star,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::LockedView1x2");
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
elemental::DistMatrixBase<T,Star,MR>::View2x1
( DistMatrixBase<T,Star,MR>& AT,
  DistMatrixBase<T,Star,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::View2x1");
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
    this->_localMatrix.View2x1
    ( AT.LocalMatrix(),
      AB.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::LockedView2x1
( const DistMatrixBase<T,Star,MR>& AT,
  const DistMatrixBase<T,Star,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::LockedView2x1");
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
elemental::DistMatrixBase<T,Star,MR>::View2x2
( DistMatrixBase<T,Star,MR>& ATL,
  DistMatrixBase<T,Star,MR>& ATR,
  DistMatrixBase<T,Star,MR>& ABL,
  DistMatrixBase<T,Star,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::View2x2");
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
elemental::DistMatrixBase<T,Star,MR>::LockedView2x2
( const DistMatrixBase<T,Star,MR>& ATL,
  const DistMatrixBase<T,Star,MR>& ATR,
  const DistMatrixBase<T,Star,MR>& ABL,
  const DistMatrixBase<T,Star,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::LockedView2x2");
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
elemental::DistMatrixBase<T,Star,MR>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    this->_height = height;
    this->_width = width;
    this->_localMatrix.ResizeTo
    ( height, LocalLength(width,this->RowShift(),this->Grid().Width()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,Star,MR>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // column within each process row
    const Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    T u;
    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = this->GetLocalEntry(i,jLoc);
    }
    Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Set");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        this->SetLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::Update
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Update");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
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
void
elemental::DistMatrixBase<T,Star,MR>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->Grid();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int c = g.Width();
    const int rowShift = this->RowShift();

    if( shape == Lower )
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            int lastZeroRow = ( side==Left ? j-offset-1
                                           : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                int boundary = min( lastZeroRow+1, height );
                T* thisCol = this->LocalBuffer(0,jLoc);
                memset( thisCol, 0, boundary*sizeof(T) );
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
            int j = rowShift + jLoc*c;
            int firstZeroRow = ( side==Left ? max(j-offset+1,0)
                                            : max(j-offset+height-width+1,0) );
            if( firstZeroRow < height )
            {
                T* thisCol = this->LocalBuffer(firstZeroRow,jLoc);
                memset( thisCol, 0, (height-firstZeroRow)*sizeof(T) );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->Grid();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int c = g.Width();
    const int rowShift = this->RowShift();

    if( shape == Upper )
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
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
            int j = rowShift + jLoc*c;
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
elemental::DistMatrixBase<T,Star,MR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int c = this->Grid().Width();
    const int rowShift = this->RowShift();

    this->SetToZero();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*c;
        if( j < height )
            this->SetLocalEntry(j,jLoc,1);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->Grid();
    const int height     = this->Height();
    const int localWidth = this->LocalWidth();
    const int bufSize    = height*localWidth;

    this->_auxMemory.Require( bufSize );

    // Create random matrix on process row 0, then broadcast
    T* buffer = this->_auxMemory.Buffer();
    if( g.MCRank() == 0 )
    {
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                buffer[i+j*height] = Random<T>();
    }
    Broadcast( buffer, bufSize, 0, g.MCComm() );

    // Unpack
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* bufferCol = &buffer[j*height];
        T* thisCol = this->LocalBuffer(0,j);
        memcpy( thisCol, bufferCol, height*sizeof(T) );
    }
    this->_auxMemory.Release();

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> 
void
elemental::DistMatrixBase<T,Star,MR>::SumOverCol()
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SumOverCol");
    this->AssertNotLockedView();
#endif
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int localSize = max( localHeight*localWidth, MinCollectContrib );

    this->_auxMemory.Require( 2*localSize );
    T* buffer = this->_auxMemory.Buffer();
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[localSize];

    // Pack
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* thisCol = this->LockedLocalBuffer(0,j);
        T* sendBufCol = &sendBuf[j*localHeight];
        memcpy( sendBufCol, thisCol, localHeight*sizeof(T) );
    }

    // AllReduce col
    AllReduce
    ( sendBuf, recvBuf, localSize, MPI_SUM, this->Grid().MCComm() );

    // Unpack
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* recvBufCol = &recvBuf[j*localHeight];
        T* thisCol = this->LocalBuffer(0,j);
        memcpy( thisCol, recvBufCol, localHeight*sizeof(T) );
    }
    this->_auxMemory.Release();

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::ConjugateTransposeFrom
( const DistMatrixBase<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR]::ConjugateTransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.ColAlignment() % g.Width();
            this->_rowShift = 
                Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->RowAlignment() == A.ColAlignment() % g.Width() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();

        const int width = this->Width();
        const int height = this->Height();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = max(height*maxLocalHeightOfA,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<height; ++i )
                originalData[i+j*height] = Conj( A.GetLocalEntry(j,i) );

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int colAlignmentOfA = A.ColAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShiftOfA = Shift( col+k*c, colAlignmentOfA, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidth; ++j )      
            {
                const T* dataCol = &data[j*height];
                T* thisCol = this->LocalBuffer(0,rowOffset+j*r);
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,MR].ConjugateTransposeFrom[VR,* ]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rank = g.VRRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowShift = this->RowShift();

        const int sendRank = (rank+p+rowAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-rowAlignment) % p;

        const int width = this->Width();
        const int height = this->Height();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = max(height*maxLocalHeightOfA,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<height; ++i )
                secondBuffer[i+j*height] = Conj( A.GetLocalEntry(j,i) );

        // Perform the SendRecv: puts the new data into the first buffer
        SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, MPI_ANY_TAG, g.VRComm() );

        // Use the SendRecv as input to the AllGather
        AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

        // Unpack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int colShiftOfA = Shift( col+c*k, rowAlignment, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &data[j*height];
                T* thisCol = this->LocalBuffer(0,rowOffset+j*r);
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MR>::TransposeFrom
( const DistMatrixBase<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR]::TransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.ColAlignment() % g.Width();
            this->_rowShift = 
                Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->RowAlignment() == A.ColAlignment() % g.Width() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();

        const int width = this->Width();
        const int height = this->Height();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = max(height*maxLocalHeightOfA,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<height; ++i )
                originalData[i+j*height] = A.GetLocalEntry(j,i);

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int colAlignmentOfA = A.ColAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShiftOfA = Shift( col+k*c, colAlignmentOfA, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &data[j*height];
                T* thisCol = this->LocalBuffer(0,rowOffset+j*r);
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,MR].TransposeFrom[VR,* ]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rank = g.VRRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowShift = this->RowShift();

        const int sendRank = (rank+p+rowAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-rowAlignment) % p;

        const int width = this->Width();
        const int height = this->Height();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(width,p);

        const int portionSize = max(height*maxLocalHeightOfA,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<height; ++i )
                secondBuffer[i+j*height] = A.GetLocalEntry(j,i);

        // Perform the SendRecv: puts the new data into the first buffer
        SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, MPI_ANY_TAG, g.VRComm() );

        // Use the SendRecv as input to the AllGather
        AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

        // Unpack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int colShiftOfA = Shift( col+c*k, rowAlignment, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &data[j*height];
                T* thisCol = this->LocalBuffer(0,rowOffset+j*r);
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.VCRank() == 0 )
    {
        cerr << 
          "The matrix redistribution [* ,MR] <- [MC,MR] potentially causes a "
          "large amount of cache-thrashing. If possible, avoid it by "
          "performing the redistribution with a (conjugate)-transpose:"
          << endl <<
          "  [MR,* ].(Conjugate)TransposeFrom([MC,MR])" << endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_rowShift = 
                Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        if( A.Height() == 1 )
        {
            const int localWidth = this->LocalWidth();

            this->_auxMemory.Require( localWidth );
            T* bcastBuf = this->_auxMemory.Buffer();

            if( g.MCRank() == A.ColAlignment() )
            {
                this->_localMatrix = A.LockedLocalMatrix();

                // Pack
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( int j=0; j<localWidth; ++j )
                    bcastBuf[j] = this->GetLocalEntry(0,j);
            }

            // Communicate
            Broadcast
            ( bcastBuf, localWidth, A.ColAlignment(), g.MCComm() );

            // Unpack
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidth; ++j )
                this->SetLocalEntry(0,j,bcastBuf[j]);

            this->_auxMemory.Release();
        }
        else
        {
            const int r = g.Height();
            const int height = this->Height();
            const int localWidth = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int portionSize = 
                max(maxLocalHeight*localWidth,MinCollectContrib);

            this->_auxMemory.Require( (r+1)*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* originalData = &buffer[0];
            T* gatheredData = &buffer[portionSize];

            // Pack 
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,j);
                T* originalDataCol = &originalData[j*localHeightOfA];
                memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
            }

            // Communicate
            AllGather
            ( originalData, portionSize,
              gatheredData, portionSize, g.MCComm() );

            // Unpack
            const int colAlignmentOfA = A.ColAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                const T* data = &gatheredData[k*portionSize];

                const int colShift = Shift( k, colAlignmentOfA, r );
                const int localHeight = LocalLength( height, colShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int j=0; j<localWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        this->SetLocalEntry
                            (colShift+i*r,j,data[i+j*localHeight]);
            }

            this->_auxMemory.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,MR] <- [MC,MR]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

        if( A.Height() == 1 )
        {
            const int localWidth = this->LocalWidth();
            T* bcastBuf;

            if( g.MCRank() == A.ColAlignment() )
            {
                const int localWidthOfA = A.LocalWidth();

                this->_auxMemory.Require( localWidth+localWidthOfA );
                T* buffer = this->_auxMemory.Buffer();
                T* sendBuf = &buffer[0];
                bcastBuf   = &buffer[localWidthOfA];

                // Pack
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for( int j=0; j<localWidthOfA; ++j )
                    sendBuf[j] = A.GetLocalEntry(0,j);

                // Communicate
                SendRecv
                ( sendBuf,  localWidthOfA, sendCol, 0,
                  bcastBuf, localWidth,    recvCol, MPI_ANY_TAG,
                  g.MRComm() );
            }
            else
            {
                this->_auxMemory.Require( localWidth );
                bcastBuf = this->_auxMemory.Buffer();
            }

            // Communicate
            Broadcast
            ( bcastBuf, localWidth, A.ColAlignment(), g.MCComm() );

            // Unpack
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidth; ++j )
                this->SetLocalEntry(0,j,bcastBuf[j]);

            this->_auxMemory.Release();
        }
        else
        {
            const int height = this->Height();
            const int localWidth  = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int localWidthOfA  = A.LocalWidth();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int portionSize = 
                max(maxLocalHeight*localWidth,MinCollectContrib);

            this->_auxMemory.Require( (r+1)*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[portionSize];

            // Pack
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidthOfA; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,j);
                T* secondBufferCol = &secondBuffer[j*localHeightOfA];
                memcpy( secondBufferCol, ACol, localHeightOfA*sizeof(T) );
            }

            // Perform the SendRecv: puts the new data into the first buffer
            SendRecv
            ( secondBuffer, portionSize, sendCol, 0,
              firstBuffer,  portionSize, recvCol, MPI_ANY_TAG, g.MRComm() );

            // Use the output of the SendRecv as input to the AllGather
            AllGather
            ( firstBuffer,  portionSize,
              secondBuffer, portionSize, g.MCComm() );

            // Unpack the contents of each member of the process col
            const int colAlignmentOfA = A.ColAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                const T* data = &secondBuffer[k*portionSize];

                const int colShift = Shift( k, colAlignmentOfA, r );
                const int localHeight = LocalLength( height, colShift, r );
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int j=0; j<localWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        this->SetLocalEntry
                            (colShift+i*r,j,data[i+j*localHeight]);
            }

            this->_auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(false,true,0,this->RowAlignment(),g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
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
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,MR] <- [* ,MR]." << endl;
#endif
        const int rank = g.MRRank();
        const int c = g.Width();

        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvRank = (rank+c+rowAlignmentOfA-rowAlignment) % c;

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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* sendBufferCol = &sendBuffer[j*height];
            memcpy( sendBufferCol, ACol, height*sizeof(T) );
        }

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.MRComm() );

        // Unpack
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*height];
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, height*sizeof(T) );
        }
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MR] = [MD,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MR] = [* ,MD] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
    ( new DistMatrix<T,Star,VC>(g) );
    *A_Star_VC = A;

    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
    ( new DistMatrix<T,Star,VR>(true,this->RowAlignment(),g) );
    *A_Star_VR = *A_Star_VC;
    delete A_Star_VC.release(); // lowers memory highwater

    *this = *A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
    ( new DistMatrix<T,VR,Star>(g) );
    *A_VR_Star = A;

    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(g) );
    *A_VC_Star = *A_VR_Star;
    delete A_VR_Star.release(); // lowers memory highwater

    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(false,true,0,this->RowAlignment(),g) );
    *A_MC_MR = *A_VC_Star;
    delete A_VC_Star.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
    ( new DistMatrix<T,Star,VC>(g) );
    *A_Star_VC = A;

    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
    ( new DistMatrix<T,Star,VR>(true,this->RowAlignment(),g) );
    *A_Star_VR = *A_Star_VC;
    delete A_Star_VC.release(); // lowers memory highwater

    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR = *A_Star_VR;
    delete A_Star_VR.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    DistMatrix<T,MC,MR> A_MC_MR(false,true,0,this->RowAlignment(),g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    DistMatrix<T,Star,VR> A_Star_VR(true,this->RowAlignment(),g);

    A_Star_VR = A;
    *this = A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(g) );
    *A_VC_Star = A;

    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(false,true,0,this->RowAlignment(),g) );
    *A_MC_MR = *A_VC_Star;
    delete A_VC_Star.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment() % g.Width();
            this->_rowShift = 
                Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() % g.Width() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();

        const int width = this->Width();
        const int height = this->Height();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidthOfA = MaxLocalLength(width,p);

        const int portionSize = max(height*maxLocalWidthOfA,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &originalData[j*height];
            memcpy( originalDataCol, ACol, height*sizeof(T) );
        }

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int rowAlignmentOfA = A.RowAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShiftOfA = Shift( col+k*c, rowAlignmentOfA, p );
            const int rowOffset = (rowShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, rowShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &data[j*height];
                T* thisCol = this->LocalBuffer(0,rowOffset+j*r);
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [* ,MR] <- [* ,VR]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rank = g.VRRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int rowShift = this->RowShift();

        const int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const int width = this->Width();
        const int height = this->Height();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidthOfA = MaxLocalLength(width,p);

        const int portionSize = max(height*maxLocalWidthOfA,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* secondBufferCol = &secondBuffer[j*height];
            memcpy( secondBufferCol, ACol, height*sizeof(T) );
        }

        // Perform the SendRecv: puts the new data into the first buffer
        SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, MPI_ANY_TAG, g.VRComm() );

        // Use the SendRecv as input to the AllGather
        AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

        // Unpack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int rowShiftOfA = Shift( col+c*k, rowAlignment, p );
            const int rowOffset = (rowShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, rowShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &data[j*height];
                T* thisCol = this->LocalBuffer(0,rowOffset+j*r);
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
        }
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MR>&
elemental::DistMatrixBase<T,Star,MR>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int c = this->Grid().Width();
    const int rowShift = this->RowShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* ACol = A.LockedLocalBuffer(0,rowShift+j*c);
        T* thisCol = this->LocalBuffer(0,j);
        memcpy( thisCol, ACol, localHeight*sizeof(T) );
    }

#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class elemental::DistMatrixBase<int,   Star,MR>;
template class elemental::DistMatrixBase<float, Star,MR>;
template class elemental::DistMatrixBase<double,Star,MR>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,Star,MR>;
template class elemental::DistMatrixBase<dcomplex,Star,MR>;
#endif

