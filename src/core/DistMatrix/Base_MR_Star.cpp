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

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::Print");
#endif
    const Grid& g = this->Grid();
    if( g.VCRank() == 0 && s != "" )
        cout << s << endl;

    const int height      = this->Height();
    const int width       = this->Width();
    const int localHeight = this->LocalHeight();
    const int c           = g.Width();
    const int colShift    = this->ColShift();

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
        for( int i=0; i<localHeight; ++i )
            for( int j=0; j<width; ++j )
                sendBuf[colShift+i*c+j*height] = this->GetLocalEntry(i,j);

        // If we are the root, allocate a receive buffer
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
elemental::DistMatrixBase<T,MR,Star>::Align
( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::Align");
    this->AssertFreeColAlignment();
#endif
    this->AlignCols( colAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignCols
( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Width() )
        throw std::runtime_error( "Invalid column alignment for [MR,* ]" );
#endif
    this->_colAlignment = colAlignment;
    this->_colShift = Shift( g.MRRank(), colAlignment, g.Width() );
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::AlignWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_colShift = A.ColShift();
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignWith
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::AlignWith([MR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_colShift = A.ColShift();
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::AlignWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_colShift = A.RowShift();
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignWith
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::AlignWith([* ,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_colShift = A.RowShift();
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignWith
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::AlignWith([VR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignWith
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::AlignWith([* ,VR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    this->_colAlignment = A.RowAlignment();
    this->_colShift = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignColsWith
( const DistMatrixBase<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignColsWith
( const DistMatrixBase<T,MR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignColsWith
( const DistMatrixBase<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignColsWith
( const DistMatrixBase<T,Star,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignColsWith
( const DistMatrixBase<T,VR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::AlignColsWith
( const DistMatrixBase<T,Star,VR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::View
( DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = A.ColShift();
    this->_localMatrix.View( A.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::LockedView
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = A.ColShift();
    this->_localMatrix.LockedView( A.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::View
( DistMatrixBase<T,MR,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::View");
    this->AssertFreeColAlignment();
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

        this->_colAlignment = (A.ColAlignment()+i) % c;
        this->_colShift = Shift( col, this->ColAlignment(), c );

        const int localHeightBefore = LocalLength( i, A.ColShift(), c );
        const int localHeight = LocalLength( height, this->ColShift(), c );

        this->_localMatrix.View
        ( A.LocalMatrix(), localHeightBefore, j, localHeight, width );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::LockedView
( const DistMatrixBase<T,MR,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::LockedView");
    this->AssertFreeColAlignment();
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

        this->_colAlignment = (A.ColAlignment()+i) % c;
        this->_colShift = Shift( col, this->ColAlignment(), c );

        const int localHeightBefore = LocalLength( i, A.ColShift(), c );
        const int localHeight = LocalLength( height, this->ColShift(), c );

        this->_localMatrix.LockedView
        ( A.LockedLocalMatrix(), localHeightBefore, j, localHeight, width );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::View1x2
( DistMatrixBase<T,MR,Star>& AL, DistMatrixBase<T,MR,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::View1x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_colAlignment = AL.ColAlignment();
    this->_colShift = AL.ColShift();
    this->_localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::LockedView1x2
( const DistMatrixBase<T,MR,Star>& AL, const DistMatrixBase<T,MR,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::LockedView1x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_colAlignment = AL.ColAlignment();
    this->_colShift = AL.ColShift();
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
elemental::DistMatrixBase<T,MR,Star>::View2x1
( DistMatrixBase<T,MR,Star>& AT,
  DistMatrixBase<T,MR,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::View2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_colAlignment = AT.ColAlignment();
    this->_colShift = AT.ColShift();
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
elemental::DistMatrixBase<T,MR,Star>::LockedView2x1
( const DistMatrixBase<T,MR,Star>& AT,
  const DistMatrixBase<T,MR,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::LockedView2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_colAlignment = AT.ColAlignment();
    this->_colShift = AT.ColShift();
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
elemental::DistMatrixBase<T,MR,Star>::View2x2
( DistMatrixBase<T,MR,Star>& ATL,
  DistMatrixBase<T,MR,Star>& ATR,
  DistMatrixBase<T,MR,Star>& ABL,
  DistMatrixBase<T,MR,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::View2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_colAlignment = ATL.ColAlignment();
    this->_colShift = ATL.ColShift();
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
elemental::DistMatrixBase<T,MR,Star>::LockedView2x2
( const DistMatrixBase<T,MR,Star>& ATL,
  const DistMatrixBase<T,MR,Star>& ATR,
  const DistMatrixBase<T,MR,Star>& ABL,
  const DistMatrixBase<T,MR,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::LockedView2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_colAlignment = ATL.ColAlignment();
    this->_colShift = ATL.ColShift();
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
elemental::DistMatrixBase<T,MR,Star>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    this->_height = height;
    this->_width = width;
    this->_localMatrix.ResizeTo
    ( LocalLength(height,this->ColShift(),this->Grid().Width()), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,MR,Star>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // columns within each process row
    const Grid& g = this->Grid();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();

    T u;
    if( g.MRRank() == ownerCol )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        u = this->GetLocalEntry(iLoc,j);
    }
    Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        this->SetLocalEntry(iLoc,j,u);
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
elemental::DistMatrixBase<T,MR,Star>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();    
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int c = this->Grid().Width();
    const int colShift = this->ColShift();

    if( shape == Lower )
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            int lastZeroRow = ( side==Left ? j-offset-1
                                           : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                int boundary = min( lastZeroRow+1, height );
                int numZeroRows = LocalLength( boundary, colShift, c );
#ifdef RELEASE
                T* thisCol = this->LocalBuffer(0,j);
                memset( thisCol, 0, numZeroRows*sizeof(T) );
#else
                for( int iLoc=0; iLoc<numZeroRows; ++iLoc )
                    this->SetLocalEntry(iLoc,j,0);
#endif
            }
        }
    }
    else
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            int firstZeroRow = ( side==Left ? max(j-offset+1,0)
                                            : max(j-offset+height-width+1,0) );
            int numNonzeroRows = LocalLength(firstZeroRow,colShift,c);
#ifdef RELEASE
            if( numNonzeroRows < localHeight )
            {
                T* thisCol = this->LocalBuffer(numNonzeroRows,j);
                memset( thisCol, 0, (localHeight-numNonzeroRows)*sizeof(T) );
            }
#else
            for( int iLoc=numNonzeroRows; iLoc<localHeight; ++iLoc )
                this->SetLocalEntry(iLoc,j,0);
#endif
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();    
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int c = this->Grid().Width();
    const int colShift = this->ColShift();

    if( shape == Upper )
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            int lastRow = ( side==Left ? j-offset : j-offset+height-width );
            int boundary = min( lastRow+1, height );
            int numRows = LocalLength( boundary, colShift, c );
#ifdef RELEASE
            T* thisCol = this->LocalBuffer(0,j);
            for( int iLoc=0; iLoc<numRows; ++iLoc )
                thisCol[iLoc] *= alpha;
#else
            for( int iLoc=0; iLoc<numRows; ++iLoc )
            {
                const T value = this->GetLocalEntry(iLoc,j);
                this->SetLocalEntry(iLoc,j,alpha*value);
            }
#endif
        }
    }
    else
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<width; ++j )
        {
            int firstRow = ( side==Left ? max(j-offset,0)
                                        : max(j-offset+height-width,0) );
            int numZeroRows = LocalLength( firstRow, colShift, c );
#ifdef RELEASE
            T* thisCol = this->LocalBuffer(numZeroRows,j);
            for( int iLoc=0; iLoc<(localHeight-numZeroRows); ++iLoc )
                thisCol[iLoc] *= alpha;
#else
            for( int iLoc=numZeroRows; iLoc<localHeight; ++iLoc )
            {
                const T value = this->GetLocalEntry(iLoc,j);
                this->SetLocalEntry(iLoc,j,alpha*value);
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
elemental::DistMatrixBase<T,MR,Star>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int c = this->Grid().Width();
    const int colShift = this->ColShift();

    this->SetToZero();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*c;
        if( i < width )
            this->SetLocalEntry(iLoc,i,1);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->Grid();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int bufSize = localHeight*width;

    this->_auxMemory.Require( bufSize );

    // Create random matrix on process row 0, then broadcast
    T* buffer = this->_auxMemory.Buffer();
    if( g.MCRank() == 0 )
    {
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                buffer[i+j*localHeight] = Random<T>();
    }
    Broadcast( buffer, bufSize, 0, g.MCComm() );

    // Unpack
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeight; ++i )
            this->SetLocalEntry(i,j,buffer[i+j*localHeight]);

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::SumOverCol()
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SumOverCol");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->Grid();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int localSize = max( localHeight*localWidth, MinCollectContrib );

    this->_auxMemory.Require( 2*localSize );
    T* buffer = this->_auxMemory.Buffer();
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[localSize];

    // Pack
#ifdef RELEASE
# ifdef _OPENMP
    #pragma omp parallel for
# endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* thisCol = this->LockedLocalBuffer(0,j);
        T* sendBufCol = &(sendBuf[j*localHeight]);
        memcpy( sendBufCol, thisCol, localHeight*sizeof(T) );
    }
#else
# ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
# endif
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            sendBuf[i+j*localHeight] = this->GetLocalEntry(i,j);
#endif

    // AllReduce sum
    AllReduce
    ( sendBuf, recvBuf, localSize, MPI_SUM, g.MCComm() );

    // Unpack
#ifdef RELEASE
# ifdef _OPENMP
    #pragma omp parallel for
# endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* recvBufCol = &(recvBuf[j*localHeight]);
        T* thisCol = this->LocalBuffer(0,j);
        memcpy( thisCol, recvBufCol, localHeight*sizeof(T) );
    }
#else
# ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
# endif
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->SetLocalEntry(i,j,recvBuf[i+j*localHeight]);
#endif

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::ConjugateTransposeFrom
( const DistMatrixBase<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ]::ConjugateTransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.RowAlignment();
            this->_colShift = 
                Shift( g.MRRank(), this->ColAlignment(), g.Width() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->ColAlignment() == A.RowAlignment() )
    {
        const int r = g.Height();

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            max(localHeight*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack 
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<localHeight; ++i )
                originalData[i+j*localHeight] = Conj( A.GetLocalEntry(j,i) );

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

            const int rowShift = Shift( k, colAlignmentOfA, r );
            const int localWidth = LocalLength( width, rowShift, r );

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &(data[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,rowShift+j*r);
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(i,rowShift+j*r,data[i+j*localHeight]);
#endif
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,* ]::ConjugateTransposeFrom" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendCol = (col+c+colAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-colAlignment) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            max(maxLocalHeight*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack the currently owned local data of A into the second buffer
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<localWidthOfA; ++i )
                secondBuffer[i+j*localWidthOfA] = Conj( A.GetLocalEntry(j,i) );

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

            const int rowShift = Shift( k, colAlignmentOfA, r );
            const int localWidth = LocalLength( width, rowShift, r );
#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &(data[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,rowShift+j*r);
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(i,rowShift+j*r,data[i+j*localHeight]);
#endif
        }

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,Star>::TransposeFrom
( const DistMatrixBase<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ]::TransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const Grid& g = this->Grid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.RowAlignment();
            this->_colShift = 
                Shift( g.MRRank(), this->ColAlignment(), g.Width() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->ColAlignment() == A.RowAlignment() )
    {
        const int r = g.Height();

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            max(localHeight*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack 
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<localHeight; ++i )
                originalData[i+j*localHeight] = A.GetLocalEntry(j,i);

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

            const int rowShift = Shift( k, colAlignmentOfA, r );
            const int localWidth = LocalLength( width, rowShift, r );

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &(data[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,rowShift+j*r);
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(i,rowShift+j*r,data[i+j*localHeight]);
#endif
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,* ]::TransposeFrom" << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendCol = (col+c+colAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-colAlignment) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            max(maxLocalHeight*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack the currently owned local data of A into the second buffer
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<localWidthOfA; ++i )
                secondBuffer[i+j*localWidthOfA] = A.GetLocalEntry(j,i);

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

            const int rowShift = Shift( k, colAlignmentOfA, r );
            const int localWidth = LocalLength( width, rowShift, r );
#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &(data[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,rowShift+j*r);
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(i,rowShift+j*r,data[i+j*localHeight]);
#endif
        }

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(g) );
    *A_VC_Star = A;

    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
    ( new DistMatrix<T,VR,Star>(true,this->ColAlignment(),g) );
    *A_VR_Star = *A_VC_Star;
    delete A_VC_Star.release(); // lowers memory highwater

    *this = *A_VR_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    if( A.Width() == 1 )
    {
        if( !this->Viewing() )
            this->ResizeTo( A.Height(), 1 );

        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int myCol = g.MRRank();
        const int rankCM = g.VCRank();
        const int rankRM = g.VRRank();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int colShiftOfA = A.ColShift();

        const int height = this->Height();
        const int maxLocalVectorHeight = MaxLocalLength(height,p);
        const int portionSize = max(maxLocalVectorHeight,MinCollectContrib);

        const int colShiftVR = Shift(rankRM,colAlignment,p);
        const int colShiftVCOfA = Shift(rankCM,colAlignmentOfA,p);
        const int sendRankRM = (rankRM+(p+colShiftVCOfA-colShiftVR)) % p;
        const int recvRankCM = (rankCM+(p+colShiftVR-colShiftVCOfA)) % p;
        const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

        this->_auxMemory.Require( (r+1)*portionSize );
        T* buffer = this->_auxMemory.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        // A[VC,* ] <- A[MC,* ]
        {
            const int shift = Shift(rankCM,colAlignmentOfA,p);
            const int offset = (shift-colShiftOfA) / r;
            const int thisLocalHeight = LocalLength(height,shift,p);

#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int i=0; i<thisLocalHeight; ++i )
                sendBuf[i] = A.GetLocalEntry(offset+i*c,0);
        }

        // A[VR,* ] <- A[VC,* ]
        SendRecv
        ( sendBuf, portionSize, sendRankRM, 0,
          recvBuf, portionSize, recvRankRM, MPI_ANY_TAG, g.VRComm() );

        // A[MR,* ] <- A[VR,* ]
        AllGather
        ( recvBuf, portionSize,
          sendBuf, portionSize, g.MCComm() );

        // Unpack
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &sendBuf[k*portionSize];

            const int shift = Shift(myCol+c*k,colAlignment,p);
            const int offset = (shift-this->ColShift()) / c;
            const int thisLocalHeight = LocalLength(height,shift,p);

            for( int i=0; i<thisLocalHeight; ++i )
                this->SetLocalEntry(offset+i*r,0,data[i]);
        }
            
        this->_auxMemory.Release();
    }
    else
    {
        auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
        ( new DistMatrix<T,VC,Star>(g) );
        *A_VC_Star = A;

        auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
        ( new DistMatrix<T,VR,Star>(true,this->ColAlignment(),g) );
        *A_VR_Star = *A_VC_Star;
        delete A_VC_Star.release(); // lowers memory highwater

        *this = *A_VR_Star;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR   = A;

    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(g) );
    *A_VC_Star = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
    ( new DistMatrix<T,VR,Star>(true,this->ColAlignment(),g) );
    *A_VR_Star = *A_VC_Star;
    delete A_VC_Star.release(); // lowers memory highwater

    *this = *A_VR_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MR,* ] = [MD,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MR,* ] = [* ,MD] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [MR,MC]");
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
            this->_colAlignment = A.ColAlignment();
            this->_colShift = 
                Shift( g.MRRank(), this->ColAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() )
    {
        const int r = g.Height();

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            max(localHeight*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack 
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &(originalData[j*localHeight]);
            memcpy( originalDataCol, ACol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeight; ++i )
                originalData[i+j*localHeight] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowAlignmentOfA = A.RowAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShift = Shift( k, rowAlignmentOfA, r );
            const int localWidth = LocalLength( width, rowShift, r );

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &(data[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,rowShift+j*r);
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(i,rowShift+j*r,data[i+j*localHeight]);
#endif
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,* ] <- [MR,MC]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            max(maxLocalHeight*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack the currently owned local data of A into the second buffer
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* secondBufferCol = &(secondBuffer[j*localHeightOfA]);
            memcpy( secondBufferCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                secondBuffer[i+j*localHeightOfA] = A.GetLocalEntry(i,j);
#endif

        // Perform the SendRecv: puts the new data into the first buffer
        SendRecv
        ( secondBuffer, portionSize, sendCol, 0,
          firstBuffer,  portionSize, recvCol, MPI_ANY_TAG, g.MRComm() );

        // Use the output of the SendRecv as input to the AllGather
        AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

        // Unpack the contents of each member of the process col
        const int rowAlignmentOfA = A.RowAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int rowShift = Shift( k, rowAlignmentOfA, r );
            const int localWidth = LocalLength( width, rowShift, r );
#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &(data[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,rowShift+j*r);
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(i,rowShift+j*r,data[i+j*localHeight]);
#endif
        }

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [MR,* ]");
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
            this->_colAlignment = A.ColAlignment();
            this->_colShift = A.ColShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() )
    {
        this->_localMatrix = A.LockedLocalMatrix();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,* ] <- [MR,* ]." << endl;
#endif
        const int rank = g.MRRank();
        const int c = g.Width();

        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (rank+c+colAlignment-colAlignmentOfA) % c;
        const int recvRank = (rank+c+colAlignmentOfA-colAlignment) % c;

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();

        const int sendSize = localHeightOfA * width;
        const int recvSize = localHeight * width;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<width; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* sendBufferCol = &(sendBuffer[j*localHeightOfA]);
            memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i+j*localHeightOfA] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.MRComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<width; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*localHeight]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    DistMatrix<T,MR,MC> A_MR_MC(g);

    A_MR_MC = A;
    *this   = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    DistMatrix<T,VR,Star> A_VR_Star(true,this->ColAlignment(),g);

    A_VR_Star = A;
    *this     = A_VR_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    DistMatrix<T,MR,MC> A_MR_MC(g);

    A_MR_MC = A;
    *this   = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.VCRank() == 0 )
    {
        cerr << 
          "[MR,* ] <- [VR,* ] potentially causes a large amount of cache-"
          "thrashing. If possible avoid it by performing the redistribution "
          "with a (conjugate)-transpose: " << endl <<
          "  [* ,MR].(Conjugate)TransposeFrom([VR,* ])" << endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.ColAlignment() % g.Width();
            this->_colShift = 
                Shift( g.MRRank(), this->ColAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() % g.Width() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(height,p);

        const int portionSize = 
            max(maxLocalHeightOfA*width,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<width; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &(originalData[j*localHeightOfA]);
            memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                originalData[i+j*localHeightOfA] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int colShift = this->ColShift();
        const int colAlignmentOfA = A.ColAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShiftOfA = Shift( col+c*k, colAlignmentOfA, p );
            const int colOffset = (colShiftOfA-colShift) / c;
            const int localHeight = LocalLength( height, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(colOffset+i*r,j,data[i+j*localHeight]);
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,* ] <- [VR,* ]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int rank = g.VRRank();

        // Perform the SendRecv to make A have the same colAlignment
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int colShift = this->ColShift();

        const int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(height,p);

        const int portionSize = 
            max(maxLocalHeightOfA*width,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<width; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* secondBufferCol = &(secondBuffer[j*localHeightOfA]);
            memcpy( secondBufferCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                secondBuffer[i+j*localHeightOfA] = A.GetLocalEntry(i,j);
#endif

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

            const int colShiftOfA = Shift( col+c*k, colAlignment, p );
            const int colOffset = (colShiftOfA-colShift) / c;
            const int localHeight = LocalLength( height, colShiftOfA, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(colOffset+i*r,j,data[i+j*localHeight]);
        }

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
    ( new DistMatrix<T,Star,VC>(g) );
    *A_Star_VC = A;

    auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
    ( new DistMatrix<T,MR,MC>(true,false,this->ColAlignment(),0,g) );
    *A_MR_MC = *A_Star_VC;
    delete A_Star_VC.release(); // lowers memory highwater

    *this = *A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,Star>&
elemental::DistMatrixBase<T,MR,Star>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int c = g.Width();
    const int colShift = this->ColShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->SetLocalEntry(i,j,A.GetLocalEntry(colShift+i*c,j));
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class elemental::DistMatrixBase<int,   MR,Star>;
template class elemental::DistMatrixBase<float, MR,Star>;
template class elemental::DistMatrixBase<double,MR,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,MR,Star>;
template class elemental::DistMatrixBase<dcomplex,MR,Star>;
#endif

