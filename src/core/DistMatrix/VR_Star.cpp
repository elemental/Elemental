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
elemental::DistMatrixBase<T,VR,Star>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Print");
#endif
    const Grid& g = this->GetGrid();
    if( g.VRRank() == 0 && s != "" )
        cout << s << endl;

    const int height      = this->Height();
    const int width       = this->Width();
    const int localHeight = this->LocalHeight();
    const int p           = g.Size();
    const int colShift    = this->ColShift();

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
    for( int i=0; i<localHeight; ++i )
        for( int j=0; j<width; ++j )
            sendBuf[colShift+i*p+j*height] = this->LocalEntry(i,j);

    // If we are the root, allocate a receive buffer
    vector<T> recvBuf;
    if( g.VRRank() == 0 )
        recvBuf.resize( height*width );

    // Sum the contributions and send to the root
    Reduce( &sendBuf[0], &recvBuf[0], height*width, MPI_SUM, 0, g.VCComm() );

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
    Barrier( g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.RowAlignment();
    this->_colShift = Shift( g.VRRank(), this->ColAlignment(), g.Size() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = Shift( g.VRRank(), this->ColAlignment(), g.Size() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([MR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = Shift( g.VRRank(), this->ColAlignment(), g.Size() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([* ,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.RowAlignment();
    this->_colShift = Shift( g.VRRank(), this->ColAlignment(), g.Size() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith(DistMatrix[VR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_colShift = A.ColShift();
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith(DistMatrix[* ,VR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_colShift = A.RowShift();
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,MR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,Star,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,VR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,Star,VR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::View
( DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View(A)");
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
elemental::DistMatrixBase<T,VR,Star>::LockedView
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView(A)");
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
elemental::DistMatrixBase<T,VR,Star>::View
( DistMatrixBase<T,VR,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View(A,i,j,height,width)");
    this->AssertFreeColAlignment();
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

        this->_colAlignment = (A.ColAlignment()+i) % size;
        this->_colShift = Shift( rowMajorRank, this->ColAlignment(), size );

        const int localHeightBefore = LocalLength( i, A.ColShift(), size );
        const int localHeight = LocalLength( height, this->ColShift(), size );

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
elemental::DistMatrixBase<T,VR,Star>::LockedView
( const DistMatrixBase<T,VR,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView");
    this->AssertFreeColAlignment();
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

        this->_colAlignment = (A.ColAlignment()+i) % size;
        this->_colShift = Shift( rowMajorRank, this->ColAlignment(), size );

        const int localHeightBefore = LocalLength( i, A.ColShift(), size );
        const int localHeight = LocalLength( height, this->ColShift(), size );

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
elemental::DistMatrixBase<T,VR,Star>::View1x2
( DistMatrixBase<T,VR,Star>& AL,
  DistMatrixBase<T,VR,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View1x2");
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
elemental::DistMatrixBase<T,VR,Star>::LockedView1x2
( const DistMatrixBase<T,VR,Star>& AL,
  const DistMatrixBase<T,VR,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView1x2");
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
elemental::DistMatrixBase<T,VR,Star>::View2x1
( DistMatrixBase<T,VR,Star>& AT,
  DistMatrixBase<T,VR,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View2x1");
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
    this->_localMatrix.View2x1( AT.LocalMatrix(), AB.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::LockedView2x1
( const DistMatrixBase<T,VR,Star>& AT,
  const DistMatrixBase<T,VR,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView2x1");
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
elemental::DistMatrixBase<T,VR,Star>::View2x2
( DistMatrixBase<T,VR,Star>& ATL,
  DistMatrixBase<T,VR,Star>& ATR,
  DistMatrixBase<T,VR,Star>& ABL,
  DistMatrixBase<T,VR,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ABL.Width() + ABR.Width();
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
elemental::DistMatrixBase<T,VR,Star>::LockedView2x2
( const DistMatrixBase<T,VR,Star>& ATL,
  const DistMatrixBase<T,VR,Star>& ATR,
  const DistMatrixBase<T,VR,Star>& ABL,
  const DistMatrixBase<T,VR,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ABL.Width() + ABR.Width();
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
elemental::DistMatrixBase<T,VR,Star>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    const Grid& g = this->GetGrid();
    this->_height = height;
    this->_width  = width;
    this->_localMatrix.ResizeTo
    ( LocalLength(height,this->ColShift(),g.Size()) ,width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,VR,Star>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const Grid& g = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    T u;
    if( g.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        u = this->LocalEntry(iLoc,j);
    }
    Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        this->LocalEntry(iLoc,j) = u;
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
elemental::DistMatrixBase<T,VR,Star>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->GetGrid().Size();
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
                int numZeroRows = LocalLength( boundary, colShift, p );
#ifdef RELEASE
                T* thisCol = &(this->LocalEntry(0,j));
                memset( thisCol, 0, numZeroRows*sizeof(T) );
#else
                for( int iLoc=0; iLoc<numZeroRows; ++iLoc )
                    this->LocalEntry(iLoc,j) = (T)0;
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
            int numNonzeroRows = LocalLength(firstZeroRow,colShift,p);
#ifdef RELEASE
            T* thisCol = &(this->LocalEntry(numNonzeroRows,j));
            memset( thisCol, 0, (localHeight-numNonzeroRows)*sizeof(T) );
#else
            for( int iLoc=numNonzeroRows; iLoc<localHeight; ++iLoc )
                this->LocalEntry(iLoc,j) = (T)0;
#endif
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->GetGrid().Size();
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
            int numRows = LocalLength( boundary, colShift, p );
#ifdef RELEASE
            T* thisCol = &(this->LocalEntry(0,j));
            for( int iLoc=0; iLoc<numRows; ++iLoc )
                thisCol[iLoc] *= alpha;
#else
            for( int iLoc=0; iLoc<numRows; ++iLoc )
                this->LocalEntry(iLoc,j) *= alpha;
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
            int numZeroRows = LocalLength(firstRow,colShift,p);
#ifdef RELEASE
            T* thisCol = &(this->LocalEntry(numZeroRows,j));
            for( int iLoc=0; iLoc<(localHeight-numZeroRows); ++iLoc )
                thisCol[iLoc] *= alpha;
#else
            for( int iLoc=numZeroRows; iLoc<localHeight; ++iLoc )
                this->LocalEntry(iLoc,j) *= alpha;
#endif
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::SetToIdentity()
{   
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->GetGrid().Size();
    const int colShift = this->ColShift();

    this->SetToZero();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*p;
        if( i < width )
            this->LocalEntry(iLoc,i) = (T)1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = Random<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,VC,Star> A_VC_Star(g);

    A_VC_Star = A;
    *this = A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,VC,Star> A_VC_Star(g);

    A_VC_Star = A;
    *this = A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR = A;

    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(g) );
    *A_VC_Star = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[VR,* ] = [MD,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[VR,* ] = [* ,MD] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.ColAlignment();
            this->_colShift = 
                Shift( g.VRRank(), this->ColAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,r);
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
            const int thisColShift = Shift(thisRank,colAlignment,p);
            const int thisColOffset = (thisColShift-colShiftOfA) / c;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                        A.LocalEntry(thisColOffset+i*r,j);
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

            const int thisRowShift = Shift(k,rowAlignmentOfA,r);
            const int thisLocalWidth = LocalLength(width,thisRowShift,r);

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* dataCol = &(data[j*localHeight]);
                T* thisCol = &(this->LocalEntry(0,thisRowShift+j*r));
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,thisRowShift+j*r) = 
                        data[i+j*localHeight];
#endif
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [VR,* ] <- [MR,MC]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int col = g.MRRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendCol = (col+c+(colAlignment%c)-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-(colAlignment%c)) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,r);
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
            const int thisColShift = Shift(thisRank,colAlignment,p);
            const int thisColOffset = (thisColShift-colShiftOfA) / c;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                        A.LocalEntry(thisColOffset+i*r,j);
        }

        // AllToAll to gather all of the unaligned [VR,*] data into firstBuffer
        AllToAll
        ( secondBuffer, portionSize,
          firstBuffer,  portionSize, g.MCComm() );

        // SendRecv: properly align the [VR,*] via a trade in the row
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

            const int thisRowShift = Shift(k,rowAlignmentOfA,r);
            const int thisLocalWidth = LocalLength(width,thisRowShift,r);

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* dataCol = &(data[j*localHeight]);
                T* thisCol = &(this->LocalEntry(0,thisRowShift+j*r));
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,thisRowShift+j*r) = 
                        data[i+j*localHeight];
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
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.ColAlignment();
            this->_colShift = 
                Shift( g.VRRank(), this->ColAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int colShift = this->ColShift();
        const int colShiftOfA = A.ColShift();
        const int colOffset = (colShift-colShiftOfA) / c;

        const int width = this->Width();
        const int localHeight = this->LocalHeight();

#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = A.LocalEntry(colOffset+i*r,j);
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [VR,* ] <- [MR,* ]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int row = g.MCRank();
        const int col = g.MRRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[VR,*] within our process row to fix alignments.
        const int sendCol = (col+c+(colAlignment%c)-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-(colAlignment%c)) % c;
        const int sendRank = sendCol + c*row;

        const int sendColShift = Shift( sendRank, colAlignment, p );
        const int sendColOffset = (sendColShift-colShiftOfA) / c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfSend = LocalLength(height,sendColShift,p);
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = maxLocalHeight * width;

        this->_auxMemory.Require( 2*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[portionSize];

        // Pack
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfSend; ++i )
                sendBuffer[i+j*localHeightOfSend] =
                    A.LocalEntry(sendColOffset+i*r,j);

        // Communicate
        SendRecv
        ( sendBuffer, portionSize, sendCol, 0,
          recvBuffer, portionSize, recvCol, MPI_ANY_TAG, g.MRComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<width; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = &(this->LocalEntry(0,j));
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,MR,MC> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    
    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();
    const int rankCM = g.VCRank();
    const int rankRM = g.VRRank();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localHeightOfA = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength(height,p);

    const int portionSize = maxLocalHeight * width;

    const int colShift = this->ColShift();
    const int colShiftOfA = A.ColShift();

    // Compute which rowmajor rank has the colShift equal to our colShiftOfA
    const int sendRankRM = (rankRM+(p+colShiftOfA-colShift)) % p;

    // Compute which rowmajor rank has the A colShift that we need
    const int recvRankCM = (rankCM+(p+colShift-colShiftOfA)) % p;
    const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

    this->_auxMemory.Require( 2*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[portionSize];

    // Pack
#ifdef RELEASE
# ifdef _OPENMP
    #pragma omp parallel for
# endif
    for( int j=0; j<width; ++j )
    {
        const T* ACol = &(A.LocalEntry(0,j));
        T* sendBufferCol = &(sendBuffer[j*localHeightOfA]);
        memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
    }
#else
# ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
# endif
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            sendBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);
#endif

    // Communicate
    SendRecv
    ( sendBuffer, portionSize, sendRankRM, 0,
      recvBuffer, portionSize, recvRankRM, MPI_ANY_TAG, g.VRComm() );

    // Unpack
#ifdef RELEASE
# ifdef _OPENMP
    #pragma omp parallel for
# endif
    for( int j=0; j<width; ++j )
    {
        const T* recvBufferCol = &(recvBuffer[j*localHeight]);
        T* thisCol = &(this->LocalEntry(0,j));
        memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
    }
#else
# ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
# endif
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];
#endif

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,MR,MC> A_MR_MC(g);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
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
        const Grid& g = this->GetGrid();
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [VR,* ] <- [VR,* ]." << endl;
#endif
        const int rank = g.VRRank();
        const int p = g.Size();

        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

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
            const T* ACol = &(A.LocalEntry(0,j));
            T* sendBufferCol = &(sendBuffer[j*localHeightOfA]);
            memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);
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
        for( int j=0; j<width; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = &(this->LocalEntry(0,j));
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(g) );
    *A_MC_MR = A;

    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(g) );
    *A_VC_Star = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int p = this->GetGrid().Size();
    const int colShift = this->ColShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = A.LocalEntry(colShift+i*p,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::SumScatterFrom
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SumScatterFrom( [MR,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.VCRank() == 0 )
    {
        cerr <<
          "[VR,* ]::SumScatterFrom([MR,* ]) potentially causes a large amount "
          "of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MR,* ] matrix instead." << endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.ColAlignment();
            this->_colShift = 
                Shift( g.VRRank(), this->ColAlignment(), g.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();
        const int colAlignment = this->ColAlignment();
        const int colShiftOfA = A.ColShift();

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalHeight = MaxLocalLength( height, p );

        const int recvSize = max(maxLocalHeight*localWidth,MinCollectContrib);
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
            const int thisColShift = Shift( thisRank, colAlignment, p );
            const int thisColOffset = (thisColShift-colShiftOfA) / r;
            const int thisLocalHeight = LocalLength( height, thisColShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                        A.LocalEntry(thisColOffset+i*r,j);
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
            T* thisCol = &(this->LocalEntry(0,j));
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];
#endif

        this->_auxMemory.Release();
    }
    else
    {
        throw logic_error
              ( "Unaligned [VR,* ]::ReduceScatterFrom( [MR,* ] ) is not "
                "yet implemented." );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::SumScatterUpdate
( T alpha, const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SumScatterUpdate( [MR,* ] )");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.VCRank() == 0 )
    {
        cerr <<
          "[VR,* ]::SumScatterUpdate([MR,* ]) potentially causes a large amount"
          " of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [MR,* ] matrix instead." << endl;
    }
#endif
    if( this->ColAlignment() % g.Width() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();
        const int colAlignment = this->ColAlignment();
        const int colShiftOfA = A.ColShift();

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalHeight = MaxLocalLength( height, p );

        const int recvSize = max(maxLocalHeight*localWidth,MinCollectContrib);
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
            const int thisColShift = Shift( thisRank, colAlignment, p );
            const int thisColOffset = (thisColShift-colShiftOfA) / c;
            const int thisLocalHeight = LocalLength( height, thisColShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                        A.LocalEntry(thisColOffset+i*r,j);
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
            T* thisCol = &(this->LocalEntry(0,j));
            for( int i=0; i<localHeight; ++i )
                thisCol[i] += alpha*recvBufferCol[i];
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) += alpha*recvBuffer[i+j*localHeight];
#endif

        this->_auxMemory.Release();
    }
    else
    {
        throw logic_error
              ( "Unaligned [VR,* ]::ReduceScatterUpdate( [MR,* ] ) is not "
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
elemental::DistMatrix<R,VR,Star>::SetToRandomHPD()
{   
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->GetGrid().Size();
    const int colShift = this->ColShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*p;
        if( i < width )
            this->LocalEntry(iLoc,i) += (R)this->Width();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,VR,Star>::SetToRandomHPD()
{   
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->GetGrid().Size();
    const int colShift = this->ColShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*p;
        if( i < width )
            this->LocalEntry(iLoc,i) = 
                real(this->LocalEntry(iLoc,i)) + (R)this->Width();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,VR,Star>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const Grid& g = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    R u;
    if( g.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        u = real(this->LocalEntry(iLoc,j));
    }
    Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,VR,Star>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const Grid& g = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    R u;
    if( g.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        u = imag(this->LocalEntry(iLoc,j));
    }
    Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,VR,Star>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        real(this->LocalEntry(iLoc,j)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,VR,Star>::SetImag
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        imag(this->LocalEntry(iLoc,j)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrixBase<int,   VR,Star>;
template class elemental::DistMatrixBase<float, VR,Star>;
template class elemental::DistMatrixBase<double,VR,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,VR,Star>;
template class elemental::DistMatrixBase<dcomplex,VR,Star>;
#endif

template class elemental::DistMatrix<int,   VR,Star>;
template class elemental::DistMatrix<float, VR,Star>;
template class elemental::DistMatrix<double,VR,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,VR,Star>;
template class elemental::DistMatrix<dcomplex,VR,Star>;
#endif

