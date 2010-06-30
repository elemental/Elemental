/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
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
elemental::DistMatrixBase<T,Star,MR>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::Print");
#endif
    const Grid& g = this->GetGrid();
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
        T* sendBuf = new T[height*width];
        for( int i=0; i<height*width; ++i )
            sendBuf[i] = (T)0;
        for( int i=0; i<height; ++i )
            for( int j=0; j<localWidth; ++j )
                sendBuf[i+(rowShift+j*c)*height] = this->LocalEntry(i,j);

        // If we are the root, fill the receive buffer
        T* recvBuf = 0;
        if( g.MRRank() == 0 )
        {
            recvBuf = new T[height*width];
            for( int i=0; i<height*width; ++i )
                recvBuf[i] = (T)0;
        }

        // Sum the contributions and send to the root
        Reduce
        ( sendBuf, recvBuf, height*width, MPI_SUM, 0, g.MRComm() );
        delete[] sendBuf;

        if( g.MRRank() == 0 )
        {
            // Print the data
            for( int i=0; i<height; ++i )
            {
                for( int j=0; j<width; ++j )
                    cout << recvBuf[i+j*height] << " ";
                cout << endl;
            }
            cout << endl;
            delete[] recvBuf;
        }
    }
    Barrier( g.VCComm() );

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
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.ColAlignment() % g.Width();
    this->_rowShift = 
        Shift( g.MRRank(), this->RowAlignment(), g.Width() );
    this->_constrainedRowAlignment = true;
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
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.RowAlignment() % g.Width();
    this->_rowShift = 
        Shift( g.MRRank(), this->RowAlignment(), g.Width() );
    this->_constrainedRowAlignment = true;
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
        const Grid& g = this->GetGrid();
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
        const Grid& g = this->GetGrid();
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
    ( height, LocalLength(width,this->RowShift(),this->GetGrid().Width()) );
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
    const Grid& g = this->GetGrid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    T u;
    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = this->LocalEntry(i,jLoc);
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
    const Grid& g = this->GetGrid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        this->LocalEntry(i,jLoc) = u;
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
    const Grid& g = this->GetGrid();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int c = g.Width();
    const int rowShift = this->RowShift();

    if( shape == Lower )
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*c;
            int firstNonzero_i;
            if( side == Left )
                firstNonzero_i = max(j-offset,0);
            else
                firstNonzero_i = max(j-offset+height-width,0);

            const int boundary = min(height,firstNonzero_i);
#ifdef RELEASE
            T* thisCol = &(this->LocalEntry(0,jLoc));
            memset( thisCol, 0, boundary*sizeof(T) );
#else
            for( int i=0; i<boundary; ++i )
                this->LocalEntry(i,jLoc) = (T)0;
#endif
        }
    }
    else
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*c;
            int firstZero_i;
            if( side == Left )
                firstZero_i = max(j-offset+1,0);
            else
                firstZero_i = max(j-offset+height-width+1,0);
#ifdef RELEASE
            T* thisCol = &(this->LocalEntry(firstZero_i,jLoc));
            memset( thisCol, 0, (height-firstZero_i)*sizeof(T) );
#else
            for( int i=firstZero_i; i<height; ++i )
                this->LocalEntry(i,jLoc) = (T)0;
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
    const int c = this->GetGrid().Width();
    const int rowShift = this->RowShift();

    this->SetToZero();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*c;
        if( j < height )
            this->LocalEntry(j,jLoc) = (T)1;
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
    const Grid& g = this->GetGrid();
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
#ifdef RELEASE
    for( int j=0; j<localWidth; ++j )
    {
        const T* bufferCol = &(buffer[j*height]);
        T* thisCol = &(this->LocalEntry(0,j));
        memcpy( thisCol, bufferCol, height*sizeof(T) );
    }
#else
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<height; ++i )
            this->LocalEntry(i,j) = buffer[i+j*height];
#endif

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
#ifdef RELEASE
    for( int j=0; j<localWidth; ++j )
    {
        const T* thisCol = &(this->LocalEntry(0,j));
        T* sendBufCol = &(sendBuf[j*localHeight]);
        memcpy( sendBufCol, thisCol, localHeight*sizeof(T) );
    }
#else
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            sendBuf[i+j*localHeight] = this->LocalEntry(i,j);
#endif

    // AllReduce col
    AllReduce
    ( sendBuf, recvBuf, localSize, MPI_SUM, this->GetGrid().MCComm() );

    // Unpack
#ifdef RELEASE
    for( int j=0; j<localWidth; ++j )
    {
        const T* recvBufCol = &(recvBuf[j*localHeight]);
        T* thisCol = &(this->LocalEntry(0,j));
        memcpy( thisCol, recvBufCol, localHeight*sizeof(T) );
    }
#else
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = recvBuf[i+j*localHeight];
#endif

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
    const Grid& g = this->GetGrid();
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
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<height; ++i )
                originalData[i+j*height] = Conj( A.LocalEntry(j,i) );

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int colAlignmentOfA = A.ColAlignment();
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShiftOfA = Shift( col+k*c, colAlignmentOfA, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, colShiftOfA, p );

#ifdef RELEASE
            for( int j=0; j<localWidth; ++j )      
            {
                const T* dataCol = &(data[j*height]);
                T* thisCol = &(this->LocalEntry(0,rowOffset+j*r));
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
#else
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    this->LocalEntry(i,rowOffset+j*r) = data[i+j*height];
#endif
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [* ,MR].ConjugateTransposeFrom[VR,* ]." << endl;
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
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<height; ++i )
                secondBuffer[i+j*height] = Conj( A.LocalEntry(j,i) );

        // Perform the SendRecv: puts the new data into the first buffer
        SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, MPI_ANY_TAG, g.VRComm() );

        // Use the SendRecv as input to the AllGather
        AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

        // Unpack
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int colShiftOfA = Shift( col+c*k, rowAlignment, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, colShiftOfA, p );

#ifdef RELEASE
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &(data[j*height]);
                T* thisCol = &(this->LocalEntry(0,rowOffset+j*r));
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
#else
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    this->LocalEntry(i,rowOffset+j*r) = data[i+j*height];
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
    const Grid& g = this->GetGrid();
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
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<height; ++i )
                originalData[i+j*height] = A.LocalEntry(j,i);

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int colAlignmentOfA = A.ColAlignment();
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShiftOfA = Shift( col+k*c, colAlignmentOfA, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, colShiftOfA, p );

#ifdef RELEASE
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &(data[j*height]);
                T* thisCol = &(this->LocalEntry(0,rowOffset+j*r));
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
#else
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    this->LocalEntry(i,rowOffset+j*r) = data[i+j*height];
#endif
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [* ,MR].TransposeFrom[VR,* ]." << endl;
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
        for( int j=0; j<localHeightOfA; ++j )
            for( int i=0; i<height; ++i )
                secondBuffer[i+j*height] = A.LocalEntry(j,i);

        // Perform the SendRecv: puts the new data into the first buffer
        SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, MPI_ANY_TAG, g.VRComm() );

        // Use the SendRecv as input to the AllGather
        AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

        // Unpack
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int colShiftOfA = Shift( col+c*k, rowAlignment, p );
            const int rowOffset = (colShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, colShiftOfA, p );

#ifdef RELEASE
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &(data[j*height]);
                T* thisCol = &(this->LocalEntry(0,rowOffset+j*r));
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
#else
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    this->LocalEntry(i,rowOffset+j*r) = data[i+j*height];
#endif
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
    if( A.Height() != 1 && A.GetGrid().VCRank() == 0 )
    {
        cout << 
          "The matrix redistribution [* ,MR] <- [MC,MR] potentially causes a "
          "large amount of cache-thrashing. If possible, avoid it by "
          "performing the redistribution with a (conjugate)-transpose:"
          << endl <<
          "  [MR,* ].(Conjugate)TransposeFrom([MC,MR])" << endl;
    }
#endif
    const Grid& g = this->GetGrid();
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
                for( int j=0; j<localWidth; ++j )
                    bcastBuf[j] = this->LocalEntry(0,j);
            }

            // Communicate
            Broadcast
            ( bcastBuf, localWidth, A.ColAlignment(), g.MCComm() );

            // Unpack
            for( int j=0; j<localWidth; ++j )
                this->LocalEntry(0,j) = bcastBuf[j];

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
#ifdef RELEASE
            for( int j=0; j<localWidth; ++j )
            {
                const T* ACol = &(A.LocalEntry(0,j));
                T* originalDataCol = &(originalData[j*localHeightOfA]);
                memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
            }
#else
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);
#endif

            // Communicate
            AllGather
            ( originalData, portionSize,
              gatheredData, portionSize, g.MCComm() );

            // Unpack
            const int colAlignmentOfA = A.ColAlignment();
            for( int k=0; k<r; ++k )
            {
                const T* data = &gatheredData[k*portionSize];

                const int colShift = Shift( k, colAlignmentOfA, r );
                const int localHeight = LocalLength( height, colShift, r );

                for( int j=0; j<localWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        this->LocalEntry(colShift+i*r,j) = 
                            data[i+j*localHeight];
            }

            this->_auxMemory.Release();
        }
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [* ,MR] <- [MC,MR]." << endl;
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
                for( int j=0; j<localWidthOfA; ++j )
                    sendBuf[j] = A.LocalEntry(0,j);

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
            for( int j=0; j<localWidth; ++j )
                this->LocalEntry(0,j) = bcastBuf[j];

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
#ifdef RELEASE
            for( int j=0; j<localWidthOfA; ++j )
            {
                const T* ACol = &(A.LocalEntry(0,j));
                T* secondBufferCol = &(secondBuffer[j*localHeightOfA]);
                memcpy( secondBufferCol, ACol, localHeightOfA*sizeof(T) );
            }
#else
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    secondBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);
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
            const int colAlignmentOfA = A.ColAlignment();
            for( int k=0; k<r; ++k )
            {
                const T* data = &secondBuffer[k*portionSize];

                const int colShift = Shift( k, colAlignmentOfA, r );
                const int localHeight = LocalLength( height, colShift, r );
                for( int j=0; j<localWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        this->LocalEntry(colShift+i*r,j) = 
                            data[i+j*localHeight];
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
    const Grid& g = this->GetGrid();
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
    const Grid& g = this->GetGrid();
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
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [* ,MR] <- [* ,MR]." << endl;
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
#ifdef RELEASE
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = &(A.LocalEntry(0,j));
            T* sendBufferCol = &(sendBuffer[j*height]);
            memcpy( sendBufferCol, ACol, height*sizeof(T) );
        }
#else
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<height; ++i )
                sendBuffer[i+j*height] = A.LocalEntry(i,j);
#endif

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.MRComm() );

        // Unpack
#ifdef RELEASE
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*height]);
            T* thisCol = &(this->LocalEntry(0,j));
            memcpy( thisCol, recvBufferCol, height*sizeof(T) );
        }
#else
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->LocalEntry(i,j) = recvBuffer[i+j*height];
#endif

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
    const Grid& g = this->GetGrid();
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
    const Grid& g = this->GetGrid();
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
    const Grid& g = this->GetGrid();
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
    const Grid& g = this->GetGrid();
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
    const Grid& g = this->GetGrid();
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
    const Grid& g = this->GetGrid();
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
    const Grid& g = this->GetGrid();
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
#ifdef RELEASE
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = &(A.LocalEntry(0,j));
            T* originalDataCol = &(originalData[j*height]);
            memcpy( originalDataCol, ACol, height*sizeof(T) );
        }
#else
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<height; ++i )
                originalData[i+j*height] = A.LocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowShift = this->RowShift();
        const int rowAlignmentOfA = A.RowAlignment();
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShiftOfA = Shift( col+k*c, rowAlignmentOfA, p );
            const int rowOffset = (rowShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, rowShiftOfA, p );

#ifdef RELEASE
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &(data[j*height]);
                T* thisCol = &(this->LocalEntry(0,rowOffset+j*r));
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
#else
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    this->LocalEntry(i,rowOffset+j*r) = data[i+j*height];
#endif
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [* ,MR] <- [* ,VR]." << endl;
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
#ifdef RELEASE
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = &(A.LocalEntry(0,j));
            T* secondBufferCol = &(secondBuffer[j*height]);
            memcpy( secondBufferCol, ACol, height*sizeof(T) );
        }
#else
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<height; ++i )
                secondBuffer[i+j*height] = A.LocalEntry(i,j);
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
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int rowShiftOfA = Shift( col+c*k, rowAlignment, p );
            const int rowOffset = (rowShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, rowShiftOfA, p );

#ifdef RELEASE
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &(data[j*height]);
                T* thisCol = &(this->LocalEntry(0,rowOffset+j*r));
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
#else
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    this->LocalEntry(i,rowOffset+j*r) = data[i+j*height];
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

    const int c = this->GetGrid().Width();
    const int rowShift = this->RowShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
#ifdef RELEASE
    for( int j=0; j<localWidth; ++j )
    {
        const T* ACol = &(A.LocalEntry(0,rowShift+j*c));
        T* thisCol = &(this->LocalEntry(0,j));
        memcpy( thisCol, ACol, localHeight*sizeof(T) );
    }
#else
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = A.LocalEntry(i,rowShift+j*c);
#endif

#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

//----------------------------------------------------------------------------//
// DistMatrix                                                                 //
//----------------------------------------------------------------------------//

template<typename R>
void
elemental::DistMatrix<R,Star,MR>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int height     = this->Height();
    const int localWidth = this->LocalWidth();
    const int c          = this->GetGrid().Width();
    const int rowShift   = this->RowShift();

    this->SetToRandom();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*c;
        if( j < height )
            this->LocalEntry(j,jLoc) += (R)this->Width();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,Star,MR>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int height     = this->Height();
    const int localWidth = this->LocalWidth();
    const int c          = this->GetGrid().Width();
    const int rowShift   = this->RowShift();

    this->SetToRandom();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*c;
        if( j < height )
            this->LocalEntry(j,jLoc) = 
                real(this->LocalEntry(j,jLoc)) + (R)this->Width();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,MR>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // column within each process row
    const Grid& g = this->GetGrid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    R u;
    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = real(this->LocalEntry(i,jLoc));
    }
    Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,MR>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // column within each process row
    const Grid& g = this->GetGrid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    R u;
    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = imag(this->LocalEntry(i,jLoc));
    }
    Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,MR>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        real(this->LocalEntry(i,jLoc)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,MR>::SetImag
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        imag(this->LocalEntry(i,jLoc)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrixBase<int,   Star,MR>;
template class elemental::DistMatrixBase<float, Star,MR>;
template class elemental::DistMatrixBase<double,Star,MR>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,Star,MR>;
template class elemental::DistMatrixBase<dcomplex,Star,MR>;
#endif

template class elemental::DistMatrix<int,   Star,MR>;
template class elemental::DistMatrix<float, Star,MR>;
template class elemental::DistMatrix<double,Star,MR>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,Star,MR>;
template class elemental::DistMatrix<dcomplex,Star,MR>;
#endif

