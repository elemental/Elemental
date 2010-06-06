/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
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
elemental::DistMatrixBase<T,Star,Star>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Print");
#endif
    const Grid& grid = this->GetGrid();
    if( grid.VCRank() == 0 && s != "" )
        cout << s << endl;

    const int height = this->Height();
    const int width  = this->Width();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( grid.VCRank() == 0 )
    {
        for( int i=0; i<height; ++i )
        {
            for( int j=0; j<width; ++j )
                cout << this->LocalEntry(i,j) << " ";
            cout << endl;
        }
        cout << endl;
    }
    Barrier( grid.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::View
( DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View");
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_localMatrix.View( A.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::LockedView
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView");
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_localMatrix.LockedView( A.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::View
( DistMatrixBase<T,Star,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View");
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    this->_localMatrix.View( A.LocalMatrix(), i, j, height, width );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::LockedView
( const DistMatrixBase<T,Star,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView");
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    this->_localMatrix.LockedView( A.LockedLocalMatrix(), i, j, height, width );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::View1x2
( DistMatrixBase<T,Star,Star>& AL,
  DistMatrixBase<T,Star,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View1x2");
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::LockedView1x2
( const DistMatrixBase<T,Star,Star>& AL,
  const DistMatrixBase<T,Star,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView1x2");
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_localMatrix.LockedView1x2
    ( AL.LockedLocalMatrix(), 
      AR.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::View2x1
( DistMatrixBase<T,Star,Star>& AT,
  DistMatrixBase<T,Star,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View2x1");
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
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
elemental::DistMatrixBase<T,Star,Star>::LockedView2x1
( const DistMatrixBase<T,Star,Star>& AT,
  const DistMatrixBase<T,Star,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView2x1");
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
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
elemental::DistMatrixBase<T,Star,Star>::View2x2
( DistMatrixBase<T,Star,Star>& ATL, 
  DistMatrixBase<T,Star,Star>& ATR,
  DistMatrixBase<T,Star,Star>& ABL,
  DistMatrixBase<T,Star,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View2x2");
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
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
elemental::DistMatrixBase<T,Star,Star>::LockedView2x2
( const DistMatrixBase<T,Star,Star>& ATL, 
  const DistMatrixBase<T,Star,Star>& ATR,
  const DistMatrixBase<T,Star,Star>& ABL,
  const DistMatrixBase<T,Star,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView2x2");
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
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
elemental::DistMatrixBase<T,Star,Star>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    this->_height = height;
    this->_width = width;
    this->_localMatrix.ResizeTo( height, width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,Star,Star>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    T u = this->LocalEntry(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    this->LocalEntry(i,j) = u;
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., SetToIdentity and MakeTrapezoidal
//

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();

    if( shape == Lower )
    {
        for( int j=0; j<width; ++j )
        {
            int firstNonzero_i;
            if( side == Left )
                firstNonzero_i = max(j-offset,0);
            else
                firstNonzero_i = max(j-offset+height-width,0);

            const int boundary = min(height,firstNonzero_i);
            for( int i=0; i<boundary; ++i )
                this->LocalEntry(i,j) = (T)0;
        }
    }
    else
    {
        for( int j=0; j<width; ++j )
        {
            int firstZero_i;
            if( side == Left )
                firstZero_i = max(j-offset+1,0);
            else
                firstZero_i = max(j-offset+height-width+1,0);
            for( int i=firstZero_i; i<height; ++i )
                this->LocalEntry(i,j) = (T)0;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();

    this->SetToZero();
    for( int j=0; j<min(height,width); ++j )
        this->LocalEntry(j,j) = (T)1;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    // Create random matrix on process 0 and then broadcast
    const Grid& grid = this->GetGrid();
    const int height = this->Height();
    const int width = this->Width();
    const int bufSize = height*width;

    this->_auxMemory.Require( bufSize );

    T* buffer = this->_auxMemory.Buffer();
    if( grid.VCRank() == 0 )
    {
        for( int j=0; j<width; ++j )
            for( int i=0; i<height; ++i )
                buffer[i+j*height] = Random<T>();
    }
    Broadcast( buffer, bufSize, 0, grid.VCComm() );

    // Unpack
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            this->LocalEntry(i,j) = buffer[i+j*height];

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,MC,MR>& A )
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

    const Grid& grid = this->GetGrid();
    const int r = grid.Height();
    const int c = grid.Width(); 
    const int p = grid.Size();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeightOfA = A.LocalHeight();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalHeight = MaxLocalLength(height,r);
    const int maxLocalWidth = MaxLocalLength(width,c);

    const int portionSize = 
        max(maxLocalHeight*maxLocalWidth,MinCollectContrib);

    this->_auxMemory.Require( (p+1)*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, grid.VCComm() );

    // Unpack
    const int colAlignmentOfA = A.ColAlignment();
    const int rowAlignmentOfA = A.RowAlignment();
    for( int l=0; l<c; ++l )
    {
        const int rowShift = Shift( l, rowAlignmentOfA, c );
        const int localWidth = LocalLength( width, rowShift, c );

        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[(k+l*r)*portionSize];

            const int colShift = Shift( k, colAlignmentOfA, r );
            const int localHeight = LocalLength( height, colShift, r );

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(colShift+i*r,rowShift+j*c) = 
                        data[i+j*localHeight];
        }
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,MC,Star>& A )
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

    const Grid& grid = this->GetGrid();
    const int r = grid.Height();
    const int height = this->Height();
    const int width = this->Width();
    const int localHeightOfA = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength(height,r);

    const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

    this->_auxMemory.Require( (r+1)*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, grid.MCComm() );

    // Unpack
    const int colAlignmentOfA = A.ColAlignment();
    for( int k=0; k<r; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int colShift = Shift( k, colAlignmentOfA, r );
        const int localHeight = LocalLength( height, colShift, r );

        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(colShift+i*r,j) = data[i+j*localHeight];
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,MR>& A )
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

    const Grid& grid = this->GetGrid();
    const int c = grid.Width();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalWidth = MaxLocalLength(width,c);

    const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

    this->_auxMemory.Require( (c+1)*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<height; ++i )
            originalData[i+j*height] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, grid.MRComm() );

    // Unpack
    const int rowAlignmentOfA = A.RowAlignment();
    for( int k=0; k<c; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int rowShift = Shift( k, rowAlignmentOfA, c );
        const int localWidth = LocalLength( width, rowShift, c );

        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->LocalEntry(i,rowShift+j*c) = data[i+j*height];
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,MD,Star>& A )
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

    const Grid& grid = this->GetGrid();
    const int p = grid.Size();
    const int lcm = grid.LCM();
    const int ownerPath = grid.DiagPath( A.ColAlignment() );
    const int ownerPathRank = grid.DiagPathRank( A.ColAlignment() );

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength( height, lcm );
    const int portionSize = max( maxLocalHeight*width, MinCollectContrib );

    // Since a MD communicator has not been implemented, we will take
    // the suboptimal route of 'rounding up' everyone's contribution over the
    // VC communicator.
    this->_auxMemory.Require( (p+1)*portionSize );
    T* buffer = this->_auxMemory.Buffer();
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    if( A.InDiagonal() )
    {
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                sendBuf[i+j*localHeight] = A.LocalEntry(i,j);
    }

    // Communicate
    AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, grid.VCComm() );

    // Unpack
    for( int k=0; k<p; ++k )
    {
        if( grid.DiagPath( k ) == ownerPath )
        {
            const T* data = &recvBuf[k*portionSize];

            const int thisPathRank = grid.DiagPathRank( k );
            const int thisColShift = Shift( thisPathRank, ownerPathRank, lcm );
            const int thisLocalHeight = 
                LocalLength( height, thisColShift, lcm );

            for( int j=0; j<width; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    this->LocalEntry(thisColShift+i*lcm,j) = 
                        data[i+j*thisLocalHeight];
        }
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,MD>& A )
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

    const Grid& grid = this->GetGrid();
    const int p = grid.Size();
    const int lcm = grid.LCM();
    const int ownerPath = grid.DiagPath( A.RowAlignment() );
    const int ownerPathRank = grid.DiagPathRank( A.RowAlignment() );

    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = A.LocalWidth();
    const int maxLocalWidth = MaxLocalLength( width, lcm );
    const int portionSize = max( height*maxLocalWidth, MinCollectContrib );

    // Since a MD communicator has not been implemented, we will take
    // the suboptimal route of 'rounding up' everyone's contribution over the
    // VC communicator.
    this->_auxMemory.Require( (p+1)*portionSize );
    T* buffer = this->_auxMemory.Buffer();
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[portionSize];

    // Pack
    if( A.InDiagonal() )
    {
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                sendBuf[i+j*height] = A.LocalEntry(i,j);
    }

    // Communicate
    AllGather
    ( sendBuf, portionSize,
      recvBuf, portionSize, grid.VCComm() );

    // Unpack
    for( int k=0; k<p; ++k )
    {
        if( grid.DiagPath( k ) == ownerPath )
        {
            const T* data = &recvBuf[k*portionSize];

            const int thisPathRank = grid.DiagPathRank( k );
            const int thisRowShift = Shift( thisPathRank, ownerPathRank, lcm );
            const int thisLocalWidth =  LocalLength( width, thisRowShift, lcm );

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<height; ++i )
                    this->LocalEntry(i,thisRowShift+j*lcm) = data[i+j*height];
        }
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,MR,MC>& A )
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

    const Grid& grid = this->GetGrid();
    const int r = grid.Height();
    const int c = grid.Width();
    const int p = grid.Size();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeightOfA = A.LocalHeight();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalHeight = MaxLocalLength(height,c);
    const int maxLocalWidth = MaxLocalLength(width,r);

    const int portionSize = max(maxLocalHeight*maxLocalWidth,MinCollectContrib);

    this->_auxMemory.Require( (p+1)*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, grid.VRComm() );

    // Unpack
    const int colAlignmentOfA = A.ColAlignment();
    const int rowAlignmentOfA = A.RowAlignment();
    for( int l=0; l<r; ++l )
    {
        const int rowShift = Shift( l, rowAlignmentOfA, r );
        const int localWidth = LocalLength( width, rowShift, r );

        for( int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[(k+l*c)*portionSize];

            const int colShift = Shift( k, colAlignmentOfA, c );
            const int localHeight = LocalLength( height, colShift, c );

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(colShift+i*c,rowShift+j*r) =
                        data[i+j*localHeight];
        }
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,MR,Star>& A )
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

    const Grid& grid = this->GetGrid();
    const int c = grid.Width();
    const int height = this->Height();
    const int width = this->Width();
    const int localHeightOfA = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength(height,c);

    const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

    this->_auxMemory.Require( (c+1)*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, grid.MRComm() );

    // Unpack
    const int colAlignmentOfA = A.ColAlignment();
    for( int k=0; k<c; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int colShift = Shift( k, colAlignmentOfA, c );
        const int localHeight = LocalLength( height, colShift, c );

        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(colShift+i*c,j) = data[i+j*localHeight];
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,MC>& A )
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

    const Grid& grid = this->GetGrid();
    const int r = grid.Height();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalWidth = MaxLocalLength(width,r);

    const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

    this->_auxMemory.Require( (r+1)*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<height; ++i )
            originalData[i+j*height] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, grid.MCComm() );

    // Unpack
    const int rowAlignmentOfA = A.RowAlignment();
    for( int k=0; k<r; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int rowShift = Shift( k, rowAlignmentOfA, r );
        const int localWidth = LocalLength( width, rowShift, r );

        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->LocalEntry(i,rowShift+j*r) = data[i+j*height];
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,VC,Star>& A )
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

    const Grid& grid = this->GetGrid();
    const int p = grid.Size();
    const int height = this->Height();
    const int width = this->Width();
    const int localHeightOfA = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength(height,p);

    const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

    this->_auxMemory.Require( (p+1)*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, grid.VCComm() );

    // Unpack
    const int colAlignmentOfA = A.ColAlignment();
    for( int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int colShift = Shift( k, colAlignmentOfA, p );
        const int localHeight = LocalLength( height, colShift, p );

        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(colShift+i*p,j) = data[i+j*localHeight];
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,VC>& A )
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

    const Grid& grid = this->GetGrid();
    const int p = grid.Size();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalWidth = MaxLocalLength(width,p);

    const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

    this->_auxMemory.Require( (p+1)*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<height; ++i )
            originalData[i+j*height] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, grid.VCComm() );

    // Unpack
    const int rowAlignmentOfA = A.RowAlignment();
    for( int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int rowShift = Shift( k, rowAlignmentOfA, p );
        const int localWidth = LocalLength( width, rowShift, p );

        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->LocalEntry(i,rowShift+j*p) = data[i+j*height];
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,VR,Star>& A )
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

    const Grid& grid = this->GetGrid();
    const int p = grid.Size();
    const int height = this->Height();
    const int width = this->Width();
    const int localHeightOfA = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength(height,p);

    const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

    this->_auxMemory.Require( (p+1)*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, grid.VRComm() );

    // Unpack
    const int colAlignmentOfA = A.ColAlignment();
    for( int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int colShift = Shift( k, colAlignmentOfA, p );
        const int localHeight = LocalLength( height, colShift, p );

        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(colShift+i*p,j) = data[i+j*localHeight];
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,VR>& A )
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

    const Grid& grid = this->GetGrid();
    const int p = grid.Size();
    const int height = this->Height();
    const int width = this->Width();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalWidth = MaxLocalLength(width,p);

    const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

    this->_auxMemory.Require( (p+1)*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<height; ++i )
            originalData[i+j*height] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, grid.VRComm() );

    // Unpack
    const int rowAlignmentOfA = A.RowAlignment();
    for( int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int rowShift = Shift( k, rowAlignmentOfA, p );
        const int localWidth = LocalLength( width, rowShift, p );

        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->LocalEntry(i,rowShift+j*p) = data[i+j*height];
    }

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,Star>& A )
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
    this->_localMatrix = A.LockedLocalMatrix();
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
elemental::DistMatrix<R,Star,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int height = this->Height();
    const int width = this->Width();

    this->SetToRandom();
    for( int j=0; j<min(height,width); ++j )
        this->LocalEntry(j,j) += (R)this->Width();
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,Star,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int height = this->Height();
    const int width = this->Width();

    this->SetToRandom();
    for( int j=0; j<min(height,width); ++j )
        this->LocalEntry(j,j) = real(this->LocalEntry(j,j)) + (R)this->Width();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,Star>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    R u = real(this->LocalEntry(i,j));
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,Star>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    R u = imag(this->LocalEntry(i,j));
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,Star>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    real(this->LocalEntry(i,j)) = u;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,Star>::SetImag
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    imag(this->LocalEntry(i,j)) = u;
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrixBase<int,   Star,Star>;
template class elemental::DistMatrixBase<float, Star,Star>;
template class elemental::DistMatrixBase<double,Star,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,Star,Star>;
template class elemental::DistMatrixBase<dcomplex,Star,Star>;
#endif

template class elemental::DistMatrix<int,     Star,Star>;
template class elemental::DistMatrix<float,   Star,Star>;
template class elemental::DistMatrix<double,  Star,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,Star,Star>;
template class elemental::DistMatrix<dcomplex,Star,Star>;
#endif

