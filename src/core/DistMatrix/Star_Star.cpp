/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "Elemental/DistMatrix.hpp"
#include "./DistMatrixMacros.hpp"
using namespace std;
using namespace Elemental;
using namespace Elemental::utilities;
using namespace Elemental::wrappers::MPI;

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::Print( const string& msg ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::Print");
#endif
    if( _grid->VCRank() == 0 && msg != "" )
        cout << msg << endl;

    const int height = Height();
    const int width  = Width();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( _grid->VCRank() == 0 )
    {
        for( int i=0; i<height; ++i )
        {
            for( int j=0; j<width; ++j )
                cout << _localMatrix(i,j) << " ";
            cout << endl;
        }
        cout << endl;
    }
    Barrier( _grid->VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::View
( DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::View(A)");
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _localMatrix.View( A.LocalMatrix() );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::LockedView
( const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::LockedView(A)");
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _localMatrix.LockedView( A.LockedLocalMatrix() );
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::View
( DistMatrix<T,Star,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::View(A,i,j,height,width)");
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    _localMatrix.View( A.LocalMatrix(), i, j, height, width );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::LockedView
( const DistMatrix<T,Star,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::LockedView(A,i,j,height,width)");
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    _localMatrix.LockedView( A.LockedLocalMatrix(), i, j, height, width );
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::View1x2
( DistMatrix<T,Star,Star>& AL,
  DistMatrix<T,Star,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::View1x2");
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::LockedView1x2
( const DistMatrix<T,Star,Star>& AL,
  const DistMatrix<T,Star,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::LockedView1x2");
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _localMatrix.LockedView1x2( AL.LockedLocalMatrix(), 
                                AR.LockedLocalMatrix() );
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::View2x1
( DistMatrix<T,Star,Star>& AT,
  DistMatrix<T,Star,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::View2x1");
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, 
                             AB );
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _localMatrix.View2x1( AT.LocalMatrix(),
                          AB.LocalMatrix() );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::LockedView2x1
( const DistMatrix<T,Star,Star>& AT,
  const DistMatrix<T,Star,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::LockedView2x1");
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, 
                             AB );
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _localMatrix.LockedView2x1( AT.LockedLocalMatrix(),
                                AB.LockedLocalMatrix() );
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::View2x2
( DistMatrix<T,Star,Star>& ATL, 
  DistMatrix<T,Star,Star>& ATR,
  DistMatrix<T,Star,Star>& ABL,
  DistMatrix<T,Star,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::View2x2");
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _localMatrix.View2x2( ATL.LocalMatrix(), ATR.LocalMatrix(),
                          ABL.LocalMatrix(), ABR.LocalMatrix() );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::LockedView2x2
( const DistMatrix<T,Star,Star>& ATL, 
  const DistMatrix<T,Star,Star>& ATR,
  const DistMatrix<T,Star,Star>& ABL,
  const DistMatrix<T,Star,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::LockedView2x2");
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _localMatrix.LockedView2x2
    ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
      ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::ResizeTo");
    CHECK_IF_LOCKED_VIEW;
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _height = height;
    _width  = width;
    _localMatrix.ResizeTo( height, width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
Elemental::DistMatrix<T,Star,Star>::Get
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::Get");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    T u = _localMatrix(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::Set");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    _localMatrix(i,j) = u;
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Utility functions, e.g., SetToIdentity and MakeTrapezoidal
//----------------------------------------------------------------------------//

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::MakeTrapezoidal");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int width = Width();

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
                _localMatrix(i,j) = (T)0;
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
                _localMatrix(i,j) = (T)0;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::SetToIdentity");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int width = Width();

    _localMatrix.SetToZero();
    for( int j=0; j<min(height,width); ++j )
        _localMatrix(j,j) = (T)1;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::SetToRandom");
    CHECK_IF_LOCKED_VIEW;
#endif
    // Create random matrix on process 0 and then broadcast
    const int height = Height();
    const int width = Width();
    const int bufSize = height*width;

    _auxMemory.Require( bufSize );

    T* buffer = _auxMemory.Buffer();
    if( _grid->VCRank() == 0 )
    {
        for( int j=0; j<width; ++j )
            for( int i=0; i<height; ++i )
                buffer[i+j*height] = Random<T>();
    }
    Broadcast( buffer, bufSize, 0, _grid->VCComm() );

    // Unpack
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            _localMatrix(i,j) = buffer[i+j*height];

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::SetToRandomDiagDominant()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::SetToRandomDiagDominant");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int width = Width();

    SetToRandom();
    for( int j=0; j<min(height,width); ++j )
        _localMatrix(j,j) += (T)max(Height(),Width());
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,Star>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::SetToZero");
    CHECK_IF_LOCKED_VIEW;
#endif
    _localMatrix.SetToZero();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix<MC,MR>");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int r = _grid->Height();
    const int c = _grid->Width(); 
    const int p = _grid->Size();

    const int height = Height();
    const int width = Width();
    const int localHeightOfA = A.LocalHeight();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalHeight = MaxLocalLength(height,r);
    const int maxLocalWidth = MaxLocalLength(width,c);

    const int portionSize = 
        max(maxLocalHeight*maxLocalWidth,MinCollectContrib);

    _auxMemory.Require( (p+1)*portionSize );

    T* buffer = _auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, _grid->VCComm() );

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
                    _localMatrix(colShift+i*r,rowShift+j*c) = 
                          data[i+j*localHeight];
        }
    }

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix<MC,*>");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int r = _grid->Height();
    const int height = Height();
    const int width = Width();
    const int localHeightOfA = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength(height,r);

    const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

    _auxMemory.Require( (r+1)*portionSize );

    T* buffer = _auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, _grid->MCComm() );

    // Unpack
    const int colAlignmentOfA = A.ColAlignment();
    for( int k=0; k<r; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int colShift = Shift( k, colAlignmentOfA, r );
        const int localHeight = LocalLength( height, colShift, r );

        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(colShift+i*r,j) = data[i+j*localHeight];
    }

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix[* ,MR]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int c = _grid->Width();
    const int height = Height();
    const int width = Width();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalWidth = MaxLocalLength(width,c);

    const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

    _auxMemory.Require( (c+1)*portionSize );

    T* buffer = _auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<height; ++i )
            originalData[i+j*height] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, _grid->MRComm() );

    // Unpack
    const int rowAlignmentOfA = A.RowAlignment();
    for( int k=0; k<c; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int rowShift = Shift( k, rowAlignmentOfA, c );
        const int localWidth = LocalLength( width, rowShift, c );

        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                _localMatrix(i,rowShift+j*c) = data[i+j*height];
    }

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix<MD,* >");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int p = _grid->Size();
    const int lcm = _grid->LCM();
    const int ownerPath = _grid->DiagPath( A.ColAlignment() );
    const int ownerPathRank = _grid->DiagPathRank( A.ColAlignment() );

    const int height = Height();
    const int width = Width();
    const int localHeight = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength( height, lcm );
    const int portionSize = max( maxLocalHeight*width, MinCollectContrib );

    // Since a MD communicator has not been implemented, we will take
    // the suboptimal route of 'rounding up' everyone's contribution over the
    // VC communicator.
    _auxMemory.Require( (p+1)*portionSize );
    T* buffer = _auxMemory.Buffer();
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
    AllGather( sendBuf, portionSize,
               recvBuf, portionSize, _grid->VCComm() );

    // Unpack
    for( int k=0; k<p; ++k )
    {
        if( _grid->DiagPath( k ) == ownerPath )
        {
            const T* data = &recvBuf[k*portionSize];

            const int thisPathRank = _grid->DiagPathRank( k );
            const int thisColShift = Shift( thisPathRank, ownerPathRank, lcm );
            const int thisLocalHeight = 
                LocalLength( height, thisColShift, lcm );

            for( int j=0; j<width; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    _localMatrix(thisColShift+i*lcm,j) = 
                        data[i+j*thisLocalHeight];
        }
    }

    _auxMemory.Release();

#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix[* ,MD]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A ); 
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int p = _grid->Size();
    const int lcm = _grid->LCM();
    const int ownerPath = _grid->DiagPath( A.RowAlignment() );
    const int ownerPathRank = _grid->DiagPathRank( A.RowAlignment() );

    const int height = Height();
    const int width = Width();
    const int localWidth = A.LocalWidth();
    const int maxLocalWidth = MaxLocalLength( width, lcm );
    const int portionSize = max( height*maxLocalWidth, MinCollectContrib );

    // Since a MD communicator has not been implemented, we will take
    // the suboptimal route of 'rounding up' everyone's contribution over the
    // VC communicator.
    _auxMemory.Require( (p+1)*portionSize );
    T* buffer = _auxMemory.Buffer();
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
    AllGather( sendBuf, portionSize,
               recvBuf, portionSize, _grid->VCComm() );

    // Unpack
    for( int k=0; k<p; ++k )
    {
        if( _grid->DiagPath( k ) == ownerPath )
        {
            const T* data = &recvBuf[k*portionSize];

            const int thisPathRank = _grid->DiagPathRank( k );
            const int thisRowShift = Shift( thisPathRank, ownerPathRank, lcm );
            const int thisLocalWidth =  LocalLength( width, thisRowShift, lcm );

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<height; ++i )
                    _localMatrix(i,thisRowShift+j*lcm) = data[i+j*height];
        }
    }

    _auxMemory.Release();

#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix<MR,MC>");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int r = _grid->Height();
    const int c = _grid->Width();
    const int p = _grid->Size();

    const int height = Height();
    const int width = Width();
    const int localHeightOfA = A.LocalHeight();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalHeight = MaxLocalLength(height,c);
    const int maxLocalWidth = MaxLocalLength(width,r);

    const int portionSize = max(maxLocalHeight*maxLocalWidth,MinCollectContrib);

    _auxMemory.Require( (p+1)*portionSize );

    T* buffer = _auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, _grid->VRComm() );

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
                    _localMatrix(colShift+i*c,rowShift+j*r) =
                          data[i+j*localHeight];
        }
    }

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix<MR,*>");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int c = _grid->Width();
    const int height = Height();
    const int width = Width();
    const int localHeightOfA = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength(height,c);

    const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

    _auxMemory.Require( (c+1)*portionSize );

    T* buffer = _auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, _grid->MRComm() );

    // Unpack
    const int colAlignmentOfA = A.ColAlignment();
    for( int k=0; k<c; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int colShift = Shift( k, colAlignmentOfA, c );
        const int localHeight = LocalLength( height, colShift, c );

        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(colShift+i*c,j) = data[i+j*localHeight];
    }

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix[* ,MC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int r = _grid->Height();
    const int height = Height();
    const int width = Width();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalWidth = MaxLocalLength(width,r);

    const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

    _auxMemory.Require( (r+1)*portionSize );

    T* buffer = _auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<height; ++i )
            originalData[i+j*height] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, _grid->MCComm() );

    // Unpack
    const int rowAlignmentOfA = A.RowAlignment();
    for( int k=0; k<r; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int rowShift = Shift( k, rowAlignmentOfA, r );
        const int localWidth = LocalLength( width, rowShift, r );

        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                _localMatrix(i,rowShift+j*r) = data[i+j*height];
    }

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix<VC,*>");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int p = _grid->Size();
    const int height = Height();
    const int width = Width();
    const int localHeightOfA = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength(height,p);

    const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

    _auxMemory.Require( (p+1)*portionSize );

    T* buffer = _auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, _grid->VCComm() );

    // Unpack
    const int colAlignmentOfA = A.ColAlignment();
    for( int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int colShift = Shift( k, colAlignmentOfA, p );
        const int localHeight = LocalLength( height, colShift, p );

        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(colShift+i*p,j) = data[i+j*localHeight];
    }

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix[* ,VC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int p = _grid->Size();
    const int height = Height();
    const int width = Width();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalWidth = MaxLocalLength(width,p);

    const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

    _auxMemory.Require( (p+1)*portionSize );

    T* buffer = _auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<height; ++i )
            originalData[i+j*height] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, _grid->VCComm() );

    // Unpack
    const int rowAlignmentOfA = A.RowAlignment();
    for( int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int rowShift = Shift( k, rowAlignmentOfA, p );
        const int localWidth = LocalLength( width, rowShift, p );

        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                _localMatrix(i,rowShift+j*p) = data[i+j*height];
    }

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix<VR,*>");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int p = _grid->Size();
    const int height = Height();
    const int width = Width();
    const int localHeightOfA = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength(height,p);

    const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

    _auxMemory.Require( (p+1)*portionSize );

    T* buffer = _auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, _grid->VRComm() );

    // Unpack
    const int colAlignmentOfA = A.ColAlignment();
    for( int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int colShift = Shift( k, colAlignmentOfA, p );
        const int localHeight = LocalLength( height, colShift, p );

        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(colShift+i*p,j) = data[i+j*localHeight];
    }

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix[* ,VR]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int p = _grid->Size();
    const int height = Height();
    const int width = Width();
    const int localWidthOfA = A.LocalWidth();
    const int maxLocalWidth = MaxLocalLength(width,p);

    const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

    _auxMemory.Require( (p+1)*portionSize );

    T* buffer = _auxMemory.Buffer();
    T* originalData = &buffer[0];
    T* gatheredData = &buffer[portionSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<height; ++i )
            originalData[i+j*height] = A.LocalEntry(i,j);

    // Communicate
    AllGather
    ( originalData, portionSize,
      gatheredData, portionSize, _grid->VRComm() );

    // Unpack
    const int rowAlignmentOfA = A.RowAlignment();
    for( int k=0; k<p; ++k )
    {
        const T* data = &gatheredData[k*portionSize];

        const int rowShift = Shift( k, rowAlignmentOfA, p );
        const int localWidth = LocalLength( width, rowShift, p );

        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                _localMatrix(i,rowShift+j*p) = data[i+j*height];
    }

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,Star>&
Elemental::DistMatrix<T,Star,Star>::operator=
( const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ] = DistMatrix[* ,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );
    _localMatrix = A.LockedLocalMatrix();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class Elemental::DistMatrix<int,     Star,Star>;
template class Elemental::DistMatrix<float,   Star,Star>;
template class Elemental::DistMatrix<double,  Star,Star>;
#ifndef WITHOUT_COMPLEX
template class Elemental::DistMatrix<scomplex,Star,Star>;
template class Elemental::DistMatrix<dcomplex,Star,Star>;
#endif

