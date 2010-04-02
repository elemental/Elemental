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
Elemental::DistMatrix<T,Star,MC>::Print( const string msg ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::Print");
#endif
    if( _grid->VCRank() == 0 && msg != "" )
        cout << msg << endl;

    const int height     = Height();
    const int width      = Width();
    const int localWidth = LocalWidth();
    const int r          = _grid->Height();
    const int rowShift   = RowShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Only one process col needs to participate
    if( _grid->MRRank() == 0 )
    {
        T* sendBuf = new T[height*width];
        for( int i=0; i<height*width; ++i )
            sendBuf[i] = (T)0;
        for( int i=0; i<height; ++i )
            for( int j=0; j<localWidth; ++j )
                sendBuf[i+(rowShift+j*r)*height] = _localMatrix(i,j);

        // If we are the root, fill the receive buffer
        T* recvBuf = 0;
        if( _grid->MCRank() == 0 )
        {
            recvBuf = new T[height*width];
            for( int i=0; i<height*width; ++i )
                recvBuf[i] = (T)0;
        }

        // Sum the contributions and send to the root
        Reduce
        ( sendBuf, recvBuf, height*width, MPI_SUM, 0, _grid->MCComm() );
        delete[] sendBuf;

        if( _grid->MCRank() == 0 )
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
    Barrier( _grid->VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::AlignWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::AlignWith(DistMatrix[MR,MC])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _rowShift     = A.RowShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::AlignWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::AlignWith(DistMatrix[* ,MC])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _rowShift     = A.RowShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::AlignWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::AlignWith(DistMatrix[MC,MR])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _rowShift     = A.ColShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::AlignWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::AlignWith(DistMatrix[MC,* ])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _rowShift     = A.ColShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::AlignRowsWith
( const DistMatrix<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::AlignRowsWith
( const DistMatrix<T,MC,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::AlignRowsWith
( const DistMatrix<T,Star,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::AlignRowsWith
( const DistMatrix<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::ConformWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::ConformWith(DistMatrix[MC,MR])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_CONFORMING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _rowShift     = A.ColShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::ConformWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::ConformWith(DistMatrix[MC,* ])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_CONFORMING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _rowShift     = A.ColShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::ConformWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::ConformWith(DistMatrix[MR,MC])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_CONFORMING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _rowShift     = A.RowShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::ConformWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::ConformWith(DistMatrix[* ,MC])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_CONFORMING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _rowShift     = A.RowShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::FreeConstraints()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::FreeConstraints");
#endif
    _constrainedRowDist = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::View
( DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::View(A)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _rowAlignment = A.RowAlignment();
    _rowShift     = A.RowShift();
    _localMatrix.View( A.LocalMatrix() );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::LockedView
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE 
    PushCallStack("DistMatrix[* ,MC]::LockedView(A)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _rowAlignment = A.RowAlignment();
    _rowShift     = A.RowShift();
    _localMatrix.LockedView( A.LockedLocalMatrix() );
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::View
( DistMatrix<T,Star,MC>& A,
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::View(A,i,j,height,width)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int r   = _grid->Height();
        const int row = _grid->MCRank();

        _rowAlignment = (A.RowAlignment()+j) % r;
        _rowShift = Shift( row, _rowAlignment, r );

        const int localWidthBefore = LocalLength( j, A.RowShift(), r );
        const int localWidth = LocalLength( width, _rowShift, r );

        _localMatrix.View
        ( A.LocalMatrix(), i, localWidthBefore, height, localWidth );
    }
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::LockedView
( const DistMatrix<T,Star,MC>& A,
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::LockedView(A,i,j,height,width)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int r   = _grid->Height();
        const int row = _grid->MCRank();

        _rowAlignment = (A.RowAlignment()+j) % r;
        _rowShift = Shift( row, _rowAlignment, r );

        const int localWidthBefore = LocalLength( j, A.RowShift(), r );
        const int localWidth = LocalLength( width, _rowShift, r );

        _localMatrix.LockedView
        ( A.LockedLocalMatrix(), i, localWidthBefore, height, localWidth );
    }
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::View1x2
( DistMatrix<T,Star,MC>& AL,
  DistMatrix<T,Star,MC>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::View1x2");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _rowAlignment = AL.RowAlignment();
    _rowShift     = AL.RowShift();
    _localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::LockedView1x2
( const DistMatrix<T,Star,MC>& AL,
  const DistMatrix<T,Star,MC>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::LockedView1x2");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _rowAlignment = AL.RowAlignment();
    _rowShift     = AL.RowShift();
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
Elemental::DistMatrix<T,Star,MC>::View2x1
( DistMatrix<T,Star,MC>& AT,
  DistMatrix<T,Star,MC>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::View2x1");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
    if( AT.RowAlignment() != AB.RowAlignment() )
        throw "Cannot combine misaligned 2x1 array of matrices.";
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _rowAlignment = AT.RowAlignment();
    _rowShift     = AT.RowShift();
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
Elemental::DistMatrix<T,Star,MC>::LockedView2x1
( const DistMatrix<T,Star,MC>& AT,
  const DistMatrix<T,Star,MC>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::LockedView2x1");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
    if( AT.RowAlignment() != AB.RowAlignment() )
        throw "Cannot combine misaligned 2x1 array of matrices.";
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _rowAlignment = AT.RowAlignment();
    _rowShift     = AT.RowShift();
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
Elemental::DistMatrix<T,Star,MC>::View2x2
( DistMatrix<T,Star,MC>& ATL,
  DistMatrix<T,Star,MC>& ATR,
  DistMatrix<T,Star,MC>& ABL,
  DistMatrix<T,Star,MC>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::View2x2");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
    if( ATL.RowAlignment() != ABL.RowAlignment() ||
        ATR.RowAlignment() != ABR.RowAlignment()    )
    {
        throw "Cannot combine misaligned 2x2 grid of matrices.";
    }
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _rowAlignment = ATL.RowAlignment();
    _rowShift     = ATL.RowShift();
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
Elemental::DistMatrix<T,Star,MC>::LockedView2x2
( const DistMatrix<T,Star,MC>& ATL,
  const DistMatrix<T,Star,MC>& ATR,
  const DistMatrix<T,Star,MC>& ABL,
  const DistMatrix<T,Star,MC>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::LockedView2x2");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
    if( ATL.RowAlignment() != ABL.RowAlignment() ||
        ATR.RowAlignment() != ABR.RowAlignment()    )
    {
        throw "Cannot combine misaligned 2x2 grid of matrices.";
    }
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _rowAlignment = ATL.RowAlignment();
    _rowShift     = ATL.RowShift();
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
Elemental::DistMatrix<T,Star,MC>::ResizeTo
( const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::ResizeTo");
    CHECK_IF_LOCKED_VIEW;
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _height = height;
    _width  = width;
    _localMatrix.ResizeTo(height,LocalLength(width,_rowShift,_grid->Height()));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
Elemental::DistMatrix<T,Star,MC>::Get
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::Get");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        throw msg.str();
    }
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that
    // row within each process column
    const int ownerRow = (j + RowAlignment()) % _grid->Height();

    T u;
    if( _grid->MCRank() == ownerRow )
    {
        const int jLoc = (j-RowShift()) / _grid->Height();
        u = _localMatrix(i,jLoc);
    }
    Broadcast( &u, 1, ownerRow, _grid->MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::Set
( const int i, const int j, const T u )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::Set");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        throw msg.str();
    }
#endif
    const int ownerRow = (j + RowAlignment()) % _grid->Height();

    if( _grid->MCRank() == ownerRow )
    {
        const int jLoc = (j-RowShift()) / _grid->Height();
        _localMatrix(i,jLoc) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Utility functions, e.g., SetToIdentity and MakeTrapezoidal
//----------------------------------------------------------------------------//

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::MakeTrapezoidal
( const Side side, const Shape shape, const int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::MakeTrapezoidal");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int width = Width();
    const int localWidth = LocalWidth();
    const int r = _grid->Height();
    const int rowShift = RowShift();

    if( shape == Lower )
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*r;
            int firstNonzero_i;
            if( side == Left )
                firstNonzero_i = max(j-offset,0);
            else
                firstNonzero_i = max(j-offset+height-width,0);

            const int boundary = min(height,firstNonzero_i);
            for( int i=0; i<boundary; ++i )
                _localMatrix(i,jLoc) = (T)0;
        }
    }
    else
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*r;
            int firstZero_i;
            if( side == Left )
                firstZero_i = max(j-offset+1,0);
            else
                firstZero_i = max(j-offset+height-width+1,0);
            for( int i=firstZero_i; i<height; ++i )
                _localMatrix(i,jLoc) = (T)0;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::SetToIdentity");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int localWidth = LocalWidth();
    const int r = _grid->Height();
    const int rowShift = RowShift();

    _localMatrix.SetToZero();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*r;
        if( j < height )
            _localMatrix(j,jLoc) = (T)1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::SetToRandom");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height     = Height();    
    const int localWidth = LocalWidth();
    const int bufSize    = height*localWidth;

    _auxMemory.Require( bufSize );

    // Create a random matrix on process column 0, then broadcast
    T* buffer = _auxMemory.Buffer();
    if( _grid->MRRank() == 0 )
    {
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                buffer[i+j*height] = Random<T>();
    }
    Broadcast( buffer, bufSize, 0, _grid->MRComm() );

    // Unpack
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<height; ++i )
            _localMatrix(i,j) = buffer[i+j*height];

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::SetToRandomDiagDominant()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::SetToRandomDiagDominant");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int localWidth = LocalWidth();
    const int r = _grid->Height();
    const int rowShift = RowShift();

    SetToRandom();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*r;
        if( j < height )
            _localMatrix(j,jLoc) += (T)max(Height(),Width());
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::SetToZero");
    CHECK_IF_LOCKED_VIEW;
#endif
    _localMatrix.SetToZero();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MC>::AllReduce()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::AllReduce");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    const int localSize = max( localHeight*localWidth, MinCollectContrib );

    _auxMemory.Require( 2*localSize );
    T* buffer = _auxMemory.Buffer();
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[localSize];

    // Pack
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            sendBuf[i+j*localHeight] = _localMatrix(i,j);

    // AllReduce sum
    wrappers::MPI::AllReduce
    ( sendBuf, recvBuf, localSize, MPI_SUM, _grid->MRComm() );

    // Unpack
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = recvBuf[i+j*localHeight];

    _auxMemory.Release();

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[MC,MR]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
        ( new DistMatrix<T,Star,VR>(*_grid) );
    *A_Star_VR = A;

    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
        ( new DistMatrix<T,Star,VC>(true,RowAlignment(),*_grid) );
    *A_Star_VC = *A_Star_VR;
    delete A_Star_VR.release(); // lowers memory highwater

    *this = *A_Star_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[MC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
        ( new DistMatrix<T,MC,MR>(*_grid) );
    *A_MC_MR   = A;

    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
        ( new DistMatrix<T,Star,VR>(*_grid) );
    *A_Star_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
        ( new DistMatrix<T,Star,VC>(true,RowAlignment(),*_grid) );
    *A_Star_VC = *A_Star_VR;
    delete A_Star_VR.release(); // lowers memory highwater

    *this = *A_Star_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[* ,MR]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( A.Height() == 1 )
    {
        if( !_viewing )
            ResizeTo( 1, A.Width() );

        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int myRow = _grid->MCRank();
        const int rankCM = _grid->VCRank();
        const int rankRM = _grid->VRRank();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int rowShiftOfA = A.RowShift();

        const int width = Width();
        const int maxLocalVectorWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxLocalVectorWidth,MinCollectContrib);

        const int rowShiftVC = Shift(rankCM,rowAlignment,p);
        const int rowShiftVROfA = Shift(rankRM,rowAlignmentOfA,p);
        const int sendRankCM = (rankCM+(p+rowShiftVROfA-rowShiftVC)) % p;
        const int recvRankRM = (rankRM+(p+rowShiftVC-rowShiftVROfA)) % p;
        const int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        _auxMemory.Require( (c+1)*portionSize );
        T* buffer = _auxMemory.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        // A[* ,VR] <- A[* ,MR]
        {
            const int shift = Shift(rankRM,rowAlignmentOfA,p);
            const int offset = (shift-rowShiftOfA) / c;
            const int thisLocalWidth = LocalLength(width,shift,p);

            for( int j=0; j<thisLocalWidth; ++j )
                sendBuf[j] = A.LocalEntry(0,offset+j*r);
        }

        // A[* ,VC] <- A[* ,VR]
        SendRecv
        ( sendBuf, portionSize, sendRankCM, 0,
          recvBuf, portionSize, recvRankCM, MPI_ANY_TAG, _grid->VCComm() );

        // A[* ,MC] <- A[* ,VC]
        AllGather
        ( recvBuf, portionSize,
          sendBuf, portionSize, _grid->MRComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &sendBuf[k*portionSize];

            const int shift = Shift(myRow+r*k,rowAlignment,p);
            const int offset = (shift-RowShift()) / r;
            const int thisLocalWidth = LocalLength(width,shift,p);

            for( int j=0; j<thisLocalWidth; ++j )
                _localMatrix(0,offset+j*c) = data[j];
        }

        _auxMemory.Release();
    }
    else
    {
        auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
            ( new DistMatrix<T,Star,VR>(*_grid) );
        *A_Star_VR = A;

        auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
            ( new DistMatrix<T,Star,VC>(true,RowAlignment(),*_grid) );
        *A_Star_VC = *A_Star_VR;
        delete A_Star_VR.release(); // lowers memory highwater

        *this = *A_Star_VC;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[MD,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    REPORT_UNIMPLEMENTED_FEATURE;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[* ,MD]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A ); 
#endif
    REPORT_UNIMPLEMENTED_FEATURE;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[MR,MC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        if( ! ConstrainedRowDist() )
        {
            _rowAlignment = A.RowAlignment();
            _rowShift = Shift
                        ( _grid->MCRank(), _rowAlignment, _grid->Height() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( RowAlignment() == A.RowAlignment() )
    {
        const int c = _grid->Width();
        const int height = Height();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(height,c);

        const int portionSize = 
            max(maxLocalHeightOfA*localWidth,MinCollectContrib);

        _auxMemory.Require( (c+1)*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
        for( int j=0; j<localWidth; ++j )
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

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(colShift+i*c,j) = data[i+j*localHeight];
        }

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [* ,MC] <- [MR,MC]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int row = _grid->MCRank();

        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const int height = Height();
        const int width = Width();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeightOfA = MaxLocalLength(height,c);
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            max(maxLocalHeightOfA*maxLocalWidth,MinCollectContrib);

        _auxMemory.Require( (c+1)*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                secondBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);

        // Perform the SendRecv: puts the new data into the first buffer
        SendRecv
        ( secondBuffer, portionSize, sendRow, 0,
          firstBuffer,  portionSize, recvRow, MPI_ANY_TAG, 
          _grid->MCComm() );

        // Use the output of the SendRecv as input to the AllGather
        AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, _grid->MRComm() );

        // Unpack the contents of each member of the process row
        const int colAlignmentOfA = A.ColAlignment();
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int colShift = Shift( k, colAlignmentOfA, c );
            const int localHeight = LocalLength( height, colShift, c );
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(colShift+i*c,j) = data[i+j*localHeight];
        }

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[MR,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,MR,MC> A_MR_MC(*_grid);

    A_MR_MC = A;
    *this   = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[* ,MC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        if( ! ConstrainedRowDist() )
        {
            _rowAlignment = A.RowAlignment();
            _rowShift = A.RowShift();
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( RowAlignment() == A.RowAlignment() )
    {
        _localMatrix = A.LockedLocalMatrix();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [* ,MC] <- [* ,MC]." << endl;
#endif
        const int rank = _grid->MCRank();
        const int r = _grid->Height();

        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRank = (rank+r+rowAlignmentOfA-rowAlignment) % r;

        const int height = Height();
        const int localWidth = LocalWidth();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = height * localWidthOfA;
        const int recvSize = height * localWidth;

        _auxMemory.Require( sendSize + recvSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<height; ++i )
                sendBuffer[i+j*height] = A.LocalEntry(i,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, _grid->MCComm() );

        // Unpack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                _localMatrix(i,j) = recvBuffer[i+j*height];

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[VC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
        ( new DistMatrix<T,VR,Star>(*_grid) );
    *A_VR_Star = A;

    auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
        ( new DistMatrix<T,MR,MC>(false,0,true,RowAlignment(),*_grid) );
    *A_MR_MC = *A_VR_Star;
    delete A_VR_Star.release();

    *this = *A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[* ,VC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        if( ! ConstrainedRowDist() )
        {
            _rowAlignment = A.RowAlignment() % _grid->Height();
            _rowShift = Shift
                        ( _grid->MCRank(), _rowAlignment, _grid->Height() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( RowAlignment() == A.RowAlignment() % _grid->Height() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int row = _grid->MCRank();

        const int height = Height();
        const int width = Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidthOfA = MaxLocalLength(width,p);

        const int portionSize = 
            max(height*maxLocalWidthOfA,MinCollectContrib);

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
        const int rowShift = RowShift();
        const int rowAlignmentOfA = A.RowAlignment();
        for( int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShiftOfA = Shift( row+k*r, rowAlignmentOfA, p );
            const int rowOffset = (rowShiftOfA-rowShift) / r;
            const int localWidth = LocalLength( width, rowShiftOfA, p );

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    _localMatrix(i,rowOffset+j*c) = data[i+j*height];
        }

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [* ,MC] <- [* ,VC]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int row = _grid->MCRank();
        const int rank = _grid->VCRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int rowShift = RowShift();

        const int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const int height = Height();
        const int width = Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidthOfA = MaxLocalLength(width,p);

        const int portionSize = 
            max(height*maxLocalWidthOfA,MinCollectContrib);

        _auxMemory.Require( (c+1)*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<height; ++i )
                secondBuffer[i+j*height] = A.LocalEntry(i,j);

        // Perform the SendRecv: puts the new data into the first buffer
        SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, 0, _grid->VCComm() );

        // Use the SendRecv as input to the AllGather
        AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, _grid->MRComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int rowShiftOfA = Shift(row+r*k,rowAlignment,p);
            const int rowOffset = (rowShiftOfA-rowShift) / r;
            const int localWidth = LocalLength( width, rowShiftOfA, p );
            
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    _localMatrix(i,rowOffset+j*c) = data[i+j*height];
        }

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[VR,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,MR,MC> A_MR_MC(*_grid);

    A_MR_MC = A;
    *this   = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[* ,VR]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,Star,VC> A_Star_VC(*_grid);

    A_Star_VC = A;
    *this     = A_Star_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MC>&
Elemental::DistMatrix<T,Star,MC>::operator=
( const DistMatrix<T,Star,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC] = DistMatrix[* ,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int r = _grid->Height();
    const int rowShift = RowShift();

    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = A.LocalEntry(i,rowShift+j*r);
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class Elemental::DistMatrix<int,     Star,MC>;
template class Elemental::DistMatrix<float,   Star,MC>;
template class Elemental::DistMatrix<double,  Star,MC>;
#ifndef WITHOUT_COMPLEX
template class Elemental::DistMatrix<scomplex,Star,MC>;
template class Elemental::DistMatrix<dcomplex,Star,MC>;
#endif

