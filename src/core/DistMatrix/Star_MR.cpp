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
Elemental::DistMatrix<T,Star,MR>::Print( const string msg ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::Print");
#endif
    if( _grid->VCRank() == 0 && msg != "" )
        cout << msg << endl;

    const int height     = Height();
    const int width      = Width();
    const int localWidth = LocalWidth();
    const int c          = _grid->Width();
    const int rowShift   = RowShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Only one process row needs to participate
    if( _grid->MCRank() == 0 )
    {
        T* sendBuf = new T[height*width];
        for( int i=0; i<height*width; ++i )
            sendBuf[i] = (T)0;
        for( int i=0; i<height; ++i )
            for( int j=0; j<localWidth; ++j )
                sendBuf[i+(rowShift+j*c)*height] = _localMatrix(i,j);

        // If we are the root, fill the receive buffer
        T* recvBuf = 0;
        if( _grid->MRRank() == 0 )
        {
            recvBuf = new T[height*width];
            for( int i=0; i<height*width; ++i )
                recvBuf[i] = (T)0;
        }

        // Sum the contributions and send to the root
        Reduce
        ( sendBuf, recvBuf, height*width, MPI_SUM, 0, _grid->MRComm() );
        delete[] sendBuf;

        if( _grid->MRRank() == 0 )
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
Elemental::DistMatrix<T,Star,MR>::AlignWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::AlignWith(DistMatrix[MC,MR])");
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
Elemental::DistMatrix<T,Star,MR>::AlignWith
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::AlignWith(DistMatrix[* ,MR])");
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
Elemental::DistMatrix<T,Star,MR>::AlignWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::AlignWith(DistMatrix[MR,MC])");
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
Elemental::DistMatrix<T,Star,MR>::AlignWith
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::AlignWith(DistMatrix[MR,* ])");
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
Elemental::DistMatrix<T,Star,MR>::AlignWith
( const DistMatrix<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::AlignWith(DistMatrix[VR,* ])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment() % _grid->Width();
    _rowShift     = Shift( _grid->MRRank(), _rowAlignment, _grid->Width() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::AlignWith
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::AlignWith(DistMatrix[* ,VR])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment() % _grid->Width();
    _rowShift     = Shift( _grid->MRRank(), _rowAlignment, _grid->Width() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::AlignRowsWith
( const DistMatrix<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::AlignRowsWith
( const DistMatrix<T,Star,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::AlignRowsWith
( const DistMatrix<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::AlignRowsWith
( const DistMatrix<T,MR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::AlignRowsWith
( const DistMatrix<T,VR,Star>& A )
{ AlignWith( A ); } 

template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::AlignRowsWith
( const DistMatrix<T,Star,VR>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::ConformWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::ConformWith(DistMatrix[MC,MR])");
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
Elemental::DistMatrix<T,Star,MR>::ConformWith
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::ConformWith(DistMatrix[* ,MR])");
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
Elemental::DistMatrix<T,Star,MR>::ConformWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::ConformWith(DistMatrix[MR,MC])");
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
Elemental::DistMatrix<T,Star,MR>::ConformWith
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::ConformWith(DistMatrix[MR,* ])");
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
Elemental::DistMatrix<T,Star,MR>::FreeConstraints()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::FreeConstraints");
#endif
    _constrainedRowDist = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::View
( DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::View(A)");
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
Elemental::DistMatrix<T,Star,MR>::LockedView
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[*, MR]::LockedView(A)");
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
Elemental::DistMatrix<T,Star,MR>::View
( DistMatrix<T,Star,MR>& A,
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::View(A,i,j,height,width)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int c   = _grid->Width();
        const int col = _grid->MRRank();

        _rowAlignment = (A.RowAlignment()+j) % c;
        _rowShift = Shift( col, _rowAlignment, c );

        const int localWidthBefore = LocalLength( j, A.RowShift(), c );
        const int localWidth = LocalLength( width, _rowShift, c );

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
Elemental::DistMatrix<T,Star,MR>::LockedView
( const DistMatrix<T,Star,MR>& A,
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::LockedView(A,i,j,height,width)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int c   = _grid->Width();
        const int col = _grid->MRRank();

        _rowAlignment = (A.RowAlignment()+j) % c;
        _rowShift = Shift( col, _rowAlignment, c );

        const int localWidthBefore = LocalLength( j, A.RowShift(), c );
        const int localWidth = LocalLength( width, _rowShift, c );

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
Elemental::DistMatrix<T,Star,MR>::View1x2
( DistMatrix<T,Star,MR>& AL,
  DistMatrix<T,Star,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::View1x2");
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
Elemental::DistMatrix<T,Star,MR>::LockedView1x2
( const DistMatrix<T,Star,MR>& AL,
  const DistMatrix<T,Star,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::LockedView1x2");
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
Elemental::DistMatrix<T,Star,MR>::View2x1
( DistMatrix<T,Star,MR>& AT,
  DistMatrix<T,Star,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::View2x1");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
    if( AT.RowAlignment() != AB.RowAlignment() )
        throw "2x1 misaligned, cannot combine.";
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
Elemental::DistMatrix<T,Star,MR>::LockedView2x1
( const DistMatrix<T,Star,MR>& AT,
  const DistMatrix<T,Star,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::LockedView2x1");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
    if( AT.RowAlignment() != AB.RowAlignment() )
        throw "2x1 misaligned, cannot combine.";
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
Elemental::DistMatrix<T,Star,MR>::View2x2
( DistMatrix<T,Star,MR>& ATL,
  DistMatrix<T,Star,MR>& ATR,
  DistMatrix<T,Star,MR>& ABL,
  DistMatrix<T,Star,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::View2x2");
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
        throw "2x2 is misaligned, cannot combine.";
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
Elemental::DistMatrix<T,Star,MR>::LockedView2x2
( const DistMatrix<T,Star,MR>& ATL,
  const DistMatrix<T,Star,MR>& ATR,
  const DistMatrix<T,Star,MR>& ABL,
  const DistMatrix<T,Star,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::LockedView2x2");
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
        throw "2x2 is misaligned, cannot combine.";
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
Elemental::DistMatrix<T,Star,MR>::ResizeTo
( const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::ResizeTo");
    CHECK_IF_LOCKED_VIEW;
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _height = height;
    _width  = width;
    _localMatrix.ResizeTo(height,LocalLength(width,_rowShift,_grid->Width()));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
Elemental::DistMatrix<T,Star,MR>::Get
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
    // We will determine the owner column of entry (i,j) and broadcast from that
    // column within each process row
    const int ownerCol = (j + RowAlignment()) % _grid->Width();

    T u;
    if( _grid->MRRank() == ownerCol )
    {
        const int jLoc = (j-RowShift()) / _grid->Width();
        u = _localMatrix(i,jLoc);
    }
    Broadcast( &u, 1, ownerCol, _grid->MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::Set
( const int i, const int j, const T u )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::Set");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        throw msg.str();
    }
#endif
    const int ownerCol = (j + RowAlignment()) % _grid->Width();

    if( _grid->MRRank() == ownerCol )
    {
        const int jLoc = (j-RowShift()) / _grid->Width();
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
Elemental::DistMatrix<T,Star,MR>::MakeTrapezoidal
( const Side side, const Shape shape, const int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::MakeTrapezoidal");
    CHECK_IF_LOCKED_VIEW;
#endif

    const int height = Height();
    const int width = Width();
    const int localWidth = LocalWidth();
    const int c = _grid->Width();
    const int rowShift = RowShift();

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
            for( int i=0; i<boundary; ++i )
                _localMatrix(i,jLoc) = (T)0;
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
Elemental::DistMatrix<T,Star,MR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::SetToIdentity");
    CHECK_IF_LOCKED_VIEW;
#endif

    const int height = Height();
    const int localWidth = LocalWidth();
    const int c = _grid->Width();
    const int rowShift = RowShift();

    _localMatrix.SetToZero();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*c;
        if( j < height )
            _localMatrix(j,jLoc) = (T)1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::SetToRandom");
    CHECK_IF_LOCKED_VIEW;
#endif

    const int height     = Height();
    const int localWidth = LocalWidth();
    const int bufSize    = height*localWidth;

    _auxMemory.Require( bufSize );

    // Create random matrix on process row 0, then broadcast
    T* buffer = _auxMemory.Buffer();
    if( _grid->MCRank() == 0 )
    {
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                buffer[i+j*height] = Random<T>();
    }
    Broadcast( buffer, bufSize, 0, _grid->MCComm() );

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
Elemental::DistMatrix<T,Star,MR>::SetToRandomDiagDominant()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::SetToRandomDiagDominant");
    CHECK_IF_LOCKED_VIEW;
#endif

    const int height     = Height();
    const int localWidth = LocalWidth();
    const int c          = _grid->Width();
    const int rowShift   = RowShift();

    SetToRandom();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*c;
        if( j < height )
            _localMatrix(j,jLoc) += (T)max(Height(),Width());
    }
#ifndef RELEASE
    PopCallStack();
#endif
}


template<typename T>
void
Elemental::DistMatrix<T,Star,MR>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::SetToZero");
    CHECK_IF_LOCKED_VIEW;
#endif
    _localMatrix.SetToZero();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> 
void
Elemental::DistMatrix<T,Star,MR>::AllReduce()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::AllReduce");
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

    // AllReduce col
    wrappers::MPI::AllReduce
    ( sendBuf, recvBuf, localSize, MPI_SUM, _grid->MCComm() );

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
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[MC,MR]");
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
                        ( _grid->MRRank(), _rowAlignment, _grid->Width() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( RowAlignment() == A.RowAlignment() )
    {
        if( A.Height() == 1 )
        {
            const int localWidth = LocalWidth();

            _auxMemory.Require( localWidth );
            T* bcastBuf = _auxMemory.Buffer();

            if( _grid->MCRank() == A.ColAlignment() )
            {
                _localMatrix = A.LockedLocalMatrix();

                // Pack
                for( int j=0; j<localWidth; ++j )
                    bcastBuf[j] = _localMatrix(0,j);
            }

            // Communicate
            Broadcast
            ( bcastBuf, localWidth, A.ColAlignment(), _grid->MCComm() );

            // Unpack
            for( int j=0; j<localWidth; ++j )
                _localMatrix(0,j) = bcastBuf[j];

            _auxMemory.Release();
        }
        else
        {
            const int r = _grid->Height();
            const int height = Height();
            const int localWidth = LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int portionSize = 
                max(maxLocalHeight*localWidth,MinCollectContrib);

            _auxMemory.Require( (r+1)*portionSize );

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
              gatheredData, portionSize, _grid->MCComm() );

            // Unpack
            const int colAlignmentOfA = A.ColAlignment();
            for( int k=0; k<r; ++k )
            {
                const T* data = &gatheredData[k*portionSize];

                const int colShift = Shift( k, colAlignmentOfA, r );
                const int localHeight = LocalLength( height, colShift, r );

                for( int j=0; j<localWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        _localMatrix(colShift+i*r,j) = data[i+j*localHeight];
            }

            _auxMemory.Release();
        }
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [* ,MR] <- [MC,MR]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int col = _grid->MRRank();

        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

        if( A.Height() == 1 )
        {
            const int localWidth = LocalWidth();
            T* bcastBuf;

            if( _grid->MCRank() == A.ColAlignment() )
            {
                const int localWidthOfA = A.LocalWidth();

                _auxMemory.Require( localWidth+localWidthOfA );
                T* buffer = _auxMemory.Buffer();
                T* sendBuf = &buffer[0];
                bcastBuf   = &buffer[localWidthOfA];

                // Pack
                for( int j=0; j<localWidthOfA; ++j )
                    sendBuf[j] = A.LocalEntry(0,j);

                // Communicate
                SendRecv
                ( sendBuf,  localWidthOfA, sendCol, 0,
                  bcastBuf, localWidth,    recvCol, MPI_ANY_TAG,
                  _grid->MRComm()                          );
            }
            else
            {
                _auxMemory.Require( localWidth );
                bcastBuf = _auxMemory.Buffer();
            }

            // Communicate
            Broadcast
            ( bcastBuf, localWidth, A.ColAlignment(), _grid->MCComm() );

            // Unpack
            for( int j=0; j<localWidth; ++j )
                _localMatrix(0,j) = bcastBuf[j];

            _auxMemory.Release();
        }
        else
        {
            const int height = Height();
            const int localWidth  = LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int localWidthOfA  = A.LocalWidth();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int portionSize = 
                max(maxLocalHeight*localWidth,MinCollectContrib);

            _auxMemory.Require( (r+1)*portionSize );

            T* buffer = _auxMemory.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[portionSize];

            // Pack
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    secondBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);

            // Perform the SendRecv: puts the new data into the first buffer
            SendRecv
            ( secondBuffer, portionSize, sendCol, 0,
              firstBuffer,  portionSize, recvCol, MPI_ANY_TAG, 
              _grid->MRComm()                            );

            // Use the output of the SendRecv as input to the AllGather
            AllGather
            ( firstBuffer,  portionSize,
              secondBuffer, portionSize, _grid->MCComm() );

            // Unpack the contents of each member of the process col
            const int colAlignmentOfA = A.ColAlignment();
            for( int k=0; k<r; ++k )
            {
                const T* data = &secondBuffer[k*portionSize];

                const int colShift = Shift( k, colAlignmentOfA, r );
                const int localHeight = LocalLength( height, colShift, r );
                for( int j=0; j<localWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        _localMatrix(colShift+i*r,j) = data[i+j*localHeight];
            }

            _auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[MC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,MC,MR> 
        A_MC_MR(false,0,true,RowAlignment(),*_grid);

    A_MC_MR = A;
    *this   = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[* ,MR]");
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
            cout << "Unaligned [* ,MR] <- [* ,MR]." << endl;
#endif
        const int rank = _grid->MRRank();
        const int c = _grid->Width();

        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvRank = (rank+c+rowAlignmentOfA-rowAlignment) % c;

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
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, _grid->MRComm() );

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
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[MD,* ]");
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
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[* ,MD]");
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
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[MR,MC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
        ( new DistMatrix<T,Star,VC>(*_grid) );
    *A_Star_VC = A;

    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
        ( new DistMatrix<T,Star,VR>(true,RowAlignment(),*_grid) );
    *A_Star_VR = *A_Star_VC;
    delete A_Star_VC.release(); // lowers memory highwater

    *this = *A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[MR,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
        ( new DistMatrix<T,VR,Star>(*_grid) );
    *A_VR_Star = A;

    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
        ( new DistMatrix<T,VC,Star>(*_grid) );
    *A_VC_Star = *A_VR_Star;
    delete A_VR_Star.release(); // lowers memory highwater

    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
        ( new DistMatrix<T,MC,MR>(false,0,true,RowAlignment(),*_grid) );
    *A_MC_MR = *A_VC_Star;
    delete A_VC_Star.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[* ,MC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
        ( new DistMatrix<T,Star,VC>(*_grid) );
    *A_Star_VC = A;

    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
        ( new DistMatrix<T,Star,VR>(true,RowAlignment(),*_grid) );
    *A_Star_VR = *A_Star_VC;
    delete A_Star_VC.release(); // lowers memory highwater

    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
        ( new DistMatrix<T,MC,MR>(*_grid) );
    *A_MC_MR = *A_Star_VR;
    delete A_Star_VR.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[VC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,MC,MR> A_MC_MR(false,0,true,RowAlignment(),*_grid);

    A_MC_MR = A;
    *this   = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[* ,VC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,Star,VR> A_Star_VR(true,RowAlignment(),*_grid);

    A_Star_VR = A;
    *this     = A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[VR,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
        ( new DistMatrix<T,VC,Star>(*_grid) );
    *A_VC_Star = A;

    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
        ( new DistMatrix<T,MC,MR>(false,0,true,RowAlignment(),*_grid) );
    *A_MC_MR = *A_VC_Star;
    delete A_VC_Star.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[* ,VR]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif

    if( !_viewing )
    {
        if( ! ConstrainedRowDist() )
        {
            _rowAlignment = A.RowAlignment() % _grid->Width();
            _rowShift = Shift
                        ( _grid->MRRank(), _rowAlignment, _grid->Width() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( RowAlignment() == A.RowAlignment() % _grid->Width() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int col = _grid->MRRank();

        const int width = Width();
        const int height = Height();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidthOfA = MaxLocalLength(width,p);

        const int portionSize = max(height*maxLocalWidthOfA,MinCollectContrib);

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
        const int rowShift = RowShift();
        const int rowAlignmentOfA = A.RowAlignment();
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShiftOfA = Shift( col+k*c, rowAlignmentOfA, p );
            const int rowOffset = (rowShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, rowShiftOfA, p );

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    _localMatrix(i,rowOffset+j*r) = data[i+j*height];
        }

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [* ,MR] <- [* ,VR]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int col = _grid->MRRank();
        const int rank = _grid->VRRank();

        // Perform the SendRecv to make A have the same rowAlignment
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int rowShift = RowShift();

        const int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const int width = Width();
        const int height = Height();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidthOfA = MaxLocalLength(width,p);

        const int portionSize = max(height*maxLocalWidthOfA,MinCollectContrib);

        _auxMemory.Require( (r+1)*portionSize );

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
          firstBuffer,  portionSize, recvRank, MPI_ANY_TAG, _grid->VRComm() );

        // Use the SendRecv as input to the AllGather
        AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, _grid->MCComm() );

        // Unpack
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int rowShiftOfA = Shift( col+c*k, rowAlignment, p );
            const int rowOffset = (rowShiftOfA-rowShift) / c;
            const int localWidth = LocalLength( width, rowShiftOfA, p );

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    _localMatrix(i,rowOffset+j*r) = data[i+j*height];
        }

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MR>&
Elemental::DistMatrix<T,Star,MR>::operator=
( const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR] = DistMatrix[* ,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int c = _grid->Width();
    const int rowShift = RowShift();

    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = A.LocalEntry(i,rowShift+j*c);
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class Elemental::DistMatrix<int,     Star,MR>;
template class Elemental::DistMatrix<float,   Star,MR>;
template class Elemental::DistMatrix<double,  Star,MR>;
#ifndef WITHOUT_COMPLEX
template class Elemental::DistMatrix<scomplex,Star,MR>;
template class Elemental::DistMatrix<dcomplex,Star,MR>;
#endif

