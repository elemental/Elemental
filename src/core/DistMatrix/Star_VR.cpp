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
#include "ElementalDistMatrix.h"
#include "./DistMatrixMacros.h"
using namespace std;
using namespace Elemental;
using namespace Elemental::utilities;
using namespace Elemental::wrappers::MPI;

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::Print( const string msg ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::Print");
#endif
    if( _grid->VRRank() == 0 && msg != "" )
        cout << msg << endl;

    const int height     = Height();
    const int width      = Width();
    const int localWidth = LocalWidth();
    const int p          = _grid->Size();
    const int rowShift   = RowShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    T* sendBuf = new T[height*width];
    for( int i=0; i<height*width; ++i )
        sendBuf[i] = (T)0;
    for( int i=0; i<height; ++i )
        for( int j=0; j<localWidth; ++j )
            sendBuf[i+(rowShift+j*p)*height] = _localMatrix(i,j);

    // If we are the root, fill the receive buffer
    T* recvBuf = 0;
    if( _grid->VRRank() == 0 )
    {
        recvBuf = new T[height*width];
        for( int i=0; i<height*width; ++i )
            recvBuf[i] = (T)0;
    }

    // Sum the contributions and send to the root
    Reduce( sendBuf, recvBuf, height*width, MPI_SUM, 0, _grid->VRComm() );
    delete[] sendBuf;

    if( _grid->VRRank() == 0 )
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
    Barrier( _grid->VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::AlignWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::AlignWith(DistMatrix[MC,MR])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _rowShift = Shift( _grid->VRRank(), _rowAlignment, _grid->Size() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::AlignWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::AlignWith(DistMatrix[MR,MC])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _rowShift = Shift( _grid->VRRank(), _rowAlignment, _grid->Size() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::AlignWith
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::AlignWith(Distmatrix[MR,* ])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _rowShift = Shift( _grid->VRRank(), _rowAlignment, _grid->Size() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::AlignWith
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::AlignWith(Distmatrix[* ,MR])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _rowShift = Shift( _grid->VRRank(), _rowAlignment, _grid->Size() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::AlignWith
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::AlignWith(DistMatrix[* ,VR])");
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
Elemental::DistMatrix<T,Star,VR>::AlignWith
( const DistMatrix<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::AlignWith(DistMatrix[VR,* ])");
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
Elemental::DistMatrix<T,Star,VR>::AlignRowsWith
( const DistMatrix<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::AlignRowsWith
( const DistMatrix<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::AlignRowsWith
( const DistMatrix<T,MR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::AlignRowsWith
( const DistMatrix<T,Star,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::AlignRowsWith
( const DistMatrix<T,Star,VR>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::AlignRowsWith
( const DistMatrix<T,VR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::ConformWith
( const DistMatrix<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::ConformWith(DistMatrix[VR,* ])");
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
Elemental::DistMatrix<T,Star,VR>::ConformWith
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::ConformWith(DistMatrix[* ,VR])");
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
Elemental::DistMatrix<T,Star,VR>::FreeConstraints()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::FreeConstraints");
    CHECK_IF_NONEXISTANT_CONSTRAINTS;
#endif
    _constrainedRowDist = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::View
( DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::View(A)");
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
Elemental::DistMatrix<T,Star,VR>::LockedView
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::LockedView(A)");
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
Elemental::DistMatrix<T,Star,VR>::View
( DistMatrix<T,Star,VR>& A,
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::View(A,i,j,height,width)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int rowMajorRank = _grid->VRRank();
        const int size = _grid->Size();

        _rowAlignment = (A.RowAlignment()+j) % size;
        _rowShift = Shift( rowMajorRank, _rowAlignment, size );

        const int localWidthBefore = LocalLength(j,A.RowShift(),size);
        const int localWidth = LocalLength(width,_rowShift,size);

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
Elemental::DistMatrix<T,Star,VR>::LockedView
( const DistMatrix<T,Star,VR>& A,
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::LockedView(A,i,j,height,width)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int rowMajorRank = _grid->VRRank();
        const int size = _grid->Size();

        _rowAlignment = (A.RowAlignment()+j) % size;
        _rowShift = Shift( rowMajorRank, _rowAlignment, size );

        const int localWidthBefore = LocalLength(j,A.RowShift(),size);
        const int localWidth = LocalLength(width,_rowShift,size);

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
Elemental::DistMatrix<T,Star,VR>::View1x2
( DistMatrix<T,Star,VR>& AL,
  DistMatrix<T,Star,VR>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::View1x2");
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
Elemental::DistMatrix<T,Star,VR>::LockedView1x2
( const DistMatrix<T,Star,VR>& AL,
  const DistMatrix<T,Star,VR>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::LockedView1x2");
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
Elemental::DistMatrix<T,Star,VR>::View2x1
( DistMatrix<T,Star,VR>& AT,
  DistMatrix<T,Star,VR>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::View2x1");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
    if( AT.RowAlignment() != AB.RowAlignment() )
    {
        if( _grid->VCRank() == 0 )
            cerr << "2x1 is misaligned, cannot combine." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _rowAlignment = AT.RowAlignment();
    _rowShift     = AT.RowShift();
    _localMatrix.View2x1( AT.LocalMatrix(), AB.LocalMatrix() );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::LockedView2x1
( const DistMatrix<T,Star,VR>& AT,
  const DistMatrix<T,Star,VR>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::LockedView2x1");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
    if( AT.RowAlignment() != AB.RowAlignment() )
    {
        if( _grid->VCRank() == 0 )
            cerr << "2x1 is misaligned, cannot combine." << endl;
        DumpCallStack();
        throw exception();
    }
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
Elemental::DistMatrix<T,Star,VR>::View2x2
( DistMatrix<T,Star,VR>& ATL,
  DistMatrix<T,Star,VR>& ATR,
  DistMatrix<T,Star,VR>& ABL,
  DistMatrix<T,Star,VR>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::View2x2");
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
        if( _grid->VCRank() == 0 )
            cerr << "2x2 is misaligned, cannot combine." << endl;
        DumpCallStack();
        throw exception();
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
Elemental::DistMatrix<T,Star,VR>::LockedView2x2
( const DistMatrix<T,Star,VR>& ATL,
  const DistMatrix<T,Star,VR>& ATR,
  const DistMatrix<T,Star,VR>& ABL,
  const DistMatrix<T,Star,VR>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::LockedView2x2");
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
        if( _grid->VCRank() == 0 )
            cerr << "2x2 is misaligned, cannot combine." << endl;
        DumpCallStack();
        throw exception();
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
Elemental::DistMatrix<T,Star,VR>::ResizeTo
( const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::ResizeTo");
    CHECK_IF_LOCKED_VIEW;
    if( height < 0 || width < 0 )
    {
        if( _grid->VCRank() == 0 )
            cerr << "Height and width must be non-negative." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height = height;
    _width  = width;
    _localMatrix.ResizeTo(height,LocalLength(width,_rowShift,_grid->Size()));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
Elemental::DistMatrix<T,Star,VR>::Get
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::Get");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        if( _grid->VRRank() == 0 )
        {
            cerr << "Entry (" << i << "," << j << ") is out of bounds of "
                 << Height() << " x " << Width() << " matrix." << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire grid
    const int ownerRank = (j + RowAlignment()) % _grid->Size();

    T u;
    if( _grid->VRRank() == ownerRank )
    {
        const int jLoc = (j-RowShift()) / _grid->Size();
        u = _localMatrix(i,jLoc);
    }
    Broadcast( &u, 1, ownerRank, _grid->VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::Set
( const int i, const int j, const T u )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::Set");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        if( _grid->VRRank() == 0 )
        {
            cerr << "Entry (" << i << "," << j << ") is out of bounds of "
                 << Height() << " x " << Width() << " matrix." << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    const int ownerRank = (j + RowAlignment()) % _grid->Size();

    if( _grid->VRRank() == ownerRank )
    {
        const int jLoc = (j-RowShift()) / _grid->Size();
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
Elemental::DistMatrix<T,Star,VR>::MakeTrapezoidal
( const Side side, const Shape shape, const int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::MakeTrapezoidal");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int width = Width();
    const int localWidth = LocalWidth();
    const int p = _grid->Size();
    const int rowShift = RowShift();

    if( shape == Lower )
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*p;
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
            const int j = rowShift + jLoc*p;
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
Elemental::DistMatrix<T,Star,VR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::SetToIdentity");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int localWidth = LocalWidth();
    const int p = _grid->Size();
    const int rowShift = RowShift();

    _localMatrix.SetToZero();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*p;
        if( j < height )
            _localMatrix(j,jLoc) = (T)1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::SetToRandom");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int localWidth = LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<height; ++i )
            _localMatrix(i,j) = Random<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::SetToRandomDiagDominant()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::SetToRandomDiagDominant");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int localWidth = LocalWidth();
    const int p = _grid->Size();
    const int rowShift = RowShift();

    SetToRandom();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*p;
        if( j < height )
            _localMatrix(j,jLoc) += (T)max(Height(),Width());
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VR>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::SetToZero");
    CHECK_IF_LOCKED_VIEW;
#endif
    _localMatrix.SetToZero();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[MC,MR]");
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
                        ( _grid->VRRank(), _rowAlignment, _grid->Size() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( RowAlignment() % _grid->Width() == A.RowAlignment() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int col = _grid->MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = Height();
        const int width = Width();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,r);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        _auxMemory.Require( 2*r*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[r*portionSize];

        // Pack
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const int thisRank = col+k*c;
            const int thisRowShift = Shift(thisRank,rowAlignment,p);
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = 
                          A.LocalEntry(i,thisRowOffset+j*r);
        }

        // Communicate
        AllToAll( sendBuffer, portionSize,
                  recvBuffer, portionSize, _grid->MCComm() );

        // Unpack
        for( int k=0; k<r; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignmentOfA,r);
            const int thisLocalHeight = LocalLength(height,thisColShift,r);

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    _localMatrix(thisColShift+i*r,j) = 
                          data[i+j*thisLocalHeight];
        }

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [* ,VR] <- [MC,MR]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int col = _grid->MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendCol = (col+c+(rowAlignment%c)-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-(rowAlignment%c)) % c;

        const int height = Height();
        const int width = Width();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,r);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        _auxMemory.Require( 2*r*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[r*portionSize];

        // Pack
        for( int k=0; k<r; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const int thisRank = sendCol+k*c;
            const int thisRowShift = Shift(thisRank,rowAlignment,p);
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / c;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = 
                          A.LocalEntry(i,thisRowOffset+j*r);
        }

        // AllToAll to gather all of the unaligned [*,VR] data into firstBuffer
        AllToAll
        ( secondBuffer, portionSize,
          firstBuffer,  portionSize, _grid->MCComm() );

        // SendRecv: properly align the [*,VR] via a trade in the column
        SendRecv
        ( firstBuffer,  portionSize, sendCol, 0,
          secondBuffer, portionSize, recvCol, MPI_ANY_TAG, 
          _grid->MRComm() );

        // Unpack
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignmentOfA,r);
            const int thisLocalHeight = LocalLength(height,thisColShift,r);

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    _localMatrix(thisColShift+i*r,j) = 
                          data[i+j*thisLocalHeight];
        }

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[MC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,MC,MR> A_MC_MR(*_grid);

    A_MC_MR = A;
    *this   = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[* ,MR]");
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
                        ( _grid->VRRank(), _rowAlignment, _grid->Size() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( RowAlignment() % _grid->Width() == A.RowAlignment() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int rowShift = RowShift();
        const int rowShiftOfA = A.RowShift();
        const int rowOffset = (rowShift-rowShiftOfA) / c;

        const int height = Height();
        const int localWidth = LocalWidth();

        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                _localMatrix(i,j) = A.LocalEntry(i,rowOffset+j*r);
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [* ,VR] <- [* ,MR]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int row = _grid->MCRank();
        const int col = _grid->MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        // We will SendRecv A[*,VR] within our process row to fix alignments.
        const int sendCol = (col+c+(rowAlignment%c)-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-(rowAlignment%c)) % c;
        const int sendRank = sendCol + c*row;

        const int sendRowShift = Shift( sendRank, rowAlignment, p );
        const int sendRowOffset = (sendRowShift-rowShiftOfA) / c;

        const int height = Height();
        const int width = Width();
        const int localWidth = LocalWidth();
        const int localWidthOfSend = LocalLength(width,sendRowShift,p);

        const int sendSize = height * localWidthOfSend;
        const int recvSize = height * localWidth;

        _auxMemory.Require( sendSize + recvSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidthOfSend; ++j )
            for( int i=0; i<height; ++i )
                sendBuffer[i+j*height] = A.LocalEntry(i,sendRowOffset+j*r);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, MPI_ANY_TAG, _grid->MRComm() );

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
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[MD,* ]");
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
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[* ,MD]");
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
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[MR,MC]");
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
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[MR,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
        ( new DistMatrix<T,MR,MC>(*_grid) );
    *A_MR_MC = A;

    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
        ( new DistMatrix<T,Star,VC>(*_grid) );
    *A_Star_VC = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_Star_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[* ,MC]");
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
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[VC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,MC,MR> A_MC_MR(*_grid);

    A_MC_MR = A;
    *this   = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[* ,VC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );
    
    const int height = Height();
    const int localWidth = LocalWidth();
    const int localWidthOfA = A.LocalWidth();

    const int sendSize = height * localWidthOfA;
    const int recvSize = height * localWidth;

    const int r = _grid->Height();
    const int c = _grid->Width();
    const int p = _grid->Size();
    const int rankCM = _grid->VCRank();
    const int rankRM = _grid->VRRank();

    const int rowShift = RowShift();
    const int rowShiftOfA = A.RowShift();

    // Compute which rowmajor rank has the rowShift equal to our rowShiftOfA
    const int sendRankRM = (rankRM+(p+rowShiftOfA-rowShift)) % p;

    // Compute which rowmajor rank has the A rowShift that we need
    const int recvRankCM = (rankCM+(p+rowShift-rowShiftOfA)) % p;
    const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

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
    ( sendBuffer, sendSize, sendRankRM, 0,
      recvBuffer, recvSize, recvRankRM, MPI_ANY_TAG, _grid->VRComm() );

    // Unpack
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<height; ++i )
            _localMatrix(i,j) = recvBuffer[i+j*height];

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[VR,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
        ( new DistMatrix<T,MR,MC>(*_grid) );
    *A_MR_MC = A;

    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
        ( new DistMatrix<T,Star,VC>(*_grid) );
    *A_Star_VC = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_Star_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[* ,VR]");
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
            cout << "Unaligned [* ,VR] <- [* ,VR]." << endl;
#endif
        const int rank = _grid->VRRank();
        const int p = _grid->Size();

        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

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
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, _grid->VRComm() );

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
const DistMatrix<T,Star,VR>&
Elemental::DistMatrix<T,Star,VR>::operator=
( const DistMatrix<T,Star,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR] = DistMatrix[* ,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int p = _grid->Size();
    const int rowShift = RowShift();

    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = A.LocalEntry(i,rowShift+j*p);
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class Elemental::DistMatrix<int,     Star,VR>;
template class Elemental::DistMatrix<float,   Star,VR>;
template class Elemental::DistMatrix<double,  Star,VR>;
#ifndef WITHOUT_COMPLEX
template class Elemental::DistMatrix<scomplex,Star,VR>;
template class Elemental::DistMatrix<dcomplex,Star,VR>;
#endif

