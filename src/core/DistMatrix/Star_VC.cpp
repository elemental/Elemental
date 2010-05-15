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
Elemental::DistMatrix<T,Star,VC>::Print( const string& msg ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::Print");
#endif
    if( _grid->VCRank() == 0 && msg != "" )
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
    if( _grid->VCRank() == 0 )
    {
        recvBuf = new T[height*width];
        for( int i=0; i<height*width; ++i )
            recvBuf[i] = (T)0;
    }

    // Sum the contributions and send to the root
    Reduce( sendBuf, recvBuf, height*width, MPI_SUM, 0, _grid->VCComm() );
    delete[] sendBuf;

    if( _grid->VCRank() == 0 )
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
    Barrier( _grid->VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::AlignWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::AlignWith(DistMatrix[MC,MR])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _rowShift = Shift( _grid->VCRank(), _rowAlignment, _grid->Size() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::AlignWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::AlignWith(DistMatrix[MR,MC])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _rowShift = Shift( _grid->VCRank(), _rowAlignment, _grid->Size() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::AlignWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::AlignWith(DistMatrix[MC,* ])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _rowShift = Shift( _grid->VCRank(), _rowAlignment, _grid->Size() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::AlignWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::AlignWith(DistMatrix[* ,MC])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _rowShift = Shift( _grid->VCRank(), _rowAlignment, _grid->Size() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::AlignWith
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::AlignWith(DistMatrix[* ,VC])");
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
Elemental::DistMatrix<T,Star,VC>::AlignWith
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::AlignWith(DistMatrix[VC,* ])");
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
Elemental::DistMatrix<T,Star,VC>::AlignRowsWith
( const DistMatrix<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::AlignRowsWith
( const DistMatrix<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::AlignRowsWith
( const DistMatrix<T,MC,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::AlignRowsWith
( const DistMatrix<T,Star,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::AlignRowsWith
( const DistMatrix<T,Star,VC>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::AlignRowsWith
( const DistMatrix<T,VC,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::ConformWith
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::ConformWith(DistMatrix[VC,* ])");
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
Elemental::DistMatrix<T,Star,VC>::ConformWith
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::ConformWith(DistMatrix[* ,VC])");
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
Elemental::DistMatrix<T,Star,VC>::FreeConstraints()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::FreeConstraints");
#endif
    _constrainedRowDist = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::View
( DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::View(A)");
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
Elemental::DistMatrix<T,Star,VC>::LockedView
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::LockedView(A)");
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
Elemental::DistMatrix<T,Star,VC>::View
( DistMatrix<T,Star,VC>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::View(A,i,j,height,width)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int colMajorRank = _grid->VCRank();
        const int size = _grid->Size();

        _rowAlignment = (A.RowAlignment()+j) % size;
        _rowShift = Shift( colMajorRank, _rowAlignment, size );

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
Elemental::DistMatrix<T,Star,VC>::LockedView
( const DistMatrix<T,Star,VC>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::LockedView(A,i,j,height,width)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int colMajorRank = _grid->VCRank();
        const int size = _grid->Size();

        _rowAlignment = (A.RowAlignment()+j) % size;
        _rowShift = Shift( colMajorRank, _rowAlignment, size );

        const int localWidth = LocalLength(width,_rowShift,size);

        _localMatrix.LockedView
        ( A.LockedLocalMatrix(), i, A.RowShift(), height, localWidth );
    }
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::View1x2
( DistMatrix<T,Star,VC>& AL,
  DistMatrix<T,Star,VC>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::View1x2");
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
Elemental::DistMatrix<T,Star,VC>::LockedView1x2
( const DistMatrix<T,Star,VC>& AL,
  const DistMatrix<T,Star,VC>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::LockedView1x2");
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
Elemental::DistMatrix<T,Star,VC>::View2x1
( DistMatrix<T,Star,VC>& AT,
  DistMatrix<T,Star,VC>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::View2x1");
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
    _localMatrix.View2x1( AT.LocalMatrix(), AB.LocalMatrix() );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::LockedView2x1
( const DistMatrix<T,Star,VC>& AT,
  const DistMatrix<T,Star,VC>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::LockedView2x1");
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
Elemental::DistMatrix<T,Star,VC>::View2x2
( DistMatrix<T,Star,VC>& ATL, 
  DistMatrix<T,Star,VC>& ATR,
  DistMatrix<T,Star,VC>& ABL,
  DistMatrix<T,Star,VC>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::View2x2");
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
Elemental::DistMatrix<T,Star,VC>::LockedView2x2
( const DistMatrix<T,Star,VC>& ATL, 
  const DistMatrix<T,Star,VC>& ATR,
  const DistMatrix<T,Star,VC>& ABL,
  const DistMatrix<T,Star,VC>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::LockedView2x2");
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
Elemental::DistMatrix<T,Star,VC>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::ResizeTo");
    CHECK_IF_LOCKED_VIEW;
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
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
Elemental::DistMatrix<T,Star,VC>::Get
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::Get");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire grid
    const int ownerRank = (j + RowAlignment()) % _grid->Size();

    T u;
    if( _grid->VCRank() == ownerRank )
    {
        const int jLoc = (j-RowShift()) / _grid->Size();
        u = _localMatrix(i,jLoc);
    }
    Broadcast( &u, 1, ownerRank, _grid->VCComm() );
    
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
Elemental::DistMatrix<T,Star,VC>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::Set");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const int ownerRank = (j + RowAlignment()) % _grid->Size();

    if( _grid->VCRank() == ownerRank )
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
Elemental::DistMatrix<T,Star,VC>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::MakeTrapezoidal");
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
Elemental::DistMatrix<T,Star,VC>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::SetToIdentity");
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
Elemental::DistMatrix<T,Star,VC>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::SetToRandom");
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
Elemental::DistMatrix<T,Star,VC>::SetToRandomDiagDominant()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::SetToRandomDiagDominant");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height     = Height();
    const int localWidth = LocalWidth();
    const int p          = _grid->Size();
    const int rowShift   = RowShift();

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
Elemental::DistMatrix<T,Star,VC>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::SetToZero");
    CHECK_IF_LOCKED_VIEW;
#endif
    _localMatrix.SetToZero();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[MC,MR]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,Star,VR> A_Star_VR(*_grid);

    A_Star_VR = A;
    *this     = A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[MC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
        ( new DistMatrix<T,MC,MR>(*_grid) );
    *A_MC_MR = A;

    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
        ( new DistMatrix<T,Star,VR>(*_grid) );
    *A_Star_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[* ,MR]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,Star,VR> A_Star_VR(*_grid);

    A_Star_VR = A;
    *this     = A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[MD,* ]");
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
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[* ,MD]");
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
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[MR,MC]");
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
                        ( _grid->VCRank(), _rowAlignment, _grid->Size() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( RowAlignment() % _grid->Height() == A.RowAlignment() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int row = _grid->MCRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = Height();
        const int width = Width();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,c);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        _auxMemory.Require( 2*c*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[c*portionSize];

        // Pack
        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const int thisRank = row+k*r;
            const int thisRowShift = Shift(thisRank,rowAlignment,p);
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = 
                          A.LocalEntry(i,thisRowOffset+j*c);
        }

        // Communicate
        AllToAll( sendBuffer, portionSize,
                  recvBuffer, portionSize, _grid->MRComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignmentOfA,c);
            const int thisLocalHeight = LocalLength(height,thisColShift,c);

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    _localMatrix(thisColShift+i*c,j) = 
                          data[i+j*thisLocalHeight];
        }

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [* ,VC] <- [MR,MC]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int row = _grid->MCRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+(rowAlignment%r)-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-(rowAlignment%r)) % r;

        const int height = Height();
        const int width = Width();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,c);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        _auxMemory.Require( 2*c*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[c*portionSize];

        // Pack
        for( int k=0; k<c; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const int thisRank = sendRow+k*r;
            const int thisRowShift = Shift(thisRank,rowAlignment,p);
            const int thisRowOffset = (thisRowShift-rowShiftOfA) / r;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = 
                          A.LocalEntry(i,thisRowOffset+j*c);
        }

        // AllToAll to gather all of the unaligned [*,VC] data into firstBuffer
        AllToAll
        ( secondBuffer, portionSize,
          firstBuffer,  portionSize, _grid->MRComm() );

        // SendRecv: properly align the [*,VC] via a trade in the column
        SendRecv
        ( firstBuffer,  portionSize, sendRow, 0,
          secondBuffer, portionSize, recvRow, MPI_ANY_TAG, 
          _grid->MCComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignmentOfA,c);
            const int thisLocalHeight = LocalLength(height,thisColShift,c);

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    _localMatrix(thisColShift+i*c,j) = 
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
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[MR,* ]");
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
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[* ,MC]");
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
                        ( _grid->VCRank(), _rowAlignment, _grid->Size() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( RowAlignment() % _grid->Height() == A.RowAlignment() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int rowShift = RowShift();
        const int rowShiftOfA = A.RowShift();
        const int rowOffset = (rowShift-rowShiftOfA) / r;

        const int height = Height();
        const int localWidth = LocalWidth();

        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                _localMatrix(i,j) = A.LocalEntry(i,rowOffset+j*c);
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [* ,VC] <- [* ,MC]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int row = _grid->MCRank();
        const int col = _grid->MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        // We will SendRecv A[*,VC] within our process column to fix alignments.
        const int sendRow = (row+r+(rowAlignment%r)-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-(rowAlignment%r)) % r;
        const int sendRank = sendRow + r*col;

        const int sendRowShift = Shift( sendRank, rowAlignment, p );
        const int sendRowOffset = (sendRowShift-rowShiftOfA) / r;

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
                sendBuffer[i+j*height] = A.LocalEntry(i,sendRowOffset+j*c);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRow, 0,
          recvBuffer, recvSize, recvRow, MPI_ANY_TAG, _grid->MCComm() );

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
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[VC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
        ( new DistMatrix<T,MC,MR>(*_grid) );
    *A_MC_MR = A;

    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
        ( new DistMatrix<T,Star,VR>(*_grid) );
    *A_Star_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[* ,VC]");
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
            cout << "Unaligned [* ,VC] <- [* ,VC]." << endl;
#endif
        const int rank = _grid->VCRank();
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
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, _grid->VCComm() );

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
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[VR,* ]");
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
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[* ,VR]");
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

    // Compute which colmajor rank has the rowShift equal to our rowShiftOfA
    const int sendRankCM = (rankCM+(p+rowShiftOfA-rowShift)) % p;

    // Compute which colmajor rank has the A rowShift that we need
    const int recvRankRM = (rankRM+(p+rowShift-rowShiftOfA)) % p;
    const int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

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
    ( sendBuffer, sendSize, sendRankCM, 0,
      recvBuffer, recvSize, recvRankCM, MPI_ANY_TAG, _grid->VCComm() );

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
const DistMatrix<T,Star,VC>&
Elemental::DistMatrix<T,Star,VC>::operator=
( const DistMatrix<T,Star,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC] = DistMatrix[* ,* ]");
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

template class Elemental::DistMatrix<int,     Star,VC>;
template class Elemental::DistMatrix<float,   Star,VC>;
template class Elemental::DistMatrix<double,  Star,VC>;
#ifndef WITHOUT_COMPLEX
template class Elemental::DistMatrix<scomplex,Star,VC>;
template class Elemental::DistMatrix<dcomplex,Star,VC>;
#endif

