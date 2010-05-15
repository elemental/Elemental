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
Elemental::DistMatrix<T,MC,MR>::Print( const string& msg ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::Print");
#endif
    const int r = _grid->Height();
    const int c = _grid->Width();

    if( _grid->VCRank() == 0 && msg != "" )
        cout << msg << endl;

    const int height = Height();
    const int width  = Width();
    const int localHeight = LocalHeight();
    const int localWidth  = LocalWidth();
    const int colShift = ColShift();
    const int rowShift = RowShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Fill the send buffer: zero it then place our entries into their 
    // appropriate locations
    T* sendBuf = new T[height*width];
    for( int i=0; i<height*width; ++i )
        sendBuf[i] = (T)0;
    for( int i=0; i<localHeight; ++i )
        for( int j=0; j<localWidth; ++j )
            sendBuf[colShift+i*r + (rowShift+j*c)*height] = _localMatrix(i,j);

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
        delete recvBuf;
    }
    Barrier( _grid->VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignWith(DistMatrix[MC,MR])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment();
    _rowAlignment = A.RowAlignment();
    _colShift     = A.ColShift();
    _rowShift     = A.RowShift();
    _constrainedColDist = true;
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignWith(DistMatrix[MC,* ])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment();
    _colShift     = A.ColShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignWith
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignWith(DistMatrix[* ,MR])");
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
Elemental::DistMatrix<T,MC,MR>::AlignWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignWith(DistMatrix[MR,MC])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment();
    _rowAlignment = A.ColAlignment();
    _colShift     = A.RowShift();
    _rowShift     = A.ColShift();
    _constrainedColDist = true;
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignWith
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignWith(DistMatrix[MR,* ])");
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
Elemental::DistMatrix<T,MC,MR>::AlignWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignWith(DistMatrix[* ,MC])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment();
    _colShift     = A.RowShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignWith
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignWith(DistMatrix[VC,* ])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment() % _grid->Height();
    _colShift     = Shift( _grid->MCRank(), _colAlignment, _grid->Height() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignWith
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignWith(DistMatrix[* ,VC])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment() % _grid->Height();
    _colShift     = Shift( _grid->MCRank(), _colAlignment, _grid->Height() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignWith
( const DistMatrix<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignWith(DistMatrix[VR,* ])");
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
Elemental::DistMatrix<T,MC,MR>::AlignWith
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignWith(DistMatrix[* ,VR])");
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
Elemental::DistMatrix<T,MC,MR>::AlignColsWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignColsWith(DistMatrix[MC,MR])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment();
    _colShift     = A.ColShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignColsWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignColsWith(DistMatrix[MC,* ])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment();
    _colShift     = A.ColShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignColsWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignColsWith(DistMatrix[MR,MC])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment();
    _colShift     = A.RowShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignColsWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignColsWith(DistMatrix[* ,MC])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment();
    _colShift     = A.RowShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignColsWith
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignColsWith(DistMatrix[VC,* ])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment() % _grid->Height();
    _colShift     = Shift( _grid->MCRank(), _colAlignment, _grid->Height() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignColsWith
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignColsWith(DistMatrix[* ,VC])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment() % _grid->Height();
    _colShift     = Shift( _grid->MCRank(), _colAlignment, _grid->Height() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::AlignRowsWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignRowsWith(DistMatrix[MC,MR])");
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
Elemental::DistMatrix<T,MC,MR>::AlignRowsWith
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignRowsWith(DistMatrix[* ,MR])");
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
Elemental::DistMatrix<T,MC,MR>::AlignRowsWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignRowsWith(DistMatrix[MR,MC])");
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
Elemental::DistMatrix<T,MC,MR>::AlignRowsWith
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignRowsWith(DistMatrix[MR,* ])");
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
Elemental::DistMatrix<T,MC,MR>::AlignRowsWith
( const DistMatrix<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignRowsWith(DistMatrix[VR,* ])");
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
Elemental::DistMatrix<T,MC,MR>::AlignRowsWith
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::AlignRowsWith(DistMatrix[* ,VR])");
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
Elemental::DistMatrix<T,MC,MR>::ConformWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ConformWith(DistMatrix[MC,* ])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_CONFORMING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment();
    _colShift     = A.ColShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::ConformWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ConformWith(DistMatrix[* ,MC])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_CONFORMING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment();
    _colShift     = A.RowShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::ConformWith
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ConformWith(DistMatrix[MR,* ])");
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
Elemental::DistMatrix<T,MC,MR>::ConformWith
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ConformWith(DistMatrix[* ,MR])");
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
Elemental::DistMatrix<T,MC,MR>::FreeConstraints()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::FreeConstraints");
#endif
    _constrainedColDist = false;
    _constrainedRowDist = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::View
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::View(A)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _colAlignment = A.ColAlignment();
    _rowAlignment = A.RowAlignment();
    _colShift     = A.ColShift();
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
Elemental::DistMatrix<T,MC,MR>::LockedView
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::LockedView(A)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _colAlignment = A.ColAlignment();
    _rowAlignment = A.RowAlignment();
    _colShift     = A.ColShift();
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
Elemental::DistMatrix<T,MC,MR>::View
( DistMatrix<T,MC,MR>& A, 
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::View(A,i,j,height,width)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int r   = _grid->Height();
        const int c   = _grid->Width();
        const int row = _grid->MCRank();
        const int col = _grid->MRRank();

        _colAlignment = (A.ColAlignment()+i) % r;
        _rowAlignment = (A.RowAlignment()+j) % c;
  
        _colShift = Shift( row, _colAlignment, r );
        _rowShift = Shift( col, _rowAlignment, c );

        const int localHeightBehind = LocalLength(i,A.ColShift(),r);
        const int localWidthBehind  = LocalLength(j,A.RowShift(),c);

        const int localHeight = LocalLength( height, _colShift, r );
        const int localWidth  = LocalLength( width,  _rowShift, c );

        _localMatrix.View
        ( A.LocalMatrix(), localHeightBehind, localWidthBehind,
                           localHeight,       localWidth       );
    }
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::LockedView
( const DistMatrix<T,MC,MR>& A, 
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::LockedView(A,i,j,height,width)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int r   = _grid->Height();
        const int c   = _grid->Width();
        const int row = _grid->MCRank();
        const int col = _grid->MRRank();

        _colAlignment = (A.ColAlignment()+i) % r;
        _rowAlignment = (A.RowAlignment()+j) % c;
  
        _colShift = Shift( row, _colAlignment, r );
        _rowShift = Shift( col, _rowAlignment, c );

        const int localHeightBehind = LocalLength(i,A.ColShift(),r);
        const int localWidthBehind  = LocalLength(j,A.RowShift(),c);

        const int localHeight = LocalLength( height, _colShift, r );
        const int localWidth  = LocalLength( width,  _rowShift, c );

        _localMatrix.LockedView
        ( A.LockedLocalMatrix(), localHeightBehind, localWidthBehind,
                                 localHeight,       localWidth       );
    }
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::View1x2
( DistMatrix<T,MC,MR>& AL, 
  DistMatrix<T,MC,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::View1x2");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
    if( AL.ColAlignment() != AR.ColAlignment() )
        throw "1x2 is misaligned, cannot combine.";
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _colAlignment = AL.ColAlignment();
    _rowAlignment = AL.RowAlignment();
    _colShift     = AL.ColShift();
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
Elemental::DistMatrix<T,MC,MR>::LockedView1x2
( const DistMatrix<T,MC,MR>& AL, 
  const DistMatrix<T,MC,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::LockedView1x2");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
    if( AL.ColAlignment() != AR.ColAlignment() )
        throw "1x2 is misaligned, cannot combine.";
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _colAlignment = AL.ColAlignment();
    _rowAlignment = AL.RowAlignment();
    _colShift     = AL.ColShift();
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
Elemental::DistMatrix<T,MC,MR>::View2x1
( DistMatrix<T,MC,MR>& AT,
  DistMatrix<T,MC,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::View2x1");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
    if( AT.RowAlignment() != AB.RowAlignment() )
        throw "2x1 is misaligned, cannot combine.";
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _colAlignment = AT.ColAlignment();
    _rowAlignment = AT.RowAlignment();
    _colShift     = AT.ColShift();
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
Elemental::DistMatrix<T,MC,MR>::LockedView2x1
( const DistMatrix<T,MC,MR>& AT,
  const DistMatrix<T,MC,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::LockedView2x1");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
    if( AT.RowAlignment() != AB.RowAlignment() )
        throw "2x1 is misaligned, cannot combine.";
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _colAlignment = AT.ColAlignment();
    _rowAlignment = AT.RowAlignment();
    _colShift     = AT.ColShift();
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
Elemental::DistMatrix<T,MC,MR>::View2x2
( DistMatrix<T,MC,MR>& ATL, 
  DistMatrix<T,MC,MR>& ATR,
  DistMatrix<T,MC,MR>& ABL,
  DistMatrix<T,MC,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::View2x2");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
    if( ATL.ColAlignment() != ATR.ColAlignment() ||
        ABL.ColAlignment() != ABR.ColAlignment() ||
        ATL.RowAlignment() != ABL.RowAlignment() ||
        ATR.RowAlignment() != ABR.RowAlignment()   )
    {
        throw "2x2 set of matrices must aligned to combine.";
    }
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _colAlignment = ATL.ColAlignment();
    _rowAlignment = ATL.RowAlignment();
    _colShift     = ATL.ColShift();
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
Elemental::DistMatrix<T,MC,MR>::LockedView2x2
( const DistMatrix<T,MC,MR>& ATL, 
  const DistMatrix<T,MC,MR>& ATR,
  const DistMatrix<T,MC,MR>& ABL,
  const DistMatrix<T,MC,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::LockedView2x2");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
    if( ATL.ColAlignment() != ATR.ColAlignment() ||
        ABL.ColAlignment() != ABR.ColAlignment() ||
        ATL.RowAlignment() != ABL.RowAlignment() ||
        ATR.RowAlignment() != ABR.RowAlignment()   )
    {
        throw "2x2 set of matrices must align to combine.";
    }
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _colAlignment = ATL.ColAlignment();
    _rowAlignment = ATL.RowAlignment();
    _colShift     = ATL.ColShift();
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
Elemental::DistMatrix<T,MC,MR>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ResizeTo");
    CHECK_IF_LOCKED_VIEW;
#endif
    _height = height;
    _width  = width;
    _localMatrix.ResizeTo( LocalLength(height,_colShift,_grid->Height()),
                           LocalLength(width, _rowShift,_grid->Width())  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
Elemental::DistMatrix<T,MC,MR>::Get
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::Get");
    if( i<0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process _grid
    const int ownerRow = (i + ColAlignment()) % _grid->Height();
    const int ownerCol = (j + RowAlignment()) % _grid->Width();
    const int ownerRank = ownerRow + ownerCol*_grid->Height();

    T u;
    if( _grid->VCRank() == ownerRank )
    {
        const int iLoc = (i-ColShift()) / _grid->Height();
        const int jLoc = (j-RowShift()) / _grid->Width();
        u = _localMatrix(iLoc,jLoc);
    }
    Broadcast( &u, 1, ownerRank, _grid->VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::Set");
    if( i < 0 || i >= Height() || j < 0 || j >=Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const int ownerRow = (i + ColAlignment()) % _grid->Height();
    const int ownerCol = (j + RowAlignment()) % _grid->Width();
    const int ownerRank = ownerRow + ownerCol*_grid->Height();

    if( _grid->VCRank() == ownerRank )
    {
        const int iLoc = (i-ColShift()) / _grid->Height();
        const int jLoc = (j-RowShift()) / _grid->Width();
        _localMatrix(iLoc,jLoc) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::GetDiagonal
( DistMatrix<T,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::GetDiagonal");
    CHECK_IF_LOCKED_VIEW;
#endif
    int width = Width();
    int height = Height();
    int length;
    if( offset > 0 )
    {
        const int remainingWidth = max(width-offset,0);
        length = min(height,remainingWidth);
    }
    else
    {
        const int remainingHeight = max(height+offset,0);
        length = min(remainingHeight,width);
    }
#ifndef RELEASE
    if( d.Viewing()  && length != d.Height() )
        throw "d is not of the correct length.";
#endif

    if( ! d.Viewing() )
    {
        if( ! d.ConstrainedColDist() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( length, 1 );
    }

    if( d.InDiagonal() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int lcm = _grid->LCM();
        const int colShift = ColShift();
        const int rowShift = RowShift();
        const int diagShift = d.ColShift();

        int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
        for( int kLoc=0; kLoc<localDiagLength; ++kLoc )
            d.LocalEntry(kLoc,0) = LocalEntry(iLocStart+kLoc*(lcm/r),
                                              jLocStart+kLoc*(lcm/c));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::GetDiagonal
( DistMatrix<T,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::GetDiagonal");
    CHECK_IF_LOCKED_VIEW;
#endif
    int height = Height();
    int width = Width();
    int length;
    if( offset > 0 )
    {
        const int remainingWidth = max(width-offset,0);
        length = min(height,remainingWidth);
    }
    else
    {
        const int remainingHeight = max(height+offset,0);
        length = min(remainingHeight,width);
    }
#ifndef RELEASE
    if( d.Viewing()  && length != d.Width() )
        throw "d is not of the correct length.";
#endif

    if( ! d.Viewing() )
    {
        if( ! d.ConstrainedRowDist() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( 1, length );
    }

    if( d.InDiagonal() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int lcm = _grid->LCM();
        const int colShift = ColShift();
        const int rowShift = RowShift();
        const int diagShift = d.RowShift();

        int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
        for( int kLoc=0; kLoc<localDiagLength; ++kLoc )
            d.LocalEntry(0,kLoc) = LocalEntry(iLocStart+kLoc*(lcm/r),
                                              jLocStart+kLoc*(lcm/c));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::SetDiagonal
( const DistMatrix<T,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::SetDiagonal");
    if( d.Width() != 1 )
        throw "d must be a column vector.";
    {
        int height = Height();
        int width = Width();
        int length;
        if( offset >= 0 )
        {
            const int remainingWidth = max(width-offset,0);
            length = min(remainingWidth,height);
        }
        else
        {
            const int remainingHeight = max(height+offset,0);
            length = min(remainingHeight,width);
        }
        if( length != d.Height() )
        {
            ostringstream msg;
            msg << "d is not of the same length as the diagonal:" << endl
                << "  A ~ " << Height() << " x " << Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            const string& s = msg.str();
            throw s.c_str();
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int lcm = _grid->LCM();
        const int colShift = ColShift();
        const int rowShift = RowShift();
        const int diagShift = d.ColShift();

        int iStart,jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
        for( int kLoc=0; kLoc<localDiagLength; ++kLoc )
            LocalEntry(iLocStart+kLoc*(lcm/r),
                       jLocStart+kLoc*(lcm/c)) = d.LocalEntry(kLoc,0);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::SetDiagonal
( const DistMatrix<T,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::SetDiagonal");
    if( d.Height() != 1 )
        throw "d must be a row vector.";
    {
        int height = Height();
        int width = Width();
        int length;
        if( offset >= 0 )
        {
            const int remainingWidth = max(width-offset,0);
            length = min(remainingWidth,height);
        }
        else
        {
            const int remainingHeight = max(height+offset,0);
            length = min(remainingHeight,width);
        }
        if( length != d.Width() )
        {
            ostringstream msg;
            msg << "d is not of the same length as the diagonal:" << endl
                << "  A ~ " << Height() << " x " << Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            const string& s = msg.str();
            throw s.c_str();
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int lcm = _grid->LCM();
        const int colShift = ColShift();
        const int rowShift = RowShift();
        const int diagShift = d.RowShift();

        int iStart,jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
        for( int kLoc=0; kLoc<localDiagLength; ++kLoc )
            LocalEntry(iLocStart+kLoc*(lcm/r),
                       jLocStart+kLoc*(lcm/c)) = d.LocalEntry(0,kLoc);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Utility functions, e.g., SetToIdentity and MakeTrapezoidal                 //
//----------------------------------------------------------------------------//

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::MakeTrapezoidal");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int width = Width();
    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    const int r = _grid->Height();
    const int c = _grid->Width();
    const int colShift = ColShift();
    const int rowShift = RowShift();

    if( shape == Lower )
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*c;
            int lastZero_i;
            if( side == Left )
                lastZero_i = j-offset-1;
            else
                lastZero_i = j-offset+height-width-1;
            if( lastZero_i >= 0 )
            {
                const int boundary = min( lastZero_i+1, height );
                const int numZeros = LocalLength( boundary, colShift, r );
                for( int iLoc=0; iLoc<numZeros; ++iLoc )
                    _localMatrix(iLoc,jLoc) = (T)0;
            }
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
                firstZero_i = max(j+height-width-offset+1,0);
            const int nonzeroLength = LocalLength(firstZero_i,colShift,r);
            for( int iLoc=nonzeroLength; iLoc<localHeight; ++iLoc )
                _localMatrix(iLoc,jLoc) = (T)0;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::SetToIdentity");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    const int r = _grid->Height();
    const int c = _grid->Width();
    const int colShift = ColShift();
    const int rowShift = RowShift();

    _localMatrix.SetToZero();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;                
        if( i % c == rowShift )
        {
            const int jLoc = (i-rowShift) / c;
            if( jLoc < localWidth )
                _localMatrix(iLoc,jLoc) = (T)1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::SetToRandom");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    for( int i=0; i<localHeight; ++i )
        for( int j=0; j<localWidth; ++j )
            _localMatrix(i,j) = Random<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::SetToRandomDiagDominant()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::SetToRandomDiagDominant");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int r = _grid->Height();
    const int c = _grid->Width();

    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    const int colShift = ColShift();
    const int rowShift = RowShift();

    SetToRandom();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;                
        if( i % c == rowShift )
        {
            const int jLoc = (i-rowShift) / c;
            if( jLoc < localWidth )
                _localMatrix(iLoc,jLoc) += (T)max(Height(),Width());
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::SetToZero");
    CHECK_IF_LOCKED_VIEW;
#endif
    _localMatrix.SetToZero();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::ConjugateTransposeFrom
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ConjugateTransposeFrom");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    if( _viewing && ( Height() != A.Width() || Width() != A.Height() ) )
        throw "Cannot resize a view.";
#endif
    if( !_viewing )
    {
        if( ! ConstrainedColDist() )
        {
            _colAlignment = A.RowAlignment();
            _colShift = Shift
                        ( _grid->MCRank(), _colAlignment, _grid->Height() );
        }
        ResizeTo( A.Width(), A.Height() );
    }

    if( ColAlignment() == A.RowAlignment() )
    {
        const int c = _grid->Width();
        const int rowShift = RowShift();

        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = Conj( A.LocalEntry(rowShift+j*c,i) );
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MC,MR]::ConjugateTransposeFrom." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int rank = _grid->MCRank();
        const int rowShift = RowShift();
        const int colAlignment = ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+r+colAlignment-rowAlignmentOfA) % r;
        const int recvRank = (rank+r+rowAlignmentOfA-colAlignment) % r;

        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = localWidthOfA * localWidth;
        const int recvSize = localHeight * localWidth;

        _auxMemory.Require( sendSize + recvSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localWidthOfA; ++i )
                sendBuffer[i+j*localWidth] = 
                    Conj( A.LocalEntry(rowShift+j*c,i) );

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, _grid->MCComm() );

        // Unpack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = recvBuffer[i+j*localHeight];

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::TransposeFrom
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::TransposeFrom");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    if( _viewing && ( Height() != A.Width() || Width() != A.Height() ) )
        throw "Cannot resize a view.";
#endif
    if( !_viewing )
    {
        if( ! ConstrainedColDist() )
        {
            _colAlignment = A.RowAlignment();
            _colShift = Shift
                        ( _grid->MCRank(), _colAlignment, _grid->Height() );
        }
        ResizeTo( A.Width(), A.Height() );
    }

    if( ColAlignment() == A.RowAlignment() )
    {
        const int c = _grid->Width();
        const int rowShift = RowShift();

        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = A.LocalEntry(rowShift+j*c,i);
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MC,MR]::TransposeFrom." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int rank = _grid->MCRank();
        const int rowShift = RowShift();
        const int colAlignment = ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+r+colAlignment-rowAlignmentOfA) % r;
        const int recvRank = (rank+r+rowAlignmentOfA-colAlignment) % r;

        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = localWidthOfA * localWidth;
        const int recvSize = localHeight * localWidth;

        _auxMemory.Require( sendSize + recvSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localWidthOfA; ++i )
                sendBuffer[i+j*localWidth] = A.LocalEntry(rowShift+j*c,i);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, _grid->MCComm() );

        // Unpack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = recvBuffer[i+j*localHeight];

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[MC,MR]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        if( ! ConstrainedColDist() )
        {
            _colAlignment = A.ColAlignment();
            _colShift = A.ColShift();
        }
        if( ! ConstrainedRowDist() )
        {
            _rowAlignment = A.RowAlignment();
            _rowShift = A.RowShift();
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( ColAlignment() == A.ColAlignment() &&
        RowAlignment() == A.RowAlignment()    )
    {
        _localMatrix = A.LockedLocalMatrix();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MC,MR] <- [MC,MR]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int row = _grid->MCRank();
        const int col = _grid->MRRank();

        const int colAlignment = ColAlignment();
        const int rowAlignment = RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
        const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;
        const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;
        const int sendRank = sendRow + sendCol*r;
        const int recvRank = recvRow + recvCol*r;

        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = localHeightOfA * localWidthOfA;
        const int recvSize = localHeight * localWidth;

        _auxMemory.Require( sendSize + recvSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, _grid->VCComm() );

        // Unpack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = recvBuffer[i+j*localHeight];

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[MC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        if( ! ConstrainedColDist() )
        {
            _colAlignment = A.ColAlignment();
            _colShift = Shift
                        ( _grid->MCRank(), _colAlignment, _grid->Height() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( ColAlignment() == A.ColAlignment() )
    {
        const int c = _grid->Width();
        const int rowShift = RowShift();

        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = A.LocalEntry(i,rowShift+j*c);
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MC,MR] <- [MC,* ]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int rank = _grid->MCRank();
        const int rowShift = RowShift();
        const int colAlignment = ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (rank+r+colAlignment-colAlignmentOfA) % r;
        const int recvRank = (rank+r+colAlignmentOfA-colAlignment) % r;

        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int sendSize = localHeightOfA * localWidth;
        const int recvSize = localHeight * localWidth;

        _auxMemory.Require( sendSize + recvSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i+j*localWidth] = A.LocalEntry(i,rowShift+j*c);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, _grid->MCComm() );

        // Unpack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = recvBuffer[i+j*localHeight];

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[* ,MR]");
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
        const int r = _grid->Height();
        const int colShift = ColShift();

        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = A.LocalEntry(colShift+i*r,j);
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MC,MR] <- [* ,MR]." << endl;
#endif
        const int r = _grid->Height(); 
        const int c = _grid->Width();
        const int col = _grid->MRRank();
        const int colShift = ColShift();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = localHeight * localWidthOfA;
        const int recvSize = localHeight * localWidth;

        _auxMemory.Require( sendSize + recvSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeight; ++i )
                sendBuffer[i+j*localHeight] = A.LocalEntry(colShift+i*r,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, MPI_ANY_TAG, _grid->MRComm() );

        // Unpack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = recvBuffer[i+j*localHeight];

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[MD,* ]");
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
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[* ,MD]");
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
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[MR,MC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( A.Width() == 1 )
    {
        if( !_viewing )
            ResizeTo( A.Height(), 1 );

        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int myRow = _grid->MCRank();
        const int myCol = _grid->MRRank();
        const int rankCM = _grid->VCRank();
        const int rankRM = _grid->VRRank();
        const int ownerCol = RowAlignment();
        const int ownerRow = A.RowAlignment();
        const int colAlignment = ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = A.Height();
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = max(maxLocalHeight,MinCollectContrib);

        const int colShiftVC = Shift(rankCM,colAlignment,p);
        const int colShiftVROfA = Shift(rankRM,colAlignmentOfA,p);
        const int sendRankCM = (rankCM+(p+colShiftVROfA-colShiftVC)) % p;
        const int recvRankRM = (rankRM+(p+colShiftVC-colShiftVROfA)) % p;
        const int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        _auxMemory.Require( (r+c)*portionSize );
        T* buffer = _auxMemory.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        if( myRow == ownerRow )
        {
            // Pack
            for( int k=0; k<r; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const int shift = Shift(myCol+c*k,colAlignmentOfA,p);
                const int offset = (shift-A.ColShift()) / c;
                const int thisLocalHeight = LocalLength(height,shift,p);

                for( int i=0; i<thisLocalHeight; ++i )
                    data[i] = A.LocalEntry(offset+i*r,0);
            }
        }

        // A[VR,* ] <- A[MR,MC]
        Scatter
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerRow, _grid->MCComm() );

        // A[VC,* ] <- A[VR,* ]
        SendRecv
        ( sendBuf, portionSize, sendRankCM, 0,
          recvBuf, portionSize, recvRankCM, MPI_ANY_TAG, _grid->VCComm() );

        // A[MC,MR] <- A[VC,* ]
        Gather
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerCol, _grid->MRComm() );

        if( myCol == ownerCol )
        {
            // Unpack
            for( int k=0; k<c; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const int shift = Shift(myRow+r*k,colAlignment,p);
                const int offset = (shift-ColShift()) / r;
                const int thisLocalHeight = LocalLength(height,shift,p);

                for( int i=0; i<thisLocalHeight; ++i )
                    _localMatrix(offset+i*c,0) = data[i];
            }
        }

        _auxMemory.Release();
    }
    else if( A.Height() == 1 )
    {
        if( !_viewing )
            ResizeTo( 1, A.Width() );

        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int myRow = _grid->MCRank();
        const int myCol = _grid->MRRank();
        const int rankCM = _grid->VCRank();
        const int rankRM = _grid->VRRank();
        const int ownerRow = ColAlignment();
        const int ownerCol = A.ColAlignment();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int width = A.Width();
        const int maxLocalWidth = MaxLocalLength(width,p);

        const int portionSize = max(maxLocalWidth,MinCollectContrib);

        const int rowShiftVR = Shift(rankRM,rowAlignment,p);
        const int rowShiftVCOfA = Shift(rankCM,rowAlignmentOfA,p);
        const int sendRankRM = (rankRM+(p+rowShiftVCOfA-rowShiftVR)) % p;
        const int recvRankCM = (rankCM+(p+rowShiftVR-rowShiftVCOfA)) % p;
        const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

        _auxMemory.Require( (r+c)*portionSize );
        T* buffer = _auxMemory.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        if( myCol == ownerCol )
        {
            // Pack
            for( int k=0; k<c; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const int shift = Shift(myRow+r*k,rowAlignmentOfA,p);
                const int offset = (shift-A.RowShift()) / r;
                const int thisLocalWidth = LocalLength(width,shift,p);

                for( int j=0; j<thisLocalWidth; ++j )
                    data[j] = A.LocalEntry(0,offset+j*c);
            }
        }

        // A[* ,VC] <- A[MR,MC]
        Scatter
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerCol, _grid->MRComm() );

        // A[* ,VR] <- A[* ,VC]
        SendRecv
        ( sendBuf, portionSize, sendRankRM, 0,
          recvBuf, portionSize, recvRankRM, MPI_ANY_TAG, _grid->VRComm() );

        // A[MC,MR] <- A[* ,VR]
        Gather
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerRow, _grid->MCComm() );

        if( myRow == ownerRow )
        {
            // Unpack
            for( int k=0; k<r; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const int shift = Shift(myCol+c*k,rowAlignment,p);
                const int offset = (shift-RowShift()) / c;
                const int thisLocalWidth = LocalLength(width,shift,p);

                for( int j=0; j<thisLocalWidth; ++j )
                    _localMatrix(0,offset+j*r) = data[j];
            }
        }

        _auxMemory.Release();
    }
    else
    {
        if( A.Height() >= A.Width() )
        {
            auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
                ( new DistMatrix<T,VR,Star>(*_grid) );

            *A_VR_Star = A;

            auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
                ( new DistMatrix<T,VC,Star>(true,ColAlignment(),*_grid) );
            *A_VC_Star = *A_VR_Star;
            delete A_VR_Star.release(); // lowers memory highwater

            *this = *A_VC_Star;
        }
        else
        {
            auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
                ( new DistMatrix<T,Star,VC>(*_grid) );
            *A_Star_VC = A;

            auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
                ( new DistMatrix<T,Star,VR>(true,RowAlignment(),*_grid) );
            *A_Star_VR = *A_Star_VC;
            delete A_Star_VC.release(); // lowers memory highwater

            *this = *A_Star_VR;
            this->ResizeTo( A_Star_VR->Height(), A_Star_VR->Width() );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[MR,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
        ( new DistMatrix<T,VR,Star>(*_grid) );
    *A_VR_Star = A;

    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
        ( new DistMatrix<T,VC,Star>(true,ColAlignment(),*_grid) );
    *A_VC_Star = *A_VR_Star;
    delete A_VR_Star.release(); // lowers memory highwater

    *this = *A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[* ,MC]");
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
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[VC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        if( ! ConstrainedColDist() )
        {
            _colAlignment = A.ColAlignment() % _grid->Height();
            _colShift = Shift
                        ( _grid->MCRank(), _colAlignment, _grid->Height() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( ColAlignment() == A.ColAlignment() % _grid->Height() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = r * c;
        const int row = _grid->MCRank();
        const int colShift = ColShift();
        const int rowAlignment = RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = Height();
        const int width = Width();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,c);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        _auxMemory.Require( 2*c*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[c*portionSize];

        // Pack
        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignment,c);
            const int thisLocalWidth = LocalLength(width,thisRowShift,c);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = A.LocalEntry(i,thisRowShift+j*c);
        }

        // Communicate
        AllToAll( sendBuffer, portionSize,
                  recvBuffer, portionSize, _grid->MRComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisRank = row+k*r;
            const int thisColShift = Shift(thisRank,colAlignmentOfA,p);
            const int thisColOffset = (thisColShift-colShift) / r;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    _localMatrix(thisColOffset+i*c,j) = 
                          data[i+j*thisLocalHeight];
        }

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MC,MR] <- [VC,* ]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = r * c;
        const int row = _grid->MCRank();
        const int colShift = ColShift();
        const int colAlignment = ColAlignment();
        const int rowAlignment = RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRow = (row+r+colAlignment-(colAlignmentOfA%r)) % r;
        const int recvRow = (row+r+(colAlignmentOfA%r)-colAlignment) % r;

        const int height = Height();
        const int width = Width();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,c);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        _auxMemory.Require( 2*c*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[c*portionSize];

        // Pack
        for( int k=0; k<c; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignment,c);
            const int thisLocalWidth = LocalLength(width,thisRowShift,c);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = A.LocalEntry(i,thisRowShift+j*c);
        }

        // SendRecv: properly align A[VC,*] via a trade in the column
        SendRecv
        ( secondBuffer, c*portionSize, sendRow, 0,
          firstBuffer,  c*portionSize, recvRow, 0, _grid->MCComm() );

        // AllToAll to gather all of the aligned A[VC,*] data into secondBuff.
        AllToAll
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, _grid->MRComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisRank = recvRow+k*r;
            const int thisColShift = Shift(thisRank,colAlignmentOfA,p);
            const int thisColOffset = (thisColShift-colShift) / r;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    _localMatrix(thisColOffset+i*c,j) = 
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
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[* ,VC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,Star,VR> A_Star_VR(true,RowAlignment(),*_grid);

    A_Star_VR = A;
    *this = A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[VR,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,VC,Star> A_VC_Star(true,ColAlignment(),*_grid);

    A_VC_Star = A;
    *this     = A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[* ,VR]");
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
        const int p = r * c;
        const int col = _grid->MRRank();
        const int rowShift = RowShift();
        const int colAlignment = ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int height = Height();
        const int width = Width();
        const int localHeight = LocalHeight();
        const int localWidthOfA = A.LocalWidth();

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

            const int thisColShift = Shift(k,colAlignment,r);
            const int thisLocalHeight = LocalLength(height,thisColShift,r);

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColShift+i*r,j);
        }

        // Communicate
        AllToAll( sendBuffer, portionSize,
                  recvBuffer, portionSize, _grid->MCComm() );

        // Unpack
        for( int k=0; k<r; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisRank = col+k*c;
            const int thisRowShift = Shift(thisRank,rowAlignmentOfA,p);
            const int thisRowOffset = (thisRowShift-rowShift) / c;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,thisRowOffset+j*r) = data[i+j*localHeight];
        }

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MC,MR] <- [* ,VR]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = r * c;
        const int col = _grid->MRRank();
        const int rowShift = RowShift();
        const int colAlignment = ColAlignment();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendCol = (col+c+rowAlignment-(rowAlignmentOfA%c)) % c;
        const int recvCol = (col+c+(rowAlignmentOfA%c)-rowAlignment) % c;

        const int height = Height();
        const int width = Width();
        const int localHeight = LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        
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

            const int thisColShift = Shift(k,colAlignment,r);
            const int thisLocalHeight = LocalLength(height,thisColShift,r);

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColShift+i*r,j);
        }

        // SendRecv: properly align A[*,VR] via a trade in the column
        SendRecv
        ( secondBuffer, r*portionSize, sendCol, 0,
          firstBuffer,  r*portionSize, recvCol, 0, _grid->MRComm() );

        // AllToAll to gather all of the aligned [*,VR] data into secondBuffer
        AllToAll
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, _grid->MCComm() );

        // Unpack
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisRank = recvCol+k*c;
            const int thisRowShift = Shift(thisRank,rowAlignmentOfA,p);
            const int thisRowOffset = (thisRowShift-rowShift) / c;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,thisRowOffset+j*r) = data[i+j*localHeight];
        }

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,MR>&
Elemental::DistMatrix<T,MC,MR>::operator=
( const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR] = DistMatrix[* ,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int r = _grid->Height();
    const int c = _grid->Width();
    const int colShift = ColShift();
    const int rowShift = RowShift();

    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = A.LocalEntry(colShift+i*r,rowShift+j*c);
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::ReduceScatterFrom
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ReduceScatterFrom([MC,* ])");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        if( ! ConstrainedColDist() )
        {
            _colAlignment = A.ColAlignment();
            _colShift = Shift
                        ( _grid->MCRank(), _colAlignment, _grid->Height() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( ColAlignment() == A.ColAlignment() )
    {
        if( Width() == 1 )
        {
            const int rowAlignment = RowAlignment();
            const int myCol = _grid->MRRank();

            const int localHeight = LocalHeight();

            const int recvSize = max(localHeight,MinCollectContrib);
            const int sendSize = recvSize;

            _auxMemory.Require( sendSize + recvSize );

            T* buffer = _auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];

            // Pack 
            for( int i=0; i<localHeight; ++i )
                sendBuffer[i] = A.LocalEntry(i,0);

            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, sendSize,
              MPI_SUM, rowAlignment, _grid->MRComm() );

            if( myCol == rowAlignment )
            {
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,0) = recvBuffer[i];
            }

            _auxMemory.Release();
        }
        else
        {
            const int c = _grid->Width();
            const int rowAlignment = RowAlignment();
        
            const int width = Width();
            const int localHeight = LocalHeight();
            const int localWidth = LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,c);

            const int recvSize = max(localHeight*maxLocalWidth,
                                     MinCollectContrib         );
            const int sendSize = c * recvSize;

            _auxMemory.Require( sendSize + recvSize );

            T* buffer = _auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];
        
            // Pack 
            int* recvSizes = new int[c];
            for( int k=0; k<c; ++k )
            {
                T* data = &sendBuffer[k*recvSize];
                recvSizes[k] = recvSize;

                const int thisRowShift = Shift( k, rowAlignment, c );
                const int thisLocalWidth = LocalLength(width,thisRowShift,c);

                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        data[i+j*localHeight] = 
                            A.LocalEntry(i,thisRowShift+j*c);
            }

            // Reduce-scatter over each process row
            ReduceScatter
            ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, 
              _grid->MRComm()                       );
            delete[] recvSizes;

            // Unpack our received data
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,j) = recvBuffer[i+j*localHeight];

            _auxMemory.Release();
        }
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned ReduceScatterFrom [MC,MR] <- [MC,* ]." << endl;
#endif
        if( Width() == 1 )
        {
            const int r = _grid->Height();
            const int rowAlignment = RowAlignment();
            const int myRow = _grid->MCRank();
            const int myCol = _grid->MRRank();

            const int height = Height();
            const int localHeight = LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int portionSize = max(maxLocalHeight,MinCollectContrib);

            const int colAlignment = ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendRow = (myRow+r+colAlignment-colAlignmentOfA) % r;
            const int recvRow = (myRow+r+colAlignmentOfA-colAlignment) % r;

            _auxMemory.Require( 2*portionSize );

            T* buffer = _auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack 
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i] = A.LocalEntry(i,0);
        
            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize, 
              MPI_SUM, rowAlignment, _grid->MRComm() );

            if( myCol == rowAlignment )
            {
                // Perform the realignment
                SendRecv
                ( recvBuffer, portionSize, sendRow, 0,
                  sendBuffer, portionSize, recvRow, 0, _grid->MCComm() );

                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,0) = sendBuffer[i];
            }

            _auxMemory.Release();
        }
        else
        {
            const int r = _grid->Height();
            const int c = _grid->Width();
            const int row = _grid->MCRank();

            const int colAlignment = ColAlignment();
            const int rowAlignment = RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
            const int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

            const int width = Width();
            const int localHeight = LocalHeight();
            const int localWidth = LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalWidth = MaxLocalLength(width,c);

            const int recvSize_RS = max(localHeightOfA * maxLocalWidth,
                                        MinCollectContrib              );
            const int sendSize_RS = c * recvSize_RS;
            const int recvSize_SR = localHeight * localWidth;

            _auxMemory.Require( recvSize_RS + max(sendSize_RS,recvSize_SR) );

            T* buffer = _auxMemory.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[recvSize_RS];

            // Pack 
            int* recvSizes = new int[c];
            for( int k=0; k<c; ++k )
            {
                T* data = &secondBuffer[k*recvSize_RS];
                recvSizes[k] = recvSize_RS;

                const int thisRowShift = Shift( k, rowAlignment, c );
                const int thisLocalWidth = LocalLength(width,thisRowShift,c);

                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeightOfA; ++i )
                        data[i+j*localHeightOfA] = 
                            A.LocalEntry(i,thisRowShift+j*c);
            }

            // Reduce-scatter over each process row
            ReduceScatter
            ( secondBuffer, firstBuffer, recvSizes, MPI_SUM, 
              _grid->MRComm()                          );
            delete[] recvSizes;

            // Trade reduced data with the appropriate process row
            SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendRow, 0,
              secondBuffer, localHeight*localWidth,    recvRow, 0, 
              _grid->MCComm()                                );

            // Unpack the received data
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,j) = secondBuffer[i+j*localHeight];

            _auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::ReduceScatterFrom
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ReduceScatterFrom([* ,MR])");
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
        const int r = _grid->Height();
        const int colAlignment = ColAlignment();

        const int height = Height();
        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);

        const int recvSize = max(maxLocalHeight * localWidth,
                                 MinCollectContrib           );
        const int sendSize = r * recvSize;

        _auxMemory.Require( sendSize + recvSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack 
        int* recvSizes = new int[r];
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColShift+i*r,j);
        }

        // Reduce-scatter over each process col
        ReduceScatter
        ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, _grid->MCComm() );
        delete[] recvSizes;

        // Unpack our received data
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = recvBuffer[i+j*localHeight];

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned ReduceScatterFrom [MC,MR] <- [* ,MR]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int col = _grid->MRRank();

        const int colAlignment = ColAlignment();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

        const int height = Height();
        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);

        const int recvSize_RS = max(maxLocalHeight * localWidthOfA,
                                    MinCollectContrib              );
        const int sendSize_RS = r * recvSize_RS;
        const int recvSize_SR = localHeight * localWidth;

        _auxMemory.Require( recvSize_RS + max(sendSize_RS,recvSize_SR) );

        T* buffer = _auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[recvSize_RS];

        // Pack 
        int* recvSizes = new int[r];
        for( int k=0; k<r; ++k )
        {
            T* data = &secondBuffer[k*recvSize_RS];
            recvSizes[k] = recvSize_RS;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColShift+i*r,j);
        }

        // Reduce-scatter over each process col
        ReduceScatter
        ( secondBuffer, firstBuffer, recvSizes, MPI_SUM, _grid->MCComm() );
        delete[] recvSizes;

        // Trade reduced data with the appropriate process col
        SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendCol, 0,
          secondBuffer, localHeight*localWidth,    recvCol, 0,
          _grid->MRComm()                                );

        // Unpack the received data
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = secondBuffer[i+j*localHeight];
        
        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::ReduceScatterFrom
( const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ReduceScatterFrom([* ,* ])");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        ResizeTo( A.Height(), A.Width() );
    }

    const int r = _grid->Height();
    const int c = _grid->Width();
    const int colAlignment = ColAlignment();
    const int rowAlignment = RowAlignment();

    const int height = Height();
    const int width = Width();
    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    const int maxLocalHeight = MaxLocalLength(height,r);
    const int maxLocalWidth = MaxLocalLength(width,c);

    const int recvSize = max(maxLocalHeight*maxLocalWidth,MinCollectContrib);
    const int sendSize = r * c * recvSize;

    _auxMemory.Require( sendSize + recvSize );

    T* buffer = _auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack 
    int* recvSizes = new int[r*c];
    for( int l=0; l<c; ++l )
    {
        const int thisRowShift = Shift( l, rowAlignment, c );
        const int thisLocalWidth = LocalLength( width, thisRowShift, c );

        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[(k+l*r)*recvSize];
            recvSizes[k+l*r] = recvSize;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColShift+i*r,thisRowShift+j*c);
        }
    }

    // Reduce-scatter over each process col
    ReduceScatter
    ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, _grid->VCComm() );
    delete[] recvSizes;

    // Unpack our received data
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = recvBuffer[i+j*localHeight];

    _auxMemory.Release();

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::ReduceScatterUpdate
( T alpha, const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ReduceScatterUpdate([MC,* ])");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_DIFF_SIZE( A );
    CHECK_IF_REDIST_DIFF_GRID( A );
#endif
    if( ColAlignment() == A.ColAlignment() )
    {
        if( Width() == 1 )
        {
            const int rowAlignment = RowAlignment();
            const int myCol = _grid->MRRank();

            const int localHeight = LocalHeight();

            const int portionSize = max(localHeight,MinCollectContrib);

            _auxMemory.Require( 2*portionSize );

            T* buffer = _auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack 
            for( int i=0; i<localHeight; ++i )
                sendBuffer[i] = A.LocalEntry(i,0);
        
            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize, 
              MPI_SUM, rowAlignment, _grid->MRComm() );

            if( myCol == rowAlignment )
            {
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,0) += alpha*recvBuffer[i];
            }

            _auxMemory.Release();
        }
        else
        {
            const int c = _grid->Width();
            const int rowAlignment = RowAlignment();

            const int width = Width();
            const int localHeight = LocalHeight();
            const int localWidth = LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,c);

            const int portionSize = max(localHeight * maxLocalWidth,
                                        MinCollectContrib           );

            _auxMemory.Require( (c+1)*portionSize );

            T* buffer = _auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[c*portionSize];

            // Pack 
            int* recvSizes = new int[c];
            for( int k=0; k<c; ++k )
            {
                T* data = &sendBuffer[k*portionSize];
                recvSizes[k] = portionSize;

                const int thisRowShift = Shift( k, rowAlignment, c );
                const int thisLocalWidth = LocalLength(width,thisRowShift,c);

                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        data[i+j*localHeight] = 
                            A.LocalEntry(i,thisRowShift+j*c);
            }
        
            // Reduce-scatter over each process row
            ReduceScatter
            ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, 
              _grid->MRComm()                       );
            delete[] recvSizes;

            // Update with our received data
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,j) += alpha*recvBuffer[i+j*localHeight];

            _auxMemory.Release();
        }
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned ReduceScatterUpdate [MC,MR] <- [MC,* ]." << endl;
#endif
        if( Width() == 1 )
        {
            const int r = _grid->Height();
            const int rowAlignment = RowAlignment();
            const int myRow = _grid->MCRank();
            const int myCol = _grid->MRRank();

            const int height = Height();
            const int localHeight = LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int portionSize = max(maxLocalHeight,MinCollectContrib);

            const int colAlignment = ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendRow = (myRow+r+colAlignment-colAlignmentOfA) % r;
            const int recvRow = (myRow+r+colAlignmentOfA-colAlignment) % r;

            _auxMemory.Require( 2*portionSize );

            T* buffer = _auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack 
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i] = A.LocalEntry(i,0);
        
            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize, 
              MPI_SUM, rowAlignment, _grid->MRComm() );

            if( myCol == rowAlignment )
            {
                // Perform the realignment
                SendRecv
                ( recvBuffer, portionSize, sendRow, 0,
                  sendBuffer, portionSize, recvRow, 0, _grid->MCComm() );

                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,0) += alpha*sendBuffer[i];
            }

            _auxMemory.Release();
        }
        else
        {
            const int r = _grid->Height();
            const int c = _grid->Width();
            const int row = _grid->MCRank();

            const int colAlignment = ColAlignment();
            const int rowAlignment = RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
            const int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

            const int width = Width();
            const int localHeight = LocalHeight();
            const int localWidth = LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalWidth = MaxLocalLength(width,c);

            const int recvSize_RS = max(localHeightOfA * maxLocalWidth,
                                        MinCollectContrib              );
            const int sendSize_RS = c * recvSize_RS;
            const int recvSize_SR = localHeight * localWidth;

            _auxMemory.Require( recvSize_RS + max(sendSize_RS,recvSize_SR) );

            T* buffer = _auxMemory.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[recvSize_RS];

            // Pack 
            int* recvSizes = new int[c];
            for( int k=0; k<c; ++k )
            {
                T* data = &secondBuffer[k*recvSize_RS];
                recvSizes[k] = recvSize_RS;

                const int thisRowShift = Shift( k, rowAlignment, c );
                const int thisLocalWidth = LocalLength(width,thisRowShift,c);

                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeightOfA; ++i )
                        data[i+j*localHeightOfA] = 
                            A.LocalEntry(i,thisRowShift+j*c);
            }

            // Reduce-scatter over each process row
            ReduceScatter
            ( secondBuffer, firstBuffer, recvSizes, MPI_SUM, 
              _grid->MRComm()                          );
            delete[] recvSizes;

            // Trade reduced data with the appropriate process row
            SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendRow, 0,
              secondBuffer, localHeight*localWidth,    recvRow, 0, 
              _grid->MCComm()                                );

            // Update with our received data
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,j) += alpha*secondBuffer[i+j*localHeight];

            _auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::ReduceScatterUpdate
( T alpha, const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ReduceScatterUpdate([* ,MR])");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_DIFF_SIZE( A );
    CHECK_IF_REDIST_DIFF_GRID( A );
#endif
    if( RowAlignment() == A.RowAlignment() )
    {
        const int r = _grid->Height();
        const int colAlignment = ColAlignment();

        const int height = Height();
        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);

        const int recvSize = max(maxLocalHeight*localWidth,MinCollectContrib);
        const int sendSize = r * recvSize;

        _auxMemory.Require( sendSize + recvSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack 
        int* recvSizes = new int[r];
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColShift+i*r,j);
        }

        // Reduce-scatter over each process col
        ReduceScatter
        ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, _grid->MCComm() );
        delete[] recvSizes;

        // Update with our received data
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) += alpha*recvBuffer[i+j*localHeight];

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned ReduceScatterUpdate [MC,MR] <- [* ,MR]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int col = _grid->MRRank();

        const int colAlignment = ColAlignment();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

        const int height = Height();
        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);

        const int recvSize_RS = max(maxLocalHeight * localWidthOfA,
                                    MinCollectContrib              );
        const int sendSize_RS = r * recvSize_RS;
        const int recvSize_SR = localHeight * localWidth;

        _auxMemory.Require( recvSize_RS + max(sendSize_RS,recvSize_SR) );

        T* buffer = _auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[recvSize_RS];

        // Pack
        int* recvSizes = new int[r];
        for( int k=0; k<r; ++k )
        {
            T* data = &secondBuffer[k*recvSize_RS];
            recvSizes[k] = recvSize_RS;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColShift+i*r,j);
        }

        // Reduce-scatter over each process col
        ReduceScatter
        ( secondBuffer, firstBuffer, recvSizes, MPI_SUM, _grid->MCComm() );
        delete[] recvSizes;

        // Trade reduced data with the appropriate process col
        SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendCol, 0,
          secondBuffer, localHeight*localWidth,    recvCol, MPI_ANY_TAG,
          _grid->MRComm() );

        // Update with our received data
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) += alpha*secondBuffer[i+j*localHeight];

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,MR>::ReduceScatterUpdate
( T alpha, const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::ReduceScatterUpdate([* ,* ])");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        ResizeTo( A.Height(), A.Width() );
    }

    const int r = _grid->Height();
    const int c = _grid->Width();
    const int colAlignment = ColAlignment();
    const int rowAlignment = RowAlignment();

    const int height = Height();
    const int width = Width();
    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    const int maxLocalHeight = MaxLocalLength(height,r);
    const int maxLocalWidth = MaxLocalLength(width,c);

    const int recvSize = max(maxLocalHeight*maxLocalWidth,MinCollectContrib);
    const int sendSize = r * c * recvSize;

    _auxMemory.Require( sendSize + recvSize );

    T* buffer = _auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack 
    int* recvSizes = new int[r*c];
    for( int l=0; l<c; ++l )
    {
        const int thisRowShift = Shift( l, rowAlignment, c );
        const int thisLocalWidth = LocalLength( width, thisRowShift, c );

        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[(k+l*r)*recvSize];
            recvSizes[k+l*r] = recvSize;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColShift+i*r,thisRowShift+j*c);
        }
    }

    // Reduce-scatter over each process col
    ReduceScatter
    ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, _grid->VCComm() );
    delete[] recvSizes;

    // Unpack our received data
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) += alpha*recvBuffer[i+j*localHeight];

    _auxMemory.Release();

#ifndef RELEASE
    PopCallStack();
#endif
}

template class Elemental::DistMatrix<int,     MC,MR>;
template class Elemental::DistMatrix<float,   MC,MR>;
template class Elemental::DistMatrix<double,  MC,MR>;
#ifndef WITHOUT_COMPLEX
template class Elemental::DistMatrix<scomplex,MC,MR>;
template class Elemental::DistMatrix<dcomplex,MC,MR>;
#endif

