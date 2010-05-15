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
Elemental::DistMatrix<T,MR,MC>::Print( const string& msg ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::Print");
#endif
    if( _grid->VCRank() == 0 && msg != "" )
        cout << msg << endl;

    const int height = Height();
    const int width  = Width();
    const int localHeight = LocalHeight();
    const int localWidth  = LocalWidth();
    const int r = _grid->Height();
    const int c = _grid->Width();
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
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            sendBuf[colShift+i*c + (rowShift+j*r)*height] = _localMatrix(i,j);

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
Elemental::DistMatrix<T,MR,MC>::AlignWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignWith(DistMatrix[MR,MC])");
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
Elemental::DistMatrix<T,MR,MC>::AlignWith
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignWith(DistMatrix[MR,* ])");
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
Elemental::DistMatrix<T,MR,MC>::AlignWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignWith(DistMatrix[* ,MC])");
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
Elemental::DistMatrix<T,MR,MC>::AlignWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignWith(DistMatrix[MC,MR])");
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
Elemental::DistMatrix<T,MR,MC>::AlignWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignWith(DistMatrix[MC,*])");
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
Elemental::DistMatrix<T,MR,MC>::AlignWith
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignWith(DistMatrix[* ,MR])");
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
Elemental::DistMatrix<T,MR,MC>::AlignWith
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignWith(DistMatrix[VC,* ])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _rowShift     = Shift( _grid->MCRank(), _rowAlignment, _grid->Height() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::AlignWith
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]:AlignWith(DistMatrix[* ,VC])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _rowShift     = Shift( _grid->MCRank(), _rowAlignment, _grid->Height() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::AlignWith
( const DistMatrix<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignWith(DistMatrix[VR,* ])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment();
    _colShift     = Shift( _grid->MRRank(), _colAlignment, _grid->Width() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::AlignWith
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]:AlignWith(DistMatrix[* ,VR])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment();
    _colShift     = Shift( _grid->MRRank(), _colAlignment, _grid->Width() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::AlignColsWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignColsWith(DistMatrix[MR,MC])");
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
Elemental::DistMatrix<T,MR,MC>::AlignColsWith
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignColsWith(DistMatrix[MR,* ])");
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
Elemental::DistMatrix<T,MR,MC>::AlignColsWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignColsWith(DistMatrix[MC,MR])");
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
Elemental::DistMatrix<T,MR,MC>::AlignColsWith
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignColsWith(DistMatrix[* ,MR])");
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
Elemental::DistMatrix<T,MR,MC>::AlignColsWith
( const DistMatrix<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignColsWith(DistMatrix[VR,* ])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment();
    _colShift     = Shift( _grid->MRRank(), _colAlignment, _grid->Width() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::AlignColsWith
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]:AlignColsWith(DistMatrix[* ,VR])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment();
    _colShift     = Shift( _grid->MRRank(), _colAlignment, _grid->Width() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::AlignRowsWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignRowsWith(DistMatrix[MR,MC])");
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
Elemental::DistMatrix<T,MR,MC>::AlignRowsWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignRowsWith(DistMatrix[* ,MC])");
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
Elemental::DistMatrix<T,MR,MC>::AlignRowsWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignRowsWith(DistMatrix[MC,MR])");
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
Elemental::DistMatrix<T,MR,MC>::AlignRowsWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignRowsWith(DistMatrix[MC,* ])");
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
Elemental::DistMatrix<T,MR,MC>::AlignRowsWith
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::AlignRowsWith(DistMatrix[VC,* ])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _rowShift     = Shift( _grid->MCRank(), _rowAlignment, _grid->Height() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::AlignRowsWith
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]:AlignRowsWith(DistMatrix[* ,VC])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _rowShift     = Shift( _grid->MCRank(), _rowAlignment, _grid->Height() );
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::ConformWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::ConformWith(DistMatrix[MC,* ])");
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
Elemental::DistMatrix<T,MR,MC>::ConformWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::ConformWith(DistMatrix[* ,MC])");
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
Elemental::DistMatrix<T,MR,MC>::ConformWith
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::ConformWith(DistMatrix[MR,* ])");
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
Elemental::DistMatrix<T,MR,MC>::ConformWith
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::ConformWith(DistMatrix[* ,MR])");
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
Elemental::DistMatrix<T,MR,MC>::FreeConstraints()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::FreeConstraints");
#endif
    _constrainedColDist = false;
    _constrainedRowDist = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::View
( DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::View(A)");
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
Elemental::DistMatrix<T,MR,MC>::LockedView
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::LockedView(A)");
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
Elemental::DistMatrix<T,MR,MC>::View
( DistMatrix<T,MR,MC>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::View(A,i,j,height,width)");
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

        _colAlignment = (A.ColAlignment()+i) % c;
        _rowAlignment = (A.RowAlignment()+j) % r;
        
        _colShift = Shift( col, _colAlignment, c );
        _rowShift = Shift( row, _rowAlignment, r );

        const int localHeightBefore = LocalLength(i,A.ColShift(),c);
        const int localWidthBefore  = LocalLength(j,A.RowShift(),r);

        const int localHeight = LocalLength( height, _colShift, c );
        const int localWidth  = LocalLength( width,  _rowShift, r );

        _localMatrix.View( A.LocalMatrix(),
                           localHeightBefore, localWidthBefore, 
                           localHeight, localWidth             );
    }
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::LockedView
( const DistMatrix<T,MR,MC>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::LockedView(A,i,j,height,width)");
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

        _colAlignment = (A.ColAlignment()+i) % c;
        _rowAlignment = (A.RowAlignment()+j) % r;
        
        _colShift = Shift( col, _colAlignment, c );
        _rowShift = Shift( row, _rowAlignment, r );

        const int localHeightBefore = LocalLength(i,A.ColShift(),c);
        const int localWidthBefore  = LocalLength(j,A.RowShift(),r);

        const int localHeight = LocalLength( height, _colShift, c );
        const int localWidth  = LocalLength( width,  _rowShift, r );

        _localMatrix.LockedView
        ( A.LockedLocalMatrix(),
          localHeightBefore, localWidthBefore, localHeight, localWidth );
    }
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::View1x2
( DistMatrix<T,MR,MC>& AL,
  DistMatrix<T,MR,MC>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::View1x2");
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
Elemental::DistMatrix<T,MR,MC>::LockedView1x2
( const DistMatrix<T,MR,MC>& AL,
  const DistMatrix<T,MR,MC>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::LockedView1x2");
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
Elemental::DistMatrix<T,MR,MC>::View2x1
( DistMatrix<T,MR,MC>& AT,
  DistMatrix<T,MR,MC>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::View2x1");
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
Elemental::DistMatrix<T,MR,MC>::LockedView2x1
( const DistMatrix<T,MR,MC>& AT,
  const DistMatrix<T,MR,MC>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::LockedView2x1");
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
Elemental::DistMatrix<T,MR,MC>::View2x2
( DistMatrix<T,MR,MC>& ATL,
  DistMatrix<T,MR,MC>& ATR,
  DistMatrix<T,MR,MC>& ABL,
  DistMatrix<T,MR,MC>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::View2x2");
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
        ATR.RowAlignment() != ABR.RowAlignment()    )
    {
        throw "2x2 must align to combine.";
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
Elemental::DistMatrix<T,MR,MC>::LockedView2x2
( const DistMatrix<T,MR,MC>& ATL,
  const DistMatrix<T,MR,MC>& ATR,
  const DistMatrix<T,MR,MC>& ABL,
  const DistMatrix<T,MR,MC>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::LockedView2x2");
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
        ATR.RowAlignment() != ABR.RowAlignment()    )
    {
        throw "2x2 must align to combine.";
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
Elemental::DistMatrix<T,MR,MC>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::ResizeTo");
    CHECK_IF_LOCKED_VIEW;
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _height = height;
    _width  = width;
    _localMatrix.ResizeTo(LocalLength(height,_colShift,_grid->Width()),
                          LocalLength(width, _rowShift,_grid->Height()));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
Elemental::DistMatrix<T,MR,MC>::Get
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::Get");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const int ownerRow = (j + RowAlignment()) % _grid->Height();
    const int ownerCol = (i + ColAlignment()) % _grid->Width();
    const int ownerRank = ownerRow + ownerCol*_grid->Height();

    T u;
    if( _grid->VCRank() == ownerRank )
    {
        const int iLoc = (i-ColShift()) / _grid->Width();
        const int jLoc = (j-RowShift()) / _grid->Height();
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
Elemental::DistMatrix<T,MR,MC>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::Set");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const int ownerRow = (j + RowAlignment()) % _grid->Height();
    const int ownerCol = (i + ColAlignment()) % _grid->Width();
    const int ownerRank = ownerRow + ownerCol*_grid->Height();

    if( _grid->VCRank() == ownerRank )
    {
        const int iLoc = (i-ColShift()) / _grid->Width();
        const int jLoc = (j-RowShift()) / _grid->Height();
        _localMatrix(iLoc,jLoc) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::GetDiagonal
( DistMatrix<T,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::GetDiagonal");
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
    if( d.Viewing() && length != d.Height() )
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

        const int iLocStart = (iStart-colShift) / c;
        const int jLocStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalHeight();
        for( int kLoc=0; kLoc<localDiagLength; ++kLoc )
            d.LocalEntry(kLoc,0) = LocalEntry(iLocStart+kLoc*(lcm/c),
                                              jLocStart+kLoc*(lcm/r));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::GetDiagonal
( DistMatrix<T,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::GetDiagonal");
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
    if( d.Viewing() && length != d.Width() )
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

        const int iLocStart = (iStart-colShift) / c;
        const int jLocStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalWidth();
        for( int kLoc=0; kLoc<localDiagLength; ++kLoc )
            d.LocalEntry(0,kLoc) = LocalEntry(iLocStart+kLoc*(lcm/c),
                                              jLocStart+kLoc*(lcm/r));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::SetDiagonal
( const DistMatrix<T,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::SetDiagonal");
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

        const int iLocStart = (iStart-colShift) / c;
        const int jLocStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalHeight();
        for( int kLoc=0; kLoc<localDiagLength; ++kLoc )
            LocalEntry(iLocStart+kLoc*(lcm/c),
                       jLocStart+kLoc*(lcm/r)) = d.LocalEntry(kLoc,0);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::SetDiagonal
( const DistMatrix<T,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::SetDiagonal");
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

        const int iLocStart = (iStart-colShift) / c;
        const int jLocStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalWidth();
        for( int kLoc=0; kLoc<localDiagLength; ++kLoc )
            LocalEntry(iLocStart+kLoc*(lcm/c),
                       jLocStart+kLoc*(lcm/r)) = d.LocalEntry(0,kLoc);
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
Elemental::DistMatrix<T,MR,MC>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::MakeTrapezoidal");
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
            const int j = rowShift + jLoc*r;
            int lastZero_i;
            if( side == Left )
                lastZero_i = j-offset-1;
            else
                lastZero_i = j-offset+height-width-1;
            if( lastZero_i >= 0 )
            {
                const int boundary = min( lastZero_i+1, height );
                const int numZeros = LocalLength( boundary, colShift, c );
                for( int iLoc=0; iLoc<numZeros; ++iLoc )
                    _localMatrix(iLoc,jLoc) = (T)0;
            }
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
            const int nonzeroLength = LocalLength(firstZero_i,colShift,c);
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
Elemental::DistMatrix<T,MR,MC>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::SetToIdentity");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int localHeight = LocalHeight();
    const int localWidth  = LocalWidth();
    const int r = _grid->Height();
    const int c = _grid->Width();
    const int colShift = ColShift();
    const int rowShift = RowShift();

    _localMatrix.SetToZero();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*c;
        if( i % r == rowShift )
        {
            const int jLoc = (i-rowShift) / r;
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
Elemental::DistMatrix<T,MR,MC>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::SetToRandom");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int localHeight = LocalHeight();
    const int localWidth  = LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = Random<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::SetToRandomDiagDominant()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::SetToRandomDiagDominant");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int localHeight = LocalHeight();
    const int localWidth  = LocalWidth();
    const int r = _grid->Height();
    const int c = _grid->Width();
    const int colShift = ColShift();
    const int rowShift = RowShift();

    _localMatrix.SetToZero();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*c;
        if( i % r == rowShift )
        {
            const int jLoc = (i-rowShift) / r;
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
Elemental::DistMatrix<T,MR,MC>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::SetToZero");
    CHECK_IF_LOCKED_VIEW;
#endif
    _localMatrix.SetToZero();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[MC,MR]");
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
        const int ownerRow = RowAlignment();
        const int ownerCol = A.RowAlignment();
        const int colAlignment = ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = A.Height();
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = max(maxLocalHeight,MinCollectContrib);

        const int colShiftVR = Shift(rankRM,colAlignment,p);
        const int colShiftVCOfA = Shift(rankCM,colAlignmentOfA,p);
        const int sendRankRM = (rankRM+(p+colShiftVCOfA-colShiftVR)) % p;
        const int recvRankCM = (rankCM+(p+colShiftVR-colShiftVCOfA)) % p;
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

                const int shift = Shift(myRow+r*k,colAlignmentOfA,p);
                const int offset = (shift-A.ColShift()) / r;
                const int thisLocalHeight = LocalLength(height,shift,p);

                for( int i=0; i<thisLocalHeight; ++i )
                    data[i] = A.LocalEntry(offset+i*c,0);
            }
        }

        // A[VC,* ] <- A[MC,MR]
        Scatter
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerCol, _grid->MRComm() );

        // A[VR,* ] <- A[VC,* ]
        SendRecv
        ( sendBuf, portionSize, sendRankRM, 0,
          recvBuf, portionSize, recvRankRM, MPI_ANY_TAG, 
          _grid->VRComm() );

        // A[MR,MC] <- A[VR,* ]
        Gather
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerRow, _grid->MCComm() );

        if( myRow == ownerRow )
        {
            // Unpack
            for( int k=0; k<r; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const int shift = Shift(myCol+c*k,colAlignment,p);
                const int offset = (shift-ColShift()) / c;
                const int thisLocalHeight = LocalLength(height,shift,p);

                for( int i=0; i<thisLocalHeight; ++i )
                    _localMatrix(offset+i*r,0) = data[i];
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
        const int ownerCol = ColAlignment();
        const int ownerRow = A.ColAlignment();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int width = A.Width();
        const int maxLocalWidth = MaxLocalLength(width,p);

        const int portionSize = max(maxLocalWidth,MinCollectContrib);

        const int rowShiftVC = Shift(rankCM,rowAlignment,p);
        const int rowShiftVROfA = Shift(rankRM,rowAlignmentOfA,p);
        const int sendRankCM = (rankCM+(p+rowShiftVROfA-rowShiftVC)) % p;
        const int recvRankRM = (rankRM+(p+rowShiftVC-rowShiftVROfA)) % p;
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

                const int shift = Shift(myCol+c*k,rowAlignmentOfA,p);
                const int offset = (shift-A.RowShift()) / c;
                const int thisLocalWidth = LocalLength(width,shift,p);

                for( int j=0; j<thisLocalWidth; ++j )
                    data[j] = A.LocalEntry(0,offset+j*r);
            }
        }

        // A[* ,VR] <- A[MC,MR]
        Scatter
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerRow, _grid->MCComm() );

        // A[* ,VC] <- A[* ,VR]
        SendRecv
        ( sendBuf, portionSize, sendRankCM, 0,
          recvBuf, portionSize, recvRankCM, MPI_ANY_TAG, _grid->VCComm() );

        // A[MR,MC] <- A[* ,VC]
        Gather
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerCol, _grid->MRComm() );

        if( myCol == ownerCol )
        {
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
        }

        _auxMemory.Release();
    }
    else
    {
        if( A.Height() >= A.Width() )
        {
            auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
                ( new DistMatrix<T,VC,Star>(*_grid) );
            *A_VC_Star = A;

            auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
                ( new DistMatrix<T,VR,Star>(true,ColAlignment(),*_grid) );
            *A_VR_Star = *A_VC_Star;
            delete A_VC_Star.release(); // lowers memory highwater

            *this = *A_VR_Star;
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
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[MC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
        ( new DistMatrix<T,VC,Star>(*_grid) );
    *A_VC_Star = A;

    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
        ( new DistMatrix<T,VR,Star>(true,ColAlignment(),*_grid) );
    *A_VR_Star = *A_VC_Star;
    delete A_VC_Star.release(); // lowers memory highwater

    *this = *A_VR_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[* ,MR]");
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
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[MD,* ]");
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
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[* ,MD]");
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
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[MR,MC]");
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
            cout << "Unaligned [MR,MC] <- [MR,MC]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int row = _grid->MCRank();
        const int col = _grid->MRRank();

        const int colAlignment = ColAlignment();
        const int rowAlignment = RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;
        const int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;
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
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[MR,* ]");
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
                        ( _grid->MRRank(), _colAlignment, _grid->Width() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( ColAlignment() == A.ColAlignment() )
    {
        const int r = _grid->Height();
        const int rowShift = RowShift();

        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = A.LocalEntry(i,rowShift+j*r);
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MR,MC] <- [MR,* ]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int col = _grid->MRRank();

        const int rowShift = RowShift();
        const int colAlignment = ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (col+c+colAlignment-colAlignmentOfA) % c;
        const int recvRank = (col+c+colAlignmentOfA-colAlignment) % c;

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
                sendBuffer[i+j*localHeightOfA] = A.LocalEntry(i,rowShift+j*r);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, _grid->MRComm() );

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
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[* ,MC]");
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
        const int colShift = ColShift();

        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = A.LocalEntry(colShift+i*c,j);
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MR,MC] <- [* ,MC]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int row = _grid->MCRank();

        const int colShift = ColShift();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

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
                sendBuffer[i+j*localHeight] = A.LocalEntry(colShift+i*c,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRow, 0,
          recvBuffer, recvSize, recvRow, MPI_ANY_TAG, _grid->MCComm() );

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
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[VC,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,VR,Star> A_VR_Star(*_grid);

    A_VR_Star = A;
    *this     = A_VR_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[* ,VC]");
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
        const int p = r * c;
        const int row = _grid->MCRank();

        const int rowShift = RowShift();
        const int colAlignment = ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int height = Height();
        const int width = Width();
        const int localHeight = LocalHeight();
        const int localWidthOfA = A.LocalWidth();

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

            const int thisColShift = Shift(k,colAlignment,c);
            const int thisLocalHeight = LocalLength(height,thisColShift,c);

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColShift+i*c,j);
        }

        // Communicate
        AllToAll( sendBuffer, portionSize,
                  recvBuffer, portionSize, _grid->MRComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisRank = row+k*r;
            const int thisRowShift = Shift(thisRank,rowAlignmentOfA,p);
            const int thisRowOffset = (thisRowShift-rowShift) / r;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,thisRowOffset+j*c) = data[i+j*localHeight];
        }

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MR,MC] <- [* ,VC]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = r * c;
        const int row = _grid->MCRank();

        const int rowShift = RowShift();
        const int colAlignment = ColAlignment();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+rowAlignment-(rowAlignmentOfA%r)) % r;
        const int recvRow = (row+r+(rowAlignmentOfA%r)-rowAlignment) % r;

        const int height = Height();
        const int width = Width();
        const int localHeight = LocalHeight();
        const int localWidthOfA = A.LocalWidth();

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

            const int thisColShift = Shift(k,colAlignment,c);
            const int thisLocalHeight = LocalLength(height,thisColShift,c);

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColShift+i*c,j);
        }

        // SendRecv to align A[* ,VC] via a trade in the column
        SendRecv
        ( secondBuffer, c*portionSize, sendRow, 0,
          firstBuffer,  c*portionSize, recvRow, MPI_ANY_TAG, 
          _grid->MCComm() );

        // AllToAll to gather all of the aligned [* ,VC] into secondBuffer
        AllToAll
        ( firstBuffer,  portionSize, 
          secondBuffer, portionSize, _grid->MRComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisRank = recvRow+k*r;
            const int thisRowShift = Shift(thisRank,rowAlignmentOfA,p);
            const int thisRowOffset = (thisRowShift-rowShift) / r;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,thisRowOffset+j*c) = data[i+j*localHeight];
        }

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[VR,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        if( ! ConstrainedColDist() )
        {
            _colAlignment = A.ColAlignment() % _grid->Width();
            _colShift = Shift
                        ( _grid->MRRank(), _colAlignment, _grid->Width() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( ColAlignment() == A.ColAlignment() % _grid->Width() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = r * c;
        const int col = _grid->MRRank();

        const int colShift = ColShift();
        const int rowAlignment = RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = Height();
        const int width = Width();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,r);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        _auxMemory.Require( 2*r*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[r*portionSize];

        // Pack
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignment,r);
            const int thisLocalWidth = LocalLength(width,thisRowShift,r);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = 
                          A.LocalEntry(i,thisRowShift+j*r);
        }

        // Communicate
        AllToAll( sendBuffer, portionSize,
                  recvBuffer, portionSize, _grid->MCComm() );

        // Unpack
        for( int k=0; k<r; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisRank = col+k*c;
            const int thisColShift = Shift(thisRank,colAlignmentOfA,p);
            const int thisColOffset = (thisColShift-colShift) / c;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);
            
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    _localMatrix(thisColOffset+i*r,j) = 
                          data[i+j*thisLocalHeight];
        }

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MR,MC] <- [* ,VC]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = r * c;
        const int col = _grid->MRRank();

        const int colShift = ColShift();
        const int colAlignment = ColAlignment();
        const int rowAlignment = RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendCol = (col+c+colAlignment-(colAlignmentOfA%c)) % c;
        const int recvCol = (col+c+(colAlignmentOfA%c)-colAlignment) % c;

        const int height = Height();
        const int width = Width();
        const int localWidth = LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,r);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        _auxMemory.Require( 2*r*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[r*portionSize];

        // Pack
        for( int k=0; k<r; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignment,r);
            const int thisLocalWidth = LocalLength(width,thisRowShift,r);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = A.LocalEntry(i,thisRowShift+j*r);
        }

        // SendRecv to align A[VR,* ] via a trade in the row
        SendRecv
        ( secondBuffer, r*portionSize, sendCol, 0,
          firstBuffer,  r*portionSize, recvCol, MPI_ANY_TAG, 
          _grid->MRComm() );

        // AllToAll to gather all of the aligned [VR,* ] data into secondBuffer
        AllToAll
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, _grid->MCComm() );

        // Unpack
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisRank = recvCol+k*c;
            const int thisColShift = Shift(thisRank,colAlignmentOfA,p);
            const int thisColOffset = (thisColShift-colShift) / c;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    _localMatrix(thisColOffset+i*r,j) = 
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
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[* ,VR]");
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
const DistMatrix<T,MR,MC>&
Elemental::DistMatrix<T,MR,MC>::operator=
( const DistMatrix<T,Star,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC] = DistMatrix[* ,* ]");
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
            _localMatrix(i,j) = A.LocalEntry(colShift+i*c,rowShift+j*r);
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
void
Elemental::DistMatrix<T,MR,MC>::ReduceScatterFrom
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::ReduceScatterFrom([MR,* ])");
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
                        ( _grid->MRRank(), _colAlignment, _grid->Width() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( ColAlignment() == A.ColAlignment() )
    {
        if( Width() == 1 )
        {
            const int rowAlignment = RowAlignment();
            const int myRow = _grid->MCRank();

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
              MPI_SUM, rowAlignment, _grid->MCComm() );

            if( myRow == rowAlignment )
            {
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,0) = recvBuffer[i];
            }

            _auxMemory.Release();
        }
        else
        {
            const int r = _grid->Height();
            const int rowAlignment = RowAlignment();

            const int width = Width();
            const int localHeight = LocalHeight();
            const int localWidth = LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int recvSize = max(localHeight * maxLocalWidth,
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

                const int thisRowShift = Shift( k, rowAlignment, r );
                const int thisLocalWidth = 
                      LocalLength( width, thisRowShift, r );

                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        data[i+j*localHeight] = 
                            A.LocalEntry(i,thisRowShift+j*r);
            }

            // Reduce-scatter over each process column
            ReduceScatter
            ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, 
              _grid->MCComm() );
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
            cout << "Unaligned ReduceScatterFrom [MR,MC] <- [MR,* ]." << endl;
#endif
        if( Width() == 1 )
        {
            const int c = _grid->Width();
            const int rowAlignment = RowAlignment();
            const int myRow = _grid->MCRank();
            const int myCol = _grid->MRRank();

            const int height = Height();
            const int localHeight = LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,c);

            const int portionSize = max(maxLocalHeight,MinCollectContrib);

            const int colAlignment = ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (myCol+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (myCol+c+colAlignmentOfA-colAlignment) % c;

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
              MPI_SUM, rowAlignment, _grid->MCComm() );

            if( myRow == rowAlignment )
            {
                // Perform the realignment
                SendRecv
                ( recvBuffer, portionSize, sendCol, 0,
                  sendBuffer, portionSize, recvCol, 0, _grid->MRComm() );

                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,0) = sendBuffer[i];
            }

            _auxMemory.Release();
        }
        else
        {
            const int r = _grid->Height();
            const int c = _grid->Width();
            const int col = _grid->MRRank();

            const int colAlignment = ColAlignment();
            const int rowAlignment = RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

            const int width = Width();
            const int localHeight = LocalHeight();
            const int localWidth = LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int recvSize_RS = max(localHeightOfA * maxLocalWidth,
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

                const int thisRowShift = Shift( k, rowAlignment, r );
                const int thisLocalWidth = LocalLength(width,thisRowShift,r);

                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeightOfA; ++i )
                        data[i+j*localHeightOfA] = 
                            A.LocalEntry(i,thisRowShift+j*r);
            }

            // Reduce-scatter over each process col
            ReduceScatter
            ( secondBuffer, firstBuffer, recvSizes, MPI_SUM, 
              _grid->MCComm() );
            delete[] recvSizes;

            // Trade reduced data with the appropriate process col
            SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendCol, 0,
              secondBuffer, localHeight*localWidth,    recvCol, MPI_ANY_TAG,
              _grid->MRComm() );

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
Elemental::DistMatrix<T,MR,MC>::ReduceScatterFrom
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::ReduceScatterFrom([* ,MC])");
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
        const int colAlignment = ColAlignment();

        const int height = Height();
        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);

        const int recvSize = max(maxLocalHeight*localWidth,MinCollectContrib);
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

            const int thisColShift = Shift( k, colAlignment, c );
            const int thisLocalHeight = LocalLength( height, thisColShift, c );

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] =
                          A.LocalEntry(thisColShift+i*c,j);
        }

        // Reduce-scatter over each process row
        ReduceScatter
        ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, _grid->MRComm() );
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
            cout << "Unaligned ReduceScatterFrom [MR,MC] <- [* ,MC]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int row = _grid->MCRank();

        const int colAlignment = ColAlignment();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const int height = Height();
        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);
        
        const int recvSize_RS = max(maxLocalHeight * localWidthOfA,
                                    MinCollectContrib              );
        const int sendSize_RS = c* recvSize_RS;
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

            const int thisColShift = Shift( k, colAlignment, c );
            const int thisLocalHeight = LocalLength( height, thisColShift, c );

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColShift+i*c,j);
        }

        // Reduce-scatter over each process row
        ReduceScatter
        ( secondBuffer, firstBuffer, recvSizes, MPI_SUM, _grid->MRComm() );
        delete[] recvSizes;

        // Trade reduced data with the appropriate process row
        SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendRow, 0,
          secondBuffer, localHeight*localWidth,    recvRow, MPI_ANY_TAG, 
          _grid->MCComm() );

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
Elemental::DistMatrix<T,MR,MC>::ReduceScatterFrom
( const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::ReduceScatterFrom([* ,* ])");
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
    const int maxLocalHeight = MaxLocalLength(height,c);
    const int maxLocalWidth = MaxLocalLength(width,r);

    const int recvSize = max(maxLocalHeight*maxLocalWidth,MinCollectContrib);
    const int sendSize = r * c * recvSize;

    _auxMemory.Require( sendSize + recvSize );

    T* buffer = _auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack 
    int* recvSizes = new int[r*c];
    for( int l=0; l<r; ++l )
    {
        const int thisRowShift = Shift( l, rowAlignment, r );
        const int thisLocalWidth = LocalLength( width, thisRowShift, r );

        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[(k+l*c)*recvSize];
            recvSizes[k+l*c] = recvSize;

            const int thisColShift = Shift( k, colAlignment, c );
            const int thisLocalHeight = LocalLength( height, thisColShift, c );

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] =
                          A.LocalEntry(thisColShift+i*c,thisRowShift+j*r);
        }
    }

    // Reduce-scatter over each process col
    ReduceScatter
    ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, _grid->VRComm() );
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
Elemental::DistMatrix<T,MR,MC>::ReduceScatterUpdate
( T alpha, const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::ReduceScatterUpdate([MR,* ])");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_DIFF_SIZE( A );
    CHECK_IF_REDIST_DIFF_GRID( A );
#endif
    if( ColAlignment() == A.ColAlignment() )
    {
        if( Width() == 1 )
        {
            const int rowAlignment = RowAlignment();
            const int myRow = _grid->MCRank();

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
              MPI_SUM, rowAlignment, _grid->MCComm() );

            if( myRow == rowAlignment )
            {
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,0) += alpha*recvBuffer[i];
            }

            _auxMemory.Release();
        }
        else
        {
            const int r = _grid->Height();
            const int rowAlignment = RowAlignment();

            const int width = Width();
            const int localHeight = LocalHeight();
            const int localWidth = LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int portionSize = max(localHeight*maxLocalWidth,
                                        MinCollectContrib         );

            _auxMemory.Require( (r+1)*portionSize );

            T* buffer = _auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack 
            int* recvSizes = new int[r];
            for( int k=0; k<r; ++k )
            {
                T* data = &sendBuffer[k*portionSize];
                recvSizes[k] = portionSize;

                const int thisRowShift = Shift( k, rowAlignment, r );
                const int thisLocalWidth = LocalLength(width,thisRowShift,r);

                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        data[i+j*localHeight] = 
                            A.LocalEntry(i,thisRowShift+j*r);
            }

            // Reduce-scatter over each process column
            ReduceScatter
            ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, 
              _grid->MCComm() );
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
            cout << "Unaligned ReduceScatterUpdate [MR,MC] <- [MR,* ]." << endl;
#endif
        if( Width() == 1 )
        {
            const int c = _grid->Width();
            const int rowAlignment = RowAlignment();
            const int myRow = _grid->MCRank();
            const int myCol = _grid->MRRank();

            const int height = Height();
            const int localHeight = LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,c);

            const int portionSize = max(maxLocalHeight,MinCollectContrib);

            const int colAlignment = ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (myCol+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (myCol+c+colAlignmentOfA-colAlignment) % c;

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
              MPI_SUM, rowAlignment, _grid->MCComm() );

            if( myRow == rowAlignment )
            {
                // Perform the realignment
                SendRecv
                ( recvBuffer, portionSize, sendCol, 0,
                  sendBuffer, portionSize, recvCol, 0, _grid->MRComm() );

                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,0) += alpha*sendBuffer[i];
            }

            _auxMemory.Release();
        }
        else
        {
            const int r = _grid->Height();
            const int c = _grid->Width();
            const int col = _grid->MRRank();

            const int colAlignment = ColAlignment();
            const int rowAlignment = RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

            const int width = Width();
            const int localHeight = LocalHeight();
            const int localWidth = LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int recvSize_RS = max(localHeightOfA * maxLocalWidth,
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

                const int thisRowShift = Shift( k, rowAlignment, r );
                const int thisLocalWidth = LocalLength(width,thisRowShift,r);

                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeightOfA; ++i )
                        data[i+j*localHeightOfA] =
                              A.LocalEntry(i,thisRowShift+j*r);
            }

            // Reduce-scatter over each process col
            ReduceScatter
            ( secondBuffer, firstBuffer, recvSizes, MPI_SUM, 
              _grid->MCComm() );
            delete[] recvSizes;

            // Trade reduced data with the appropriate process col
            SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendCol, 0, 
              secondBuffer, localHeight*localWidth,    recvCol, MPI_ANY_TAG,
              _grid->MRComm() );

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
Elemental::DistMatrix<T,MR,MC>::ReduceScatterUpdate
( T alpha, const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::ReduceScatterUpdate([* ,MC])");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_DIFF_SIZE( A );
    CHECK_IF_REDIST_DIFF_GRID( A );
#endif
    if( RowAlignment() == A.RowAlignment() )
    {
        const int c = _grid->Width();
        const int colAlignment = ColAlignment();

        const int height = Height();
        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);

        const int recvSize = max(maxLocalHeight*localWidth,MinCollectContrib);
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

            const int thisColShift = Shift( k, colAlignment, c );
            const int thisLocalHeight = LocalLength( height, thisColShift, c );

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] =
                          A.LocalEntry(thisColShift+i*c,j);
        }

        // Reduce-scatter over each process row
        ReduceScatter
        ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, _grid->MRComm() );
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
            cout << "Unaligned ReduceScatterUpdate [MR,MC] <- [* ,MC]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int row = _grid->MCRank();

        const int colAlignment = ColAlignment();
        const int rowAlignment = RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const int height = Height();
        const int localHeight = LocalHeight();
        const int localWidth = LocalWidth();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);

        const int recvSize_RS = max(maxLocalHeight * localWidthOfA,
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

            const int thisColShift = Shift( k, colAlignment, c );
            const int thisLocalHeight = LocalLength( height, thisColShift, c );

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] =
                          A.LocalEntry(thisColShift+i*c,j);
        }

        // Reduce-scatter over each process row
        ReduceScatter
        ( secondBuffer, firstBuffer, recvSizes, MPI_SUM, _grid->MRComm() );
        delete[] recvSizes;

        // Trade reduced data with the appropriate process row
        SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendRow, 0,
          secondBuffer, localHeight*localWidth,    recvRow, MPI_ANY_TAG,
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
Elemental::DistMatrix<T,MR,MC>::ReduceScatterUpdate
( T alpha, const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::ReduceScatterUpdate([* ,* ])");
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
    const int maxLocalHeight = MaxLocalLength(height,c);
    const int maxLocalWidth = MaxLocalLength(width,r);

    const int recvSize = max(maxLocalHeight*maxLocalWidth,MinCollectContrib);
    const int sendSize = r * c * recvSize;

    _auxMemory.Require( sendSize + recvSize );

    T* buffer = _auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack 
    int* recvSizes = new int[r*c];
    for( int l=0; l<r; ++l )
    {
        const int thisRowShift = Shift( l, rowAlignment, r );
        const int thisLocalWidth = LocalLength( width, thisRowShift, r );

        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[(k+l*c)*recvSize];
            recvSizes[k+l*c] = recvSize;

            const int thisColShift = Shift( k, colAlignment, c );
            const int thisLocalHeight = LocalLength( height, thisColShift, c );

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] =
                          A.LocalEntry(thisColShift+i*c,thisRowShift+j*r);
        }
    }

    // Reduce-scatter over each process col
    ReduceScatter
    ( sendBuffer, recvBuffer, recvSizes, MPI_SUM, _grid->VRComm() );
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

template class Elemental::DistMatrix<int,    MR,MC>;
template class Elemental::DistMatrix<float,  MR,MC>;
template class Elemental::DistMatrix<double, MR,MC>;
#ifndef WITHOUT_COMPLEX
template class Elemental::DistMatrix<scomplex,MR,MC>;
template class Elemental::DistMatrix<dcomplex,MR,MC>;
#endif

