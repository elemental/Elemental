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
Elemental::DistMatrix<T,Star,MD>::Print( const string msg ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::Print");
#endif
    if( _grid->VCRank() == 0 && msg != "" )
        cout << msg << endl;
        
    const int height     = Height();
    const int width      = Width();
    const int localWidth = LocalWidth();
    const int lcm        = _grid->LCM();
    const int inDiagonal = InDiagonal();

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
    if( inDiagonal )
    {
        const int colShift = ColShift();
        for( int i=0; i<height; ++i )
            for( int j=0; j<localWidth; ++j )
                sendBuf[colShift+i+j*lcm*height] = _localMatrix(i,j);
    }

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

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::AlignWith
( const DistMatrix<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::AlignWith(DistMatrix[* ,MD])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
        _rowShift = A.RowShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::AlignWith
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::AlignWith(DistMatrix[MD,* ])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
        _rowShift = A.ColShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::AlignRowsWith
( const DistMatrix<T,Star,MD>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::AlignRowsWith
( const DistMatrix<T,MD,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::AlignWithDiag
( const DistMatrix<T,MC,MR>& A, const int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::AlignWithDiag");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    const int r = _grid->Height();
    const int c = _grid->Width();
    const int lcm = _grid->LCM();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = colAlignment;
        const int ownerCol = (rowAlignment + offset) % c;
        _rowAlignment = ownerRow + r*ownerCol;
        _inDiagonal = ( _grid->DiagPath() == _grid->DiagPath( _rowAlignment ) );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        _rowAlignment = ownerRow + r*ownerCol;
        _inDiagonal = ( _grid->DiagPath() == _grid->DiagPath( _rowAlignment ) );
    }
    if( _inDiagonal )
    {
        _rowShift = (_grid->DiagPathRank()+lcm-
                     _grid->DiagPathRank(_rowAlignment)) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::AlignWithDiag
( const DistMatrix<T,MR,MC>& A, const int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::AlignWithDiag");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    const int r = _grid->Height();
    const int c = _grid->Width();
    const int lcm = _grid->LCM();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = rowAlignment;
        const int ownerCol = (colAlignment + offset) % c;
        _rowAlignment = ownerRow + r*ownerCol;
        _inDiagonal = ( _grid->DiagPath() == _grid->DiagPath( _rowAlignment ) );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        _rowAlignment = ownerRow + r*ownerCol;
        _inDiagonal = ( _grid->DiagPath() == _grid->DiagPath( _rowAlignment ) );
    }
    if( _inDiagonal )
    {
        _rowShift = (_grid->DiagPathRank()+lcm-
                     _grid->DiagPathRank(_rowAlignment)) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::ConformWith
( const DistMatrix<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::ConformWith(DistMatrix[* ,MD])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_CONFORMING_DIFF_GRID( A );
#endif
    _rowAlignment = A.RowAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
        _rowShift = A.RowShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::ConformWith
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::ConformWith(DistMatrix[MD,* ])");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_CONFORMING_DIFF_GRID( A );
#endif
    _rowAlignment = A.ColAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
        _rowShift = A.ColShift();
    _constrainedRowDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::FreeConstraints()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::FreeConstraints");
#endif
    _constrainedRowDist = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::View
( DistMatrix<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::View(A)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _rowAlignment = A.RowAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
    {
        _rowShift = A.RowShift();
        _localMatrix.View( A.LocalMatrix() );
    }
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::LockedView
( const DistMatrix<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::LockedView(A)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _rowAlignment = A.RowAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
    {
        _rowShift = A.RowShift();
        _localMatrix.LockedView( A.LockedLocalMatrix() );
    }
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::View
( DistMatrix<T,Star,MD>& A,
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::View(A,i,j,height,width)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int lcm = _grid->LCM();
        const int diagPathRank = _grid->DiagPathRank(); 
        const int alignmentRank = A.RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        _rowAlignment = newAlignmentRank;
        _inDiagonal   = A.InDiagonal();

        if( _inDiagonal )
        {
            _rowShift = Shift
                        ( diagPathRank, _grid->DiagPathRank(_rowAlignment), lcm );
            int localWidthBefore = LocalLength( j, A.RowShift(), lcm );
            int localWidth = LocalLength( width, _rowShift, lcm );
        
            _localMatrix.View( A.LocalMatrix(),
                               i, localWidthBefore, height, localWidth );
        }

    }
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::LockedView
( const DistMatrix<T,Star,MD>& A,
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::LockedView(A,i,j,height,width)");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int lcm = _grid->LCM();
        const int diagPathRank = _grid->DiagPathRank();
        const int alignmentRank = A.RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        _rowAlignment = newAlignmentRank;
        _inDiagonal   = A.InDiagonal();

        if( _inDiagonal )
        {
            _rowShift = Shift
                        ( diagPathRank, _grid->DiagPathRank(_rowAlignment), lcm );
            int localWidthBefore = LocalLength( j, A.RowShift(), lcm);
            int localWidth = LocalLength( width, _rowShift, lcm );
        
            _localMatrix.LockedView
            ( A.LockedLocalMatrix(), i, localWidthBefore, height, localWidth );
        }
    }
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::View1x2
( DistMatrix<T,Star,MD>& AL,
  DistMatrix<T,Star,MD>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::View1x2");    
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _rowAlignment = AL.RowAlignment();
    _inDiagonal   = AL.InDiagonal();
    if( _inDiagonal )
    {
        _rowShift = AL.RowShift();
        _localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    }
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::LockedView1x2
( const DistMatrix<T,Star,MD>& AL,
  const DistMatrix<T,Star,MD>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::LockedView1x2");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _rowAlignment = AL.RowAlignment();
    _inDiagonal   = AL.InDiagonal();
    if( _inDiagonal )
    {
        _rowShift = AL.RowShift();
        _localMatrix.LockedView1x2( AL.LockedLocalMatrix(), 
                                    AR.LockedLocalMatrix() );
    }
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::View2x1
( DistMatrix<T,Star,MD>& AT,
  DistMatrix<T,Star,MD>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::View2x1");
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
    _inDiagonal   = AT.InDiagonal();
    if( _inDiagonal )
    {
        _rowShift = AT.RowShift();
        _localMatrix.View2x1( AT.LocalMatrix(),
                              AB.LocalMatrix() );
    }
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::LockedView2x1
( const DistMatrix<T,Star,MD>& AT,
  const DistMatrix<T,Star,MD>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::LockedView2x1");
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
    _inDiagonal   = AT.InDiagonal();
    if( _inDiagonal )
    {
        _rowShift = AT.RowShift();
        _localMatrix.LockedView2x1( AT.LockedLocalMatrix(),
                                    AB.LockedLocalMatrix() );
    }
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::View2x2
( DistMatrix<T,Star,MD>& ATL,
  DistMatrix<T,Star,MD>& ATR,
  DistMatrix<T,Star,MD>& ABL,
  DistMatrix<T,Star,MD>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::View2x2");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
    if( ATL.RowAlignment() != ABL.RowAlignment() ||
        ATR.RowAlignment() != ABR.RowAlignment()   )
    {
        throw "Cannot combine misaligned 2x2 grid of matrices.";
    }
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _rowAlignment = ATL.RowAlignment();
    _inDiagonal   = ATL.InDiagonal();
    if( _inDiagonal )
    {
        _rowShift = ATL.RowShift();
        _localMatrix.View2x2( ATL.LocalMatrix(), ATR.LocalMatrix(),
                              ABL.LocalMatrix(), ABR.LocalMatrix() );
    }
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::LockedView2x2
( const DistMatrix<T,Star,MD>& ATL,
  const DistMatrix<T,Star,MD>& ATR,
  const DistMatrix<T,Star,MD>& ABL,
  const DistMatrix<T,Star,MD>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::LockedView2x2");
    CHECK_IF_UNFREED_ROW_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
    if( ATL.RowAlignment() != ABL.RowAlignment() ||
        ATR.RowAlignment() != ABR.RowAlignment()   )
    {
        throw "Cannot combine misaligned 2x2 grid of matrices.";
    }
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _rowAlignment = ATL.RowAlignment();
    _inDiagonal   = ATL.InDiagonal();
    if( _inDiagonal )
    {
        _rowShift = ATL.RowShift();
        _localMatrix.LockedView2x2
        ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
          ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    }
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::ResizeTo
( const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::ResizeTo");
    CHECK_IF_LOCKED_VIEW;
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _height = height;
    _width  = width;
    if( _inDiagonal )
    {
        const int lcm = _grid->LCM();
        _localMatrix.ResizeTo
        ( height, LocalLength(width,_rowShift,lcm) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
Elemental::DistMatrix<T,Star,MD>::Get
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::Get");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        throw msg.str();
    }
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int alignmentRank = RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    T u;
    if( _grid->VCRank() == ownerRank )
    {
        const int jLoc = (j-RowShift()) / _grid->LCM();
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
Elemental::DistMatrix<T,Star,MD>::Set
( const int i, const int j, const T u )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::Set");
    if( i < 0 || i >= Height() || j < 0 || j >=Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        throw msg.str();
    }
#endif
    int ownerRank;
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int alignmentRank = RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( _grid->VCRank() == ownerRank )
    {
        const int jLoc = (j-RowShift()) / _grid->LCM();
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
Elemental::DistMatrix<T,Star,MD>::MakeTrapezoidal
( const Side side, const Shape shape, const int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::MakeTrapezoidal");
    CHECK_IF_LOCKED_VIEW;
#endif
    if( InDiagonal() )
    {
        const int height = Height();
        const int width = Width();
        const int localWidth = LocalWidth();
        const int lcm = _grid->LCM();
        const int rowShift = RowShift();

        if( shape == Lower )
        {
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*lcm;
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
                const int j = rowShift + jLoc*lcm;
                int firstZero_i;
                if( side == Left )
                    firstZero_i = max(j-offset+1,0);
                else
                    firstZero_i = max(j-offset+height-width+1,0);
                for( int i=firstZero_i; i<height; ++i )
                    _localMatrix(i,jLoc) = (T)0;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star, MD>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::SetToIdentity");
    CHECK_IF_LOCKED_VIEW;
#endif
    if( InDiagonal() )
    {
        const int lcm = _grid->LCM();
        const int height = Height();
        const int localWidth = LocalWidth();
        const int rowShift = RowShift();

        _localMatrix.SetToZero();
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*lcm;
            if( j < height )
                _localMatrix(j,jLoc) = (T)1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::SetToRandom");
    CHECK_IF_LOCKED_VIEW;
#endif
    if( InDiagonal() )
    {
        const int height = Height();
        const int localWidth = LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                _localMatrix(i,j) = Random<T>();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::SetToRandomDiagDominant()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::SetToRandomDiagDominant");
    CHECK_IF_LOCKED_VIEW;
#endif
    SetToRandom();

    if( InDiagonal() )
    {
        const int height = Height();
        const int localWidth = LocalWidth();
        const int lcm = _grid->LCM();
        const int rowShift = RowShift();

        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*lcm;
            if( j < height )
                _localMatrix(j,jLoc) += (T)max(Height(),Width());
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,Star,MD>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::SetToZero");
#endif
    if( InDiagonal() )
        _localMatrix.SetToZero();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[MC,MR]");
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
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[MC,* ]");
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
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[* ,MR]");
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
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[MD,* ]");
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
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[* ,MD]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        if( ! ConstrainedRowDist() )
        {
            _rowAlignment = A.RowAlignment();
            _inDiagonal   = A.InDiagonal();
            if( _inDiagonal )
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
            cout << "Unaligned [* ,MD] <- [* ,MD]." << endl;
#endif
        REPORT_UNIMPLEMENTED_FEATURE;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[MR,MC]");
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
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[MR,* ]");
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
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[* ,MC]");
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
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[VC,* ]");
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
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[* ,VC]");
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
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[VR,* ]");
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
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[* ,VR]");
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
const DistMatrix<T,Star,MD>&
Elemental::DistMatrix<T,Star,MD>::operator=
( const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD] = DistMatrix[* ,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    if( InDiagonal() )
    {
        const int lcm = _grid->LCM();
        const int rowShift = RowShift();

        const int height = Height();
        const int localWidth = LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                _localMatrix(i,j) = A.LocalEntry(i,rowShift+j*lcm);
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class Elemental::DistMatrix<int,     Star,MD>;
template class Elemental::DistMatrix<float,   Star,MD>;
template class Elemental::DistMatrix<double,  Star,MD>;
#ifndef WITHOUT_COMPLEX
template class Elemental::DistMatrix<scomplex,Star,MD>;
template class Elemental::DistMatrix<dcomplex,Star,MD>;
#endif

