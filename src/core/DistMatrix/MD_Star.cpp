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
Elemental::DistMatrix<T,MD,Star>::Print( const string& msg ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::Print");
#endif
    if( _grid->VCRank() == 0 && msg != "" )
        cout << msg << endl;
        
    const int height      = Height();
    const int width       = Width();
    const int localHeight = LocalHeight();
    const int lcm         = _grid->LCM();
    const int inDiagonal  = InDiagonal();

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
        for( int i=0; i<localHeight; ++i )
            for( int j=0; j<width; ++j )
                sendBuf[colShift+i*lcm+j*height] = _localMatrix(i,j);
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
Elemental::DistMatrix<T,MD,Star>::AlignWith
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::AlignWith(DistMatrix[MD,* ])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
        _colShift = A.ColShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::AlignWith
( const DistMatrix<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::AlignWith(DistMatrix[* ,MD])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
        _colShift = A.RowShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::AlignColsWith
( const DistMatrix<T,MD,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::AlignColsWith
( const DistMatrix<T,Star,MD>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::AlignWithDiag
( const DistMatrix<T,MC,MR>& A, int offset )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::AlignWithDiag");
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
        _colAlignment = ownerRow + r*ownerCol;
        _inDiagonal = ( _grid->DiagPath() == _grid->DiagPath( _colAlignment ) );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        _colAlignment = ownerRow + r*ownerCol;
        _inDiagonal = ( _grid->DiagPath() == _grid->DiagPath( _colAlignment ) );
    }
    if( _inDiagonal )
    {
        _colShift = (_grid->DiagPathRank()+lcm-
                     _grid->DiagPathRank(_colAlignment)) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::AlignWithDiag
( const DistMatrix<T,MR,MC>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::AlignWithDiag");
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
        _colAlignment = ownerRow + r*ownerCol;
        _inDiagonal = ( _grid->DiagPath() == _grid->DiagPath( _colAlignment ) );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        _colAlignment = ownerRow + r*ownerCol;
        _inDiagonal = ( _grid->DiagPath() == _grid->DiagPath( _colAlignment ) );
    }
    if( _inDiagonal )
    {
        _colShift = (_grid->DiagPathRank()+lcm-
                     _grid->DiagPathRank(_colAlignment)) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::ConformWith
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::ConformWith(DistMatrix[MD,* ])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_CONFORMING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
        _colShift = A.ColShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::ConformWith
( const DistMatrix<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::ConformWith(DistMatrix[* ,MD])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_CONFORMING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
        _colShift = A.RowShift();
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::FreeConstraints()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::FreeConstraints");
#endif
    _constrainedColDist = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::View
( DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::View(A)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _colAlignment = A.ColAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
    {
        _colShift = A.ColShift();
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
Elemental::DistMatrix<T,MD,Star>::LockedView
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::LockedView(A)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _colAlignment = A.ColAlignment();
    _inDiagonal   = A.InDiagonal();
    if( _inDiagonal )
    {
        _colShift = A.ColShift();
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
Elemental::DistMatrix<T,MD,Star>::View
( DistMatrix<T,MD,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::View(A,i,j,height,width)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
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
        const int alignmentRank = A.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        _colAlignment = newAlignmentRank;
        _inDiagonal   = A.InDiagonal();

        if( _inDiagonal )
        {
            _colShift = Shift
                        ( diagPathRank, _grid->DiagPathRank(_colAlignment), lcm );
            int localHeightBefore = LocalLength( i, A.ColShift(), lcm);
            int localHeight = LocalLength( height, _colShift, lcm );

            _localMatrix.View
            ( A.LocalMatrix(), localHeightBefore, j, localHeight, width );
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
Elemental::DistMatrix<T,MD,Star>::LockedView
( const DistMatrix<T,MD,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::LockedView(A,i,j,height,width)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
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
        const int alignmentRank = A.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        _colAlignment = newAlignmentRank;
        _inDiagonal   = A.InDiagonal();

        if( _inDiagonal )
        {
            _colShift = Shift
                        ( diagPathRank, _grid->DiagPathRank(_colAlignment), lcm );
            int localHeightBefore = LocalLength( i, A.ColShift(), lcm);
            int localHeight = LocalLength( height, _colShift, lcm );
        
            _localMatrix.LockedView
            ( A.LockedLocalMatrix(), localHeightBefore, j, localHeight, width );
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
Elemental::DistMatrix<T,MD,Star>::View1x2
( DistMatrix<T,MD,Star>& AL,
  DistMatrix<T,MD,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::View1x2");    
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
    if( AL.ColAlignment() != AR.ColAlignment() )
        throw "Cannot combine misaligned 1x2 array of matrices.";
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _colAlignment = AL.ColAlignment();
    _inDiagonal   = AL.InDiagonal();
    if( _inDiagonal )
    {
        _colShift = AL.ColShift();
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
Elemental::DistMatrix<T,MD,Star>::LockedView1x2
( const DistMatrix<T,MD,Star>& AL,
  const DistMatrix<T,MD,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::LockedView1x2");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
    if( AL.ColAlignment() != AR.ColAlignment() )
        throw "Cannot combine misaligned 1x2 array of matrices.";
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _colAlignment = AL.ColAlignment();
    _inDiagonal   = AL.InDiagonal();
    if( _inDiagonal )
    {
        _colShift = AL.ColShift();
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
Elemental::DistMatrix<T,MD,Star>::View2x1
( DistMatrix<T,MD,Star>& AT,
  DistMatrix<T,MD,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::View2x1");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _colAlignment = AT.ColAlignment();
    _inDiagonal   = AT.InDiagonal();
    if( _inDiagonal )
    {
        _colShift = AT.ColShift();
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
Elemental::DistMatrix<T,MD,Star>::LockedView2x1
( const DistMatrix<T,MD,Star>& AT,
  const DistMatrix<T,MD,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::LockedView2x1");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _colAlignment = AT.ColAlignment();
    _inDiagonal   = AT.InDiagonal();
    if( _inDiagonal )
    {
        _colShift = AT.ColShift();
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
Elemental::DistMatrix<T,MD,Star>::View2x2
( DistMatrix<T,MD,Star>& ATL,
  DistMatrix<T,MD,Star>& ATR,
  DistMatrix<T,MD,Star>& ABL,
  DistMatrix<T,MD,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::View2x2");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
    if( ATL.ColAlignment() != ATR.ColAlignment() ||
        ABL.ColAlignment() != ABR.ColAlignment()   )
    {
        throw "Cannot combine misaligned 2x2 grid of matrices.";
    }
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _colAlignment = ATL.ColAlignment();
    _inDiagonal   = ATL.InDiagonal();
    if( _inDiagonal )
    {
        _colShift = ATL.ColShift();
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
Elemental::DistMatrix<T,MD,Star>::LockedView2x2
( const DistMatrix<T,MD,Star>& ATL,
  const DistMatrix<T,MD,Star>& ATR,
  const DistMatrix<T,MD,Star>& ABL,
  const DistMatrix<T,MD,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::LockedView2x2");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
    if( ATL.ColAlignment() != ATR.ColAlignment() ||
        ABL.ColAlignment() != ABR.ColAlignment()   )
    {
        throw "Cannot combine misaligned 2x2 grid of matrices.";
    }
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _colAlignment = ATL.ColAlignment();
    _inDiagonal   = ATL.InDiagonal();
    if( _inDiagonal )
    {
        _colShift = ATL.ColShift();
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
Elemental::DistMatrix<T,MD,Star>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::ResizeTo");
    CHECK_IF_LOCKED_VIEW;
    if( height < 0 || width < 0 )
    {
        if( _grid->VCRank() == 0 )
            cerr << "Height and width must be non-negative." << endl;
        DumpCallStack();
    }
#endif
    _height = height;
    _width  = width;
    if( _inDiagonal )
    {
        const int lcm = _grid->LCM();
        _localMatrix.ResizeTo
        ( LocalLength(height,_colShift,lcm), width );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
Elemental::DistMatrix<T,MD,Star>::Get
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::Get");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int alignmentRank = ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    T u;
    if( _grid->VCRank() == ownerRank )
    {
        const int iLoc = (i-ColShift()) / _grid->LCM();
        u = _localMatrix(iLoc,j);
    }
    Broadcast( &u, 1, ownerRank, _grid->VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::Set");
    if( i < 0 || i >= Height() || j < 0 || j >=Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    int ownerRank;
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int alignmentRank = ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( _grid->VCRank() == ownerRank )
    {
        const int iLoc = (i-ColShift()) / _grid->LCM();
        _localMatrix(iLoc,j) = u;
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
Elemental::DistMatrix<T,MD,Star>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::MakeTrapezoidal");
    CHECK_IF_LOCKED_VIEW;
#endif
    if( InDiagonal() )
    {
        const int height = Height();
        const int width = Width();
        const int localHeight = LocalHeight();
        const int lcm = _grid->LCM();
        const int colShift = ColShift();

        if( shape == Lower )
        {
            for( int j=0; j<width; ++j )
            {
                int lastZero_i;
                if( side == Left )
                    lastZero_i = j-offset-1;
                else
                    lastZero_i = j-offset+height-width-1;
                if( lastZero_i >= 0 )
                {
                    const int boundary = min( lastZero_i+1, height );
                    const int numZeros = LocalLength( boundary, colShift, lcm );
                    for( int iLoc=0; iLoc<numZeros; ++iLoc )
                        _localMatrix(iLoc,j) = (T)0;
                }
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
                const int nonzeroLength = LocalLength(firstZero_i,colShift,lcm);
                for( int iLoc=nonzeroLength; iLoc<localHeight; ++iLoc )
                    _localMatrix(iLoc,j) = (T)0;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::SetToIdentity");
    CHECK_IF_LOCKED_VIEW;
#endif
    if( InDiagonal() )
    {
        const int lcm = _grid->LCM();
        const int width = Width();
        const int localHeight = LocalHeight();
        const int colShift = ColShift();

        _localMatrix.SetToZero();
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*lcm;
            if( i < width )
                _localMatrix(iLoc,i) = (T)1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::SetToRandom");
    CHECK_IF_LOCKED_VIEW;
#endif
    if( InDiagonal() )
    {
        const int width       = Width();
        const int localHeight = LocalHeight();
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = Random<T>();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::SetToRandomDiagDominant()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::SetToRandomDiagDominant");
    CHECK_IF_LOCKED_VIEW;
#endif
    SetToRandom();

    if( InDiagonal() )
    {
        const int width = Width();
        const int localHeight = LocalHeight();
        const int lcm = _grid->LCM();
        const int colShift = ColShift();

        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*lcm;
            if( i < width )
                _localMatrix(iLoc,i) += (T)max(Height(),Width());
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MD,Star>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::SetToZero");
#endif
    if( InDiagonal() )
        _localMatrix.SetToZero();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[MC,MR]");
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
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[MC,* ]");
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
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[* ,MR]");
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
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[MD,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
    {
        if( ! ConstrainedColDist() )
        {
            _colAlignment = A.ColAlignment();
            _inDiagonal   = A.InDiagonal();
            if( _inDiagonal )
                _colShift = A.ColShift();
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( ColAlignment() == A.ColAlignment() )
    {
        _localMatrix = A.LockedLocalMatrix();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MD,* ] <- [MD,* ]." << endl;
#endif
        REPORT_UNIMPLEMENTED_FEATURE;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[* ,MD]");
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
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[MR,MC]");
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
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[MR,* ]");
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
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[* ,MC]");
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
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[VC,* ]");
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
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[* ,VC]");
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
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[VR,* ]");
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
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[* ,VR]");
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
const DistMatrix<T,MD,Star>&
Elemental::DistMatrix<T,MD,Star>::operator=
( const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ] = DistMatrix[* ,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    if( InDiagonal() )
    {
        const int lcm = _grid->LCM();
        const int colShift = ColShift();

        const int width = Width();
        const int localHeight = LocalHeight();
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = A.LocalEntry(colShift+i*lcm,j);
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class Elemental::DistMatrix<int,     MD,Star>;
template class Elemental::DistMatrix<float,   MD,Star>;
template class Elemental::DistMatrix<double,  MD,Star>;
#ifndef WITHOUT_COMPLEX
template class Elemental::DistMatrix<scomplex,MD,Star>;
template class Elemental::DistMatrix<dcomplex,MD,Star>;
#endif

