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
Elemental::DistMatrix<T,MC,Star>::Print( const string msg ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::Print");
#endif
    if( _grid->VCRank() == 0 && msg != "" )
        cout << msg << endl;
        
    const int height      = Height();
    const int width       = Width();
    const int localHeight = LocalHeight();
    const int r           = _grid->Height();
    const int colShift    = ColShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Only one process column needs to participate
    if( _grid->MRRank() == 0 )
    {
        T* sendBuf = new T[height*width];
        for( int i=0; i<height*width; ++i )
            sendBuf[i] = (T)0;
        for( int i=0; i<localHeight; ++i )
            for( int j=0; j<width; ++j )
                sendBuf[colShift+i*r+j*height] = _localMatrix(i,j);

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
Elemental::DistMatrix<T,MC,Star>::AlignWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::AlignWith(DistMatrix[MC,MR])");
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
Elemental::DistMatrix<T,MC,Star>::AlignWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::AlignWith(DistMatrix[MC,* ])");
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
Elemental::DistMatrix<T,MC,Star>::AlignWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::AlignWith(DistMatrix[MR,MC])");
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
Elemental::DistMatrix<T,MC,Star>::AlignWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::AlignWith(DistMatrix[* ,MC])");
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
Elemental::DistMatrix<T,MC,Star>::AlignWith
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::AlignWith(DistMatrix[VC,* ])");
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
Elemental::DistMatrix<T,MC,Star>::AlignWith
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::AlignWith(DistMatrix[* ,VC])");
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
Elemental::DistMatrix<T,MC,Star>::AlignColsWith
( const DistMatrix<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::AlignColsWith
( const DistMatrix<T,MC,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::AlignColsWith
( const DistMatrix<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::AlignColsWith
( const DistMatrix<T,Star,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::AlignColsWith
( const DistMatrix<T,VC,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::AlignColsWith
( const DistMatrix<T,Star,VC>& A )
{ AlignWith( A ); }

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::ConformWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::ConformWith(DistMatrix[MC,MR])");
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
Elemental::DistMatrix<T,MC,Star>::ConformWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::ConformWith(DistMatrix[MC,* ])");
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
Elemental::DistMatrix<T,MC,Star>::ConformWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::ConformWith(DistMatrix[MR,MC])");
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
Elemental::DistMatrix<T,MC,Star>::ConformWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::ConformWith(DistMatrix[* ,MC])");
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
Elemental::DistMatrix<T,MC,Star>::FreeConstraints()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::FreeConstraints");
    CHECK_IF_NONEXISTANT_CONSTRAINTS;
#endif
    _constrainedColDist = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::View
( DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::View(A)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _colAlignment = A.ColAlignment();
    _colShift     = A.ColShift();
    _localMatrix.View( A.LocalMatrix() );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::LockedView
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::LockedView(A)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
#endif
    _height = A.Height();
    _width  = A.Width();
    _colAlignment = A.ColAlignment();
    _colShift     = A.ColShift();
    _localMatrix.LockedView( A.LockedLocalMatrix() );
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::View
( DistMatrix<T,MC,Star>& A,
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::View(A,i,j,height,width)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int r   = _grid->Height();
        const int row = _grid->MCRank();

        _colAlignment = (A.ColAlignment()+i) % r;
        _colShift = Shift( row, _colAlignment, r );

        const int localHeightBefore = LocalLength( i, A.ColShift(), r );
        const int localHeight = LocalLength( height, _colShift, r );

        _localMatrix.View
        ( A.LocalMatrix(), localHeightBefore, j, localHeight, width );
    }
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::LockedView
( const DistMatrix<T,MC,Star>& A,
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::LockedView(A,i,j,height,width)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int r   = _grid->Height();
        const int row = _grid->MCRank();

        _colAlignment = (A.ColAlignment()+i) % r;
        _colShift = Shift( row, _colAlignment, r );

        const int localHeightBefore = LocalLength( i, A.ColShift(), r );
        const int localHeight = LocalLength( height, _colShift, r);

        _localMatrix.LockedView
        ( A.LockedLocalMatrix(), localHeightBefore, j, localHeight, width );
    }
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::View1x2
( DistMatrix<T,MC,Star>& AL,
  DistMatrix<T,MC,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::View1x2");    
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
    if( AL.ColAlignment() != AR.ColAlignment() )
    {
        if( _grid->VCRank() == 0 )
            cerr << "1x2 is misaligned, cannot horizontally combine." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _colAlignment = AL.ColAlignment();
    _colShift     = AL.ColShift();
    _localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::LockedView1x2
( const DistMatrix<T,MC,Star>& AL,
  const DistMatrix<T,MC,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::LockedView1x2");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AL );
    CHECK_IF_VIEWING_DIFF_GRID( AR );
    CHECK_IF_CONFORMING_1x2( AL, AR );
    if( AL.ColAlignment() != AR.ColAlignment() )
    {
        if( _grid->VCRank() == 0 )
            cerr << "1x2 is misaligned, cannot horizontally combine." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _colAlignment = AL.ColAlignment();
    _colShift     = AL.ColShift();
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
Elemental::DistMatrix<T,MC,Star>::View2x1
( DistMatrix<T,MC,Star>& AT,
  DistMatrix<T,MC,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::View2x1");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _colAlignment = AT.ColAlignment();
    _colShift     = AT.ColShift();
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
Elemental::DistMatrix<T,MC,Star>::LockedView2x1
( const DistMatrix<T,MC,Star>& AT,
  const DistMatrix<T,MC,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::LockedView2x1");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( AT );
    CHECK_IF_VIEWING_DIFF_GRID( AB );
    CHECK_IF_CONFORMING_2x1( AT, AB );
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _colAlignment = AT.ColAlignment();
    _colShift     = AT.ColShift();
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
Elemental::DistMatrix<T,MC,Star>::View2x2
( DistMatrix<T,MC,Star>& ATL,
  DistMatrix<T,MC,Star>& ATR,
  DistMatrix<T,MC,Star>& ABL,
  DistMatrix<T,MC,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::View2x2");
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
        if( _grid->VCRank() == 0 )
            cerr << "2x2 must align to combine." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _colAlignment = ATL.ColAlignment();
    _colShift     = ATL.ColShift();
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
Elemental::DistMatrix<T,MC,Star>::LockedView2x2
( const DistMatrix<T,MC,Star>& ATL,
  const DistMatrix<T,MC,Star>& ATR,
  const DistMatrix<T,MC,Star>& ABL,
  const DistMatrix<T,MC,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::LockedView2x2");
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
        if( _grid->VCRank() == 0 )
            cerr << "2x2 must align to combine." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _colAlignment = ATL.ColAlignment();
    _colShift     = ATL.ColShift();
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
Elemental::DistMatrix<T,MC,Star>::ResizeTo
( const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::ResizeTo");
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
    _localMatrix.ResizeTo( LocalLength(height,_colShift,_grid->Height()), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
Elemental::DistMatrix<T,MC,Star>::Get
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix<MC,* >::Get");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        if( _grid->VCRank() == 0 )
        {
            cerr << "Entry (" << i << "," << j << ") is out of bounds of "
                 << Height() << " x " << Width() << " matrix." << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that 
    // row within each process column
    const int ownerRow = (i + ColAlignment()) % _grid->Height();

    T u;
    if( _grid->MCRank() == ownerRow )
    {
        const int iLoc = (i-ColShift()) / _grid->Height();
        u = _localMatrix(iLoc,j);
    }
    Broadcast( &u, 1, ownerRow, _grid->MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::Set
( const int i, const int j, const T u )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::Set");
    if( i < 0 || i >= Height() || j < 0 || j >=Width() )
    {
        if( _grid->VCRank() == 0 )
        {
            cerr << "Entry (" << i << "," << j << ") is out of bounds of "
                 << Height() << " x " << Width() << " matrix." << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    const int ownerRow = (i + ColAlignment()) % _grid->Height();

    if( _grid->MCRank() == ownerRow )
    {
        const int iLoc = (i-ColShift()) / _grid->Height();
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
Elemental::DistMatrix<T,MC,Star>::MakeTrapezoidal
( const Side side, const Shape shape, const int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix<MC,* >::MakeTrapezoidal");
    CHECK_IF_LOCKED_VIEW;
#endif

    const int height = Height();
    const int width = Width();
    const int localHeight = LocalHeight();
    const int r = _grid->Height();
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
                const int numZeros = LocalLength( boundary, colShift, r );
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
            const int nonzeroLength = LocalLength(firstZero_i,colShift,r);
            for( int iLoc=nonzeroLength; iLoc<localHeight; ++iLoc )
                _localMatrix(iLoc,j) = (T)0;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("DistMatrix<MC,* >::SetToIdentity");
    CHECK_IF_LOCKED_VIEW;
#endif

    const int width = Width();
    const int localHeight = LocalHeight();
    const int r = _grid->Height();
    const int colShift = ColShift();

    _localMatrix.SetToZero();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;
        if( i < width )
            _localMatrix(iLoc,i) = (T)1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("DistMatrix<MC,* >::SetToRandom");
    CHECK_IF_LOCKED_VIEW;
#endif

    const int width       = Width();
    const int localHeight = LocalHeight();
    const int bufSize     = localHeight*width;

    _auxMemory.Require( bufSize );

    // Create random matrix on process column 0, then broadcast
    T* buffer = _auxMemory.Buffer();
    if( _grid->MRRank() == 0 )
    {
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                buffer[i+j*localHeight] = Random<T>();
    }
    Broadcast( buffer, bufSize, 0, _grid->MRComm() );

    // Unpack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = buffer[i+j*localHeight];

    _auxMemory.Release();

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::SetToRandomDiagDominant()
{
#ifndef RELEASE
    PushCallStack("DistMatrix<MC,* >::SetToRandomDiagDominant");
    CHECK_IF_LOCKED_VIEW;
#endif

    const int width = Width();
    const int localHeight = LocalHeight();
    const int r = _grid->Height();
    const int colShift = ColShift();

    SetToRandom();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;
        if( i < width )
            _localMatrix(iLoc,i) += (T)max(Height(),Width());
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::DistMatrix<T,MC,Star>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("DistMatrix<MC,* >::SetToZero");
#endif
    _localMatrix.SetToZero();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void 
Elemental::DistMatrix<T,MC,Star>::AllReduce()
{
#ifndef RELEASE
    PushCallStack("DistMatrix<MC,* >::AllReduce");
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
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix<MC,* > = DistMatrix<MC,MR>");
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
        if( A.Width() == 1 )
        {
            if( _grid->MRRank() == A.RowAlignment() )
                _localMatrix = A.LockedLocalMatrix();

            // Communicate
            Broadcast
            ( _localMatrix.Buffer(), LocalHeight(),
              A.RowAlignment(), _grid->MRComm() );
        }
        else
        {
            const int c = _grid->Width();

            const int width = Width();
            const int localHeight = LocalHeight();
            const int localWidthOfA = A.LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,c);

            const int portionSize = 
                max(localHeight*maxLocalWidth,MinCollectContrib);

            _auxMemory.Require( (c+1)*portionSize );

            T* buffer = _auxMemory.Buffer();
            T* originalData = &buffer[0];
            T* gatheredData = &buffer[portionSize];

            // Pack
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<localHeight; ++i )
                    originalData[i+j*localHeight] = A.LocalEntry(i,j);

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
                    for( int i=0; i<localHeight; ++i )
                        _localMatrix(i,rowShift+j*c) = data[i+j*localHeight];
            }

            _auxMemory.Release();
        }
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MC,* ] <- [MC,MR]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int row = _grid->MCRank();

        const int colAlignment = ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
        const int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

        if( A.Width()==1 )
        {
            const int localHeight = LocalHeight();

            if( _grid->MRRank() == A.RowAlignment() )
            {
                const int localHeightOfA = A.LocalHeight();

                _auxMemory.Require( localHeightOfA );
                T* buffer = _auxMemory.Buffer();

                // Pack
                for( int i=0; i<localHeightOfA; ++i )
                    buffer[i] = A.LocalEntry(i,0);

                // Communicate
                SendRecv
                ( buffer, localHeightOfA, sendRow, 0,
                  _localMatrix.Buffer(), localHeight, recvRow, MPI_ANY_TAG, 
                  _grid->MCComm()                                     );

                _auxMemory.Release();
            }

            // Communicate
            Broadcast
            ( _localMatrix.Buffer(), localHeight, 
              A.RowAlignment(), _grid->MRComm() );
        }
        else
        {
            const int height = Height();
            const int width = Width();
            const int localHeight = LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int localWidthOfA  = A.LocalWidth();
            const int maxLocalHeight = MaxLocalLength(height,r);
            const int maxLocalWidth  = MaxLocalLength(width,c);

            const int portionSize = 
                max(maxLocalHeight*maxLocalWidth,MinCollectContrib);

            _auxMemory.Require( (c+1)*portionSize );

            T* buffer = _auxMemory.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[portionSize];

            // Pack the currently owned local data of A into the second buffer
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    secondBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);

            // Perform the SendRecv: puts the new data into the first buffer
            SendRecv
            ( secondBuffer, portionSize, sendRow, 0,
              firstBuffer,  portionSize, recvRow, MPI_ANY_TAG, 
              _grid->MCComm()                            );

            // Use the output of the SendRecv as the input to the AllGather
            AllGather
            ( firstBuffer,  portionSize, 
              secondBuffer, portionSize, _grid->MRComm() );

            // Unpack the contents of each member of the process row
            const int rowAlignmentOfA = A.RowAlignment();
            for( int k=0; k<c; ++k )
            {
                const T* data = &secondBuffer[k*portionSize];

                const int rowShift = Shift( k, rowAlignmentOfA, c ); 
                const int localWidth = LocalLength( width, rowShift, c );
                for( int j=0; j<localWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        _localMatrix(i,rowShift+j*c) = data[i+j*localHeight];
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
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[MC,* ]");
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
            cout << "Unaligned [MC,* ] <- [MC,* ]." << endl;
#endif
        const int rank = _grid->MCRank();
        const int r = _grid->Height();

        const int colAlignment = ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (rank+r+colAlignment-colAlignmentOfA) % r;
        const int recvRank = (rank+r+colAlignmentOfA-colAlignment) % r;

        const int width = Width();
        const int localHeight = LocalHeight();
        const int localHeightOfA = A.LocalHeight();

        const int sendSize = localHeightOfA * width;
        const int recvSize = localHeight * width;

        _auxMemory.Require( sendSize + recvSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, _grid->MCComm() );

        // Unpack
        for( int j=0; j<width; ++j )
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
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[* ,MR]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,MC,MR> 
        A_MC_MR(true,ColAlignment(),false,0,*_grid);

    A_MC_MR = A;
    *this   = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[MD,* ]");
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
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[* ,MD]");
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
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[MR,MC]");
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
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[MR,* ]");
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
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[* ,MC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
        ( new DistMatrix<T,MR,MC>(*_grid) );
    *A_MR_MC = A;

    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
        ( new DistMatrix<T,VR,Star>(*_grid) );
    *A_VR_Star = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

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
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[VC,* ]");
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

        const int height = Height();
        const int width = Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(height,p);

        const int portionSize = max(maxLocalHeightOfA*width,MinCollectContrib);

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
        const int colShift = ColShift();
        const int colAlignmentOfA = A.ColAlignment();
        for( int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];    

            const int colShiftOfA = Shift( row+r*k, colAlignmentOfA, p );
            const int colOffset = (colShiftOfA-colShift) / r;
            const int localHeight = LocalLength( height, colShiftOfA, p );

            for( int j=0; j<width; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(colOffset+i*c,j) = data[i+j*localHeight];
        }

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [MC,* ] <- [VC,* ]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = r * c;
        const int row = _grid->MCRank();
        const int rank = _grid->VCRank();

        // Perform the SendRecv to make A have the same colAlignment
        const int colAlignment = ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int colShift = ColShift();

        const int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

        const int height = Height();
        const int width = Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(height,p);

        const int portionSize = max(maxLocalHeightOfA*width,MinCollectContrib);

        _auxMemory.Require( (c+1)*portionSize );

        T* buffer = _auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                secondBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);

        // Perform the SendRecv: puts the new data into the first buffer
        SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, MPI_ANY_TAG, 
          _grid->VCComm() );

        // Use the SendRecv as input to the AllGather
        AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, _grid->MRComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int colShiftOfA = Shift( row+r*k, colAlignment, p );
            const int colOffset = (colShiftOfA-colShift) / r;
            const int localHeight = LocalLength( height, colShiftOfA, p );

            for( int j=0; j<width; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(colOffset+i*c,j) = data[i+j*localHeight];
        }

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[* ,VC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
        ( new DistMatrix<T,Star,VR>(*_grid) );
    *A_Star_VR = A;

    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
        ( new DistMatrix<T,MC,MR>(true,ColAlignment(),false,0,*_grid) );
    *A_MC_MR = *A_Star_VR;
    delete A_Star_VR.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[VR,* ]");
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
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[* ,VR]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    DistMatrix<T,MC,MR> A_MC_MR(true,ColAlignment(),false,0,*_grid);

    A_MC_MR = A;
    *this   = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,MC,Star>&
Elemental::DistMatrix<T,MC,Star>::operator=
( const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ] = DistMatrix[* ,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int r = _grid->Height(); 
    const int colShift = ColShift();

    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = A.LocalEntry(colShift+i*r,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class Elemental::DistMatrix<int,     MC,Star>;
template class Elemental::DistMatrix<float,   MC,Star>;
template class Elemental::DistMatrix<double,  MC,Star>;
#ifndef WITHOUT_COMPLEX
template class Elemental::DistMatrix<scomplex,MC,Star>;
template class Elemental::DistMatrix<dcomplex,MC,Star>;
#endif

