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
#include "elemental/dist_matrix.hpp"
#include "./dist_matrix_macros.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::utilities;
using namespace elemental::wrappers::mpi;

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::Print( const string& msg ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::Print");
#endif
    if( _grid->VCRank() == 0 && msg != "" )
        cout << msg << endl;

    const int height      = Height();
    const int width       = Width();
    const int localHeight = LocalHeight();
    const int p           = _grid->Size();
    const int colShift    = ColShift();

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
    for( int i=0; i<localHeight; ++i )
        for( int j=0; j<width; ++j )
            sendBuf[colShift+i*p+j*height] = _localMatrix(i,j);

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
elemental::DistMatrix<T,VC,Star>::AlignWith
( const DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::AlignWith(DistMatrix[MC,MR])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment();
    _colShift = Shift( _grid->VCRank(), _colAlignment, _grid->Size() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::AlignWith
( const DistMatrix<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::AlignWith(DistMatrix[MR,MC])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment();
    _colShift = Shift( _grid->VCRank(), _colAlignment, _grid->Size() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::AlignWith
( const DistMatrix<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::AlignWith(DistMatrix[MC,* ])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.ColAlignment();
    _colShift = Shift( _grid->VCRank(), _colAlignment, _grid->Size() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::AlignWith
( const DistMatrix<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::AlignWith(DistMatrix[* ,MC])");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_ALIGNING_DIFF_GRID( A );
#endif
    _colAlignment = A.RowAlignment();
    _colShift = Shift( _grid->VCRank(), _colAlignment, _grid->Size() );
    _constrainedColDist = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::AlignWith
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::AlignWith(DistMatrix[VC,* ])");
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
elemental::DistMatrix<T,VC,Star>::AlignWith
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::AlignWith(DistMatrix[* ,VC])");
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
elemental::DistMatrix<T,VC,Star>::AlignColsWith
( const DistMatrix<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::AlignColsWith
( const DistMatrix<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::AlignColsWith
( const DistMatrix<T,MC,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::AlignColsWith
( const DistMatrix<T,Star,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::AlignColsWith
( const DistMatrix<T,VC,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::AlignColsWith
( const DistMatrix<T,Star,VC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::ConformWith
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::ConformWith(DistMatrix[VC,* ])");
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
elemental::DistMatrix<T,VC,Star>::ConformWith
( const DistMatrix<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::ConformWith(DistMatrix[* ,VC])");
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
elemental::DistMatrix<T,VC,Star>::FreeConstraints()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::FreeConstraints");
#endif
    _constrainedColDist = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::View
( DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::View(A)");
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
elemental::DistMatrix<T,VC,Star>::LockedView
( const DistMatrix<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::LockedView(A)");
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
elemental::DistMatrix<T,VC,Star>::View
( DistMatrix<T,VC,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::View(A,i,j,height,width)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int colMajorRank = _grid->VCRank();
        const int size = _grid->Size();

        _colAlignment = (A.ColAlignment()+i) % size;
        _colShift = Shift( colMajorRank, _colAlignment, size );

        const int localHeightBefore = LocalLength(i,A.ColShift(),size);
        const int localHeight = LocalLength(height,_colShift,size);

        _localMatrix.View( A.LocalMatrix(),
                           localHeightBefore, j, localHeight, width );
    }
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::LockedView
( const DistMatrix<T,VC,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::LockedView(A,i,j,height,width)");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( A );
    CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width );
#endif
    _height = height;
    _width  = width;
    {
        const int colMajorRank = _grid->VCRank();
        const int size = _grid->Size();

        _colAlignment = (A.ColAlignment()+i) % size;
        _colShift = Shift( colMajorRank, _colAlignment, size );

        const int localHeightBefore = LocalLength(i,A.ColShift(),size);
        const int localHeight = LocalLength(height,_colShift,size);

        _localMatrix.LockedView( A.LockedLocalMatrix(),
                                 localHeightBefore, j, localHeight, width );
    }
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::View1x2
( DistMatrix<T,VC,Star>& AL,
  DistMatrix<T,VC,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::View1x2");
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
elemental::DistMatrix<T,VC,Star>::LockedView1x2
( const DistMatrix<T,VC,Star>& AL,
  const DistMatrix<T,VC,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::LockedView1x2");
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
elemental::DistMatrix<T,VC,Star>::View2x1
( DistMatrix<T,VC,Star>& AT,
  DistMatrix<T,VC,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::View2x1");
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
    _localMatrix.View2x1( AT.LocalMatrix(), AB.LocalMatrix() );
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::LockedView2x1
( const DistMatrix<T,VC,Star>& AT,
  const DistMatrix<T,VC,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::LockedView2x1");
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
elemental::DistMatrix<T,VC,Star>::View2x2
( DistMatrix<T,VC,Star>& ATL,
  DistMatrix<T,VC,Star>& ATR,
  DistMatrix<T,VC,Star>& ABL,
  DistMatrix<T,VC,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::View2x2");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
    if( ATL.ColAlignment() != ATR.ColAlignment() ||
        ABL.ColAlignment() != ABR.ColAlignment()    )
    {
        throw "Cannot combine misaligned 2x2 grid of matrices.";
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
elemental::DistMatrix<T,VC,Star>::LockedView2x2
( const DistMatrix<T,VC,Star>& ATL,
  const DistMatrix<T,VC,Star>& ATR,
  const DistMatrix<T,VC,Star>& ABL,
  const DistMatrix<T,VC,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::LockedView2x2");
    CHECK_IF_UNFREED_COL_CONSTRAINT;
    CHECK_IF_VIEWING_AND_STORING;
    CHECK_IF_VIEWING_DIFF_GRID( ATL );
    CHECK_IF_VIEWING_DIFF_GRID( ATR );
    CHECK_IF_VIEWING_DIFF_GRID( ABL );
    CHECK_IF_VIEWING_DIFF_GRID( ABR );
    CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR );
    if( ATL.ColAlignment() != ATR.ColAlignment() ||
        ABL.ColAlignment() != ABR.ColAlignment()    )
    {
        throw "Cannot combine misaligned 2x2 grid of matrices.";
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
elemental::DistMatrix<T,VC,Star>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::ResizeTo");
    CHECK_IF_LOCKED_VIEW;
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _height = height;
    _width  = width;
    _localMatrix.ResizeTo(LocalLength(height,_colShift,_grid->Size()),width);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrix<T,VC,Star>::Get
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::Get");
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
    const int ownerRank = (i + ColAlignment()) % _grid->Size();

    T u;
    if( _grid->VCRank() == ownerRank )
    {
        const int iLoc = (i-ColShift()) / _grid->Size();
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
elemental::DistMatrix<T,VC,Star>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::Set");
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix." << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const int ownerRank = (i + ColAlignment()) % _grid->Size();

    if( _grid->VCRank() == ownerRank )
    {
        const int iLoc = (i-ColShift()) / _grid->Size();
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
elemental::DistMatrix<T,VC,Star>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::MakeTrapezoidal");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int height = Height();
    const int width = Width();
    const int localHeight = LocalHeight();
    const int p = _grid->Size();
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
                const int numZeros = LocalLength( boundary, colShift, p );
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
            const int nonzeroLength = LocalLength(firstZero_i,colShift,p);
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
elemental::DistMatrix<T,VC,Star>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::SetToIdentity");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int width       = Width();
    const int localHeight = LocalHeight();
    const int p           = _grid->Size();
    const int colShift    = ColShift();

    _localMatrix.SetToZero();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*p;
        if( i < width )
            _localMatrix(iLoc,i) = (T)1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::SetToRandom");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int width       = Width();
    const int localHeight = LocalHeight();
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = Random<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::SetToRandomDiagDominant()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::SetToRandomDiagDominant");
    CHECK_IF_LOCKED_VIEW;
#endif
    const int width       = Width();
    const int localHeight = LocalHeight();
    const int p           = _grid->Size();
    const int colShift    = ColShift();

    SetToRandom();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*p;
        if( i < width )
            _localMatrix(iLoc,i) += (T)max(Height(),Width());
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrix<T,VC,Star>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::SetToZero");
    CHECK_IF_LOCKED_VIEW;
#endif
    _localMatrix.SetToZero();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[MC,MR]");
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
                        ( _grid->VCRank(), _colAlignment, _grid->Size() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( ColAlignment() % _grid->Height() == A.ColAlignment() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int row = _grid->MCRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int height = Height();
        const int width = Width();
        const int localHeight = LocalHeight();
        const int localWidthOfA = A.LocalWidth();

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

            const int thisRank = row+k*r;
            const int thisColShift = Shift(thisRank,colAlignment,p); 
            const int thisColOffset = (thisColShift-colShiftOfA) / r;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColOffset+i*c,j);
        }

        // Communicate
        AllToAll( sendBuffer, portionSize,
                  recvBuffer, portionSize, _grid->MRComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignmentOfA,c);
            const int thisLocalWidth = LocalLength(width,thisRowShift,c);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,thisRowShift+j*c) = data[i+j*localHeight];
        }

        _auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [VC,* ] <- [MC,MR]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int row = _grid->MCRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        
        const int sendRow = (row+r+(colAlignment%r)-colAlignmentOfA) % r;
        const int recvRow = (row+r+colAlignmentOfA-(colAlignment%r)) % r;

        const int height = Height();
        const int width = Width();
        const int localHeight = LocalHeight();
        const int localWidthOfA = A.LocalWidth();

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

            const int thisRank = sendRow+k*r;
            const int thisColShift = Shift(thisRank,colAlignment,p);
            const int thisColOffset = (thisColShift-colShiftOfA) / r;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.LocalEntry(thisColOffset+i*c,j);
        }

        // AllToAll to gather all of the unaligned [VC,*] data into firstBuffer
        AllToAll( secondBuffer, portionSize, 
                  firstBuffer,  portionSize, _grid->MRComm() );

        // SendRecv: properly align the [VC,*] via a trade in the column
        SendRecv
        ( firstBuffer,  portionSize, sendRow, 0,
          secondBuffer, portionSize, recvRow, MPI_ANY_TAG, 
          _grid->MCComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignmentOfA,c);
            const int thisLocalWidth = LocalLength(width,thisRowShift,c);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    _localMatrix(i,thisRowShift+j*c) = data[i+j*localHeight];
        }

        _auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[MC,* ]");
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
                        ( _grid->VCRank(), _colAlignment, _grid->Size() );
        }
        ResizeTo( A.Height(), A.Width() );
    }

    if( ColAlignment() % _grid->Height() == A.ColAlignment() )
    {
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int colShift = ColShift();
        const int colShiftOfA = A.ColShift();
        const int colOffset = (colShift-colShiftOfA) / r;
        
        const int width = Width();
        const int localHeight = LocalHeight();

        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                _localMatrix(i,j) = A.LocalEntry(colOffset+i*c,j);
    }
    else
    {
#ifndef RELEASE
        if( _grid->VCRank() == 0 )
            cout << "Unaligned [VC,* ] <- [MC,* ]." << endl;
#endif
        const int r = _grid->Height();
        const int c = _grid->Width();
        const int p = _grid->Size();
        const int row = _grid->MCRank();
        const int col = _grid->MRRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[VC,*] within our process column to fix alignments.
        const int sendRow = (row+r+(colAlignment%r)-colAlignmentOfA) % r;
        const int recvRow = (row+r+colAlignmentOfA-(colAlignment%r)) % r;
        const int sendRank = sendRow + r*col;

        const int sendColShift = Shift( sendRank, colAlignment, p );
        const int sendColOffset = (sendColShift-colShiftOfA) / r;

        const int height = Height();
        const int width = Width();
        const int localHeight = LocalHeight();
        const int localHeightOfSend = LocalLength(height,sendColShift,p);

        const int sendSize = localHeightOfSend * width;
        const int recvSize = localHeight * width;

        _auxMemory.Require( sendSize + recvSize );

        T* buffer = _auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfSend; ++i )
                sendBuffer[i+j*localHeightOfSend] = 
                      A.LocalEntry(sendColOffset+i*c,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRow, 0,
          recvBuffer, recvSize, recvRow, MPI_ANY_TAG, _grid->MCComm() );

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
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[* ,MR]");
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
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[MD,* ]");
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
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[* ,MD]");
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
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[MR,MC]");
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
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[MR,* ]");
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
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[* ,MC]");
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

    *this = *A_VR_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[VC,* ]");
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
            cout << "Unaligned [VC,* ] <- [VC,* ]." << endl;
#endif
        const int rank = _grid->VCRank();
        const int p = _grid->Size();

        const int colAlignment = ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

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
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, _grid->VCComm() );

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
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[* ,VC]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC
        ( new DistMatrix<T,MR,MC>(*_grid) );
    *A_MR_MC = A;

    auto_ptr< DistMatrix<T,VC,Star> > A_VR_Star
        ( new DistMatrix<T,VC,Star>(*_grid) );
    *A_VR_Star = *A_MR_MC; 
    delete A_MR_MC.release(); // lowers memory highwater

    *this = *A_VR_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[VR,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );
    
    const int width = Width();
    const int localHeight = LocalHeight();
    const int localHeightOfA = A.LocalHeight();

    const int sendSize = localHeightOfA * width;
    const int recvSize = localHeight * width;

    const int r = _grid->Height();
    const int c = _grid->Width();
    const int p = _grid->Size();
    const int rankCM = _grid->VCRank();
    const int rankRM = _grid->VRRank();

    const int colShift = ColShift();
    const int colShiftOfA = A.ColShift();

    // Compute which colmajor rank has the colShift equal to our colShiftOfA
    const int sendRankCM = (rankCM+(p+colShiftOfA-colShift)) % p;

    // Compute which colmajor rank has the A colShift that we need
    const int recvRankRM = (rankRM+(p+colShift-colShiftOfA)) % p;
    const int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

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
    ( sendBuffer, sendSize, sendRankCM, 0,
      recvBuffer, recvSize, recvRankCM, MPI_ANY_TAG, _grid->VCComm() );

    // Unpack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = recvBuffer[i+j*localHeight];

    _auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[* ,VR]");
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
const DistMatrix<T,VC,Star>&
elemental::DistMatrix<T,VC,Star>::operator=
( const DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ] = DistMatrix[* ,* ]");
    CHECK_IF_LOCKED_VIEW;
    CHECK_IF_REDIST_DIFF_GRID( A );
    CHECK_IF_VIEWING_DIFF_SIZE( A );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int p = _grid->Size();
    const int colShift = ColShift();

    const int localHeight = LocalHeight();
    const int localWidth = LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            _localMatrix(i,j) = A.LocalEntry(colShift+i*p,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class elemental::DistMatrix<int,     VC,Star>;
template class elemental::DistMatrix<float,   VC,Star>;
template class elemental::DistMatrix<double,  VC,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,VC,Star>;
template class elemental::DistMatrix<dcomplex,VC,Star>;
#endif

