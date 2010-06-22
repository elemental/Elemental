/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/dist_matrix.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::utilities;
using namespace elemental::wrappers::mpi;

//----------------------------------------------------------------------------//
// DistMatrixBase                                                             //
//----------------------------------------------------------------------------//

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Print");
#endif
    const Grid& g = this->GetGrid();
    if( g.VCRank() == 0 && s != "" )
        cout << s << endl;

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

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
            sendBuf[colShift+i*c + (rowShift+j*r)*height] = 
                this->LocalEntry(i,j);

    // If we are the root, fill the receive buffer
    T* recvBuf = 0;
    if( g.VCRank() == 0 )
    {
        recvBuf = new T[height*width];
        for( int i=0; i<height*width; ++i )
            recvBuf[i] = (T)0;
    }

    // Sum the contributions and send to the root
    Reduce( sendBuf, recvBuf, height*width, MPI_SUM, 0, g.VCComm() );
    delete[] sendBuf;

    if( g.VCRank() == 0 )
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
    Barrier( g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_rowAlignment = A.RowAlignment();
    this->_colShift     = A.ColShift();
    this->_rowShift     = A.RowShift();
    this->_constrainedColAlignment = true; 
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_colShift = A.ColShift();
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([* ,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = A.RowShift();
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_rowAlignment = A.ColAlignment();
    this->_colShift     = A.RowShift();
    this->_rowShift     = A.ColShift();
    this->_constrainedColAlignment = true;
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MC,*])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = A.ColShift();
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([* ,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_colShift = A.RowShift();
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([VC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = 
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignWith([* ,VC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = 
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([VR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignWith([* ,VR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.RowAlignment();
    this->_colShift = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignColsWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_colShift = A.ColShift();
    this->_constrainedColAlignment = true; 
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignColsWith
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([MR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_colShift = A.ColShift();
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignColsWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_colShift = A.RowShift();
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignColsWith
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([* ,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_colShift = A.RowShift();
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignColsWith
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([VR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignColsWith
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignColsWith([* ,VR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.RowAlignment();
    this->_colShift = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignRowsWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = A.RowShift();
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignRowsWith
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([* ,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = A.RowShift();
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignRowsWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = A.ColShift();
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignRowsWith
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([MC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = A.ColShift();
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignRowsWith
( const DistMatrixBase<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([VC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = 
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::AlignRowsWith
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignRowsWith([* ,VC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = 
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::View
( DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::View");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width  = A.Width();
    this->_colAlignment = A.ColAlignment();
    this->_rowAlignment = A.RowAlignment();
    this->_colShift     = A.ColShift();
    this->_rowShift     = A.RowShift();
    this->_localMatrix.View( A.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::LockedView
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width  = A.Width();
    this->_colAlignment = A.ColAlignment();
    this->_rowAlignment = A.RowAlignment();
    this->_colShift     = A.ColShift();
    this->_rowShift     = A.RowShift();
    this->_localMatrix.LockedView( A.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::View
( DistMatrixBase<T,MR,MC>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::View");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const Grid& g = this->GetGrid();
        const int r   = g.Height();
        const int c   = g.Width();
        const int row = g.MCRank();
        const int col = g.MRRank();

        this->_colAlignment = (A.ColAlignment()+i) % c;
        this->_rowAlignment = (A.RowAlignment()+j) % r;
        
        this->_colShift = Shift( col, this->ColAlignment(), c );
        this->_rowShift = Shift( row, this->RowAlignment(), r );

        const int localHeightBefore = LocalLength( i, A.ColShift(), c );
        const int localWidthBefore  = LocalLength( j, A.RowShift(), r );

        const int localHeight = LocalLength( height, this->ColShift(), c );
        const int localWidth  = LocalLength( width,  this->RowShift(), r );

        this->_localMatrix.View
        ( A.LocalMatrix(),
          localHeightBefore, localWidthBefore, localHeight, localWidth );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::LockedView
( const DistMatrixBase<T,MR,MC>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const Grid& g = this->GetGrid();
        const int r   = g.Height();
        const int c   = g.Width();
        const int row = g.MCRank();
        const int col = g.MRRank();

        this->_colAlignment = (A.ColAlignment()+i) % c;
        this->_rowAlignment = (A.RowAlignment()+j) % r;
        
        this->_colShift = Shift( col, this->ColAlignment(), c );
        this->_rowShift = Shift( row, this->RowAlignment(), r );

        const int localHeightBefore = LocalLength( i, A.ColShift(), c );
        const int localWidthBefore  = LocalLength( j, A.RowShift(), r );

        const int localHeight = LocalLength( height, this->ColShift(), c );
        const int localWidth  = LocalLength( width,  this->RowShift(), r );

        this->_localMatrix.LockedView
        ( A.LockedLocalMatrix(),
          localHeightBefore, localWidthBefore, localHeight, localWidth );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::View1x2
( DistMatrixBase<T,MR,MC>& AL, DistMatrixBase<T,MR,MC>& AR )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::View1x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width  = AL.Width() + AR.Width();
    this->_colAlignment = AL.ColAlignment();
    this->_rowAlignment = AL.RowAlignment();
    this->_colShift     = AL.ColShift();
    this->_rowShift     = AL.RowShift();
    this->_localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::LockedView1x2
( const DistMatrixBase<T,MR,MC>& AL,
  const DistMatrixBase<T,MR,MC>& AR )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::LockedView1x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width  = AL.Width() + AR.Width();
    this->_colAlignment = AL.ColAlignment();
    this->_rowAlignment = AL.RowAlignment();
    this->_colShift     = AL.ColShift();
    this->_rowShift     = AL.RowShift();
    this->_localMatrix.LockedView1x2
    ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::View2x1
( DistMatrixBase<T,MR,MC>& AT,
  DistMatrixBase<T,MR,MC>& AB )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::View2x1");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width  = AT.Width();
    this->_colAlignment = AT.ColAlignment();
    this->_rowAlignment = AT.RowAlignment();
    this->_colShift     = AT.ColShift();
    this->_rowShift     = AT.RowShift();
    this->_localMatrix.View2x1
    ( AT.LocalMatrix(),
      AB.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::LockedView2x1
( const DistMatrixBase<T,MR,MC>& AT,
  const DistMatrixBase<T,MR,MC>& AB )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::LockedView2x1");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width  = AT.Width();
    this->_colAlignment = AT.ColAlignment();
    this->_rowAlignment = AT.RowAlignment();
    this->_colShift     = AT.ColShift();
    this->_rowShift     = AT.RowShift();
    this->_localMatrix.LockedView2x1
    ( AT.LockedLocalMatrix(),
      AB.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::View2x2
( DistMatrixBase<T,MR,MC>& ATL,
  DistMatrixBase<T,MR,MC>& ATR,
  DistMatrixBase<T,MR,MC>& ABL,
  DistMatrixBase<T,MR,MC>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::View2x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width  = ATL.Width() + ATR.Width();
    this->_colAlignment = ATL.ColAlignment();
    this->_rowAlignment = ATL.RowAlignment();
    this->_colShift     = ATL.ColShift();
    this->_rowShift     = ATL.RowShift();
    this->_localMatrix.View2x2
    ( ATL.LocalMatrix(), ATR.LocalMatrix(),
      ABL.LocalMatrix(), ABR.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::LockedView2x2
( const DistMatrixBase<T,MR,MC>& ATL,
  const DistMatrixBase<T,MR,MC>& ATR,
  const DistMatrixBase<T,MR,MC>& ABL,
  const DistMatrixBase<T,MR,MC>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::LockedView2x2");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width  = ATL.Width() + ATR.Width();
    this->_colAlignment = ATL.ColAlignment();
    this->_rowAlignment = ATL.RowAlignment();
    this->_colShift     = ATL.ColShift();
    this->_rowShift     = ATL.RowShift();
    this->_localMatrix.LockedView2x2
    ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
      ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    const Grid& g = this->GetGrid();
    this->_height = height;
    this->_width = width;
    this->_localMatrix.ResizeTo
    ( LocalLength( height, this->ColShift(), g.Width() ),
      LocalLength( width,  this->RowShift(), g.Height() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,MR,MC>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const Grid& g = this->GetGrid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        u = this->LocalEntry(iLoc,jLoc);
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        this->LocalEntry(iLoc,jLoc) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::GetDiagonal
( DistMatrixBase<T,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetDiagonal([MD,* ])");
    this->AssertNotLockedView();
#endif
    int height = this->Height();
    int width = this->Width();
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
        throw logic_error( "d is not of the correct length." );
#endif
    if( !d.Viewing() )
    {
        if( !d.ConstrainedColAlignment() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( length, 1 );
    }

    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            d.LocalEntry(k,0) = 
                this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::GetDiagonal
( DistMatrixBase<T,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetDiagonal([* ,MD])");
    this->AssertNotLockedView();
#endif
    int height = this->Height();
    int width = this->Width();
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
        throw logic_error( "d is not of the correct length." );
#endif

    if( !d.Viewing() )
    {
        if( !d.ConstrainedRowAlignment() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( 1, length );
    }

    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            d.LocalEntry(0,k) = 
                this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::SetDiagonal
( const DistMatrixBase<T,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetDiagonal([MD,* ])");
    if( d.Width() != 1 )
        throw logic_error( "d must be a column vector." );
    {
        int height = this->Height();
        int width = this->Width();
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
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r)) = 
                d.LocalEntry(k,0);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::SetDiagonal
( const DistMatrixBase<T,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetDiagonal([* ,MD])");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        int height = this->Height();
        int width = this->Width();
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
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r)) = 
                d.LocalEntry(0,k);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., SetToIdentity and MakeTrapezoidal
//

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::MakeTrapezoidal");
    this->AssertNotLockedView(); 
#endif
    const Grid& g = this->GetGrid();
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

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
                    this->LocalEntry(iLoc,jLoc) = (T)0;
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
                this->LocalEntry(iLoc,jLoc) = (T)0;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->GetGrid();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToZero();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*c;
        if( i % r == rowShift )
        {
            const int jLoc = (i-rowShift) / r;
            if( jLoc < localWidth )
                this->LocalEntry(iLoc,jLoc) = (T)1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = Random<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( A.Width() == 1 )
    {
        if( !this->Viewing() )
            ResizeTo( A.Height(), 1 );

        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int myRow = g.MCRank();
        const int myCol = g.MRRank();
        const int rankCM = g.VCRank();
        const int rankRM = g.VRRank();
        const int ownerRow = this->RowAlignment();
        const int ownerCol = A.RowAlignment();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = A.Height();
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = max(maxLocalHeight,MinCollectContrib);

        const int colShiftVR = Shift(rankRM,colAlignment,p);
        const int colShiftVCOfA = Shift(rankCM,colAlignmentOfA,p);
        const int sendRankRM = (rankRM+(p+colShiftVCOfA-colShiftVR)) % p;
        const int recvRankCM = (rankCM+(p+colShiftVR-colShiftVCOfA)) % p;
        const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

        this->_auxMemory.Require( (r+c)*portionSize );
        T* buffer = this->_auxMemory.Buffer();
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
          sendBuf, portionSize, ownerCol, g.MRComm() );

        // A[VR,* ] <- A[VC,* ]
        SendRecv
        ( sendBuf, portionSize, sendRankRM, 0,
          recvBuf, portionSize, recvRankRM, MPI_ANY_TAG, g.VRComm() );

        // A[MR,MC] <- A[VR,* ]
        Gather
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerRow, g.MCComm() );

        if( myRow == ownerRow )
        {
            // Unpack
            for( int k=0; k<r; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const int shift = Shift(myCol+c*k,colAlignment,p);
                const int offset = (shift-this->ColShift()) / c;
                const int thisLocalHeight = LocalLength(height,shift,p);

                for( int i=0; i<thisLocalHeight; ++i )
                    this->LocalEntry(offset+i*r,0) = data[i];
            }
        }

        this->_auxMemory.Release();
    }
    else if( A.Height() == 1 )
    {
        if( !this->Viewing() )
            ResizeTo( 1, A.Width() );

        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int myRow = g.MCRank();
        const int myCol = g.MRRank();
        const int rankCM = g.VCRank();
        const int rankRM = g.VRRank();
        const int ownerCol = this->ColAlignment();
        const int ownerRow = A.ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int width = A.Width();
        const int maxLocalWidth = MaxLocalLength(width,p);

        const int portionSize = max(maxLocalWidth,MinCollectContrib);

        const int rowShiftVC = Shift(rankCM,rowAlignment,p);
        const int rowShiftVROfA = Shift(rankRM,rowAlignmentOfA,p);
        const int sendRankCM = (rankCM+(p+rowShiftVROfA-rowShiftVC)) % p;
        const int recvRankRM = (rankRM+(p+rowShiftVC-rowShiftVROfA)) % p;
        const int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        this->_auxMemory.Require( (r+c)*portionSize );
        T* buffer = this->_auxMemory.Buffer();
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
          sendBuf, portionSize, ownerRow, g.MCComm() );

        // A[* ,VC] <- A[* ,VR]
        SendRecv
        ( sendBuf, portionSize, sendRankCM, 0,
          recvBuf, portionSize, recvRankCM, MPI_ANY_TAG, g.VCComm() );

        // A[MR,MC] <- A[* ,VC]
        Gather
        ( recvBuf, portionSize,
          sendBuf, portionSize, ownerCol, g.MRComm() );

        if( myCol == ownerCol )
        {
            // Unpack
            for( int k=0; k<c; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const int shift = Shift(myRow+r*k,rowAlignment,p);
                const int offset = (shift-this->RowShift()) / r;
                const int thisLocalWidth = LocalLength(width,shift,p);

                for( int j=0; j<thisLocalWidth; ++j )
                    this->LocalEntry(0,offset+j*c) = data[j];
            }
        }

        this->_auxMemory.Release();
    }
    else
    {
        if( A.Height() >= A.Width() )
        {
            auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
            ( new DistMatrix<T,VC,Star>(g) );
            *A_VC_Star = A;

            auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
            ( new DistMatrix<T,VR,Star>(true,this->ColAlignment(),g) );
            *A_VR_Star = *A_VC_Star;
            delete A_VC_Star.release(); // lowers memory highwater

            *this = *A_VR_Star;
        }
        else
        {
            auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
            ( new DistMatrix<T,Star,VR>(g) );
            *A_Star_VR = A;

            auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
            ( new DistMatrix<T,Star,VC>(true,this->RowAlignment(),g) );
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
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(g) );
    *A_VC_Star = A;

    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
    ( new DistMatrix<T,VR,Star>(true,this->ColAlignment(),g) );
    *A_VR_Star = *A_VC_Star;
    delete A_VC_Star.release(); // lowers memory highwater

    *this = *A_VR_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
    ( new DistMatrix<T,Star,VR>(g) );
    *A_Star_VR = A;

    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
    ( new DistMatrix<T,Star,VC>(true,this->RowAlignment(),g) );
    *A_Star_VC = *A_Star_VR;
    delete A_Star_VR.release(); // lowers memory highwater

    *this = *A_Star_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MR,MC] = [MD,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MR,MC] = [* ,MD] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.ColAlignment();
            this->_colShift = A.ColShift();
        }
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_rowShift = A.RowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() &&
        this->RowAlignment() == A.RowAlignment() )
    {
        this->_localMatrix = A.LockedLocalMatrix();
    }
    else
    {
        const Grid& g = this->GetGrid();
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [MR,MC] <- [MR,MC]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int row = g.MCRank();
        const int col = g.MRRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;
        const int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;
        const int sendRank = sendRow + sendCol*r;
        const int recvRank = recvRow + recvCol*r;

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = localHeightOfA * localWidthOfA;
        const int recvSize = localHeight * localWidth;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.VCComm() );

        // Unpack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.ColAlignment();
            this->_colShift = 
                Shift( g.MRRank(), this->ColAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() )
    {
        const int r = g.Height();
        const int rowShift = this->RowShift();

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = A.LocalEntry(i,rowShift+j*r);
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [MR,MC] <- [MR,* ]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int rowShift = this->RowShift();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (col+c+colAlignment-colAlignmentOfA) % c;
        const int recvRank = (col+c+colAlignmentOfA-colAlignment) % c;

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int sendSize = localHeightOfA * localWidth;
        const int recvSize = localHeight * localWidth;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i+j*localHeightOfA] = A.LocalEntry(i,rowShift+j*r);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.MRComm() );

        // Unpack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_rowShift = 
                Shift( g.MCRank(), this->RowAlignment(), g.Height() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const int c = g.Width();
        const int colShift = this->ColShift();

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = A.LocalEntry(colShift+i*c,j);
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [MR,MC] <- [* ,MC]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int row = g.MCRank(); 

        const int colShift = this->ColShift();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = localHeight * localWidthOfA;
        const int recvSize = localHeight * localWidth;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeight; ++i )
                sendBuffer[i+j*localHeight] = A.LocalEntry(colShift+i*c,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRow, 0,
          recvBuffer, recvSize, recvRow, MPI_ANY_TAG, g.MCComm() );

        // Unpack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,VR,Star> A_VR_Star(g);

    A_VR_Star = A;
    *this     = A_VR_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment() % g.Height();
            this->_rowShift = 
                Shift( g.MCRank(), this->RowAlignment(), g.Height() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() % g.Height() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int row = g.MCRank();

        const int rowShift = this->RowShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int maxHeight = MaxLocalLength(height,c);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        this->_auxMemory.Require( 2*c*portionSize );

        T* buffer = this->_auxMemory.Buffer();
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
        AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, g.MRComm() );

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
                    this->LocalEntry(i,thisRowOffset+j*c) = 
                        data[i+j*localHeight];
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [MR,MC] <- [* ,VC]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int row = g.MCRank();

        const int rowShift = this->RowShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+rowAlignment-(rowAlignmentOfA%r)) % r;
        const int recvRow = (row+r+(rowAlignmentOfA%r)-rowAlignment) % r;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int maxHeight = MaxLocalLength(height,c);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        this->_auxMemory.Require( 2*c*portionSize );

        T* buffer = this->_auxMemory.Buffer();
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
          firstBuffer,  c*portionSize, recvRow, MPI_ANY_TAG, g.MCComm() );

        // AllToAll to gather all of the aligned [* ,VC] into secondBuffer
        AllToAll
        ( firstBuffer,  portionSize, 
          secondBuffer, portionSize, g.MRComm() );

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
                    this->LocalEntry(i,thisRowOffset+j*c) = 
                        data[i+j*localHeight];
        }

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.ColAlignment() % g.Width();
            this->_colShift = 
                Shift( g.MRRank(), this->ColAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() % g.Width() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();

        const int colShift = this->ColShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,r);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        this->_auxMemory.Require( 2*r*portionSize );

        T* buffer = this->_auxMemory.Buffer();
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
        AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, g.MCComm() );

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
                    this->LocalEntry(thisColOffset+i*r,j) = 
                          data[i+j*thisLocalHeight];
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [MR,MC] <- [* ,VC]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();

        const int colShift = this->ColShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendCol = (col+c+colAlignment-(colAlignmentOfA%c)) % c;
        const int recvCol = (col+c+(colAlignmentOfA%c)-colAlignment) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,r);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        this->_auxMemory.Require( 2*r*portionSize );

        T* buffer = this->_auxMemory.Buffer();
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
          firstBuffer,  r*portionSize, recvCol, MPI_ANY_TAG, g.MRComm() );

        // AllToAll to gather all of the aligned [VR,* ] data into secondBuffer
        AllToAll
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

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
                    this->LocalEntry(thisColOffset+i*r,j) = 
                        data[i+j*thisLocalHeight];
        }

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,Star,VC> A_Star_VC(g);

    A_Star_VC = A;
    *this = A_Star_VC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MR,MC>&
elemental::DistMatrixBase<T,MR,MC>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = A.LocalEntry(colShift+i*c,rowShift+j*r);
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::SumScatterFrom
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterFrom([MR,* ])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.ColAlignment();
            this->_colShift = 
                Shift( g.MRRank(), this->ColAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() )
    {
        if( this->Width() == 1 )
        {
            const int rowAlignment = this->RowAlignment();
            const int myRow = g.MCRank();

            const int localHeight = this->LocalHeight();

            const int portionSize = max(localHeight,MinCollectContrib);

            this->_auxMemory.Require( 2*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack 
            for( int i=0; i<localHeight; ++i )
                sendBuffer[i] = A.LocalEntry(i,0);

            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize,
              MPI_SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
            {
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,0) = recvBuffer[i];
            }

            this->_auxMemory.Release();
        }
        else
        {
            const int r = g.Height();
            const int rowAlignment = this->RowAlignment();

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int recvSize = 
                max(localHeight*maxLocalWidth,MinCollectContrib);
            const int sendSize = r * recvSize;

            this->_auxMemory.Require( sendSize + recvSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];

            // Pack 
            vector<int> recvSizes(r);
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
            ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

            // Unpack our received data
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];

            this->_auxMemory.Release();
        }
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned SumScatterFrom [MR,MC] <- [MR,* ]." << endl;
#endif
        if( this->Width() == 1 )
        {
            const int c = g.Width();
            const int rowAlignment = this->RowAlignment();
            const int myRow = g.MCRank();
            const int myCol = g.MRRank();

            const int height = this->Height();
            const int localHeight = this->LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,c);

            const int portionSize = max(maxLocalHeight,MinCollectContrib);

            const int colAlignment = this->ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (myCol+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (myCol+c+colAlignmentOfA-colAlignment) % c;

            this->_auxMemory.Require( 2*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i] = A.LocalEntry(i,0);

            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize,
              MPI_SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
            {
                // Perform the realignment
                SendRecv
                ( recvBuffer, portionSize, sendCol, 0,
                  sendBuffer, portionSize, recvCol, 0, g.MRComm() );

                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,0) = sendBuffer[i];
            }

            this->_auxMemory.Release();
        }
        else
        {
            const int r = g.Height();
            const int c = g.Width();
            const int col = g.MRRank();

            const int colAlignment = this->ColAlignment();
            const int rowAlignment = this->RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int recvSize_RS = 
                max(localHeightOfA*maxLocalWidth,MinCollectContrib);
            const int sendSize_RS = r * recvSize_RS;
            const int recvSize_SR = localHeight * localWidth;

            this->_auxMemory.Require
            ( recvSize_RS + max(sendSize_RS,recvSize_SR) );

            T* buffer = this->_auxMemory.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[recvSize_RS];

            // Pack
            vector<int> recvSizes(r);
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
            ( secondBuffer, firstBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

            // Trade reduced data with the appropriate process col
            SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendCol, 0,
              secondBuffer, localHeight*localWidth,    recvCol, MPI_ANY_TAG,
              g.MRComm() );

            // Unpack the received data
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,j) = secondBuffer[i+j*localHeight];

            this->_auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::SumScatterFrom
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterFrom([* ,MC])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
    if( A.GetGrid().VCRank() == 0 )
    {
        if( A.Height() == 1 )
        {
            cout <<    
              "The vector version of [MR,MC].SumScatterFrom([* ,MC]) is not yet"
              " written, but it only requires a modification of the vector "
              "version of [MR,MC].SumScatterFrom([MR,* ])." << endl;
        }
        else
        {
            cout << 
              "[MR,MC]::SumScatterFrom([* ,MC]) necessarily causes a large "
              "amount of cache-thrashing. If possible, avoid it by forming the "
              "(conjugate-)transpose of the [* ,MC] matrix instead." << endl;
        }
    }
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_rowShift = 
                Shift( g.MCRank(), this->RowAlignment(), g.Height() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const int c = g.Width();
        const int colAlignment = this->ColAlignment();

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);

        const int recvSize = max(maxLocalHeight*localWidth,MinCollectContrib);
        const int sendSize = c * recvSize;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack 
        vector<int> recvSizes(c);
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
        ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MRComm() );

        // Unpack our received data
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];

        this->_auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned SumScatterFrom [MR,MC] <- [* ,MC]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int row = g.MCRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);
        
        const int recvSize_RS = 
            max(maxLocalHeight*localWidthOfA,MinCollectContrib);
        const int sendSize_RS = c* recvSize_RS;
        const int recvSize_SR = localHeight * localWidth;

        this->_auxMemory.Require( recvSize_RS + max(sendSize_RS,recvSize_SR) );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[recvSize_RS];

        // Pack 
        vector<int> recvSizes(c);
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
        ( secondBuffer, firstBuffer, &recvSizes[0], MPI_SUM, g.MRComm() );

        // Trade reduced data with the appropriate process row
        SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendRow, 0,
          secondBuffer, localHeight*localWidth,    recvRow, MPI_ANY_TAG, 
          g.MCComm() );

        // Unpack the received data
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = secondBuffer[i+j*localHeight];

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::SumScatterFrom
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterFrom([* ,* ])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int colAlignment = this->ColAlignment();
    const int rowAlignment = this->RowAlignment();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int maxLocalHeight = MaxLocalLength(height,c);
    const int maxLocalWidth = MaxLocalLength(width,r);

    const int recvSize = max(maxLocalHeight*maxLocalWidth,MinCollectContrib);
    const int sendSize = r * c * recvSize;

    this->_auxMemory.Require( sendSize + recvSize );

    T* buffer = this->_auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack 
    vector<int> recvSizes(r*c);
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
    ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.VRComm() );

    // Unpack our received data
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::SumScatterUpdate
( T alpha, const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterUpdate([MR,* ])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A ); 
#endif
    const Grid& g = this->GetGrid();
    if( this->ColAlignment() == A.ColAlignment() )
    {
        if( this->Width() == 1 )
        {
            const int rowAlignment = this->RowAlignment();
            const int myRow = g.MCRank();

            const int localHeight = this->LocalHeight();

            const int portionSize = max(localHeight,MinCollectContrib);

            this->_auxMemory.Require( 2*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack
            for( int i=0; i<localHeight; ++i )
                sendBuffer[i] = A.LocalEntry(i,0);

            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize,
              MPI_SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
            {
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,0) += alpha*recvBuffer[i];
            }

            this->_auxMemory.Release();
        }
        else
        {
            const int r = g.Height();
            const int rowAlignment = this->RowAlignment();

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int portionSize = 
                max(localHeight*maxLocalWidth,MinCollectContrib);

            this->_auxMemory.Require( (r+1)*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack 
            vector<int> recvSizes(r);
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
            ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

            // Update with our received data
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,j) += alpha*recvBuffer[i+j*localHeight];

            this->_auxMemory.Release();
        }
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned SumScatterUpdate [MR,MC] <- [MR,* ]." << endl;
#endif
        if( this->Width() == 1 )
        {
            const int c = g.Width();
            const int rowAlignment = this->RowAlignment();
            const int myRow = g.MCRank();
            const int myCol = g.MRRank();

            const int height = this->Height();
            const int localHeight = this->LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,c);

            const int portionSize = max(maxLocalHeight,MinCollectContrib);

            const int colAlignment = this->ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (myCol+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (myCol+c+colAlignmentOfA-colAlignment) % c;

            this->_auxMemory.Require( 2*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i] = A.LocalEntry(i,0);

            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize,
              MPI_SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
            {
                // Perform the realignment
                SendRecv
                ( recvBuffer, portionSize, sendCol, 0,
                  sendBuffer, portionSize, recvCol, 0, g.MRComm() );

                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,0) += alpha*sendBuffer[i];
            }

            this->_auxMemory.Release();
        }
        else
        {
            const int r = g.Height();
            const int c = g.Width();
            const int col = g.MRRank();

            const int colAlignment = this->ColAlignment();
            const int rowAlignment = this->RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendCol = (col+c+colAlignment-colAlignmentOfA) % c;
            const int recvCol = (col+c+colAlignmentOfA-colAlignment) % c;

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalWidth = MaxLocalLength(width,r);

            const int recvSize_RS = 
                max(localHeightOfA*maxLocalWidth,MinCollectContrib);
            const int sendSize_RS = r * recvSize_RS;
            const int recvSize_SR = localHeight * localWidth;

            this->_auxMemory.Require
            ( recvSize_RS + max(sendSize_RS,recvSize_SR) );

            T* buffer = this->_auxMemory.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[recvSize_RS];

            // Pack 
            vector<int> recvSizes(r);
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
            ( secondBuffer, firstBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

            // Trade reduced data with the appropriate process col
            SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendCol, 0, 
              secondBuffer, localHeight*localWidth,    recvCol, MPI_ANY_TAG,
              g.MRComm() );

            // Update with our received data
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,j) += 
                        alpha*secondBuffer[i+j*localHeight];

            this->_auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::SumScatterUpdate
( T alpha, const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterUpdate([* ,MC])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
    if( A.GetGrid().VCRank() == 0 )
    {
        if( A.Height() == 1 )
        {
            cout <<    
              "The vector version of [MR,MC].SumScatterUpdate([* ,MC]) is not "
              "yet written, but it only requires a modification of the vector "
              "version of [MR,MC].SumScatterUpdate([MR,* ])." << endl;
        }
        else
        {
            cout << 
              "[MR,MC]::SumScatterUpdate([* ,MC]) necessarily causes a large "
              "amount of cache-thrashing. If possible, avoid it by forming the "
              "(conjugate-)transpose of the [* ,MC] matrix instead." << endl;
        }
    }
#endif
    const Grid& g = this->GetGrid();
    if( this->RowAlignment() == A.RowAlignment() )
    {
        const int c = g.Width();
        const int colAlignment = this->ColAlignment();

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);

        const int recvSize = max(maxLocalHeight*localWidth,MinCollectContrib);
        const int sendSize = c * recvSize;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        vector<int> recvSizes(c);
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
        ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MRComm() );

        // Update with our received data
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) += alpha*recvBuffer[i+j*localHeight];

        this->_auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned SumScatterUpdate [MR,MC] <- [* ,MC]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int row = g.MCRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendRow = (row+r+rowAlignment-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-rowAlignment) % r;

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);

        const int recvSize_RS = 
            max(maxLocalHeight*localWidthOfA,MinCollectContrib);
        const int sendSize_RS = c * recvSize_RS;
        const int recvSize_SR = localHeight * localWidth;

        this->_auxMemory.Require( recvSize_RS + max(sendSize_RS,recvSize_SR) );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[recvSize_RS];

        // Pack 
        vector<int> recvSizes(c);
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
        ( secondBuffer, firstBuffer, &recvSizes[0], MPI_SUM, g.MRComm() );

        // Trade reduced data with the appropriate process row
        SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendRow, 0,
          secondBuffer, localHeight*localWidth,    recvRow, MPI_ANY_TAG,
          g.MRComm() );

        // Update with our received data
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) += alpha*secondBuffer[i+j*localHeight];

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::SumScatterUpdate
( T alpha, const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterUpdate([* ,* ])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int colAlignment = this->ColAlignment();
    const int rowAlignment = this->RowAlignment();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int maxLocalHeight = MaxLocalLength(height,c);
    const int maxLocalWidth = MaxLocalLength(width,r);

    const int recvSize = max(maxLocalHeight*maxLocalWidth,MinCollectContrib);
    const int sendSize = r * c * recvSize;

    this->_auxMemory.Require( sendSize + recvSize );

    T* buffer = this->_auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack 
    vector<int> recvSizes(r*c);
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
    ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.VRComm() );

    // Unpack our received data
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) += alpha*recvBuffer[i+j*localHeight];

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// DistMatrix                                                                 //
//----------------------------------------------------------------------------//

template<typename R>
void
elemental::DistMatrix<R,MR,MC>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const Grid& g = this->GetGrid();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToRandom();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*c;
        if( i % r == rowShift )
        {
            const int jLoc = (i-rowShift) / r;
            if( jLoc < localWidth )
                this->LocalEntry(iLoc,jLoc) += (R)this->Width();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const Grid& g = this->GetGrid();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToRandom();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*c;
        if( i % r == rowShift )
        {
            const int jLoc = (i-rowShift) / r;
            if( jLoc < localWidth )
                this->LocalEntry(iLoc,jLoc) = 
                    real(this->LocalEntry(iLoc,jLoc)) + (R)this->Width();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,MR,MC>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const Grid& g = this->GetGrid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        u = real(this->LocalEntry(iLoc,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,MR,MC>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const Grid& g = this->GetGrid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        u = imag(this->LocalEntry(iLoc,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        real(this->LocalEntry(iLoc,jLoc)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::SetImag
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        imag(this->LocalEntry(iLoc,jLoc)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::GetRealDiagonal
( DistMatrix<R,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetRealDiagonal([MD,* ])");
    this->AssertNotLockedView();
#endif
    int height = this->Height();
    int width = this->Width();
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
        throw logic_error( "d is not of the correct length." );
#endif
    if( !d.Viewing() )
    {
        if( !d.ConstrainedColAlignment() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( length, 1 );
    }

    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            d.LocalEntry(k,0) = 
                real(this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r)));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::GetImagDiagonal
( DistMatrix<R,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetImagDiagonal([MD,* ])");
    this->AssertNotLockedView();
#endif
    int height = this->Height();
    int width = this->Width();
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
        throw logic_error( "d is not of the correct length." );
#endif
    if( !d.Viewing() )
    {
        if( !d.ConstrainedColAlignment() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( length, 1 );
    }

    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            d.LocalEntry(k,0) = 
                imag(this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r)));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::GetRealDiagonal
( DistMatrix<R,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetRealDiagonal([* ,MD])");
    this->AssertNotLockedView();
#endif
    int height = this->Height();
    int width = this->Width();
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
        throw logic_error( "d is not of the correct length." );
#endif

    if( !d.Viewing() )
    {
        if( !d.ConstrainedRowAlignment() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( 1, length );
    }

    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            d.LocalEntry(0,k) = 
                real(this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r)));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::GetImagDiagonal
( DistMatrix<R,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetImagDiagonal([* ,MD])");
    this->AssertNotLockedView();
#endif
    int height = this->Height();
    int width = this->Width();
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
        throw logic_error( "d is not of the correct length." );
#endif

    if( !d.Viewing() )
    {
        if( !d.ConstrainedRowAlignment() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( 1, length );
    }

    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            d.LocalEntry(0,k) = 
                imag(this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r)));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::SetDiagonal
( const DistMatrixBase<R,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetDiagonal([MD,* ])");
    if( d.Width() != 1 )
        throw logic_error( "d must be a column vector." );
    {
        int height = this->Height();
        int width = this->Width();
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
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r)) = 
                d.LocalEntry(k,0);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::SetRealDiagonal
( const DistMatrixBase<R,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetRealDiagonal([MD,* ])");
    if( d.Width() != 1 )
        throw logic_error( "d must be a column vector." );
    {
        int height = this->Height();
        int width = this->Width();
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
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            real(this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r))) = 
                d.LocalEntry(k,0);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::SetImagDiagonal
( const DistMatrixBase<R,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetImagDiagonal([MD,* ])");
    if( d.Width() != 1 )
        throw logic_error( "d must be a column vector." );
    {
        int height = this->Height();
        int width = this->Width();
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
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            imag(this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r))) = 
                d.LocalEntry(k,0);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::SetDiagonal
( const DistMatrixBase<R,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetDiagonal([* ,MD])");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        int height = this->Height();
        int width = this->Width();
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
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r)) = 
                d.LocalEntry(0,k);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::SetRealDiagonal
( const DistMatrixBase<R,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetRealDiagonal([* ,MD])");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        int height = this->Height();
        int width = this->Width();
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
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            real(this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r))) = 
                d.LocalEntry(0,k);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MR,MC>::SetImagDiagonal
( const DistMatrixBase<R,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetImagDiagonal([* ,MD])");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        int height = this->Height();
        int width = this->Width();
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
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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
        for( int k=0; k<localDiagLength; ++k )
            imag(this->LocalEntry(iLocStart+k*(lcm/c),jLocStart+k*(lcm/r))) = 
                d.LocalEntry(0,k);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrixBase<int,   MR,MC>;
template class elemental::DistMatrixBase<float, MR,MC>;
template class elemental::DistMatrixBase<double,MR,MC>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,MR,MC>;
template class elemental::DistMatrixBase<dcomplex,MR,MC>;
#endif

template class elemental::DistMatrix<int,   MR,MC>;
template class elemental::DistMatrix<float, MR,MC>;
template class elemental::DistMatrix<double,MR,MC>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,MR,MC>;
template class elemental::DistMatrix<dcomplex,MR,MC>;
#endif

