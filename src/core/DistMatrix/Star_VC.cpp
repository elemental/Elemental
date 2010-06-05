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
elemental::DistMatrixBase<T,Star,VC>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Print");
#endif
    const Grid& grid = this->GetGrid();
    if( grid.VCRank() == 0 && s != "" )
        cout << s << endl;

    const int height     = this->Height();
    const int width      = this->Width();
    const int localWidth = this->LocalWidth();
    const int p          = grid.Size();
    const int rowShift   = this->RowShift();

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
            sendBuf[i+(rowShift+j*p)*height] = this->LocalEntry(i,j);

    // If we are the root, fill the receive buffer
    T* recvBuf = 0;
    if( grid.VCRank() == 0 )
    {
        recvBuf = new T[height*width];
        for( int i=0; i<height*width; ++i )
            recvBuf[i] = (T)0;
    }

    // Sum the contributions and send to the root
    Reduce( sendBuf, recvBuf, height*width, MPI_SUM, 0, grid.VCComm() );
    delete[] sendBuf;

    if( grid.VCRank() == 0 )
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
    Barrier( grid.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::AlignWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& grid = this->GetGrid();
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = Shift( grid.VCRank(), this->RowAlignment(), grid.Size() );
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::AlignWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& grid = this->GetGrid();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = Shift( grid.VCRank(), this->RowAlignment(), grid.Size() );
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::AlignWith
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([MC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& grid = this->GetGrid();
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = Shift( grid.VCRank(), this->RowAlignment(), grid.Size() );
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::AlignWith
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([* ,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& grid = this->GetGrid();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = Shift( grid.VCRank(), this->RowAlignment(), grid.Size() );
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::AlignWith
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([* ,VC])");
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
elemental::DistMatrixBase<T,Star,VC>::AlignWith
( const DistMatrixBase<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::AlignWith([VC,* ])");
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
elemental::DistMatrixBase<T,Star,VC>::AlignRowsWith
( const DistMatrixBase<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::AlignRowsWith
( const DistMatrixBase<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::AlignRowsWith
( const DistMatrixBase<T,MC,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::AlignRowsWith
( const DistMatrixBase<T,Star,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::AlignRowsWith
( const DistMatrixBase<T,Star,VC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::AlignRowsWith
( const DistMatrixBase<T,VC,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::View
( DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = A.RowShift();
    this->_localMatrix.View( A.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::LockedView
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = A.RowShift();
    this->_localMatrix.LockedView( A.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::View
( DistMatrixBase<T,Star,VC>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const Grid& grid = this->GetGrid();
        const int colMajorRank = grid.VCRank();
        const int size = grid.Size();

        this->_rowAlignment = (A.RowAlignment()+j) % size;
        this->_rowShift = Shift( colMajorRank, this->RowAlignment(), size );

        const int localWidthBefore = LocalLength( j, A.RowShift(), size );
        const int localWidth = LocalLength( width, this->RowShift(), size );

        this->_localMatrix.View
        ( A.LocalMatrix(), i, localWidthBefore, height, localWidth );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::LockedView
( const DistMatrixBase<T,Star,VC>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const Grid& grid = this->GetGrid();
        const int colMajorRank = grid.VCRank();
        const int size = grid.Size();

        this->_rowAlignment = (A.RowAlignment()+j) % size;
        this->_rowShift = Shift( colMajorRank, this->RowAlignment(), size );

        const int localWidth = LocalLength( width, this->RowShift(), size );

        this->_localMatrix.LockedView
        ( A.LockedLocalMatrix(), i, A.RowShift(), height, localWidth );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::View1x2
( DistMatrixBase<T,Star,VC>& AL,
  DistMatrixBase<T,Star,VC>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View1x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_rowAlignment = AL.RowAlignment();
    this->_rowShift = AL.RowShift();
    this->_localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::LockedView1x2
( const DistMatrixBase<T,Star,VC>& AL,
  const DistMatrixBase<T,Star,VC>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView1x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_rowAlignment = AL.RowAlignment();
    this->_rowShift = AL.RowShift();
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
elemental::DistMatrixBase<T,Star,VC>::View2x1
( DistMatrixBase<T,Star,VC>& AT,
  DistMatrixBase<T,Star,VC>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View2x1");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_rowAlignment = AT.RowAlignment();
    this->_rowShift = AT.RowShift();
    this->_localMatrix.View2x1( AT.LocalMatrix(), AB.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::LockedView2x1
( const DistMatrixBase<T,Star,VC>& AT,
  const DistMatrixBase<T,Star,VC>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView2x1");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_rowAlignment = AT.RowAlignment();
    this->_rowShift = AT.RowShift();
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
elemental::DistMatrixBase<T,Star,VC>::View2x2
( DistMatrixBase<T,Star,VC>& ATL, 
  DistMatrixBase<T,Star,VC>& ATR,
  DistMatrixBase<T,Star,VC>& ABL,
  DistMatrixBase<T,Star,VC>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::View2x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_rowAlignment = ATL.RowAlignment();
    this->_rowShift = ATL.RowShift();
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
elemental::DistMatrixBase<T,Star,VC>::LockedView2x2
( const DistMatrixBase<T,Star,VC>& ATL, 
  const DistMatrixBase<T,Star,VC>& ATR,
  const DistMatrixBase<T,Star,VC>& ABL,
  const DistMatrixBase<T,Star,VC>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::LockedView2x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_rowAlignment = ATL.RowAlignment();
    this->_rowShift = ATL.RowShift();
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
elemental::DistMatrixBase<T,Star,VC>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    const Grid& grid = this->GetGrid();
    this->_height = height;
    this->_width  = width;
    this->_localMatrix.ResizeTo
    ( height, LocalLength( width, this->RowShift(), grid.Size() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,Star,VC>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire grid
    const Grid& grid = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % grid.Size();

    T u;
    if( grid.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / grid.Size();
        u = this->LocalEntry(i,jLoc);
    }
    Broadcast( &u, 1, ownerRank, grid.VCComm() );
    
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const Grid& grid = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % grid.Size();

    if( grid.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / grid.Size();
        this->LocalEntry(i,jLoc) = u;
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
elemental::DistMatrixBase<T,Star,VC>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localWidth = this->LocalWidth();
    const int p = this->GetGrid().Size();
    const int rowShift = this->RowShift();

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
                this->LocalEntry(i,jLoc) = (T)0;
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
                this->LocalEntry(i,jLoc) = (T)0;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int p = this->GetGrid().Size();
    const int rowShift = this->RowShift();

    this->SetToZero();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*p;
        if( j < height )
            this->LocalEntry(j,jLoc) = (T)1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,VC>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<height; ++i )
            this->LocalEntry(i,j) = Random<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    DistMatrix<T,Star,VR> A_Star_VR(grid);

    A_Star_VR = A;
    *this = A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(grid) );
    *A_MC_MR = A;

    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
    ( new DistMatrix<T,Star,VR>(grid) );
    *A_Star_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    DistMatrix<T,Star,VR> A_Star_VR(grid);

    A_Star_VR = A;
    *this = A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw "[* ,VC] = [MD,* ] not yet implemented.";
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw "[* ,VC] = [* ,MD] not yet implemented.";
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_rowShift = 
                Shift( grid.VCRank(), this->RowAlignment(), grid.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() % grid.Height() == A.RowAlignment() )
    {
        const int r = grid.Height();
        const int c = grid.Width();
        const int p = grid.Size();
        const int row = grid.MCRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

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
        AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, grid.MRComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignmentOfA,c);
            const int thisLocalHeight = LocalLength(height,thisColShift,c);

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    this->LocalEntry(thisColShift+i*c,j) = 
                        data[i+j*thisLocalHeight];
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( grid.VCRank() == 0 )
            cout << "Unaligned [* ,VC] <- [MR,MC]." << endl;
#endif
        const int r = grid.Height();
        const int c = grid.Width();
        const int p = grid.Size();
        const int row = grid.MCRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+(rowAlignment%r)-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-(rowAlignment%r)) % r;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

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
          firstBuffer,  portionSize, grid.MRComm() );

        // SendRecv: properly align the [*,VC] via a trade in the column
        SendRecv
        ( firstBuffer,  portionSize, sendRow, 0,
          secondBuffer, portionSize, recvRow, MPI_ANY_TAG, grid.MCComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignmentOfA,c);
            const int thisLocalHeight = LocalLength(height,thisColShift,c);

            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    this->LocalEntry(thisColShift+i*c,j) = 
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
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    DistMatrix<T,MR,MC> A_MR_MC(grid);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_rowShift = 
                Shift( grid.VCRank(), this->RowAlignment(), grid.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() % grid.Height() == A.RowAlignment() )
    {
        const int r = grid.Height();
        const int c = grid.Width();
        const int rowShift = this->RowShift();
        const int rowShiftOfA = A.RowShift();
        const int rowOffset = (rowShift-rowShiftOfA) / r;

        const int height = this->Height();
        const int localWidth = this->LocalWidth();

        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->LocalEntry(i,j) = A.LocalEntry(i,rowOffset+j*c);
    }
    else
    {
#ifndef RELEASE
        if( grid.VCRank() == 0 )
            cout << "Unaligned [* ,VC] <- [* ,MC]." << endl;
#endif
        const int r = grid.Height();
        const int c = grid.Width();
        const int p = grid.Size();
        const int row = grid.MCRank();
        const int col = grid.MRRank();
        const int rowShiftOfA = A.RowShift();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        // We will SendRecv A[*,VC] within our process column to fix alignments.
        const int sendRow = (row+r+(rowAlignment%r)-rowAlignmentOfA) % r;
        const int recvRow = (row+r+rowAlignmentOfA-(rowAlignment%r)) % r;
        const int sendRank = sendRow + r*col;

        const int sendRowShift = Shift( sendRank, rowAlignment, p );
        const int sendRowOffset = (sendRowShift-rowShiftOfA) / r;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localWidthOfSend = LocalLength(width,sendRowShift,p);

        const int sendSize = height * localWidthOfSend;
        const int recvSize = height * localWidth;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidthOfSend; ++j )
            for( int i=0; i<height; ++i )
                sendBuffer[i+j*height] = A.LocalEntry(i,sendRowOffset+j*c);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRow, 0,
          recvBuffer, recvSize, recvRow, MPI_ANY_TAG, grid.MCComm() );

        // Unpack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->LocalEntry(i,j) = recvBuffer[i+j*height];

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(grid) );
    *A_MC_MR = A;

    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
    ( new DistMatrix<T,Star,VR>(grid) );
    *A_Star_VR = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid(); 
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_rowShift = A.RowShift();
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        this->_localMatrix = A.LockedLocalMatrix();
    }
    else
    {
#ifndef RELEASE
        if( grid.VCRank() == 0 )
            cout << "Unaligned [* ,VC] <- [* ,VC]." << endl;
#endif
        const int rank = grid.VCRank();
        const int p = grid.Size();

        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+p+rowAlignment-rowAlignmentOfA) % p;
        const int recvRank = (rank+p+rowAlignmentOfA-rowAlignment) % p;

        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = height * localWidthOfA;
        const int recvSize = height * localWidth;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<height; ++i )
                sendBuffer[i+j*height] = A.LocalEntry(i,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, grid.VCComm() );

        // Unpack
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->LocalEntry(i,j) = recvBuffer[i+j*height];

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    DistMatrix<T,MR,MC> A_MR_MC(grid);

    A_MR_MC = A;
    *this = A_MR_MC;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    
    const int height = this->Height();
    const int localWidth = this->LocalWidth();
    const int localWidthOfA = A.LocalWidth();
    
    const int sendSize = height * localWidthOfA;
    const int recvSize = height * localWidth;

    const Grid& grid = this->GetGrid();
    const int r = grid.Height();
    const int c = grid.Width();
    const int p = grid.Size();
    const int rankCM = grid.VCRank();
    const int rankRM = grid.VRRank(); 

    const int rowShift = this->RowShift();
    const int rowShiftOfA = A.RowShift();

    // Compute which colmajor rank has the rowShift equal to our rowShiftOfA
    const int sendRankCM = (rankCM+(p+rowShiftOfA-rowShift)) % p;

    // Compute which colmajor rank has the A rowShift that we need
    const int recvRankRM = (rankRM+(p+rowShift-rowShiftOfA)) % p;
    const int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

    this->_auxMemory.Require( sendSize + recvSize );

    T* buffer = this->_auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack
    for( int j=0; j<localWidthOfA; ++j )
        for( int i=0; i<height; ++i )
            sendBuffer[i+j*height] = A.LocalEntry(i,j);

    // Communicate
    SendRecv
    ( sendBuffer, sendSize, sendRankCM, 0,
      recvBuffer, recvSize, recvRankCM, MPI_ANY_TAG, grid.VCComm() );

    // Unpack
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<height; ++i )
            this->LocalEntry(i,j) = recvBuffer[i+j*height];

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,VC>&
elemental::DistMatrixBase<T,Star,VC>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,VC] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int p = this->GetGrid().Size();
    const int rowShift = this->RowShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = A.LocalEntry(i,rowShift+j*p);
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

//----------------------------------------------------------------------------//
// DistMatrix                                                                 //
//----------------------------------------------------------------------------//

template<typename R>
void
elemental::DistMatrix<R,Star,VC>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw "Positive-definite matrices must be square.";
#endif
    const int height     = this->Height();
    const int localWidth = this->LocalWidth();
    const int p          = this->GetGrid().Size();
    const int rowShift   = this->RowShift();

    this->SetToRandom();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*p;
        if( j < height )
            this->LocalEntry(j,jLoc) += (R)this->Width();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,Star,VC>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw "Positive-definite matrices must be square.";
#endif
    const int height     = this->Height();
    const int localWidth = this->LocalWidth();
    const int p          = this->GetGrid().Size();
    const int rowShift   = this->RowShift();

    this->SetToRandom();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*p;
        if( j < height )
            this->LocalEntry(j,jLoc) = 
                real(this->LocalEntry(j,jLoc)) + (R)this->Width();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,VC>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire grid
    const Grid& grid = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % grid.Size();

    R u;
    if( grid.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / grid.Size();
        u = real(this->LocalEntry(i,jLoc));
    }
    Broadcast( &u, 1, ownerRank, grid.VCComm() );
    
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,VC>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire grid
    const Grid& grid = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % grid.Size();

    R u;
    if( grid.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / grid.Size();
        u = imag(this->LocalEntry(i,jLoc));
    }
    Broadcast( &u, 1, ownerRank, grid.VCComm() );
    
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,VC>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& grid = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % grid.Size();

    if( grid.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / grid.Size();
        real(this->LocalEntry(i,jLoc)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,VC>::SetImag
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& grid = this->GetGrid();
    const int ownerRank = (j + this->RowAlignment()) % grid.Size();

    if( grid.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / grid.Size();
        imag(this->LocalEntry(i,jLoc)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrixBase<int,   Star,VC>;
template class elemental::DistMatrixBase<float, Star,VC>;
template class elemental::DistMatrixBase<double,Star,VC>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,Star,VC>;
template class elemental::DistMatrixBase<dcomplex,Star,VC>;
#endif

template class elemental::DistMatrix<int,     Star,VC>;
template class elemental::DistMatrix<float,   Star,VC>;
template class elemental::DistMatrix<double,  Star,VC>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,Star,VC>;
template class elemental::DistMatrix<dcomplex,Star,VC>;
#endif

