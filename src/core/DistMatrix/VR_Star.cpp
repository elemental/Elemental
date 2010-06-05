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
elemental::DistMatrixBase<T,VR,Star>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Print");
#endif
    const Grid& grid = this->GetGrid();
    if( grid.VRRank() == 0 && s != "" )
        cout << s << endl;

    const int height      = this->Height();
    const int width       = this->Width();
    const int localHeight = this->LocalHeight();
    const int p           = grid.Size();
    const int colShift    = this->ColShift();

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
            sendBuf[colShift+i*p+j*height] = this->LocalEntry(i,j);

    // If we are the root, fill the receive buffer
    T* recvBuf = 0;
    if( grid.VRRank() == 0 )
    {
        recvBuf = new T[height*width];
        for( int i=0; i<height*width; ++i )
            recvBuf[i] = (T)0;
    }

    // Sum the contributions and send to the root
    Reduce( sendBuf, recvBuf, height*width, MPI_SUM, 0, grid.VCComm() );
    delete[] sendBuf;

    if( grid.VRRank() == 0 )
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
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& grid = this->GetGrid();
    this->_colAlignment = A.RowAlignment();
    this->_colShift = Shift( grid.VRRank(), this->ColAlignment(), grid.Size() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& grid = this->GetGrid();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = Shift( grid.VRRank(), this->ColAlignment(), grid.Size() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([MR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& grid = this->GetGrid();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = Shift( grid.VRRank(), this->ColAlignment(), grid.Size() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith([* ,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& grid = this->GetGrid();
    this->_colAlignment = A.RowAlignment();
    this->_colShift = Shift( grid.VRRank(), this->ColAlignment(), grid.Size() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith(DistMatrix[VR,* ])");
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
elemental::DistMatrixBase<T,VR,Star>::AlignWith
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::AlignWith(DistMatrix[* ,VR])");
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
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,MR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,Star,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,VR,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::AlignColsWith
( const DistMatrixBase<T,Star,VR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::View
( DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View(A)");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = A.ColShift();
    this->_localMatrix.View( A.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::LockedView
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView(A)");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = A.ColShift();
    this->_localMatrix.LockedView( A.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::View
( DistMatrixBase<T,VR,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View(A,i,j,height,width)");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const Grid& grid = this->GetGrid();
        const int rowMajorRank = grid.VRRank();
        const int size = grid.Size();

        this->_colAlignment = (A.ColAlignment()+i) % size;
        this->_colShift = Shift( rowMajorRank, this->ColAlignment(), size );

        const int localHeightBefore = LocalLength( i, A.ColShift(), size );
        const int localHeight = LocalLength( height, this->ColShift(), size );

        this->_localMatrix.View
        ( A.LocalMatrix(), localHeightBefore, j, localHeight, width );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::LockedView
( const DistMatrixBase<T,VR,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const Grid& grid = this->GetGrid();
        const int rowMajorRank = grid.VRRank();
        const int size = grid.Size();

        this->_colAlignment = (A.ColAlignment()+i) % size;
        this->_colShift = Shift( rowMajorRank, this->ColAlignment(), size );

        const int localHeightBefore = LocalLength( i, A.ColShift(), size );
        const int localHeight = LocalLength( height, this->ColShift(), size );

        this->_localMatrix.LockedView
        ( A.LockedLocalMatrix(), localHeightBefore, j, localHeight, width );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::View1x2
( DistMatrixBase<T,VR,Star>& AL,
  DistMatrixBase<T,VR,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View1x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_colAlignment = AL.ColAlignment();
    this->_colShift = AL.ColShift();
    this->_localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::LockedView1x2
( const DistMatrixBase<T,VR,Star>& AL,
  const DistMatrixBase<T,VR,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView1x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_colAlignment = AL.ColAlignment();
    this->_colShift = AL.ColShift();
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
elemental::DistMatrixBase<T,VR,Star>::View2x1
( DistMatrixBase<T,VR,Star>& AT,
  DistMatrixBase<T,VR,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_colAlignment = AT.ColAlignment();
    this->_colShift = AT.ColShift();
    this->_localMatrix.View2x1( AT.LocalMatrix(), AB.LocalMatrix() );
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::LockedView2x1
( const DistMatrixBase<T,VR,Star>& AT,
  const DistMatrixBase<T,VR,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_colAlignment = AT.ColAlignment();
    this->_colShift = AT.ColShift();
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
elemental::DistMatrixBase<T,VR,Star>::View2x2
( DistMatrixBase<T,VR,Star>& ATL,
  DistMatrixBase<T,VR,Star>& ATR,
  DistMatrixBase<T,VR,Star>& ABL,
  DistMatrixBase<T,VR,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::View2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ABL.Width() + ABR.Width();
    this->_colAlignment = ATL.ColAlignment();
    this->_colShift = ATL.ColShift();
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
elemental::DistMatrixBase<T,VR,Star>::LockedView2x2
( const DistMatrixBase<T,VR,Star>& ATL,
  const DistMatrixBase<T,VR,Star>& ATR,
  const DistMatrixBase<T,VR,Star>& ABL,
  const DistMatrixBase<T,VR,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::LockedView2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ABL.Width() + ABR.Width();
    this->_colAlignment = ATL.ColAlignment();
    this->_colShift = ATL.ColShift();
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
elemental::DistMatrixBase<T,VR,Star>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    const Grid& grid = this->GetGrid();
    this->_height = height;
    this->_width  = width;
    this->_localMatrix.ResizeTo
    ( LocalLength(height,this->ColShift(),grid.Size()) ,width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,VR,Star>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire grid
    const Grid& grid = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % grid.Size();

    T u;
    if( grid.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / grid.Size();
        u = this->LocalEntry(iLoc,j);
    }
    Broadcast( &u, 1, ownerRank, grid.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const Grid& grid = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % grid.Size();

    if( grid.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / grid.Size();
        this->LocalEntry(iLoc,j) = u;
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
elemental::DistMatrixBase<T,VR,Star>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->GetGrid().Size();
    const int colShift = this->ColShift();

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
                    this->LocalEntry(iLoc,j) = (T)0;
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
                this->LocalEntry(iLoc,j) = (T)0;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::SetToIdentity()
{   
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->GetGrid().Size();
    const int colShift = this->ColShift();

    this->SetToZero();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*p;
        if( i < width )
            this->LocalEntry(iLoc,i) = (T)1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,VR,Star>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = Random<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    DistMatrix<T,VC,Star> A_VC_Star(grid);

    A_VC_Star = A;
    *this = A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    DistMatrix<T,VC,Star> A_VC_Star(grid);

    A_VC_Star = A;
    *this = A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(grid) );
    *A_MC_MR = A;

    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(grid) );
    *A_VC_Star = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw "[VR,* ] = [MD,* ] not yet implemented.";
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw "[VR,* ] = [* ,MD] not yet implemented.";
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.ColAlignment();
            this->_colShift = 
                Shift( grid.VRRank(), this->ColAlignment(), grid.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() % grid.Width() == A.ColAlignment() )
    {
        const int r = grid.Height();
        const int c = grid.Width();
        const int p = grid.Size();
        const int col = grid.MRRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

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

            const int thisRank = col+k*c;
            const int thisColShift = Shift(thisRank,colAlignment,p);
            const int thisColOffset = (thisColShift-colShiftOfA) / c;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                        A.LocalEntry(thisColOffset+i*r,j);
        }

        // Communicate
        AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, grid.MCComm() );

        // Unpack
        for( int k=0; k<r; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignmentOfA,r);
            const int thisLocalWidth = LocalLength(width,thisRowShift,r);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,thisRowShift+j*r) = 
                        data[i+j*localHeight];
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( grid.VCRank() == 0 )
            cout << "Unaligned [VR,* ] <- [MR,MC]." << endl;
#endif
        const int r = grid.Height();
        const int c = grid.Width();
        const int p = grid.Size();
        const int col = grid.MRRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendCol = (col+c+(colAlignment%c)-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-(colAlignment%c)) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

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

            const int thisRank = sendCol+k*c;
            const int thisColShift = Shift(thisRank,colAlignment,p);
            const int thisColOffset = (thisColShift-colShiftOfA) / c;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                        A.LocalEntry(thisColOffset+i*r,j);
        }

        // AllToAll to gather all of the unaligned [VR,*] data into firstBuffer
        AllToAll
        ( secondBuffer, portionSize,
          firstBuffer,  portionSize, grid.MCComm() );

        // SendRecv: properly align the [VR,*] via a trade in the row
        SendRecv
        ( firstBuffer,  portionSize, sendCol, 0,
          secondBuffer, portionSize, recvCol, MPI_ANY_TAG, grid.MRComm() );

        // Unpack
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignmentOfA,r);
            const int thisLocalWidth = LocalLength(width,thisRowShift,r);

            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(i,thisRowShift+j*r) = 
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
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.ColAlignment();
            this->_colShift = 
                Shift( grid.VRRank(), this->ColAlignment(), grid.Size() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() % grid.Width() == A.ColAlignment() )
    {
        const int r = grid.Height();
        const int c = grid.Width();
        const int colShift = this->ColShift();
        const int colShiftOfA = A.ColShift();
        const int colOffset = (colShift-colShiftOfA) / c;

        const int width = this->Width();
        const int localHeight = this->LocalHeight();

        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = A.LocalEntry(colOffset+i*r,j);
    }
    else
    {
#ifndef RELEASE
        if( grid.VCRank() == 0 )
            cout << "Unaligned [VR,* ] <- [MR,* ]." << endl;
#endif
        const int r = grid.Height();
        const int c = grid.Width();
        const int p = grid.Size();
        const int row = grid.MCRank();
        const int col = grid.MRRank();
        const int colShiftOfA = A.ColShift();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        // We will SendRecv A[VR,*] within our process row to fix alignments.
        const int sendCol = (col+c+(colAlignment%c)-colAlignmentOfA) % c;
        const int recvCol = (col+c+colAlignmentOfA-(colAlignment%c)) % c;
        const int sendRank = sendCol + c*row;

        const int sendColShift = Shift( sendRank, colAlignment, p );
        const int sendColOffset = (sendColShift-colShiftOfA) / c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfSend = LocalLength(height,sendColShift,p);
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = maxLocalHeight * width;

        this->_auxMemory.Require( 2*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[portionSize];

        // Pack
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfSend; ++i )
                sendBuffer[i+j*localHeightOfSend] =
                    A.LocalEntry(sendColOffset+i*r,j);

        // Communicate
        SendRecv
        ( sendBuffer, portionSize, sendCol, 0,
          recvBuffer, portionSize, recvCol, MPI_ANY_TAG, grid.MRComm() );

        // Unpack
        for( int j=0; j<width; ++j )
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
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,MC]");
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
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    
    const Grid& grid = this->GetGrid();
    const int r = grid.Height();
    const int c = grid.Width();
    const int p = grid.Size();
    const int rankCM = grid.VCRank();
    const int rankRM = grid.VRRank();

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localHeightOfA = A.LocalHeight();
    const int maxLocalHeight = MaxLocalLength(height,p);

    const int portionSize = maxLocalHeight * width;

    const int colShift = this->ColShift();
    const int colShiftOfA = A.ColShift();

    // Compute which rowmajor rank has the colShift equal to our colShiftOfA
    const int sendRankRM = (rankRM+(p+colShiftOfA-colShift)) % p;

    // Compute which rowmajor rank has the A colShift that we need
    const int recvRankCM = (rankCM+(p+colShift-colShiftOfA)) % p;
    const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

    this->_auxMemory.Require( 2*portionSize );

    T* buffer = this->_auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[portionSize];

    // Pack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeightOfA; ++i )
            sendBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);

    // Communicate
    SendRecv
    ( sendBuffer, portionSize, sendRankRM, 0,
      recvBuffer, portionSize, recvRankRM, MPI_ANY_TAG, grid.VRComm() );

    // Unpack
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,VC]");
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
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [VR,* ]");
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
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() )
    {
        this->_localMatrix = A.LockedLocalMatrix();
    }
    else
    {
        const Grid& grid = this->GetGrid();
#ifndef RELEASE
        if( grid.VCRank() == 0 )
            cout << "Unaligned [VR,* ] <- [VR,* ]." << endl;
#endif
        const int rank = grid.VRRank();
        const int p = grid.Size();

        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localHeightOfA = A.LocalHeight();

        const int sendSize = localHeightOfA * width;
        const int recvSize = localHeight * width;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, grid.VRComm() );

        // Unpack
        for( int j=0; j<width; ++j )
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
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& grid = this->GetGrid();
    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(grid) );
    *A_MC_MR = A;

    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(grid) );
    *A_VC_Star = *A_MC_MR;
    delete A_MC_MR.release(); // lowers memory highwater

    *this = *A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,VR,Star>&
elemental::DistMatrixBase<T,VR,Star>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[VR,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int p = this->GetGrid().Size();
    const int colShift = this->ColShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = A.LocalEntry(colShift+i*p,j);
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
elemental::DistMatrix<R,VR,Star>::SetToRandomHPD()
{   
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw "Positive-definite matrices must be square.";
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->GetGrid().Size();
    const int colShift = this->ColShift();

    this->SetToRandom();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*p;
        if( i < width )
            this->LocalEntry(iLoc,i) += (R)this->Width();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,VR,Star>::SetToRandomHPD()
{   
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw "Positive-definite matrices must be square.";
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int p = this->GetGrid().Size();
    const int colShift = this->ColShift();

    this->SetToRandom();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*p;
        if( i < width )
            this->LocalEntry(iLoc,i) = 
                real(this->LocalEntry(iLoc,i)) + (R)this->Width();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,VR,Star>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire grid
    const Grid& grid = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % grid.Size();

    R u;
    if( grid.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / grid.Size();
        u = real(this->LocalEntry(iLoc,j));
    }
    Broadcast( &u, 1, ownerRank, grid.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,VR,Star>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire grid
    const Grid& grid = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % grid.Size();

    R u;
    if( grid.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / grid.Size();
        u = imag(this->LocalEntry(iLoc,j));
    }
    Broadcast( &u, 1, ownerRank, grid.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,VR,Star>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& grid = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % grid.Size();

    if( grid.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / grid.Size();
        real(this->LocalEntry(iLoc,j)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,VR,Star>::SetImag
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& grid = this->GetGrid();
    const int ownerRank = (i + this->ColAlignment()) % grid.Size();

    if( grid.VRRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / grid.Size();
        imag(this->LocalEntry(iLoc,j)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrixBase<int,   VR,Star>;
template class elemental::DistMatrixBase<float, VR,Star>;
template class elemental::DistMatrixBase<double,VR,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,VR,Star>;
template class elemental::DistMatrixBase<dcomplex,VR,Star>;
#endif

template class elemental::DistMatrix<int,     VR,Star>;
template class elemental::DistMatrix<float,   VR,Star>;
template class elemental::DistMatrix<double,  VR,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,VR,Star>;
template class elemental::DistMatrix<dcomplex,VR,Star>;
#endif

