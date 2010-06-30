/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
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
elemental::DistMatrixBase<T,MC,Star>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Print");
#endif
    const Grid& g = this->GetGrid();
    if( g.VCRank() == 0 && s != "" )
        cout << s << endl;
        
    const int height      = this->Height();
    const int width       = this->Width();
    const int localHeight = this->LocalHeight();
    const int r           = g.Height();
    const int colShift    = this->ColShift();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Only one process column needs to participate
    if( g.MRRank() == 0 )
    {
        T* sendBuf = new T[height*width];
        for( int i=0; i<height*width; ++i )
            sendBuf[i] = (T)0;
        for( int i=0; i<localHeight; ++i )
            for( int j=0; j<width; ++j )
                sendBuf[colShift+i*r+j*height] = this->LocalEntry(i,j);

        // If we are the root, fill the receive buffer
        T* recvBuf = 0;
        if( g.MCRank() == 0 )
        {
            recvBuf = new T[height*width];     
            for( int i=0; i<height*width; ++i )
                recvBuf[i] = (T)0;
        }

        // Sum the contributions and send to the root
        Reduce
        ( sendBuf, recvBuf, height*width, MPI_SUM, 0, g.MCComm() );
        delete[] sendBuf;

        if( g.MCRank() == 0 )
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
    Barrier( g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::AlignWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([MC,MR])");
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
elemental::DistMatrixBase<T,MC,Star>::AlignWith
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([MC,* ])");
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
elemental::DistMatrixBase<T,MC,Star>::AlignWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([MR,MC])");
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
elemental::DistMatrixBase<T,MC,Star>::AlignWith
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([* ,MC])");
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
elemental::DistMatrixBase<T,MC,Star>::AlignWith
( const DistMatrixBase<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([VC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.ColAlignment() % g.Height();
    this->_colShift = 
        Shift( g.MCRank(), this->ColAlignment(), g.Height() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::AlignWith
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([* ,VC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.RowAlignment() % g.Height();
    this->_colShift = 
        Shift( g.MCRank(), this->ColAlignment(), g.Height() );
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::AlignColsWith
( const DistMatrixBase<T,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::AlignColsWith
( const DistMatrixBase<T,MC,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::AlignColsWith
( const DistMatrixBase<T,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::AlignColsWith
( const DistMatrixBase<T,Star,MC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::AlignColsWith
( const DistMatrixBase<T,VC,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::AlignColsWith
( const DistMatrixBase<T,Star,VC>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::View
( DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::View");
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
elemental::DistMatrixBase<T,MC,Star>::LockedView
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedView");
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
elemental::DistMatrixBase<T,MC,Star>::View
( DistMatrixBase<T,MC,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const int r   = this->GetGrid().Height();
        const int row = this->GetGrid().MCRank();

        this->_colAlignment = (A.ColAlignment()+i) % r;
        this->_colShift = Shift( row, this->ColAlignment(), r );

        const int localHeightBefore = LocalLength( i, A.ColShift(), r );
        const int localHeight = LocalLength( height, this->ColShift(), r );

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
elemental::DistMatrixBase<T,MC,Star>::LockedView
( const DistMatrixBase<T,MC,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const int r   = this->GetGrid().Height();
        const int row = this->GetGrid().MCRank();

        this->_colAlignment = (A.ColAlignment()+i) % r;
        this->_colShift = Shift( row, this->ColAlignment(), r );

        const int localHeightBefore = LocalLength( i, A.ColShift(), r );
        const int localHeight = LocalLength( height, this->ColShift(), r );

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
elemental::DistMatrixBase<T,MC,Star>::View1x2
( DistMatrixBase<T,MC,Star>& AL,
  DistMatrixBase<T,MC,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::View1x2");    
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
elemental::DistMatrixBase<T,MC,Star>::LockedView1x2
( const DistMatrixBase<T,MC,Star>& AL,
  const DistMatrixBase<T,MC,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedView1x2");
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
elemental::DistMatrixBase<T,MC,Star>::View2x1
( DistMatrixBase<T,MC,Star>& AT,
  DistMatrixBase<T,MC,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::View2x1");
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
elemental::DistMatrixBase<T,MC,Star>::LockedView2x1
( const DistMatrixBase<T,MC,Star>& AT,
  const DistMatrixBase<T,MC,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedView2x1");
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
elemental::DistMatrixBase<T,MC,Star>::View2x2
( DistMatrixBase<T,MC,Star>& ATL,
  DistMatrixBase<T,MC,Star>& ATR,
  DistMatrixBase<T,MC,Star>& ABL,
  DistMatrixBase<T,MC,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::View2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
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
elemental::DistMatrixBase<T,MC,Star>::LockedView2x2
( const DistMatrixBase<T,MC,Star>& ATL,
  const DistMatrixBase<T,MC,Star>& ATR,
  const DistMatrixBase<T,MC,Star>& ABL,
  const DistMatrixBase<T,MC,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::LockedView2x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
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
elemental::DistMatrixBase<T,MC,Star>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    this->_height = height;
    this->_width = width;
    this->_localMatrix.ResizeTo
    ( LocalLength(height,this->ColShift(),this->GetGrid().Height()), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,MC,Star>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that 
    // row within each process column
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();

    T u;
    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        u = this->LocalEntry(iLoc,j);
    }
    Broadcast( &u, 1, ownerRow, g.MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
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
elemental::DistMatrixBase<T,MC,Star>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif

    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int r = this->GetGrid().Height();
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
                const int numZeros = LocalLength( boundary, colShift, r );
#ifdef RELEASE
                T* thisCol = &(this->LocalEntry(0,j));
                memset( thisCol, 0, numZeros*sizeof(T) );
#else
                for( int iLoc=0; iLoc<numZeros; ++iLoc )
                    this->LocalEntry(iLoc,j) = (T)0;
#endif
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
#ifdef RELEASE
            T* thisCol = &(this->LocalEntry(nonzeroLength,j));
            memset( thisCol, 0, (localHeight-nonzeroLength)*sizeof(T) );
#else
            for( int iLoc=nonzeroLength; iLoc<localHeight; ++iLoc )
                this->LocalEntry(iLoc,j) = (T)0;
#endif
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif

    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int r = this->GetGrid().Height();
    const int colShift = this->ColShift();

    this->_localMatrix.SetToZero();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;
        if( i < width )
            this->LocalEntry(iLoc,i) = (T)1;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,Star>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif

    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int bufSize = localHeight*width;

    this->_auxMemory.Require( bufSize );

    // Create random matrix on process column 0, then broadcast
    T* buffer = this->_auxMemory.Buffer();
    if( this->_g->MRRank() == 0 )
    {
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                buffer[i+j*localHeight] = Random<T>();
    }
    Broadcast( buffer, bufSize, 0, this->GetGrid().MRComm() );

    // Unpack
#ifdef RELEASE
    for( int j=0; j<width; ++j )
    {
        const T* bufferCol = &(buffer[j*localHeight]);
        T* thisCol = &(this->LocalEntry(0,j));
        memcpy( thisCol, bufferCol, localHeight*sizeof(T) );
    }
#else
    for( int j=0; j<width; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = buffer[i+j*localHeight];
#endif

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void 
elemental::DistMatrixBase<T,MC,Star>::SumOverRow()
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SumOverRow");
    this->AssertNotLockedView();
#endif

    const int localHeight = this->LocalHeight(); 
    const int localWidth = this->LocalWidth();
    const int localSize = max( localHeight*localWidth, MinCollectContrib );

    this->_auxMemory.Require( 2*localSize );
    T* buffer = this->_auxMemory.Buffer();
    T* sendBuf = &buffer[0];
    T* recvBuf = &buffer[localSize];

    // Pack
#ifdef RELEASE
    for( int j=0; j<localWidth; ++j )
    {
        const T* thisCol = &(this->LocalEntry(0,j));
        T* sendBufCol = &(sendBuf[j*localHeight]);
        memcpy( sendBufCol, thisCol, localHeight*sizeof(T) );
    }
#else
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            sendBuf[i+j*localHeight] = this->LocalEntry(i,j);
#endif

    // AllReduce sum
    AllReduce
    ( sendBuf, recvBuf, localSize, MPI_SUM, this->GetGrid().MRComm() );

    // Unpack
#ifdef RELEASE
    for( int j=0; j<localWidth; ++j )
    {
        const T* recvBufCol = &(recvBuf[j*localHeight]);
        T* thisCol = &(this->LocalEntry(0,j));
        memcpy( thisCol, recvBufCol, localHeight*sizeof(T) );
    }
#else
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = recvBuf[i+j*localHeight];
#endif

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MC,MR]");
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
                Shift( g.MCRank(), this->ColAlignment(), g.Height() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() )
    {
        if( A.Width() == 1 )
        {
            if( g.MRRank() == A.RowAlignment() )
                this->_localMatrix = A.LockedLocalMatrix();

            // Communicate
            Broadcast
            ( this->_localMatrix.Buffer(), this->LocalHeight(),
              A.RowAlignment(), g.MRComm() );
        }
        else
        {
            const int c = g.Width();

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidthOfA = A.LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,c);

            const int portionSize = 
                max(localHeight*maxLocalWidth,MinCollectContrib);

            this->_auxMemory.Require( (c+1)*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* originalData = &buffer[0];
            T* gatheredData = &buffer[portionSize];

            // Pack
#ifdef RELEASE
            for( int j=0; j<localWidthOfA; ++j )
            {
                const T* ACol = &(A.LocalEntry(0,j));
                T* originalDataCol = &(originalData[j*localHeight]);
                memcpy( originalDataCol, ACol, localHeight*sizeof(T) );
            }
#else
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<localHeight; ++i )
                    originalData[i+j*localHeight] = A.LocalEntry(i,j);
#endif

            // Communicate
            AllGather
            ( originalData, portionSize,
              gatheredData, portionSize, g.MRComm() );

            // Unpack
            const int rowAlignmentOfA = A.RowAlignment();
            for( int k=0; k<c; ++k )
            {
                const T* data = &gatheredData[k*portionSize];

                const int rowShift = Shift( k, rowAlignmentOfA, c );
                const int localWidth = LocalLength( width, rowShift, c );

#ifdef RELEASE
                for( int j=0; j<localWidth; ++j )
                {
                    const T* dataCol = &(data[j*localHeight]);
                    T* thisCol = &(this->LocalEntry(0,rowShift+j*c));
                    memcpy( thisCol, dataCol, localHeight*sizeof(T) );
                }
#else
                for( int j=0; j<localWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        this->LocalEntry(i,rowShift+j*c) = 
                            data[i+j*localHeight];
#endif
            }

            this->_auxMemory.Release();
        }
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [MC,* ] <- [MC,MR]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int row = g.MCRank();

        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
        const int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

        if( A.Width()==1 )
        {
            const int localHeight = this->LocalHeight();

            if( this->_g->MRRank() == A.RowAlignment() )
            {
                const int localHeightOfA = A.LocalHeight();

                this->_auxMemory.Require( localHeightOfA );
                T* buffer = this->_auxMemory.Buffer();

                // Pack
#ifdef RELEASE
                const T* ACol = &(A.LocalEntry(0,0));
                memcpy( buffer, ACol, localHeightOfA*sizeof(T) );
#else
                for( int i=0; i<localHeightOfA; ++i )
                    buffer[i] = A.LocalEntry(i,0);
#endif

                // Communicate
                SendRecv
                ( buffer, localHeightOfA, sendRow, 0,
                  this->_localMatrix.Buffer(), localHeight, recvRow, 
                  MPI_ANY_TAG, g.MCComm() );

                this->_auxMemory.Release();
            }

            // Communicate
            Broadcast
            ( this->_localMatrix.Buffer(), localHeight, A.RowAlignment(),
              g.MRComm() );
        }
        else
        {
            const int height = this->Height();
            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int localWidthOfA  = A.LocalWidth();
            const int maxLocalHeight = MaxLocalLength(height,r);
            const int maxLocalWidth  = MaxLocalLength(width,c);

            const int portionSize = 
                max(maxLocalHeight*maxLocalWidth,MinCollectContrib);

            this->_auxMemory.Require( (c+1)*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[portionSize];

            // Pack the currently owned local data of A into the second buffer
#ifdef RELEASE
            for( int j=0; j<localWidthOfA; ++j )
            {
                const T* ACol = &(A.LocalEntry(0,j));
                T* secondBufferCol = &(secondBuffer[j*localHeightOfA]);
                memcpy( secondBufferCol, ACol, localHeightOfA*sizeof(T) );
            }
#else
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    secondBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);
#endif

            // Perform the SendRecv: puts the new data into the first buffer
            SendRecv
            ( secondBuffer, portionSize, sendRow, 0,
              firstBuffer,  portionSize, recvRow, MPI_ANY_TAG, 
              g.MCComm() );

            // Use the output of the SendRecv as the input to the AllGather
            AllGather
            ( firstBuffer,  portionSize, 
              secondBuffer, portionSize, g.MRComm() );

            // Unpack the contents of each member of the process row
            const int rowAlignmentOfA = A.RowAlignment();
            for( int k=0; k<c; ++k )
            {
                const T* data = &secondBuffer[k*portionSize];

                const int rowShift = Shift( k, rowAlignmentOfA, c ); 
                const int localWidth = LocalLength( width, rowShift, c );
#ifdef RELEASE
                for( int j=0; j<localWidth; ++j )
                {
                    const T* dataCol = &(data[j*localHeight]);
                    T* thisCol = &(this->LocalEntry(0,rowShift+j*c));
                    memcpy( thisCol, dataCol, localHeight*sizeof(T) );
                }
#else
                for( int j=0; j<localWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        this->LocalEntry(i,rowShift+j*c) = 
                            data[i+j*localHeight];
#endif
            }

            this->_auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MC,* ]");
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
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [MC,* ] <- [MC,* ]." << endl;
#endif
        const int rank = g.MCRank();
        const int r = g.Height();

        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (rank+r+colAlignment-colAlignmentOfA) % r;
        const int recvRank = (rank+r+colAlignmentOfA-colAlignment) % r;

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
#ifdef RELEASE
        for( int j=0; j<width; ++j )
        {
            const T* ACol = &(A.LocalEntry(0,j));
            T* sendBufferCol = &(sendBuffer[j*localHeightOfA]);
            memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);
#endif

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.MCComm() );

        // Unpack
#ifdef RELEASE
        for( int j=0; j<width; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = &(this->LocalEntry(0,j));
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = recvBuffer[i+j*localHeight];
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,MC,MR> A_MC_MR(true,false,this->ColAlignment(),0,g);

    A_MC_MR = A;
    *this   = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MC,* ] = [MD,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MC,* ] = [* ,MD] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
    ( new DistMatrix<T,VR,Star>(g) );
    *A_VR_Star = A;

    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(true,this->ColAlignment(),g) );
    *A_VC_Star = *A_VR_Star;
    delete A_VR_Star.release(); // lowers memory highwater

    *this = *A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
    ( new DistMatrix<T,VR,Star>(g) );
    *A_VR_Star = A;

    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(true,this->ColAlignment(),g) );
    *A_VC_Star = *A_VR_Star;
    delete A_VR_Star.release(); // lowers memory highwater

    *this = *A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    auto_ptr< DistMatrix<T,MR,MC> > A_MR_MC( new DistMatrix<T,MR,MC>(g) );
    *A_MR_MC = A;

    auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
    ( new DistMatrix<T,VR,Star>(g) );
    *A_VR_Star = *A_MR_MC;
    delete A_MR_MC.release(); // lowers memory highwater

    auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
    ( new DistMatrix<T,VC,Star>(true,this->ColAlignment(),g) );
    *A_VC_Star = *A_VR_Star;
    delete A_VR_Star.release(); // lowers memory highwater

    *this = *A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
    if( A.GetGrid().VCRank() == 0 )
    {
        if( A.Width() == 1 )
        {
            cout << 
              "The vector version of [MC,* ] <- [VC,* ] is not yet written, but"
              " it only requires a modification of the vector version of "
              "[* ,MR] <- [* ,VR]" << endl;
        }
        else
        {
            cout << 
              "[MC,* ] <- [VC,* ] potentially causes a large amount of cache-"
              "thrashing. If possible avoid it by performing the redistribution"
              " with a (conjugate)-transpose: " << endl << 
              "  [* ,MC].(Conjugate)TransposeFrom([VC,* ])" << endl;
        }
    }
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.ColAlignment() % g.Height();
            this->_colShift = 
                Shift( g.MCRank(), this->ColAlignment(), g.Height() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->ColAlignment() == A.ColAlignment() % g.Height() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int row = g.MCRank();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(height,p);

        const int portionSize = max(maxLocalHeightOfA*width,MinCollectContrib);

        this->_auxMemory.Require( (c+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack 
#ifdef RELEASE
        for( int j=0; j<width; ++j )
        {
            const T* ACol = &(A.LocalEntry(0,j));
            T* originalDataCol = &(originalData[j*localHeightOfA]);
            memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                originalData[i+j*localHeightOfA] = A.LocalEntry(i,j);
#endif

        // Communicate 
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MRComm() );

        // Unpack
        const int colShift = this->ColShift();
        const int colAlignmentOfA = A.ColAlignment();
        for( int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];    

            const int colShiftOfA = Shift( row+r*k, colAlignmentOfA, p );
            const int colOffset = (colShiftOfA-colShift) / r;
            const int localHeight = LocalLength( height, colShiftOfA, p );

            for( int j=0; j<width; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(colOffset+i*c,j) = data[i+j*localHeight];
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifndef RELEASE
        if( g.VCRank() == 0 )
            cout << "Unaligned [MC,* ] <- [VC,* ]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int row = g.MCRank();
        const int rank = g.VCRank();

        // Perform the SendRecv to make A have the same colAlignment
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int colShift = this->ColShift();

        const int sendRank = (rank+p+colAlignment-colAlignmentOfA) % p;
        const int recvRank = (rank+p+colAlignmentOfA-colAlignment) % p;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeightOfA = MaxLocalLength(height,p);

        const int portionSize = max(maxLocalHeightOfA*width,MinCollectContrib);

        this->_auxMemory.Require( (c+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[portionSize];

        // Pack
#ifdef RELEASE
        for( int j=0; j<width; ++j )
        {
            const T* ACol = &(A.LocalEntry(0,j));
            T* secondBufferCol = &(secondBuffer[j*localHeightOfA]);
            memcpy( secondBufferCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                secondBuffer[i+j*localHeightOfA] = A.LocalEntry(i,j);
#endif

        // Perform the SendRecv: puts the new data into the first buffer
        SendRecv
        ( secondBuffer, portionSize, sendRank, 0,
          firstBuffer,  portionSize, recvRank, MPI_ANY_TAG, g.VCComm() );

        // Use the SendRecv as input to the AllGather
        AllGather
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MRComm() );

        // Unpack
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int colShiftOfA = Shift( row+r*k, colAlignment, p );
            const int colOffset = (colShiftOfA-colShift) / r;
            const int localHeight = LocalLength( height, colShiftOfA, p );

            for( int j=0; j<width; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->LocalEntry(colOffset+i*c,j) = data[i+j*localHeight];
        }

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
    ( new DistMatrix<T,Star,VR>(g) );
    *A_Star_VR = A;

    auto_ptr< DistMatrix<T,MC,MR> > A_MC_MR
    ( new DistMatrix<T,MC,MR>(true,false,this->ColAlignment(),0,g) );
    *A_MC_MR = *A_Star_VR;
    delete A_Star_VR.release(); // lowers memory highwater

    *this = *A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,VC,Star> A_VC_Star(true,this->ColAlignment(),g);

    A_VC_Star = A;
    *this = A_VC_Star;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,MC,MR> A_MC_MR(true,false,this->ColAlignment(),0,g);

    A_MC_MR = A;
    *this = A_MC_MR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,Star>&
elemental::DistMatrixBase<T,MC,Star>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int r = this->GetGrid().Height(); 
    const int colShift = this->ColShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->LocalEntry(i,j) = A.LocalEntry(colShift+i*r,j);
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
elemental::DistMatrix<R,MC,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int r = this->GetGrid().Height();
    const int colShift = this->ColShift();

    this->SetToRandom();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;
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
elemental::DistMatrix<complex<R>,MC,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int r = this->GetGrid().Height();
    const int colShift = this->ColShift();

    this->SetToRandom();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;
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
elemental::DistMatrix<complex<R>,MC,Star>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that 
    // row within each process column
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();

    R u;
    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        u = real(this->LocalEntry(iLoc,j));
    }
    Broadcast( &u, 1, ownerRow, g.MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,MC,Star>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that 
    // row within each process column
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();

    R u;
    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        u = imag(this->LocalEntry(iLoc,j));
    }
    Broadcast( &u, 1, ownerRow, g.MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}


template<typename R>
void
elemental::DistMatrix<complex<R>,MC,Star>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        real(this->LocalEntry(iLoc,j)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,Star>::SetImag
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        imag(this->LocalEntry(iLoc,j)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrixBase<int,   MC,Star>;
template class elemental::DistMatrixBase<float, MC,Star>;
template class elemental::DistMatrixBase<double,MC,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,MC,Star>;
template class elemental::DistMatrixBase<dcomplex,MC,Star>;
#endif

template class elemental::DistMatrix<int,   MC,Star>;
template class elemental::DistMatrix<float, MC,Star>;
template class elemental::DistMatrix<double,MC,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,MC,Star>;
template class elemental::DistMatrix<dcomplex,MC,Star>;
#endif

