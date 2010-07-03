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
elemental::DistMatrixBase<T,MD,Star>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Print");
#endif
    const Grid& g = this->GetGrid();
    if( g.VCRank() == 0 && s != "" )
        cout << s << endl;
        
    const int height      = this->Height();
    const int width       = this->Width();
    const int localHeight = this->LocalHeight();
    const int inDiagonal  = this->InDiagonal();
    const int lcm         = g.LCM();

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
        const int colShift = this->ColShift();
        for( int i=0; i<localHeight; ++i )
            for( int j=0; j<width; ++j )
                sendBuf[colShift+i*lcm+j*height] = this->LocalEntry(i,j);
    }

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

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::AlignWith
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWith([MD,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_inDiagonal = A.InDiagonal();
    if( this->InDiagonal() )
        this->_colShift = A.ColShift();
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::AlignWith
( const DistMatrixBase<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWith([* ,MD])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_inDiagonal = A.InDiagonal();
    if( this->InDiagonal() )
        this->_colShift = A.RowShift();
    this->_constrainedColAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::AlignColsWith
( const DistMatrixBase<T,MD,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::AlignColsWith
( const DistMatrixBase<T,Star,MD>& A )
{ AlignWith( A ); }

template<typename T>
bool
elemental::DistMatrixBase<T,MD,Star>::AlignedWithDiag
( const DistMatrixBase<T,MC,MR>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignedWithDiag([MC,MR])");
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const int ownerRow = colAlignment;
        const int ownerCol = (rowAlignment + offset) % c;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::AlignWithDiag
( const DistMatrixBase<T,MC,MR>& A, int offset )
{ 
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWithDiag([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int lcm = g.LCM();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = colAlignment;
        const int ownerCol = (rowAlignment + offset) % c;
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_colShift = 
            ( g.DiagPathRank() + lcm - 
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
bool
elemental::DistMatrixBase<T,MD,Star>::AlignedWithDiag
( const DistMatrixBase<T,MR,MC>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignedWithDiag([MR,MC])");
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const int ownerRow = rowAlignment;
        const int ownerCol = (colAlignment + offset) % c;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::AlignWithDiag
( const DistMatrixBase<T,MR,MC>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWithDiag([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int lcm = g.LCM();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = rowAlignment;
        const int ownerCol = (colAlignment + offset) % c;
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_colShift = 
            ( g.DiagPathRank() + lcm - 
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::View
( DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_colAlignment = A.ColAlignment();
    this->_inDiagonal = A.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_colShift = A.ColShift();
        this->_localMatrix.View( A.LocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::LockedView
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_colAlignment = A.ColAlignment();
    this->_inDiagonal = A.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_colShift = A.ColShift();
        this->_localMatrix.LockedView( A.LockedLocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::View
( DistMatrixBase<T,MD,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int diagPathRank = g.DiagPathRank();
        const int alignmentRank = A.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        this->_colAlignment = newAlignmentRank;
        this->_inDiagonal   = A.InDiagonal();

        if( this->_inDiagonal )
        {
            this->_colShift = 
                Shift( diagPathRank,
                       g.DiagPathRank( this->ColAlignment() ),
                       lcm );
            int localHeightBefore = LocalLength( i, A.ColShift(), lcm);
            int localHeight = LocalLength( height, this->ColShift(), lcm );

            this->_localMatrix.View
            ( A.LocalMatrix(), localHeightBefore, j, localHeight, width );
        }

    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::LockedView
( const DistMatrixBase<T,MD,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width  = width;
    {
        const Grid& g = this->GetGrid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int diagPathRank = g.DiagPathRank();
        const int alignmentRank = A.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        this->_colAlignment = newAlignmentRank;
        this->_inDiagonal = A.InDiagonal();

        if( this->InDiagonal() )
        {
            this->_colShift = 
                Shift( diagPathRank,
                       g.DiagPathRank( this->ColAlignment() ),
                       lcm );
            int localHeightBefore = LocalLength( i, A.ColShift(), lcm);
            int localHeight = LocalLength( height, this->ColShift(), lcm );
        
            this->_localMatrix.LockedView
            ( A.LockedLocalMatrix(), localHeightBefore, j, localHeight, width );
        }
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::View1x2
( DistMatrixBase<T,MD,Star>& AL,
  DistMatrixBase<T,MD,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View1x2");    
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_colAlignment = AL.ColAlignment();
    this->_inDiagonal = AL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_colShift = AL.ColShift();
        this->_localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::LockedView1x2
( const DistMatrixBase<T,MD,Star>& AL,
  const DistMatrixBase<T,MD,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView1x2");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_colAlignment = AL.ColAlignment();
    this->_inDiagonal = AL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_colShift = AL.ColShift();
        this->_localMatrix.LockedView1x2
        ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::View2x1
( DistMatrixBase<T,MD,Star>& AT,
  DistMatrixBase<T,MD,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_colAlignment = AT.ColAlignment();
    this->_inDiagonal = AT.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_colShift = AT.ColShift();
        this->_localMatrix.View2x1
        ( AT.LocalMatrix(),
          AB.LocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::LockedView2x1
( const DistMatrixBase<T,MD,Star>& AT,
  const DistMatrixBase<T,MD,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView2x1");
    this->AssertFreeColAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_colAlignment = AT.ColAlignment();
    this->_inDiagonal = AT.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_colShift = AT.ColShift();
        this->_localMatrix.LockedView2x1
        ( AT.LockedLocalMatrix(),
          AB.LockedLocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::View2x2
( DistMatrixBase<T,MD,Star>& ATL,
  DistMatrixBase<T,MD,Star>& ATR,
  DistMatrixBase<T,MD,Star>& ABL,
  DistMatrixBase<T,MD,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::View2x2");
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
    this->_inDiagonal = ATL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_colShift = ATL.ColShift();
        this->_localMatrix.View2x2
        ( ATL.LocalMatrix(), ATR.LocalMatrix(),
          ABL.LocalMatrix(), ABR.LocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::LockedView2x2
( const DistMatrixBase<T,MD,Star>& ATL,
  const DistMatrixBase<T,MD,Star>& ATR,
  const DistMatrixBase<T,MD,Star>& ABL,
  const DistMatrixBase<T,MD,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::LockedView2x2");
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
    this->_inDiagonal = ATL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_colShift = ATL.ColShift();
        this->_localMatrix.LockedView2x2
        ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
          ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    this->_height = height;
    this->_width = width;
    if( this->InDiagonal() )
    {
        const int lcm = this->GetGrid().LCM();
        this->_localMatrix.ResizeTo
        ( LocalLength(height,this->ColShift(),lcm), width );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,MD,Star>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    T u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        u = this->LocalEntry(iLoc,j);
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
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
elemental::DistMatrixBase<T,MD,Star>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int lcm = this->GetGrid().LCM();
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
                    const int numZeros = LocalLength( boundary, colShift, lcm );
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
                const int nonzeroLength = LocalLength(firstZero_i,colShift,lcm);
#ifdef RELEASE
                T* thisCol = &(this->LocalEntry(nonzeroLength,j));
                memset( thisCol, 0, (localHeight-nonzeroLength)*sizeof(T) );
#else
                for( int iLoc=nonzeroLength; iLoc<localHeight; ++iLoc )
                    this->LocalEntry(iLoc,j) = (T)0;
#endif
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int lcm = this->GetGrid().LCM();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int colShift = this->ColShift();

        this->_localMatrix.SetToZero();
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*lcm;
            if( i < width )
                this->LocalEntry(iLoc,i) = (T)1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = Random<T>();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MD,* ] = [MC,MR] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MD,* ] = [MC,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MD,* ] = [* ,MR] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [MD,* ]");
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
            this->_inDiagonal = A.InDiagonal();
            if( this->InDiagonal() )
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
        if( this->GetGrid().VCRank() == 0 )
            cout << "Unaligned [MD,* ] <- [MD,* ]." << endl;
#endif
        throw logic_error( "Unaligned [MD,* ] = [MD,* ] not yet implemented." );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MD,* ] = [* ,MD] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A ); 
#endif
    throw logic_error( "[MD,* ] = [MR,MC] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MD,* ] = [MR,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MD,* ] = [* ,MC] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MD,* ] = [VC,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MD,* ] = [* ,VC] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MD,* ] = [VR,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MD,* ] = [* ,VR] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MD,Star>&
elemental::DistMatrixBase<T,MD,Star>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    if( this->InDiagonal() )
    {
        const int lcm = this->_g->LCM();
        const int colShift = this->ColShift();

        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->LocalEntry(i,j) = A.LocalEntry(colShift+i*lcm,j);
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

//----------------------------------------------------------------------------//
// DistMatrix                                                                 //
//----------------------------------------------------------------------------//

template<typename R>
bool
elemental::DistMatrix<R,MD,Star>::AlignedWithDiag
( const DistMatrixBase<R,MC,MR>& A, int offset ) const
{ return DMB::AlignedWithDiag( A, offset ); }

template<typename R>
void
elemental::DistMatrix<R,MD,Star>::AlignWithDiag
( const DistMatrixBase<R,MC,MR>& A, int offset )
{ DMB::AlignWithDiag( A, offset ); }

template<typename R>
bool
elemental::DistMatrix<R,MD,Star>::AlignedWithDiag
( const DistMatrixBase<R,MR,MC>& A, int offset ) const
{ return DMB::AlignedWithDiag( A, offset ); }

template<typename R>
void
elemental::DistMatrix<R,MD,Star>::AlignWithDiag
( const DistMatrixBase<R,MR,MC>& A, int offset )
{ DMB::AlignWithDiag( A, offset ); }

template<typename R>
void
elemental::DistMatrix<R,MD,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    this->SetToRandom();

    if( this->InDiagonal() )
    {
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int lcm = this->GetGrid().LCM();
        const int colShift = this->ColShift();

        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*lcm;
            if( i < width )
                this->LocalEntry(iLoc,i) += (R)this->Width();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,MD,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    this->SetToRandom();

    if( this->InDiagonal() )
    {
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int lcm = this->GetGrid().LCM();
        const int colShift = this->ColShift();

        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*lcm;
            if( i < width )
            {
                this->LocalEntry(iLoc,i) = 
                    real(this->LocalEntry(iLoc,i)) + (R)this->Width();
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,MD,Star>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        u = real(this->LocalEntry(iLoc,j));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,MD,Star>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        u = imag(this->LocalEntry(iLoc,j));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MD,Star>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        real(this->LocalEntry(iLoc,j)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MD,Star>::SetImag
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        imag(this->LocalEntry(iLoc,j)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
bool
elemental::DistMatrix<R,MD,Star>::AlignedWithDiag
( const DistMatrixBase<complex<R>,MC,MR>& A, int offset ) const
{ 
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignedWithDiag([MC,MR])");
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const int ownerRow = colAlignment;
        const int ownerCol = (rowAlignment + offset) % c;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename R>
void
elemental::DistMatrix<R,MD,Star>::AlignWithDiag
( const DistMatrixBase<complex<R>,MC,MR>& A, int offset )
{ 
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWithDiag([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int lcm = g.LCM();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = colAlignment;
        const int ownerCol = (rowAlignment + offset) % c;
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_colShift = 
            ( g.DiagPathRank() + lcm - 
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
bool
elemental::DistMatrix<R,MD,Star>::AlignedWithDiag
( const DistMatrixBase<complex<R>,MR,MC>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignedWithDiag([MR,MC])");
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const int ownerRow = rowAlignment;
        const int ownerCol = (colAlignment + offset) % c;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}
template<typename R>
void
elemental::DistMatrix<R,MD,Star>::AlignWithDiag
( const DistMatrixBase<complex<R>,MR,MC>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWithDiag([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();
    const int lcm = g.LCM();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = rowAlignment;
        const int ownerCol = (colAlignment + offset) % c;
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_colShift = 
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrixBase<int,   MD,Star>;
template class elemental::DistMatrixBase<float, MD,Star>;
template class elemental::DistMatrixBase<double,MD,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,MD,Star>;
template class elemental::DistMatrixBase<dcomplex,MD,Star>;
#endif

template class elemental::DistMatrix<int,   MD,Star>;
template class elemental::DistMatrix<float, MD,Star>;
template class elemental::DistMatrix<double,MD,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,MD,Star>;
template class elemental::DistMatrix<dcomplex,MD,Star>;
#endif

