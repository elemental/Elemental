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
elemental::DistMatrixBase<T,Star,MD>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::Print");
#endif
    if( this->GetGrid().VCRank() == 0 && s != "" )
        cout << s << endl;
        
    const int height     = this->Height();
    const int width      = this->Width();
    const int localWidth = this->LocalWidth();
    const int lcm        = this->GetGrid().LCM();
    const int inDiagonal = this->InDiagonal();

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
        for( int i=0; i<height; ++i )
            for( int j=0; j<localWidth; ++j )
                sendBuf[colShift+i+j*lcm*height] = this->_localMatrix(i,j);
    }

    // If we are the root, fill the receive buffer
    T* recvBuf = 0;
    if( this->GetGrid().VCRank() == 0 )
    {
        recvBuf = new T[height*width];     
        for( int i=0; i<height*width; ++i )
            recvBuf[i] = (T)0;
    }

    // Sum the contributions and send to the root
    Reduce
    ( sendBuf, recvBuf, height*width, MPI_SUM, 0, this->GetGrid().VCComm() );
    delete[] sendBuf;

    if( this->GetGrid().VCRank() == 0 )
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
elemental::DistMatrixBase<T,Star,MD>::AlignWith
( const DistMatrixBase<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWith([* ,MD])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.RowAlignment();
    this->_inDiagonal   = A.InDiagonal();
    if( this->InDiagonal() )
        this->_rowShift = A.RowShift();
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MD>::AlignWith
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWith([MD,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.ColAlignment();
    this->_inDiagonal   = A.InDiagonal();
    if( this->InDiagonal() )
        this->_rowShift = A.ColShift();
    this->_constrainedRowAlignment = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MD>::AlignRowsWith
( const DistMatrixBase<T,Star,MD>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,MD>::AlignRowsWith
( const DistMatrixBase<T,MD,Star>& A )
{ AlignWith( A ); }

template<typename T>
void
elemental::DistMatrixBase<T,Star,MD>::AlignWithDiag
( const DistMatrixBase<T,MC,MR>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiag([MC,MR])");
    this->AssertFreeRowAlignment();
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
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_rowShift = 
            ( g.DiagPathRank() + lcm - 
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MD>::AlignWithDiag
( const DistMatrixBase<T,MR,MC>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiag([MR,MC])");
    this->AssertFreeRowAlignment();
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
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_rowShift = 
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MD>::View
( DistMatrixBase<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width  = A.Width();
    this->_rowAlignment = A.RowAlignment();
    this->_inDiagonal   = A.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = A.RowShift();
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
elemental::DistMatrixBase<T,Star,MD>::LockedView
( const DistMatrixBase<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width  = A.Width();
    this->_rowAlignment = A.RowAlignment();
    this->_inDiagonal   = A.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = A.RowShift();
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
elemental::DistMatrixBase<T,Star,MD>::View
( DistMatrixBase<T,Star,MD>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View");
    this->AssertFreeRowAlignment();
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
        const int alignmentRank = A.RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        this->_rowAlignment = newAlignmentRank;
        this->_inDiagonal = A.InDiagonal();

        if( this->InDiagonal() )
        {
            this->_rowShift = 
                Shift( diagPathRank,
                       g.DiagPathRank(this->RowAlignment()),
                       lcm );
            int localWidthBefore = LocalLength( j, A.RowShift(), lcm );
            int localWidth = LocalLength( width, this->RowShift(), lcm );
        
            this->_localMatrix.View
            ( A.LocalMatrix(), i, localWidthBefore, height, localWidth );
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
elemental::DistMatrixBase<T,Star,MD>::LockedView
( const DistMatrixBase<T,Star,MD>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView");
    this->AssertFreeRowAlignment();
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
        const int alignmentRank = A.RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int newAlignmentRow = (alignmentRow + i) % r;
        const int newAlignmentCol = (alignmentCol + i) % c;
        const int newAlignmentRank = newAlignmentRow + r*newAlignmentCol;

        this->_rowAlignment = newAlignmentRank;
        this->_inDiagonal = A.InDiagonal();

        if( this->InDiagonal() )
        {
            this->_rowShift = 
                Shift( diagPathRank,
                       g.DiagPathRank( this->RowAlignment() ),
                       lcm );
            int localWidthBefore = LocalLength( j, A.RowShift(), lcm);
            int localWidth = LocalLength( width, this->RowShift(), lcm );
        
            this->_localMatrix.LockedView
            ( A.LockedLocalMatrix(), i, localWidthBefore, height, localWidth );
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
elemental::DistMatrixBase<T,Star,MD>::View1x2
( DistMatrixBase<T,Star,MD>& AL,
  DistMatrixBase<T,Star,MD>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View1x2");    
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width  = AL.Width() + AR.Width();
    this->_rowAlignment = AL.RowAlignment();
    this->_inDiagonal = AL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = AL.RowShift();
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
elemental::DistMatrixBase<T,Star,MD>::LockedView1x2
( const DistMatrixBase<T,Star,MD>& AL,
  const DistMatrixBase<T,Star,MD>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView1x2");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_rowAlignment = AL.RowAlignment();
    this->_inDiagonal = AL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = AL.RowShift();
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
elemental::DistMatrixBase<T,Star,MD>::View2x1
( DistMatrixBase<T,Star,MD>& AT,
  DistMatrixBase<T,Star,MD>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View2x1");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_rowAlignment = AT.RowAlignment();
    this->_inDiagonal = AT.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = AT.RowShift();
        this->_localMatrix.View2x1
        ( AT.LocalMatrix(), AB.LocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MD>::LockedView2x1
( const DistMatrixBase<T,Star,MD>& AT,
  const DistMatrixBase<T,Star,MD>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView2x1");
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_rowAlignment = AT.RowAlignment();
    this->_inDiagonal = AT.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = AT.RowShift();
        this->_localMatrix.LockedView2x1
        ( AT.LockedLocalMatrix(), AB.LockedLocalMatrix() );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MD>::View2x2
( DistMatrixBase<T,Star,MD>& ATL,
  DistMatrixBase<T,Star,MD>& ATR,
  DistMatrixBase<T,Star,MD>& ABL,
  DistMatrixBase<T,Star,MD>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::View2x2");
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
    this->_inDiagonal = ATL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = ATL.RowShift();
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
elemental::DistMatrixBase<T,Star,MD>::LockedView2x2
( const DistMatrixBase<T,Star,MD>& ATL,
  const DistMatrixBase<T,Star,MD>& ATR,
  const DistMatrixBase<T,Star,MD>& ABL,
  const DistMatrixBase<T,Star,MD>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::LockedView2x2");
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
    this->_inDiagonal = ATL.InDiagonal();
    if( this->InDiagonal() )
    {
        this->_rowShift = ATL.RowShift();
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
elemental::DistMatrixBase<T,Star,MD>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::ResizeTo");
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
        ( height, LocalLength(width,this->RowShift(),lcm) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,Star,MD>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    T u;
    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        u = this->_localMatrix(i,jLoc);
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MD>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::Set");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        this->_localMatrix(i,jLoc) = u;
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
elemental::DistMatrixBase<T,Star,MD>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int lcm = this->GetGrid().LCM();
        const int rowShift = this->RowShift();

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
                    this->_localMatrix(i,jLoc) = (T)0;
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
                    this->_localMatrix(i,jLoc) = (T)0;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star, MD>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int lcm = this->GetGrid().LCM();
        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        const int rowShift = this->RowShift();

        this->_localMatrix.SetToZero();
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*lcm;
            if( j < height )
                this->_localMatrix(j,jLoc) = (T)1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,MD>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToRandom");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->_localMatrix(i,j) = Random<T>();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [MC,MR] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [MC,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [* ,MR] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [MD,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_inDiagonal = A.InDiagonal();
            if( this->InDiagonal() )
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
        if( this->GetGrid().VCRank() == 0 )
            cout << "Unaligned [* ,MD] <- [* ,MD]." << endl;
#endif
        throw logic_error( "Unaligned [* ,MD] = [* ,MD] not yet implemented." );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [MR,MC] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [MR,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [* ,MC] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [VC,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [* ,VC] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [VR,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[* ,MD] = [* ,VR] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,MD>&
elemental::DistMatrixBase<T,Star,MD>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    if( this->InDiagonal() )
    {
        const int lcm = this->GetGrid().LCM();
        const int rowShift = this->RowShift();

        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<height; ++i )
                this->_localMatrix(i,j) = A.LocalEntry(i,rowShift+j*lcm);
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
void
elemental::DistMatrix<R,Star,MD>::AlignWithDiag
( const DistMatrixBase<R,MC,MR>& A, int offset )
{ DMB::AlignWithDiag( A, offset ); }

template<typename R>
void
elemental::DistMatrix<R,Star,MD>::AlignWithDiag
( const DistMatrixBase<R,MR,MC>& A, int offset )
{ DMB::AlignWithDiag( A, offset ); }

template<typename R>
void
elemental::DistMatrix<R,Star,MD>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    this->SetToRandom();

    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        const int lcm = this->GetGrid().LCM();
        const int rowShift = this->RowShift();

        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*lcm;
            if( j < height )
                this->LocalEntry(j,jLoc) += (R)this->Width();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,Star,MD>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    this->SetToRandom();

    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        const int lcm = this->GetGrid().LCM();
        const int rowShift = this->RowShift();

        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*lcm;
            if( j < height )
            {
                this->LocalEntry(j,jLoc) = 
                    real(this->LocalEntry(j,jLoc)) + (R)this->Width();
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,MD>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        u = real(this->_localMatrix(i,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,MD>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        u = imag(this->_localMatrix(i,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,MD>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        real(this->_localMatrix(i,jLoc)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,MD>::SetImag
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->GetGrid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        imag(this->_localMatrix(i,jLoc)) = u;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<R,Star,MD>::AlignWithDiag
( const DistMatrixBase<complex<R>,MC,MR>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiag([MC,MR])");
    this->AssertFreeRowAlignment();
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
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_rowShift = 
            ( g.DiagPathRank() + lcm - 
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<R,Star,MD>::AlignWithDiag
( const DistMatrixBase<complex<R>,MR,MC>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiag([MR,MC])");
    this->AssertFreeRowAlignment();
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
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_rowShift = 
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrixBase<int,   Star,MD>;
template class elemental::DistMatrixBase<float, Star,MD>;
template class elemental::DistMatrixBase<double,Star,MD>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,Star,MD>;
template class elemental::DistMatrixBase<dcomplex,Star,MD>;
#endif

template class elemental::DistMatrix<int,   Star,MD>;
template class elemental::DistMatrix<float, Star,MD>;
template class elemental::DistMatrix<double,Star,MD>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,Star,MD>;
template class elemental::DistMatrix<dcomplex,Star,MD>;
#endif

