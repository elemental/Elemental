/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/dist_matrix.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::utilities;
using namespace elemental::imports::mpi;

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Print");
#endif
    const Grid& g = this->Grid();
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

    if( g.InGrid() )
    {
        vector<T> sendBuf(height*width,0);
        if( inDiagonal )
        {
            const int colShift = this->ColShift();
#ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int i=0; i<localHeight; ++i )
                for( int j=0; j<width; ++j )
                    sendBuf[colShift+i*lcm+j*height] = this->GetLocalEntry(i,j);
        }

        // If we are the root, allocate a receive buffer
        vector<T> recvBuf;
        if( g.VCRank() == 0 )
            recvBuf.resize( height*width );

        // Sum the contributions and send to the root
        Reduce
        ( &sendBuf[0], &recvBuf[0], height*width, MPI_SUM, 0, g.VCComm() );

        if( g.VCRank() == 0 )
        {
            // Print the data
            for( int i=0; i<height; ++i )
            {
                for( int j=0; j<width; ++j )
                    cout << recvBuf[i+j*height] << " ";
                cout << "\n";
            }
            cout << endl;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::Align
( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MD,Star]::Align");
    this->AssertFreeColAlignment();
#endif
    this->AlignCols( colAlignment );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::AlignCols
( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MD,Star]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Size() )
        throw std::runtime_error( "Invalid column alignment for [MD,Star]" );
#endif
    this->_colAlignment = colAlignment;
    this->_inDiagonal = ( g.DiagPath() == g.DiagPath(colAlignment) );
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    if( g.InGrid() )
    {
        if( this->_inDiagonal )
            this->_colShift = Shift( g.DiagPathRank(), colAlignment, g.Size() );
        this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    if( this->Grid().InGrid() )
    {
        if( this->InDiagonal() )
            this->_colShift = A.ColShift();
        this->_localMatrix.ResizeTo( 0, 0 );
    }
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
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    if( this->Grid().InGrid() )
    {
        if( this->InDiagonal() )
            this->_colShift = A.RowShift();
        this->_localMatrix.ResizeTo( 0, 0 );
    }
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
    const Grid& g = this->Grid();
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
    const Grid& g = this->Grid();
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
        if( g.InGrid() )
        {
            this->_inDiagonal = 
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->_inDiagonal = false;
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        this->_colAlignment = ownerRow + r*ownerCol;
        if( g.InGrid() )
        {
            this->_inDiagonal = 
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->_inDiagonal = false;
    }
    if( this->InDiagonal() )
    {
        this->_colShift = 
            ( g.DiagPathRank() + lcm - 
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    const Grid& g = this->Grid();
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
    const Grid& g = this->Grid();
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
        if( g.InGrid() )
        {
            this->_inDiagonal = 
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->_inDiagonal = false;
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        this->_colAlignment = ownerRow + r*ownerCol;
        if( g.InGrid() )
        {
            this->_inDiagonal = 
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->_inDiagonal = false;
    }
    if( this->InDiagonal() )
    {
        this->_colShift = 
            ( g.DiagPathRank() + lcm - 
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
        const Grid& g = this->Grid();
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
        const Grid& g = this->Grid();
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
        const int lcm = this->Grid().LCM();
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
    const Grid& g = this->Grid();
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
        u = this->GetLocalEntry(iLoc,j);
    }
    Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );

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
    const Grid& g = this->Grid();
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
        this->SetLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::Update
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::Update");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->Grid();
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
        this->UpdateLocalEntry(iLoc,j,u);
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
        const int lcm = this->Grid().LCM();
        const int colShift = this->ColShift();

        if( shape == Lower )
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int lastZeroRow = ( side==Left ? j-offset-1
                                               : j-offset+height-width-1 );
                if( lastZeroRow >= 0 )
                {
                    int boundary = min( lastZeroRow+1, height );
                    int numZeroRows = LocalLength( boundary, colShift, lcm );
                    T* thisCol = this->LocalBuffer(0,j);
                    memset( thisCol, 0, numZeroRows*sizeof(T) );
                }
            }
        }
        else
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int firstZeroRow = 
                    ( side==Left ? max(j-offset+1,0)
                                 : max(j-offset+height-width+1,0) );
                int numNonzeroRows = LocalLength(firstZeroRow,colShift,lcm);
                if( numNonzeroRows < localHeight )
                {
                    T* thisCol = this->LocalBuffer(numNonzeroRows,j);
                    memset
                    ( thisCol, 0, (localHeight-numNonzeroRows)*sizeof(T) );
                }
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MD,Star>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int lcm = this->Grid().LCM();
        const int colShift = this->ColShift();

        if( shape == Upper )
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int lastRow = ( side==Left ? j-offset : j-offset+height-width );
                int boundary = min( lastRow+1, height );
                int numRows = LocalLength( boundary, colShift, lcm );
#ifdef RELEASE
                T* thisCol = this->LocalBuffer(0,j);
                for( int iLoc=0; iLoc<numRows; ++iLoc )
                    thisCol[iLoc] *= alpha;
#else
                for( int iLoc=0; iLoc<numRows; ++iLoc )
                {
                    const T value = this->GetLocalEntry(iLoc,j);
                    this->SetLocalEntry(iLoc,j,alpha*value);
                }
#endif
            }
        }
        else
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int firstRow = ( side==Left ? max(j-offset,0)
                                            : max(j+height-width-offset,0) );
                int numZeroRows = LocalLength( firstRow, colShift, lcm );
#ifdef RELEASE
                T* thisCol = this->LocalBuffer(numZeroRows,j);
                for( int iLoc=0; iLoc<(localHeight-numZeroRows); ++iLoc )
                    thisCol[iLoc] *= alpha;
#else
                for( int iLoc=numZeroRows; iLoc<localHeight; ++iLoc )
                {
                    const T value = this->GetLocalEntry(iLoc,j);
                    this->SetLocalEntry(iLoc,j,alpha*value);
                }
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
        const int lcm = this->Grid().LCM();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int colShift = this->ColShift();

        this->_localMatrix.SetToZero();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*lcm;
            if( i < width )
                this->SetLocalEntry(iLoc,i,1);
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
                this->SetLocalEntry(i,j,Random<T>());
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
#ifdef UNALIGNED_WARNINGS
        if( this->Grid().VCRank() == 0 )
            cerr << "Unaligned [MD,* ] <- [MD,* ]." << endl;
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
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,A.GetLocalEntry(colShift+i*lcm,j));
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class elemental::DistMatrixBase<int,   MD,Star>;
template class elemental::DistMatrixBase<float, MD,Star>;
template class elemental::DistMatrixBase<double,MD,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,MD,Star>;
template class elemental::DistMatrixBase<dcomplex,MD,Star>;
#endif

