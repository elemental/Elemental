/*
   Copyright (c) 2009-2010, Jack Poulson
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
using namespace elemental::wrappers::mpi;

//----------------------------------------------------------------------------//
// DistMatrixBase                                                             //
//----------------------------------------------------------------------------//

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::Print");
#endif
    const Grid& g = this->GetGrid();
    const int r = g.Height();
    const int c = g.Width();

    if( g.VCRank() == 0 && s != "" )
        cout << s << endl;

    const int height = this->Height();
    const int width  = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth  = this->LocalWidth();
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
    vector<T> sendBuf(height*width,0);
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            sendBuf[colShift+i*r+(rowShift+j*c)*height] = 
                this->GetLocalEntry(i,j);

    // If we are the root, allocate a receive buffer
    vector<T> recvBuf;
    if( g.VCRank() == 0 )
        recvBuf.resize( height*width );

    // Sum the contributions and send to the root
    Reduce( &sendBuf[0], &recvBuf[0], height*width, MPI_SUM, 0, g.VCComm() );

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
    Barrier( g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::Align
( int colAlignment, int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::Align");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
#endif
    const Grid& g = this->GetGrid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Height() )
        throw std::runtime_error( "Invalid column alignment for [MC,MR]" );
    if( rowAlignment < 0 || rowAlignment >= g.Width() )
        throw std::runtime_error( "Invalid row alignment for [MC,MR]" );
#endif
    this->_colAlignment = colAlignment;
    this->_rowAlignment = rowAlignment;
    this->_colShift = Shift( g.MCRank(), colAlignment, g.Height() );
    this->_rowShift = Shift( g.MRRank(), rowAlignment, g.Width() );
    this->_constrainedColAlignment = true;
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignCols
( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const Grid& g = this->GetGrid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Height() )
        throw std::runtime_error( "Invalid column alignment for [MC,MR]" );
#endif
    this->_colAlignment = colAlignment;
    this->_colShift = Shift( g.MCRank(), colAlignment, g.Height() );
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
elemental::DistMatrixBase<T,MC,MR>::AlignRows
( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const Grid& g = this->GetGrid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Width() )
        throw std::runtime_error( "Invalid row alignment for [MC,MR]" );
#endif
    this->_rowAlignment = rowAlignment;
    this->_rowShift = Shift( g.MRRank(), rowAlignment, g.Width() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MC,MR])");
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignWith
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_colShift     = A.ColShift();
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
elemental::DistMatrixBase<T,MC,MR>::AlignWith
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift     = A.RowShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MR,MC])");
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignWith
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift     = A.ColShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignWith
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_colShift     = A.RowShift();
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
elemental::DistMatrixBase<T,MC,MR>::AlignWith
( const DistMatrixBase<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([VC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.ColAlignment() % g.Height();
    this->_colShift = 
        Shift( g.MCRank(), this->ColAlignment(), g.Height() );
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
elemental::DistMatrixBase<T,MC,MR>::AlignWith
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,VC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.RowAlignment() % g.Height();
    this->_colShift = 
        Shift( g.MCRank(), this->ColAlignment(), g.Height() );
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
elemental::DistMatrixBase<T,MC,MR>::AlignWith
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([VR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.ColAlignment() % g.Width();
    this->_rowShift = 
        Shift( g.MRRank(), this->RowAlignment(), g.Width() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignWith
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,VR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.RowAlignment() % g.Width();
    this->_rowShift = 
        Shift( g.MRRank(), this->RowAlignment(), g.Width() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignColsWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_colShift     = A.ColShift();
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
elemental::DistMatrixBase<T,MC,MR>::AlignColsWith
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([MC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_colShift     = A.ColShift();
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
elemental::DistMatrixBase<T,MC,MR>::AlignColsWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_colShift     = A.RowShift();
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
elemental::DistMatrixBase<T,MC,MR>::AlignColsWith
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([* ,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_colShift     = A.RowShift();
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
elemental::DistMatrixBase<T,MC,MR>::AlignColsWith
( const DistMatrixBase<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([VC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.ColAlignment() % g.Height();
    this->_colShift = 
        Shift( g.MCRank(), this->ColAlignment(), g.Height() );
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
elemental::DistMatrixBase<T,MC,MR>::AlignColsWith
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([* ,VC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_colAlignment = A.RowAlignment() % g.Height();
    this->_colShift = 
        Shift( g.MCRank(), this->ColAlignment(), g.Height() );
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
elemental::DistMatrixBase<T,MC,MR>::AlignRowsWith
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift     = A.RowShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignRowsWith
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([* ,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift     = A.RowShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignRowsWith
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift     = A.ColShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignRowsWith
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([MR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift     = A.ColShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignRowsWith
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([VR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.ColAlignment() % g.Width();
    this->_rowShift = 
        Shift( g.MRRank(), this->RowAlignment(), g.Width() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::AlignRowsWith
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([* ,VR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->GetGrid();
    this->_rowAlignment = A.RowAlignment() % g.Width();
    this->_rowShift = 
        Shift( g.MRRank(), this->RowAlignment(), g.Width() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::View
( DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View");
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
elemental::DistMatrixBase<T,MC,MR>::LockedView
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView");
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
elemental::DistMatrixBase<T,MC,MR>::View
( DistMatrixBase<T,MC,MR>& A, 
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width  = width;
    {
        const Grid& g = this->GetGrid();
        const int r   = g.Height();
        const int c   = g.Width();
        const int row = g.MCRank();
        const int col = g.MRRank();

        this->_colAlignment = (A.ColAlignment()+i) % r;
        this->_rowAlignment = (A.RowAlignment()+j) % c;
  
        this->_colShift = Shift( row, this->ColAlignment(), r );
        this->_rowShift = Shift( col, this->RowAlignment(), c );

        const int localHeightBehind = LocalLength(i,A.ColShift(),r);
        const int localWidthBehind  = LocalLength(j,A.RowShift(),c);

        const int localHeight = LocalLength( height, this->ColShift(), r );
        const int localWidth  = LocalLength( width,  this->RowShift(), c );

        this->_localMatrix.View
        ( A.LocalMatrix(), localHeightBehind, localWidthBehind,
                           localHeight,       localWidth );
    }
    this->_viewing = true;
    this->_lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::LockedView
( const DistMatrixBase<T,MC,MR>& A, 
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width  = width;
    {
        const Grid& g = this->GetGrid();
        const int r   = g.Height();
        const int c   = g.Width();
        const int row = g.MCRank();
        const int col = g.MRRank();

        this->_colAlignment = (A.ColAlignment()+i) % r;
        this->_rowAlignment = (A.RowAlignment()+j) % c;
  
        this->_colShift = Shift( row, this->ColAlignment(), r );
        this->_rowShift = Shift( col, this->RowAlignment(), c );

        const int localHeightBehind = LocalLength(i,A.ColShift(),r);
        const int localWidthBehind  = LocalLength(j,A.RowShift(),c);

        const int localHeight = LocalLength( height, this->ColShift(), r );
        const int localWidth  = LocalLength( width,  this->RowShift(), c );

        this->_localMatrix.LockedView
        ( A.LockedLocalMatrix(), localHeightBehind, localWidthBehind,
                                 localHeight,       localWidth );
    }
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::View1x2
( DistMatrixBase<T,MC,MR>& AL, 
  DistMatrixBase<T,MC,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View1x2");
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
elemental::DistMatrixBase<T,MC,MR>::LockedView1x2
( const DistMatrixBase<T,MC,MR>& AL, 
  const DistMatrixBase<T,MC,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView1x2");
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
    ( AL.LockedLocalMatrix(), 
      AR.LockedLocalMatrix() );
    this->_viewing = true;
    this->_lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::View2x1
( DistMatrixBase<T,MC,MR>& AT,
  DistMatrixBase<T,MC,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View2x1");
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
elemental::DistMatrixBase<T,MC,MR>::LockedView2x1
( const DistMatrixBase<T,MC,MR>& AT,
  const DistMatrixBase<T,MC,MR>& AB )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView2x1");
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
elemental::DistMatrixBase<T,MC,MR>::View2x2
( DistMatrixBase<T,MC,MR>& ATL, 
  DistMatrixBase<T,MC,MR>& ATR,
  DistMatrixBase<T,MC,MR>& ABL,
  DistMatrixBase<T,MC,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::View2x2");
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
elemental::DistMatrixBase<T,MC,MR>::LockedView2x2
( const DistMatrixBase<T,MC,MR>& ATL, 
  const DistMatrixBase<T,MC,MR>& ATR,
  const DistMatrixBase<T,MC,MR>& ABL,
  const DistMatrixBase<T,MC,MR>& ABR )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::LockedView2x2");
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
elemental::DistMatrixBase<T,MC,MR>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::ResizeTo");
    this->AssertNotLockedView();
#endif
    this->_height = height;
    this->_width  = width;
    this->_localMatrix.ResizeTo
    ( LocalLength(height,this->ColShift(),this->GetGrid().Height()),
      LocalLength(width, this->RowShift(),this->GetGrid().Width()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,MC,MR>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::Get");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = this->GetLocalEntry(iLoc,jLoc);
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::Set");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const int jLoc = (j-this->RowShift()) / g.Width();
        this->SetLocalEntry(iLoc,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::GetDiagonal
( DistMatrixBase<T,MD,Star>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetDiagonal");
    this->AssertNotLockedView();
#endif
    int length = this->DiagonalLength( offset );
#ifndef RELEASE
    if( d.Viewing() && length != d.Height() )
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:" << endl
            << "  A ~ " << this->Height() << " x " << this->Width() << endl
            << "  d ~ " << d.Height() << " x " << d.Width() << endl
            << "  A diag length: " << length << endl;
        throw logic_error( msg.str() );
    }
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiag( *this, offset ) )
        throw logic_error( "d must be aligned with the 'offset' diagonal." );
#endif
    if( !d.Viewing() )
    {
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiag( *this, offset );
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const T value = 
                this->GetLocalEntry(iLocStart+k*(lcm/r),jLocStart+k*(lcm/c));
            d.SetLocalEntry(k,0,value);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::GetDiagonal
( DistMatrixBase<T,Star,MD>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetDiagonal");
    this->AssertNotLockedView();
#endif
    int length = this->DiagonalLength( offset );
#ifndef RELEASE
    if( d.Viewing() && length != d.Width() )
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:" << endl
            << "  A ~ " << this->Height() << " x " << this->Width() << endl
            << "  d ~ " << d.Height() << " x " << d.Width() << endl
            << "  A diag length: " << length << endl;
        throw logic_error( msg.str() );
    }
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiag( *this, offset ) )
        throw logic_error( "d must be aligned with the 'offset' diagonal." );
#endif
    if( !d.Viewing() )
    {
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiag( *this, offset );
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const T value = 
                this->GetLocalEntry(iLocStart+k*(lcm/r),jLocStart+k*(lcm/c));
            d.SetLocalEntry(0,k,value);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::SetDiagonal
( const DistMatrixBase<T,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetDiagonal");
    if( d.Width() != 1 )
        throw logic_error( "d must be a column vector." );
    {
        int length = this->DiagonalLength( offset );
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
    if( !d.AlignedWithDiag( *this, offset ) )
        throw logic_error( "d must be aligned with the 'offset' diagonal." );
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const T value = d.GetLocalEntry(k,0);
            this->SetLocalEntry(iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),value);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::SetDiagonal
( const DistMatrixBase<T,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetDiagonal");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        int length = this->DiagonalLength( offset );
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
    if( !d.AlignedWithDiag( *this, offset ) )
        throw logic_error( "d must be aligned with the 'offset' diagonal." );
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const T value = d.GetLocalEntry(0,k);
            this->SetLocalEntry(iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),value);
        }
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
elemental::DistMatrixBase<T,MC,MR>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = this->GetGrid().Height();
    const int c = this->GetGrid().Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    if( shape == Lower )
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            int lastZeroRow = ( side==Left ? j-offset-1
                                           : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                int boundary = min( lastZeroRow+1, height );
                int numZeroRows = LocalLength( boundary, colShift, r );
#ifdef RELEASE
                T* thisCol = this->LocalBuffer(0,jLoc);
                memset( thisCol, 0, numZeroRows*sizeof(T) );
#else
                for( int iLoc=0; iLoc<numZeroRows; ++iLoc )
                    this->SetLocalEntry(iLoc,jLoc,0);
#endif
            }
        }
    }
    else
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            int firstZeroRow = ( side==Left ? max(j-offset+1,0)
                                            : max(j-offset+height-width+1,0) );
            int numNonzeroRows = LocalLength(firstZeroRow,colShift,r);
#ifdef RELEASE
            if( numNonzeroRows < localHeight )
            {
                T* thisCol = this->LocalBuffer(numNonzeroRows,jLoc);
                memset( thisCol, 0, (localHeight-numNonzeroRows)*sizeof(T) );
            }
#else
            for( int iLoc=numNonzeroRows; iLoc<localHeight; ++iLoc )
                this->SetLocalEntry(iLoc,jLoc,0);
#endif
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = this->GetGrid().Height();
    const int c = this->GetGrid().Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    if( shape == Upper )
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            int lastRow = ( side==Left ? j-offset : j-offset+height-width );
            int boundary = min( lastRow+1, height );
            int numRows = LocalLength( boundary, colShift, r );
#ifdef RELEASE
            T* thisCol = this->LocalBuffer(0,jLoc);
            for( int iLoc=0; iLoc<numRows; ++iLoc )
                thisCol[iLoc] *= alpha;
#else
            for( int iLoc=0; iLoc<numRows; ++iLoc )
            {
                const T value = this->GetLocalEntry(iLoc,jLoc);
                this->SetLocalEntry(iLoc,jLoc,alpha*value);
            }
#endif
        }
    }
    else
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            int firstRow = ( side==Left ? max(j-offset,0) 
                                        : max(j-offset+height-width,0) );
            int numZeroRows = LocalLength( firstRow, colShift, r );
#ifdef RELEASE
            T* thisCol = this->LocalBuffer(numZeroRows,jLoc);
            for( int iLoc=0; iLoc<(localHeight-numZeroRows); ++iLoc )
                thisCol[iLoc] *= alpha;
#else
            for( int iLoc=numZeroRows; iLoc<localHeight; ++iLoc )
            {
                const T value = this->GetLocalEntry(iLoc,jLoc);
                this->SetLocalEntry(iLoc,jLoc,alpha*value);
            }
#endif
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = this->GetGrid().Height();
    const int c = this->GetGrid().Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->_localMatrix.SetToZero();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;                
        if( i % c == rowShift )
        {
            const int jLoc = (i-rowShift) / c;
            if( jLoc < localWidth )
                this->SetLocalEntry(iLoc,jLoc,1);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToRandom");
    this->AssertNotLockedView();
#endif
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    for( int i=0; i<localHeight; ++i )
        for( int j=0; j<localWidth; ++j )
            this->SetLocalEntry(i,j,Random<T>());
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::ConjugateTransposeFrom
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::ConjugateTransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.RowAlignment();
            this->_colShift = 
                Shift( g.MCRank(), this->ColAlignment(), g.Height() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->ColAlignment() == A.RowAlignment() )
    {
        const int c = g.Width();
        const int rowShift = this->RowShift();

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidth; ++j )
        {
            for( int i=0; i<localHeight; ++i )
            {
                const T value = Conj( A.GetLocalEntry(rowShift+j*c,i) );
                this->SetLocalEntry(i,j,value);
            }
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MC,MR]::ConjugateTransposeFrom." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int rank = g.MCRank();
        const int rowShift = this->RowShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+r+colAlignment-rowAlignmentOfA) % r;
        const int recvRank = (rank+r+rowAlignmentOfA-colAlignment) % r;

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = localWidthOfA * localWidth;
        const int recvSize = localHeight * localWidth;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localWidthOfA; ++i )
                sendBuffer[i+j*localWidth] = 
                    Conj( A.GetLocalEntry(rowShift+j*c,i) );

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.MCComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*localHeight]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::TransposeFrom
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::TransposeFrom");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSizeAsTranspose( A );
#endif
    const Grid& g = this->GetGrid();
    if( !this->Viewing() )
    {
        if( !this->ConstrainedColAlignment() )
        {
            this->_colAlignment = A.RowAlignment();
            this->_colShift = 
                Shift( g.MCRank(), this->ColAlignment(), g.Height() );
        }
        this->ResizeTo( A.Width(), A.Height() );
    }

    if( this->ColAlignment() == A.RowAlignment() )
    {
        const int c = g.Width();
        const int rowShift = this->RowShift();

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,A.GetLocalEntry(rowShift+j*c,i));
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MC,MR]::TransposeFrom." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int rank = g.MCRank();
        const int rowShift = this->RowShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRank = (rank+r+colAlignment-rowAlignmentOfA) % r;
        const int recvRank = (rank+r+rowAlignmentOfA-colAlignment) % r;

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();

        const int sendSize = localWidthOfA * localWidth;
        const int recvSize = localHeight * localWidth;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localWidthOfA; ++i )
                sendBuffer[i+j*localWidth] = A.GetLocalEntry(rowShift+j*c,i);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.MCComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*localHeight]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR] = [MC,MR]");
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
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MC,MR] <- [MC,MR]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int row = g.MCRank();
        const int col = g.MRRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
        const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;
        const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;
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
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* sendBufferCol = &(sendBuffer[j*localHeightOfA]);
            memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i+j*localHeightOfA] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.VCComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*localHeight]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR] = [MC,* ]");
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
        const int c = g.Width();
        const int rowShift = this->RowShift();

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,rowShift+j*c);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, ACol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,A.GetLocalEntry(i,rowShift+j*c));
#endif
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MC,MR] <- [MC,* ]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int rank = g.MCRank();
        const int rowShift = this->RowShift();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRank = (rank+r+colAlignment-colAlignmentOfA) % r;
        const int recvRank = (rank+r+colAlignmentOfA-colAlignment) % r;

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
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,rowShift+j*c);
            T* sendBufferCol = &(sendBuffer[j*localHeightOfA]);
            memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i+j*localHeightOfA] = 
                    A.GetLocalEntry(i,rowShift+j*c);
#endif

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.MCComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*localHeight]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,MR]");
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
                Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int colShift = this->ColShift();

        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,A.GetLocalEntry(colShift+i*r,j));
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MC,MR] <- [* ,MR]." << endl;
#endif
        const int r = g.Height(); 
        const int c = g.Width();
        const int col = g.MRRank();
        const int colShift = this->ColShift();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

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
#ifdef _OPENMP
        #pragma omp parallel for  COLLAPSE(2)
#endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeight; ++i )
                sendBuffer[i+j*localHeight] = A.GetLocalEntry(colShift+i*r,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendCol, 0,
          recvBuffer, recvSize, recvCol, MPI_ANY_TAG, g.MRComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for  
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for  COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*localHeight]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MC,MR] = [MD,* ] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    throw logic_error( "[MC,MR] = [* ,MD] not yet implemented." );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    if( A.Width() == 1 )
    {
        if( !this->Viewing() )
            this->ResizeTo( A.Height(), 1 );

        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int myRow = g.MCRank();
        const int myCol = g.MRRank();
        const int rankCM = g.VCRank();
        const int rankRM = g.VRRank();
        const int ownerCol = this->RowAlignment();
        const int ownerRow = A.RowAlignment();
        const int colAlignment = this->ColAlignment();
        const int colAlignmentOfA = A.ColAlignment();
        const int colShift = this->ColShift();
        const int colShiftOfA = A.ColShift();

        const int height = A.Height();
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = max(maxLocalHeight,MinCollectContrib);

        const int colShiftVC = Shift(rankCM,colAlignment,p);
        const int colShiftVROfA = Shift(rankRM,colAlignmentOfA,p);
        const int sendRankCM = (rankCM+(p+colShiftVROfA-colShiftVC)) % p;
        const int recvRankRM = (rankRM+(p+colShiftVC-colShiftVROfA)) % p;
        const int recvRankCM = (recvRankRM/c)+r*(recvRankRM%c);

        this->_auxMemory.Require( (r+c)*portionSize );
        T* buffer = this->_auxMemory.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[c*portionSize];

        if( myRow == ownerRow )
        {
            // Pack
#ifdef _OPENMP
            #pragma omp parallel for  
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const int shift = Shift(myCol+c*k,colAlignmentOfA,p);
                const int offset = (shift-colShiftOfA) / c;
                const int thisLocalHeight = LocalLength(height,shift,p);

                for( int i=0; i<thisLocalHeight; ++i )
                    data[i] = A.GetLocalEntry(offset+i*r,0);
            }
        }

        // A[VR,* ] <- A[MR,MC]
        Scatter
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerRow, g.MCComm() );

        // A[VC,* ] <- A[VR,* ]
        SendRecv
        ( sendBuf, portionSize, sendRankCM, 0,
          recvBuf, portionSize, recvRankCM, MPI_ANY_TAG, g.VCComm() );

        // A[MC,MR] <- A[VC,* ]
        Gather
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerCol, g.MRComm() );

        if( myCol == ownerCol )
        {
            // Unpack
#ifdef _OPENMP
            #pragma omp parallel for  
#endif
            for( int k=0; k<c; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const int shift = Shift(myRow+r*k,colAlignment,p);
                const int offset = (shift-colShift) / r;
                const int thisLocalHeight = LocalLength(height,shift,p);

                for( int i=0; i<thisLocalHeight; ++i )
                    this->SetLocalEntry(offset+i*c,0,data[i]);
            }
        }

        this->_auxMemory.Release();
    }
    else if( A.Height() == 1 )
    {
        if( !this->Viewing() )
            this->ResizeTo( 1, A.Width() );

        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();
        const int myRow = g.MCRank();
        const int myCol = g.MRRank();
        const int rankCM = g.VCRank();
        const int rankRM = g.VRRank();
        const int ownerRow = this->ColAlignment();
        const int ownerCol = A.ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int rowShift = this->RowShift();
        const int rowShiftOfA = A.RowShift();

        const int width = A.Width();
        const int maxLocalWidth = MaxLocalLength(width,p);

        const int portionSize = max(maxLocalWidth,MinCollectContrib);

        const int rowShiftVR = Shift(rankRM,rowAlignment,p);
        const int rowShiftVCOfA = Shift(rankCM,rowAlignmentOfA,p);
        const int sendRankRM = (rankRM+(p+rowShiftVCOfA-rowShiftVR)) % p;
        const int recvRankCM = (rankCM+(p+rowShiftVR-rowShiftVCOfA)) % p;
        const int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

        this->_auxMemory.Require( (r+c)*portionSize );
        T* buffer = this->_auxMemory.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[r*portionSize];

        if( myCol == ownerCol )
        {
            // Pack
#ifdef _OPENMP
            #pragma omp parallel for  
#endif
            for( int k=0; k<c; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const int shift = Shift(myRow+r*k,rowAlignmentOfA,p);
                const int offset = (shift-rowShiftOfA) / r;
                const int thisLocalWidth = LocalLength(width,shift,p);

                for( int j=0; j<thisLocalWidth; ++j )
                    data[j] = A.GetLocalEntry(0,offset+j*c);
            }
        }

        // A[* ,VC] <- A[MR,MC]
        Scatter
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerCol, g.MRComm() );

        // A[* ,VR] <- A[* ,VC]
        SendRecv
        ( sendBuf, portionSize, sendRankRM, 0,
          recvBuf, portionSize, recvRankRM, MPI_ANY_TAG, g.VRComm() );

        // A[MC,MR] <- A[* ,VR]
        Gather
        ( recvBuf, portionSize, 
          sendBuf, portionSize, ownerRow, g.MCComm() );

        if( myRow == ownerRow )
        {
            // Unpack
#ifdef _OPENMP
            #pragma omp parallel for  
#endif
            for( int k=0; k<r; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const int shift = Shift(myCol+c*k,rowAlignment,p);
                const int offset = (shift-rowShift) / c;
                const int thisLocalWidth = LocalLength(width,shift,p);

                for( int j=0; j<thisLocalWidth; ++j )
                    this->SetLocalEntry(0,offset+j*r,data[j]);
            }
        }

        this->_auxMemory.Release();
    }
    else
    {
        if( A.Height() >= A.Width() )
        {
            auto_ptr< DistMatrix<T,VR,Star> > A_VR_Star
            ( new DistMatrix<T,VR,Star>(g) );

            *A_VR_Star = A;

            auto_ptr< DistMatrix<T,VC,Star> > A_VC_Star
            ( new DistMatrix<T,VC,Star>(true,this->ColAlignment(),g) );
            *A_VC_Star = *A_VR_Star;
            delete A_VR_Star.release(); // lowers memory highwater

            *this = *A_VC_Star;
        }
        else
        {
            auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
            ( new DistMatrix<T,Star,VC>(g) );
            *A_Star_VC = A;

            auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
            ( new DistMatrix<T,Star,VR>(true,this->RowAlignment(),g) );
            *A_Star_VR = *A_Star_VC;
            delete A_Star_VC.release(); // lowers memory highwater

            *this = *A_Star_VR;
            this->ResizeTo( A_Star_VR->Height(), A_Star_VR->Width() );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [MR,* ]");
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
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();

    auto_ptr< DistMatrix<T,Star,VC> > A_Star_VC
    ( new DistMatrix<T,Star,VC>(g) );
    *A_Star_VC = A;

    auto_ptr< DistMatrix<T,Star,VR> > A_Star_VR
    ( new DistMatrix<T,Star,VR>(true,this->RowAlignment(),g) );
    *A_Star_VR = *A_Star_VC;
    delete A_Star_VC.release(); // lowers memory highwater

    *this = *A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [VC,* ]");
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
        const int colShift = this->ColShift();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,c);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        this->_auxMemory.Require( 2*c*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[c*portionSize];

        // Pack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for  
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignment,c);
            const int thisLocalWidth = LocalLength(width,thisRowShift,c);

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for 
# endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*c);
                T* dataCol = &(data[j*localHeightOfA]);
                memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for  COLLAPSE(2)
# endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = 
                        A.GetLocalEntry(i,thisRowShift+j*c);
#endif
        }

        // Communicate
        AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, g.MRComm() );

        // Unpack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for 
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisRank = row+k*r;
            const int thisColShift = Shift(thisRank,colAlignmentOfA,p);
            const int thisColOffset = (thisColShift-colShift) / r;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    this->SetLocalEntry
                        (thisColOffset+i*c,j,data[i+j*thisLocalHeight]);
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MC,MR] <- [VC,* ]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int row = g.MCRank();
        const int colShift = this->ColShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int colAlignmentOfA = A.ColAlignment();

        const int sendRow = (row+r+colAlignment-(colAlignmentOfA%r)) % r;
        const int recvRow = (row+r+(colAlignmentOfA%r)-colAlignment) % r;

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = this->LocalWidth();
        const int localHeightOfA = A.LocalHeight();

        const int maxHeight = MaxLocalLength(height,p);
        const int maxWidth = MaxLocalLength(width,c);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        this->_auxMemory.Require( 2*c*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[c*portionSize];

        // Pack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignment,c);
            const int thisLocalWidth = LocalLength(width,thisRowShift,c);

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*c);
                T* dataCol = &(data[j*localHeightOfA]);
                memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeightOfA; ++i )
                    data[i+j*localHeightOfA] = 
                        A.GetLocalEntry(i,thisRowShift+j*c);
#endif
        }

        // SendRecv: properly align A[VC,*] via a trade in the column
        SendRecv
        ( secondBuffer, c*portionSize, sendRow, 0,
          firstBuffer,  c*portionSize, recvRow, 0, g.MCComm() );

        // AllToAll to gather all of the aligned A[VC,*] data into secondBuff.
        AllToAll
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MRComm() );

        // Unpack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisRank = recvRow+k*r;
            const int thisColShift = Shift(thisRank,colAlignmentOfA,p);
            const int thisColOffset = (thisColShift-colShift) / r;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    this->SetLocalEntry
                        (thisColOffset+i*c,j,data[i+j*thisLocalHeight]);
        }

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,VC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
    DistMatrix<T,Star,VR> A_Star_VR(true,this->RowAlignment(),g);

    A_Star_VR = A;
    *this = A_Star_VR;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [VR,* ]");
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
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,VR]");
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
            this->_rowAlignment = A.RowAlignment() % g.Width();
            this->_rowShift = 
                Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() % g.Width() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();
        const int rowShift = this->RowShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();

        const int maxHeight = MaxLocalLength(height,r);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        this->_auxMemory.Require( 2*r*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[r*portionSize];

        // Pack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignment,r);
            const int thisLocalHeight = LocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.GetLocalEntry(thisColShift+i*r,j);
        }

        // Communicate
        AllToAll
        ( sendBuffer, portionSize,
          recvBuffer, portionSize, g.MCComm() );

        // Unpack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &recvBuffer[k*portionSize];

            const int thisRank = col+k*c;
            const int thisRowShift = Shift(thisRank,rowAlignmentOfA,p);
            const int thisRowOffset = (thisRowShift-rowShift) / c;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* dataCol = &(data[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,thisRowOffset+j*r);
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry
                        (i,thisRowOffset+j*r,data[i+j*localHeight]);
#endif
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MC,MR] <- [* ,VR]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int p = r * c;
        const int col = g.MRRank();
        const int rowShift = this->RowShift();
        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();

        const int sendCol = (col+c+rowAlignment-(rowAlignmentOfA%c)) % c;
        const int recvCol = (col+c+(rowAlignmentOfA%c)-rowAlignment) % c;

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        
        const int maxHeight = MaxLocalLength(height,r);
        const int maxWidth = MaxLocalLength(width,p);
        const int portionSize = max(maxHeight*maxWidth,MinCollectContrib);

        this->_auxMemory.Require( 2*r*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[r*portionSize];

        // Pack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignment,r);
            const int thisLocalHeight = LocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.GetLocalEntry(thisColShift+i*r,j);
        }

        // SendRecv: properly align A[*,VR] via a trade in the column
        SendRecv
        ( secondBuffer, r*portionSize, sendCol, 0,
          firstBuffer,  r*portionSize, recvCol, 0, g.MRComm() );

        // AllToAll to gather all of the aligned [*,VR] data into secondBuffer
        AllToAll
        ( firstBuffer,  portionSize,
          secondBuffer, portionSize, g.MCComm() );

        // Unpack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisRank = recvCol+k*c;
            const int thisRowShift = Shift(thisRank,rowAlignmentOfA,p);
            const int thisRowOffset = (thisRowShift-rowShift) / c;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* dataCol = &(data[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,thisRowOffset+j*r);
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry
                        (i,thisRowOffset+j*r,data[i+j*localHeight]);
#endif
        }

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,MC,MR>&
elemental::DistMatrixBase<T,MC,MR>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const int r = this->GetGrid().Height();
    const int c = this->GetGrid().Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->SetLocalEntry(i,j,A.GetLocalEntry(colShift+i*r,rowShift+j*c));
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::SumScatterFrom
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterFrom([MC,* ])");
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
        if( this->Width() == 1 )
        {
            const int rowAlignment = this->RowAlignment();
            const int myCol = g.MRRank();

            const int localHeight = this->LocalHeight();

            const int recvSize = max(localHeight,MinCollectContrib);
            const int sendSize = recvSize;

            this->_auxMemory.Require( sendSize + recvSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];

            // Pack 
#ifdef RELEASE
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeight*sizeof(T) );
#else
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int i=0; i<localHeight; ++i )
                sendBuffer[i] = A.GetLocalEntry(i,0);
#endif

            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, sendSize, 
              MPI_SUM, rowAlignment, g.MRComm() );

            if( myCol == rowAlignment )
            {
#ifdef RELEASE
                T* thisCol = this->LocalBuffer(0,0);
                memcpy( thisCol, recvBuffer, localHeight*sizeof(T) );
#else
# ifdef _OPENMP
                #pragma omp parallel for
# endif
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(i,0,recvBuffer[i]);
#endif
            }

            this->_auxMemory.Release();
        }
        else
        {
            const int c = g.Width();
            const int rowAlignment = this->RowAlignment();
        
            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,c);

            const int recvSize = 
                max(localHeight*maxLocalWidth,MinCollectContrib);
            const int sendSize = c * recvSize;

            this->_auxMemory.Require( sendSize + recvSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[sendSize];
        
            // Pack 
            vector<int> recvSizes(c);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<c; ++k )
            {
                T* data = &sendBuffer[k*recvSize];
                recvSizes[k] = recvSize;

                const int thisRowShift = Shift( k, rowAlignment, c );
                const int thisLocalWidth = LocalLength(width,thisRowShift,c);

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
# endif
                for( int j=0; j<thisLocalWidth; ++j )
                {
                    const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*c);
                    T* dataCol = &(data[j*localHeight]);
                    memcpy( dataCol, ACol, localHeight*sizeof(T) );
                }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
# endif
                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        data[i+j*localHeight] = 
                            A.GetLocalEntry(i,thisRowShift+j*c);
#endif
            }

            // Reduce-scatter over each process row
            ReduceScatter
            ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MRComm() );

            // Unpack our received data
#ifdef RELEASE
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* recvBufferCol = &(recvBuffer[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,j);
                memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
#else
# ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(i,j,recvBuffer[i+j*localHeight]);
#endif

            this->_auxMemory.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned SumScatterFrom [MC,MR] <- [MC,* ]." << endl;
#endif
        if( this->Width() == 1 )
        {
            const int r = g.Height();
            const int rowAlignment = this->RowAlignment();
            const int myRow = g.MCRank();
            const int myCol = g.MRRank();

            const int height = this->Height();
            const int localHeight = this->LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int portionSize = max(maxLocalHeight,MinCollectContrib);

            const int colAlignment = this->ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendRow = (myRow+r+colAlignment-colAlignmentOfA) % r;
            const int recvRow = (myRow+r+colAlignmentOfA-colAlignment) % r;

            this->_auxMemory.Require( 2*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack 
#ifdef RELEASE
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeightOfA*sizeof(T) );
#else
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i] = A.GetLocalEntry(i,0);
#endif
        
            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize, 
              MPI_SUM, rowAlignment, g.MRComm() );

            if( myCol == rowAlignment )
            {
                // Perform the realignment
                SendRecv
                ( recvBuffer, portionSize, sendRow, 0,
                  sendBuffer, portionSize, recvRow, 0, g.MCComm() );

#ifdef RELEASE
                T* thisCol = this->LocalBuffer(0,0);
                memcpy( thisCol, sendBuffer, localHeight*sizeof(T) );
#else
# ifdef _OPENMP
                #pragma omp parallel for
# endif
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(i,0,sendBuffer[i]);
#endif
            }

            this->_auxMemory.Release();
        }
        else
        {
            const int r = g.Height();
            const int c = g.Width();
            const int row = g.MCRank();

            const int colAlignment = this->ColAlignment();
            const int rowAlignment = this->RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
            const int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalWidth = MaxLocalLength(width,c);

            const int recvSize_RS = 
                max(localHeightOfA*maxLocalWidth,MinCollectContrib);
            const int sendSize_RS = c * recvSize_RS;
            const int recvSize_SR = localHeight * localWidth;

            this->_auxMemory.Require
            ( recvSize_RS + max(sendSize_RS,recvSize_SR) );

            T* buffer = this->_auxMemory.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[recvSize_RS];

            // Pack 
            vector<int> recvSizes(c);
            for( int k=0; k<c; ++k )
            {
                T* data = &secondBuffer[k*recvSize_RS];
                recvSizes[k] = recvSize_RS;

                const int thisRowShift = Shift( k, rowAlignment, c );
                const int thisLocalWidth = LocalLength(width,thisRowShift,c);

#ifdef RELEASE
# ifdef _OPENMP
                #pragma omp parallel for
# endif
                for( int j=0; j<thisLocalWidth; ++j )
                {
                    const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*c);
                    T* dataCol = &(data[j*localHeightOfA]);
                    memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                }
#else
# ifdef _OPENMP
                #pragma omp parallel for COLLAPSE(2)
# endif
                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeightOfA; ++i )
                        data[i+j*localHeightOfA] = 
                            A.GetLocalEntry(i,thisRowShift+j*c);
#endif
            }

            // Reduce-scatter over each process row
            ReduceScatter
            ( secondBuffer, firstBuffer, &recvSizes[0], MPI_SUM, g.MRComm() );

            // Trade reduced data with the appropriate process row
            SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendRow, 0,
              secondBuffer, localHeight*localWidth,    recvRow, 0, g.MCComm() );

            // Unpack the received data
#ifdef RELEASE
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* secondBufferCol = &(secondBuffer[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,j);
                memcpy( thisCol, secondBufferCol, localHeight*sizeof(T) );
            }
#else
# ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(i,j,secondBuffer[i+j*localHeight]);
#endif

            this->_auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::SumScatterFrom
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterFrom([* ,MR])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
#ifdef VECTOR_WARNINGS
    if( A.Width() == 1 && g.VCRank() == 0 )
    {
        cerr <<
          "The vector version of [MC,MR].SumScatterFrom([* ,MR]) does not "
          "yet have a vector version implemented, but it would only require"
          " a modification of the vector version of "
          "[MC,MR].SumScatterFrom([MC,* ])" << endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.VCRank() == 0 )
    {
        cerr << 
          "[MC,MR]::SumScatterFrom([* ,MR]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MR] matrix instead." << endl;
    }
#endif
    if( !this->Viewing() )
    {
        if( !this->ConstrainedRowAlignment() )
        {
            this->_rowAlignment = A.RowAlignment();
            this->_rowShift = 
                Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        }
        this->ResizeTo( A.Height(), A.Width() );
    }

    if( this->RowAlignment() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int colAlignment = this->ColAlignment();

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);

        const int recvSize = 
            max(maxLocalHeight*localWidth,MinCollectContrib);
        const int sendSize = r * recvSize;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack 
        vector<int> recvSizes(r);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.GetLocalEntry(thisColShift+i*r,j);
        }

        // Reduce-scatter over each process col
        ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

        // Unpack our received data
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,recvBuffer[i+j*localHeight]);
#endif

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned SumScatterFrom [MC,MR] <- [* ,MR]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);

        const int recvSize_RS = 
            max(maxLocalHeight*localWidthOfA,MinCollectContrib);
        const int sendSize_RS = r * recvSize_RS;
        const int recvSize_SR = localHeight * localWidth;

        this->_auxMemory.Require( recvSize_RS + max(sendSize_RS,recvSize_SR) );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[recvSize_RS];

        // Pack 
        vector<int> recvSizes(r);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &secondBuffer[k*recvSize_RS];
            recvSizes[k] = recvSize_RS;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.GetLocalEntry(thisColShift+i*r,j);
        }

        // Reduce-scatter over each process col
        ReduceScatter
        ( secondBuffer, firstBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

        // Trade reduced data with the appropriate process col
        SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendCol, 0,
          secondBuffer, localHeight*localWidth,    recvCol, 0, g.MRComm() );

        // Unpack the received data
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* secondBufferCol = &(secondBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, secondBufferCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,secondBuffer[i+j*localHeight]);
#endif
        
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::SumScatterFrom
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterFrom([* ,* ])");
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
    const int maxLocalHeight = MaxLocalLength(height,r);
    const int maxLocalWidth = MaxLocalLength(width,c);

    const int recvSize = max(maxLocalHeight*maxLocalWidth,MinCollectContrib);
    const int sendSize = r * c * recvSize;

    this->_auxMemory.Require( sendSize + recvSize );

    T* buffer = this->_auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack 
    vector<int> recvSizes(r*c);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int l=0; l<c; ++l )
    {
        const int thisRowShift = Shift( l, rowAlignment, c );
        const int thisLocalWidth = LocalLength( width, thisRowShift, c );

        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[(k+l*r)*recvSize];
            recvSizes[k+l*r] = recvSize;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.GetLocalEntry(thisColShift+i*r,thisRowShift+j*c);
        }
    }

    // Reduce-scatter over each process col
    ReduceScatter
    ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.VCComm() );

    // Unpack our received data
#ifdef RELEASE
# ifdef _OPENMP
    #pragma omp parallel for
# endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* recvBufferCol = &(recvBuffer[j*localHeight]);
        T* thisCol = this->LocalBuffer(0,j);
        memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
    }
#else
# ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
# endif
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->SetLocalEntry(i,j,recvBuffer[i+j*localHeight]);
#endif

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::SumScatterUpdate
( T alpha, const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterUpdate([MC,* ])");
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
            const int myCol = g.MRRank();

            const int localHeight = this->LocalHeight();

            const int portionSize = max(localHeight,MinCollectContrib);

            this->_auxMemory.Require( 2*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack 
#ifdef RELEASE
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeight*sizeof(T) );
#else
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int i=0; i<localHeight; ++i )
                sendBuffer[i] = A.GetLocalEntry(i,0);
#endif
        
            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize, 
              MPI_SUM, rowAlignment, g.MRComm() );

            if( myCol == rowAlignment )
            {
#ifdef RELEASE
                T* thisCol = this->LocalBuffer(0,0);
# ifdef _OPENMP
                #pragma omp parallel for
# endif
                for( int i=0; i<localHeight; ++i )
                    thisCol[i] += alpha*recvBuffer[i];
#else
# ifdef _OPENMP
                #pragma omp parallel for
# endif
                for( int i=0; i<localHeight; ++i )
                {
                    const T value = this->GetLocalEntry(i,0);
                    this->SetLocalEntry(i,0,value+alpha*recvBuffer[i]);
                }
#endif
            }

            this->_auxMemory.Release();
        }
        else
        {
            const int c = g.Width();
            const int rowAlignment = this->RowAlignment();

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int maxLocalWidth = MaxLocalLength(width,c);

            const int portionSize = 
                max(localHeight*maxLocalWidth,MinCollectContrib);

            this->_auxMemory.Require( (c+1)*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[c*portionSize];

            // Pack 
            vector<int> recvSizes(c);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<c; ++k )
            {
                T* data = &sendBuffer[k*portionSize];
                recvSizes[k] = portionSize;

                const int thisRowShift = Shift( k, rowAlignment, c );
                const int thisLocalWidth = LocalLength(width,thisRowShift,c);

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
# endif
                for( int j=0; j<thisLocalWidth; ++j )
                {
                    const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*c);
                    T* dataCol = &(data[j*localHeight]);
                    memcpy( dataCol, ACol, localHeight*sizeof(T) );
                }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
# endif
                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        data[i+j*localHeight] = 
                            A.GetLocalEntry(i,thisRowShift+j*c);
#endif
            }
        
            // Reduce-scatter over each process row
            ReduceScatter
            ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MRComm() );

            // Update with our received data
#ifdef RELEASE
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* recvBufferCol = &(recvBuffer[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,j);
                wrappers::blas::Axpy
                ( localHeight, alpha, recvBufferCol, 1, thisCol, 1 );
            }
#else
# ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
            {
                for( int i=0; i<localHeight; ++i )
                {
                    const T value = this->GetLocalEntry(i,j);
                    this->SetLocalEntry
                        (i,j,value+alpha*recvBuffer[i+j*localHeight]);
                }
            }
#endif

            this->_auxMemory.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned SumScatterUpdate [MC,MR] <- [MC,* ]." << endl;
#endif
        if( this->Width() == 1 )
        {
            const int r = g.Height();
            const int rowAlignment = this->RowAlignment();
            const int myRow = g.MCRank();
            const int myCol = g.MRRank();

            const int height = this->Height();
            const int localHeight = this->LocalHeight();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalHeight = MaxLocalLength(height,r);

            const int portionSize = max(maxLocalHeight,MinCollectContrib);

            const int colAlignment = this->ColAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendRow = (myRow+r+colAlignment-colAlignmentOfA) % r;
            const int recvRow = (myRow+r+colAlignmentOfA-colAlignment) % r;

            this->_auxMemory.Require( 2*portionSize );

            T* buffer = this->_auxMemory.Buffer();
            T* sendBuffer = &buffer[0];
            T* recvBuffer = &buffer[portionSize];

            // Pack 
#ifdef RELEASE
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeightOfA*sizeof(T) );
#else
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int i=0; i<localHeightOfA; ++i )
                sendBuffer[i] = A.GetLocalEntry(i,0);
#endif
        
            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize, 
              MPI_SUM, rowAlignment, g.MRComm() );

            if( myCol == rowAlignment )
            {
                // Perform the realignment
                SendRecv
                ( recvBuffer, portionSize, sendRow, 0,
                  sendBuffer, portionSize, recvRow, 0, g.MCComm() );

#ifdef RELEASE
                T* thisCol = this->LocalBuffer(0,0);
# ifdef _OPENMP
                #pragma omp parallel for
# endif
                for( int i=0; i<localHeight; ++i )
                    thisCol[i] += alpha*sendBuffer[i];
#else
                for( int i=0; i<localHeight; ++i )
                {
                    const T value = this->GetLocalEntry(i,0);
                    this->SetLocalEntry(i,0,value+alpha*sendBuffer[i]);
                }
#endif
            }

            this->_auxMemory.Release();
        }
        else
        {
            const int r = g.Height();
            const int c = g.Width();
            const int row = g.MCRank();

            const int colAlignment = this->ColAlignment();
            const int rowAlignment = this->RowAlignment();
            const int colAlignmentOfA = A.ColAlignment();
            const int sendRow = (row+r+colAlignment-colAlignmentOfA) % r;
            const int recvRow = (row+r+colAlignmentOfA-colAlignment) % r;

            const int width = this->Width();
            const int localHeight = this->LocalHeight();
            const int localWidth = this->LocalWidth();
            const int localHeightOfA = A.LocalHeight();
            const int maxLocalWidth = MaxLocalLength(width,c);

            const int recvSize_RS = 
                max(localHeightOfA*maxLocalWidth,MinCollectContrib);
            const int sendSize_RS = c * recvSize_RS;
            const int recvSize_SR = localHeight * localWidth;

            this->_auxMemory.Require
            ( recvSize_RS + max(sendSize_RS,recvSize_SR) );

            T* buffer = this->_auxMemory.Buffer();
            T* firstBuffer = &buffer[0];
            T* secondBuffer = &buffer[recvSize_RS];

            // Pack 
            vector<int> recvSizes(c);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<c; ++k )
            {
                T* data = &secondBuffer[k*recvSize_RS];
                recvSizes[k] = recvSize_RS;

                const int thisRowShift = Shift( k, rowAlignment, c );
                const int thisLocalWidth = LocalLength(width,thisRowShift,c);

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
# endif
                for( int j=0; j<thisLocalWidth; ++j )
                {
                    const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*c);
                    T* dataCol = &(data[j*localHeightOfA]);
                    memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
# endif
                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<localHeightOfA; ++i )
                        data[i+j*localHeightOfA] = 
                            A.GetLocalEntry(i,thisRowShift+j*c);
#endif
            }

            // Reduce-scatter over each process row
            ReduceScatter
            ( secondBuffer, firstBuffer, &recvSizes[0], MPI_SUM, g.MRComm() );

            // Trade reduced data with the appropriate process row
            SendRecv
            ( firstBuffer,  localHeightOfA*localWidth, sendRow, 0,
              secondBuffer, localHeight*localWidth,    recvRow, 0, 
              g.MCComm() );

            // Update with our received data
#ifdef RELEASE
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* secondBufferCol = &(secondBuffer[j*localHeight]);
                T* thisCol = this->LocalBuffer(0,j);
                for( int i=0; i<localHeight; ++i )
                    thisCol[i] += alpha*secondBufferCol[i];
            }
#else
# ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
            {
                for( int i=0; i<localHeight; ++i )
                {
                    const T value = this->GetLocalEntry(i,j);
                    this->SetLocalEntry
                        (i,j,value+alpha*secondBuffer[i+j*localHeight]);
                }
            }
#endif

            this->_auxMemory.Release();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::SumScatterUpdate
( T alpha, const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterUpdate([* ,MR])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const Grid& g = this->GetGrid();
#ifdef VECTOR_WARNINGS
    if( A.Width() == 1 && g.VCRank() == 0 )
    {
        cerr <<
          "The vector version of [MC,MR].SumScatterUpdate([* ,MR]) does not"
          " yet have a vector version implemented, but it would only "
          "require a modification of the vector version of "
          "[MC,MR].SumScatterUpdate([MC,* ])" << endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Width() != 1 && g.VCRank() == 0 )
    {
        cerr << 
          "[MC,MR]::SumScatterUpdate([* ,MR]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MR] matrix instead." << endl;
    }
#endif
    if( this->RowAlignment() == A.RowAlignment() )
    {
        const int r = g.Height();
        const int colAlignment = this->ColAlignment();

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);

        const int recvSize = max(maxLocalHeight*localWidth,MinCollectContrib);
        const int sendSize = r * recvSize;

        this->_auxMemory.Require( sendSize + recvSize );

        T* buffer = this->_auxMemory.Buffer();
        T* sendBuffer = &buffer[0];
        T* recvBuffer = &buffer[sendSize];

        // Pack 
        vector<int> recvSizes(r);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.GetLocalEntry(thisColShift+i*r,j);
        }

        // Reduce-scatter over each process col
        ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

        // Update with our received data
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &(recvBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            for( int i=0; i<localHeight; ++i )
                thisCol[i] += alpha*recvBufferCol[i];
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
        {
            for( int i=0; i<localHeight; ++i )
            {
                const T value = this->GetLocalEntry(i,j);
                this->SetLocalEntry
                    (i,j,value+alpha*recvBuffer[i+j*localHeight]);
            }
        }
#endif

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned SumScatterUpdate [MC,MR] <- [* ,MR]." << endl;
#endif
        const int r = g.Height();
        const int c = g.Width();
        const int col = g.MRRank();

        const int colAlignment = this->ColAlignment();
        const int rowAlignment = this->RowAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
        const int sendCol = (col+c+rowAlignment-rowAlignmentOfA) % c;
        const int recvCol = (col+c+rowAlignmentOfA-rowAlignment) % c;

        const int height = this->Height();
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);

        const int recvSize_RS = 
            max(maxLocalHeight*localWidthOfA,MinCollectContrib);
        const int sendSize_RS = r * recvSize_RS;
        const int recvSize_SR = localHeight * localWidth;

        this->_auxMemory.Require( recvSize_RS + max(sendSize_RS,recvSize_SR) );

        T* buffer = this->_auxMemory.Buffer();
        T* firstBuffer = &buffer[0];
        T* secondBuffer = &buffer[recvSize_RS];

        // Pack
        vector<int> recvSizes(r);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &secondBuffer[k*recvSize_RS];
            recvSizes[k] = recvSize_RS;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.GetLocalEntry(thisColShift+i*r,j);
        }

        // Reduce-scatter over each process col
        ReduceScatter
        ( secondBuffer, firstBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

        // Trade reduced data with the appropriate process col
        SendRecv
        ( firstBuffer,  localHeight*localWidthOfA, sendCol, 0,
          secondBuffer, localHeight*localWidth,    recvCol, MPI_ANY_TAG,
          g.MRComm() );

        // Update with our received data
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* secondBufferCol = &(secondBuffer[j*localHeight]);
            T* thisCol = this->LocalBuffer(0,j);
            for( int i=0; i<localHeight; ++i )
                thisCol[i] += alpha*secondBufferCol[i];
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
        {
            for( int i=0; i<localHeight; ++i )
            {
                const T value = this->GetLocalEntry(i,j);
                this->SetLocalEntry
                    (i,j,value+alpha*secondBuffer[i+j*localHeight]);
            }
        }
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MC,MR>::SumScatterUpdate
( T alpha, const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SumScatterUpdate([* ,* ])");
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
    const int maxLocalHeight = MaxLocalLength(height,r);
    const int maxLocalWidth = MaxLocalLength(width,c);

    const int recvSize = max(maxLocalHeight*maxLocalWidth,MinCollectContrib);
    const int sendSize = r * c * recvSize;

    this->_auxMemory.Require( sendSize + recvSize );

    T* buffer = this->_auxMemory.Buffer();
    T* sendBuffer = &buffer[0];
    T* recvBuffer = &buffer[sendSize];

    // Pack 
    vector<int> recvSizes(r*c);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int l=0; l<c; ++l )
    {
        const int thisRowShift = Shift( l, rowAlignment, c );
        const int thisLocalWidth = LocalLength( width, thisRowShift, c );

        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[(k+l*r)*recvSize];
            recvSizes[k+l*r] = recvSize;

            const int thisColShift = Shift( k, colAlignment, r );
            const int thisLocalHeight = LocalLength( height, thisColShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.GetLocalEntry(thisColShift+i*r,thisRowShift+j*c);
        }
    }

    // Reduce-scatter over each process col
    ReduceScatter
    ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.VCComm() );

    // Unpack our received data
#ifdef RELEASE
# ifdef _OPENMP
    #pragma omp parallel for
# endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* recvBufferCol = &(recvBuffer[j*localHeight]);
        T* thisCol = this->LocalBuffer(0,j);
        for( int i=0; i<localHeight; ++i )
            thisCol[i] += alpha*recvBufferCol[i];
    }
#else
# ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
# endif
    for( int j=0; j<localWidth; ++j )
    {
        for( int i=0; i<localHeight; ++i )
        {
            const T value = this->GetLocalEntry(i,j);
            this->SetLocalEntry(i,j,value+alpha*recvBuffer[i+j*localHeight]);
        }
    }
#endif

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
elemental::DistMatrix<R,MC,MR>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int r = this->GetGrid().Height();
    const int c = this->GetGrid().Width();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;                
        if( i % c == rowShift )
        {
            const int jLoc = (i-rowShift) / c;
            if( jLoc < localWidth )
            {
                const R value = this->GetLocalEntry(iLoc,jLoc);
                this->SetLocalEntry(iLoc,jLoc,value+this->Width());
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int r = this->GetGrid().Height();
    const int c = this->GetGrid().Width();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;                
        if( i % c == rowShift )
        {
            const int jLoc = (i-rowShift) / c;
            if( jLoc < localWidth )
            {
                const R value = real(this->GetLocalEntry(iLoc,jLoc));
                this->SetLocalEntry(iLoc,jLoc,value+this->Width());
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,MC,MR>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = real(this->GetLocalEntry(iLoc,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,MC,MR>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = imag(this->GetLocalEntry(iLoc,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const int jLoc = (j-this->RowShift()) / g.Width();
        const R v = imag(this->GetLocalEntry(iLoc,jLoc));
        this->SetLocalEntry(iLoc,jLoc,complex<R>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetImag
( int i, int j, R v )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const int jLoc = (j-this->RowShift()) / g.Width();
        R u = real(this->GetLocalEntry(iLoc,jLoc));
        this->SetLocalEntry(iLoc,jLoc,complex<R>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::GetRealDiagonal
( DistMatrix<R,MD,Star>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetRealDiagonal");
    this->AssertNotLockedView();
#endif
    int width = this->Width();
    int height = this->Height();
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
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:" << endl
            << "  A ~ " << this->Height() << " x " << this->Width() << endl
            << "  d ~ " << d.Height() << " x " << d.Width() << endl
            << "  A diag length: " << length << endl;
        throw logic_error( msg.str() );
    }
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R u = real(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            d.SetLocalEntry(k,0,u);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::GetImagDiagonal
( DistMatrix<R,MD,Star>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImagDiagonal");
    this->AssertNotLockedView();
#endif
    int width = this->Width();
    int height = this->Height();
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
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:" << endl
            << "  A ~ " << this->Height() << " x " << this->Width() << endl
            << "  d ~ " << d.Height() << " x " << d.Width() << endl
            << "  A diag length: " << length << endl;
        throw logic_error( msg.str() );
    }
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R v = imag(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            d.SetLocalEntry(k,0,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::GetRealDiagonal
( DistMatrix<R,Star,MD>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetRealDiagonal");
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
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:" << endl
            << "  A ~ " << this->Height() << " x " << this->Width() << endl
            << "  d ~ " << d.Height() << " x " << d.Width() << endl
            << "  A diag length: " << length << endl;
        throw logic_error( msg.str() );
    }
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R u = real(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            d.SetLocalEntry(0,k,u);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::GetImagDiagonal
( DistMatrix<R,Star,MD>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImagDiagonal");
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
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:" << endl
            << "  A ~ " << this->Height() << " x " << this->Width() << endl
            << "  d ~ " << d.Height() << " x " << d.Width() << endl
            << "  A diag length: " << length << endl;
        throw logic_error( msg.str() );
    }
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R v = imag(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            d.SetLocalEntry(0,k,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetDiagonal
( const DistMatrixBase<R,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetDiagonal");
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),d.GetLocalEntry(k,0));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetRealDiagonal
( const DistMatrixBase<R,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetRealDiagonal");
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R v = imag(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),
                 complex<R>(d.GetLocalEntry(k,0),v));
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetImagDiagonal
( const DistMatrixBase<R,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImagDiagonal");
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R u = real(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),
                 complex<R>(u,d.GetLocalEntry(k,0)));
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetDiagonal
( const DistMatrixBase<R,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetDiagonal");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        const int height = this->Height();
        const int width = this->Width();

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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),d.GetLocalEntry(0,k));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetRealDiagonal
( const DistMatrixBase<R,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetRealDiagonal");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        const int height = this->Height();
        const int width = this->Width();

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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R v = imag(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),
                 complex<R>(d.GetLocalEntry(0,k),v));
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetImagDiagonal
( const DistMatrixBase<R,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImagDiagonal");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        const int height = this->Height();
        const int width = this->Width();

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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R u = real(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),
                 complex<R>(u,d.GetLocalEntry(0,k)));
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrixBase<int,   MC,MR>;
template class elemental::DistMatrixBase<float, MC,MR>;
template class elemental::DistMatrixBase<double,MC,MR>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,MC,MR>;
template class elemental::DistMatrixBase<dcomplex,MC,MR>;
#endif

template class elemental::DistMatrix<int,   MC,MR>;
template class elemental::DistMatrix<float, MC,MR>;
template class elemental::DistMatrix<double,MC,MR>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,MC,MR>;
template class elemental::DistMatrix<dcomplex,MC,MR>;
#endif

