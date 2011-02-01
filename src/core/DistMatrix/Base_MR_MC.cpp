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
elemental::DistMatrixBase<T,MR,MC>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Print");
#endif
    const Grid& g = this->Grid();
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
    vector<T> sendBuf(height*width,0);
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            sendBuf[colShift+i*c + (rowShift+j*r)*height] = 
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
elemental::DistMatrixBase<T,MR,MC>::Align
( int colAlignment, int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Align");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
#endif
    const Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Width() )
        throw std::runtime_error( "Invalid column alignment for [MR,MC]" );
    if( rowAlignment < 0 || rowAlignment >= g.Height() )
        throw std::runtime_error( "Invalid row alignment for [MR,MC]" );
#endif
    this->_colAlignment = colAlignment;
    this->_rowAlignment = rowAlignment;
    this->_colShift = Shift( g.MRRank(), colAlignment, g.Width() );
    this->_rowShift = Shift( g.MCRank(), rowAlignment, g.Height() );
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
elemental::DistMatrixBase<T,MR,MC>::AlignCols
( int colAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignCols");
    this->AssertFreeColAlignment();
#endif
    const Grid& g = this->Grid();
#ifndef RELEASE
    if( colAlignment < 0 || colAlignment >= g.Width() )
        throw std::runtime_error( "Invalid column alignment for [MR,MC]" );
#endif
    this->_colAlignment = colAlignment;
    this->_colShift = Shift( g.MRRank(), colAlignment, g.Width() );
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
elemental::DistMatrixBase<T,MR,MC>::AlignRows
( int rowAlignment )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRows");
    this->AssertFreeRowAlignment();
#endif
    const Grid& g = this->Grid();
#ifndef RELEASE
    if( rowAlignment < 0 || rowAlignment >= g.Height() )
        throw std::runtime_error( "Invalid row alignment for [MR,MC]" );
#endif
    this->_rowAlignment = rowAlignment;
    this->_rowShift = Shift( g.MCRank(), rowAlignment, g.Height() );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    const Grid& g = this->Grid();
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = 
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
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
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignWith([* ,VC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = 
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
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
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([VR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
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
elemental::DistMatrixBase<T,MR,MC>::AlignWith
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignWith([* ,VR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    this->_colAlignment = A.RowAlignment();
    this->_colShift = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    const Grid& g = this->Grid();
    this->_colAlignment = A.ColAlignment();
    this->_colShift = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
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
elemental::DistMatrixBase<T,MR,MC>::AlignColsWith
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignColsWith([* ,VR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    this->_colAlignment = A.RowAlignment();
    this->_colShift = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
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
    const Grid& g = this->Grid();
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = 
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
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
elemental::DistMatrixBase<T,MR,MC>::AlignRowsWith
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignRowsWith([* ,VC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = 
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
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
        const Grid& g = this->Grid();
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
        const Grid& g = this->Grid();
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
    const Grid& g = this->Grid();
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
    const Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    T u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
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
elemental::DistMatrixBase<T,MR,MC>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Set");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        this->SetLocalEntry(iLoc,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::Update
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::Update");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateLocalEntry(iLoc,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::GetDiagonal
( DistMatrixBase<T,MD,Star>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetDiagonal([MD,* ])");
    this->AssertNotLockedView();
#endif
    int length = this->DiagonalLength( offset );
#ifndef RELEASE
    if( d.Viewing() && length != d.Height() )
        throw logic_error( "d is not of the correct length." );
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
        const Grid& g = this->Grid();
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const T value = this->GetLocalEntry(iLocStart+k*(lcm/c),
                                                jLocStart+k*(lcm/r));
            d.SetLocalEntry(k,0,value);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::GetDiagonal
( DistMatrixBase<T,Star,MD>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetDiagonal([* ,MD])");
    this->AssertNotLockedView();
#endif
    int length = this->DiagonalLength( offset );
#ifndef RELEASE
    if( d.Viewing() && length != d.Width() )
        throw logic_error( "d is not of the correct length." );
    if( ( d.Viewing() && d.ConstrainedRowAlignment() ) &&
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
        const Grid& g = this->Grid();
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const T value = this->GetLocalEntry(iLocStart+k*(lcm/c),
                                                jLocStart+k*(lcm/r));
            d.SetLocalEntry(0,k,value);
        }
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
        const Grid& g = this->Grid();
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
            this->SetLocalEntry
                (iLocStart+k*(lcm/c),jLocStart+k*(lcm/r),d.GetLocalEntry(k,0));
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
        const Grid& g = this->Grid();
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
            this->SetLocalEntry
                (iLocStart+k*(lcm/c),jLocStart+k*(lcm/r),d.GetLocalEntry(0,k));
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
    const Grid& g = this->Grid();
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*r;
            int lastZeroRow = ( side==Left ? j-offset-1
                                           : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                int boundary = min( lastZeroRow+1, height );
                int numZeroRows = LocalLength( boundary, colShift, c );
                T* thisCol = this->LocalBuffer(0,jLoc);
                memset( thisCol, 0, numZeroRows*sizeof(T) );
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
            int j = rowShift + jLoc*r;
            int firstZeroRow = ( side==Left ? max(j-offset+1,0)
                                            : max(j-offset+height-width+1,0) );
            int numNonzeroRows = LocalLength(firstZeroRow,colShift,c);
            if( numNonzeroRows < localHeight )
            {
                T* thisCol = this->LocalBuffer(numNonzeroRows,jLoc);
                memset( thisCol, 0, (localHeight-numNonzeroRows)*sizeof(T) );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,MR,MC>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::ScaleTrapezoidal");
    this->AssertNotLockedView(); 
#endif
    const Grid& g = this->Grid();
    const int height = this->Height();
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    if( shape == Upper )
    {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*r;
            int lastRow = ( side==Left ? j-offset : j-offset+height-width );
            int boundary = min( lastRow+1, height );
            int numRows = LocalLength( boundary, colShift, c );
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
            int j = rowShift + jLoc*r;
            int firstRow = ( side==Left ? max(j-offset,0)
                                        : max(j-offset+height-width,0) );
            int numZeroRows = LocalLength( firstRow, colShift, c );
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
elemental::DistMatrixBase<T,MR,MC>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const Grid& g = this->Grid();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToZero();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*c;
        if( i % r == rowShift )
        {
            const int jLoc = (i-rowShift) / r;
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
            this->SetLocalEntry(i,j,Random<T>());
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
    const Grid& g = this->Grid();
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
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int k=0; k<c; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const int shift = Shift(myRow+r*k,colAlignmentOfA,p);
                const int offset = (shift-A.ColShift()) / r;
                const int thisLocalHeight = LocalLength(height,shift,p);

                for( int i=0; i<thisLocalHeight; ++i )
                    data[i] = A.GetLocalEntry(offset+i*c,0);
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
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const int shift = Shift(myCol+c*k,colAlignment,p);
                const int offset = (shift-this->ColShift()) / c;
                const int thisLocalHeight = LocalLength(height,shift,p);

                for( int i=0; i<thisLocalHeight; ++i )
                    this->SetLocalEntry(offset+i*r,0,data[i]);
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
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const int shift = Shift(myCol+c*k,rowAlignmentOfA,p);
                const int offset = (shift-A.RowShift()) / c;
                const int thisLocalWidth = LocalLength(width,shift,p);

                for( int j=0; j<thisLocalWidth; ++j )
                    data[j] = A.GetLocalEntry(0,offset+j*r);
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
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int k=0; k<c; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const int shift = Shift(myRow+r*k,rowAlignment,p);
                const int offset = (shift-this->RowShift()) / r;
                const int thisLocalWidth = LocalLength(width,shift,p);

                for( int j=0; j<thisLocalWidth; ++j )
                    this->SetLocalEntry(0,offset+j*c,data[j]);
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
    const Grid& g = this->Grid();
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
    const Grid& g = this->Grid();
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
        const Grid& g = this->Grid();
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,MC] <- [MR,MC]." << endl;
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* sendBufferCol = &sendBuffer[j*localHeightOfA];
            memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.VCComm() );

        // Unpack
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*localHeight];
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
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
( const DistMatrixBase<T,MR,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,rowShift+j*r);
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, ACol, localHeight*sizeof(T) );
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,MC] <- [MR,* ]." << endl;
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,rowShift+j*r);
            T* sendBufferCol = &sendBuffer[j*localHeightOfA];
            memcpy( sendBufferCol, ACol, localHeightOfA*sizeof(T) );
        }

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRank, 0,
          recvBuffer, recvSize, recvRank, MPI_ANY_TAG, g.MRComm() );

        // Unpack
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*localHeight];
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
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
( const DistMatrixBase<T,Star,MC>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
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
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,A.GetLocalEntry(colShift+i*c,j));
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,MC] <- [* ,MC]." << endl;
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
#ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeight; ++i )
                sendBuffer[i+j*localHeight] = A.GetLocalEntry(colShift+i*c,j);

        // Communicate
        SendRecv
        ( sendBuffer, sendSize, sendRow, 0,
          recvBuffer, recvSize, recvRow, MPI_ANY_TAG, g.MCComm() );

        // Unpack
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*localHeight];
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
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
( const DistMatrixBase<T,VC,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("[MR,MC] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
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
    const Grid& g = this->Grid();
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignment,c);
            const int thisLocalHeight = LocalLength(height,thisColShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.GetLocalEntry(thisColShift+i*c,j);
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
            const int thisRowShift = Shift(thisRank,rowAlignmentOfA,p);
            const int thisRowOffset = (thisRowShift-rowShift) / r;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* dataCol = &data[j*localHeight];
                T* thisCol = this->LocalBuffer(0,thisRowOffset+j*c);
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,MC] <- [* ,VC]." << endl;
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const int thisColShift = Shift(k,colAlignment,c);
            const int thisLocalHeight = LocalLength(height,thisColShift,c);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.GetLocalEntry(thisColShift+i*c,j);
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisRank = recvRow+k*r;
            const int thisRowShift = Shift(thisRank,rowAlignmentOfA,p);
            const int thisRowOffset = (thisRowShift-rowShift) / r;
            const int thisLocalWidth = LocalLength(width,thisRowShift,p);

#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* dataCol = &data[j*localHeight];
                T* thisCol = this->LocalBuffer(0,thisRowOffset+j*c);
                memcpy( thisCol, dataCol, localHeight*sizeof(T) );
            }
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
    const Grid& g = this->Grid();
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &sendBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignment,r);
            const int thisLocalWidth = LocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*r);
                T* dataCol = &data[j*localHeightOfA];
                memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
            }
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
            const int thisColShift = Shift(thisRank,colAlignmentOfA,p);
            const int thisColOffset = (thisColShift-colShift) / c;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);
            
#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    this->SetLocalEntry
                        (thisColOffset+i*r,j,data[i+j*thisLocalHeight]);
        }

        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned [MR,MC] <- [* ,VC]." << endl;
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            T* data = &secondBuffer[k*portionSize];

            const int thisRowShift = Shift(k,rowAlignment,r);
            const int thisLocalWidth = LocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int j=0; j<thisLocalWidth; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*r);
                T* dataCol = &data[j*localHeightOfA];
                memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
            }
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &secondBuffer[k*portionSize];

            const int thisRank = recvCol+k*c;
            const int thisColShift = Shift(thisRank,colAlignmentOfA,p);
            const int thisColOffset = (thisColShift-colShift) / c;
            const int thisLocalHeight = LocalLength(height,thisColShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    this->SetLocalEntry
                        (thisColOffset+i*r,j,data[i+j*thisLocalHeight]);
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
    const Grid& g = this->Grid();
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

    const Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
#ifdef _OPENMP
    #pragma omp parallel for COLLAPSE(2)
#endif
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            this->SetLocalEntry(i,j,A.GetLocalEntry(colShift+i*c,rowShift+j*r));
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
    const Grid& g = this->Grid();
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
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeight*sizeof(T) );

            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize,
              MPI_SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
            {
                T* thisCol = this->LocalBuffer(0,0);
                memcpy( thisCol, recvBuffer, localHeight*sizeof(T) );
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &sendBuffer[k*recvSize];
                recvSizes[k] = recvSize;

                const int thisRowShift = Shift( k, rowAlignment, r );
                const int thisLocalWidth = 
                      LocalLength( width, thisRowShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int j=0; j<thisLocalWidth; ++j )
                {
                    const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*r);
                    T* dataCol = &data[j*localHeight];
                    memcpy( dataCol, ACol, localHeight*sizeof(T) );
                }
            }

            // Reduce-scatter over each process column
            ReduceScatter
            ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

            // Unpack our received data
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* recvBufferCol = &recvBuffer[j*localHeight];
                T* thisCol = this->LocalBuffer(0,j);
                memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
            }
            this->_auxMemory.Release();
        }
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned SumScatterFrom [MR,MC] <- [MR,* ]." << endl;
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
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeightOfA*sizeof(T) );

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

                T* thisCol = this->LocalBuffer(0,0);
                memcpy( thisCol, sendBuffer, localHeight*sizeof(T) );
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &secondBuffer[k*recvSize_RS];
                recvSizes[k] = recvSize_RS;

                const int thisRowShift = Shift( k, rowAlignment, r );
                const int thisLocalWidth = LocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int j=0; j<thisLocalWidth; ++j )
                {
                    const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*r);
                    T* dataCol = &data[j*localHeightOfA];
                    memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                }
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
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* secondBufferCol = &secondBuffer[j*localHeight];
                T* thisCol = this->LocalBuffer(0,j);
                memcpy( thisCol, secondBufferCol, localHeight*sizeof(T) );
            }
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
#endif
    const Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Height() == 1 && g.VCRank() == 0 )
    {
        cerr <<    
          "The vector version of [MR,MC].SumScatterFrom([* ,MC]) is not yet"
          " written, but it only requires a modification of the vector "
          "version of [MR,MC].SumScatterFrom([MR,* ])." << endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.VCRank() == 0 )
    {
        cerr << 
          "[MR,MC]::SumScatterFrom([* ,MC]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MC] matrix instead." << endl;
    }
#endif
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisColShift = Shift( k, colAlignment, c );
            const int thisLocalHeight = LocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] =
                          A.GetLocalEntry(thisColShift+i*c,j);
        }

        // Reduce-scatter over each process row
        ReduceScatter
        ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MRComm() );

        // Unpack our received data
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufferCol = &recvBuffer[j*localHeight];
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
        }
        this->_auxMemory.Release();
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned SumScatterFrom [MR,MC] <- [* ,MC]." << endl;
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &secondBuffer[k*recvSize_RS];
            recvSizes[k] = recvSize_RS;

            const int thisColShift = Shift( k, colAlignment, c );
            const int thisLocalHeight = LocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] = 
                          A.GetLocalEntry(thisColShift+i*c,j);
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
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* secondBufferCol = &secondBuffer[j*localHeight];
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, secondBufferCol, localHeight*sizeof(T) );
        }
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

    const Grid& g = this->Grid();
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
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

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] =
                          A.GetLocalEntry(thisColShift+i*c,thisRowShift+j*r);
        }
    }

    // Reduce-scatter over each process col
    ReduceScatter
    ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.VRComm() );

    // Unpack our received data
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* recvBufferCol = &recvBuffer[j*localHeight];
        T* thisCol = this->LocalBuffer(0,j);
        memcpy( thisCol, recvBufferCol, localHeight*sizeof(T) );
    }
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
    const Grid& g = this->Grid();
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
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeight*sizeof(T) );

            // Reduce to rowAlignment
            Reduce
            ( sendBuffer, recvBuffer, portionSize,
              MPI_SUM, rowAlignment, g.MCComm() );

            if( myRow == rowAlignment )
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
            T* recvBuffer = &buffer[0];
            T* sendBuffer = &buffer[portionSize];

            // Pack 
            vector<int> recvSizes(r);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &sendBuffer[k*portionSize];
                recvSizes[k] = portionSize;

                const int thisRowShift = Shift( k, rowAlignment, r );
                const int thisLocalWidth = LocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int j=0; j<thisLocalWidth; ++j )
                {
                    const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*r);
                    T* dataCol = &data[j*localHeight];
                    memcpy( dataCol, ACol, localHeight*sizeof(T) );
                }
            }

            // Reduce-scatter over each process column
            ReduceScatter
            ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.MCComm() );

            // Update with our received data
#ifdef RELEASE
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* recvBufferCol = &recvBuffer[j*localHeight];
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
    }
    else
    {
#ifdef UNALIGNED_WARNINGS
        if( g.VCRank() == 0 )
            cerr << "Unaligned SumScatterUpdate [MR,MC] <- [MR,* ]." << endl;
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
            const T* ACol = A.LockedLocalBuffer(0,0);
            memcpy( sendBuffer, ACol, localHeightOfA*sizeof(T) );

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

#ifdef RELEASE
                T* thisCol = this->LocalBuffer(0,0);
# ifdef _OPENMP
                #pragma omp parallel for
# endif
                for( int i=0; i<localHeight; ++i )
                    thisCol[i] += alpha*sendBuffer[i];
#else
# ifdef _OPENMP
                #pragma omp parallel for
# endif
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
#endif
            for( int k=0; k<r; ++k )
            {
                T* data = &secondBuffer[k*recvSize_RS];
                recvSizes[k] = recvSize_RS;

                const int thisRowShift = Shift( k, rowAlignment, r );
                const int thisLocalWidth = LocalLength(width,thisRowShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
#endif
                for( int j=0; j<thisLocalWidth; ++j )
                {
                    const T* ACol = A.LockedLocalBuffer(0,thisRowShift+j*r);
                    T* dataCol = &data[j*localHeightOfA];
                    memcpy( dataCol, ACol, localHeightOfA*sizeof(T) );
                }
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
#ifdef RELEASE
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* secondBufferCol = &secondBuffer[j*localHeight];
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
elemental::DistMatrixBase<T,MR,MC>::SumScatterUpdate
( T alpha, const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SumScatterUpdate([* ,MC])");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    this->AssertSameSize( A );
#endif
    const Grid& g = this->Grid();
#ifdef VECTOR_WARNINGS
    if( A.Height() == 1 && g.VCRank() == 0 )
    {
        cerr <<    
          "The vector version of [MR,MC].SumScatterUpdate([* ,MC]) is not "
          "yet written, but it only requires a modification of the vector "
          "version of [MR,MC].SumScatterUpdate([MR,* ])." << endl;
    }
#endif
#ifdef CACHE_WARNINGS
    if( A.Height() != 1 && g.VCRank() == 0 )
    {
        cerr <<
          "[MR,MC]::SumScatterUpdate([* ,MC]) potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,MC] matrix instead." << endl;
    }
#endif
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &sendBuffer[k*recvSize];
            recvSizes[k] = recvSize;

            const int thisColShift = Shift( k, colAlignment, c );
            const int thisLocalHeight = LocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] =
                          A.GetLocalEntry(thisColShift+i*c,j);
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
            const T* recvBufferCol = &recvBuffer[j*localHeight];
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
            cerr << "Unaligned SumScatterUpdate [MR,MC] <- [* ,MC]." << endl;
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            T* data = &secondBuffer[k*recvSize_RS];
            recvSizes[k] = recvSize_RS;

            const int thisColShift = Shift( k, colAlignment, c );
            const int thisLocalHeight = LocalLength( height, thisColShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<localWidthOfA; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] =
                          A.GetLocalEntry(thisColShift+i*c,j);
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
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* secondBufferCol = &secondBuffer[j*localHeight];
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

    const Grid& g = this->Grid();
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
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
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

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<thisLocalWidth; ++j )
                for( int i=0; i<thisLocalHeight; ++i )
                    data[i+j*thisLocalHeight] =
                          A.GetLocalEntry(thisColShift+i*c,thisRowShift+j*r);
        }
    }

    // Reduce-scatter over each process col
    ReduceScatter
    ( sendBuffer, recvBuffer, &recvSizes[0], MPI_SUM, g.VRComm() );

    // Unpack our received data
#ifdef RELEASE
# ifdef _OPENMP
    #pragma omp parallel for
# endif
    for( int j=0; j<localWidth; ++j )
    {
        const T* recvBufferCol = &recvBuffer[j*localHeight];
        T* thisCol = this->LocalBuffer(0,j);
        imports::blas::Axpy
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
            this->SetLocalEntry(i,j,value+alpha*recvBuffer[i+j*localHeight]);
        }
    }
#endif

    this->_auxMemory.Release();
#ifndef RELEASE
    PopCallStack();
#endif
}

template class elemental::DistMatrixBase<int,   MR,MC>;
template class elemental::DistMatrixBase<float, MR,MC>;
template class elemental::DistMatrixBase<double,MR,MC>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,MR,MC>;
template class elemental::DistMatrixBase<dcomplex,MR,MC>;
#endif

