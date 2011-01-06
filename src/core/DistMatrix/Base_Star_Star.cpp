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
using namespace elemental::wrappers::mpi;

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
elemental::DistMatrixBase<T,Star,Star>::Print( const string& s ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Print");
#endif
    const Grid& g = this->Grid();
    if( g.VCRank() == 0 && s != "" )
        cout << s << endl;

    const int height = this->Height();
    const int width  = this->Width();

    if( height == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    if( g.InGrid() )
    {
        if( g.VCRank() == 0 )
        {
            for( int i=0; i<height; ++i )
            {
                for( int j=0; j<width; ++j )
                    cout << this->GetLocalEntry(i,j) << " ";
                cout << "\n";
            }
            cout << endl;
        }
        Barrier( g.VCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::View
( DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View");
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_viewing = true;
    this->_lockedView = false;
    if( this->Grid().InGrid() )
        this->_localMatrix.View( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::LockedView
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView");
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
#endif
    this->_height = A.Height();
    this->_width = A.Width();
    this->_viewing = true;
    this->_lockedView = true;
    if( this->Grid().InGrid() )
        this->_localMatrix.LockedView( A.LockedLocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::View
( DistMatrixBase<T,Star,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View");
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    this->_viewing = true;
    this->_lockedView = false;
    if( this->Grid().InGrid() )
        this->_localMatrix.View( A.LocalMatrix(), i, j, height, width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::LockedView
( const DistMatrixBase<T,Star,Star>& A,
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView");
    this->AssertNotStoringData();
    this->AssertSameGrid( A );
    this->AssertValidSubmatrix( A, i, j, height, width );
#endif
    this->_height = height;
    this->_width = width;
    this->_viewing = true;
    this->_lockedView = true;
    if( this->Grid().InGrid() )
    {
        this->_localMatrix.LockedView
        ( A.LockedLocalMatrix(), i, j, height, width );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::View1x2
( DistMatrixBase<T,Star,Star>& AL,
  DistMatrixBase<T,Star,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View1x2");
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_viewing = true;
    this->_lockedView = false;
    if( this->Grid().InGrid() )
        this->_localMatrix.View1x2( AL.LocalMatrix(), AR.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::LockedView1x2
( const DistMatrixBase<T,Star,Star>& AL,
  const DistMatrixBase<T,Star,Star>& AR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView1x2");
    this->AssertNotStoringData();
    this->AssertSameGrid( AL );
    this->AssertSameGrid( AR );
    this->AssertConforming1x2( AL, AR );
#endif
    this->_height = AL.Height();
    this->_width = AL.Width() + AR.Width();
    this->_viewing = true;
    this->_lockedView = true;
    if( this->Grid().InGrid() )
    {
        this->_localMatrix.LockedView1x2
        ( AL.LockedLocalMatrix(), AR.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::View2x1
( DistMatrixBase<T,Star,Star>& AT,
  DistMatrixBase<T,Star,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View2x1");
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_viewing = true;
    this->_lockedView = false;
    if( this->Grid().InGrid() )
    {
        this->_localMatrix.View2x1
        ( AT.LocalMatrix(),
          AB.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::LockedView2x1
( const DistMatrixBase<T,Star,Star>& AT,
  const DistMatrixBase<T,Star,Star>& AB )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView2x1");
    this->AssertNotStoringData();
    this->AssertSameGrid( AT );
    this->AssertSameGrid( AB );
    this->AssertConforming2x1( AT, AB );
#endif
    this->_height = AT.Height() + AB.Height();
    this->_width = AT.Width();
    this->_viewing = true;
    this->_lockedView = true;
    if( this->Grid().InGrid() )
    {
        this->_localMatrix.LockedView2x1
        ( AT.LockedLocalMatrix(),
          AB.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::View2x2
( DistMatrixBase<T,Star,Star>& ATL, 
  DistMatrixBase<T,Star,Star>& ATR,
  DistMatrixBase<T,Star,Star>& ABL,
  DistMatrixBase<T,Star,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::View2x2");
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_viewing = true;
    this->_lockedView = false;
    if( this->Grid().InGrid() )
    {
        this->_localMatrix.View2x2
        ( ATL.LocalMatrix(), ATR.LocalMatrix(),
          ABL.LocalMatrix(), ABR.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::LockedView2x2
( const DistMatrixBase<T,Star,Star>& ATL, 
  const DistMatrixBase<T,Star,Star>& ATR,
  const DistMatrixBase<T,Star,Star>& ABL,
  const DistMatrixBase<T,Star,Star>& ABR )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::LockedView2x2");
    this->AssertNotStoringData();
    this->AssertSameGrid( ATL );
    this->AssertSameGrid( ATR );
    this->AssertSameGrid( ABL );
    this->AssertSameGrid( ABR );
    this->AssertConforming2x2( ATL, ATR, ABL, ABR );
#endif
    this->_height = ATL.Height() + ABL.Height();
    this->_width = ATL.Width() + ATR.Width();
    this->_viewing = true;
    this->_lockedView = true;
    if( this->Grid().InGrid() )
    {
        this->_localMatrix.LockedView2x2
        ( ATL.LockedLocalMatrix(), ATR.LockedLocalMatrix(),
          ABL.LockedLocalMatrix(), ABR.LockedLocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::ResizeTo");
    this->AssertNotLockedView();
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
#endif
    this->_height = height;
    this->_width = width;
    if( this->Grid().InGrid() )
        this->_localMatrix.ResizeTo( height, width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
T
elemental::DistMatrixBase<T,Star,Star>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Get");
    this->AssertValidEntry( i, j );
#endif
    int viewingSize, owningSize;
    MPI_Comm_size( this->Grid().ViewingComm(), &viewingSize );
    MPI_Group_size( this->Grid().OwningGroup(), &owningSize );
    T u;
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        u = this->GetLocalEntry(i,j);
    }
    else
    {
        // Have the root broadcast its data
        if( this->Grid().VCRank() == 0 )
            u = this->GetLocalEntry(i,j);
        Broadcast
        ( &u, 1, this->Grid().OwningToViewingMap(0), 
          this->Grid().ViewingComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::Set
( int i, int j, T u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::Set");
    this->AssertValidEntry( i, j );
#endif
    if( this->Grid().InGrid() )
        this->SetLocalEntry(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utility functions, e.g., SetToIdentity and MakeTrapezoidal
//

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::MakeTrapezoidal
( Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::MakeTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();

    if( this->Grid().InGrid() )
    {
        if( shape == Lower )
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int lastZeroRow = 
                    ( side==Left ? j-offset-1 : j-offset+height-width-1 );
                if( lastZeroRow >= 0 )
                {
                    int boundary = min( lastZeroRow+1, height );
#ifdef RELEASE
                    T* thisCol = this->LocalBuffer(0,j);
                    memset( thisCol, 0, boundary*sizeof(T) );
#else
                    for( int i=0; i<boundary; ++i )
                        this->SetLocalEntry(i,j,0);
#endif
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
#ifdef RELEASE
                if( firstZeroRow < height )
                {
                    T* thisCol = this->LocalBuffer(firstZeroRow,j);
                    memset( thisCol, 0, (height-firstZeroRow)*sizeof(T) );
                }
#else
                for( int i=firstZeroRow; i<height; ++i )
                    this->SetLocalEntry(i,j,0);
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
elemental::DistMatrixBase<T,Star,Star>::ScaleTrapezoidal
( T alpha, Side side, Shape shape, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::ScaleTrapezoidal");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();

    if( this->Grid().InGrid() )
    {
        if( shape == Upper )
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                int lastRow = 
                    ( side==Left ? j-offset : j-offset+height-width );
                int boundary = min( lastRow+1, height );
#ifdef RELEASE
                T* thisCol = this->LocalBuffer(0,j);
                for( int i=0; i<boundary; ++i )
                    thisCol[i] *= alpha;
#else
                for( int i=0; i<boundary; ++i )
                {
                    const T value = this->GetLocalEntry(i,j);
                    this->SetLocalEntry(i,j,alpha*value);
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
                                            : max(j-offset+height-width,0) );
#ifdef RELEASE
                T* thisCol = this->LocalBuffer(firstRow,j);
                for( int i=0; i<(height-firstRow); ++i )
                    thisCol[i] *= alpha;
#else
                for( int i=firstRow; i<height; ++i )
                {
                    const T value = this->GetLocalEntry(i,j);
                    this->SetLocalEntry(i,j,alpha*value);
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
elemental::DistMatrixBase<T,Star,Star>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToIdentity");
    this->AssertNotLockedView();
#endif
    const int height = this->Height();
    const int width = this->Width();

    if( this->Grid().InGrid() )
    {
        this->SetToZero();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<min(height,width); ++j )
            this->SetLocalEntry(j,j,1);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandom");
    this->AssertNotLockedView();
#endif
    // Create random matrix on process 0 and then broadcast
    const Grid& g = this->Grid();
    const int height = this->Height();
    const int width = this->Width();
    const int bufSize = height*width;

    if( g.InGrid() )
    {
        this->_auxMemory.Require( bufSize );

        T* buffer = this->_auxMemory.Buffer();
        if( g.VCRank() == 0 )
        {
            for( int j=0; j<width; ++j )
                for( int i=0; i<height; ++i )
                    buffer[i+j*height] = Random<T>();
        }
        Broadcast( buffer, bufSize, 0, g.VCComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<width; ++j )
        {
            const T* bufferCol = &buffer[j*height];
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, bufferCol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<height; ++i )
                this->SetLocalEntry(i,j,buffer[i+j*height]);
#endif
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MC,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int r = g.Height();
        const int c = g.Width(); 
        const int p = g.Size();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,r);
        const int maxLocalWidth = MaxLocalLength(width,c);

        const int portionSize = 
            max(maxLocalHeight*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (p+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &originalData[j*localHeightOfA];
            memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                originalData[i+j*localHeightOfA] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.VCComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int l=0; l<c; ++l )
        {
            const int rowShift = Shift( l, rowAlignmentOfA, c );
            const int localWidth = LocalLength( width, rowShift, c );

            for( int k=0; k<r; ++k )
            {
                const T* data = &gatheredData[(k+l*r)*portionSize];

                const int colShift = Shift( k, colAlignmentOfA, r );
                const int localHeight = LocalLength( height, colShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int j=0; j<localWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        this->SetLocalEntry
                            (colShift+i*r,rowShift+j*c,data[i+j*localHeight]);
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
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,MC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int r = g.Height();
        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeight = MaxLocalLength(height,r);

        const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<width; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &originalData[j*localHeightOfA];
            memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                originalData[i+j*localHeightOfA] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShift = Shift( k, colAlignmentOfA, r );
            const int localHeight = LocalLength( height, colShift, r );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(colShift+i*r,j,data[i+j*localHeight]);
        }
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,MR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int c = g.Width();
        const int height = this->Height();
        const int width = this->Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength(width,c);

        const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (c+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &originalData[j*height];
            memcpy( originalDataCol, ACol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<height; ++i )
                originalData[i+j*height] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MRComm() );

        // Unpack
        const int rowAlignmentOfA = A.RowAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShift = Shift( k, rowAlignmentOfA, c );
            const int localWidth = LocalLength( width, rowShift, c );

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &data[j*height];
                T* thisCol = this->LocalBuffer(0,rowShift+j*c);
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    this->SetLocalEntry(i,rowShift+j*c,data[i+j*height]);
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
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,MD,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MD,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int lcm = g.LCM();
        const int ownerPath = g.DiagPath( A.ColAlignment() );
        const int ownerPathRank = g.DiagPathRank( A.ColAlignment() );

        const int height = this->Height();
        const int width = this->Width();
        const int localHeight = A.LocalHeight();
        const int maxLocalHeight = MaxLocalLength( height, lcm );
        const int portionSize = max( maxLocalHeight*width, MinCollectContrib );

        // Since a MD communicator has not been implemented, we will take
        // the suboptimal route of 'rounding up' everyone's contribution over 
        // the VC communicator.
        this->_auxMemory.Require( (p+1)*portionSize );
        T* buffer = this->_auxMemory.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        if( A.InDiagonal() )
        {
#ifdef RELEASE
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int j=0; j<width; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,j);
                T* sendBufCol = &sendBuf[j*localHeight];
                memcpy( sendBufCol, ACol, localHeight*sizeof(T) );
            }
#else
# ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<width; ++j )
                for( int i=0; i<localHeight; ++i )
                    sendBuf[i+j*localHeight] = A.GetLocalEntry(i,j);
#endif
        }

        // Communicate
        AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.VCComm() );

        // Unpack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            if( g.DiagPath( k ) == ownerPath )
            {
                const T* data = &recvBuf[k*portionSize];

                const int thisPathRank = g.DiagPathRank( k );
                const int thisColShift = 
                    Shift( thisPathRank, ownerPathRank, lcm );
                const int thisLocalHeight = 
                    LocalLength( height, thisColShift, lcm );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int j=0; j<width; ++j )
                    for( int i=0; i<thisLocalHeight; ++i )
                        this->SetLocalEntry
                            (thisColShift+i*lcm,j,data[i+j*thisLocalHeight]);
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
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,MD>& A )
{ 
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,MD]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int lcm = g.LCM();
        const int ownerPath = g.DiagPath( A.RowAlignment() );
        const int ownerPathRank = g.DiagPathRank( A.RowAlignment() );

        const int height = this->Height();
        const int width = this->Width();
        const int localWidth = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength( width, lcm );
        const int portionSize = max( height*maxLocalWidth, MinCollectContrib );

        // Since a MD communicator has not been implemented, we will take
        // the suboptimal route of 'rounding up' everyone's contribution over 
        // the VC communicator.
        this->_auxMemory.Require( (p+1)*portionSize );
        T* buffer = this->_auxMemory.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[portionSize];

        // Pack
        if( A.InDiagonal() )
        {
#ifdef RELEASE
# ifdef _OPENMP
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* ACol = A.LockedLocalBuffer(0,j);
                T* sendBufCol = &sendBuf[j*height];
                memcpy( sendBufCol, ACol, height*sizeof(T) );
            }
#else
# ifdef _OPENMP
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    sendBuf[i+j*height] = A.GetLocalEntry(i,j);
#endif
        }

        // Communicate
        AllGather
        ( sendBuf, portionSize,
          recvBuf, portionSize, g.VCComm() );

        // Unpack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            if( g.DiagPath( k ) == ownerPath )
            {
                const T* data = &recvBuf[k*portionSize];

                const int thisPathRank = g.DiagPathRank( k );
                const int thisRowShift = 
                    Shift( thisPathRank, ownerPathRank, lcm );
                const int thisLocalWidth = 
                    LocalLength( width, thisRowShift, lcm );

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for
# endif
                for( int j=0; j<thisLocalWidth; ++j )
                {
                    const T* dataCol = &data[j*height];
                    T* thisCol = this->LocalBuffer(0,thisRowShift+j*lcm);
                    memcpy( thisCol, dataCol, height*sizeof(T) );
                }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
# endif
                for( int j=0; j<thisLocalWidth; ++j )
                    for( int i=0; i<height; ++i )
                        this->SetLocalEntry
                            (i,thisRowShift+j*lcm,data[i+j*height]);
#endif
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
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MR,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int p = g.Size();

        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalHeight = MaxLocalLength(height,c);
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = 
            max(maxLocalHeight*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (p+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &originalData[j*localHeightOfA];
            memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                originalData[i+j*localHeightOfA] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.VRComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
        const int rowAlignmentOfA = A.RowAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int l=0; l<r; ++l )
        {
            const int rowShift = Shift( l, rowAlignmentOfA, r );
            const int localWidth = LocalLength( width, rowShift, r );

            for( int k=0; k<c; ++k )
            {
                const T* data = &gatheredData[(k+l*c)*portionSize];

                const int colShift = Shift( k, colAlignmentOfA, c );
                const int localHeight = LocalLength( height, colShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
                #pragma omp parallel for COLLAPSE(2)
#endif
                for( int j=0; j<localWidth; ++j )
                    for( int i=0; i<localHeight; ++i )
                        this->SetLocalEntry
                            (colShift+i*c,rowShift+j*r,data[i+j*localHeight]);
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
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,MR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [MR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int c = g.Width();
        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeight = MaxLocalLength(height,c);

        const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

        this->_auxMemory.Require( (c+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<width; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &originalData[j*localHeightOfA];
            memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                originalData[i+j*localHeightOfA] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MRComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<c; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShift = Shift( k, colAlignmentOfA, c );
            const int localHeight = LocalLength( height, colShift, c );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(colShift+i*c,j,data[i+j*localHeight]);
        }
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,MC]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int r = g.Height();
        const int height = this->Height();
        const int width = this->Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength(width,r);

        const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (r+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &originalData[j*height];
            memcpy( originalDataCol, ACol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<height; ++i )
                originalData[i+j*height] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.MCComm() );

        // Unpack
        const int rowAlignmentOfA = A.RowAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<r; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShift = Shift( k, rowAlignmentOfA, r );
            const int localWidth = LocalLength( width, rowShift, r );

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &data[j*height];
                T* thisCol = this->LocalBuffer(0,rowShift+j*r);
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    this->SetLocalEntry(i,rowShift+j*r,data[i+j*height]);
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
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,VC,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [VC,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

        this->_auxMemory.Require( (p+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<width; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &originalData[j*localHeightOfA];
            memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                originalData[i+j*localHeightOfA] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.VCComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShift = Shift( k, colAlignmentOfA, p );
            const int localHeight = LocalLength( height, colShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(colShift+i*p,j,data[i+j*localHeight]);
        }
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int height = this->Height();
        const int width = this->Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength(width,p);

        const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (p+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &originalData[j*height];
            memcpy( originalDataCol, ACol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<height; ++i )
                originalData[i+j*height] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.VCComm() );

        // Unpack
        const int rowAlignmentOfA = A.RowAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShift = Shift( k, rowAlignmentOfA, p );
            const int localWidth = LocalLength( width, rowShift, p );

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &data[j*height];
                T* thisCol = this->LocalBuffer(0,rowShift+j*p);
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    this->SetLocalEntry(i,rowShift+j*p,data[i+j*height]);
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
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,VR,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [VR,* ]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int height = this->Height();
        const int width = this->Width();
        const int localHeightOfA = A.LocalHeight();
        const int maxLocalHeight = MaxLocalLength(height,p);

        const int portionSize = max(maxLocalHeight*width,MinCollectContrib);

        this->_auxMemory.Require( (p+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<width; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &originalData[j*localHeightOfA];
            memcpy( originalDataCol, ACol, localHeightOfA*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<width; ++j )
            for( int i=0; i<localHeightOfA; ++i )
                originalData[i+j*localHeightOfA] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.VRComm() );

        // Unpack
        const int colAlignmentOfA = A.ColAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int colShift = Shift( k, colAlignmentOfA, p );
            const int localHeight = LocalLength( height, colShift, p );

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
            for( int j=0; j<width; ++j )
                for( int i=0; i<localHeight; ++i )
                    this->SetLocalEntry(colShift+i*p,j,data[i+j*localHeight]);
        }
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,VR]");
    this->AssertNotLockedView();
    this->AssertSameGrid( A );
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );

    const Grid& g = this->Grid();
    if( g.InGrid() )
    {
        const int p = g.Size();
        const int height = this->Height();
        const int width = this->Width();
        const int localWidthOfA = A.LocalWidth();
        const int maxLocalWidth = MaxLocalLength(width,p);

        const int portionSize = max(height*maxLocalWidth,MinCollectContrib);

        this->_auxMemory.Require( (p+1)*portionSize );

        T* buffer = this->_auxMemory.Buffer();
        T* originalData = &buffer[0];
        T* gatheredData = &buffer[portionSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidthOfA; ++j )
        {
            const T* ACol = A.LockedLocalBuffer(0,j);
            T* originalDataCol = &originalData[j*height];
            memcpy( originalDataCol, ACol, height*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidthOfA; ++j )
            for( int i=0; i<height; ++i )
                originalData[i+j*height] = A.GetLocalEntry(i,j);
#endif

        // Communicate
        AllGather
        ( originalData, portionSize,
          gatheredData, portionSize, g.VRComm() );

        // Unpack
        const int rowAlignmentOfA = A.RowAlignment();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int k=0; k<p; ++k )
        {
            const T* data = &gatheredData[k*portionSize];

            const int rowShift = Shift( k, rowAlignmentOfA, p );
            const int localWidth = LocalLength( width, rowShift, p );

#ifdef RELEASE
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for
# endif
            for( int j=0; j<localWidth; ++j )
            {
                const T* dataCol = &data[j*height];
                T* thisCol = this->LocalBuffer(0,rowShift+j*p);
                memcpy( thisCol, dataCol, height*sizeof(T) );
            }
#else
# if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
# endif
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<height; ++i )
                    this->SetLocalEntry(i,rowShift+j*p,data[i+j*height]);
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
const DistMatrixBase<T,Star,Star>&
elemental::DistMatrixBase<T,Star,Star>::operator=
( const DistMatrixBase<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,* ] = [* ,* ]");
    this->AssertNotLockedView();
    if( this->Viewing() )
        this->AssertSameSize( A );
#endif
    if( !this->Viewing() )
        this->ResizeTo( A.Height(), A.Width() );
    if( this->Grid() == A.Grid() )
    {
        this->_localMatrix = A.LockedLocalMatrix();
    }
    else
    {
        if( !CongruentComms( A.Grid().ViewingComm(),
                             this->Grid().ViewingComm() ) )
            throw logic_error
                ("Redistributing between nonmatching grids currently requires"
                 " the viewing communicators to match.");

        // Compute and allocate the amount of required memory
        int requiredMemory = 0;
        if( A.Grid().VCRank() == 0 )
            requiredMemory += A.Height()*A.Width();
        if( this->Grid().InGrid() )
            requiredMemory += A.Height()*A.Width();
        this->_auxMemory.Require( requiredMemory );
        T* buffer = this->_auxMemory.Buffer();
        int offset = 0;
        T* sendBuffer = &buffer[offset];
        if( A.Grid().VCRank() == 0 )
            offset += A.Height()*A.Width();
        T* bcastBuffer = &buffer[offset];

        // Send from the root of A to the root of this matrix's grid
        MPI_Request sendRequest;
        if( A.Grid().VCRank() == 0 )
        {
            for( int j=0; j<A.Width(); ++j ) 
                for( int i=0; i<A.Height(); ++i )
                    sendBuffer[i+j*A.Height()] = A.GetLocalEntry(i,j);
            const int recvViewingRank = this->Grid().OwningToViewingMap(0);
            ISend
            ( sendBuffer, A.Height()*A.Width(), recvViewingRank, 0,
              this->Grid().ViewingComm(), sendRequest );
        }

        // Receive on the root of this matrix's grid and then broadcast
        // over this matrix's owning communicator
        if( this->Grid().InGrid() )
        {
            if( this->Grid().VCRank() == 0 )
            {
                const int sendViewingRank = A.Grid().OwningToViewingMap(0);
                Recv
                ( bcastBuffer, A.Height()*A.Width(), sendViewingRank, 0,
                  this->Grid().ViewingComm() );
            }

            Broadcast
            ( bcastBuffer, A.Height()*A.Width(), 0, this->Grid().VCComm() );

            for( int j=0; j<A.Width(); ++j )
                for( int i=0; i<A.Height(); ++i )
                    this->SetLocalEntry(i,j,bcastBuffer[i+j*A.Height()]);
        }

        if( A.Grid().VCRank() == 0 )
            Wait( sendRequest );
        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::SumOverCol()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SumOverCol");
    this->AssertNotLockedView();
#endif
    if( this->Grid().InGrid() )
    {
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localSize = localHeight*localWidth;
        this->_auxMemory.Require( 2*localSize );
        T* buffer = this->_auxMemory.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[localSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* thisCol = this->LockedLocalBuffer(0,j);
            T* sendBufCol = &sendBuf[j*localHeight];
            memcpy( sendBufCol, thisCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                sendBuf[i+j*localHeight] = this->GetLocalEntry(i,j);
#endif

        // Sum
        AllReduce
        ( sendBuf, recvBuf, localSize, MPI_SUM, this->Grid().MCComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufCol = &recvBuf[j*localHeight];
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,recvBuf[i+j*localHeight]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::SumOverRow()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SumOverRow");
    this->AssertNotLockedView();
#endif
    if( this->Grid().InGrid() )
    {
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localSize = localHeight*localWidth;
        this->_auxMemory.Require( 2*localSize );
        T* buffer = this->_auxMemory.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[localSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* thisCol = this->LockedLocalBuffer(0,j);
            T* sendBufCol = &sendBuf[j*localHeight];
            memcpy( sendBufCol, thisCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                sendBuf[i+j*localHeight] = this->GetLocalEntry(i,j);
#endif

        // Sum
        AllReduce
        ( sendBuf, recvBuf, localSize, MPI_SUM, this->Grid().MRComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufCol = &recvBuf[j*localHeight];
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,recvBuf[i+j*localHeight]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::DistMatrixBase<T,Star,Star>::SumOverGrid()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SumOverGrid");
    this->AssertNotLockedView();
#endif
    if( this->Grid().InGrid() )
    {
        const int localHeight = this->LocalHeight();
        const int localWidth = this->LocalWidth();
        const int localSize = localHeight*localWidth;
        this->_auxMemory.Require( 2*localSize );
        T* buffer = this->_auxMemory.Buffer();
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[localSize];

        // Pack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* thisCol = this->LockedLocalBuffer(0,j);
            T* sendBufCol = &sendBuf[j*localHeight];
            memcpy( sendBufCol, thisCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                sendBuf[i+j*localHeight] = this->GetLocalEntry(i,j);
#endif

        // Sum
        AllReduce
        ( sendBuf, recvBuf, localSize, MPI_SUM, this->Grid().VCComm() );

        // Unpack
#ifdef RELEASE
# ifdef _OPENMP
        #pragma omp parallel for
# endif
        for( int j=0; j<localWidth; ++j )
        {
            const T* recvBufCol = &recvBuf[j*localHeight];
            T* thisCol = this->LocalBuffer(0,j);
            memcpy( thisCol, recvBufCol, localHeight*sizeof(T) );
        }
#else
# ifdef _OPENMP
        #pragma omp parallel for COLLAPSE(2)
# endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<localHeight; ++i )
                this->SetLocalEntry(i,j,recvBuf[i+j*localHeight]);
#endif

        this->_auxMemory.Release();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template class elemental::DistMatrixBase<int,   Star,Star>;
template class elemental::DistMatrixBase<float, Star,Star>;
template class elemental::DistMatrixBase<double,Star,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrixBase<scomplex,Star,Star>;
template class elemental::DistMatrixBase<dcomplex,Star,Star>;
#endif

