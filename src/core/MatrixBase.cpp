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
#include "elemental/environment.hpp"
#include "elemental/matrix.hpp"
using namespace std;
using namespace elemental;

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
elemental::MatrixBase<T>::Print( ostream& os, const string msg ) const
{
#ifndef RELEASE
    PushCallStack("MatrixBase::Print");
#endif
    if( msg != "" )
        os << msg << endl;

    const int height = Height();
    const int width = Width();

    for( int i=0; i<height; ++i )
    {
        for( int j=0; j<width; ++j )
            os << Get(i,j) << " ";
        os << endl;
    }
    os << endl;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::Print( const string msg ) const
{
    Print( cout, msg );
}

template<typename T>
void
elemental::MatrixBase<T>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::ResizeTo(height,width)");
    if( height < 0 || width < 0 )
        throw logic_error("Height and width must be non-negative");
    if( _viewing && (height>_height || width>_width) )
        throw logic_error("Cannot increase the size of a view");
#endif
    // Only change the ldim when necessary. Simply 'shrink' our view if 
    // possible.
    const int minLDim = 1;
    if( height > _height || width > _width )
        _ldim = max( height, minLDim );

    _height = height;
    _width = width;

    _memory.Require(_ldim*width);
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::ResizeTo
( int height, int width, int ldim )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::ResizeTo(height,width,ldim)");
    if( height < 0 || width < 0 )
        throw logic_error("Height and width must be non-negative");
    if( _viewing && (height > _height || width > _width || ldim != _ldim) )
        throw logic_error("Illogical ResizeTo on viewed data");
    if( ldim < height )
    {
        ostringstream msg;
        msg << "Tried to set ldim(" << ldim << ") < height (" << height << ")";
        throw logic_error( msg.str() );
    }
#endif
    _height = height;
    _width = width;
    _ldim = ldim;

    _memory.Require(ldim*width);
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::View( MatrixBase<T>& A )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::View(A)");
    if( _memory.Size() > 0 )
        throw logic_error("Viewing with MatrixBase after allocating memory");
#endif
    _height = A.Height();
    _width  = A.Width();
    _ldim   = A.LDim();
    _data   = A.Buffer();
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::LockedView( const MatrixBase<T>& A )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::LockedView(A)");
    if( _memory.Size() > 0 )
        throw logic_error("Viewing with MatrixBase after allocating memory");
#endif
    _height     = A.Height();
    _width      = A.Width();
    _ldim       = A.LDim();
    _lockedData = A.LockedBuffer();
    _viewing    = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::View
( MatrixBase<T>& A, 
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::View(A,i,j,height,width)");
    if( i < 0 || j < 0 )
        throw logic_error("Indices must be non-negative");
    if( height < 0 || width < 0 )
        throw logic_error("Height and width must be non-negative");
    if( _memory.Size() > 0 )
        throw logic_error("Viewing with MatrixBase after allocating memory");
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        ostringstream msg;
        msg << "Trying to view outside of a MatrixBase: " 
            << "up to (" << i+height-1 << "," << j+width-1 << ") "
            << "of " << A.Height() << " x " << A.Width() << " MatrixBase.";
        throw logic_error( msg.str() );
    }
#endif
    _height     = height;
    _width      = width;
    _ldim       = A.LDim();
    _data       = A.Buffer(i,j,height,width);
    _viewing    = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::LockedView
( const MatrixBase<T>& A, 
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::LockedView(A,i,j,height,width)");
    if( i < 0 || j < 0 )
        throw logic_error("Indices must be non-negative");
    if( height < 0 || width < 0 )
        throw logic_error("Height and width must be non-negative");
    if( _memory.Size() > 0 )
        throw logic_error("Viewing with MatrixBase after allocating memory");
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        ostringstream msg;
        msg << "Trying to view outside of a MatrixBase: " 
            << "up to (" << i+height-1 << "," << j+width-1 << ") "
            << "of " << A.Height() << " x " << A.Width() << " MatrixBase."; 
        throw logic_error( msg.str() );
    }
#endif
    _height     = height;
    _width      = width;
    _ldim       = A.LDim();
    _lockedData = A.LockedBuffer(i,j,height,width);
    _viewing    = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::View1x2
( MatrixBase<T>& AL, MatrixBase<T>& AR )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::View1x2");
    if( _memory.Size() > 0 )
        throw logic_error("Viewing with MatrixBase after allocating memory");
    if( AL.Height() != AR.Height() )
        throw logic_error("1x2 must have consistent height to combine");
    if( AL.LDim() != AR.LDim() )
        throw logic_error("1x2 must have consistent ldims to combine");
    if( AR.Buffer() != (AL.Buffer()+AL.LDim()*AL.Width()) )
        throw logic_error("1x2 must have contiguous memory");
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _ldim   = AL.LDim();
    _data   = AL.Buffer();
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::LockedView1x2
( const MatrixBase<T>& AL, const MatrixBase<T>& AR )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::LockedView1x2");
    if( _memory.Size() > 0 )
        throw logic_error("Viewing with MatrixBase after allocating memory");
    if( AL.Height() != AR.Height() )
        throw logic_error("1x2 must have consistent height to combine");
    if( AL.LDim() != AR.LDim() )
        throw logic_error("1x2 must have consistent ldims to combine");
    if( AR.LockedBuffer() != (AL.LockedBuffer()+AL.LDim()*AL.Width()) )
        throw logic_error("1x2 must have contiguous memory");
#endif
    _height     = AL.Height();
    _width      = AL.Width() + AR.Width();
    _ldim       = AL.LDim();
    _lockedData = AL.LockedBuffer();
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::View2x1
( MatrixBase<T>& AT, 
  MatrixBase<T>& AB )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::View2x1");
    if( _memory.Size() > 0 )
        throw logic_error("Viewing with matrix after allocating memory");
    if( AT.Width() != AB.Width() )
        throw logic_error("2x1 must have consistent width to combine");
    if( AT.LDim() != AB.LDim() )
        throw logic_error("2x1 must have consistent ldim to combine");
    if( AB.Buffer() != (AT.Buffer() + AT.Height()) )
        throw logic_error("2x1 must have contiguous memory");
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _ldim   = AT.LDim();
    _data   = AT.Buffer();
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::LockedView2x1
( const MatrixBase<T>& AT, 
  const MatrixBase<T>& AB )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::LockedView2x1");
    if( _memory.Size() > 0 )
        throw logic_error("Viewing with matrix after allocating memory");
    if( AT.Width() != AB.Width() )
        throw logic_error( "2x1 must have consistent width to combine");
    if( AT.LDim() != AB.LDim() )
        throw logic_error("2x1 must have consistent ldim to combine");
    if( AB.LockedBuffer() != (AT.LockedBuffer()+AT.Height()) )
        throw logic_error("2x1 must have contiguous memory");
#endif
    _height     = AT.Height() + AB.Height();
    _width      = AT.Width();
    _ldim       = AT.LDim();
    _lockedData = AT.LockedBuffer();
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::View2x2
( MatrixBase<T>& ATL, MatrixBase<T>& ATR,
  MatrixBase<T>& ABL, MatrixBase<T>& ABR )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::View2x2");
    if( _memory.Size() > 0 )
        throw logic_error("Viewing a matrix after allocating memory");
    if( ATL.Width() != ABL.Width()   ||
        ATR.Width() != ABR.Width()   ||
        ATL.Height() != ATR.Height() ||
        ABL.Height() != ABR.Height()   )
        throw logic_error("2x2 must conform to combine");
    if( ATL.LDim() != ATR.LDim() ||
        ATR.LDim() != ABL.LDim() ||
        ABL.LDim() != ABR.LDim()   )
        throw logic_error("2x2 must have consistent ldims to combine");
    if( ABL.Buffer() != (ATL.Buffer() + ATL.Height()) ||
        ABR.Buffer() != (ATR.Buffer() + ATR.Height()) ||
        ATR.Buffer() != (ATL.Buffer() + ATL.LDim()*ATL.Width()) )
        throw logic_error("2x2 must have contiguous memory");
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _ldim   = ATL.LDim();
    _data   = ATL.Buffer();
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::LockedView2x2
( const MatrixBase<T>& ATL, const MatrixBase<T>& ATR,
  const MatrixBase<T>& ABL, const MatrixBase<T>& ABR )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::LockedView2x2");
    if( _memory.Size() > 0 )
        throw logic_error("Viewing a matrix after allocating memory");
    if( ATL.Width() != ABL.Width()   ||
        ATR.Width() != ABR.Width()   ||
        ATL.Height() != ATR.Height() ||
        ABL.Height() != ABR.Height()   )
        throw logic_error("2x2 must conform to combine");
    if( ATL.LDim() != ATR.LDim() ||
        ATR.LDim() != ABL.LDim() ||
        ABL.LDim() != ABR.LDim()   )
        throw logic_error("2x2 must have consistent ldims to combine");
    if( ABL.LockedBuffer() != (ATL.LockedBuffer() + ATL.Height()) ||
        ABR.LockedBuffer() != (ATR.LockedBuffer() + ATR.Height()) ||
        ATR.LockedBuffer() != (ATL.LockedBuffer() + ATL.LDim()*ATL.Width()) )
        throw logic_error("2x2 must have contiguous memory");
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _ldim   = ATL.LDim();
    _lockedData = ATL.LockedBuffer();
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("MatrixBase::SetToIdentity");
    if( _lockedView )
        throw logic_error("Cannot set a locked view to identity");
#endif
    const int height = Height();
    const int width = Width();

    SetToZero();
    for( int j=0; j<min(height,width); ++j )
        _data[j+j*_ldim] = (T)1;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("MatrixBase::SetToZero");
    if( _lockedView )
        throw logic_error("Cannot set a locked view to zero");
#endif
    const int height = Height();
    const int width = Width();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<width; ++j )
    {
        T* _dataCol = &(_data[j*_ldim]);
        memset( _dataCol, 0, height*sizeof(T) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::MatrixBase<T>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("MatrixBase::SetToRandom");
    if( _lockedView )
        throw logic_error("Cannot change the data of a locked view");
#endif
    const int height = Height();
    const int width = Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            _data[i+j*_ldim] = SampleUnitBall<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const MatrixBase<T>&
elemental::MatrixBase<T>::operator=
( const MatrixBase<T>& A )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::operator=");
    if( _lockedView )
        throw logic_error("Cannot assign to a locked view");
    if( _viewing && ( A.Height() != Height() || A.Width() != Width() ) )
        throw logic_error("Cannot assign to a view of different dimensions");
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int height = Height();
    const int width = Width();
    const int ldim = LDim();
    const int ldimOfA = A.LDim();
    const T* data = A.LockedBuffer();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<width; ++j )
    {
        const T* dataCol = &(data[j*ldimOfA]);
        T* _dataCol = &(_data[j*ldim]);
        memcpy( _dataCol, dataCol, height*sizeof(T) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class elemental::MatrixBase<int>;
template class elemental::MatrixBase<float>;
template class elemental::MatrixBase<double>;
#ifndef WITHOUT_COMPLEX
template class elemental::MatrixBase<scomplex>;
template class elemental::MatrixBase<dcomplex>;
#endif

