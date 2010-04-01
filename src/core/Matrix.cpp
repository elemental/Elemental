/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "ElementalEnvironment.h"
#include "ElementalMatrix.h"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::Matrix<T>::Print( const string msg ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::Print");
#endif
    if( msg != "" )
        cout << msg << endl;

    const int height = Height();
    const int width = Width();

    for( int i=0; i<height; ++i )
    {
        for( int j=0; j<width; ++j )
        {
            cout << operator()(i,j) << " ";
        }
        cout << endl;
    }
    cout << endl;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::ResizeTo
( const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::ResizeTo(height,width)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
    if( _viewing )
        throw "Does not make sense to resize matrix when viewing other data.";
#endif
    const int minLDim = 1;
    _height = height;
    _width  = width;
    _ldim   = max( height, minLDim );

    _memory.Require(_ldim*width);
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::ResizeTo
( const int height, const int width, const int ldim )
{
#ifndef RELEASE
    PushCallStack("Matrix::ResizeTo(height,width,ldim)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
    if( _viewing )
        throw "Does not make sense to resize matrix when viewing other data.";
    if( ldim < height )
    {
        ostringstream msg;
        msg << "Tried to set ldim(" << ldim << ") < height (" << height << ")"
            << endl;
        throw msg.str();
    }
#endif
    _height = height;
    _width  = width;
    _ldim   = ldim;

    _memory.Require(ldim*width);
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::View( Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(A)");
    if( _memory.Size() > 0 )
        throw "Viewing with Matrix after allocating memory.";
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
Elemental::Matrix<T>::LockedView( const Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(A)");
    if( _memory.Size() > 0 )
        throw "Viewing with Matrix after allocating memory.";
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
Elemental::Matrix<T>::View
( Matrix<T>& A, 
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(A,i,j,height,width)");
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
    if( _memory.Size() > 0 )
        throw "Viewing with Matrix after allocating memory.";
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        ostringstream msg;
        msg << "Trying to view outside of a Matrix: " 
            << "up to (" << i+height-1 << "," << j+width-1 << ") "
            << "of " << A.Height() << " x " << A.Width() << " Matrix." << endl;
        throw msg.str();
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
Elemental::Matrix<T>::LockedView
( const Matrix<T>& A, 
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(A,i,j,height,width)");
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
    if( _memory.Size() > 0 )
        throw "Viewing with Matrix after allocating memory.";
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        ostringstream msg;
        msg << "Trying to view outside of a Matrix: " 
            << "up to (" << i+height-1 << "," << j+width-1 << ") "
            << "of " << A.Height() << " x " << A.Width() << " Matrix." << endl;
        throw msg.str();
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
Elemental::Matrix<T>::View1x2
( Matrix<T>& AL, Matrix<T>& AR )
{
#ifndef RELEASE
    PushCallStack("Matrix::View1x2");
    if( _memory.Size() > 0 )
        throw "Viewing with Matrix after allocating memory.";
    if( AL.Height() != AR.Height() )
        throw "1x2 must have consistent height to combine.";
    if( AL.LDim() != AR.LDim() )
        throw "1x2 must have consistent ldims to combine.";
    if( AR.Buffer() != (AL.Buffer()+AL.LDim()*AL.Width()) )
        throw "1x2 must have contiguous memory.";
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
Elemental::Matrix<T>::LockedView1x2
( const Matrix<T>& AL, const Matrix<T>& AR )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView1x2");
    if( _memory.Size() > 0 )
        throw "Viewing with Matrix after allocating memory.";
    if( AL.Height() != AR.Height() )
        throw "1x2 must have consistent height to combine.";
    if( AL.LDim() != AR.LDim() )
        throw "1x2 must have consistent ldims to combine.";
    if( AR.LockedBuffer() != (AL.LockedBuffer()+AL.LDim()*AL.Width()) )
        throw "1x2 must have contiguous memory.";
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
Elemental::Matrix<T>::View2x1( Matrix<T>& AT, 
                           Matrix<T>& AB )
{
#ifndef RELEASE
    PushCallStack("Matrix::View2x1");
    if( _memory.Size() > 0 )
        throw "Viewing with matrix after allocating memory.";
    if( AT.Width() != AB.Width() )
        throw "2x1 must have consistent width to combine.";
    if( AT.LDim() != AB.LDim() )
        throw "2x1 must have consistent ldim to combine.";
    if( AB.Buffer() != (AT.Buffer() + AT.Height()) )
        throw "2x1 must have contiguous memory.";
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
Elemental::Matrix<T>::LockedView2x1( const Matrix<T>& AT, 
                                 const Matrix<T>& AB )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView2x1");
    if( _memory.Size() > 0 )
        throw "Viewing with matrix after allocating memory.";
    if( AT.Width() != AB.Width() )
        throw "2x1 must have consistent width to combine.";
    if( AT.LDim() != AB.LDim() )
        throw "2x1 must have consistent ldim to combine.";
    if( AB.LockedBuffer() != (AT.LockedBuffer()+AT.Height()) )
        throw "2x1 must have contiguous memory.";
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
Elemental::Matrix<T>::View2x2( Matrix<T>& ATL, Matrix<T>& ATR,
                           Matrix<T>& ABL, Matrix<T>& ABR )
{
#ifndef RELEASE
    PushCallStack("Matrix::View2x2");
    if( _memory.Size() > 0 )
        throw "Viewing a matrix after allocating memory.";
    if( ATL.Width() != ABL.Width()   ||
        ATR.Width() != ABR.Width()   ||
        ATL.Height() != ATR.Height() ||
        ABL.Height() != ABR.Height()   )
    {
        throw "2x2 must conform to combine.";
    }
    if( ATL.LDim() != ATR.LDim() ||
        ATR.LDim() != ABL.LDim() ||
        ABL.LDim() != ABR.LDim()   )
    {
        throw "2x2 must have consistent ldims to combine.";
    }
    if( ABL.Buffer() != (ATL.Buffer() + ATL.Height()) ||
        ABR.Buffer() != (ATR.Buffer() + ATR.Height()) ||
        ATR.Buffer() != (ATL.Buffer() + ATL.LDim()*ATL.Width()) )
    {
        throw "2x2 must have contiguous memory.";
    }
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
Elemental::Matrix<T>::LockedView2x2
( const Matrix<T>& ATL, const Matrix<T>& ATR,
  const Matrix<T>& ABL, const Matrix<T>& ABR )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView2x2");
    if( _memory.Size() > 0 )
        throw "Viewing a matrix after allocating memory.";
    if( ATL.Width() != ABL.Width()   ||
        ATR.Width() != ABR.Width()   ||
        ATL.Height() != ATR.Height() ||
        ABL.Height() != ABR.Height()   )
    {
        throw "2x2 must conform to combine.";
    }
    if( ATL.LDim() != ATR.LDim() ||
        ATR.LDim() != ABL.LDim() ||
        ABL.LDim() != ABR.LDim()   )
    {
        throw "2x2 must have consistent ldims to combine.";
    }
    if( ABL.LockedBuffer() != (ATL.LockedBuffer() + ATL.Height()) ||
        ABR.LockedBuffer() != (ATR.LockedBuffer() + ATR.Height()) ||
        ATR.LockedBuffer() != (ATL.LockedBuffer() + ATL.LDim()*ATL.Width()) )
    {
        throw "2x2 must have contiguous memory.";
    }
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
Elemental::Matrix<T>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToIdentity");
    if( _lockedView )
        throw "Cannot set a locked view to identity.";
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
Elemental::Matrix<T>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToZero");
    if( _lockedView )
        throw "Cannot set a locked view to zero.";
#endif
    const int height = Height();
    const int width = Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            _data[i+j*_ldim] = (T)0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToRandom");
    if( _lockedView )
        throw "Cannot change the data of a locked view.";
#endif
    const int height = Height();
    const int width = Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            _data[i+j*_ldim] = Random<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const Matrix<T>&
Elemental::Matrix<T>::operator=
( const Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::operator=");
    if( _lockedView )
        throw "Cannot assign to a locked view.";
    if( _viewing && ( A.Height() != Height() || A.Width() != Width() ) )
        throw "Cannot assign to a view of different dimensions.";
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int height = Height();
    const int width = Width();
    const int ldim = LDim();
    const int ldimOfA = A.LDim();
    const T* data = A.LockedBuffer();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            _data[i+j*ldim] = data[i+j*ldimOfA];
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class Elemental::Matrix<int>;
template class Elemental::Matrix<float>;
template class Elemental::Matrix<double>;
#ifndef WITHOUT_COMPLEX
template class Elemental::Matrix<scomplex>;
template class Elemental::Matrix<dcomplex>;
#endif

