/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/environment.hpp"
#include "elemental/matrix.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::Matrix<T>::Print( string msg ) const
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
            cout << operator()(i,j) << " ";
        cout << endl;
    }
    cout << endl;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::Matrix<T>::ResizeTo
( int height, int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::ResizeTo(height,width)");
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
    if( _viewing )
        throw logic_error
        ( "Does not make sense to resize matrix when viewing other data." );
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
elemental::Matrix<T>::ResizeTo
( int height, int width, int ldim )
{
#ifndef RELEASE
    PushCallStack("Matrix::ResizeTo(height,width,ldim)");
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
    if( _viewing )
        throw logic_error
        ( "Does not make sense to resize matrix when viewing other data." );
    if( ldim < height )
    {
        ostringstream msg;
        msg << "Tried to set ldim(" << ldim << ") < height (" << height << ")"
            << endl;
        throw logic_error( msg.str() );
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
elemental::Matrix<T>::View( Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(A)");
    if( _memory.Size() > 0 )
        throw logic_error( "Viewing with Matrix after allocating memory." );
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
elemental::Matrix<T>::LockedView( const Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(A)");
    if( _memory.Size() > 0 )
        throw logic_error( "Viewing with Matrix after allocating memory." );
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
elemental::Matrix<T>::View
( Matrix<T>& A, 
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(A,i,j,height,width)");
    if( i < 0 || j < 0 )
        throw logic_error( "Indices must be non-negative." );
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
    if( _memory.Size() > 0 )
        throw logic_error( "Viewing with Matrix after allocating memory." );
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        ostringstream msg;
        msg << "Trying to view outside of a Matrix: " 
            << "up to (" << i+height-1 << "," << j+width-1 << ") "
            << "of " << A.Height() << " x " << A.Width() << " Matrix." << endl;
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
elemental::Matrix<T>::LockedView
( const Matrix<T>& A, 
  int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(A,i,j,height,width)");
    if( i < 0 || j < 0 )
        throw logic_error( "Indices must be non-negative." );
    if( height < 0 || width < 0 )
        throw logic_error( "Height and width must be non-negative." );
    if( _memory.Size() > 0 )
        throw logic_error( "Viewing with Matrix after allocating memory." );
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        ostringstream msg;
        msg << "Trying to view outside of a Matrix: " 
            << "up to (" << i+height-1 << "," << j+width-1 << ") "
            << "of " << A.Height() << " x " << A.Width() << " Matrix." << endl;
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
elemental::Matrix<T>::View1x2
( Matrix<T>& AL, Matrix<T>& AR )
{
#ifndef RELEASE
    PushCallStack("Matrix::View1x2");
    if( _memory.Size() > 0 )
        throw logic_error( "Viewing with Matrix after allocating memory." );
    if( AL.Height() != AR.Height() )
        throw logic_error( "1x2 must have consistent height to combine." );
    if( AL.LDim() != AR.LDim() )
        throw logic_error( "1x2 must have consistent ldims to combine." );
    if( AR.Buffer() != (AL.Buffer()+AL.LDim()*AL.Width()) )
        throw logic_error( "1x2 must have contiguous memory." );
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
elemental::Matrix<T>::LockedView1x2
( const Matrix<T>& AL, const Matrix<T>& AR )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView1x2");
    if( _memory.Size() > 0 )
        throw logic_error( "Viewing with Matrix after allocating memory." );
    if( AL.Height() != AR.Height() )
        throw logic_error( "1x2 must have consistent height to combine." );
    if( AL.LDim() != AR.LDim() )
        throw logic_error( "1x2 must have consistent ldims to combine." );
    if( AR.LockedBuffer() != (AL.LockedBuffer()+AL.LDim()*AL.Width()) )
        throw logic_error( "1x2 must have contiguous memory." );
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
elemental::Matrix<T>::View2x1( Matrix<T>& AT, 
                           Matrix<T>& AB )
{
#ifndef RELEASE
    PushCallStack("Matrix::View2x1");
    if( _memory.Size() > 0 )
        throw logic_error( "Viewing with matrix after allocating memory." );
    if( AT.Width() != AB.Width() )
        throw logic_error( "2x1 must have consistent width to combine." );
    if( AT.LDim() != AB.LDim() )
        throw logic_error( "2x1 must have consistent ldim to combine." );
    if( AB.Buffer() != (AT.Buffer() + AT.Height()) )
        throw logic_error( "2x1 must have contiguous memory." );
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
elemental::Matrix<T>::LockedView2x1
( const Matrix<T>& AT, 
  const Matrix<T>& AB )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView2x1");
    if( _memory.Size() > 0 )
        throw logic_error( "Viewing with matrix after allocating memory." );
    if( AT.Width() != AB.Width() )
        throw logic_error( "2x1 must have consistent width to combine." );
    if( AT.LDim() != AB.LDim() )
        throw logic_error( "2x1 must have consistent ldim to combine." );
    if( AB.LockedBuffer() != (AT.LockedBuffer()+AT.Height()) )
        throw logic_error( "2x1 must have contiguous memory." );
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
elemental::Matrix<T>::View2x2( Matrix<T>& ATL, Matrix<T>& ATR,
                           Matrix<T>& ABL, Matrix<T>& ABR )
{
#ifndef RELEASE
    PushCallStack("Matrix::View2x2");
    if( _memory.Size() > 0 )
        throw logic_error( "Viewing a matrix after allocating memory." );
    if( ATL.Width() != ABL.Width()   ||
        ATR.Width() != ABR.Width()   ||
        ATL.Height() != ATR.Height() ||
        ABL.Height() != ABR.Height()   )
        throw logic_error( "2x2 must conform to combine." );
    if( ATL.LDim() != ATR.LDim() ||
        ATR.LDim() != ABL.LDim() ||
        ABL.LDim() != ABR.LDim()   )
        throw logic_error( "2x2 must have consistent ldims to combine." );
    if( ABL.Buffer() != (ATL.Buffer() + ATL.Height()) ||
        ABR.Buffer() != (ATR.Buffer() + ATR.Height()) ||
        ATR.Buffer() != (ATL.Buffer() + ATL.LDim()*ATL.Width()) )
        throw logic_error( "2x2 must have contiguous memory." );
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
elemental::Matrix<T>::LockedView2x2
( const Matrix<T>& ATL, const Matrix<T>& ATR,
  const Matrix<T>& ABL, const Matrix<T>& ABR )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView2x2");
    if( _memory.Size() > 0 )
        throw logic_error( "Viewing a matrix after allocating memory." );
    if( ATL.Width() != ABL.Width()   ||
        ATR.Width() != ABR.Width()   ||
        ATL.Height() != ATR.Height() ||
        ABL.Height() != ABR.Height()   )
        throw logic_error( "2x2 must conform to combine." );
    if( ATL.LDim() != ATR.LDim() ||
        ATR.LDim() != ABL.LDim() ||
        ABL.LDim() != ABR.LDim()   )
        throw logic_error( "2x2 must have consistent ldims to combine." );
    if( ABL.LockedBuffer() != (ATL.LockedBuffer() + ATL.Height()) ||
        ABR.LockedBuffer() != (ATR.LockedBuffer() + ATR.Height()) ||
        ATR.LockedBuffer() != (ATL.LockedBuffer() + ATL.LDim()*ATL.Width()) )
        throw logic_error( "2x2 must have contiguous memory." );
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
elemental::Matrix<T>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToIdentity");
    if( _lockedView )
        throw logic_error( "Cannot set a locked view to identity." );
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
elemental::Matrix<T>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToZero");
    if( _lockedView )
        throw logic_error( "Cannot set a locked view to zero." );
#endif
    const int height = Height();
    const int width = Width();
#ifdef RELEASE
    for( int j=0; j<width; ++j )
    {
        T* _dataCol = &(_data[j*_ldim]);
        memset( _dataCol, 0, height*sizeof(T) );
    }
#else
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            _data[i+j*_ldim] = (T)0;
#endif

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::Matrix<T>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToRandom");
    if( _lockedView )
        throw logic_error( "Cannot change the data of a locked view." );
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
elemental::Matrix<T>::operator=
( const Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::operator=");
    if( _lockedView )
        throw logic_error( "Cannot assign to a locked view." );
    if( _viewing && ( A.Height() != Height() || A.Width() != Width() ) )
        throw logic_error( "Cannot assign to a view of different dimensions." );
#endif
    if( !_viewing )
        ResizeTo( A.Height(), A.Width() );

    const int height = Height();
    const int width = Width();
    const int ldim = LDim();
    const int ldimOfA = A.LDim();
    const T* data = A.LockedBuffer();
#ifdef RELEASE
    for( int j=0; j<width; ++j )
    {
        const T* dataCol = &(data[j*ldimOfA]);
        T* _dataCol = &(_data[j*ldim]);
        memcpy( _dataCol, dataCol, height*sizeof(T) );
    }
#else
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            _data[i+j*ldim] = data[i+j*ldimOfA];
#endif

#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class elemental::Matrix<int>;
template class elemental::Matrix<float>;
template class elemental::Matrix<double>;
#ifndef WITHOUT_COMPLEX
template class elemental::Matrix<scomplex>;
template class elemental::Matrix<dcomplex>;
#endif

