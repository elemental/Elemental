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
#ifndef ELEMENTAL_MATRIX_HPP
#define ELEMENTAL_MATRIX_HPP 1

#include "Elemental/Environment.hpp"

namespace Elemental
{
    template<typename T>
    class Matrix
    {
        bool      _viewing;
        bool      _lockedView;
        int       _height;
        int       _width;
        int       _ldim;
        T*        _data;
        const T*  _lockedData;
        Memory<T> _memory;

    public:    

        Matrix();                

        Matrix( const int height, const int width );

        Matrix( const int height, const int width, const int ldim );

        Matrix( const Matrix<T>& A );

        ~Matrix();
        
        void Print() const;
        void Print(std::string msg) const;

        void SetToRandom();

        T& operator() ( const int i, const int j );
        T  operator() ( const int i, const int j ) const;
        
        int Height() const;
        int Width() const;
        int LDim() const;
        int MemorySize() const;

        T* Buffer();
        T* Buffer
        ( const int i, const int j );
        T* Buffer
        ( const int i, const int j, const int height, const int width );

        T* Pointer();
        T* Pointer( const int i, const int j );

        const T* LockedBuffer() const;
        const T* LockedBuffer
        ( const int i, const int j ) const;
        const T* LockedBuffer
        ( const int i, const int j, const int height, const int width ) const;

        const T* LockedPointer() const;
        const T* LockedPointer( const int i, const int j ) const;

        // Resize the matrix
        void ResizeTo( const int height, const int width );
        void ResizeTo( const int height, const int width, const int ldim );

        void View
        ( Matrix<T>& A);

        void View
        ( Matrix<T>& A, 
          const int i, const int j, const int height, const int width );

        void View1x2( Matrix<T>& AL, Matrix<T>& AR );

        void View2x1( Matrix<T>& AT, 
                      Matrix<T>& AB );

        void View2x2( Matrix<T>& ATL, Matrix<T>& ATR,
                      Matrix<T>& ABL, Matrix<T>& ABR );
        
        void LockedView( const Matrix<T>& A );

        void LockedView
        ( const Matrix<T>& A, 
          const int i, const int j, const int height, const int width );

        void LockedView1x2( const Matrix<T>& AL, const Matrix<T>& AR );

        void LockedView2x1( const Matrix<T>& AT, 
                            const Matrix<T>& AB );

        void LockedView2x2( const Matrix<T>& ATL, const Matrix<T>& ATR,
                            const Matrix<T>& ABL, const Matrix<T>& ABR );

        void SetToIdentity();

        void SetToZero();

        const Matrix<T>& operator=( const Matrix<T>& A );
    };
}

/*----------------------------------------------------------------------------*/

template<typename T>
inline
Elemental::Matrix<T>::Matrix()
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _ldim(0), _data(NULL), _lockedData(NULL),
  _memory()
{ }

template<typename T>
inline
Elemental::Matrix<T>::Matrix( const int height, const int width )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _ldim(std::max(height,1)), 
  _lockedData(NULL)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix(height,width)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _memory.Require( _ldim*width );
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline 
Elemental::Matrix<T>::Matrix
( const int height, const int width, const int ldim )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _ldim(ldim), _lockedData(NULL)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix(height,width,ldim)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
    if( ldim < height )
    {
        std::ostringstream msg;
        msg << "Initialized with ldim(" << ldim << ") < "
            << "height(" << height << ")." << std::endl;
        throw msg.str();
    }
    if( ldim == 0 )
        throw "Leading dimensions cannot be zero (for BLAS compat.).";
#endif
    _memory.Require( ldim*width );
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::Matrix<T>::Matrix
( const Matrix<T>& A )
: _viewing(false), _lockedView(false), _lockedData(false)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix( const Matrix& )");
#endif
    if( &A != this )
        *this = A;
    else
        throw "You just tried to construct a Matrix with itself!";
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::Matrix<T>::~Matrix()
{ }

template<typename T>
inline T& 
Elemental::Matrix<T>::operator()
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::operator");
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
    if( i > Height() || j > Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << Height() 
            << " x " << Width() << " Matrix." << std::endl;
        throw msg.str();
    }
    PopCallStack();
#endif
    return _data[i+j*_ldim];
}

template<typename T>
inline T
Elemental::Matrix<T>::operator()
( const int i, const int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::operator");
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
    if( i > Height() || j > Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << Height()
            << " x " << Width() << " Matrix." << std::endl;
        throw msg.str();
    }
    PopCallStack();
#endif
    if( _lockedData )
        return _lockedData[i+j*_ldim];
    else
        return _data[i+j*_ldim];
}

template<typename T>
inline int
Elemental::Matrix<T>::Height() const
{ return _height; }

template<typename T>
inline int
Elemental::Matrix<T>::Width() const
{ return _width; }

template<typename T>
inline int
Elemental::Matrix<T>::LDim() const
{ return _ldim; }

template<typename T>
inline int
Elemental::Matrix<T>::MemorySize() const
{ return _memory.Size(); }

template<typename T>
inline T*
Elemental::Matrix<T>::Buffer()
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( _lockedView )
        throw "Cannot return non-const buffer of locked Matrix.";
    PopCallStack();
#endif
    return _data;
}

template<typename T>
inline const T*
Elemental::Matrix<T>::LockedBuffer() const
{
    if( _lockedView )
        return _lockedData;
    else
        return _data;
}

template<typename T>
inline T*
Elemental::Matrix<T>::Buffer
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
    if( _lockedView )
        throw "Cannot return non-const buffer of locked Matrix.";
    // The height or width of the buffer could be zero, so we 
    // use strict inequalities for flexibility. Pointer() does not.
    if( i>_height || j>_width )
        throw "Requested out-of-bounds buffer of Matrix.";
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
Elemental::Matrix<T>::LockedBuffer
( const int i, const int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedBuffer");
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
    // The height or width of the buffer could be zero, so we 
    // use strict inequalities for flexibility. LockedPointer() does not.
    if( i>_height || j>_width )
        throw "Requested out-of-bounds buffer of Matrix.";
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

template<typename T>
inline T*
Elemental::Matrix<T>::Buffer
( const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
    if( _lockedView )
        throw "Cannot return non-const buffer of locked Matrix.";
    if( (i+height)>_height || (j+width)>_width )
        throw "Requested out-of-bounds buffer of Matrix.";
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
Elemental::Matrix<T>::LockedBuffer
( const int i, const int j, const int height, const int width ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedBuffer");
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
    if( (i+height)>_height || (j+width)>_width )
        throw "Requested out-of-bounds buffer of Matrix.";
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

template<typename T>
inline T*
Elemental::Matrix<T>::Pointer()
{
#ifndef RELEASE
    PushCallStack("Matrix::Pointer");
    if( _lockedView )
        throw "Cannot return non-const pointer to locked Matrix.";
    if( _height == 0 || _width == 0 )
        throw "Requested pointer to empty Matrix.";
    PopCallStack();
#endif
    return _data;
}

template<typename T>
inline const T*
Elemental::Matrix<T>::LockedPointer() const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedPointer");
    if( _height == 0 || _width == 0 )
        throw "Requested locked pointer to empty Matrix.";
    PopCallStack();
#endif
    if( _lockedView )
        return _lockedData;
    else
        return _data;
}

template<typename T>
inline T*
Elemental::Matrix<T>::Pointer
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::Pointer");
    if( i < 0 || j < 0 )
        throw "indices must be non-negative.";
    if( _lockedView )
        throw "Cannot return non-const pointer to locked Matrix.";
    if( i >= _height || j >= _width )
        throw "Requested out-of-bounds pointer to Matrix.";
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
Elemental::Matrix<T>::LockedPointer
( const int i, const int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedPointer");
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
    if( i >= _height || j >= _width )
        throw "Requested out-of-bounds locked pointer to Matrix.";
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

#endif /* ELEMENTAL_MATRIX_HPP */

