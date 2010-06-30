/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#ifndef ELEMENTAL_MATRIX_HPP
#define ELEMENTAL_MATRIX_HPP 1

#include "elemental/environment.hpp"

namespace elemental {

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

    Matrix( int height, int width );

    Matrix( int height, int width, int ldim );

    Matrix( const Matrix<T>& A );

    ~Matrix();
    
    void Print() const;
    void Print(std::string msg) const;

    void SetToRandom();

          T& operator() ( int i, int j );
    const T& operator() ( int i, int j ) const;
    
    int Height() const;
    int Width() const;
    int LDim() const;
    int MemorySize() const;

    T* Buffer();
    T* Buffer( int i, int j );
    T* Buffer( int i, int j, int height, int width );

    T* Pointer();
    T* Pointer( int i, int j );

    const T* LockedBuffer() const;
    const T* LockedBuffer( int i, int j ) const;
    const T* LockedBuffer( int i, int j, int height, int width ) const;

    const T* LockedPointer() const;
    const T* LockedPointer( int i, int j ) const;

    // Resize the matrix
    void ResizeTo( int height, int width );
    void ResizeTo( int height, int width, int ldim );

    void View( Matrix<T>& A);

    void View( Matrix<T>& A, int i, int j, int height, int width );

    void View1x2( Matrix<T>& AL, Matrix<T>& AR );

    void View2x1( Matrix<T>& AT, 
                  Matrix<T>& AB );

    void View2x2( Matrix<T>& ATL, Matrix<T>& ATR,
                  Matrix<T>& ABL, Matrix<T>& ABR );
        
    void LockedView( const Matrix<T>& A );

    void LockedView
    ( const Matrix<T>& A, int i, int j, int height, int width );

    void LockedView1x2( const Matrix<T>& AL, const Matrix<T>& AR );

    void LockedView2x1( const Matrix<T>& AT, 
                        const Matrix<T>& AB );

    void LockedView2x2( const Matrix<T>& ATL, const Matrix<T>& ATR,
                        const Matrix<T>& ABL, const Matrix<T>& ABR );

    void SetToIdentity();

    void SetToZero();

    T operator=( T alpha );
    const Matrix<T>& operator=( const Matrix<T>& A );
};

} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
inline
elemental::Matrix<T>::Matrix()
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _ldim(0), _data(NULL), _lockedData(NULL),
  _memory()
{ }

template<typename T>
inline
elemental::Matrix<T>::Matrix( int height, int width )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _ldim(std::max(height,1)), 
  _lockedData(NULL)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix(height,width)");
    if( height < 0 || width < 0 )
        throw std::logic_error( "Height and width must be non-negative." );
#endif
    _memory.Require( _ldim*width );
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline 
elemental::Matrix<T>::Matrix
( int height, int width, int ldim )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _ldim(ldim), _lockedData(NULL)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix(height,width,ldim)");
    if( height < 0 || width < 0 )
        throw std::logic_error( "Height and width must be non-negative." );
    if( ldim < height )
    {
        std::ostringstream msg;
        msg << "Initialized with ldim(" << ldim << ") < "
            << "height(" << height << ")." << std::endl;
        throw std::logic_error( msg.str() );
    }
    if( ldim == 0 )
        throw std::logic_error
        ( "Leading dimensions cannot be zero (for BLAS compat.)." );
#endif
    _memory.Require( ldim*width );
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
elemental::Matrix<T>::Matrix
( const Matrix<T>& A )
: _viewing(false), _lockedView(false), _lockedData(false)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix( const Matrix& )");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "You just tried to construct a Matrix with itself!" );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
elemental::Matrix<T>::~Matrix()
{ }

template<typename T>
inline T& 
elemental::Matrix<T>::operator()
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::operator");
    if( i < 0 || j < 0 )
        throw std::logic_error( "Indices must be non-negative." );
    if( i > Height() || j > Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << Height() 
            << " x " << Width() << " Matrix." << std::endl;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
    return _data[i+j*_ldim];
}

template<typename T>
inline const T&
elemental::Matrix<T>::operator()
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::operator");
    if( i < 0 || j < 0 )
        throw std::logic_error( "Indices must be non-negative." );
    if( i > Height() || j > Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << Height()
            << " x " << Width() << " Matrix." << std::endl;
        throw std::logic_error( msg.str() );
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
elemental::Matrix<T>::Height() const
{ return _height; }

template<typename T>
inline int
elemental::Matrix<T>::Width() const
{ return _width; }

template<typename T>
inline int
elemental::Matrix<T>::LDim() const
{ return _ldim; }

template<typename T>
inline int
elemental::Matrix<T>::MemorySize() const
{ return _memory.Size(); }

template<typename T>
inline T*
elemental::Matrix<T>::Buffer()
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( _lockedView )
        throw std::logic_error
        ( "Cannot return non-const buffer of locked Matrix." );
    PopCallStack();
#endif
    return _data;
}

template<typename T>
inline const T*
elemental::Matrix<T>::LockedBuffer() const
{
    if( _lockedView )
        return _lockedData;
    else
        return _data;
}

template<typename T>
inline T*
elemental::Matrix<T>::Buffer
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( i < 0 || j < 0 )
        throw std::logic_error( "Indices must be non-negative." );
    if( _lockedView )
        throw std::logic_error
        ( "Cannot return non-const buffer of locked Matrix." );
    // The height or width of the buffer could be zero, so we 
    // use strict inequalities for flexibility. Pointer() does not.
    if( i>_height || j>_width )
        throw std::logic_error( "Requested out-of-bounds buffer of Matrix." );
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
elemental::Matrix<T>::LockedBuffer
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedBuffer");
    if( i < 0 || j < 0 )
        throw std::logic_error( "Indices must be non-negative." );
    // The height or width of the buffer could be zero, so we 
    // use strict inequalities for flexibility. LockedPointer() does not.
    if( i>_height || j>_width )
        throw std::logic_error( "Requested out-of-bounds buffer of Matrix." );
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

template<typename T>
inline T*
elemental::Matrix<T>::Buffer
( int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( i < 0 || j < 0 )
        throw std::logic_error( "Indices must be non-negative." );
    if( _lockedView )
        throw std::logic_error
        ( "Cannot return non-const buffer of locked Matrix." );
    if( (i+height)>_height || (j+width)>_width )
        throw std::logic_error( "Requested out-of-bounds buffer of Matrix." );
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
elemental::Matrix<T>::LockedBuffer
( int i, int j, int height, int width ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedBuffer");
    if( i < 0 || j < 0 )
        throw std::logic_error( "Indices must be non-negative." );
    if( height < 0 || width < 0 )
        throw std::logic_error( "Height and width must be non-negative." );
    if( (i+height)>_height || (j+width)>_width )
        throw std::logic_error( "Requested out-of-bounds buffer of Matrix." );
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

template<typename T>
inline T*
elemental::Matrix<T>::Pointer()
{
#ifndef RELEASE
    PushCallStack("Matrix::Pointer");
    if( _lockedView )
        throw std::logic_error
        ( "Cannot return non-const pointer to locked Matrix." );
    if( _height == 0 || _width == 0 )
        throw std::logic_error( "Requested pointer to empty Matrix." );
    PopCallStack();
#endif
    return _data;
}

template<typename T>
inline const T*
elemental::Matrix<T>::LockedPointer() const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedPointer");
    if( _height == 0 || _width == 0 )
        throw std::logic_error( "Requested locked pointer to empty Matrix." );
    PopCallStack();
#endif
    if( _lockedView )
        return _lockedData;
    else
        return _data;
}

template<typename T>
inline T*
elemental::Matrix<T>::Pointer
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::Pointer");
    if( i < 0 || j < 0 )
        throw std::logic_error( "indices must be non-negative." );
    if( _lockedView )
        throw std::logic_error
        ( "Cannot return non-const pointer to locked Matrix." );
    if( i >= _height || j >= _width )
        throw std::logic_error( "Requested out-of-bounds pointer to Matrix." );
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
elemental::Matrix<T>::LockedPointer
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedPointer");
    if( i < 0 || j < 0 )
        throw std::logic_error( "Indices must be non-negative." );
    if( i >= _height || j >= _width )
        throw std::logic_error
        ( "Requested out-of-bounds locked pointer to Matrix." );
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

template<typename T>
inline T
elemental::Matrix<T>::operator=( T alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::operator=");
#endif
    if( this->Height() == 1 && this->Width() == 1 )
        _data[0] = alpha;
    else
        throw std::logic_error("Scalars can only be assigned to 1x1 matrices.");
#ifndef RELEASE
    PopCallStack();
#endif
    return alpha;
}

#endif /* ELEMENTAL_MATRIX_HPP */

