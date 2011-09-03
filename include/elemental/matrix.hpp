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
#ifndef ELEMENTAL_MATRIX_HPP
#define ELEMENTAL_MATRIX_HPP 1

#include "elemental/environment.hpp"

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

namespace elemental {

// Matrix base for arbitrary rings
template<typename T>
class Matrix
{
public:    
    Matrix(); 
    Matrix( int height, int width );
    Matrix( int height, int width, int ldim );
    Matrix( int height, int width, const T* buffer, int ldim );
    Matrix( int height, int width, T* buffer, int ldim );
    Matrix( const Matrix<T>& A );
    ~Matrix();

    const Matrix<T>& operator=( const Matrix<T>& A );
    
    void Print( const std::string msg="" ) const;
    void Print( std::ostream& os, const std::string msg="" ) const;

    void SetToRandom();

    // Return the value of entry (i,j)
    T Get( int i, int j ) const;
    // Set the new value of entry (i,j)
    void Set( int i, int j, T alpha );
    // Update entry (i,j), i.e., A(i,j) += alpha
    void Update( int i, int j, T alpha );

    bool Viewing() const;
    bool LockedView() const;
    
    int Height() const;
    int Width() const;
    int LDim() const;
    int MemorySize() const;

    T* Buffer();
    T* Buffer( int i, int j );
    T* Buffer( int i, int j, int height, int width );

    const T* LockedBuffer() const;
    const T* LockedBuffer( int i, int j ) const;
    const T* LockedBuffer( int i, int j, int height, int width ) const;

    // Resize the matrix
    void ResizeTo( int height, int width );
    void ResizeTo( int height, int width, int ldim );

    // Empty the contents (frees all memory)
    void Empty();

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
    
    //
    // Only valid for complex datatypes
    //

    // Return the real part of entry (i,j)
    typename RealBase<T>::type GetReal( int i, int j ) const;
    // Return the imag part of entry (i,j)
    typename RealBase<T>::type GetImag( int i, int j ) const;
    // Set the new value of the real part of entry (i,j)
    void SetReal( int i, int j, typename RealBase<T>::type alpha );
    // Set the new value of the real part of entry (i,j)
    void SetImag( int i, int j, typename RealBase<T>::type alpha );
    // Update the real part of entry (i,j), i.e., real(A(i,j)) += alpha
    void UpdateReal( int i, int j, typename RealBase<T>::type alpha );
    // Update the imag part of entry (i,j), i.e., imag(A(i,j)) += alpha
    void UpdateImag( int i, int j, typename RealBase<T>::type alpha );

private:
    bool      _viewing;
    bool      _lockedView;
    int       _height;
    int       _width;
    T*        _data;
    const T*  _lockedData;
    int       _ldim;
    Memory<T> _memory;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func( const Matrix<Z>& parent, int i, int j );
    };
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func( const Matrix<std::complex<Z> >& parent, int i, int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const Matrix<Z>& parent, int i, int j );
    };
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func( const Matrix<std::complex<Z> >& parent, int i, int j );
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func( Matrix<Z>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func( Matrix<Z>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func( Matrix<Z>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func( Matrix<Z>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagHelper;
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
inline
Matrix<T>::Matrix()
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _data(0), _lockedData(0), _ldim(0),
  _memory()
{ }

template<typename T>
inline
Matrix<T>::Matrix( int height, int width )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _lockedData(0), _ldim(std::max(height,1))
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    _memory.Require( _ldim*width );
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline 
Matrix<T>::Matrix
( int height, int width, int ldim )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _lockedData(0), _ldim(ldim)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( ldim < height )
    {
        std::ostringstream msg;
        msg << "Initialized with ldim(" << ldim << ") < "
            << "height(" << height << ").";
        throw std::logic_error( msg.str() );
    }
    if( ldim == 0 )
        throw std::logic_error
        ("Leading dimensions cannot be zero (for BLAS compatibility)");
#endif
    _memory.Require( ldim*width );
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Matrix<T>::Matrix
( int height, int width, const T* buffer, int ldim )
: _viewing(true), _lockedView(true),
  _height(height), _width(width), _data(0), _lockedData(buffer), _ldim(ldim)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( ldim < height )
    {
        std::ostringstream msg;
        msg << "Initialized with ldim(" << ldim << ") < "
            << "height(" << height << ").";
        throw std::logic_error( msg.str() );
    }
    if( ldim == 0 )
        throw std::logic_error
        ( "Leading dimensions cannot be zero (for BLAS compatibility)." );
    PopCallStack();
#endif
}

template<typename T>
inline
Matrix<T>::Matrix
( int height, int width, T* buffer, int ldim )
: _viewing(true), _lockedView(false),
  _height(height), _width(width), _data(buffer), _lockedData(0), _ldim(ldim)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( ldim < height )
    {
        std::ostringstream msg;
        msg << "Initialized with ldim(" << ldim << ") < "
            << "height(" << height << ").";
        throw std::logic_error( msg.str() );
    }
    if( ldim == 0 )
        throw std::logic_error
        ("Leading dimensions cannot be zero (for BLAS compatibility)");
    PopCallStack();
#endif
}

template<typename T>
inline
Matrix<T>::Matrix
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
        ("You just tried to construct a Matrix with itself!");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Matrix<T>::~Matrix()
{ }

template<typename T>
inline void
Matrix<T>::Empty()
{
    _memory.Empty();
    _height = 0;
    _width = 0;
    _ldim = 0;
    _data = 0;
    _lockedData = 0;
    _viewing = false;
    _lockedView = false;
}

template<typename T>
inline T
Matrix<T>::Get( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::Get");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > Height() || j > Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << Height()
            << " x " << Width() << " Matrix.";
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
inline void
Matrix<T>::Set
( int i, int j, T alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::Set");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > Height() || j > Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << Height()
            << " x " << Width() << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( _lockedData )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    _data[i+j*_ldim] = alpha;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::Update
( int i, int j, T alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::Update");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > Height() || j > Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << Height()
            << " x " << Width() << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( _lockedData )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    _data[i+j*_ldim] += alpha;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline bool
Matrix<T>::Viewing() const
{ return _viewing; }

template<typename T>
inline bool
Matrix<T>::LockedView() const
{ return _lockedView; }

template<typename T>
inline int
Matrix<T>::Height() const
{ return _height; }

template<typename T>
inline int
Matrix<T>::Width() const
{ return _width; }

template<typename T>
inline int
Matrix<T>::LDim() const
{ return _ldim; }

template<typename T>
inline int
Matrix<T>::MemorySize() const
{ return _memory.Size(); }

template<typename T>
inline T*
Matrix<T>::Buffer()
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( _lockedView )
        throw std::logic_error
        ("Cannot return non-const buffer of locked Matrix");
    PopCallStack();
#endif
    return _data;
}

template<typename T>
inline const T*
Matrix<T>::LockedBuffer() const
{
    if( _lockedView )
        return _lockedData;
    else
        return _data;
}

template<typename T>
inline T*
Matrix<T>::Buffer
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( i < 0 || j < 0 )
        throw std::logic_error( "Indices must be non-negative." );
    if( _lockedView )
        throw std::logic_error
        ("Cannot return non-const buffer of locked Matrix");
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
Matrix<T>::LockedBuffer
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedBuffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

template<typename T>
inline T*
Matrix<T>::Buffer
( int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( _lockedView )
        throw std::logic_error
        ("Cannot return non-const buffer of locked Matrix");
    if( (height>0 && (i+height)>_height) || (width>0 && (j+width)>_width) )
        throw std::logic_error("Requested out-of-bounds buffer of Matrix");
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
Matrix<T>::LockedBuffer
( int i, int j, int height, int width ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedBuffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( (height>0 && (i+height)>_height) || (width>0 && (j+width)>_width) )
        throw std::logic_error("Requested out-of-bounds buffer of Matrix");
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

template<typename T>
inline typename RealBase<T>::type
Matrix<T>::GetReal( int i, int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
Matrix<T>::GetRealHelper<Z>::Func
( const Matrix<Z>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::GetRealHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T>
template<typename Z>
inline Z
Matrix<T>::GetRealHelper<std::complex<Z> >::Func
( const Matrix<std::complex<Z> >& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::GetRealHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent._height || j > parent._width )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent._height
            << " x " << parent._width << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
    if( parent._lockedData )
        return std::real(parent._lockedData[i+j*parent._ldim]);
    else
        return std::real(parent._data[i+j*parent._ldim]);
}

template<typename T>
inline typename RealBase<T>::type
Matrix<T>::GetImag( int i, int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
Matrix<T>::GetImagHelper<Z>::Func
( const Matrix<Z>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::GetImagHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T>
template<typename Z>
inline Z
Matrix<T>::GetImagHelper<std::complex<Z> >::Func
( const Matrix<std::complex<Z> >& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::GetImagHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent._height || j > parent._width )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent._height
            << " x " << parent._width << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
    if( parent._lockedData )
        return std::imag(parent._lockedData[i+j*parent._ldim]);
    else
        return std::imag(parent._data[i+j*parent._ldim]);
}

template<typename T>
inline void
Matrix<T>::SetReal( int i, int j, typename RealBase<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
Matrix<T>::SetRealHelper<Z>::Func
( Matrix<Z>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::SetRealHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T>
template<typename Z>
inline void
Matrix<T>::SetRealHelper<std::complex<Z> >::Func
( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::SetRealHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent._height || j > parent._width )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent._height
            << " x " << parent._width << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( parent._lockedData )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const Z beta = std::imag(parent._data[i+j*parent._ldim]);
    parent._data[i+j*parent._ldim] = std::complex<Z>( alpha, beta );
}

template<typename T>
inline void
Matrix<T>::SetImag( int i, int j, typename RealBase<T>::type alpha ) 
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
Matrix<T>::SetImagHelper<Z>::Func
( Matrix<Z>& parent, int i, int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::SetImagHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T>
template<typename Z>
inline void
Matrix<T>::SetImagHelper<std::complex<Z> >::Func
( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::SetImagHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent._height || j > parent._width )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent._height
            << " x " << parent._width << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( parent._lockedData )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const Z beta = std::real(parent._data[i+j*parent._ldim]);
    parent._data[i+j*parent._ldim] = std::complex<Z>( beta, alpha );
}

template<typename T>
inline void
Matrix<T>::UpdateReal( int i, int j, typename RealBase<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
Matrix<T>::UpdateRealHelper<Z>::Func
( Matrix<Z>& parent, int i, int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateRealHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T>
template<typename Z>
inline void
Matrix<T>::UpdateRealHelper<std::complex<Z> >::Func
( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateRealHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent._height || j > parent._width )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent._height
            << " x " << parent._width << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( parent._lockedData )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const std::complex<Z> beta = parent._data[i+j*parent._ldim];
    parent._data[i+j*parent._ldim] = 
        std::complex<Z>( std::real(beta)+alpha, std::imag(beta) );
}

template<typename T>
inline void
Matrix<T>::UpdateImag( int i, int j, typename RealBase<T>::type alpha ) 
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
Matrix<T>::UpdateImagHelper<Z>::Func
( Matrix<Z>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateImagHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T>
template<typename Z>
inline void
Matrix<T>::UpdateImagHelper<std::complex<Z> >::Func
( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateImagHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent._height || j > parent._width )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent._height
            << " x " << parent._width << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( parent._lockedData )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const std::complex<Z> beta = parent._data[i+j*parent._ldim];
    parent._data[i+j*parent._ldim] = 
        std::complex<Z>( std::real(beta), std::imag(beta)+alpha );
}

} // namespace elemental

#endif /* ELEMENTAL_MATRIX_HPP */

