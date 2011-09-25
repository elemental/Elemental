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

#include "elemental/core/environment.hpp"

namespace elemental {

// Matrix base for arbitrary rings
template<typename T>
class Matrix
{
public:    
    //
    // Constructors
    // 

    Matrix(); 
    Matrix( int height, int width );
    Matrix( int height, int width, int ldim );
    Matrix( int height, int width, const T* buffer, int ldim );
    Matrix( int height, int width, T* buffer, int ldim );
    Matrix( const Matrix<T>& A );

    //
    // Destructor
    //

    ~Matrix();

    //
    // Basic information
    //

    int Height() const;
    int Width() const;
    int DiagonalLength( int offset=0 ) const;
    int LDim() const;
    int MemorySize() const;

    T* Buffer();
    T* Buffer( int i, int j );
    T* Buffer( int i, int j, int height, int width );

    const T* LockedBuffer() const;
    const T* LockedBuffer( int i, int j ) const;
    const T* LockedBuffer( int i, int j, int height, int width ) const;

    //
    // I/O
    //

    void Print( const std::string msg="" ) const;
    void Print( std::ostream& os, const std::string msg="" ) const;

    //
    // Entry manipulation
    //

    T Get( int i, int j ) const;
    void Set( int i, int j, T alpha );
    void Update( int i, int j, T alpha );

    void GetDiagonal( Matrix<T>& d, int offset=0 ) const;
    void SetDiagonal( const Matrix<T>& d, int offset=0 );
    void UpdateDiagonal( const Matrix<T>& d, int offset=0 );

    // Only valid for complex datatypes

    typename RealBase<T>::type GetReal( int i, int j ) const;
    typename RealBase<T>::type GetImag( int i, int j ) const;
    void SetReal( int i, int j, typename RealBase<T>::type alpha );
    void SetImag( int i, int j, typename RealBase<T>::type alpha );
    void UpdateReal( int i, int j, typename RealBase<T>::type alpha );
    void UpdateImag( int i, int j, typename RealBase<T>::type alpha );

    void GetRealDiagonal
    ( Matrix<typename RealBase<T>::type>& d, int offset=0 ) const;
    void GetImagDiagonal
    ( Matrix<typename RealBase<T>::type>& d, int offset=0 ) const;
    void SetRealDiagonal
    ( const Matrix<typename RealBase<T>::type>& d, int offset=0 );
    void SetImagDiagonal
    ( const Matrix<typename RealBase<T>::type>& d, int offset=0 );
    void UpdateRealDiagonal
    ( const Matrix<typename RealBase<T>::type>& d, int offset=0 );
    void UpdateImagDiagonal
    ( const Matrix<typename RealBase<T>::type>& d, int offset=0 );

    //
    // Viewing other matrix instances (or buffers)
    //

    bool Viewing() const;
    bool LockedView() const;

    void View( int height, int width, T* buffer, int ldim );
    void View( Matrix<T>& A);
    void View( Matrix<T>& A, int i, int j, int height, int width );
    void View1x2( Matrix<T>& AL, Matrix<T>& AR );
    void View2x1( Matrix<T>& AT, 
                  Matrix<T>& AB );
    void View2x2( Matrix<T>& ATL, Matrix<T>& ATR,
                  Matrix<T>& ABL, Matrix<T>& ABR );

    void LockedView( int height, int width, const T* buffer, int ldim );
    void LockedView( const Matrix<T>& A );
    void LockedView
    ( const Matrix<T>& A, int i, int j, int height, int width );
    void LockedView1x2( const Matrix<T>& AL, const Matrix<T>& AR );
    void LockedView2x1( const Matrix<T>& AT, 
                        const Matrix<T>& AB );
    void LockedView2x2( const Matrix<T>& ATL, const Matrix<T>& ATR,
                        const Matrix<T>& ABL, const Matrix<T>& ABR );

    //
    // Utilities
    //

    const Matrix<T>& operator=( const Matrix<T>& A );

    void Empty();

    void ResizeTo( int height, int width );
    void ResizeTo( int height, int width, int ldim );

    void SetToIdentity();
    void SetToRandom();
    void SetToZero();

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
#ifndef WITHOUT_COMPLEX
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func( const Matrix<std::complex<Z> >& parent, int i, int j );
    };
#endif
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const Matrix<Z>& parent, int i, int j );
    };
#ifndef WITHOUT_COMPLEX
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func( const Matrix<std::complex<Z> >& parent, int i, int j );
    };
#endif
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func( Matrix<Z>& parent, int i, int j, Z alpha );
    };
#ifndef WITHOUT_COMPLEX
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha );
    };
#endif
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func( Matrix<Z>& parent, int i, int j, Z alpha );
    };
#ifndef WITHOUT_COMPLEX
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha );
    };
#endif
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func( Matrix<Z>& parent, int i, int j, Z alpha );
    };
#ifndef WITHOUT_COMPLEX
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha );
    };
#endif
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func( Matrix<Z>& parent, int i, int j, Z alpha );
    };
#ifndef WITHOUT_COMPLEX
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, int i, int j, Z alpha );
    };
#endif
    template<typename Z> friend struct UpdateImagHelper;
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// Constructors
//

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
: _viewing(false), _lockedView(false), 
  _height(0), _width(0), _data(0), _lockedData(0), _ldim(1)
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

//
// Destructor
//

template<typename T>
inline
Matrix<T>::~Matrix()
{ }

//
// Basic information
//

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
Matrix<T>::DiagonalLength( int offset ) const
{ return elemental::DiagonalLength(_height,_width,offset); }

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
Matrix<T>::Buffer( int i, int j )
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
Matrix<T>::LockedBuffer( int i, int j ) const
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
Matrix<T>::Buffer( int i, int j, int height, int width )
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
Matrix<T>::LockedBuffer( int i, int j, int height, int width ) const
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

//
// I/O
//

template<typename T>
inline void
Matrix<T>::Print( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::Print");
#endif
    if( msg != "" )
        os << msg << std::endl;

    const int height = Height();
    const int width = Width();

    for( int i=0; i<height; ++i )
    {
        for( int j=0; j<width; ++j )
            os << WrapScalar(Get(i,j)) << " ";
        os << std::endl;
    }
    os << std::endl;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::Print( const std::string msg ) const
{ Print( std::cout, msg ); }

//
// Entry manipulation
//

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
Matrix<T>::Set( int i, int j, T alpha ) 
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
Matrix<T>::Update( int i, int j, T alpha ) 
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
inline void
Matrix<T>::GetDiagonal( Matrix<T>& d, int offset ) const
{ 
#ifndef RELEASE
    PushCallStack("Matrix::GetDiagonal");
    if( d.LockedView() )
        throw std::logic_error("d must not be a locked view");
    if( d.Viewing() && 
        (d.Height() != DiagonalLength(offset) || d.Width() != 1 ))
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const int diagLength = DiagonalLength(offset);
    if( !d.Viewing() )    
        d.ResizeTo( diagLength, 1 );
    if( offset >= 0 )
        for( int j=0; j<diagLength; ++j )
            d.Set( j, 0, Get(j,j+offset) );
    else
        for( int j=0; j<diagLength; ++j )
            d.Set( j, 0, Get(j-offset,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::SetDiagonal( const Matrix<T>& d, int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::SetDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( int j=0; j<diagLength; ++j )
            Set( j, j+offset, d.Get(j,0) );
    else
        for( int j=0; j<diagLength; ++j )
            Set( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::UpdateDiagonal( const Matrix<T>& d, int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::UpdateDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( int j=0; j<diagLength; ++j )
            Update( j, j+offset, d.Get(j,0) );
    else
        for( int j=0; j<diagLength; ++j )
            Update( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
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
    
#ifndef WITHOUT_COMPLEX
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
#endif // WITHOUT_COMPLEX

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
    
#ifndef WITHOUT_COMPLEX
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
#endif // WITHOUT_COMPLEX

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
    
#ifndef WITHOUT_COMPLEX
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
#endif // WITHOUT_COMPLEX

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
    
#ifndef WITHOUT_COMPLEX
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
#endif // WITHOUT_COMPLEX

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
    
#ifndef WITHOUT_COMPLEX
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
#endif // WITHOUT_COMPLEX

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
    
#ifndef WITHOUT_COMPLEX
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
#endif // WITHOUT_COMPLEX

template<typename T>
inline void
Matrix<T>::GetRealDiagonal
( Matrix<typename RealBase<T>::type>& d, int offset ) const
{ 
#ifndef RELEASE
    PushCallStack("Matrix::GetRealDiagonal");
    if( d.LockedView() )
        throw std::logic_error("d must not be a locked view");
    if( d.Viewing() && 
        (d.Height() != DiagonalLength(offset) || d.Width() != 1))
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const int diagLength = DiagonalLength(offset);
    if( !d.Viewing() )    
        d.ResizeTo( diagLength, 1 );
    if( offset >= 0 )
        for( int j=0; j<diagLength; ++j )
            d.Set( j, 0, GetReal(j,j+offset) );
    else
        for( int j=0; j<diagLength; ++j )
            d.Set( j, 0, GetReal(j-offset,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::GetImagDiagonal
( Matrix<typename RealBase<T>::type>& d, int offset ) const
{ 
#ifndef RELEASE
    PushCallStack("Matrix::GetImagDiagonal");
    if( d.LockedView() )
        throw std::logic_error("d must not be a locked view");
    if( d.Viewing() && 
        (d.Height() != DiagonalLength(offset) || d.Width() != 1))
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const int diagLength = DiagonalLength(offset);
    if( !d.Viewing() )    
        d.ResizeTo( diagLength, 1 );
    if( offset >= 0 )
        for( int j=0; j<diagLength; ++j )
            d.Set( j, 0, GetImag(j,j+offset) );
    else
        for( int j=0; j<diagLength; ++j )
            d.Set( j, 0, GetImag(j-offset,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::SetRealDiagonal
( const Matrix<typename RealBase<T>::type>& d, int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::SetRealDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( int j=0; j<diagLength; ++j )
            SetReal( j, j+offset, d.Get(j,0) );
    else
        for( int j=0; j<diagLength; ++j )
            SetReal( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::SetImagDiagonal
( const Matrix<typename RealBase<T>::type>& d, int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::SetImagDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( int j=0; j<diagLength; ++j )
            SetImag( j, j+offset, d.Get(j,0) );
    else
        for( int j=0; j<diagLength; ++j )
            SetImag( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::UpdateRealDiagonal
( const Matrix<typename RealBase<T>::type>& d, int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::UpdateRealDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( int j=0; j<diagLength; ++j )
            UpdateReal( j, j+offset, d.Get(j,0) );
    else
        for( int j=0; j<diagLength; ++j )
            UpdateReal( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::UpdateImagDiagonal
( const Matrix<typename RealBase<T>::type>& d, int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::UpdateImagDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( int j=0; j<diagLength; ++j )
            UpdateImag( j, j+offset, d.Get(j,0) );
    else
        for( int j=0; j<diagLength; ++j )
            UpdateImag( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Viewing other Matrix instances
//

template<typename T>
inline bool
Matrix<T>::Viewing() const
{ return _viewing; }

template<typename T>
inline bool
Matrix<T>::LockedView() const
{ return _lockedView; }

template<typename T>
inline void
Matrix<T>::View( int height, int width, T* buffer, int ldim )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(buffer)");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
#endif
    _height = height;
    _width = width;
    _ldim = ldim;
    _data = buffer;
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::View( Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(A)");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
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
inline void
Matrix<T>::LockedView( int height, int width, const T* buffer, int ldim )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(buffer)");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
#endif
    _height = height;
    _width = width;
    _ldim = ldim;
    _lockedData = buffer;
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::LockedView( const Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(A)");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
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
inline void
Matrix<T>::View
( Matrix<T>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(A,i,j,height,width)");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        std::ostringstream msg;
        msg << "Trying to view outside of a Matrix: "
            << "up to (" << i+height-1 << "," << j+width-1 << ") "
            << "of " << A.Height() << " x " << A.Width() << " Matrix.";
        throw std::logic_error( msg.str().c_str() );
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
inline void
Matrix<T>::LockedView
( const Matrix<T>& A, int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(A,i,j,height,width)");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        std::ostringstream msg;
        msg << "Trying to view outside of a Matrix: "
            << "up to (" << i+height-1 << "," << j+width-1 << ") "
            << "of " << A.Height() << " x " << A.Width() << " Matrix.";
        throw std::logic_error( msg.str().c_str() );
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
inline void
Matrix<T>::View1x2( Matrix<T>& AL, Matrix<T>& AR )
{
#ifndef RELEASE
    PushCallStack("Matrix::View1x2");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
    if( AL.Height() != AR.Height() )
        throw std::logic_error("1x2 must have consistent height to combine");
    if( AL.LDim() != AR.LDim() )
        throw std::logic_error("1x2 must have consistent ldims to combine");
    if( AR.Buffer() != (AL.Buffer()+AL.LDim()*AL.Width()) )
        throw std::logic_error("1x2 must have contiguous memory");
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
inline void
Matrix<T>::LockedView1x2( const Matrix<T>& AL, const Matrix<T>& AR )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView1x2");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
    if( AL.Height() != AR.Height() )
        throw std::logic_error("1x2 must have consistent height to combine");
    if( AL.LDim() != AR.LDim() )
        throw std::logic_error("1x2 must have consistent ldims to combine");
    if( AR.LockedBuffer() != (AL.LockedBuffer()+AL.LDim()*AL.Width()) )
        throw std::logic_error("1x2 must have contiguous memory");
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
inline void
Matrix<T>::View2x1
( Matrix<T>& AT,
  Matrix<T>& AB )
{
#ifndef RELEASE
    PushCallStack("Matrix::View2x1");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing with matrix after allocating memory");
    if( AT.Width() != AB.Width() )
        throw std::logic_error("2x1 must have consistent width to combine");
    if( AT.LDim() != AB.LDim() )
        throw std::logic_error("2x1 must have consistent ldim to combine");
    if( AB.Buffer() != (AT.Buffer() + AT.Height()) )
        throw std::logic_error("2x1 must have contiguous memory");
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
inline void
Matrix<T>::LockedView2x1
( const Matrix<T>& AT,
  const Matrix<T>& AB )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView2x1");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing with matrix after allocating memory");
    if( AT.Width() != AB.Width() )
        throw std::logic_error( "2x1 must have consistent width to combine");
    if( AT.LDim() != AB.LDim() )
        throw std::logic_error("2x1 must have consistent ldim to combine");
    if( AB.LockedBuffer() != (AT.LockedBuffer()+AT.Height()) )
        throw std::logic_error("2x1 must have contiguous memory");
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
inline void
Matrix<T>::View2x2
( Matrix<T>& ATL, Matrix<T>& ATR,
  Matrix<T>& ABL, Matrix<T>& ABR )
{
#ifndef RELEASE
    PushCallStack("Matrix::View2x2");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing a matrix after allocating memory");
    if( ATL.Width() != ABL.Width()   ||
        ATR.Width() != ABR.Width()   ||
        ATL.Height() != ATR.Height() ||
        ABL.Height() != ABR.Height()   )
        throw std::logic_error("2x2 must conform to combine");
    if( ATL.LDim() != ATR.LDim() ||
        ATR.LDim() != ABL.LDim() ||
        ABL.LDim() != ABR.LDim()   )
        throw std::logic_error("2x2 must have consistent ldims to combine");
    if( ABL.Buffer() != (ATL.Buffer() + ATL.Height()) ||
        ABR.Buffer() != (ATR.Buffer() + ATR.Height()) ||
        ATR.Buffer() != (ATL.Buffer() + ATL.LDim()*ATL.Width()) )
        throw std::logic_error("2x2 must have contiguous memory");
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
inline void
Matrix<T>::LockedView2x2
( const Matrix<T>& ATL, const Matrix<T>& ATR,
  const Matrix<T>& ABL, const Matrix<T>& ABR )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView2x2");
    if( _memory.Size() > 0 )
        throw std::logic_error("Viewing a matrix after allocating memory");
    if( ATL.Width() != ABL.Width()   ||
        ATR.Width() != ABR.Width()   ||
        ATL.Height() != ATR.Height() ||
        ABL.Height() != ABR.Height()   )
        throw std::logic_error("2x2 must conform to combine");
    if( ATL.LDim() != ATR.LDim() ||
        ATR.LDim() != ABL.LDim() ||
        ABL.LDim() != ABR.LDim()   )
        throw std::logic_error("2x2 must have consistent ldims to combine");
    if( ABL.LockedBuffer() != (ATL.LockedBuffer() + ATL.Height()) ||
        ABR.LockedBuffer() != (ATR.LockedBuffer() + ATR.Height()) ||
        ATR.LockedBuffer() != (ATL.LockedBuffer() + ATL.LDim()*ATL.Width()) )
        throw std::logic_error("2x2 must have contiguous memory");
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

//
// Utilities
//

template<typename T>
inline const Matrix<T>&
Matrix<T>::operator=( const Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::operator=");
    if( _lockedView )
        throw std::logic_error("Cannot assign to a locked view");
    if( _viewing && ( A.Height() != Height() || A.Width() != Width() ) )
        throw std::logic_error
        ("Cannot assign to a view of different dimensions");
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
        std::memcpy( _dataCol, dataCol, height*sizeof(T) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

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
inline void
Matrix<T>::ResizeTo( int height, int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::ResizeTo(height,width)");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( _viewing && (height>_height || width>_width) )
        throw std::logic_error("Cannot increase the size of a view");
#endif
    // Only change the ldim when necessary. Simply 'shrink' our view if 
    // possible.
    const int minLDim = 1;
    if( height > _height || width > _width )
        _ldim = std::max( height, minLDim );

    _height = height;
    _width = width;

    _memory.Require(_ldim*width);
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::ResizeTo( int height, int width, int ldim )
{
#ifndef RELEASE
    PushCallStack("Matrix::ResizeTo(height,width,ldim)");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( _viewing && (height > _height || width > _width || ldim != _ldim) )
        throw std::logic_error("Illogical ResizeTo on viewed data");
    if( ldim < height )
    {
        std::ostringstream msg;
        msg << "Tried to set ldim(" << ldim << ") < height (" << height << ")";
        throw std::logic_error( msg.str().c_str() );
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
inline void
Matrix<T>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToIdentity");
    if( _lockedView )
        throw std::logic_error("Cannot set a locked view to identity");
#endif
    const int height = Height();
    const int width = Width();

    SetToZero();
    for( int j=0; j<std::min(height,width); ++j )
        _data[j+j*_ldim] = (T)1;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Matrix<T>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToRandom");
    if( _lockedView )
        throw std::logic_error("Cannot change the data of a locked view");
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
inline void
Matrix<T>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToZero");
    if( _lockedView )
        throw std::logic_error("Cannot set a locked view to zero");
#endif
    const int height = Height();
    const int width = Width();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<width; ++j )
    {
        T* _dataCol = &(_data[j*_ldim]);
        std::memset( _dataCol, 0, height*sizeof(T) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental

#endif /* ELEMENTAL_MATRIX_HPP */

