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
class MatrixBase
{
protected:
    bool      _viewing;
    bool      _lockedView;
    int       _height;
    int       _width;
    T*        _data;
    const T*  _lockedData;
    int       _ldim;
    Memory<T> _memory;

    MatrixBase(); 

    MatrixBase( int height, int width );

    MatrixBase( int height, int width, int ldim );

    MatrixBase( int height, int width, const T* buffer, int ldim );

    MatrixBase( int height, int width, T* buffer, int ldim );

    MatrixBase( const MatrixBase<T>& A );

public:    
    virtual ~MatrixBase();
    
    void Print( const std::string msg="" ) const;
    void Print( std::ostream& os, const std::string msg="" ) const;

    void SetToRandom();

    // Return the value of entry (i,j)
    const T Get( int i, int j ) const;
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

    // Empty the contents (frees all memory)
    void Empty();

    void View( MatrixBase<T>& A);

    void View( MatrixBase<T>& A, int i, int j, int height, int width );

    void View1x2( MatrixBase<T>& AL, MatrixBase<T>& AR );

    void View2x1( MatrixBase<T>& AT, 
                  MatrixBase<T>& AB );

    void View2x2( MatrixBase<T>& ATL, MatrixBase<T>& ATR,
                  MatrixBase<T>& ABL, MatrixBase<T>& ABR );
        
    void LockedView( const MatrixBase<T>& A );

    void LockedView
    ( const MatrixBase<T>& A, int i, int j, int height, int width );

    void LockedView1x2( const MatrixBase<T>& AL, const MatrixBase<T>& AR );

    void LockedView2x1( const MatrixBase<T>& AT, 
                        const MatrixBase<T>& AB );

    void LockedView2x2( const MatrixBase<T>& ATL, const MatrixBase<T>& ATR,
                        const MatrixBase<T>& ABL, const MatrixBase<T>& ABR );

    void SetToIdentity();

    void SetToZero();

    const MatrixBase<T>& operator=( const MatrixBase<T>& A );
};

// Partial specialization of Matrix for real rings.
template<typename Z>
class Matrix : public MatrixBase<Z>
{
public:
    Matrix(); 

    Matrix( int height, int width );

    Matrix( int height, int width, int ldim );

    Matrix( int height, int width, const Z* buffer, int ldim );

    Matrix( int height, int width, Z* buffer, int ldim );

    Matrix( const MatrixBase<Z>& A );

    ~Matrix();

    const Matrix<Z>& operator=( const MatrixBase<Z>& A );
};

#ifndef WITHOUT_COMPLEX
// Partial specialization of Matrix for complex rings.
template<typename Z>
class Matrix< std::complex<Z> > : public MatrixBase< std::complex<Z> >
{
public:
    Matrix(); 

    Matrix( int height, int width );

    Matrix( int height, int width, int ldim );

    Matrix( int height, int width, const std::complex<Z>* buffer, int ldim );

    Matrix( int height, int width, std::complex<Z>* buffer, int ldim );

    Matrix( const MatrixBase< std::complex<Z> >& A );

    ~Matrix();

    // Return the real part of entry (i,j)
    const Z GetReal( int i, int j ) const;
    // Return the imag part of entry (i,j)
    const Z GetImag( int i, int j ) const;
    // Set the new value of the real part of entry (i,j)
    void SetReal( int i, int j, Z alpha );
    // Set the new value of the real part of entry (i,j)
    void SetImag( int i, int j, Z alpha );
    // Update the real part of entry (i,j), i.e., real(A(i,j)) += alpha
    void UpdateReal( int i, int j, Z alpha );
    // Update the imag part of entry (i,j), i.e., imag(A(i,j)) += alpha
    void UpdateImag( int i, int j, Z alpha );

    const Matrix< std::complex<Z> >& 
    operator=( const MatrixBase< std::complex<Z> >& A );
};
#endif // WITHOUT_COMPLEX

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

// MatrixBase
template<typename T>
inline
MatrixBase<T>::MatrixBase()
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _data(0), _lockedData(0), _ldim(0),
  _memory()
{ }

template<typename T>
inline
MatrixBase<T>::MatrixBase( int height, int width )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _lockedData(0), _ldim(std::max(height,1))
{
#ifndef RELEASE
    PushCallStack("MatrixBase::MatrixBase");
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
MatrixBase<T>::MatrixBase
( int height, int width, int ldim )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _lockedData(0), _ldim(ldim)
{
#ifndef RELEASE
    PushCallStack("MatrixBase::MatrixBase");
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
MatrixBase<T>::MatrixBase
( int height, int width, const T* buffer, int ldim )
: _viewing(true), _lockedView(true),
  _height(height), _width(width), _data(0), _lockedData(buffer), _ldim(ldim)
{
#ifndef RELEASE
    PushCallStack("MatrixBase::MatrixBase");
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
MatrixBase<T>::MatrixBase
( int height, int width, T* buffer, int ldim )
: _viewing(true), _lockedView(false),
  _height(height), _width(width), _data(buffer), _lockedData(0), _ldim(ldim)
{
#ifndef RELEASE
    PushCallStack("MatrixBase::MatrixBase");
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
MatrixBase<T>::MatrixBase
( const MatrixBase<T>& A )
: _viewing(false), _lockedView(false), _lockedData(false)
{
#ifndef RELEASE
    PushCallStack("MatrixBase::MatrixBase( const MatrixBase& )");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ("You just tried to construct a MatrixBase with itself!");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
MatrixBase<T>::~MatrixBase()
{ }

template<typename T>
inline void
MatrixBase<T>::Empty()
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
inline const T
MatrixBase<T>::Get
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("MatrixBase::Get");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > Height() || j > Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << Height()
            << " x " << Width() << " MatrixBase.";
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
MatrixBase<T>::Set
( int i, int j, T alpha ) 
{
#ifndef RELEASE
    PushCallStack("MatrixBase::Set");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > Height() || j > Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << Height()
            << " x " << Width() << " MatrixBase.";
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
MatrixBase<T>::Update
( int i, int j, T alpha ) 
{
#ifndef RELEASE
    PushCallStack("MatrixBase::Update");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > Height() || j > Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << Height()
            << " x " << Width() << " MatrixBase.";
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
MatrixBase<T>::Viewing() const
{ return _viewing; }

template<typename T>
inline bool
MatrixBase<T>::LockedView() const
{ return _lockedView; }

template<typename T>
inline int
MatrixBase<T>::Height() const
{ return _height; }

template<typename T>
inline int
MatrixBase<T>::Width() const
{ return _width; }

template<typename T>
inline int
MatrixBase<T>::LDim() const
{ return _ldim; }

template<typename T>
inline int
MatrixBase<T>::MemorySize() const
{ return _memory.Size(); }

template<typename T>
inline T*
MatrixBase<T>::Buffer()
{
#ifndef RELEASE
    PushCallStack("MatrixBase::Buffer");
    if( _lockedView )
        throw std::logic_error
        ("Cannot return non-const buffer of locked MatrixBase");
    PopCallStack();
#endif
    return _data;
}

template<typename T>
inline const T*
MatrixBase<T>::LockedBuffer() const
{
    if( _lockedView )
        return _lockedData;
    else
        return _data;
}

template<typename T>
inline T*
MatrixBase<T>::Buffer
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::Buffer");
    if( i < 0 || j < 0 )
        throw std::logic_error( "Indices must be non-negative." );
    if( _lockedView )
        throw std::logic_error
        ("Cannot return non-const buffer of locked MatrixBase");
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
MatrixBase<T>::LockedBuffer
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("MatrixBase::LockedBuffer");
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
MatrixBase<T>::Buffer
( int i, int j, int height, int width )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::Buffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( _lockedView )
        throw std::logic_error
        ("Cannot return non-const buffer of locked MatrixBase");
    if( (height>0 && (i+height)>_height) || (width>0 && (j+width)>_width) )
        throw std::logic_error("Requested out-of-bounds buffer of MatrixBase");
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
MatrixBase<T>::LockedBuffer
( int i, int j, int height, int width ) const
{
#ifndef RELEASE
    PushCallStack("MatrixBase::LockedBuffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( (height>0 && (i+height)>_height) || (width>0 && (j+width)>_width) )
        throw std::logic_error("Requested out-of-bounds buffer of MatrixBase");
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

template<typename T>
inline T*
MatrixBase<T>::Pointer()
{
#ifndef RELEASE
    PushCallStack("MatrixBase::Pointer");
    if( _lockedView )
        throw std::logic_error
        ("Cannot return non-const pointer to locked MatrixBase");
    if( _height == 0 || _width == 0 )
        throw std::logic_error("Requested pointer to empty MatrixBase");
    PopCallStack();
#endif
    return _data;
}

template<typename T>
inline const T*
MatrixBase<T>::LockedPointer() const
{
#ifndef RELEASE
    PushCallStack("MatrixBase::LockedPointer");
    if( _height == 0 || _width == 0 )
        throw std::logic_error("Requested locked pointer to empty MatrixBase");
    PopCallStack();
#endif
    if( _lockedView )
        return _lockedData;
    else
        return _data;
}

template<typename T>
inline T*
MatrixBase<T>::Pointer
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("MatrixBase::Pointer");
    if( i < 0 || j < 0 )
        throw std::logic_error("indices must be non-negative");
    if( _lockedView )
        throw std::logic_error
        ("Cannot return non-const pointer to locked MatrixBase");
    if( i >= _height || j >= _width )
        throw std::logic_error("Requested out-of-bounds pointer to MatrixBase");
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
MatrixBase<T>::LockedPointer
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("MatrixBase::LockedPointer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i >= _height || j >= _width )
        throw std::logic_error
        ("Requested out-of-bounds locked pointer to MatrixBase");
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

// Matrix for real rings

template<typename Z>
inline
Matrix<Z>::Matrix()
: MatrixBase<Z>()
{ }

template<typename Z>
inline
Matrix<Z>::Matrix
( int height, int width )
: MatrixBase<Z>(height,width)
{ }

template<typename Z>
inline
Matrix<Z>::Matrix
( int height, int width, int ldim )
: MatrixBase<Z>(height,width,ldim)
{ }

template<typename Z>
inline
Matrix<Z>::Matrix
( int height, int width, const Z* buffer, int ldim )
: MatrixBase<Z>(height,width,buffer,ldim)
{ }

template<typename Z>
inline
Matrix<Z>::Matrix
( int height, int width, Z* buffer, int ldim )
: MatrixBase<Z>(height,width,buffer,ldim)
{ }

template<typename Z>
inline
Matrix<Z>::Matrix
( const MatrixBase<Z>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Attempted to construct a Matrix with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
Matrix<Z>::~Matrix()
{ }

template<typename Z>
inline const Matrix<Z>&
Matrix<Z>::operator=
( const MatrixBase<Z>& A )
{ MatrixBase<Z>::operator=( A ); return *this; }

#ifndef WITHOUT_COMPLEX
// Matrix for complex rings

template<typename Z>
inline
Matrix< std::complex<Z> >::Matrix()
: MatrixBase< std::complex<Z> >()
{ }

template<typename Z>
inline
Matrix< std::complex<Z> >::Matrix
( int height, int width )
: MatrixBase< std::complex<Z> >(height,width)
{ }

template<typename Z>
inline
Matrix< std::complex<Z> >::Matrix
( int height, int width, int ldim )
: MatrixBase< std::complex<Z> >(height,width,ldim)
{ }

template<typename Z>
inline
Matrix< std::complex<Z> >::Matrix
( int height, int width, const std::complex<Z>* buffer, int ldim )
: MatrixBase< std::complex<Z> >(height,width,buffer,ldim)
{ }

template<typename Z>
inline
Matrix< std::complex<Z> >::Matrix
( int height, int width, std::complex<Z>* buffer, int ldim )
: MatrixBase< std::complex<Z> >(height,width,buffer,ldim)
{ }

template<typename Z>
inline const Z
Matrix< std::complex<Z> >::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::GetReal");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > this->Height() || j > this->Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << this->Height()
            << " x " << this->Width() << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
    if( this->_lockedData )
        return std::real(this->_lockedData[i+j*this->_ldim]);
    else
        return std::real(this->_data[i+j*this->_ldim]);
}

template<typename Z>
inline const Z
Matrix< std::complex<Z> >::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::GetImag");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > this->Height() || j > this->Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << this->Height()
            << " x " << this->Width() << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
    if( this->_lockedData )
        return std::imag(this->_lockedData[i+j*this->_ldim]);
    else
        return std::imag(this->_data[i+j*this->_ldim]);
}

template<typename Z>
void
Matrix< std::complex<Z> >::SetReal
( int i, int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::SetReal");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > this->Height() || j > this->Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << this->Height()
            << " x " << this->Width() << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( this->_lockedData )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    const Z beta = std::imag(this->_data[i+j*this->_ldim]);
    this->_data[i+j*this->_ldim] = std::complex<Z>( alpha, beta );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline void
Matrix< std::complex<Z> >::SetImag
( int i, int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::SetImag");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > this->Height() || j > this->Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << this->Height()
            << " x " << this->Width() << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( this->_lockedData )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    const Z beta = std::real(this->_data[i+j*this->_ldim]);
    this->_data[i+j*this->_ldim] = std::complex<Z>( beta, alpha );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline void
Matrix< std::complex<Z> >::UpdateReal
( int i, int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateReal");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > this->Height() || j > this->Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << this->Height()
            << " x " << this->Width() << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( this->_lockedData )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    const std::complex<Z> beta = this->_data[i+j*this->_ldim];
    this->_data[i+j*this->_ldim] = 
        std::complex<Z>( std::real(beta)+alpha, std::imag(beta) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline void
Matrix< std::complex<Z> >::UpdateImag
( int i, int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateImag");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > this->Height() || j > this->Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << this->Height()
            << " x " << this->Width() << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( this->_lockedData )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    const std::complex<Z> beta = this->_data[i+j*this->_ldim];
    this->_data[i+j*this->_ldim] = 
        std::complex<Z>( std::real(beta), std::imag(beta)+alpha );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
Matrix< std::complex<Z> >::Matrix
( const MatrixBase< std::complex<Z> >& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Attempted to construct a Matrix with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
Matrix< std::complex<Z> >::~Matrix()
{ }

template<typename Z>
inline const Matrix< std::complex<Z> >&
Matrix< std::complex<Z> >::operator=
( const MatrixBase< std::complex<Z> >& A )
{ MatrixBase< std::complex<Z> >::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // namespace elemental

#endif /* ELEMENTAL_MATRIX_HPP */

