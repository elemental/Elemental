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
template<typename T,typename Int=int>
class Matrix
{
public:    
    //
    // Constructors
    // 

    Matrix(); 
    Matrix( Int height, Int width );
    Matrix( Int height, Int width, Int ldim );
    Matrix( Int height, Int width, const T* buffer, Int ldim );
    Matrix( Int height, Int width, T* buffer, Int ldim );
    Matrix( const Matrix<T,Int>& A );

    //
    // Destructor
    //

    ~Matrix();

    //
    // Basic information
    //

    Int Height() const;
    Int Width() const;
    Int DiagonalLength( Int offset=0 ) const;
    Int LDim() const;
    Int MemorySize() const;

    T* Buffer();
    T* Buffer( Int i, Int j );
    T* Buffer( Int i, Int j, Int height, Int width );

    const T* LockedBuffer() const;
    const T* LockedBuffer( Int i, Int j ) const;
    const T* LockedBuffer( Int i, Int j, Int height, Int width ) const;

    //
    // I/O
    //

    void Print( const std::string msg="" ) const;
    void Print( std::ostream& os, const std::string msg="" ) const;

    //
    // Entry manipulation
    //

    T Get( Int i, Int j ) const;
    void Set( Int i, Int j, T alpha );
    void Update( Int i, Int j, T alpha );

    void GetDiagonal( Matrix<T,Int>& d, Int offset=0 ) const;
    void SetDiagonal( const Matrix<T,Int>& d, Int offset=0 );
    void UpdateDiagonal( const Matrix<T,Int>& d, Int offset=0 );

    // Only valid for complex datatypes

    typename RealBase<T>::type GetReal( Int i, Int j ) const;
    typename RealBase<T>::type GetImag( Int i, Int j ) const;
    void SetReal( Int i, Int j, typename RealBase<T>::type alpha );
    void SetImag( Int i, Int j, typename RealBase<T>::type alpha );
    void UpdateReal( Int i, Int j, typename RealBase<T>::type alpha );
    void UpdateImag( Int i, Int j, typename RealBase<T>::type alpha );

    void GetRealDiagonal
    ( Matrix<typename RealBase<T>::type>& d, Int offset=0 ) const;
    void GetImagDiagonal
    ( Matrix<typename RealBase<T>::type>& d, Int offset=0 ) const;
    void SetRealDiagonal
    ( const Matrix<typename RealBase<T>::type>& d, Int offset=0 );
    void SetImagDiagonal
    ( const Matrix<typename RealBase<T>::type>& d, Int offset=0 );
    void UpdateRealDiagonal
    ( const Matrix<typename RealBase<T>::type>& d, Int offset=0 );
    void UpdateImagDiagonal
    ( const Matrix<typename RealBase<T>::type>& d, Int offset=0 );

    //
    // Viewing other matrix instances (or buffers)
    //

    bool Viewing() const;
    bool LockedView() const;

    void View( Int height, Int width, T* buffer, Int ldim );
    void View( Matrix<T,Int>& A);
    void View( Matrix<T,Int>& A, Int i, Int j, Int height, Int width );
    void View1x2( Matrix<T,Int>& AL, Matrix<T,Int>& AR );
    void View2x1( Matrix<T,Int>& AT, 
                  Matrix<T,Int>& AB );
    void View2x2( Matrix<T,Int>& ATL, Matrix<T,Int>& ATR,
                  Matrix<T,Int>& ABL, Matrix<T,Int>& ABR );

    void LockedView( Int height, Int width, const T* buffer, Int ldim );
    void LockedView( const Matrix<T,Int>& A );
    void LockedView
    ( const Matrix<T,Int>& A, Int i, Int j, Int height, Int width );
    void LockedView1x2
    ( const Matrix<T,Int>& AL, const Matrix<T,Int>& AR );
    void LockedView2x1
    ( const Matrix<T,Int>& AT, 
      const Matrix<T,Int>& AB );
    void LockedView2x2
    ( const Matrix<T,Int>& ATL, const Matrix<T,Int>& ATR,
      const Matrix<T,Int>& ABL, const Matrix<T,Int>& ABR );

    //
    // Utilities
    //

    const Matrix<T,Int>& operator=( const Matrix<T,Int>& A );

    void Empty();

    void ResizeTo( Int height, Int width );
    void ResizeTo( Int height, Int width, Int ldim );

    void MakeTrapezoidal( Side side, UpperOrLower uplo, Int offset=0 );
    void ScaleTrapezoid( T alpha, Side side, UpperOrLower uplo, Int offset=0 );

    void SetToIdentity();
    void SetToRandom();
    void SetToZero();

private:
    bool      viewing_;
    bool      lockedView_;
    Int   height_;
    Int   width_;
    T*        data_;
    const T*  lockedData_;
    Int   ldim_;
    Memory<T> memory_;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func( const Matrix<Z>& parent, Int i, Int j );
    };
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func
        ( const Matrix<std::complex<Z> >& parent, Int i, Int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const Matrix<Z>& parent, Int i, Int j );
    };
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func
        ( const Matrix<std::complex<Z> >& parent, Int i, Int j );
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func( Matrix<Z>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func( Matrix<Z>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func( Matrix<Z>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func( Matrix<Z>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( Matrix<std::complex<Z> >& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagHelper;
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// Constructors
//

template<typename T,typename Int>
inline
Matrix<T,Int>::Matrix()
: viewing_(false), lockedView_(false),
  height_(0), width_(0), data_(0), lockedData_(0), ldim_(1),
  memory_()
{ }

template<typename T,typename Int>
inline
Matrix<T,Int>::Matrix( Int height, Int width )
: viewing_(false), lockedView_(false),
  height_(height), width_(width), lockedData_(0), ldim_(std::max(height,1))
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    memory_.Require( ldim_*width );
    data_ = memory_.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline 
Matrix<T,Int>::Matrix
( Int height, Int width, Int ldim )
: viewing_(false), lockedView_(false),
  height_(height), width_(width), lockedData_(0), ldim_(ldim)
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
    memory_.Require( ldim*width );
    data_ = memory_.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
Matrix<T,Int>::Matrix
( Int height, Int width, const T* buffer, Int ldim )
: viewing_(true), lockedView_(true),
  height_(height), width_(width), data_(0), lockedData_(buffer), ldim_(ldim)
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

template<typename T,typename Int>
inline
Matrix<T,Int>::Matrix
( Int height, Int width, T* buffer, Int ldim )
: viewing_(true), lockedView_(false),
  height_(height), width_(width), data_(buffer), lockedData_(0), ldim_(ldim)
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

template<typename T,typename Int>
inline
Matrix<T,Int>::Matrix
( const Matrix<T,Int>& A )
: viewing_(false), lockedView_(false), 
  height_(0), width_(0), data_(0), lockedData_(0), ldim_(1)
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

template<typename T,typename Int>
inline
Matrix<T,Int>::~Matrix()
{ }

//
// Basic information
//

template<typename T,typename Int>
inline Int 
Matrix<T,Int>::Height() const
{ return height_; }

template<typename T,typename Int>
inline Int
Matrix<T,Int>::Width() const
{ return width_; }

template<typename T,typename Int>
inline Int
Matrix<T,Int>::DiagonalLength( Int offset ) const
{ return elemental::DiagonalLength(height_,width_,offset); }

template<typename T,typename Int>
inline Int
Matrix<T,Int>::LDim() const
{ return ldim_; }

template<typename T,typename Int>
inline Int
Matrix<T,Int>::MemorySize() const
{ return memory_.Size(); }

template<typename T,typename Int>
inline T*
Matrix<T,Int>::Buffer()
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( lockedView_ )
        throw std::logic_error
        ("Cannot return non-const buffer of locked Matrix");
    PopCallStack();
#endif
    return data_;
}

template<typename T,typename Int>
inline const T*
Matrix<T,Int>::LockedBuffer() const
{
    if( lockedView_ )
        return lockedData_;
    else
        return data_;
}

template<typename T,typename Int>
inline T*
Matrix<T,Int>::Buffer( Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( lockedView_ )
        throw std::logic_error
        ("Cannot return non-const buffer of locked Matrix");
    PopCallStack();
#endif
    return &data_[i+j*ldim_];
}

template<typename T,typename Int>
inline const T*
Matrix<T,Int>::LockedBuffer( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedBuffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    PopCallStack();
#endif
    if( lockedView_ )
        return &lockedData_[i+j*ldim_];
    else
        return &data_[i+j*ldim_];
}

template<typename T,typename Int>
inline T*
Matrix<T,Int>::Buffer( Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( lockedView_ )
        throw std::logic_error
        ("Cannot return non-const buffer of locked Matrix");
    if( (height>0 && (i+height)>height_) || (width>0 && (j+width)>width_) )
        throw std::logic_error("Requested out-of-bounds buffer of Matrix");
    PopCallStack();
#endif
    return &data_[i+j*ldim_];
}

template<typename T,typename Int>
inline const T*
Matrix<T,Int>::LockedBuffer
( Int i, Int j, Int height, Int width ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedBuffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( (height>0 && (i+height)>height_) || (width>0 && (j+width)>width_) )
        throw std::logic_error("Requested out-of-bounds buffer of Matrix");
    PopCallStack();
#endif
    if( lockedView_ )
        return &lockedData_[i+j*ldim_];
    else
        return &data_[i+j*ldim_];
}

//
// I/O
//

template<typename T,typename Int>
inline void
Matrix<T,Int>::Print( std::ostream& os, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::Print");
#endif
    if( msg != "" )
        os << msg << std::endl;

    const Int height = Height();
    const Int width = Width();

    for( Int i=0; i<height; ++i )
    {
        for( Int j=0; j<width; ++j )
            os << WrapScalar(Get(i,j)) << " ";
        os << std::endl;
    }
    os << std::endl;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::Print( const std::string msg ) const
{ Print( std::cout, msg ); }

//
// Entry manipulation
//

template<typename T,typename Int>
inline T
Matrix<T,Int>::Get( Int i, Int j ) const
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
    if( lockedData_ )
        return lockedData_[i+j*ldim_];
    else
        return data_[i+j*ldim_];
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::Set( Int i, Int j, T alpha ) 
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
    if( lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    data_[i+j*ldim_] = alpha;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::Update( Int i, Int j, T alpha ) 
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
    if( lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    data_[i+j*ldim_] += alpha;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::GetDiagonal( Matrix<T,Int>& d, Int offset ) const
{ 
#ifndef RELEASE
    PushCallStack("Matrix::GetDiagonal");
    if( d.LockedView() )
        throw std::logic_error("d must not be a locked view");
    if( d.Viewing() && 
        (d.Height() != DiagonalLength(offset) || d.Width() != 1 ))
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const Int diagLength = DiagonalLength(offset);
    if( !d.Viewing() )    
        d.ResizeTo( diagLength, 1 );
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            d.Set( j, 0, Get(j,j+offset) );
    else
        for( Int j=0; j<diagLength; ++j )
            d.Set( j, 0, Get(j-offset,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetDiagonal( const Matrix<T,Int>& d, Int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::SetDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const Int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            Set( j, j+offset, d.Get(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            Set( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::UpdateDiagonal( const Matrix<T,Int>& d, Int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::UpdateDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const Int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            Update( j, j+offset, d.Get(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            Update( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline typename RealBase<T>::type
Matrix<T,Int>::GetReal( Int i, Int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
Matrix<T,Int>::GetRealHelper<Z>::Func
( const Matrix<Z>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::GetRealHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T,typename Int>
template<typename Z>
inline Z
Matrix<T,Int>::GetRealHelper<std::complex<Z> >::Func
( const Matrix<std::complex<Z> >& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::GetRealHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent.height_ || j > parent.width_ )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent.height_
            << " x " << parent.width_ << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
    if( parent.lockedData_ )
        return std::real(parent.lockedData_[i+j*parent.ldim_]);
    else
        return std::real(parent.data_[i+j*parent.ldim_]);
}

template<typename T,typename Int>
inline typename RealBase<T>::type
Matrix<T,Int>::GetImag( Int i, Int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
Matrix<T,Int>::GetImagHelper<Z>::Func
( const Matrix<Z>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::GetImagHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T,typename Int>
template<typename Z>
inline Z
Matrix<T,Int>::GetImagHelper<std::complex<Z> >::Func
( const Matrix<std::complex<Z> >& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::GetImagHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent.height_ || j > parent.width_ )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent.height_
            << " x " << parent.width_ << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
    if( parent.lockedData_ )
        return std::imag(parent.lockedData_[i+j*parent.ldim_]);
    else
        return std::imag(parent.data_[i+j*parent.ldim_]);
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetReal
( Int i, Int j, typename RealBase<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::SetRealHelper<Z>::Func
( Matrix<Z>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::SetRealHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::SetRealHelper<std::complex<Z> >::Func
( Matrix<std::complex<Z> >& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::SetRealHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent.height_ || j > parent.width_ )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent.height_
            << " x " << parent.width_ << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const Z beta = std::imag(parent.data_[i+j*parent.ldim_]);
    parent.data_[i+j*parent.ldim_] = std::complex<Z>( alpha, beta );
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetImag
( Int i, Int j, typename RealBase<T>::type alpha ) 
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::SetImagHelper<Z>::Func
( Matrix<Z>& parent, Int i, Int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::SetImagHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::SetImagHelper<std::complex<Z> >::Func
( Matrix<std::complex<Z> >& parent, Int i, Int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::SetImagHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent.height_ || j > parent.width_ )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent.height_
            << " x " << parent.width_ << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const Z beta = std::real(parent.data_[i+j*parent.ldim_]);
    parent.data_[i+j*parent.ldim_] = std::complex<Z>( beta, alpha );
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::UpdateReal
( Int i, Int j, typename RealBase<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::UpdateRealHelper<Z>::Func
( Matrix<Z>& parent, Int i, Int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateRealHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::UpdateRealHelper<std::complex<Z> >::Func
( Matrix<std::complex<Z> >& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateRealHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent.height_ || j > parent.width_ )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent.height_
            << " x " << parent.width_ << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const std::complex<Z> beta = parent.data_[i+j*parent.ldim_];
    parent.data_[i+j*parent.ldim_] = 
        std::complex<Z>( std::real(beta)+alpha, std::imag(beta) );
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::UpdateImag
( Int i, Int j, typename RealBase<T>::type alpha ) 
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::UpdateImagHelper<Z>::Func
( Matrix<Z>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateImagHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::UpdateImagHelper<std::complex<Z> >::Func
( Matrix<std::complex<Z> >& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateImagHelper::Func");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > parent.height_ || j > parent.width_ )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << parent.height_
            << " x " << parent.width_ << " Matrix.";
        throw std::logic_error( msg.str() );
    }
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const std::complex<Z> beta = parent.data_[i+j*parent.ldim_];
    parent.data_[i+j*parent.ldim_] = 
        std::complex<Z>( std::real(beta), std::imag(beta)+alpha );
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::GetRealDiagonal
( Matrix<typename RealBase<T>::type>& d, Int offset ) const
{ 
#ifndef RELEASE
    PushCallStack("Matrix::GetRealDiagonal");
    if( d.LockedView() )
        throw std::logic_error("d must not be a locked view");
    if( d.Viewing() && 
        (d.Height() != DiagonalLength(offset) || d.Width() != 1))
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const Int diagLength = DiagonalLength(offset);
    if( !d.Viewing() )    
        d.ResizeTo( diagLength, 1 );
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            d.Set( j, 0, GetReal(j,j+offset) );
    else
        for( Int j=0; j<diagLength; ++j )
            d.Set( j, 0, GetReal(j-offset,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::GetImagDiagonal
( Matrix<typename RealBase<T>::type>& d, Int offset ) const
{ 
#ifndef RELEASE
    PushCallStack("Matrix::GetImagDiagonal");
    if( d.LockedView() )
        throw std::logic_error("d must not be a locked view");
    if( d.Viewing() && 
        (d.Height() != DiagonalLength(offset) || d.Width() != 1))
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const Int diagLength = DiagonalLength(offset);
    if( !d.Viewing() )    
        d.ResizeTo( diagLength, 1 );
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            d.Set( j, 0, GetImag(j,j+offset) );
    else
        for( Int j=0; j<diagLength; ++j )
            d.Set( j, 0, GetImag(j-offset,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetRealDiagonal
( const Matrix<typename RealBase<T>::type>& d, Int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::SetRealDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const Int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            SetReal( j, j+offset, d.Get(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            SetReal( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetImagDiagonal
( const Matrix<typename RealBase<T>::type>& d, Int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::SetImagDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const Int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            SetImag( j, j+offset, d.Get(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            SetImag( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::UpdateRealDiagonal
( const Matrix<typename RealBase<T>::type>& d, Int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::UpdateRealDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const Int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            UpdateReal( j, j+offset, d.Get(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            UpdateReal( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::UpdateImagDiagonal
( const Matrix<typename RealBase<T>::type>& d, Int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::UpdateImagDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const Int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            UpdateImag( j, j+offset, d.Get(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            UpdateImag( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Viewing other Matrix instances
//

template<typename T,typename Int>
inline bool
Matrix<T,Int>::Viewing() const
{ return viewing_; }

template<typename T,typename Int>
inline bool
Matrix<T,Int>::LockedView() const
{ return lockedView_; }

template<typename T,typename Int>
inline void
Matrix<T,Int>::View
( Int height, Int width, T* buffer, Int ldim )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(buffer)");
    if( memory_.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
#endif
    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewing_ = true;
    lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::View( Matrix<T,Int>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(A)");
    if( memory_.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
#endif
    height_ = A.Height();
    width_  = A.Width();
    ldim_   = A.LDim();
    data_   = A.Buffer();
    viewing_ = true;
    lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::LockedView
( Int height, Int width, const T* buffer, Int ldim )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(buffer)");
    if( memory_.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
#endif
    height_ = height;
    width_ = width;
    ldim_ = ldim;
    lockedData_ = buffer;
    viewing_ = true;
    lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::LockedView( const Matrix<T,Int>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(A)");
    if( memory_.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
#endif
    height_     = A.Height();
    width_      = A.Width();
    ldim_       = A.LDim();
    lockedData_ = A.LockedBuffer();
    viewing_    = true;
    lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::View
( Matrix<T,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(A,i,j,height,width)");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( memory_.Size() > 0 )
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
    height_     = height;
    width_      = width;
    ldim_       = A.LDim();
    data_       = A.Buffer(i,j,height,width);
    viewing_    = true;
    lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::LockedView
( const Matrix<T,Int>& A, Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(A,i,j,height,width)");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( memory_.Size() > 0 )
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
    height_     = height;
    width_      = width;
    ldim_       = A.LDim();
    lockedData_ = A.LockedBuffer(i,j,height,width);
    viewing_    = true;
    lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::View1x2( Matrix<T,Int>& AL, Matrix<T,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("Matrix::View1x2");
    if( memory_.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
    if( AL.Height() != AR.Height() )
        throw std::logic_error("1x2 must have consistent height to combine");
    if( AL.LDim() != AR.LDim() )
        throw std::logic_error("1x2 must have consistent ldims to combine");
    if( AR.Buffer() != (AL.Buffer()+AL.LDim()*AL.Width()) )
        throw std::logic_error("1x2 must have contiguous memory");
#endif
    height_ = AL.Height();
    width_  = AL.Width() + AR.Width();
    ldim_   = AL.LDim();
    data_   = AL.Buffer();
    viewing_ = true;
    lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::LockedView1x2( const Matrix<T,Int>& AL, const Matrix<T,Int>& AR )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView1x2");
    if( memory_.Size() > 0 )
        throw std::logic_error("Viewing with Matrix after allocating memory");
    if( AL.Height() != AR.Height() )
        throw std::logic_error("1x2 must have consistent height to combine");
    if( AL.LDim() != AR.LDim() )
        throw std::logic_error("1x2 must have consistent ldims to combine");
    if( AR.LockedBuffer() != (AL.LockedBuffer()+AL.LDim()*AL.Width()) )
        throw std::logic_error("1x2 must have contiguous memory");
#endif
    height_     = AL.Height();
    width_      = AL.Width() + AR.Width();
    ldim_       = AL.LDim();
    lockedData_ = AL.LockedBuffer();
    viewing_ = true;
    lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::View2x1
( Matrix<T,Int>& AT,
  Matrix<T,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("Matrix::View2x1");
    if( memory_.Size() > 0 )
        throw std::logic_error("Viewing with matrix after allocating memory");
    if( AT.Width() != AB.Width() )
        throw std::logic_error("2x1 must have consistent width to combine");
    if( AT.LDim() != AB.LDim() )
        throw std::logic_error("2x1 must have consistent ldim to combine");
    if( AB.Buffer() != (AT.Buffer() + AT.Height()) )
        throw std::logic_error("2x1 must have contiguous memory");
#endif
    height_ = AT.Height() + AB.Height();
    width_  = AT.Width();
    ldim_   = AT.LDim();
    data_   = AT.Buffer();
    viewing_ = true;
    lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::LockedView2x1
( const Matrix<T,Int>& AT,
  const Matrix<T,Int>& AB )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView2x1");
    if( memory_.Size() > 0 )
        throw std::logic_error("Viewing with matrix after allocating memory");
    if( AT.Width() != AB.Width() )
        throw std::logic_error("2x1 must have consistent width to combine");
    if( AT.LDim() != AB.LDim() )
        throw std::logic_error("2x1 must have consistent ldim to combine");
    if( AB.LockedBuffer() != (AT.LockedBuffer()+AT.Height()) )
        throw std::logic_error("2x1 must have contiguous memory");
#endif
    height_     = AT.Height() + AB.Height();
    width_      = AT.Width();
    ldim_       = AT.LDim();
    lockedData_ = AT.LockedBuffer();
    viewing_ = true;
    lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::View2x2
( Matrix<T,Int>& ATL, Matrix<T,Int>& ATR,
  Matrix<T,Int>& ABL, Matrix<T,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("Matrix::View2x2");
    if( memory_.Size() > 0 )
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
    height_ = ATL.Height() + ABL.Height();
    width_  = ATL.Width() + ATR.Width();
    ldim_   = ATL.LDim();
    data_   = ATL.Buffer();
    viewing_ = true;
    lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::LockedView2x2
( const Matrix<T,Int>& ATL, const Matrix<T,Int>& ATR,
  const Matrix<T,Int>& ABL, const Matrix<T,Int>& ABR )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView2x2");
    if( memory_.Size() > 0 )
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
    height_ = ATL.Height() + ABL.Height();
    width_  = ATL.Width() + ATR.Width();
    ldim_   = ATL.LDim();
    lockedData_ = ATL.LockedBuffer();
    viewing_ = true;
    lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Utilities
//

template<typename T,typename Int>
inline const Matrix<T,Int>&
Matrix<T,Int>::operator=( const Matrix<T,Int>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::operator=");
    if( lockedView_ )
        throw std::logic_error("Cannot assign to a locked view");
    if( viewing_ && ( A.Height() != Height() || A.Width() != Width() ) )
        throw std::logic_error
        ("Cannot assign to a view of different dimensions");
#endif
    if( !viewing_ )
        ResizeTo( A.Height(), A.Width() );

    const Int height = Height();
    const Int width = Width();
    const Int ldim = LDim();
    const Int ldimOfA = A.LDim();
    const T* data = A.LockedBuffer();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
        std::memcpy( &data_[j*ldim], &data[j*ldimOfA], height*sizeof(T) );
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::Empty()
{
    memory_.Empty();
    height_ = 0;
    width_ = 0;
    ldim_ = 1;
    data_ = 0;
    lockedData_ = 0;
    viewing_ = false;
    lockedView_ = false;
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::ResizeTo(height,width)");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( viewing_ && (height>height_ || width>width_) )
        throw std::logic_error("Cannot increase the size of a view");
#endif
    // Only change the ldim when necessary. Simply 'shrink' our view if 
    // possible.
    const Int minLDim = 1;
    if( height > height_ || width > width_ )
        ldim_ = std::max( height, minLDim );

    height_ = height;
    width_ = width;

    memory_.Require(ldim_*width);
    data_ = memory_.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    PushCallStack("Matrix::ResizeTo(height,width,ldim)");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( viewing_ && (height > height_ || width > width_ || ldim != ldim_) )
        throw std::logic_error("Illogical ResizeTo on viewed data");
    if( ldim < height )
    {
        std::ostringstream msg;
        msg << "Tried to set ldim(" << ldim << ") < height (" << height << ")";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    height_ = height;
    width_ = width;
    ldim_ = ldim;

    memory_.Require(ldim*width);
    data_ = memory_.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::MakeTrapezoidal( Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("Matrix::MakeTrapezoidal");
    if( lockedView_ )
        throw std::logic_error("Cannot modify a locked view");
#endif
    const Int height = Height();
    const Int width = Width();
    const Int ldim = LDim();
    T* buffer = Buffer();

    if( uplo == LOWER )
    {
        if( side == LEFT )
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int j=std::max(0,offset+1); j<width; ++j )
            {
                const Int lastZeroRow = j-offset-1;
                const Int numZeroRows = std::min( lastZeroRow+1, height );
                std::memset( &buffer[j*ldim], 0, numZeroRows*sizeof(T) );
            }
        }
        else
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int j=std::max(0,offset-height+width+1); j<width; ++j )
            {
                const Int lastZeroRow = j-offset+height-width-1;
                const Int numZeroRows = std::min( lastZeroRow+1, height );
                std::memset( &buffer[j*ldim], 0, numZeroRows*sizeof(T) );
            }
        }
    }
    else
    {
        if( side == LEFT )
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                const Int firstZeroRow = std::max(j-offset+1,0);
                if( firstZeroRow < height )
                    std::memset
                    ( &buffer[firstZeroRow+j*ldim], 0, 
                      (height-firstZeroRow)*sizeof(T) );
            }
        }
        else
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                const Int firstZeroRow = std::max(j-offset+height-width+1,0);
                if( firstZeroRow < height )
                    std::memset
                    ( &buffer[firstZeroRow+j*ldim], 0,
                      (height-firstZeroRow)*sizeof(T) );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::ScaleTrapezoid
( T alpha, Side side, UpperOrLower uplo, Int offset )
{
#ifndef RELEASE
    PushCallStack("Matrix::ScaleTrapezoid");
    if( lockedView_ )
        throw std::logic_error("Cannot modify a locked view");
#endif
    const Int height = Height();
    const Int width = Width();
    const Int ldim = LDim();
    T* buffer = Buffer();

    if( uplo == UPPER )
    {
        if( side == LEFT )
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int j=std::max(0,offset-1); j<width; ++j )
            {
                const Int numRows = j-offset+1;
                for( Int i=0; i<numRows; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
        else
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int j=std::max(0,offset-height+width-1); j<width; ++j )
            {
                const Int numRows = j-offset+height-width+1;
                for( Int i=0; i<numRows; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
    }
    else
    {
        if( side == LEFT )
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                const Int numZeroRows = std::max(j-offset,0);
                for( Int i=numZeroRows; i<height; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
        else
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                const Int numZeroRows = std::max(j-offset+height-width,0);
                for( Int i=numZeroRows; i<height; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToIdentity");
    if( lockedView_ )
        throw std::logic_error("Cannot set a locked view to identity");
#endif
    const Int height = Height();
    const Int width = Width();

    SetToZero();
    for( Int j=0; j<std::min(height,width); ++j )
        data_[j+j*ldim_] = (T)1;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToRandom");
    if( lockedView_ )
        throw std::logic_error("Cannot change the data of a locked view");
#endif
    const Int height = Height();
    const Int width = Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            data_[i+j*ldim_] = SampleUnitBall<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToZero");
    if( lockedView_ )
        throw std::logic_error("Cannot set a locked view to zero");
#endif
    const Int height = Height();
    const Int width = Width();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
        std::memset( &data_[j*ldim_], 0, height*sizeof(T) );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental

#endif /* ELEMENTAL_MATRIX_HPP */

