/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {

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

    const T* LockedBuffer() const;
    const T* LockedBuffer( Int i, Int j ) const;

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

    //
    // Though the following routines are meant for complex data, all but four
    // logically apply to real data.
    //

    typename Base<T>::type GetRealPart( Int i, Int j ) const;
    typename Base<T>::type GetImagPart( Int i, Int j ) const;
    void SetRealPart( Int i, Int j, typename Base<T>::type alpha );
    // Only valid for complex data
    void SetImagPart( Int i, Int j, typename Base<T>::type alpha );
    void UpdateRealPart( Int i, Int j, typename Base<T>::type alpha );
    // Only valid for complex data
    void UpdateImagPart( Int i, Int j, typename Base<T>::type alpha );

    void GetRealPartOfDiagonal
    ( Matrix<typename Base<T>::type>& d, Int offset=0 ) const;
    void GetImagPartOfDiagonal
    ( Matrix<typename Base<T>::type>& d, Int offset=0 ) const;
    void SetRealPartOfDiagonal
    ( const Matrix<typename Base<T>::type>& d, Int offset=0 );
    // Only valid for complex data
    void SetImagPartOfDiagonal
    ( const Matrix<typename Base<T>::type>& d, Int offset=0 );
    void UpdateRealPartOfDiagonal
    ( const Matrix<typename Base<T>::type>& d, Int offset=0 );
    // Only valid for complex data
    void UpdateImagPartOfDiagonal
    ( const Matrix<typename Base<T>::type>& d, Int offset=0 );

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

private:
    bool viewing_, lockedView_;
    Int height_, width_, ldim_;
    T* data_;
    const T* lockedData_;
    Memory<T> memory_;

    void AssertValidEntry( Int i, Int j ) const;

    template<typename Z>
    struct SetRealPartHelper
    {
        static void Func( Matrix<Z,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetRealPartHelper<Complex<Z> >
    {
        static void Func
        ( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetRealPartHelper;

    template<typename Z>
    struct SetImagPartHelper
    {
        static void Func( Matrix<Z,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetImagPartHelper<Complex<Z> >
    {
        static void Func
        ( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetImagPartHelper;

    template<typename Z>
    struct UpdateRealPartHelper
    {
        static void Func( Matrix<Z,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealPartHelper<Complex<Z> >
    {
        static void Func
        ( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealPartHelper;

    template<typename Z>
    struct UpdateImagPartHelper
    {
        static void Func( Matrix<Z,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagPartHelper<Complex<Z> >
    {
        static void Func
        ( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagPartHelper;
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
  height_(0), width_(0), ldim_(1), data_(0), lockedData_(0),
  memory_()
{ }

template<typename T,typename Int>
inline
Matrix<T,Int>::Matrix( Int height, Int width )
: viewing_(false), lockedView_(false),
  height_(height), width_(width), ldim_(std::max(height,1)), lockedData_(0)
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
  height_(height), width_(width), ldim_(ldim), lockedData_(0)
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
  height_(height), width_(width), ldim_(ldim), data_(0), lockedData_(buffer)
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
  height_(height), width_(width), ldim_(ldim), data_(buffer), lockedData_(0)
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
  height_(0), width_(0), ldim_(1), data_(0), lockedData_(0)
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
{ return elem::DiagonalLength(height_,width_,offset); }

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
            os << Get(i,j) << " ";
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
    AssertValidEntry( i, j );
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
    AssertValidEntry( i, j );
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
    AssertValidEntry( i, j );
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
inline typename Base<T>::type
Matrix<T,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::GetRealPart");
    AssertValidEntry( i, j );
    PopCallStack();
#endif
    if( lockedData_ )
        return RealPart(lockedData_[i+j*ldim_]);
    else
        return RealPart(data_[i+j*ldim_]);
}

template<typename T,typename Int>
inline typename Base<T>::type
Matrix<T,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::GetImagPart");
    AssertValidEntry( i, j );
    PopCallStack();
#endif
    if( lockedData_ )
        return ImagPart(lockedData_[i+j*ldim_]);
    else
        return ImagPart(data_[i+j*ldim_]);
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetRealPart
( Int i, Int j, typename Base<T>::type alpha )
{ SetRealPartHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::SetRealPartHelper<Z>::Func
( Matrix<Z,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::SetRealPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    parent.data_[i+j*parent.ldim_] = alpha;
}
    
template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::SetRealPartHelper<Complex<Z> >::Func
( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::SetRealPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const Z beta = parent.data_[i+j*parent.ldim_].imag;
    parent.data_[i+j*parent.ldim_] = Complex<Z>( alpha, beta );
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetImagPart
( Int i, Int j, typename Base<T>::type alpha ) 
{ SetImagPartHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::SetImagPartHelper<Z>::Func
( Matrix<Z,Int>& parent, Int i, Int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::SetImagPartHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::SetImagPartHelper<Complex<Z> >::Func
( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::SetImagPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const Z beta = parent.data_[i+j*parent.ldim_].real;
    parent.data_[i+j*parent.ldim_] = Complex<Z>( beta, alpha );
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::UpdateRealPart
( Int i, Int j, typename Base<T>::type alpha )
{ UpdateRealPartHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::UpdateRealPartHelper<Z>::Func
( Matrix<Z,Int>& parent, Int i, Int j, Z alpha ) 
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateRealPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    parent.data_[i+j*parent.ldim_] += alpha;
}
    
template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::UpdateRealPartHelper<Complex<Z> >::Func
( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateRealPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const Complex<Z> beta = parent.data_[i+j*parent.ldim_];
    parent.data_[i+j*parent.ldim_] = Complex<Z>( beta.real+alpha, beta.imag );
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::UpdateImagPart
( Int i, Int j, typename Base<T>::type alpha ) 
{ UpdateImagPartHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::UpdateImagPartHelper<Z>::Func
( Matrix<Z,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateImagPartHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T,typename Int>
template<typename Z>
inline void
Matrix<T,Int>::UpdateImagPartHelper<Complex<Z> >::Func
( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("Matrix::UpdateImagPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
    PopCallStack();
#endif
    const Complex<Z> beta = parent.data_[i+j*parent.ldim_];
    parent.data_[i+j*parent.ldim_] = Complex<Z>( beta.real, beta.imag+alpha );
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::GetRealPartOfDiagonal
( Matrix<typename Base<T>::type>& d, Int offset ) const
{ 
#ifndef RELEASE
    PushCallStack("Matrix::GetRealPartOfDiagonal");
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
            d.Set( j, 0, GetRealPart(j,j+offset) );
    else
        for( Int j=0; j<diagLength; ++j )
            d.Set( j, 0, GetRealPart(j-offset,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::GetImagPartOfDiagonal
( Matrix<typename Base<T>::type>& d, Int offset ) const
{ 
#ifndef RELEASE
    PushCallStack("Matrix::GetImagPartOfDiagonal");
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
            d.Set( j, 0, GetImagPart(j,j+offset) );
    else
        for( Int j=0; j<diagLength; ++j )
            d.Set( j, 0, GetImagPart(j-offset,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetRealPartOfDiagonal
( const Matrix<typename Base<T>::type>& d, Int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::SetRealPartOfDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const Int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            SetRealPart( j, j+offset, d.Get(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            SetRealPart( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::SetImagPartOfDiagonal
( const Matrix<typename Base<T>::type>& d, Int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::SetImagPartOfDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Cannot set imaginary part of real matrix");

    const Int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            SetImagPart( j, j+offset, d.Get(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            SetImagPart( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::UpdateRealPartOfDiagonal
( const Matrix<typename Base<T>::type>& d, Int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::UpdateRealPartOfDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    const Int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            UpdateRealPart( j, j+offset, d.Get(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            UpdateRealPart( j-offset, j, d.Get(j,0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
Matrix<T,Int>::UpdateImagPartOfDiagonal
( const Matrix<typename Base<T>::type>& d, Int offset )
{ 
#ifndef RELEASE
    PushCallStack("Matrix::UpdateImagPartOfDiagonal");
    if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
        throw std::logic_error("d is not a column-vector of the right length");
#endif
    if( !IsComplex<T>::val )
        throw std::logic_error("Cannot update imaginary part of real matrix");

    const Int diagLength = DiagonalLength(offset);
    if( offset >= 0 )
        for( Int j=0; j<diagLength; ++j )
            UpdateImagPart( j, j+offset, d.Get(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            UpdateImagPart( j-offset, j, d.Get(j,0) );
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
#endif
    Empty();

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
#endif
    Empty();

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
#endif
    Empty();

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
#endif
    Empty();

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
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        std::ostringstream msg;
        msg << "Trying to view outside of a Matrix: "
            << "up to (" << i+height-1 << "," << j+width-1 << ") "
            << "of " << A.Height() << " x " << A.Width() << " Matrix.";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    Empty();

    height_     = height;
    width_      = width;
    ldim_       = A.LDim();
    data_       = A.Buffer(i,j);
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
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        std::ostringstream msg;
        msg << "Trying to view outside of a Matrix: "
            << "up to (" << i+height-1 << "," << j+width-1 << ") "
            << "of " << A.Height() << " x " << A.Width() << " Matrix.";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    Empty();

    height_     = height;
    width_      = width;
    ldim_       = A.LDim();
    lockedData_ = A.LockedBuffer(i,j);
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
    if( AL.Height() != AR.Height() )
        throw std::logic_error("1x2 must have consistent height to combine");
    if( AL.LDim() != AR.LDim() )
        throw std::logic_error("1x2 must have consistent ldims to combine");
    if( AR.Buffer() != (AL.Buffer()+AL.LDim()*AL.Width()) )
        throw std::logic_error("1x2 must have contiguous memory");
#endif
    Empty();

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
    if( AL.Height() != AR.Height() )
        throw std::logic_error("1x2 must have consistent height to combine");
    if( AL.LDim() != AR.LDim() )
        throw std::logic_error("1x2 must have consistent ldims to combine");
    if( AR.LockedBuffer() != (AL.LockedBuffer()+AL.LDim()*AL.Width()) )
        throw std::logic_error("1x2 must have contiguous memory");
#endif
    Empty();

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
    if( AT.Width() != AB.Width() )
        throw std::logic_error("2x1 must have consistent width to combine");
    if( AT.LDim() != AB.LDim() )
        throw std::logic_error("2x1 must have consistent ldim to combine");
    if( AB.Buffer() != (AT.Buffer() + AT.Height()) )
        throw std::logic_error("2x1 must have contiguous memory");
#endif
    Empty();

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
    if( AT.Width() != AB.Width() )
        throw std::logic_error("2x1 must have consistent width to combine");
    if( AT.LDim() != AB.LDim() )
        throw std::logic_error("2x1 must have consistent ldim to combine");
    if( AB.LockedBuffer() != (AT.LockedBuffer()+AT.Height()) )
        throw std::logic_error("2x1 must have contiguous memory");
#endif
    Empty();

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
    Empty();

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
    Empty();

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
        MemCopy( &data_[j*ldim], &data[j*ldimOfA], height );
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
Matrix<T,Int>::AssertValidEntry( Int i, Int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::AssertValidEntry");
#endif
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
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif /* ELEMENTAL_MATRIX_HPP */

