/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

//
// Constructors
//

template<typename T,typename Int>
Matrix<T,Int>::Matrix()
: viewing_(false), locked_(false),
  height_(0), width_(0), ldim_(1), data_(0), lockedData_(0),
  memory_()
{ }

template<typename T,typename Int>
Matrix<T,Int>::Matrix( Int height, Int width )
: viewing_(false), locked_(false),
  height_(height), width_(width), ldim_(std::max(height,1)), lockedData_(0)
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Matrix");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    memory_.Require( ldim_*width );
    data_ = memory_.Buffer();
}

template<typename T,typename Int>
Matrix<T,Int>::Matrix
( Int height, Int width, Int ldim )
: viewing_(false), locked_(false),
  height_(height), width_(width), ldim_(ldim), lockedData_(0)
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Matrix");
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
}

template<typename T,typename Int>
Matrix<T,Int>::Matrix
( Int height, Int width, const T* buffer, Int ldim )
: viewing_(true), locked_(true),
  height_(height), width_(width), ldim_(ldim), data_(0), lockedData_(buffer)
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Matrix");
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
}

template<typename T,typename Int>
Matrix<T,Int>::Matrix
( Int height, Int width, T* buffer, Int ldim )
: viewing_(true), locked_(false),
  height_(height), width_(width), ldim_(ldim), data_(buffer), lockedData_(0)
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Matrix");
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
}

template<typename T,typename Int>
Matrix<T,Int>::Matrix
( const Matrix<T,Int>& A )
: viewing_(false), locked_(false), 
  height_(0), width_(0), ldim_(1), data_(0), lockedData_(0)
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Matrix( const Matrix& )");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ("You just tried to construct a Matrix with itself!");
}

//
// Destructor
//

template<typename T,typename Int>
Matrix<T,Int>::~Matrix()
{ }

//
// Basic information
//

template<typename T,typename Int>
Int 
Matrix<T,Int>::Height() const
{ return height_; }

template<typename T,typename Int>
Int
Matrix<T,Int>::Width() const
{ return width_; }

template<typename T,typename Int>
Int
Matrix<T,Int>::DiagonalLength( Int offset ) const
{ return elem::DiagonalLength(height_,width_,offset); }

template<typename T,typename Int>
Int
Matrix<T,Int>::LDim() const
{ return ldim_; }

template<typename T,typename Int>
Int
Matrix<T,Int>::MemorySize() const
{ return memory_.Size(); }

template<typename T,typename Int>
T*
Matrix<T,Int>::Buffer()
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Buffer");
    if( locked_ )
        throw std::logic_error
        ("Cannot return non-const buffer of locked Matrix");
#endif
    return data_;
}

template<typename T,typename Int>
const T*
Matrix<T,Int>::LockedBuffer() const
{
    if( locked_ )
        return lockedData_;
    else
        return data_;
}

template<typename T,typename Int>
T*
Matrix<T,Int>::Buffer( Int i, Int j )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Buffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( locked_ )
        throw std::logic_error
        ("Cannot return non-const buffer of locked Matrix");
#endif
    return &data_[i+j*ldim_];
}

template<typename T,typename Int>
const T*
Matrix<T,Int>::LockedBuffer( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::LockedBuffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
#endif
    if( locked_ )
        return &lockedData_[i+j*ldim_];
    else
        return &data_[i+j*ldim_];
}

//
// Entry manipulation
//

template<typename T,typename Int>
T
Matrix<T,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Get");
    AssertValidEntry( i, j );
#endif
    if( lockedData_ )
        return lockedData_[i+j*ldim_];
    else
        return data_[i+j*ldim_];
}

template<typename T,typename Int>
void
Matrix<T,Int>::Set( Int i, Int j, T alpha ) 
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Set");
    AssertValidEntry( i, j );
    if( lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    data_[i+j*ldim_] = alpha;
}

template<typename T,typename Int>
void
Matrix<T,Int>::Update( Int i, Int j, T alpha ) 
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Update");
    AssertValidEntry( i, j );
    if( lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    data_[i+j*ldim_] += alpha;
}

template<typename T,typename Int>
void
Matrix<T,Int>::GetDiagonal( Matrix<T,Int>& d, Int offset ) const
{ 
#ifndef RELEASE
    CallStackEntry entry("Matrix::GetDiagonal");
    if( d.Locked() )
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
}

template<typename T,typename Int>
void
Matrix<T,Int>::SetDiagonal( const Matrix<T,Int>& d, Int offset )
{ 
#ifndef RELEASE
    CallStackEntry entry("Matrix::SetDiagonal");
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
}

template<typename T,typename Int>
void
Matrix<T,Int>::UpdateDiagonal( const Matrix<T,Int>& d, Int offset )
{ 
#ifndef RELEASE
    CallStackEntry entry("Matrix::UpdateDiagonal");
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
}

template<typename T,typename Int>
BASE(T)
Matrix<T,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::GetRealPart");
    AssertValidEntry( i, j );
#endif
    if( lockedData_ )
        return RealPart(lockedData_[i+j*ldim_]);
    else
        return RealPart(data_[i+j*ldim_]);
}

template<typename T,typename Int>
BASE(T)
Matrix<T,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::GetImagPart");
    AssertValidEntry( i, j );
#endif
    if( lockedData_ )
        return ImagPart(lockedData_[i+j*ldim_]);
    else
        return ImagPart(data_[i+j*ldim_]);
}

template<typename T,typename Int>
void
Matrix<T,Int>::SetRealPart( Int i, Int j, BASE(T) alpha )
{ SetRealPartHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
void
Matrix<T,Int>::SetRealPartHelper<Z>::Func
( Matrix<Z,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::SetRealPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    parent.data_[i+j*parent.ldim_] = alpha;
}
    
template<typename T,typename Int>
template<typename Z>
void
Matrix<T,Int>::SetRealPartHelper<Complex<Z> >::Func
( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::SetRealPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    const Z beta = parent.data_[i+j*parent.ldim_].imag;
    parent.data_[i+j*parent.ldim_] = Complex<Z>( alpha, beta );
}

template<typename T,typename Int>
void
Matrix<T,Int>::SetImagPart( Int i, Int j, BASE(T) alpha ) 
{ SetImagPartHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
void
Matrix<T,Int>::SetImagPartHelper<Z>::Func
( Matrix<Z,Int>& parent, Int i, Int j, Z alpha ) 
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::SetImagPartHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T,typename Int>
template<typename Z>
void
Matrix<T,Int>::SetImagPartHelper<Complex<Z> >::Func
( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha ) 
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::SetImagPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    const Z beta = parent.data_[i+j*parent.ldim_].real;
    parent.data_[i+j*parent.ldim_] = Complex<Z>( beta, alpha );
}

template<typename T,typename Int>
void
Matrix<T,Int>::UpdateRealPart( Int i, Int j, BASE(T) alpha )
{ UpdateRealPartHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
void
Matrix<T,Int>::UpdateRealPartHelper<Z>::Func
( Matrix<Z,Int>& parent, Int i, Int j, Z alpha ) 
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::UpdateRealPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    parent.data_[i+j*parent.ldim_] += alpha;
}
    
template<typename T,typename Int>
template<typename Z>
void
Matrix<T,Int>::UpdateRealPartHelper<Complex<Z> >::Func
( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::UpdateRealPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    const Complex<Z> beta = parent.data_[i+j*parent.ldim_];
    parent.data_[i+j*parent.ldim_] = Complex<Z>( beta.real+alpha, beta.imag );
}

template<typename T,typename Int>
void
Matrix<T,Int>::UpdateImagPart( Int i, Int j, BASE(T) alpha ) 
{ UpdateImagPartHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
void
Matrix<T,Int>::UpdateImagPartHelper<Z>::Func
( Matrix<Z,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::UpdateImagPartHelper::Func");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}
    
template<typename T,typename Int>
template<typename Z>
void
Matrix<T,Int>::UpdateImagPartHelper<Complex<Z> >::Func
( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::UpdateImagPartHelper::Func");
    parent.AssertValidEntry( i, j );
    if( parent.lockedData_ )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    const Complex<Z> beta = parent.data_[i+j*parent.ldim_];
    parent.data_[i+j*parent.ldim_] = Complex<Z>( beta.real, beta.imag+alpha );
}

template<typename T,typename Int>
void
Matrix<T,Int>::GetRealPartOfDiagonal( Matrix<BASE(T)>& d, Int offset ) const
{ 
#ifndef RELEASE
    CallStackEntry entry("Matrix::GetRealPartOfDiagonal");
    if( d.Locked() )
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
}

template<typename T,typename Int>
void
Matrix<T,Int>::GetImagPartOfDiagonal( Matrix<BASE(T)>& d, Int offset ) const
{ 
#ifndef RELEASE
    CallStackEntry entry("Matrix::GetImagPartOfDiagonal");
    if( d.Locked() )
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
}

template<typename T,typename Int>
void
Matrix<T,Int>::SetRealPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset )
{ 
#ifndef RELEASE
    CallStackEntry entry("Matrix::SetRealPartOfDiagonal");
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
}

template<typename T,typename Int>
void
Matrix<T,Int>::SetImagPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset )
{ 
#ifndef RELEASE
    CallStackEntry entry("Matrix::SetImagPartOfDiagonal");
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
}

template<typename T,typename Int>
void
Matrix<T,Int>::UpdateRealPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset )
{ 
#ifndef RELEASE
    CallStackEntry entry("Matrix::UpdateRealPartOfDiagonal");
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
}

template<typename T,typename Int>
void
Matrix<T,Int>::UpdateImagPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset )
{ 
#ifndef RELEASE
    CallStackEntry entry("Matrix::UpdateImagPartOfDiagonal");
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
}

template<typename T,typename Int>
void
Matrix<T,Int>::Control
( Int height, Int width, T* buffer, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Control");
#endif
    Empty();

    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewing_ = false;
    locked_ = false;
}

template<typename T,typename Int>
void
Matrix<T,Int>::Attach
( Int height, Int width, T* buffer, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Attach");
#endif
    Empty();

    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewing_ = true;
    locked_ = false;
}

template<typename T,typename Int>
void
Matrix<T,Int>::LockedAttach
( Int height, Int width, const T* buffer, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::LockedAttach");
#endif
    Empty();

    height_ = height;
    width_ = width;
    ldim_ = ldim;
    lockedData_ = buffer;
    viewing_ = true;
    locked_ = true;
}

//
// Viewing other Matrix instances
//

template<typename T,typename Int>
bool
Matrix<T,Int>::Viewing() const
{ return viewing_; }

template<typename T,typename Int>
bool
Matrix<T,Int>::Locked() const
{ return locked_; }

//
// Utilities
//

template<typename T,typename Int>
const Matrix<T,Int>&
Matrix<T,Int>::operator=( const Matrix<T,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::operator=");
    if( locked_ )
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
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
        MemCopy( &data_[j*ldim], &data[j*ldimOfA], height );
    return *this;
}

template<typename T,typename Int>
void
Matrix<T,Int>::Empty()
{
    memory_.Empty();
    height_ = 0;
    width_ = 0;
    ldim_ = 1;
    data_ = 0;
    lockedData_ = 0;
    viewing_ = false;
    locked_ = false;
}

template<typename T,typename Int>
void
Matrix<T,Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::ResizeTo(height,width)");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    const bool reallocate = ( height > height_ || width > width_ );
#ifndef RELEASE
    if( viewing_ && reallocate )
        throw std::logic_error("Cannot increase the size of a view");
#endif
    height_ = height;
    width_ = width;

    // Only change the ldim when necessary. Simply 'shrink' our view if 
    // possible.
    if( reallocate )
    {
        ldim_ = std::max( height, 1 );
        memory_.Require(ldim_*width);
        data_ = memory_.Buffer();
    }
}

template<typename T,typename Int>
void
Matrix<T,Int>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::ResizeTo(height,width,ldim)");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( ldim < height )
    {
        std::ostringstream msg;
        msg << "Tried to set ldim(" << ldim << ") < height (" << height << ")";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const bool reallocate = 
        ( height > height_ || width > width_ || ldim != ldim_ );
#ifndef RELEASE
    if( viewing_ && reallocate )
        throw std::logic_error("Illogical ResizeTo on viewed data");

#endif
    height_ = height;
    width_ = width;
    ldim_ = ldim;

    if( reallocate )
    {
        memory_.Require(ldim*width);
        data_ = memory_.Buffer();
    }
}

template<typename T,typename Int>
void
Matrix<T,Int>::AssertValidEntry( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::AssertValidEntry");
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
}

template class Matrix<int,int>;
#ifndef DISABLE_FLOAT
template class Matrix<float,int>;
#endif // ifndef DISABLE_FLOAT
template class Matrix<double,int>;
#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
template class Matrix<Complex<float>,int>;
#endif // ifndef DISABLE_FLOAT
template class Matrix<Complex<double>,int>;
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
