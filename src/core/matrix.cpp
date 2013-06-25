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
Matrix<T,Int>::Matrix( bool fixed )
: viewType_( fixed ? OWNER_FIXED : OWNER ),
  height_(0), width_(0), ldim_(1), 
  data_(0)
{ }

template<typename T,typename Int>
Matrix<T,Int>::Matrix( Int height, Int width, bool fixed )
: viewType_( fixed ? OWNER_FIXED : OWNER ),
  height_(height), width_(width), ldim_(std::max(height,1))
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Matrix");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
#endif
    memory_.Require( ldim_ * width );
    data_ = memory_.Buffer();
}

template<typename T,typename Int>
Matrix<T,Int>::Matrix
( Int height, Int width, Int ldim, bool fixed )
: viewType_( fixed ? OWNER_FIXED : OWNER ),
  height_(height), width_(width), ldim_(ldim)
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
( Int height, Int width, const T* buffer, Int ldim, bool fixed )
: viewType_( fixed ? LOCKED_VIEW_FIXED: LOCKED_VIEW ),
  height_(height), width_(width), ldim_(ldim), 
  data_(buffer)
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
( Int height, Int width, T* buffer, Int ldim, bool fixed )
: viewType_( fixed ? VIEW_FIXED: VIEW ),
  height_(height), width_(width), ldim_(ldim), 
  data_(buffer)
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
: viewType_( OWNER ),
  height_(0), width_(0), ldim_(1), 
  data_(0)
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
bool
Matrix<T,Int>::Owner() const
{ return IsOwner( viewType_ ); }

template<typename T,typename Int>
bool
Matrix<T,Int>::Viewing() const
{ return !IsOwner( viewType_ ); }

template<typename T,typename Int>
bool
Matrix<T,Int>::Shrinkable() const
{ return IsShrinkable( viewType_ ); }

template<typename T,typename Int>
bool
Matrix<T,Int>::FixedSize() const
{ return !IsShrinkable( viewType_ ); }

template<typename T,typename Int>
bool
Matrix<T,Int>::Locked() const
{ return IsLocked( viewType_ ); }

template<typename T,typename Int>
T*
Matrix<T,Int>::Buffer()
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Buffer");
    if( Locked() )
        throw std::logic_error
        ("Cannot return non-const buffer of locked Matrix");
#endif
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return const_cast<T*>(data_);
}

template<typename T,typename Int>
const T*
Matrix<T,Int>::LockedBuffer() const
{ return data_; }

template<typename T,typename Int>
T*
Matrix<T,Int>::Buffer( Int i, Int j )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Buffer");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( Locked() )
        throw std::logic_error
        ("Cannot return non-const buffer of locked Matrix");
#endif
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return &const_cast<T*>(data_)[i+j*ldim_];
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
    return &data_[i+j*ldim_];
}

//
// Entry manipulation
//

template<typename T,typename Int>
const T&
Matrix<T,Int>::Get_( Int i, Int j ) const
{ return data_[i+j*ldim_]; }

template<typename T,typename Int>
T&
Matrix<T,Int>::Set_( Int i, Int j ) 
{
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return (const_cast<T*>(data_))[i+j*ldim_];
}

template<typename T,typename Int>
T
Matrix<T,Int>::Get( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Get");
    AssertValidEntry( i, j );
#endif
    return Get_( i, j );
}

template<typename T,typename Int>
void
Matrix<T,Int>::Set( Int i, Int j, T alpha ) 
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Set");
    AssertValidEntry( i, j );
    if( Locked() )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    Set_( i, j ) = alpha;
}

template<typename T,typename Int>
void
Matrix<T,Int>::Update( Int i, Int j, T alpha ) 
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Update");
    AssertValidEntry( i, j );
    if( Locked() )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    Set_( i, j ) += alpha;
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
            d.Set_( j, 0 ) = Get_(j,j+offset);
    else
        for( Int j=0; j<diagLength; ++j )
            d.Set_( j, 0 ) = Get_(j-offset,j);
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
            Set_( j, j+offset ) = d.Get_(j,0);
    else
        for( Int j=0; j<diagLength; ++j )
            Set_( j-offset, j ) = d.Get_(j,0);
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
            Set_( j, j+offset ) += d.Get(j,0);
    else
        for( Int j=0; j<diagLength; ++j )
            Set_( j-offset, j ) += d.Get(j,0);
}

template<typename T,typename Int>
BASE(T)
Matrix<T,Int>::GetRealPart( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::GetRealPart");
    AssertValidEntry( i, j );
#endif
    return elem::RealPart( Get_( i, j ) );
}

template<typename T,typename Int>
BASE(T)
Matrix<T,Int>::GetImagPart( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::GetImagPart");
    AssertValidEntry( i, j );
#endif
    return elem::ImagPart( Get_( i, j ) );
}

template<typename T,typename Int>
void 
Matrix<T,Int>::SetRealPart( Int i, Int j, BASE(T) alpha )
{
#ifndef RELEASE
    CallStackEntry cse("Matrix::SetRealPart");
    AssertValidEntry( i, j );
    if( Locked() )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    elem::SetRealPart( Set_( i, j ), alpha );
}

template<typename T,typename Int>
void 
Matrix<T,Int>::SetImagPart( Int i, Int j, BASE(T) alpha )
{
#ifndef RELEASE
    CallStackEntry cse("Matrix::SetImagPart");
    AssertValidEntry( i, j );
    if( Locked() )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    elem::SetImagPart( Set_( i, j ), alpha );
}

template<typename T,typename Int>
void 
Matrix<T,Int>::UpdateRealPart( Int i, Int j, BASE(T) alpha )
{
#ifndef RELEASE
    CallStackEntry cse("Matrix::UpdateRealPart");
    AssertValidEntry( i, j );
    if( Locked() )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    elem::UpdateRealPart( Set_( i, j ), alpha );
}

template<typename T,typename Int>
void 
Matrix<T,Int>::UpdateImagPart( Int i, Int j, BASE(T) alpha )
{
#ifndef RELEASE
    CallStackEntry cse("Matrix::UpdateImagPart");
    AssertValidEntry( i, j );
    if( Locked() )
        throw std::logic_error("Cannot modify data of locked matrices");
#endif
    elem::UpdateImagPart( Set_( i, j ), alpha );
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
            d.Set_( j, 0 ) = elem::RealPart( Get_(j,j+offset) );
    else
        for( Int j=0; j<diagLength; ++j )
            d.Set_( j, 0 ) = elem::RealPart( Get_(j-offset,j) );
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
            d.Set_( j, 0 ) = elem::ImagPart( Get_(j,j+offset) );
    else
        for( Int j=0; j<diagLength; ++j )
            d.Set_( j, 0 ) = elem::ImagPart( Get_(j-offset,j) );
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
            elem::SetRealPart( Set_(j,j+offset), d.Get_(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            elem::SetRealPart( Set_(j-offset,j), d.Get_(j,0) );
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
            elem::SetImagPart( Set_(j,j+offset), d.Get_(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            elem::SetImagPart( Set_(j-offset,j), d.Get_(j,0) );
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
            elem::UpdateRealPart( Set_(j,j+offset), d.Get_(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            elem::UpdateRealPart( Set_(j-offset,j), d.Get_(j,0) );
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
            elem::UpdateImagPart( Set_(j,j+offset), d.Get_(j,0) );
    else
        for( Int j=0; j<diagLength; ++j )
            elem::UpdateImagPart( Set_(j-offset,j), d.Get_(j,0) );
}

template<typename T,typename Int>
void
Matrix<T,Int>::Control( Int height, Int width, T* buffer, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Control");
#endif
    Empty();

    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewType_ = OWNER;
}

template<typename T,typename Int>
void
Matrix<T,Int>::Attach( Int height, Int width, T* buffer, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::Attach");
    if( FixedSize() )
        throw std::logic_error( "Cannot attach a new buffer to a view with fixed size" );
#endif
    Empty();

    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewType_ = VIEW;
}

template<typename T,typename Int>
void
Matrix<T,Int>::LockedAttach
( Int height, Int width, const T* buffer, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::LockedAttach");
    if( FixedSize() )
        throw std::logic_error( "Cannot attach a new buffer to a view with fixed size" );
#endif
    Empty();

    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewType_ = LOCKED_VIEW;
}

//
// Utilities
//

template<typename T,typename Int>
const Matrix<T,Int>&
Matrix<T,Int>::operator=( const Matrix<T,Int>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::operator=");
    if( Locked() )
        throw std::logic_error("Cannot assign to a locked view");
    if( viewType_ != OWNER && (A.Height() != Height() || A.Width() != Width()) )
        throw std::logic_error
        ("Cannot assign to a view of different dimensions");
#endif
    if( viewType_ == OWNER )
        ResizeTo( A.Height(), A.Width() );
    const Int height = Height();
    const Int width = Width();
    const Int ldim = LDim();
    const Int ldimOfA = A.LDim();
    const T* src = A.LockedBuffer();
    T* dst = this->Buffer();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( Int j=0; j<width; ++j )
        MemCopy( &dst[j*ldim], &src[j*ldimOfA], height );
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
    viewType_ = OWNER;
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
    if( height == height_ && width == width_ )
        return;
    const bool reallocate = ( height > height_ || width > width_ );
#ifndef RELEASE
    if ( FixedSize() )
        throw std::logic_error("Cannot change the size of this matrix");
    if ( Viewing() && reallocate )
        throw std::logic_error("Cannot increase the size of this matrix");
#endif
    height_ = height;
    width_ = width;
    // Only change the ldim when necessary. Simply 'shrink' our view if 
    // possible.
    if( reallocate )
    {
        ldim_ = std::max( height, 1 );
        memory_.Require( ldim_ * width );
        data_ = memory_.Buffer();
    }
}

template<typename T,typename Int>
void
Matrix<T,Int>::ResizeTo( Int height, Int width, Int ldim )
{
#ifndef RELEASE
    CallStackEntry entry("Matrix::ResizeTo(height,width,ldim)");
    if ( FixedSize() )
        throw std::logic_error("Cannot change the size of this matrix");
    if( ldim < height )
    {
        std::ostringstream msg;
        msg << "Tried to set ldim(" << ldim << ") < height (" << height << ")";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const bool reallocate = ( height > height_ || width > width_ || ldim != ldim_ );
#ifndef RELEASE
    if ( Viewing() && reallocate )
        throw std::logic_error( "Cannot increase the size of this matrix" );
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
