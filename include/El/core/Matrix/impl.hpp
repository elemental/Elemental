/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_MATRIX_IMPL_HPP
#define EL_MATRIX_IMPL_HPP

namespace El {

// Public routines
// ###############

// Constructors and destructors
// ============================

template<typename T>
Matrix<T>::Matrix( bool fixed )
: viewType_( fixed ? OWNER_FIXED : OWNER ),
  height_(0), width_(0), ldim_(1), 
  data_(nullptr)
{ }

template<typename T>
Matrix<T>::Matrix( Int height, Int width, bool fixed )
: viewType_( fixed ? OWNER_FIXED : OWNER ),
  height_(height), width_(width), ldim_(Max(height,1))
{
    DEBUG_CSE
    DEBUG_ONLY(AssertValidDimensions( height, width ))
    memory_.Require( ldim_ * width );
    data_ = memory_.Buffer();
    // TODO: Consider explicitly zeroing
}

template<typename T>
Matrix<T>::Matrix
( Int height, Int width, Int ldim, bool fixed )
: viewType_( fixed ? OWNER_FIXED : OWNER ),
  height_(height), width_(width), ldim_(ldim)
{
    DEBUG_CSE
    DEBUG_ONLY(AssertValidDimensions( height, width, ldim ))
    memory_.Require( ldim*width );
    data_ = memory_.Buffer();
}

template<typename T>
Matrix<T>::Matrix
( Int height, Int width, const T* buffer, Int ldim, bool fixed )
: viewType_( fixed ? LOCKED_VIEW_FIXED: LOCKED_VIEW ),
  height_(height), width_(width), ldim_(ldim), 
  data_(const_cast<T*>(buffer))
{
    DEBUG_CSE
    DEBUG_ONLY(AssertValidDimensions( height, width, ldim ))
}

template<typename T>
Matrix<T>::Matrix
( Int height, Int width, T* buffer, Int ldim, bool fixed )
: viewType_( fixed ? VIEW_FIXED: VIEW ),
  height_(height), width_(width), ldim_(ldim), 
  data_(buffer)
{
    DEBUG_CSE
    DEBUG_ONLY(AssertValidDimensions( height, width, ldim ))
}

template<typename T>
Matrix<T>::Matrix( const Matrix<T>& A )
: viewType_( OWNER ),
  height_(0), width_(0), ldim_(1), 
  data_(nullptr)
{
    DEBUG_CSE
    if( &A != this )
        *this = A;
    else
        LogicError("You just tried to construct a Matrix with itself!");
}

template<typename T>
Matrix<T>::Matrix( Matrix<T>&& A ) EL_NO_EXCEPT
: viewType_(A.viewType_),
  height_(A.height_), width_(A.width_), ldim_(A.ldim_),
  memory_(std::move(A.memory_)), data_(nullptr)
{ std::swap( data_, A.data_ ); }

template<typename T>
Matrix<T>::~Matrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
void Matrix<T>::Empty( bool freeMemory )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( FixedSize() )
          LogicError("Cannot empty a fixed-size matrix" );
    )
    Empty_( freeMemory );
}

template<typename T>
void Matrix<T>::Resize( Int height, Int width )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertValidDimensions( height, width );
      if( FixedSize() && ( height != height_ || width != width_ ) )
      {
          DumpCallStack();
          LogicError
          ("Cannot resize this matrix from ",
           height_," x ",width_," to ",height," x ",width);
      }
      if( Viewing() && ( height > height_ || width > width_ ) )
          LogicError("Cannot increase the size of this matrix");
    )
    Resize_( height, width );
}

template<typename T>
void Matrix<T>::Resize( Int height, Int width, Int ldim )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertValidDimensions( height, width, ldim );
      if( FixedSize() && 
          ( height != height_ || width != width_ || ldim != ldim_ ) )
      {
          DumpCallStack();
          LogicError
          ("Cannot resize this matrix from ",
           height_," x ",width_," (",ldim_,") to ",
           height," x ",width," (",ldim,")");
      }
      if( Viewing() && (height > height_ || width > width_ || ldim != ldim_) )
          LogicError("Cannot increase the size of this matrix");
    )
    Resize_( height, width, ldim );
}

template<typename T>
void Matrix<T>::Attach( Int height, Int width, T* buffer, Int ldim )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( FixedSize() )
          LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    Attach_( height, width, buffer, ldim );
}

template<typename T>
void Matrix<T>::LockedAttach( Int height, Int width, const T* buffer, Int ldim )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( FixedSize() )
          LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    LockedAttach_( height, width, buffer, ldim );
}

template<typename T>
void Matrix<T>::Control( Int height, Int width, T* buffer, Int ldim )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( FixedSize() )
          LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    Control_( height, width, buffer, ldim );
}

// Operator overloading
// ====================

// Return a view
// -------------
template<typename T>
Matrix<T> Matrix<T>::operator()( Range<Int> I, Range<Int> J )
{
    DEBUG_CSE
    if( this->Locked() )
        return LockedView( *this, I, J );
    else
        return View( *this, I, J );
}

template<typename T>
const Matrix<T> Matrix<T>::operator()( Range<Int> I, Range<Int> J ) const
{
    DEBUG_CSE
    return LockedView( *this, I, J );
}

// Return a (potentially non-contiguous) subset of indices
// -------------------------------------------------------
template<typename T>
Matrix<T> Matrix<T>::operator()
( Range<Int> I, const vector<Int>& J ) const
{
    DEBUG_CSE
    Matrix<T> ASub; 
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
Matrix<T> Matrix<T>::operator()
( const vector<Int>& I, Range<Int> J ) const
{
    DEBUG_CSE
    Matrix<T> ASub;
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
Matrix<T> Matrix<T>::operator()
( const vector<Int>& I, const vector<Int>& J ) const
{
    DEBUG_CSE
    Matrix<T> ASub;
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

// Make a copy
// -----------
template<typename T>
const Matrix<T>& Matrix<T>::operator=( const Matrix<T>& A )
{
    DEBUG_CSE
    Copy( A, *this );
    return *this;
}

// Move assignment
// ---------------
template<typename T>
Matrix<T>& Matrix<T>::operator=( Matrix<T>&& A )
{
    DEBUG_CSE
    if( Viewing() || A.Viewing() )
    {
        operator=( (const Matrix<T>&)A );
    }
    else
    {
        memory_.ShallowSwap( A.memory_ );
        std::swap( data_, A.data_ );
        viewType_ = A.viewType_;
        height_ = A.height_;
        width_ = A.width_;
        ldim_ = A.ldim_;
    }
    return *this;
}

// Rescaling
// ---------
template<typename T>
const Matrix<T>& Matrix<T>::operator*=( T alpha )
{
    DEBUG_CSE
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename T>
const Matrix<T>& Matrix<T>::operator+=( const Matrix<T>& A )
{
    DEBUG_CSE
    Axpy( T(1), A, *this );
    return *this;
}

template<typename T>
const Matrix<T>& Matrix<T>::operator-=( const Matrix<T>& A )
{
    DEBUG_CSE
    Axpy( T(-1), A, *this );
    return *this;
}

// Basic queries
// =============

template<typename T>
Int Matrix<T>::Height() const EL_NO_EXCEPT { return height_; }

template<typename T>
Int Matrix<T>::Width() const EL_NO_EXCEPT { return width_; }

template<typename T>
Int Matrix<T>::LDim() const EL_NO_EXCEPT { return ldim_; }

template<typename T>
Int Matrix<T>::MemorySize() const EL_NO_EXCEPT { return memory_.Size(); }

template<typename T>
Int Matrix<T>::DiagonalLength( Int offset ) const EL_NO_EXCEPT
{ return El::DiagonalLength(height_,width_,offset); }

template<typename T>
T* Matrix<T>::Buffer() EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( Locked() )
          LogicError("Cannot return non-const buffer of locked Matrix");
    )
    return data_;
}

template<typename T>
T* Matrix<T>::Buffer( Int i, Int j ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( Locked() )
          LogicError("Cannot return non-const buffer of locked Matrix");
    )
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    return &data_[i+j*ldim_];
}

template<typename T>
const T* Matrix<T>::LockedBuffer() const EL_NO_EXCEPT { return data_; }

template<typename T>
const T* Matrix<T>::LockedBuffer( Int i, Int j ) const EL_NO_EXCEPT
{
    DEBUG_CSE
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    return &data_[i+j*ldim_];
}

template<typename T>
bool Matrix<T>::Viewing() const EL_NO_EXCEPT
{ return IsViewing( viewType_ ); }

template<typename T>
bool Matrix<T>::FixedSize() const EL_NO_EXCEPT
{ return IsFixedSize( viewType_ ); }

template<typename T>
bool Matrix<T>::Locked() const EL_NO_EXCEPT
{ return IsLocked( viewType_ ); }

template<typename T>
void Matrix<T>::SetViewType( El::ViewType viewType ) EL_NO_EXCEPT
{ viewType_ = viewType; }

template<typename T>
El::ViewType Matrix<T>::ViewType() const EL_NO_EXCEPT
{ return viewType_; }

// Single-entry manipulation
// =========================

template<typename T>
T Matrix<T>::Get( Int i, Int j ) const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(AssertValidEntry( i, j ))
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    return CRef( i, j );
}

template<typename T>
Base<T> Matrix<T>::GetRealPart( Int i, Int j ) const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(AssertValidEntry( i, j ))
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    return El::RealPart( CRef( i, j ) );
}

template<typename T>
Base<T> Matrix<T>::GetImagPart( Int i, Int j ) const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(AssertValidEntry( i, j ))
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    return El::ImagPart( CRef( i, j ) );
}

template<typename T>
void Matrix<T>::Set( Int i, Int j, T alpha ) 
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    Ref( i, j ) = alpha;
}

template<typename T>
void Matrix<T>::Set( const Entry<T>& entry )
EL_NO_RELEASE_EXCEPT
{ Set( entry.i, entry.j, entry.value ); }

template<typename T>
void Matrix<T>::SetRealPart( Int i, Int j, Base<T> alpha )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    El::SetRealPart( Ref( i, j ), alpha );
}

template<typename T>
void Matrix<T>::SetRealPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ SetRealPart( entry.i, entry.j, entry.value ); }

template<typename T>
void Matrix<T>::SetImagPart( Int i, Int j, Base<T> alpha )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    El::SetImagPart( Ref( i, j ), alpha );
}

template<typename T>
void Matrix<T>::SetImagPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ SetImagPart( entry.i, entry.j, entry.value ); }

template<typename T>
void Matrix<T>::Update( Int i, Int j, T alpha ) 
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    Ref( i, j ) += alpha;
}

template<typename T>
void Matrix<T>::Update( const Entry<T>& entry )
EL_NO_RELEASE_EXCEPT
{ Update( entry.i, entry.j, entry.value ); }

template<typename T>
void Matrix<T>::UpdateRealPart( Int i, Int j, Base<T> alpha )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    El::UpdateRealPart( Ref( i, j ), alpha );
}

template<typename T>
void Matrix<T>::UpdateRealPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ UpdateRealPart( entry.i, entry.j, entry.value ); }

template<typename T>
void Matrix<T>::UpdateImagPart( Int i, Int j, Base<T> alpha )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    El::UpdateImagPart( Ref( i, j ), alpha );
}

template<typename T>
void Matrix<T>::UpdateImagPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ UpdateImagPart( entry.i, entry.j, entry.value ); }

template<typename T>
void Matrix<T>::MakeReal( Int i, Int j )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    Set( i, j, GetRealPart(i,j) );
}

template<typename T>
void Matrix<T>::Conjugate( Int i, Int j )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    Set( i, j, El::Conj(Get(i,j)) );
}

// Private routines
// ################

// Exchange metadata with another matrix
// =====================================
template<typename T>
void Matrix<T>::ShallowSwap( Matrix<T>& A )
{
    memory_.ShallowSwap( A.memory_ );
    std::swap( data_, A.data_ );
    std::swap( viewType_, A.viewType_ );
    std::swap( height_, A.height_ );
    std::swap( width_, A.width_ );
    std::swap( ldim_, A.ldim_ );
}

// Reconfigure without error-checking
// ==================================

template<typename T>
void Matrix<T>::Empty_( bool freeMemory )
{
    if( freeMemory )
        memory_.Empty();
    height_ = 0;
    width_ = 0;
    ldim_ = 1;
    data_ = nullptr;
    viewType_ = (El::ViewType)( viewType_ & ~LOCKED_VIEW );
}

template<typename T>
void Matrix<T>::Attach_( Int height, Int width, T* buffer, Int ldim )
{
    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewType_ = (El::ViewType)( ( viewType_ & ~LOCKED_OWNER ) | VIEW );
}

template<typename T>
void Matrix<T>::LockedAttach_
( Int height, Int width, const T* buffer, Int ldim )
{
    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = const_cast<T*>(buffer);
    viewType_ = (El::ViewType)( viewType_ | LOCKED_VIEW );
}

template<typename T>
void Matrix<T>::Control_( Int height, Int width, T* buffer, Int ldim )
{
    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewType_ = (El::ViewType)( viewType_ & ~LOCKED_VIEW );
}

// Return a reference to a single entry without error-checking
// ===========================================================
template<typename T>
const T& Matrix<T>::CRef( Int i, Int j ) const 
EL_NO_RELEASE_EXCEPT
{ 
    return data_[i+j*ldim_]; 
}

template<typename T>
const T& Matrix<T>::operator()( Int i, Int j ) const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(AssertValidEntry( i, j ))
    return data_[i+j*ldim_];
}

template<typename T>
T& Matrix<T>::Ref( Int i, Int j ) 
EL_NO_RELEASE_EXCEPT
{
    return data_[i+j*ldim_];
}

template<typename T>
T& Matrix<T>::operator()( Int i, Int j ) 
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    return data_[i+j*ldim_];
}

// Assertions
// ==========

template<typename T>
void Matrix<T>::AssertValidDimensions( Int height, Int width ) const
{
    DEBUG_CSE
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
}

template<typename T>
void Matrix<T>::AssertValidDimensions( Int height, Int width, Int ldim ) const
{
    DEBUG_CSE
    AssertValidDimensions( height, width );
    if( ldim < height )
        LogicError("Leading dimension must be no less than height");
    if( ldim == 0 )
        LogicError("Leading dimension cannot be zero (for BLAS compatibility)");
}

template<typename T>
void Matrix<T>::AssertValidEntry( Int i, Int j ) const
{
    DEBUG_CSE
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    if( i < 0 || j < 0 )
        LogicError("Indices must be non-negative");
    if( i >= Height() || j >= Width() )
        LogicError
        ("Out of bounds: (",i,",",j,") of ",Height()," x ",Width()," Matrix");
}

template<typename T>
void Matrix<T>::Resize_( Int height, Int width )
{
    bool reallocate = height > ldim_ || width > width_;
    height_ = height;
    width_ = width;
    // Only change the ldim when necessary. Simply 'shrink' our view if 
    // possible.
    if( reallocate )
    {
        ldim_ = Max( height, 1 );
        memory_.Require( ldim_ * width );
        data_ = memory_.Buffer();
    }
}

template<typename T>
void Matrix<T>::Resize_( Int height, Int width, Int ldim )
{
    bool reallocate = height > ldim_ || width > width_ || ldim != ldim_;
    height_ = height;
    width_ = width;
    if( reallocate )
    {
        ldim_ = ldim;
        memory_.Require(ldim*width);
        data_ = memory_.Buffer();
    }
}

#ifdef EL_INSTANTIATE_CORE
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) EL_EXTERN template class Matrix<T>;
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_MATRIX_IMPL_HPP
