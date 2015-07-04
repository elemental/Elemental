/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Public routines
// ###############

// Constructors and destructors
// ============================

template<typename Ring>
Matrix<Ring>::Matrix( bool fixed )
: viewType_( fixed ? OWNER_FIXED : OWNER ),
  height_(0), width_(0), ldim_(1), 
  data_(nullptr)
{ }

template<typename Ring>
Matrix<Ring>::Matrix( Int height, Int width, bool fixed )
: viewType_( fixed ? OWNER_FIXED : OWNER ),
  height_(height), width_(width), ldim_(Max(height,1))
{
    DEBUG_ONLY(
      CSE cse("Matrix::Matrix");
      AssertValidDimensions( height, width );
    )
    memory_.Require( ldim_ * width );
    data_ = memory_.Buffer();
    // TODO: Consider explicitly zeroing
}

template<typename Ring>
Matrix<Ring>::Matrix
( Int height, Int width, Int ldim, bool fixed )
: viewType_( fixed ? OWNER_FIXED : OWNER ),
  height_(height), width_(width), ldim_(ldim)
{
    DEBUG_ONLY(
      CSE cse("Matrix::Matrix");
      AssertValidDimensions( height, width, ldim );
    )
    memory_.Require( ldim*width );
    data_ = memory_.Buffer();
}

template<typename Ring>
Matrix<Ring>::Matrix
( Int height, Int width, const Ring* buffer, Int ldim, bool fixed )
: viewType_( fixed ? LOCKED_VIEW_FIXED: LOCKED_VIEW ),
  height_(height), width_(width), ldim_(ldim), 
  data_(buffer)
{
    DEBUG_ONLY(
      CSE cse("Matrix::Matrix");
      AssertValidDimensions( height, width, ldim );
    )
}

template<typename Ring>
Matrix<Ring>::Matrix
( Int height, Int width, Ring* buffer, Int ldim, bool fixed )
: viewType_( fixed ? VIEW_FIXED: VIEW ),
  height_(height), width_(width), ldim_(ldim), 
  data_(buffer)
{
    DEBUG_ONLY(
      CSE cse("Matrix::Matrix");
      AssertValidDimensions( height, width, ldim );
    )
}

template<typename Ring>
Matrix<Ring>::Matrix( const Matrix<Ring>& A )
: viewType_( OWNER ),
  height_(0), width_(0), ldim_(1), 
  data_(nullptr)
{
    DEBUG_ONLY(CSE cse("Matrix::Matrix( const Matrix& )"))
    if( &A != this )
        *this = A;
    else
        LogicError("You just tried to construct a Matrix with itself!");
}

template<typename Ring>
Matrix<Ring>::Matrix( Matrix<Ring>&& A ) EL_NOEXCEPT
: viewType_(A.viewType_),
  height_(A.height_), width_(A.width_), ldim_(A.ldim_),
  data_(nullptr), memory_(std::move(A.memory_))
{ std::swap( data_, A.data_ ); }

template<typename Ring>
Matrix<Ring>::~Matrix() { }

// Assignment and reconfiguration
// ==============================

template<typename Ring>
void Matrix<Ring>::Empty()
{
    DEBUG_ONLY(
      CSE cse("Matrix::Empty()");
      if( FixedSize() )
          LogicError("Cannot empty a fixed-size matrix" );
    )
    Empty_();
}

template<typename Ring>
void Matrix<Ring>::Resize( Int height, Int width )
{
    DEBUG_ONLY(
      CSE cse("Matrix::Resize(height,width)");
      AssertValidDimensions( height, width );
      if( FixedSize() && ( height != height_ || width != width_ ) )
          LogicError("Cannot change the size of this matrix");
      if( Viewing() && ( height > height_ || width > width_ ) )
          LogicError("Cannot increase the size of this matrix");
    )
    Resize_( height, width );
}

template<typename Ring>
void Matrix<Ring>::Resize( Int height, Int width, Int ldim )
{
    DEBUG_ONLY(
      CSE cse("Matrix::Resize(height,width,ldim)");
      AssertValidDimensions( height, width, ldim );
      if( FixedSize() && 
          ( height != height_ || width != width_ || ldim != ldim_ ) )
          LogicError("Cannot change the size of this matrix");
      if( Viewing() && (height > height_ || width > width_ || ldim != ldim_) )
          LogicError("Cannot increase the size of this matrix");
    )
    Resize_( height, width, ldim );
}

template<typename Ring>
void Matrix<Ring>::Attach( Int height, Int width, Ring* buffer, Int ldim )
{
    DEBUG_ONLY(
      CSE cse("Matrix::Attach");
      if( FixedSize() )
          LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    Attach_( height, width, buffer, ldim );
}

template<typename Ring>
void Matrix<Ring>::LockedAttach
( Int height, Int width, const Ring* buffer, Int ldim )
{
    DEBUG_ONLY(
      CSE cse("Matrix::LockedAttach");
      if( FixedSize() )
          LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    LockedAttach_( height, width, buffer, ldim );
}

template<typename Ring>
void Matrix<Ring>::Control( Int height, Int width, Ring* buffer, Int ldim )
{
    DEBUG_ONLY(
      CSE cse("Matrix::Control");
      if( FixedSize() )
          LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    Control_( height, width, buffer, ldim );
}

// Operator overloading
// ====================

// Return a view
// -------------
template<typename Ring>
Matrix<Ring> Matrix<Ring>::operator()( Range<Int> I, Range<Int> J )
{
    DEBUG_ONLY(CSE cse("Matrix::operator()"))
    if( this->Locked() )
        return LockedView( *this, I, J );
    else
        return View( *this, I, J );
}

template<typename Ring>
const Matrix<Ring> Matrix<Ring>::operator()( Range<Int> I, Range<Int> J ) const
{
    DEBUG_ONLY(CSE cse("Matrix::operator()"))
    return LockedView( *this, I, J );
}

// Make a copy
// -----------
template<typename Ring>
const Matrix<Ring>& Matrix<Ring>::operator=( const Matrix<Ring>& A )
{
    DEBUG_ONLY(CSE cse("Matrix::operator="))
    Copy( A, *this );
    return *this;
}

// Move assignment
// ---------------
template<typename Ring>
Matrix<Ring>& Matrix<Ring>::operator=( Matrix<Ring>&& A )
{
    DEBUG_ONLY(CSE cse("Matrix::operator=( Matrix&& )"))
    if( Viewing() || A.Viewing() )
    {
        operator=( (const Matrix<Ring>&)A );
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
template<typename Ring>
const Matrix<Ring>& Matrix<Ring>::operator*=( Ring alpha )
{
    DEBUG_ONLY(CSE cse("Matrix::operator*=(Ring)"))
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename Ring>
const Matrix<Ring>& Matrix<Ring>::operator+=( const Matrix<Ring>& A )
{
    DEBUG_ONLY(CSE cse("Matrix::operator+=( const Matrix<Ring>& )"))
    Axpy( Ring(1), A, *this );
    return *this;
}

template<typename Ring>
const Matrix<Ring>& Matrix<Ring>::operator-=( const Matrix<Ring>& A )
{
    DEBUG_ONLY(CSE cse("Matrix::operator+=( const Matrix<Ring>& )"))
    Axpy( Ring(-1), A, *this );
    return *this;
}

// Basic queries
// =============

template<typename Ring>
Int Matrix<Ring>::Height() const { return height_; }

template<typename Ring>
Int Matrix<Ring>::Width() const { return width_; }

template<typename Ring>
Int Matrix<Ring>::LDim() const { return ldim_; }

template<typename Ring>
Int Matrix<Ring>::MemorySize() const { return memory_.Size(); }

template<typename Ring>
Int Matrix<Ring>::DiagonalLength( Int offset ) const
{ return El::DiagonalLength(height_,width_,offset); }

template<typename Ring>
Ring* Matrix<Ring>::Buffer()
{
    DEBUG_ONLY(
        CSE cse("Matrix::Buffer");
        if( Locked() )
            LogicError("Cannot return non-const buffer of locked Matrix");
    )
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return const_cast<Ring*>(data_);
}

template<typename Ring>
Ring* Matrix<Ring>::Buffer( Int i, Int j )
{
    DEBUG_ONLY(
      CSE cse("Matrix::Buffer");
      if( Locked() )
          LogicError("Cannot return non-const buffer of locked Matrix");
    )
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return &const_cast<Ring*>(data_)[i+j*ldim_];
}

template<typename Ring>
const Ring* Matrix<Ring>::LockedBuffer() const { return data_; }

template<typename Ring>
const Ring* Matrix<Ring>::LockedBuffer( Int i, Int j ) const
{
    DEBUG_ONLY(CSE cse("Matrix::LockedBuffer"))
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    return &data_[i+j*ldim_];
}

template<typename Ring>
bool Matrix<Ring>::Viewing() const { return IsViewing( viewType_ ); }

template<typename Ring>
bool Matrix<Ring>::FixedSize() const { return IsFixedSize( viewType_ ); }

template<typename Ring>
bool Matrix<Ring>::Locked() const { return IsLocked( viewType_ ); }

// Single-entry manipulation
// =========================

template<typename Ring>
Ring Matrix<Ring>::Get( Int i, Int j ) const
{
    DEBUG_ONLY(
      CSE cse("Matrix::Get");
      AssertValidEntry( i, j );
    )
    return Get_( i, j );
}

template<typename Ring>
Base<Ring> Matrix<Ring>::GetRealPart( Int i, Int j ) const
{
    DEBUG_ONLY(
      CSE cse("Matrix::GetRealPart");
      AssertValidEntry( i, j );
    )
    return El::RealPart( Get_( i, j ) );
}

template<typename Ring>
Base<Ring> Matrix<Ring>::GetImagPart( Int i, Int j ) const
{
    DEBUG_ONLY(
      CSE cse("Matrix::GetImagPart");
      AssertValidEntry( i, j );
    )
    return El::ImagPart( Get_( i, j ) );
}

template<typename Ring>
void Matrix<Ring>::Set( Int i, Int j, Ring alpha ) 
{
    DEBUG_ONLY(
      CSE cse("Matrix::Set");
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    Set_( i, j ) = alpha;
}

template<typename Ring>
void Matrix<Ring>::Set( const Entry<Ring>& entry )
{ Set( entry.i, entry.j, entry.value ); }

template<typename Ring>
void Matrix<Ring>::SetRealPart( Int i, Int j, Base<Ring> alpha )
{
    DEBUG_ONLY(
      CSE cse("Matrix::SetRealPart");
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    El::SetRealPart( Set_( i, j ), alpha );
}

template<typename Ring>
void Matrix<Ring>::SetRealPart( const Entry<Base<Ring>>& entry )
{ SetRealPart( entry.i, entry.j, entry.value ); }

template<typename Ring>
void Matrix<Ring>::SetImagPart( Int i, Int j, Base<Ring> alpha )
{
    DEBUG_ONLY(
      CSE cse("Matrix::SetImagPart");
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    ComplainIfReal();
    El::SetImagPart( Set_( i, j ), alpha );
}

template<typename Ring>
void Matrix<Ring>::SetImagPart( const Entry<Base<Ring>>& entry )
{ SetImagPart( entry.i, entry.j, entry.value ); }

template<typename Ring>
void Matrix<Ring>::Update( Int i, Int j, Ring alpha ) 
{
    DEBUG_ONLY(
      CSE cse("Matrix::Update");
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    Set_( i, j ) += alpha;
}

template<typename Ring>
void Matrix<Ring>::Update( const Entry<Ring>& entry )
{ Update( entry.i, entry.j, entry.value ); }

template<typename Ring>
void Matrix<Ring>::UpdateRealPart( Int i, Int j, Base<Ring> alpha )
{
    DEBUG_ONLY(
      CSE cse("Matrix::UpdateRealPart");
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    El::UpdateRealPart( Set_( i, j ), alpha );
}

template<typename Ring>
void Matrix<Ring>::UpdateRealPart( const Entry<Base<Ring>>& entry )
{ UpdateRealPart( entry.i, entry.j, entry.value ); }

template<typename Ring>
void Matrix<Ring>::UpdateImagPart( Int i, Int j, Base<Ring> alpha )
{
    DEBUG_ONLY(
      CSE cse("Matrix::UpdateImagPart");
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    ComplainIfReal();
    El::UpdateImagPart( Set_( i, j ), alpha );
}

template<typename Ring>
void Matrix<Ring>::UpdateImagPart( const Entry<Base<Ring>>& entry )
{ UpdateImagPart( entry.i, entry.j, entry.value ); }

template<typename Ring>
void Matrix<Ring>::MakeReal( Int i, Int j )
{
    DEBUG_ONLY(
      CSE cse("Matrix::MakeReal");
      AssertValidEntry( i, j );
      if( Locked() )
          LogicError("Cannot modify data of locked matrices");
    )
    Set( i, j, GetRealPart(i,j) );
}

template<typename Ring>
void Matrix<Ring>::Conjugate( Int i, Int j )
{
    DEBUG_ONLY(
      CSE cse("Matrix::Conjugate");
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
template<typename Ring>
void Matrix<Ring>::ShallowSwap( Matrix<Ring>& A )
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

template<typename Ring>
void Matrix<Ring>::Empty_()
{
    memory_.Empty();
    height_ = 0;
    width_ = 0;
    ldim_ = 1;
    data_ = nullptr;
    viewType_ = (ViewType)( viewType_ & ~LOCKED_VIEW );
}

template<typename Ring>
void Matrix<Ring>::Attach_( Int height, Int width, Ring* buffer, Int ldim )
{
    memory_.Empty();
    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewType_ = (ViewType)( ( viewType_ & ~LOCKED_OWNER ) | VIEW );
}

template<typename Ring>
void Matrix<Ring>::LockedAttach_
( Int height, Int width, const Ring* buffer, Int ldim )
{
    memory_.Empty();
    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewType_ = (ViewType)( viewType_ | VIEW );
}

template<typename Ring>
void Matrix<Ring>::Control_( Int height, Int width, Ring* buffer, Int ldim )
{
    memory_.Empty();
    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewType_ = (ViewType)( viewType_ & ~LOCKED_VIEW );
}

// Return a reference to a single entry without error-checking
// ===========================================================
template<typename Ring>
const Ring& Matrix<Ring>::Get_( Int i, Int j ) const 
{ 
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    return data_[i+j*ldim_]; 
}

template<typename Ring>
Ring& Matrix<Ring>::Set_( Int i, Int j ) 
{
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return (const_cast<Ring*>(data_))[i+j*ldim_];
}

// Assertions
// ==========

template<typename Ring>
void Matrix<Ring>::ComplainIfReal() const
{ 
    if( !IsComplex<Ring>::val )
        LogicError("Called complex-only routine with real data");
}

template<typename Ring>
void Matrix<Ring>::AssertValidDimensions( Int height, Int width ) const
{
    DEBUG_ONLY(CSE cse("Matrix::AssertValidDimensions"))
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
}

template<typename Ring>
void Matrix<Ring>::AssertValidDimensions( Int height, Int width, Int ldim ) const
{
    DEBUG_ONLY(CSE cse("Matrix::AssertValidDimensions"))
    AssertValidDimensions( height, width );
    if( ldim < height )
        LogicError("Leading dimension must be no less than height");
    if( ldim == 0 )
        LogicError("Leading dimension cannot be zero (for BLAS compatibility)");
}

template<typename Ring>
void Matrix<Ring>::AssertValidEntry( Int i, Int j ) const
{
    DEBUG_ONLY(CSE cse("Matrix::AssertValidEntry"))
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    if( i < 0 || j < 0 )
        LogicError("Indices must be non-negative");
    if( i >= Height() || j >= Width() )
        LogicError
        ("Out of bounds: (",i,",",j,") of ",Height()," x ",Width()," Matrix");
}

template<typename Ring>
void Matrix<Ring>::Resize_( Int height, Int width )
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

template<typename Ring>
void Matrix<Ring>::Resize_( Int height, Int width, Int ldim )
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

#define PROTO(Ring) template class Matrix<Ring>;
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
