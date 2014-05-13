/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

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
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Matrix");
        AssertValidDimensions( height, width );
    )
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
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Matrix");
        AssertValidDimensions( height, width, ldim );
    )
    memory_.Require( ldim*width );
    data_ = memory_.Buffer();
}

template<typename T>
Matrix<T>::Matrix
( Int height, Int width, const T* buffer, Int ldim, bool fixed )
: viewType_( fixed ? LOCKED_VIEW_FIXED: LOCKED_VIEW ),
  height_(height), width_(width), ldim_(ldim), 
  data_(buffer)
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Matrix");
        AssertValidDimensions( height, width, ldim );
    )
}

template<typename T>
Matrix<T>::Matrix
( Int height, Int width, T* buffer, Int ldim, bool fixed )
: viewType_( fixed ? VIEW_FIXED: VIEW ),
  height_(height), width_(width), ldim_(ldim), 
  data_(buffer)
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Matrix");
        AssertValidDimensions( height, width, ldim );
    )
}

template<typename T>
Matrix<T>::Matrix( const Matrix<T>& A )
: viewType_( OWNER ),
  height_(0), width_(0), ldim_(1), 
  data_(nullptr)
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::Matrix( const Matrix& )"))
    if( &A != this )
        *this = A;
    else
        LogicError("You just tried to construct a Matrix with itself!");
}

template<typename T>
Matrix<T>::Matrix( Matrix<T>&& A ) EL_NOEXCEPT
: viewType_(A.viewType_),
  height_(A.height_), width_(A.width_), ldim_(A.ldim_),
  data_(nullptr), memory_(std::move(A.memory_))
{ std::swap( data_, A.data_ ); }

template<typename T>
Matrix<T>::~Matrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
Matrix<T>&
Matrix<T>::operator=( const Matrix<T>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::operator=");
        if( Locked() )
            LogicError("Cannot assign to a locked view");
        if( viewType_ != OWNER && 
            (A.Height() != Height() || A.Width() != Width()) )
            LogicError("Cannot assign to a view of different dimensions");
    )
    if( viewType_ == OWNER )
        Resize( A.Height(), A.Width() );
    const Int height = Height();
    const Int width = Width();
    const Int ldim = LDim();
    const Int ldimOfA = A.LDim();
    const T* src = A.LockedBuffer();
    T* dst = Buffer();
    EL_PARALLEL_FOR
    for( Int j=0; j<width; ++j )
        MemCopy( &dst[j*ldim], &src[j*ldimOfA], height );
    return *this;
}

template<typename T>
Matrix<T>&
Matrix<T>::operator=( Matrix<T>&& A )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::operator=( Matrix&& )"))
    if( Viewing() && !A.Viewing() )
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

template<typename T>
void
Matrix<T>::Empty()
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Empty()");
        if( FixedSize() )
            LogicError("Cannot empty a fixed-size matrix" );
    )
    Empty_();
}

template<typename T>
void
Matrix<T>::Resize( Int height, Int width )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Resize(height,width)");
        AssertValidDimensions( height, width );
        if( FixedSize() && ( height != height_ || width != width_ ) )
            LogicError("Cannot change the size of this matrix");
        if( Viewing() && ( height > height_ || width > width_ ) )
            LogicError("Cannot increase the size of this matrix");
    )
    Resize_( height, width );
}

template<typename T>
void
Matrix<T>::Resize( Int height, Int width, Int ldim )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Resize(height,width,ldim)");
        AssertValidDimensions( height, width, ldim );
        if( FixedSize() && 
            ( height != height_ || width != width_ || ldim != ldim_ ) )
            LogicError("Cannot change the size of this matrix");
        if( Viewing() && (height > height_ || width > width_ || ldim != ldim_) )
            LogicError("Cannot increase the size of this matrix");
    )
    Resize_( height, width, ldim );
}

template<typename T>
void
Matrix<T>::Attach( Int height, Int width, T* buffer, Int ldim )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Attach");
        if( FixedSize() )
            LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    Attach_( height, width, buffer, ldim );
}

template<typename T>
void
Matrix<T>::LockedAttach
( Int height, Int width, const T* buffer, Int ldim )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::LockedAttach");
        if( FixedSize() )
            LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    LockedAttach_( height, width, buffer, ldim );
}

template<typename T>
void
Matrix<T>::Control( Int height, Int width, T* buffer, Int ldim )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Control");
        if( FixedSize() )
            LogicError("Cannot attach a new buffer to a view with fixed size");
    )
    Control_( height, width, buffer, ldim );
}

// Basic queries
// =============

template<typename T>
Int Matrix<T>::Height() const { return height_; }

template<typename T>
Int Matrix<T>::Width() const { return width_; }

template<typename T>
Int Matrix<T>::LDim() const { return ldim_; }

template<typename T>
Int Matrix<T>::MemorySize() const { return memory_.Size(); }

template<typename T>
Int Matrix<T>::DiagonalLength( Int offset ) const
{ return El::DiagonalLength(height_,width_,offset); }

template<typename T>
T*
Matrix<T>::Buffer()
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Buffer");
        if( Locked() )
            LogicError("Cannot return non-const buffer of locked Matrix");
    )
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return const_cast<T*>(data_);
}

template<typename T>
T*
Matrix<T>::Buffer( Int i, Int j )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Buffer");
        if( Locked() )
            LogicError("Cannot return non-const buffer of locked Matrix");
    )
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return &const_cast<T*>(data_)[i+j*ldim_];
}

template<typename T>
const T* Matrix<T>::LockedBuffer() const { return data_; }

template<typename T>
const T*
Matrix<T>::LockedBuffer( Int i, Int j ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::LockedBuffer"))
    return &data_[i+j*ldim_];
}

template<typename T>
bool Matrix<T>::Viewing() const { return IsViewing( viewType_ ); }

template<typename T>
bool Matrix<T>::FixedSize() const { return IsFixedSize( viewType_ ); }

template<typename T>
bool Matrix<T>::Locked() const { return IsLocked( viewType_ ); }

// Single-entry manipulation
// =========================

template<typename T>
T
Matrix<T>::Get( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Get");
        AssertValidEntry( i, j );
    )
    return Get_( i, j );
}

template<typename T>
Base<T>
Matrix<T>::GetRealPart( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::GetRealPart");
        AssertValidEntry( i, j );
    )
    return El::RealPart( Get_( i, j ) );
}

template<typename T>
Base<T>
Matrix<T>::GetImagPart( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::GetImagPart");
        AssertValidEntry( i, j );
    )
    return El::ImagPart( Get_( i, j ) );
}

template<typename T>
void
Matrix<T>::Set( Int i, Int j, T alpha ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Set");
        AssertValidEntry( i, j );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    Set_( i, j ) = alpha;
}

template<typename T>
void 
Matrix<T>::SetRealPart( Int i, Int j, Base<T> alpha )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::SetRealPart");
        AssertValidEntry( i, j );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    El::SetRealPart( Set_( i, j ), alpha );
}

template<typename T>
void 
Matrix<T>::SetImagPart( Int i, Int j, Base<T> alpha )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::SetImagPart");
        AssertValidEntry( i, j );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    ComplainIfReal();
    El::SetImagPart( Set_( i, j ), alpha );
}

template<typename T>
void
Matrix<T>::Update( Int i, Int j, T alpha ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Update");
        AssertValidEntry( i, j );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    Set_( i, j ) += alpha;
}

template<typename T>
void 
Matrix<T>::UpdateRealPart( Int i, Int j, Base<T> alpha )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::UpdateRealPart");
        AssertValidEntry( i, j );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    El::UpdateRealPart( Set_( i, j ), alpha );
}

template<typename T>
void 
Matrix<T>::UpdateImagPart( Int i, Int j, Base<T> alpha )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::UpdateImagPart");
        AssertValidEntry( i, j );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    ComplainIfReal();
    El::UpdateImagPart( Set_( i, j ), alpha );
}

template<typename T>
void
Matrix<T>::MakeReal( Int i, Int j )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::MakeReal");
        AssertValidEntry( i, j );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    Set( i, j, GetRealPart(i,j) );
}

template<typename T>
void
Matrix<T>::Conjugate( Int i, Int j )
{
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::Conjugate");
        AssertValidEntry( i, j );
        if( Locked() )
            LogicError("Cannot modify data of locked matrices");
    )
    Set( i, j, El::Conj(Get(i,j)) );
}

// Diagonal manipulation
// =====================

template<typename T>
void
Matrix<T>::GetDiagonal( Matrix<T>& d, Int offset ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::GetDiagonal");
        if( d.Locked() )
            LogicError("d must not be a locked view");
    )
    const Int diagLength = DiagonalLength(offset);
    d.Resize( diagLength, 1 );
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        d.Set_( k, 0 ) = Get_(k+iOff,k+jOff);
}

template<typename T>
void
Matrix<T>::GetRealPartOfDiagonal( Matrix<Base<T>>& d, Int offset ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::GetRealPartOfDiagonal");
        if( d.Locked() )
            LogicError("d must not be a locked view");
    )
    const Int diagLength = DiagonalLength(offset);
    d.Resize( diagLength, 1 );
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        d.Set_( k, 0 ) = El::RealPart( Get_(k+iOff,k+jOff) );
}

template<typename T>
void
Matrix<T>::GetImagPartOfDiagonal( Matrix<Base<T>>& d, Int offset ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::GetImagPartOfDiagonal");
        if( d.Locked() )
            LogicError("d must not be a locked view");
    )
    const Int diagLength = DiagonalLength(offset);
    d.Resize( diagLength, 1 );
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        d.Set_( k, 0 ) = El::ImagPart( Get_(k+iOff,k+jOff) );
}

template<typename T>
Matrix<T>
Matrix<T>::GetDiagonal( Int offset ) const
{ 
    Matrix<T> d;
    GetDiagonal( d, offset );
    return d;
}

template<typename T>
Matrix<Base<T>>
Matrix<T>::GetRealPartOfDiagonal( Int offset ) const
{ 
    Matrix<Base<T>> d;
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
Matrix<Base<T>>
Matrix<T>::GetImagPartOfDiagonal( Int offset ) const
{ 
    Matrix<Base<T>> d;
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T>
void
Matrix<T>::SetDiagonal( const Matrix<T>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::SetDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        Set_( k+iOff, k+jOff ) = d.Get_(k,0);
}

template<typename T>
void
Matrix<T>::SetRealPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::SetRealPartOfDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        El::SetRealPart( Set_(k+iOff,k+jOff), d.Get_(k,0) );
}

template<typename T>
void
Matrix<T>::SetImagPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::SetImagPartOfDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    ComplainIfReal();
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        El::SetImagPart( Set_(k+iOff,k+jOff), d.Get_(k,0) );
}

template<typename T>
void
Matrix<T>::UpdateDiagonal( const Matrix<T>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::UpdateDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        Set_( k+iOff, k+jOff ) += d.Get(k,0);
}

template<typename T>
void
Matrix<T>::UpdateRealPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::UpdateRealPartOfDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        El::UpdateRealPart( Set_(k+iOff,k+jOff), d.Get_(k,0) );
}

template<typename T>
void
Matrix<T>::UpdateImagPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("Matrix::UpdateImagPartOfDiagonal");
        if( d.Height() != DiagonalLength(offset) || d.Width() != 1 )
            LogicError("d is not a column-vector of the right length");
    )
    ComplainIfReal();
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        El::UpdateImagPart( Set_(k+iOff,k+jOff), d.Get_(k,0) );
}

template<typename T>
void
Matrix<T>::MakeDiagonalReal( Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::MakeDiagonalReal"))
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        Set( k+iOff, k+jOff, RealPart(Get(k+iOff,k+jOff)) );
}

template<typename T>
void
Matrix<T>::ConjugateDiagonal( Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::ConjugateDiagonal"))
    const Int diagLength = DiagonalLength(offset);
    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );
    for( Int k=0; k<diagLength; ++k )
        Set( k+iOff, k+jOff, Conj(Get(k+iOff,k+jOff)) );
}

// Arbitrary submatrix manipulation
// ================================

template<typename T>
void
Matrix<T>::GetSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Matrix<T>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::GetSubmatrix"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    ASub.Resize( m, n );
    T* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    const T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(AssertValidEntry(rowInd[i],colInd[j]))
            bufSub[i+j*ldSub] = buf[rowInd[i]+jSub*ld];
        }
    }
}

template<typename T>
void
Matrix<T>::GetRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Matrix<Base<T>>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::GetRealPartOfSubmatrix"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    ASub.Resize( m, n );
    Base<T>* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    const T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(AssertValidEntry(rowInd[i],colInd[j]))
            bufSub[i+j*ldSub] = RealPart(buf[rowInd[i]+jSub*ld]);
        }
    }
}

template<typename T>
void
Matrix<T>::GetImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Matrix<Base<T>>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::GetImagPartOfSubmatrix"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    ASub.Resize( m, n );
    Base<T>* bufSub = ASub.Buffer();
    const Int ldSub = ASub.LDim();
    const T* buf = LockedBuffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(AssertValidEntry(rowInd[i],colInd[j]))
            bufSub[i+j*ldSub] = ImagPart(buf[rowInd[i]+jSub*ld]);
        }
    }
}

template<typename T>
Matrix<T>
Matrix<T>::GetSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    Matrix<T> ASub;
    GetSubmatrix( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
Matrix<Base<T>>
Matrix<T>::GetRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    Matrix<Base<T>> ASub;
    GetRealPartOfSubmatrix( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
Matrix<Base<T>>
Matrix<T>::GetImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    Matrix<Base<T>> ASub;
    GetImagPartOfSubmatrix( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
void
Matrix<T>::SetSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::SetSubmatrix"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const T* bufSub = ASub.LockedBuffer();
    const Int ldSub = ASub.LDim();
    T* buf = Buffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(AssertValidEntry(rowInd[i],colInd[j]))
            buf[rowInd[i]+jSub*ld] = bufSub[i+j*ldSub];
        }
    }
}

template<typename T>
void
Matrix<T>::SetRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::SetRealPartOfSubmatrix"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const Base<T>* bufSub = ASub.LockedBuffer();
    const Int ldSub = ASub.LDim();
    T* buf = Buffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(AssertValidEntry(rowInd[i],colInd[j]))
            El::SetRealPart( buf[rowInd[i]+jSub*ld], bufSub[i+j*ldSub] );
        }
    }
}

template<typename T>
void
Matrix<T>::SetImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::SetImagPartOfSubmatrix"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const Base<T>* bufSub = ASub.LockedBuffer();
    const Int ldSub = ASub.LDim();
    T* buf = Buffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(AssertValidEntry(rowInd[i],colInd[j]))
            El::SetImagPart( buf[rowInd[i]+jSub*ld], bufSub[i+j*ldSub] );
        }
    }
}

template<typename T>
void
Matrix<T>::UpdateSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  T alpha, const Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::UpdateSubmatrix"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const T* bufSub = ASub.LockedBuffer();
    const Int ldSub = ASub.LDim();
    T* buf = Buffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(AssertValidEntry(rowInd[i],colInd[j]))
            buf[rowInd[i]+jSub*ld] += alpha*bufSub[i+j*ldSub];
        }
    }
}

template<typename T>
void
Matrix<T>::UpdateRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Base<T> alpha, const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::UpdateRealPartOfSubmatrix"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const Base<T>* bufSub = ASub.LockedBuffer();
    const Int ldSub = ASub.LDim();
    T* buf = Buffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(AssertValidEntry(rowInd[i],colInd[j]))
            El::UpdateRealPart
            ( buf[rowInd[i]+jSub*ld], alpha*bufSub[i+j*ldSub] );
        }
    }
}

template<typename T>
void
Matrix<T>::UpdateImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Base<T> alpha, const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::UpdateImagPartOfSubmatrix"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    const Base<T>* bufSub = ASub.LockedBuffer();
    const Int ldSub = ASub.LDim();
    T* buf = Buffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(AssertValidEntry(rowInd[i],colInd[j]))
            El:: UpdateImagPart
            ( buf[rowInd[i]+jSub*ld], alpha*bufSub[i+j*ldSub] );
        }
    }
}

template<typename T>
void
Matrix<T>::MakeSubmatrixReal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::MakeSubmatrixReal"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    T* buf = Buffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(AssertValidEntry(rowInd[i],colInd[j]))
            buf[rowInd[i]+jSub*ld] = RealPart(buf[rowInd[i]+jSub*ld]);
        }
    }
}

template<typename T>
void
Matrix<T>::ConjugateSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::ConjugateSubmatrix"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    T* buf = Buffer();
    const Int ld = LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int jSub = colInd[j];
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(AssertValidEntry(rowInd[i],colInd[j]))
            buf[rowInd[i]+jSub*ld] = Conj( buf[rowInd[i]+jSub*ld] );
        }
    }
}

// Private routines
// ################

// Exchange metadata with another matrix
// =====================================
template<typename T>
void
Matrix<T>::ShallowSwap( Matrix<T>& A )
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
void
Matrix<T>::Empty_()
{
    memory_.Empty();
    height_ = 0;
    width_ = 0;
    ldim_ = 1;
    data_ = nullptr;
    viewType_ = (ViewType)( viewType_ & ~LOCKED_VIEW );
}

template<typename T>
void
Matrix<T>::Attach_( Int height, Int width, T* buffer, Int ldim )
{
    memory_.Empty();
    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewType_ = (ViewType)( ( viewType_ & ~LOCKED_OWNER ) | VIEW );
}

template<typename T>
void
Matrix<T>::LockedAttach_( Int height, Int width, const T* buffer, Int ldim )
{
    memory_.Empty();
    height_ = height;
    width_ = width;
    ldim_ = ldim;
    data_ = buffer;
    viewType_ = (ViewType)( viewType_ | VIEW );
}

template<typename T>
void
Matrix<T>::Control_( Int height, Int width, T* buffer, Int ldim )
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
template<typename T>
const T&
Matrix<T>::Get_( Int i, Int j ) const
{ return data_[i+j*ldim_]; }

template<typename T>
T&
Matrix<T>::Set_( Int i, Int j ) 
{
    // NOTE: This const_cast has been carefully considered and should be safe
    //       since the underlying data should be non-const if this is called.
    return (const_cast<T*>(data_))[i+j*ldim_];
}

// Assertions
// ==========

template<typename T>
void
Matrix<T>::ComplainIfReal() const
{ 
    if( !IsComplex<T>::val )
        LogicError("Called complex-only routine with real data");
}

template<typename T>
void
Matrix<T>::AssertValidDimensions( Int height, Int width ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::AssertValidDimensions"))
    if( height < 0 || width < 0 )
        LogicError("Height and width must be non-negative");
}

template<typename T>
void
Matrix<T>::AssertValidDimensions( Int height, Int width, Int ldim ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::AssertValidDimensions"))
    AssertValidDimensions( height, width );
    if( ldim < height )
        LogicError("Leading dimension must be no less than height");
    if( ldim == 0 )
        LogicError("Leading dimension cannot be zero (for BLAS compatibility)");
}

template<typename T>
void
Matrix<T>::AssertValidEntry( Int i, Int j ) const
{
    DEBUG_ONLY(CallStackEntry cse("Matrix::AssertValidEntry"))
    if( i < 0 || j < 0 )
        LogicError("Indices must be non-negative");
    if( i >= Height() || j >= Width() )
        LogicError
        ("Out of bounds: (",i,",",j,") of ",Height()," x ",Width()," Matrix");
}

template<typename T>
void
Matrix<T>::Resize_( Int height, Int width )
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
void
Matrix<T>::Resize_( Int height, Int width, Int ldim )
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

// Instantiate for {Int,Real,Complex<Real>} for each Real in {float,double}
// ########################################################################
template class Matrix<Int>;
#ifndef EL_DISABLE_FLOAT
template class Matrix<float>;
#endif // ifndef EL_DISABLE_FLOAT
template class Matrix<double>;
#ifndef EL_DISABLE_COMPLEX
#ifndef EL_DISABLE_FLOAT
template class Matrix<Complex<float>>;
#endif // ifndef EL_DISABLE_FLOAT
template class Matrix<Complex<double>>;
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El
