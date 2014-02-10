/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {

template<typename T>
using ADM = AbstractDistMatrix<T>;

// NOTE: It seems that member functions cannot be defined using a 
//       fully-specified template alias, e.g., ADM<T>::AbstractDistMatrix(),
//       but DM<T> is okay if it is only partially specified, e.g., 
//       DM<T> = DistMatrix<T,MC,MR> and DM<T>::DistMatrix()

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
AbstractDistMatrix<T>::AbstractDistMatrix( ADM<T>&& A )
: viewType_(A.viewType_),
  height_(A.height_), width_(A.width_), 
  colConstrained_(A.colConstrained_), rowConstrained_(A.rowConstrained_),
  colAlign_(A.colAlign_), rowAlign_(A.rowAlign_),
  colShift_(A.colShift_), rowShift_(A.rowShift_), 
  root_(A.root_),
  grid_(A.grid_)
{ 
    matrix_.ShallowSwap( A.matrix_ );
    auxMemory_.ShallowSwap( A.auxMemory_ );
}

// Optional to override
// --------------------

template<typename T>
AbstractDistMatrix<T>::~AbstractDistMatrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
ADM<T>& 
AbstractDistMatrix<T>::operator=( ADM<T>&& A )
{
    auxMemory_.ShallowSwap( A.auxMemory_ );
    matrix_.ShallowSwap( A.matrix_ );
    viewType_ = A.viewType_;
    height_ = A.height_;
    width_ = A.width_;
    colConstrained_ = A.colConstrained_;
    rowConstrained_ = A.rowConstrained_;
    colAlign_ = A.colAlign_;
    rowAlign_ = A.rowAlign_;
    colShift_ = A.colShift_;
    rowShift_ = A.rowShift_;
    root_ = A.root_;
    grid_ = A.grid_;
    return *this;
}

template<typename T>
void
AbstractDistMatrix<T>::Empty()
{
    matrix_.Empty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
    colAlign_ = 0;
    rowAlign_ = 0;
    colConstrained_ = false;
    rowConstrained_ = false;
}

template<typename T>
void
AbstractDistMatrix<T>::EmptyData()
{
    matrix_.Empty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
}

template<typename T>
void
AbstractDistMatrix<T>::SetGrid( const elem::Grid& grid )
{
    if( grid_ != &grid )
    {
        Empty();
        grid_ = &grid; 
        SetShifts();
    }
}

template<typename T>
void
AbstractDistMatrix<T>::Resize( Int height, Int width )
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::Resize");
        AssertNotLocked();
    )
    height_ = height; 
    width_ = width;
    if( Participating() )
        matrix_.Resize_
        ( Length(height,ColShift(),ColStride()),
          Length(width,RowShift(),RowStride()) );
}

template<typename T>
void
AbstractDistMatrix<T>::Resize( Int height, Int width, Int ldim )
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::Resize");
        AssertNotLocked();
    )
    height_ = height; 
    width_ = width;
    if( Participating() )
        matrix_.Resize_
        ( Length(height,ColShift(),ColStride()),
          Length(width,RowShift(),RowStride()), ldim );
}

template<typename T>
void
AbstractDistMatrix<T>::MakeConsistent()
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeConsistent"))
    const elem::Grid& g = *grid_;
    const Int vcRoot = g.VCToViewingMap(0);
    Int message[8];
    if( g.ViewingRank() == vcRoot )
    {
        message[0] = viewType_;
        message[1] = height_;
        message[2] = width_;
        message[3] = colConstrained_;
        message[4] = rowConstrained_;
        message[5] = colAlign_;
        message[6] = rowAlign_;
        message[7] = root_;
    }
    mpi::Broadcast( message, 8, vcRoot, g.ViewingComm() );
    const ViewType newViewType = static_cast<ViewType>(message[0]);
    const Int newHeight = message[1]; 
    const Int newWidth = message[2];
    const bool newConstrainedCol = message[3];
    const bool newConstrainedRow = message[4];
    const Int newColAlign = message[5];
    const Int newRowAlign = message[6];
    const Int root = message[7];
    if( !Participating() )
    {
        SetRoot( root );
        viewType_ = newViewType;
        colConstrained_ = newConstrainedCol;
        rowConstrained_ = newConstrainedRow;
        colAlign_ = newColAlign;
        rowAlign_ = newRowAlign;
        SetShifts();
        Resize( newHeight, newWidth );
    }
    DEBUG_ONLY(
        else
        {
            if( viewType_ != newViewType )
                LogicError("Inconsistent ViewType");
            if( height_ != newHeight )
                LogicError("Inconsistent height");
            if( width_ != newWidth )
                LogicError("Inconsistent width");
            if( colConstrained_ != newConstrainedCol || 
                colAlign_ != newColAlign )
                LogicError("Inconsistent column constraint");
            if( rowConstrained_ != newConstrainedRow || 
                rowAlign_ != newRowAlign )
                LogicError("Inconsistent row constraint");
            if( root != root_ )
                LogicError("Inconsistent root");
        }
    )
}

// Realignment
// -----------

template<typename T>
void
AbstractDistMatrix<T>::Align( Int colAlign, Int rowAlign )
{ 
    DEBUG_ONLY(CallStackEntry cse("ADM::Align"))
    if( colAlign_ != colAlign || rowAlign_ != rowAlign )
        Empty();
    colConstrained_ = true;
    rowConstrained_ = true;
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::AlignCols( Int colAlign )
{ 
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignCols"))
    if( colAlign_ != colAlign )
        EmptyData();
    colConstrained_ = true;
    colAlign_ = colAlign;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::AlignRows( Int rowAlign )
{ 
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignRows"))
    if( rowAlign_ != rowAlign )
        EmptyData();
    rowConstrained_ = true;
    rowAlign_ = rowAlign;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::FreeAlignments() 
{ 
    colConstrained_ = false;
    rowConstrained_ = false;
}

template<typename T>
void
AbstractDistMatrix<T>::SetRoot( Int root )
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::SetRoot");
        if( root < 0 || root >= mpi::CommSize(CrossComm()) )
            LogicError("Invalid root");
    )
    if( root != root_ )
        Empty();
    root_ = root;
}

template<typename T>
void
AbstractDistMatrix<T>::AlignWith( const elem::DistData& data )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::AlignWith");
        if( colAlign_ != 0 || rowAlign_ != 0 )
            LogicError("Alignments should have been zero");
        if( colConstrained_ || rowConstrained_ )
            LogicError("There should not have been constraints");
    )
    SetGrid( *data.grid ); 
}

template<typename T>
void
AbstractDistMatrix<T>::AlignColsWith( const elem::DistData& data )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::AlignColsWith");
        if( colAlign_ != 0 )
            LogicError("Alignment should have been zero");
        if( colConstrained_ )
            LogicError("There should not have been a constraint");
    )
    SetGrid( *data.grid );
}

template<typename T>
void
AbstractDistMatrix<T>::AlignRowsWith( const elem::DistData& data )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::AlignRowsWith");
        if( rowAlign_ != 0 )
            LogicError("Alignment should have been zero");
        if( rowConstrained_ )
            LogicError("There should not have been a constraint");
    )
    SetGrid( *data.grid );
}

// Buffer attachment
// -----------------

template<typename T>
void
AbstractDistMatrix<T>::Attach
( Int height, Int width, Int colAlign, Int rowAlign,
  T* buffer, Int ldim, const elem::Grid& g, Int root )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Attach"))
    Empty();

    grid_ = &g;
    root_ = root;
    height_ = height;
    width_ = width;
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    viewType_ = VIEW;
    SetShifts();
    if( Participating() )
    {
        Int localHeight = Length(height,colShift_,ColStride());
        Int localWidth = Length(width,rowShift_,RowStride());
        matrix_.Attach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T>
void
AbstractDistMatrix<T>::Attach
( elem::Matrix<T>& A, Int colAlign, Int rowAlign, const elem::Grid& g, 
  Int root )
{
    Attach
    ( A.Height(), A.Width(), colAlign, rowAlign, A.Buffer(), A.LDim(),
      g, root );
}

template<typename T>
void
AbstractDistMatrix<T>::LockedAttach
( Int height, Int width, Int colAlign, Int rowAlign,
  const T* buffer, Int ldim, const elem::Grid& g, Int root )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::LockedAttach"))
    Empty();

    grid_ = &g;
    root_ = root;
    height_ = height;
    width_ = width;
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    viewType_ = LOCKED_VIEW;
    SetShifts();
    if( Participating() )
    {
        Int localHeight = Length(height,colShift_,ColStride());
        Int localWidth = Length(width,rowShift_,RowStride());
        matrix_.LockedAttach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T>
void
AbstractDistMatrix<T>::LockedAttach
( const elem::Matrix<T>& A, Int colAlign, Int rowAlign,
  const elem::Grid& g, Int root )
{
    LockedAttach
    ( A.Height(), A.Width(), colAlign, rowAlign,
      A.LockedBuffer(), A.LDim(), g, root );
}

// Basic queries
// =============

// Global matrix information
// -------------------------

template<typename T>
Int AbstractDistMatrix<T>::Height() const { return height_; }
template<typename T>
Int AbstractDistMatrix<T>::Width() const { return width_; }

template<typename T>
Int AbstractDistMatrix<T>::DiagonalLength( Int offset ) const
{ return elem::DiagonalLength(height_,width_,offset); }

template<typename T>
bool AbstractDistMatrix<T>::Viewing() const { return IsViewing( viewType_ ); }
template<typename T>
bool AbstractDistMatrix<T>::Locked() const { return IsLocked( viewType_ ); }

// Local matrix information
// ------------------------

template<typename T>
Int AbstractDistMatrix<T>::LocalHeight() const { return matrix_.Height(); }
template<typename T>
Int AbstractDistMatrix<T>::LocalWidth() const { return matrix_.Width(); }
template<typename T>
Int AbstractDistMatrix<T>::LDim() const { return matrix_.LDim(); }

template<typename T>
elem::Matrix<T>& 
AbstractDistMatrix<T>::Matrix() { return matrix_; }
template<typename T>
const elem::Matrix<T>& 
AbstractDistMatrix<T>::LockedMatrix() const { return matrix_; }

template<typename T>
size_t
AbstractDistMatrix<T>::AllocatedMemory() const { return matrix_.MemorySize(); }

template<typename T>
T*
AbstractDistMatrix<T>::Buffer() { return matrix_.Buffer(); }

template<typename T>
T*
AbstractDistMatrix<T>::Buffer( Int iLoc, Int jLoc )
{ return matrix_.Buffer(iLoc,jLoc); }

template<typename T>
const T*
AbstractDistMatrix<T>::LockedBuffer() const
{ return matrix_.LockedBuffer(); }

template<typename T>
const T*
AbstractDistMatrix<T>::LockedBuffer( Int iLoc, Int jLoc ) const
{ return matrix_.LockedBuffer(iLoc,jLoc); }

// Distribution information
// ------------------------

template<typename T>
const elem::Grid& AbstractDistMatrix<T>::Grid() const { return *grid_; }

template<typename T>
bool AbstractDistMatrix<T>::ColConstrained() const { return colConstrained_; }
template<typename T>
bool AbstractDistMatrix<T>::RowConstrained() const { return rowConstrained_; }

template<typename T>
Int AbstractDistMatrix<T>::ColAlign() const { return colAlign_; }
template<typename T>
Int AbstractDistMatrix<T>::RowAlign() const { return rowAlign_; }

template<typename T>
Int AbstractDistMatrix<T>::ColShift() const { return colShift_; }
template<typename T>
Int AbstractDistMatrix<T>::RowShift() const { return rowShift_; }

template<typename T>
Int
AbstractDistMatrix<T>::ColRank() const
{ 
    if( grid_->InGrid() )
        return mpi::CommRank(ColComm());
    else
        return mpi::UNDEFINED;
}

template<typename T>
Int
AbstractDistMatrix<T>::RowRank() const
{
    if( grid_->InGrid() )
        return mpi::CommRank(RowComm());
    else
        return mpi::UNDEFINED;
}

template<typename T>
Int
AbstractDistMatrix<T>::DistRank() const
{ return mpi::CommRank(DistComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::CrossRank() const
{ return mpi::CommRank(CrossComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::RedundantRank() const
{ return mpi::CommRank(RedundantComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::DistSize() const
{ return mpi::CommSize(DistComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::CrossSize() const
{ return mpi::CommSize(CrossComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::RedundantSize() const
{ return mpi::CommSize(RedundantComm()); }

template<typename T>
Int AbstractDistMatrix<T>::Root() const { return root_; }

template<typename T>
bool
AbstractDistMatrix<T>::Participating() const
{ return grid_->InGrid() && (CrossRank()==root_); }

template<typename T>
Int
AbstractDistMatrix<T>::RowOwner( Int i ) const
{ return (i+ColAlign()) % ColStride(); }

template<typename T>
Int
AbstractDistMatrix<T>::ColOwner( Int j ) const
{ return (j+RowAlign()) % RowStride(); }

template<typename T>
Int
AbstractDistMatrix<T>::Owner( Int i, Int j ) const
{ return RowOwner(i)+ColOwner(j)*ColStride(); }

template<typename T>
Int 
AbstractDistMatrix<T>::LocalRow( Int i ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::LocalRow");
        if( !IsLocalRow(i) )
            LogicError("Requested local index of non-local row");
    )
    return (i-ColShift()) / ColStride();
}

template<typename T>
Int
AbstractDistMatrix<T>::LocalCol( Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::LocalCol");
        if( !IsLocalCol(j) )
            LogicError("Requested local index of non-local column");
    )
    return (j-RowShift()) / RowStride();
}

template<typename T>
bool
AbstractDistMatrix<T>::IsLocalRow( Int i ) const
{ return Participating() && ((i-ColShift()) % ColStride()) == 0; }

template<typename T>
bool
AbstractDistMatrix<T>::IsLocalCol( Int j ) const
{ return Participating() && ((j-RowShift()) % RowStride()) == 0; }

template<typename T>
bool
AbstractDistMatrix<T>::IsLocal( Int i, Int j ) const
{ return IsLocalRow(i) && IsLocalCol(j); }

// Single-entry manipulation
// =========================

// Global entry manipulation
// -------------------------

template<typename T>
T
AbstractDistMatrix<T>::Get( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::Get");
        if( !grid_->InGrid() )
            LogicError("Get should only be called in-grid");
    )
    T value;
    if( CrossRank() == Root() )
    {
        const Int owner = Owner( i, j );
        if( owner == DistRank() )
        {
            const Int iLoc = (i-ColShift()) / ColStride();
            const Int jLoc = (j-RowShift()) / RowStride();
            value = GetLocal( iLoc, jLoc );
        }
        mpi::Broadcast( value, owner, DistComm() );
    }
    mpi::Broadcast( value, Root(), CrossComm() ); 
    return value;
}

template<typename T>
Base<T>
AbstractDistMatrix<T>::GetRealPart( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::GetRealPart");
        if( !grid_->InGrid() )
            LogicError("Get should only be called in-grid");
    )
    Base<T> value;
    if( CrossRank() == Root() )
    {
        const Int owner = Owner( i, j );
        if( owner == DistRank() )
        {
            const Int iLoc = (i-ColShift()) / ColStride();
            const Int jLoc = (j-RowShift()) / RowStride();
            value = GetLocalRealPart( iLoc, jLoc );
        }
        mpi::Broadcast( value, owner, DistComm() );
    }
    mpi::Broadcast( value, Root(), CrossComm() );
    return value;
}

template<typename T>
Base<T>
AbstractDistMatrix<T>::GetImagPart( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::GetImagPart");
        if( !grid_->InGrid() )
            LogicError("Get should only be called in-grid");
    )
    Base<T> value;
    if( IsComplex<T>::val )
    {
        if( CrossRank() == Root() )
        {
            const Int owner = Owner( i, j );
            if( owner == DistRank() )
            {
                const Int iLoc = (i-ColShift()) / ColStride();
                const Int jLoc = (j-RowShift()) / RowStride();
                value = GetLocalRealPart( iLoc, jLoc );
            }
            mpi::Broadcast( value, owner, DistComm() );
        }
        mpi::Broadcast( value, Root(), CrossComm() );
    }
    else
        value = 0;
    return value;
}

template<typename T>
void
AbstractDistMatrix<T>::Set( Int i, Int j, T value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Set"))
    if( CrossRank() == Root() )
    {
        const Int owner = Owner( i, j );
        if( owner == DistRank() )
        {
            const Int iLoc = (i-ColShift()) / ColStride();
            const Int jLoc = (j-RowShift()) / RowStride();
            SetLocal( iLoc, jLoc, value );
        }
    }
}

template<typename T>
void
AbstractDistMatrix<T>::SetRealPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetRealPart"))
    if( CrossRank() == Root() )
    {
        const Int owner = Owner( i, j );
        if( owner == DistRank() )
        {
            const Int iLoc = (i-ColShift()) / ColStride();
            const Int jLoc = (j-RowShift()) / RowStride();
            SetLocalRealPart( iLoc, jLoc, value );
        }
    }
}

template<typename T>
void
AbstractDistMatrix<T>::SetImagPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetImagPart"))
    if( CrossRank() == Root() )
    {
        const Int owner = Owner( i, j );
        if( owner == DistRank() )
        {
            const Int iLoc = (i-ColShift()) / ColStride();
            const Int jLoc = (j-RowShift()) / RowStride();
            SetLocalImagPart( iLoc, jLoc, value );
        }
    }
}

template<typename T>
void
AbstractDistMatrix<T>::Update( Int i, Int j, T value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Update"))
    if( CrossRank() == Root() )
    {
        const Int owner = Owner( i, j );
        if( owner == DistRank() )
        {
            const Int iLoc = (i-ColShift()) / ColStride();
            const Int jLoc = (j-RowShift()) / RowStride();
            UpdateLocal( iLoc, jLoc, value );
        }
    }
}

template<typename T>
void
AbstractDistMatrix<T>::UpdateRealPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateRealPart"))
    if( CrossRank() == Root() )
    {
        const Int owner = Owner( i, j );
        if( owner == DistRank() )
        {
            const Int iLoc = (i-ColShift()) / ColStride();
            const Int jLoc = (j-RowShift()) / RowStride();
            UpdateLocalRealPart( iLoc, jLoc, value );
        }
    }
}

template<typename T>
void
AbstractDistMatrix<T>::UpdateImagPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateImagPart"))
    if( CrossRank() == Root() )
    {
        const Int owner = Owner( i, j );
        if( owner == DistRank() )
        {
            const Int iLoc = (i-ColShift()) / ColStride();
            const Int jLoc = (j-RowShift()) / RowStride();
            UpdateLocalImagPart( iLoc, jLoc, value );
        }
    }
}

template<typename T>
void
AbstractDistMatrix<T>::MakeReal( Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeReal"))
    if( CrossRank() == Root() )
    {
        const Int owner = Owner( i, j );
        if( owner == DistRank() )
        {
            const Int iLoc = (i-ColShift()) / ColStride();
            const Int jLoc = (j-RowShift()) / RowStride();
            MakeRealLocal( iLoc, jLoc );
        }
    }
}

template<typename T>
void
AbstractDistMatrix<T>::Conjugate( Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Conjugate"))
    if( CrossRank() == Root() )
    {
        const Int owner = Owner( i, j );
        if( owner == DistRank() )
        {
            const Int iLoc = (i-ColShift()) / ColStride();
            const Int jLoc = (j-RowShift()) / RowStride();
            ConjugateLocal( iLoc, jLoc );
        }
    }
}

// Local entry manipulation
// ------------------------

template<typename T>
T
AbstractDistMatrix<T>::GetLocal( Int i, Int j ) const
{ return matrix_.Get(i,j); }

template<typename T>
Base<T>
AbstractDistMatrix<T>::GetLocalRealPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetRealPart(iLoc,jLoc); }

template<typename T>
Base<T>
AbstractDistMatrix<T>::GetLocalImagPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetImagPart(iLoc,jLoc); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Set(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocalRealPart( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.SetRealPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocalImagPart( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.SetImagPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Update(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocalRealPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.UpdateRealPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocalImagPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.UpdateImagPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::MakeRealLocal( Int iLoc, Int jLoc )
{ matrix_.MakeReal( iLoc, jLoc ); }

template<typename T>
void
AbstractDistMatrix<T>::ConjugateLocal( Int iLoc, Int jLoc )
{ matrix_.Conjugate( iLoc, jLoc ); }

// Diagonal manipulation
// =====================

template<typename T>
void
AbstractDistMatrix<T>::MakeDiagonalReal( Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeDiagonalReal"))
    const Int height = this->Height();
    const Int width = this->Width();
    const Int minDim = Min(height,width);
    const Int localWidth = this->LocalWidth();
    const Int rowShift = this->RowShift();
    const Int rowStride = this->RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        if( j < minDim && IsLocal(j,j) )
        {
            const Int iLoc = LocalRow(j);
            MakeRealLocal( iLoc, jLoc );
        }
    }
}

template<typename T>
void
AbstractDistMatrix<T>::ConjugateDiagonal( Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::ConjugateDiagonal"))
    const Int height = this->Height();
    const Int width = this->Width();
    const Int minDim = Min(height,width);
    const Int localWidth = this->LocalWidth();
    const Int rowShift = this->RowShift();
    const Int rowStride = this->RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        if( j < minDim && IsLocal(j,j) )
        {
            const Int iLoc = LocalRow(j);
            ConjugateLocal( iLoc, jLoc );
        }
    }
}

// Arbitrary submatrix manipulation
// ================================

// Global submatrix manipulation
// -----------------------------

template<typename T>
void
AbstractDistMatrix<T>::Get
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  DistMatrix<T,STAR,STAR>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Get"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    ASub.SetGrid( Grid() );
    ASub.Resize( m, n, m );
    Zeros( ASub, m, n );
    if( Participating() )
    {
        // Fill in our locally-owned entries
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = colInd[jSub];
            if( IsLocalCol(j) )
            {
                const Int jLoc = LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = rowInd[iSub];
                    if( IsLocalRow(i) )
                    {
                        const Int iLoc = LocalRow(i);
                        ASub.SetLocal( iSub, jSub, GetLocal(iLoc,jLoc) );
                    }
                }
            }
        }
        // Sum over the distribution communicator
        mpi::AllReduce( ASub.Buffer(), m*n, DistComm() ); 
    }
    // Broadcast over the cross communicator
    mpi::Broadcast( ASub.Buffer(), m*n, Root(), CrossComm() );
}

template<typename T>
void
AbstractDistMatrix<T>::GetRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  DistMatrix<BASE(T),STAR,STAR>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetRealPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    ASub.SetGrid( Grid() );
    ASub.Resize( m, n, m );
    Zeros( ASub, m, n );
    if( Participating() )
    {
        // Fill in our locally-owned entries
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = colInd[jSub];
            if( IsLocalCol(j) )
            {
                const Int jLoc = LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = rowInd[iSub];
                    if( IsLocalRow(i) )
                    {
                        const Int iLoc = LocalRow(i);
                        ASub.SetLocal
                        ( iSub, jSub, GetLocalRealPart(iLoc,jLoc) );
                    }
                }
            }
        }
        // Sum over the distribution communicator
        mpi::AllReduce( ASub.Buffer(), m*n, DistComm() ); 
    }
    // Broadcast over the cross communicator
    mpi::Broadcast( ASub.Buffer(), m*n, Root(), CrossComm() );
}

template<typename T>
void
AbstractDistMatrix<T>::GetImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  DistMatrix<BASE(T),STAR,STAR>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetImagPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    ASub.SetGrid( Grid() );
    ASub.Resize( m, n, m );
    Zeros( ASub, m, n );
    if( Participating() )
    {
        // Fill in our locally-owned entries
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = colInd[jSub];
            if( IsLocalCol(j) )
            {
                const Int jLoc = LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = rowInd[iSub];
                    if( IsLocalRow(i) )
                    {
                        const Int iLoc = LocalRow(i);
                        ASub.SetLocal
                        ( iSub, jSub, GetLocalImagPart(iLoc,jLoc) );
                    }
                }
            }
        }
        // Sum over the distribution communicator
        mpi::AllReduce( ASub.Buffer(), m*n, DistComm() ); 
    }
    // Broadcast over the cross communicator
    mpi::Broadcast( ASub.Buffer(), m*n, Root(), CrossComm() );
}

template<typename T>
DistMatrix<T,STAR,STAR>
AbstractDistMatrix<T>::Get
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Get"))
    DistMatrix<T,STAR,STAR> ASub( Grid() );
    Get( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
DistMatrix<BASE(T),STAR,STAR>
AbstractDistMatrix<T>::GetRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetRealPart"))
    DistMatrix<Base<T>,STAR,STAR> ASub( Grid() );
    GetRealPart( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
DistMatrix<BASE(T),STAR,STAR>
AbstractDistMatrix<T>::GetImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetImagPart"))
    DistMatrix<Base<T>,STAR,STAR> ASub( Grid() );
    GetImagPart( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
void 
AbstractDistMatrix<T>::Set
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const DistMatrix<T,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Set"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    if( Participating() )
    {
        // Fill in our locally-owned entries
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = colInd[jSub];
            if( IsLocalCol(j) )
            {
                const Int jLoc = LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = rowInd[iSub];
                    if( IsLocalRow(i) )
                    {
                        const Int iLoc = LocalRow(i);
                        SetLocal( iLoc, jLoc, ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

template<typename T>
void 
AbstractDistMatrix<T>::SetRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const DistMatrix<BASE(T),STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetRealPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    if( Participating() )
    {
        // Fill in our locally-owned entries
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = colInd[jSub];
            if( IsLocalCol(j) )
            {
                const Int jLoc = LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = rowInd[iSub];
                    if( IsLocalRow(i) )
                    {
                        const Int iLoc = LocalRow(i);
                        SetLocalRealPart
                        ( iLoc, jLoc, ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

template<typename T>
void 
AbstractDistMatrix<T>::SetImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const DistMatrix<BASE(T),STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetImagPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    if( Participating() )
    {
        // Fill in our locally-owned entries
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = colInd[jSub];
            if( IsLocalCol(j) )
            {
                const Int jLoc = LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = rowInd[iSub];
                    if( IsLocalRow(i) )
                    {
                        const Int iLoc = LocalRow(i);
                        SetLocalImagPart
                        ( iLoc, jLoc, ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

template<typename T>
void 
AbstractDistMatrix<T>::Update
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  T alpha, const DistMatrix<T,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Update"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    if( Participating() )
    {
        // Modify our locally-owned entries
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = colInd[jSub];
            if( IsLocalCol(j) )
            {
                const Int jLoc = LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = rowInd[iSub];
                    if( IsLocalRow(i) )
                    {
                        const Int iLoc = LocalRow(i);
                        UpdateLocal
                        ( iLoc, jLoc, alpha*ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

template<typename T>
void 
AbstractDistMatrix<T>::UpdateRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  BASE(T) alpha, const DistMatrix<BASE(T),STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateRealPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    if( Participating() )
    {
        // Modify our locally-owned entries
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = colInd[jSub];
            if( IsLocalCol(j) )
            {
                const Int jLoc = LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = rowInd[iSub];
                    if( IsLocalRow(i) )
                    {
                        const Int iLoc = LocalRow(i);
                        UpdateLocalRealPart
                        ( iLoc, jLoc, alpha*ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

template<typename T>
void 
AbstractDistMatrix<T>::UpdateImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  BASE(T) alpha, const DistMatrix<BASE(T),STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateImagPart"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    if( Participating() )
    {
        // Modify our locally-owned entries
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = colInd[jSub];
            if( IsLocalCol(j) )
            {
                const Int jLoc = LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = rowInd[iSub];
                    if( IsLocalRow(i) )
                    {
                        const Int iLoc = LocalRow(i);
                        UpdateLocalImagPart
                        ( iLoc, jLoc, alpha*ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

template<typename T>
void 
AbstractDistMatrix<T>::MakeReal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeReal"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    if( Participating() )
    {
        // Modify the locally-owned entries 
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = colInd[jSub];
            if( IsLocalCol(j) )
            {
                const Int jLoc = LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = rowInd[iSub];
                    if( IsLocalRow(i) )
                    {
                        const Int iLoc = LocalRow(i);
                        MakeRealLocal( iLoc, jLoc );
                    }
                }
            }
        }
    }
}

template<typename T>
void 
AbstractDistMatrix<T>::Conjugate
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Conjugate"))
    const Int m = rowInd.size();
    const Int n = colInd.size();
    if( Participating() )
    {
        // Modify the locally-owned entries 
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = colInd[jSub];
            if( IsLocalCol(j) )
            {
                const Int jLoc = LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = rowInd[iSub];
                    if( IsLocalRow(i) )
                    {
                        const Int iLoc = LocalRow(i);
                        ConjugateLocal( iLoc, jLoc );
                    }
                }
            }
        }
    }
}

// Local submatrix manipulation
// ----------------------------

template<typename T>
void
AbstractDistMatrix<T>::GetLocal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  elem::Matrix<T>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetLocal"))
    LockedMatrix().Get( rowInd, colInd, ASub );
}

template<typename T>
void
AbstractDistMatrix<T>::GetLocalRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  elem::Matrix<BASE(T)>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetLocalRealPart"))
    LockedMatrix().GetRealPart( rowInd, colInd, ASub );
}

template<typename T>
void
AbstractDistMatrix<T>::GetLocalImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  elem::Matrix<BASE(T)>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetLocalImagPart"))
    LockedMatrix().GetImagPart( rowInd, colInd, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::SetLocal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const elem::Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetLocal"))
    Matrix().Set( rowInd, colInd, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::SetLocalRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const elem::Matrix<BASE(T)>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetLocalRealPart"))
    Matrix().SetRealPart( rowInd, colInd, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::SetLocalImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const elem::Matrix<BASE(T)>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetLocalImagPart"))
    Matrix().SetImagPart( rowInd, colInd, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::UpdateLocal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  T alpha, const elem::Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateLocal"))
    Matrix().Update( rowInd, colInd, alpha, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::UpdateLocalRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  BASE(T) alpha, const elem::Matrix<BASE(T)>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateLocalRealPart"))
    Matrix().UpdateRealPart( rowInd, colInd, alpha, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::UpdateLocalImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  BASE(T) alpha, const elem::Matrix<BASE(T)>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateLocalImagPart"))
    Matrix().UpdateImagPart( rowInd, colInd, alpha, ASub );
}

template<typename T>
void
AbstractDistMatrix<T>::MakeRealLocal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeRealLocal"))
    Matrix().MakeReal( rowInd, colInd );
}

template<typename T>
void
AbstractDistMatrix<T>::ConjugateLocal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::ConjugateLocal"))
    Matrix().Conjugate( rowInd, colInd );
}

// Combined realignment and resize
// ===============================

template<typename T>
void
AbstractDistMatrix<T>::AlignAndResize
( Int colAlign, Int rowAlign, Int height, Int width, bool force )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignAndResize"))
    if( !Viewing() )
    {
        if( force || !ColConstrained() )
        {
            colAlign_ = colAlign;
            SetColShift(); 
        }
        if( force || !RowConstrained() )
        {
            rowAlign_ = rowAlign;
            SetRowShift();
        }
    }
    if( force && (colAlign_ != colAlign || rowAlign_ != rowAlign) )
        LogicError("Could not set alignments"); 
    Resize( height, width );
}

template<typename T>
void
AbstractDistMatrix<T>::AlignColsAndResize
( Int colAlign, Int height, Int width, bool force )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignColsAndResize"))
    if( !Viewing() && (force || !ColConstrained()) )
    {
        colAlign_ = colAlign;
        SetColShift(); 
    }
    if( force && colAlign_ != colAlign )
        LogicError("Could not set col alignment");
    Resize( height, width );
}

template<typename T>
void
AbstractDistMatrix<T>::AlignRowsAndResize
( Int rowAlign, Int height, Int width, bool force )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignRowsAndResize"))
    if( !Viewing() && (force || !RowConstrained()) )
    {
        rowAlign_ = rowAlign;
        SetRowShift(); 
    }
    if( force && rowAlign_ != rowAlign )
        LogicError("Could not set row alignment");
    Resize( height, width );
}

// Assertions
// ==========

template<typename T>
void
AbstractDistMatrix<T>::ComplainIfReal() const
{ 
    if( !IsComplex<T>::val )
        LogicError("Called complex-only routine with real data");
}

template<typename T>
void
AbstractDistMatrix<T>::AssertNotLocked() const
{
    if( Locked() )
        LogicError("Assertion that matrix not be a locked view failed");
}

template<typename T>
void
AbstractDistMatrix<T>::AssertNotStoringData() const
{
    if( matrix_.MemorySize() > 0 )
        LogicError("Assertion that matrix not be storing data failed");
}

template<typename T>
void
AbstractDistMatrix<T>::AssertValidEntry( Int i, Int j ) const
{
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
        LogicError
        ("Entry (",i,",",j,") is out of bounds of ",Height(),
         " x ",Width()," matrix");
}

template<typename T>
void
AbstractDistMatrix<T>::AssertValidSubmatrix
( Int i, Int j, Int height, Int width ) const
{
    if( i < 0 || j < 0 )
        LogicError("Indices of submatrix were negative");
    if( height < 0 || width < 0 )
        LogicError("Dimensions of submatrix were negative");
    if( (i+height) > Height() || (j+width) > Width() )
        LogicError
        ("Submatrix is out of bounds: accessing up to (",i+height-1,
         ",",j+width-1,") of ",Height()," x ",Width()," matrix");
}

template<typename T> 
void
AbstractDistMatrix<T>::AssertSameGrid( const elem::Grid& grid ) const
{
    if( Grid() != grid )
        LogicError("Assertion that grids match failed");
}

template<typename T> 
void
AbstractDistMatrix<T>::AssertSameSize( Int height, Int width ) const
{
    if( Height() != height || Width() != width )
        LogicError("Assertion that matrices be the same size failed");
}

// Private section
// ###############

// Construct using a particular process grid
// =========================================

template<typename T>
AbstractDistMatrix<T>::AbstractDistMatrix( const elem::Grid& grid )
: viewType_(OWNER),
  height_(0), width_(0), 
  auxMemory_(), 
  matrix_(0,0,true), 
  colConstrained_(false), rowConstrained_(false),
  colAlign_(0), rowAlign_(0),
  colShift_(0), rowShift_(0), 
  root_(0),
  grid_(&grid)
{ }

// Exchange metadata with another matrix
// =====================================

template<typename T>
void 
AbstractDistMatrix<T>::ShallowSwap( ADM<T>& A )
{
    matrix_.ShallowSwap( A.matrix_ );
    auxMemory_.ShallowSwap( A.auxMemory_ );
    std::swap( viewType_, A.viewType_ );
    std::swap( height_ , A.height_ );
    std::swap( width_, A.width_ );
    std::swap( colConstrained_, A.colConstrained_ );
    std::swap( rowConstrained_, A.rowConstrained_ );
    std::swap( colAlign_, A.colAlign_ );
    std::swap( rowAlign_, A.rowAlign_ );
    std::swap( colShift_, A.colShift_ );
    std::swap( rowShift_, A.rowShift_ );
    std::swap( root_, A.root_ );
    std::swap( grid_, A.grid_ );
}

// Modify the distribution metadata
// ================================

template<typename T>
void
AbstractDistMatrix<T>::SetShifts()
{
    if( Participating() )
    {
        colShift_ = Shift(ColRank(),colAlign_,ColStride());
        rowShift_ = Shift(RowRank(),rowAlign_,RowStride());
    }
    else
    {
        colShift_ = 0;
        rowShift_ = 0;
    }
}

template<typename T>
void
AbstractDistMatrix<T>::SetColShift()
{
    if( Participating() )
        colShift_ = Shift(ColRank(),colAlign_,ColStride());
    else
        colShift_ = 0;
}

template<typename T>
void
AbstractDistMatrix<T>::SetRowShift()
{
    if( Participating() )
        rowShift_ = Shift(RowRank(),rowAlign_,RowStride());
    else
        rowShift_ = 0;
}

// Outside of class
// ----------------

template<typename T> 
void
AssertConforming1x2( const ADM<T>& AL, const ADM<T>& AR )
{
    if( AL.Height() != AR.Height() )    
        LogicError
        ("1x2 not conformant:\n",
         DimsString(AL,"Left"),"\n",DimsString(AR,"Right"));
    if( AL.ColAlign() != AR.ColAlign() )
        LogicError("1x2 is misaligned");
}

template<typename T> 
void
AssertConforming2x1( const ADM<T>& AT, const ADM<T>& AB )
{
    if( AT.Width() != AB.Width() )
        LogicError
        ("2x1 is not conformant:\n",
         DimsString(AT,"Top"),"\n",DimsString(AB,"Bottom"));
    if( AT.RowAlign() != AB.RowAlign() )
        LogicError("2x1 is not aligned");
}

template<typename T> 
void
AssertConforming2x2
( const ADM<T>& ATL, const ADM<T>& ATR, 
  const ADM<T>& ABL, const ADM<T>& ABR ) 
{
    if( ATL.Width() != ABL.Width() || ATR.Width() != ABR.Width() ||
        ATL.Height() != ATR.Height() || ABL.Height() != ABR.Height() )
        LogicError
        ("2x2 is not conformant:\n",
         DimsString(ATL,"TL"),"\n",DimsString(ATR,"TR"),"\n",
         DimsString(ABL,"BL"),"\n",DimsString(ABR,"BR"));
    if( ATL.ColAlign() != ATR.ColAlign() ||
        ABL.ColAlign() != ABR.ColAlign() ||
        ATL.RowAlign() != ABL.RowAlign() ||
        ATR.RowAlign() != ABR.RowAlign() )
        LogicError("2x2 set of matrices must aligned to combine");
}

// Instantiations for {Int,Real,Complex<Real>} for each Real in {float,double}
// ###########################################################################

#ifndef RELEASE
 #define PROTO(T) \
  template class AbstractDistMatrix<T>;\
  template void AssertConforming1x2\
  ( const ADM<T>& AL, const ADM<T>& AR );\
  template void AssertConforming2x1\
  ( const ADM<T>& AT, const ADM<T>& AB );\
  template void AssertConforming2x2\
  ( const ADM<T>& ATL, const ADM<T>& ATR,\
    const ADM<T>& ABL, const ADM<T>& ABR )
#else
 #define PROTO(T) template class AbstractDistMatrix<T>
#endif
 
#ifndef DISABLE_COMPLEX
 #ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(float);
  PROTO(double);
  PROTO(Complex<float>);
  PROTO(Complex<double>);
 #else // ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(double);
  PROTO(Complex<double>);
 #endif // ifndef DISABLE_FLOAT
#else // ifndef DISABLE_COMPLEX
 #ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(float);
  PROTO(double);
 #else // ifndef DISABLE_FLOAT
  PROTO(Int);
  PROTO(double);
 #endif // ifndef DISABLE_FLOAT
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
