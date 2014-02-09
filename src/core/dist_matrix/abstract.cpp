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

template<typename T,Dist U,Dist V>
using ADM = AbstractDistMatrix<T,U,V>;

// NOTE: It seems that member functions cannot be defined using a 
//       fully-specified template alias, e.g., ADM<T,U,V>::AbstractDistMatrix(),
//       but DM<T> is okay if it is only partially specified, e.g., 
//       DM<T> = DistMatrix<T,MC,MR> and DM<T>::DistMatrix()

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T,Dist U,Dist V>
AbstractDistMatrix<T,U,V>::AbstractDistMatrix( ADM<T,U,V>&& A )
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

template<typename T,Dist U,Dist V>
AbstractDistMatrix<T,U,V>::~AbstractDistMatrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T,Dist U,Dist V>
ADM<T,U,V>& 
AbstractDistMatrix<T,U,V>::operator=( ADM<T,U,V>&& A )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::Empty()
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::EmptyData()
{
    matrix_.Empty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetGrid( const elem::Grid& grid )
{
    if( grid_ != &grid )
    {
        Empty();
        grid_ = &grid; 
        SetShifts();
    }
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::Resize( Int height, Int width )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::Resize( Int height, Int width, Int ldim )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::MakeConsistent()
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::Align( Int colAlign, Int rowAlign )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AlignCols( Int colAlign )
{ 
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignCols"))
    if( colAlign_ != colAlign )
        EmptyData();
    colConstrained_ = true;
    colAlign_ = colAlign;
    SetShifts();
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AlignRows( Int rowAlign )
{ 
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignRows"))
    if( rowAlign_ != rowAlign )
        EmptyData();
    rowConstrained_ = true;
    rowAlign_ = rowAlign;
    SetShifts();
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::FreeAlignments() 
{ 
    colConstrained_ = false;
    rowConstrained_ = false;
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetRoot( Int root )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AlignWith( const elem::DistData& data )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AlignColsWith( const elem::DistData& data )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AlignRowsWith( const elem::DistData& data )
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

// Basic queries
// =============

// Global matrix information
// -------------------------

template<typename T,Dist U,Dist V>
Int AbstractDistMatrix<T,U,V>::Height() const { return height_; }
template<typename T,Dist U,Dist V>
Int AbstractDistMatrix<T,U,V>::Width() const { return width_; }

template<typename T,Dist U,Dist V>
Int AbstractDistMatrix<T,U,V>::DiagonalLength( Int offset ) const
{ return elem::DiagonalLength(height_,width_,offset); }

template<typename T,Dist U,Dist V>
bool AbstractDistMatrix<T,U,V>::Viewing() const { return IsViewing( viewType_ ); }
template<typename T,Dist U,Dist V>
bool AbstractDistMatrix<T,U,V>::Locked() const { return IsLocked( viewType_ ); }

// Local matrix information
// ------------------------

template<typename T,Dist U,Dist V>
Int AbstractDistMatrix<T,U,V>::LocalHeight() const { return matrix_.Height(); }
template<typename T,Dist U,Dist V>
Int AbstractDistMatrix<T,U,V>::LocalWidth() const { return matrix_.Width(); }
template<typename T,Dist U,Dist V>
Int AbstractDistMatrix<T,U,V>::LDim() const { return matrix_.LDim(); }

template<typename T,Dist U,Dist V>
elem::Matrix<T>& 
AbstractDistMatrix<T,U,V>::Matrix() { return matrix_; }
template<typename T,Dist U,Dist V>
const elem::Matrix<T>& 
AbstractDistMatrix<T,U,V>::LockedMatrix() const { return matrix_; }

template<typename T,Dist U,Dist V>
size_t
AbstractDistMatrix<T,U,V>::AllocatedMemory() const { return matrix_.MemorySize(); }

template<typename T,Dist U,Dist V>
T*
AbstractDistMatrix<T,U,V>::Buffer() { return matrix_.Buffer(); }

template<typename T,Dist U,Dist V>
T*
AbstractDistMatrix<T,U,V>::Buffer( Int iLoc, Int jLoc )
{ return matrix_.Buffer(iLoc,jLoc); }

template<typename T,Dist U,Dist V>
const T*
AbstractDistMatrix<T,U,V>::LockedBuffer() const
{ return matrix_.LockedBuffer(); }

template<typename T,Dist U,Dist V>
const T*
AbstractDistMatrix<T,U,V>::LockedBuffer( Int iLoc, Int jLoc ) const
{ return matrix_.LockedBuffer(iLoc,jLoc); }

// Distribution information
// ------------------------

template<typename T,Dist U,Dist V>
const elem::Grid& AbstractDistMatrix<T,U,V>::Grid() const { return *grid_; }

template<typename T,Dist U,Dist V>
bool AbstractDistMatrix<T,U,V>::ColConstrained() const { return colConstrained_; }
template<typename T,Dist U,Dist V>
bool AbstractDistMatrix<T,U,V>::RowConstrained() const { return rowConstrained_; }

template<typename T,Dist U,Dist V>
Int AbstractDistMatrix<T,U,V>::ColAlign() const { return colAlign_; }
template<typename T,Dist U,Dist V>
Int AbstractDistMatrix<T,U,V>::RowAlign() const { return rowAlign_; }

template<typename T,Dist U,Dist V>
Int AbstractDistMatrix<T,U,V>::ColShift() const { return colShift_; }
template<typename T,Dist U,Dist V>
Int AbstractDistMatrix<T,U,V>::RowShift() const { return rowShift_; }

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::ColRank() const
{ 
    if( grid_->InGrid() )
        return mpi::CommRank(ColComm());
    else
        return mpi::UNDEFINED;
}

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::RowRank() const
{
    if( grid_->InGrid() )
        return mpi::CommRank(RowComm());
    else
        return mpi::UNDEFINED;
}

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::DistRank() const
{ return mpi::CommRank(DistComm()); }

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::CrossRank() const
{ return mpi::CommRank(CrossComm()); }

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::RedundantRank() const
{ return mpi::CommRank(RedundantComm()); }

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::DistSize() const
{ return mpi::CommSize(DistComm()); }

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::CrossSize() const
{ return mpi::CommSize(CrossComm()); }

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::RedundantSize() const
{ return mpi::CommSize(RedundantComm()); }

template<typename T,Dist U,Dist V>
Int AbstractDistMatrix<T,U,V>::Root() const { return root_; }

template<typename T,Dist U,Dist V>
bool
AbstractDistMatrix<T,U,V>::Participating() const
{ return grid_->InGrid() && (CrossRank()==root_); }

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::RowOwner( Int i ) const
{ return (i+ColAlign()) % ColStride(); }

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::ColOwner( Int j ) const
{ return (j+RowAlign()) % RowStride(); }

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::Owner( Int i, Int j ) const
{ return RowOwner(i)+ColOwner(j)*ColStride(); }

template<typename T,Dist U,Dist V>
Int 
AbstractDistMatrix<T,U,V>::LocalRow( Int i ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::LocalRow");
        if( !IsLocalRow(i) )
            LogicError("Requested local index of non-local row");
    )
    return (i-ColShift()) / ColStride();
}

template<typename T,Dist U,Dist V>
Int
AbstractDistMatrix<T,U,V>::LocalCol( Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::LocalCol");
        if( !IsLocalCol(j) )
            LogicError("Requested local index of non-local column");
    )
    return (j-RowShift()) / RowStride();
}

template<typename T,Dist U,Dist V>
bool
AbstractDistMatrix<T,U,V>::IsLocalRow( Int i ) const
{ return Participating() && ((i-ColShift()) % ColStride()) == 0; }

template<typename T,Dist U,Dist V>
bool
AbstractDistMatrix<T,U,V>::IsLocalCol( Int j ) const
{ return Participating() && ((j-RowShift()) % RowStride()) == 0; }

template<typename T,Dist U,Dist V>
bool
AbstractDistMatrix<T,U,V>::IsLocal( Int i, Int j ) const
{ return IsLocalRow(i) && IsLocalCol(j); }

// Single-entry manipulation
// =========================

// Global entry manipulation
// -------------------------

template<typename T,Dist U,Dist V>
T
AbstractDistMatrix<T,U,V>::Get( Int i, Int j ) const
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

template<typename T,Dist U,Dist V>
Base<T>
AbstractDistMatrix<T,U,V>::GetRealPart( Int i, Int j ) const
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

template<typename T,Dist U,Dist V>
Base<T>
AbstractDistMatrix<T,U,V>::GetImagPart( Int i, Int j ) const
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::Set( Int i, Int j, T value )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetRealPart( Int i, Int j, Base<T> value )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetImagPart( Int i, Int j, Base<T> value )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::Update( Int i, Int j, T value )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::UpdateRealPart( Int i, Int j, Base<T> value )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::UpdateImagPart( Int i, Int j, Base<T> value )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::MakeReal( Int i, Int j )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::Conjugate( Int i, Int j )
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

template<typename T,Dist U,Dist V>
T
AbstractDistMatrix<T,U,V>::GetLocal( Int i, Int j ) const
{ return matrix_.Get(i,j); }

template<typename T,Dist U,Dist V>
Base<T>
AbstractDistMatrix<T,U,V>::GetLocalRealPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetRealPart(iLoc,jLoc); }

template<typename T,Dist U,Dist V>
Base<T>
AbstractDistMatrix<T,U,V>::GetLocalImagPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetImagPart(iLoc,jLoc); }

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Set(iLoc,jLoc,alpha); }

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetLocalRealPart( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.SetRealPart(iLoc,jLoc,alpha); }

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetLocalImagPart( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.SetImagPart(iLoc,jLoc,alpha); }

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::UpdateLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Update(iLoc,jLoc,alpha); }

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::UpdateLocalRealPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.UpdateRealPart(iLoc,jLoc,alpha); }

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::UpdateLocalImagPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.UpdateImagPart(iLoc,jLoc,alpha); }

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::MakeRealLocal( Int iLoc, Int jLoc )
{ matrix_.MakeReal( iLoc, jLoc ); }

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::ConjugateLocal( Int iLoc, Int jLoc )
{ matrix_.Conjugate( iLoc, jLoc ); }

// Diagonal manipulation
// =====================
template<typename T,Dist U,Dist V>
template<typename S>
bool
AbstractDistMatrix<T,U,V>::DiagonalAligned
( const DistMatrix<S,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::DiagonalAligned"))
    if( this->Grid() != d.Grid() )
        return false;

    if( U == MC && V == MR )
    {
        // The result is an [MD,* ]
        const Int diagRow = d.ColAlign()            % this->ColStride();
        const Int diagCol = (d.ColAlign()+d.Root()) % this->RowStride();
        if( offset >= 0 )
        {
            const Int procRow = this->ColAlign();
            const Int procCol = (this->RowAlign()+offset) % this->RowStride();
            return procRow==diagRow && procCol==diagCol;
        }
        else
        {
            const Int procRow = (this->ColAlign()-offset) % this->ColStride(); 
            const Int procCol = this->RowAlign();
            return procRow==diagRow && procCol==diagCol;
        }
    }
    else if( U == MR && V == MC )
    {
        // The result is an [MD,* ]
        const Int diagRow = d.ColAlign()            % this->ColStride();
        const Int diagCol = (d.ColAlign()+d.Root()) % this->RowStride();
        if( offset >= 0 )
        {
            const Int procCol = this->ColAlign();
            const Int procRow = (this->RowAlign()+offset) % this->RowStride();
            return procRow==diagRow && procCol==diagCol;
        }
        else
        {
            const Int procCol = (this->ColAlign()-offset) % this->ColStride();
            const Int procRow = this->RowAlign();
            return procRow==diagRow && procCol==diagCol;
        }
    }
    else if( U == STAR )
    {
        // The result is a [V,* ]
        if( offset >= 0 )
            return this->Root()==d.Root() && 
                   ((this->RowAlign()+offset)%this->RowStride())==d.ColAlign();
        else
            return this->Root()==d.Root() && this->RowAlign()==d.ColAlign();
    }
    else
    {
        // The result is a [U,V], where V is either STAR or CIRC
        if( offset >= 0 )
            return this->Root()==d.Root() && this->ColAlign() == d.ColAlign();
        else
            return this->Root()==d.Root() && 
                   ((this->ColAlign()-offset)%this->ColStride())==d.ColAlign();
    }
}

template<typename T,Dist U,Dist V>
template<typename S>
void
AbstractDistMatrix<T,U,V>::ForceDiagonalAlign
( DistMatrix<S,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::ForceDiagonalAlign"))
    const elem::Grid& grid = this->Grid();
    d.SetGrid( grid );

    if( U == MC && V == MR )
    {
        // Result is an [MD,* ]
        Int owner;
        if( offset >= 0 )
        {
            const Int procRow = this->ColAlign();
            const Int procCol = (this->RowAlign()+offset) % this->RowStride();
            owner = procRow + this->ColStride()*procCol;
        }
        else
        {
            const Int procRow = (this->ColAlign()-offset) % this->ColStride();
            const Int procCol = this->RowAlign();
            owner = procRow + this->ColStride()*procCol;
        }
        d.SetRoot( grid.DiagPath(owner) );
        d.AlignCols( grid.DiagPathRank(owner) );
    }
    else if( U == MR && V == MC )
    {
        // Result is an [MD,* ]
        Int owner;
        if( offset >= 0 )
        {
            const Int procCol = this->ColAlign();
            const Int procRow = (this->RowAlign()+offset) % this->RowStride();
            owner = procRow + this->ColStride()*procCol;
        }
        else
        {
            const Int procCol = (this->ColAlign()-offset) % this->ColStride();
            const Int procRow = this->RowAlign();
            owner = procRow + this->ColStride()*procCol;
        }
        d.SetRoot( grid.DiagPath(owner) );
        d.AlignCols( grid.DiagPathRank(owner) );
    }
    else if( U == STAR )
    {
        // Result is a [V,* ] 
        Int colAlign;
        if( offset >= 0 )
            colAlign = (this->RowAlign()+offset) % this->RowStride();
        else
            colAlign = this->RowAlign();
        d.SetRoot( this->Root() );
        d.AlignCols( colAlign );
    }
    else
    {
        // Result is a [U,V], where V is either STAR or CIRC
        Int colAlign;
        if( offset >= 0 )
            colAlign = this->ColAlign();
        else
            colAlign = (this->ColAlign()-offset) % this->ColStride();
        d.SetRoot( this->Root() );
        d.AlignCols( colAlign );
    }
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::GetDiagonal
( DistMatrix<T,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::GetRealPartOfDiagonal
( DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetRealPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = RealPart(beta); } );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::GetImagPartOfDiagonal
( DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetImagPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = ImagPart(beta); } );
}

template<typename T,Dist U,Dist V>
auto
AbstractDistMatrix<T,U,V>::GetDiagonal( Int offset ) const
-> DistMatrix<T,UDiag,VDiag>
{
    DistMatrix<T,UDiag,VDiag> d( this->Grid() );
    GetDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
AbstractDistMatrix<T,U,V>::GetRealPartOfDiagonal( Int offset ) const
-> DistMatrix<Base<T>,UDiag,VDiag>
{
    DistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
AbstractDistMatrix<T,U,V>::GetImagPartOfDiagonal( Int offset ) const
-> DistMatrix<Base<T>,UDiag,VDiag>
{
    DistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetDiagonal
( const DistMatrix<T,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetRealPartOfDiagonal
( const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetRealPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetImagPartOfDiagonal
( const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { elem::SetImagPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::UpdateDiagonal
( T gamma, const DistMatrix<T,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, [gamma]( T& alpha, T beta ) { alpha += gamma*beta; } );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::UpdateRealPartOfDiagonal
( Base<T> gamma, const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      [gamma]( T& alpha, Base<T> beta ) 
      { elem::UpdateRealPart(alpha,gamma*beta); } );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::UpdateImagPartOfDiagonal
( Base<T> gamma, const DistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      [gamma]( T& alpha, Base<T> beta ) 
      { elem::UpdateImagPart(alpha,gamma*beta); } );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::MakeDiagonalReal( Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeDiagonalReal"))
    // TODO
    LogicError("This routine not yet written");
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::ConjugateDiagonal( Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::ConjugateDiagonal"))
    // TODO
    LogicError("This routine not yet written");
}

// Arbitrary submatrix manipulation
// ================================

// Global submatrix manipulation
// -----------------------------

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::Get
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::GetRealPart
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::GetImagPart
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

template<typename T,Dist U,Dist V>
DistMatrix<T,STAR,STAR>
AbstractDistMatrix<T,U,V>::Get
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Get"))
    DistMatrix<T,STAR,STAR> ASub( Grid() );
    Get( rowInd, colInd, ASub );
    return ASub;
}

template<typename T,Dist U,Dist V>
DistMatrix<BASE(T),STAR,STAR>
AbstractDistMatrix<T,U,V>::GetRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetRealPart"))
    DistMatrix<Base<T>,STAR,STAR> ASub( Grid() );
    GetRealPart( rowInd, colInd, ASub );
    return ASub;
}

template<typename T,Dist U,Dist V>
DistMatrix<BASE(T),STAR,STAR>
AbstractDistMatrix<T,U,V>::GetImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetImagPart"))
    DistMatrix<Base<T>,STAR,STAR> ASub( Grid() );
    GetImagPart( rowInd, colInd, ASub );
    return ASub;
}

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::Set
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

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::SetRealPart
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

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::SetImagPart
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

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::Update
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

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::UpdateRealPart
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

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::UpdateImagPart
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

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::MakeReal
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

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::Conjugate
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::GetLocal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  elem::Matrix<T>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetLocal"))
    LockedMatrix().Get( rowInd, colInd, ASub );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::GetLocalRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  elem::Matrix<BASE(T)>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetLocalRealPart"))
    LockedMatrix().GetRealPart( rowInd, colInd, ASub );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::GetLocalImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  elem::Matrix<BASE(T)>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetLocalImagPart"))
    LockedMatrix().GetImagPart( rowInd, colInd, ASub );
}

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::SetLocal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const elem::Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetLocal"))
    Matrix().Set( rowInd, colInd, ASub );
}

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::SetLocalRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const elem::Matrix<BASE(T)>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetLocalRealPart"))
    Matrix().SetRealPart( rowInd, colInd, ASub );
}

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::SetLocalImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const elem::Matrix<BASE(T)>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetLocalImagPart"))
    Matrix().SetImagPart( rowInd, colInd, ASub );
}

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::UpdateLocal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  T alpha, const elem::Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateLocal"))
    Matrix().Update( rowInd, colInd, alpha, ASub );
}

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::UpdateLocalRealPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  BASE(T) alpha, const elem::Matrix<BASE(T)>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateLocalRealPart"))
    Matrix().UpdateRealPart( rowInd, colInd, alpha, ASub );
}

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::UpdateLocalImagPart
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  BASE(T) alpha, const elem::Matrix<BASE(T)>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateLocalImagPart"))
    Matrix().UpdateImagPart( rowInd, colInd, alpha, ASub );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::MakeRealLocal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeRealLocal"))
    Matrix().MakeReal( rowInd, colInd );
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::ConjugateLocal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::ConjugateLocal"))
    Matrix().Conjugate( rowInd, colInd );
}

// Private section
// ###############

// Construct using a particular process grid
// =========================================

template<typename T,Dist U,Dist V>
AbstractDistMatrix<T,U,V>::AbstractDistMatrix( const elem::Grid& grid )
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

template<typename T,Dist U,Dist V>
void 
AbstractDistMatrix<T,U,V>::ShallowSwap( ADM<T,U,V>& A )
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetShifts()
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetColShift()
{
    if( Participating() )
        colShift_ = Shift(ColRank(),colAlign_,ColStride());
    else
        colShift_ = 0;
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::SetRowShift()
{
    if( Participating() )
        rowShift_ = Shift(RowRank(),rowAlign_,RowStride());
    else
        rowShift_ = 0;
}

// Combined realignment and resize
// ===============================

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AlignAndResize
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AlignColsAndResize
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AlignRowsAndResize
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

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::ComplainIfReal() const
{ 
    if( !IsComplex<T>::val )
        LogicError("Called complex-only routine with real data");
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AssertNotLocked() const
{
    if( Locked() )
        LogicError("Assertion that matrix not be a locked view failed");
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AssertNotStoringData() const
{
    if( matrix_.MemorySize() > 0 )
        LogicError("Assertion that matrix not be storing data failed");
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AssertValidEntry( Int i, Int j ) const
{
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
        LogicError
        ("Entry (",i,",",j,") is out of bounds of ",Height(),
         " x ",Width()," matrix");
}

template<typename T,Dist U,Dist V>
void
AbstractDistMatrix<T,U,V>::AssertValidSubmatrix
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

template<typename T,Dist U,Dist V> 
void
AbstractDistMatrix<T,U,V>::AssertSameGrid( const elem::Grid& grid ) const
{
    if( Grid() != grid )
        LogicError("Assertion that grids match failed");
}

template<typename T,Dist U,Dist V> 
void
AbstractDistMatrix<T,U,V>::AssertSameSize( Int height, Int width ) const
{
    if( Height() != height || Width() != width )
        LogicError("Assertion that matrices be the same size failed");
}

// Helper functions
// ================
template<typename T,Dist U,Dist V>
template<typename S,class Function>
void
AbstractDistMatrix<T,U,V>::GetDiagonalHelper
( DistMatrix<S,UDiag,VDiag>& d, Int offset, Function func ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetDiagonalHelper"))
    d.SetGrid( this->Grid() );
    this->ForceDiagonalAlign( d, offset );
    d.Resize( this->DiagonalLength(offset), 1 );
    if( !d.Participating() )
        return;

    const Int diagShift = d.ColShift();
    const Int diagStride = d.ColStride();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int iLocStart = (iStart-this->ColShift()) / colStride;
    const Int jLocStart = (jStart-this->RowShift()) / rowStride;

    const Int localDiagLength = d.LocalHeight();
    S* dBuf = d.Buffer();
    const T* buffer = this->LockedBuffer();
    const Int ldim = this->LDim();

    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(diagStride/colStride);
        const Int jLoc = jLocStart + k*(diagStride/rowStride);
        func( dBuf[k], buffer[iLoc+jLoc*ldim] );
    }
}

template<typename T,Dist U,Dist V>
template<typename S,class Function>
void
AbstractDistMatrix<T,U,V>::SetDiagonalHelper
( const DistMatrix<S,UDiag,VDiag>& d, Int offset, Function func ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::SetDiagonalHelper");
        if( !this->DiagonalAligned( d, offset ) )
            LogicError("Invalid diagonal alignment");
    )
    if( !d.Participating() )
        return;

    const Int diagShift = d.ColShift();
    const Int diagStride = d.ColStride();
    const Int iStart = ( offset>=0 ? diagShift        : diagShift-offset );
    const Int jStart = ( offset>=0 ? diagShift+offset : diagShift        );

    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int iLocStart = (iStart-this->ColShift()) / colStride;
    const Int jLocStart = (jStart-this->RowShift()) / rowStride;

    const Int localDiagLength = d.LocalHeight();
    const S* dBuf = d.LockedBuffer();
    T* buffer = this->Buffer();
    const Int ldim = this->LDim();

    PARALLEL_FOR
    for( Int k=0; k<localDiagLength; ++k )
    {
        const Int iLoc = iLocStart + k*(diagStride/colStride);
        const Int jLoc = jLocStart + k*(diagStride/rowStride);
        func( buffer[iLoc+jLoc*ldim], dBuf[k] );
    }
}

// Outside of class
// ----------------

template<typename T,Dist U,Dist V> 
void
AssertConforming1x2( const ADM<T,U,V>& AL, const ADM<T,U,V>& AR )
{
    if( AL.Height() != AR.Height() )    
        LogicError
        ("1x2 not conformant:\n",
         DimsString(AL,"Left"),"\n",DimsString(AR,"Right"));
    if( AL.ColAlign() != AR.ColAlign() )
        LogicError("1x2 is misaligned");
}

template<typename T,Dist U,Dist V> 
void
AssertConforming2x1( const ADM<T,U,V>& AT, const ADM<T,U,V>& AB )
{
    if( AT.Width() != AB.Width() )
        LogicError
        ("2x1 is not conformant:\n",
         DimsString(AT,"Top"),"\n",DimsString(AB,"Bottom"));
    if( AT.RowAlign() != AB.RowAlign() )
        LogicError("2x1 is not aligned");
}

template<typename T,Dist U,Dist V> 
void
AssertConforming2x2
( const ADM<T,U,V>& ATL, const ADM<T,U,V>& ATR, 
  const ADM<T,U,V>& ABL, const ADM<T,U,V>& ABR ) 
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

// TODO: Modify to instantiate for all matrix distributions

#define DIAGALIGNED(S,T,U,V)\
  template bool AbstractDistMatrix<T,U,V>::DiagonalAligned\
  ( const DistMatrix<S,UDiag,VDiag>&, Int ) const;\
  template void AbstractDistMatrix<T,U,V>::ForceDiagonalAlign\
  ( DistMatrix<S,UDiag,VDiag>&, Int ) const;

#define DISTPROTO(T,U,V)\
  template class AbstractDistMatrix<T,U,V>;\
  DIAGALIGNED(Int,T,U,V);\
  DIAGALIGNED(float,T,U,V);\
  DIAGALIGNED(double,T,U,V);\
  DIAGALIGNED(Complex<float>,T,U,V);\
  DIAGALIGNED(Complex<double>,T,U,V);
  
#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT

#define DISTPROTO(T,U,V)\
  template class AbstractDistMatrix<T,U,V>;\
  DIAGALIGNED(Int,T,U,V);\
  DIAGALIGNED(float,T,U,V);\
  DIAGALIGNED(double,T,U,V);\
  DIAGALIGNED(Complex<float>,T,U,V);\
  DIAGALIGNED(Complex<double>,T,U,V);
 
#else // ifndef DISABLE_FLOAT

#define DISTPROTO(T,U,V)\
  template class AbstractDistMatrix<T,U,V>;\
  DIAGALIGNED(Int,T,U,V);\
  DIAGALIGNED(double,T,U,V);\
  DIAGALIGNED(Complex<double>,T,U,V);

#endif // ifndef DISABLE_FLOAT
#else // ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT

#define DISTPROTO(T,U,V)\
  template class AbstractDistMatrix<T,U,V>;\
  DIAGALIGNED(Int,T,U,V);\
  DIAGALIGNED(float,T,U,V);\
  DIAGALIGNED(double,T,U,V);

#else // ifndef DISABLE_FLOAT

#define DISTPROTO(T,U,V)\
  template class AbstractDistMatrix<T,U,V>;\
  DIAGALIGNED(Int,T,U,V);\
  DIAGALIGNED(double,T,U,V);

#endif // ifndef DISABLE_FLOAT
#endif // ifndef DISABLE_COMPLEX

#define PROTO(T)\
  DISTPROTO(T,CIRC,CIRC);\
  DISTPROTO(T,MC,  MR  );\
  DISTPROTO(T,MC,  STAR);\
  DISTPROTO(T,MD,  STAR);\
  DISTPROTO(T,MR,  MC  );\
  DISTPROTO(T,MR,  STAR);\
  DISTPROTO(T,STAR,MC  );\
  DISTPROTO(T,STAR,MD  );\
  DISTPROTO(T,STAR,MR  );\
  DISTPROTO(T,STAR,STAR);\
  DISTPROTO(T,STAR,VC  );\
  DISTPROTO(T,STAR,VR  );\
  DISTPROTO(T,VC,  STAR);\
  DISTPROTO(T,VR,  STAR);

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

#ifndef RELEASE

#define DISTCONFORMING(T,U,V) \
  template void AssertConforming1x2\
  ( const ADM<T,U,V>& AL, const ADM<T,U,V>& AR ); \
  template void AssertConforming2x1\
  ( const ADM<T,U,V>& AT, const ADM<T,U,V>& AB ); \
  template void AssertConforming2x2\
  ( const ADM<T,U,V>& ATL, const ADM<T,U,V>& ATR,\
    const ADM<T,U,V>& ABL, const ADM<T,U,V>& ABR )
#define CONFORMING(T)\
  DISTCONFORMING(T,CIRC,CIRC);\
  DISTCONFORMING(T,MC,  MR  );\
  DISTCONFORMING(T,MC,  STAR);\
  DISTCONFORMING(T,MD,  STAR);\
  DISTCONFORMING(T,MR,  MC  );\
  DISTCONFORMING(T,MR,  STAR);\
  DISTCONFORMING(T,STAR,MC  );\
  DISTCONFORMING(T,STAR,MD  );\
  DISTCONFORMING(T,STAR,MR  );\
  DISTCONFORMING(T,STAR,STAR);\
  DISTCONFORMING(T,STAR,VC  );\
  DISTCONFORMING(T,STAR,VR  );\
  DISTCONFORMING(T,VC,  STAR);\
  DISTCONFORMING(T,VR,  STAR);

CONFORMING(Int);
#ifndef DISABLE_FLOAT
CONFORMING(float);
#endif // ifndef DISABLE_FLOAT
CONFORMING(double);
#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
CONFORMING(Complex<float>);
#endif // ifndef DISABLE_FLOAT
CONFORMING(Complex<double>);
#endif // ifndef DISABLE_COMPLEX
#endif // ifndef RELEASE

} // namespace elem
