/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include ELEM_ZEROS_INC

namespace elem {

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
AbstractDistMatrix<T>::AbstractDistMatrix( const elem::Grid& grid, Int root )
: viewType_(OWNER),
  height_(0), width_(0),
  auxMemory_(),
  matrix_(0,0,true),
  colConstrained_(false), rowConstrained_(false), rootConstrained_(false),
  colAlign_(0), rowAlign_(0),
  root_(root), grid_(&grid)
{ }

template<typename T>
AbstractDistMatrix<T>::AbstractDistMatrix( AbstractDistMatrix<T>&& A ) 
ELEM_NOEXCEPT
: viewType_(A.viewType_),
  height_(A.height_), width_(A.width_), 
  colConstrained_(A.colConstrained_), rowConstrained_(A.rowConstrained_),
  rootConstrained_(A.rootConstrained_),
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
AbstractDistMatrix<T>& 
AbstractDistMatrix<T>::operator=( AbstractDistMatrix<T>&& A )
{
    if( Viewing() && !A.Viewing() )
    {
        LogicError("Cannot move a non-view into a viewing AbstractDistMatrix");
    }
    else
    {
        auxMemory_.ShallowSwap( A.auxMemory_ );
        matrix_.ShallowSwap( A.matrix_ );
        viewType_ = A.viewType_;
        height_ = A.height_;
        width_ = A.width_;
        colConstrained_ = A.colConstrained_;
        rowConstrained_ = A.rowConstrained_;
        rootConstrained_ = A.rootConstrained_;
        colAlign_ = A.colAlign_;
        rowAlign_ = A.rowAlign_;
        colShift_ = A.colShift_;
        rowShift_ = A.rowShift_;
        root_ = A.root_;
        grid_ = A.grid_;
    }
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
    rootConstrained_ = false;
    SetShifts();
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
        grid_ = &grid; 
        Empty();
    }
}

template<typename T>
void
AbstractDistMatrix<T>::Resize( Int height, Int width )
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::Resize");
        AssertNotLocked();
        if( Viewing() && (height > height_ || width > width_) )
            LogicError("Tried to increase the size of a view");
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
        if( Viewing() && 
            (height > height_ || width > width_ || ldim > matrix_.LDim()) )
            LogicError("Tried to increase the size of a view");
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
AbstractDistMatrix<T>::MakeConsistent( bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeConsistent"))

    const Int msgLength = 9;
    Int message[msgLength];
    if( CrossRank() == Root() )
    {
        message[0] = viewType_;
        message[1] = height_;
        message[2] = width_;
        message[3] = colConstrained_;
        message[4] = rowConstrained_;
        message[5] = rootConstrained_;
        message[6] = colAlign_;
        message[7] = rowAlign_;
        message[8] = root_;
    }

    const elem::Grid& g = *grid_;
    if( !g.InGrid() && !includingViewers )
        LogicError("Non-participating process called MakeConsistent");
    if( g.InGrid() )
    {
        // TODO: Ensure roots are consistent within each cross communicator
        mpi::Broadcast( message, msgLength, Root(), CrossComm() );
    }
    if( includingViewers )
    {
        const Int vcRoot = g.VCToViewingMap(0);
        mpi::Broadcast( message, msgLength, vcRoot, g.ViewingComm() );
    }
    const ViewType newViewType    = static_cast<ViewType>(message[0]);
    const Int newHeight           = message[1]; 
    const Int newWidth            = message[2];
    const bool newConstrainedCol  = message[3];
    const bool newConstrainedRow  = message[4];
    const bool newConstrainedRoot = message[5];
    const Int newColAlign         = message[6];
    const Int newRowAlign         = message[7];
    const Int root                = message[8];

    root_            = root;
    viewType_        = newViewType;
    colConstrained_  = newConstrainedCol;
    rowConstrained_  = newConstrainedRow;
    rootConstrained_ = newConstrainedRoot;
    colAlign_        = newColAlign;
    rowAlign_        = newRowAlign;

    SetShifts();
    Resize( newHeight, newWidth );
}

template<typename T>
void
AbstractDistMatrix<T>::MakeSizeConsistent( bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeSizeConsistent"))

    const Int msgSize = 2;
    Int message[msgSize];
    if( CrossRank() == Root() )
    {
        message[0] = height_;
        message[1] = width_;
    }

    const elem::Grid& g = *grid_;
    if( !g.InGrid() && !includingViewers )
        LogicError("Non-participating process called MakeSizeConsistent");
    if( g.InGrid() )
        mpi::Broadcast( message, msgSize, Root(), CrossComm() );
    if( includingViewers )
    {
        const Int vcRoot = g.VCToViewingMap(0);
        mpi::Broadcast( message, msgSize, vcRoot, g.ViewingComm() );
    }
    const Int newHeight = message[0]; 
    const Int newWidth  = message[1];
    Resize( newHeight, newWidth );
}

// Realignment
// -----------

template<typename T>
void
AbstractDistMatrix<T>::Align( Int colAlign, Int rowAlign, bool constrain )
{ 
    DEBUG_ONLY(CallStackEntry cse("ADM::Align"))
    const bool requireChange = colAlign_ != colAlign || rowAlign_ != rowAlign;
    DEBUG_ONLY(
        if( Viewing() && requireChange )
            LogicError("Tried to realign a view");
    )
    if( requireChange )
        Empty();
    if( constrain )
    {
        colConstrained_ = true;
        rowConstrained_ = true;
    }
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::AlignCols( Int colAlign, bool constrain )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::AlignCols");
        if( Viewing() && colAlign_ != colAlign )
            LogicError("Tried to realign a view");
    )
    if( colAlign_ != colAlign )
        EmptyData();
    if( constrain )
        colConstrained_ = true;
    colAlign_ = colAlign;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::AlignRows( Int rowAlign, bool constrain )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::AlignRows");
        if( Viewing() && rowAlign_ != rowAlign )
            LogicError("Tried to realign a view");
    )
    if( rowAlign_ != rowAlign )
        EmptyData();
    if( constrain )
        rowConstrained_ = true;
    rowAlign_ = rowAlign;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::FreeAlignments() 
{ 
    if( !Viewing() )
    {
        colConstrained_ = false;
        rowConstrained_ = false;
        rootConstrained_ = false;
    }
    else
        LogicError("Cannot free alignments of views");
}

template<typename T>
void
AbstractDistMatrix<T>::SetRoot( Int root, bool constrain )
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::SetRoot");
        if( root < 0 || root >= mpi::Size(CrossComm()) )
            LogicError("Invalid root");
    )
    if( root != root_ )
        Empty();
    root_ = root;
    if( constrain )
        rootConstrained_ = true;
}

template<typename T>
void
AbstractDistMatrix<T>::AlignWith
( const elem::DistData& data, bool constrain )
{ 
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignWith"))
    AlignColsWith( data, constrain );
    AlignRowsWith( data, constrain );
}

template<typename T>
void
AbstractDistMatrix<T>::AlignColsWith
( const elem::DistData& data, bool constrain )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::AlignColsWith");
        if( colAlign_ != 0 )
            LogicError("Alignment should have been zero");
    )
    SetGrid( *data.grid );
}

template<typename T>
void
AbstractDistMatrix<T>::AlignRowsWith
( const elem::DistData& data, bool constrain )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::AlignRowsWith");
        if( rowAlign_ != 0 )
            LogicError("Alignment should have been zero");
    )
    SetGrid( *data.grid );
}

template<typename T>
void
AbstractDistMatrix<T>::AlignAndResize
( Int colAlign, Int rowAlign, Int height, Int width, 
  bool force, bool constrain )
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
    if( constrain )
    {
        colConstrained_ = true;
        rowConstrained_ = true;
    }
    if( force && (colAlign_ != colAlign || rowAlign_ != rowAlign) )
        LogicError("Could not set alignments"); 
    Resize( height, width );
}

template<typename T>
void
AbstractDistMatrix<T>::AlignColsAndResize
( Int colAlign, Int height, Int width, bool force, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignColsAndResize"))
    if( !Viewing() && (force || !ColConstrained()) )
    {
        colAlign_ = colAlign;
        SetColShift(); 
    }
    if( constrain )
        colConstrained_ = true;
    if( force && colAlign_ != colAlign )
        LogicError("Could not set col alignment");
    Resize( height, width );
}

template<typename T>
void
AbstractDistMatrix<T>::AlignRowsAndResize
( Int rowAlign, Int height, Int width, bool force, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignRowsAndResize"))
    if( !Viewing() && (force || !RowConstrained()) )
    {
        rowAlign_ = rowAlign;
        SetRowShift(); 
    }
    if( constrain )
        rowConstrained_ = true;
    if( force && rowAlign_ != rowAlign )
        LogicError("Could not set row alignment");
    Resize( height, width );
}

// Buffer attachment
// -----------------

template<typename T>
void
AbstractDistMatrix<T>::Attach
( Int height, Int width, const elem::Grid& g, 
  Int colAlign, Int rowAlign, T* buffer, Int ldim, Int root )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Attach"))
    Empty();

    grid_ = &g;
    root_ = root;
    height_ = height;
    width_ = width;
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    colConstrained_ = true;
    rowConstrained_ = true;
    rootConstrained_ = true;
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
( Int height, Int width, const elem::Grid& g,
  Int colAlign, Int rowAlign, elem::Matrix<T>& A, Int root )
{
    // TODO: Assert that the local dimensions are correct
    Attach( height, width, g, colAlign, rowAlign, A.Buffer(), A.LDim(), root );
}

template<typename T>
void
AbstractDistMatrix<T>::LockedAttach
( Int height, Int width, const elem::Grid& g, 
  Int colAlign, Int rowAlign, const T* buffer, Int ldim, Int root )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::LockedAttach"))
    Empty();

    grid_ = &g;
    root_ = root;
    height_ = height;
    width_ = width;
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    colConstrained_ = true;
    rowConstrained_ = true;
    rootConstrained_ = true;
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
( Int height, Int width, const elem::Grid& g, 
  Int colAlign, Int rowAlign, const elem::Matrix<T>& A, Int root )
{
    // TODO: Assert that the local dimensions are correct
    LockedAttach
    ( height, width, g, colAlign, rowAlign, A.LockedBuffer(), A.LDim(), root );
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
bool AbstractDistMatrix<T>::RootConstrained() const { return rootConstrained_; }

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
{ return mpi::Rank(ColComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::RowRank() const
{ return mpi::Rank(RowComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::PartialColRank() const
{ return mpi::Rank(PartialColComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::PartialUnionColRank() const
{ return mpi::Rank(PartialUnionColComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::PartialRowRank() const
{ return mpi::Rank(PartialRowComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::PartialUnionRowRank() const
{ return mpi::Rank(PartialUnionRowComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::DistRank() const
{ return mpi::Rank(DistComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::CrossRank() const
{ return mpi::Rank(CrossComm()); }

template<typename T>
Int
AbstractDistMatrix<T>::RedundantRank() const
{ return mpi::Rank(RedundantComm()); }

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
    return LocalRowOffset(i);
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
    return LocalColOffset(j);
}

template<typename T>
Int
AbstractDistMatrix<T>::LocalRowOffset( Int i ) const
{ return Length_(i,ColShift(),ColStride()); }

template<typename T>
Int
AbstractDistMatrix<T>::LocalColOffset( Int j ) const
{ return Length_(j,RowShift(),RowStride()); }

template<typename T>
Int
AbstractDistMatrix<T>::GlobalRow( Int iLoc ) const
{ return ColShift() + iLoc*ColStride(); }

template<typename T>
Int
AbstractDistMatrix<T>::GlobalCol( Int jLoc ) const
{ return RowShift() + jLoc*RowStride(); }

template<typename T>
bool
AbstractDistMatrix<T>::IsLocalRow( Int i ) const
{ return Participating() && RowOwner(i) == ColRank(); }

template<typename T>
bool
AbstractDistMatrix<T>::IsLocalCol( Int j ) const
{ return Participating() && ColOwner(j) == RowRank(); }

template<typename T>
bool
AbstractDistMatrix<T>::IsLocal( Int i, Int j ) const
{ return IsLocalRow(i) && IsLocalCol(j); }

template<typename T>
mpi::Comm
AbstractDistMatrix<T>::PartialColComm() const
{ return ColComm(); }

template<typename T>
mpi::Comm
AbstractDistMatrix<T>::PartialRowComm() const
{ return RowComm(); }

template<typename T>
mpi::Comm
AbstractDistMatrix<T>::PartialUnionColComm() const
{ return mpi::COMM_SELF; }

template<typename T>
mpi::Comm
AbstractDistMatrix<T>::PartialUnionRowComm() const
{ return mpi::COMM_SELF; }

template<typename T>
Int
AbstractDistMatrix<T>::PartialColStride() const
{ return ColStride(); }

template<typename T>
Int
AbstractDistMatrix<T>::PartialRowStride() const
{ return RowStride(); }

template<typename T>
Int
AbstractDistMatrix<T>::PartialUnionColStride() const
{ return 1; }

template<typename T>
Int
AbstractDistMatrix<T>::PartialUnionRowStride() const
{ return 1; }

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
            const Int iLoc = LocalRow(i);
            const Int jLoc = LocalCol(j);
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
            const Int iLoc = LocalRow(i);
            const Int jLoc = LocalCol(j);
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
                const Int iLoc = LocalRow(i);
                const Int jLoc = LocalCol(j);
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
            const Int iLoc = LocalRow(i);
            const Int jLoc = LocalCol(j);
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
            const Int iLoc = LocalRow(i);
            const Int jLoc = LocalCol(j);
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
            const Int iLoc = LocalRow(i);
            const Int jLoc = LocalCol(j);
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
            const Int iLoc = LocalRow(i);
            const Int jLoc = LocalCol(j);
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
            const Int iLoc = LocalRow(i);
            const Int jLoc = LocalCol(j);
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
            const Int iLoc = LocalRow(i);
            const Int jLoc = LocalCol(j);
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
            const Int iLoc = LocalRow(i);
            const Int jLoc = LocalCol(j);
            MakeLocalReal( iLoc, jLoc );
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
            const Int iLoc = LocalRow(i);
            const Int jLoc = LocalCol(j);
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
AbstractDistMatrix<T>::MakeLocalReal( Int iLoc, Int jLoc )
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
    const Int height = Height();
    const Int localWidth = LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = GlobalCol(jLoc);
        if( j < height && IsLocal(j,j) )
        {
            const Int iLoc = LocalRow(j);
            MakeLocalReal( iLoc, jLoc );
        }
    }
}

template<typename T>
void
AbstractDistMatrix<T>::ConjugateDiagonal( Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::ConjugateDiagonal"))
    const Int height = Height();
    const Int localWidth = LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = GlobalCol(jLoc);
        if( j < height && IsLocal(j,j) )
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
AbstractDistMatrix<T>::GetSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  DistMatrix<T,STAR,STAR>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetSubmatrix"))
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
AbstractDistMatrix<T>::GetRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  DistMatrix<Base<T>,STAR,STAR>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetRealPartOfSubmatrix"))
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
AbstractDistMatrix<T>::GetImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  DistMatrix<Base<T>,STAR,STAR>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetImagPartOfSubmatrix"))
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
AbstractDistMatrix<T>::GetSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetSubmatrix"))
    DistMatrix<T,STAR,STAR> ASub( Grid() );
    GetSubmatrix( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
DistMatrix<Base<T>,STAR,STAR>
AbstractDistMatrix<T>::GetRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetRealPartOfSubmatrix"))
    DistMatrix<Base<T>,STAR,STAR> ASub( Grid() );
    GetRealPartOfSubmatrix( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
DistMatrix<Base<T>,STAR,STAR>
AbstractDistMatrix<T>::GetImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetImagPartOfSubmatrix"))
    DistMatrix<Base<T>,STAR,STAR> ASub( Grid() );
    GetImagPartOfSubmatrix( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
void 
AbstractDistMatrix<T>::SetSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const DistMatrix<T,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetSubmatrix"))
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
AbstractDistMatrix<T>::SetRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const DistMatrix<Base<T>,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetRealPartOfSubmatrix"))
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
AbstractDistMatrix<T>::SetImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const DistMatrix<Base<T>,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetImagPartOfSubmatrix"))
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
AbstractDistMatrix<T>::UpdateSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  T alpha, const DistMatrix<T,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateSubmatrix"))
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
AbstractDistMatrix<T>::UpdateRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateRealPartOfSubmatrix"))
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
AbstractDistMatrix<T>::UpdateImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateImagPartOfSubmatrix"))
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
AbstractDistMatrix<T>::MakeSubmatrixReal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeSubmatrixReal"))
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
                        MakeLocalReal( iLoc, jLoc );
                    }
                }
            }
        }
    }
}

template<typename T>
void 
AbstractDistMatrix<T>::ConjugateSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::ConjugateSubmatrix"))
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
AbstractDistMatrix<T>::GetLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
  elem::Matrix<T>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetLocalSubmatrix"))
    LockedMatrix().GetSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void
AbstractDistMatrix<T>::GetRealPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
  elem::Matrix<Base<T>>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetRealPartOfLocalSubmatrix"))
    LockedMatrix().GetRealPartOfSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void
AbstractDistMatrix<T>::GetImagPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
  elem::Matrix<Base<T>>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::GetImagPartOfLocalSubmatrix"))
    LockedMatrix().GetImagPartOfSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::SetLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  const elem::Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetLocalSubmatrix"))
    Matrix().SetSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::SetRealPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  const elem::Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetRealPartOfLocalSubmatrix"))
    Matrix().SetRealPartOfSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::SetImagPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  const elem::Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetImagPartOfLocalSubmatrix"))
    Matrix().SetImagPartOfSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::UpdateLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  T alpha, const elem::Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateLocalSubmatrix"))
    Matrix().UpdateSubmatrix( rowIndLoc, colIndLoc, alpha, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::UpdateRealPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  Base<T> alpha, const elem::Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateRealPartOfLocalSubmatrix"))
    Matrix().UpdateRealPartOfSubmatrix( rowIndLoc, colIndLoc, alpha, ASub );
}

template<typename T>
void 
AbstractDistMatrix<T>::UpdateImagPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  Base<T> alpha, const elem::Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateImagPartOfLocalSubmatrix"))
    Matrix().UpdateImagPartOfSubmatrix( rowIndLoc, colIndLoc, alpha, ASub );
}

template<typename T>
void
AbstractDistMatrix<T>::MakeLocalSubmatrixReal
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeLocalSubmatrixReal"))
    Matrix().MakeSubmatrixReal( rowIndLoc, colIndLoc );
}

template<typename T>
void
AbstractDistMatrix<T>::ConjugateLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::ConjugateLocalSubmatrix"))
    Matrix().ConjugateSubmatrix( rowIndLoc, colIndLoc );
}

// Sum the local matrix over a particular communicator
// ===================================================
// NOTE: The matrix dimensions *must* be uniform over the communicator.

template<typename T>
void
AbstractDistMatrix<T>::SumOver( mpi::Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SumOver"))
    if( !Participating() )
        return;

    const Int localHeight = LocalHeight();
    const Int localWidth = LocalWidth();
    const Int localSize = mpi::Pad( localHeight*localWidth );
    T* sumBuf = auxMemory_.Require( localSize );   

    // Pack
    T* buf = Buffer();
    const Int ldim = LDim(); 
    ELEM_PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* thisCol = &buf[jLoc*ldim];
        T* sumCol = &sumBuf[jLoc*localHeight];
        MemCopy( sumCol, thisCol, localHeight );
    }

    // AllReduce sum
    mpi::AllReduce( sumBuf, localSize, comm );

    // Unpack
    ELEM_PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* sumCol = &sumBuf[jLoc*localHeight];
        T* thisCol = &buf[jLoc*ldim];
        MemCopy( thisCol, sumCol, localHeight );
    } 
    auxMemory_.Release();
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

// Exchange metadata with another matrix
// =====================================

template<typename T>
void 
AbstractDistMatrix<T>::ShallowSwap( AbstractDistMatrix<T>& A )
{
    matrix_.ShallowSwap( A.matrix_ );
    auxMemory_.ShallowSwap( A.auxMemory_ );
    std::swap( viewType_, A.viewType_ );
    std::swap( height_ , A.height_ );
    std::swap( width_, A.width_ );
    std::swap( colConstrained_, A.colConstrained_ );
    std::swap( rowConstrained_, A.rowConstrained_ );
    std::swap( rootConstrained_, A.rootConstrained_ );
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
AssertConforming1x2
( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR )
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
AssertConforming2x1
( const AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& AB )
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
( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR, 
  const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR ) 
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

#ifndef ELEM_RELEASE
 #define PROTO(T) \
  template class AbstractDistMatrix<T>;\
  template void AssertConforming1x2\
  ( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR );\
  template void AssertConforming2x1\
  ( const AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& AB );\
  template void AssertConforming2x2\
  ( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR,\
    const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR )
#else
 #define PROTO(T) template class AbstractDistMatrix<T>
#endif
 
#ifndef ELEM_DISABLE_COMPLEX
 #ifndef ELEM_DISABLE_FLOAT
  PROTO(Int);
  PROTO(float);
  PROTO(double);
  PROTO(Complex<float>);
  PROTO(Complex<double>);
 #else // ifndef ELEM_DISABLE_FLOAT
  PROTO(Int);
  PROTO(double);
  PROTO(Complex<double>);
 #endif // ifndef ELEM_DISABLE_FLOAT
#else // ifndef ELEM_DISABLE_COMPLEX
 #ifndef ELEM_DISABLE_FLOAT
  PROTO(Int);
  PROTO(float);
  PROTO(double);
 #else // ifndef ELEM_DISABLE_FLOAT
  PROTO(Int);
  PROTO(double);
 #endif // ifndef ELEM_DISABLE_FLOAT
#endif // ifndef ELEM_DISABLE_COMPLEX

} // namespace elem
