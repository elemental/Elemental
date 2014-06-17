/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"


namespace El {

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
AbstractBlockDistMatrix<T>::AbstractBlockDistMatrix
( const El::Grid& g, Int root )
: viewType_(OWNER),
  height_(0), width_(0),
  auxMemory_(),
  matrix_(0,0,true),
  colConstrained_(false), rowConstrained_(false), rootConstrained_(false),
  blockHeight_(DefaultBlockHeight()), blockWidth_(DefaultBlockWidth()),
  colAlign_(0), rowAlign_(0),
  colCut_(0), rowCut_(0),
  root_(root), grid_(&g)
{ }

template<typename T>
AbstractBlockDistMatrix<T>::AbstractBlockDistMatrix
( const El::Grid& g, Int blockHeight, Int blockWidth, Int root )
: viewType_(OWNER),
  height_(0), width_(0),
  auxMemory_(),
  matrix_(0,0,true),
  colConstrained_(true), rowConstrained_(true), rootConstrained_(false),
  blockHeight_(blockHeight), blockWidth_(blockWidth),
  colAlign_(0), rowAlign_(0),
  colCut_(0), rowCut_(0),
  root_(root), grid_(&g)
{ }

template<typename T>
AbstractBlockDistMatrix<T>::AbstractBlockDistMatrix
( AbstractBlockDistMatrix<T>&& A ) EL_NOEXCEPT
: viewType_(A.viewType_),
  height_(A.height_), width_(A.width_), 
  colConstrained_(A.colConstrained_), rowConstrained_(A.rowConstrained_),
  rootConstrained_(A.rootConstrained_),
  blockHeight_(A.blockHeight_), blockWidth_(A.blockWidth_),
  colAlign_(A.colAlign_), rowAlign_(A.rowAlign_),
  colCut_(A.colCut_), rowCut_(A.rowCut_),
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
AbstractBlockDistMatrix<T>::~AbstractBlockDistMatrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
AbstractBlockDistMatrix<T>& 
AbstractBlockDistMatrix<T>::operator=( AbstractBlockDistMatrix<T>&& A )
{
    if( Viewing() && !A.Viewing() )
    {
        LogicError
        ("Cannot move a non-view into a viewing AbstractBlockDistMatrix");
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
        blockHeight_ = A.blockHeight_;
        blockWidth_ = A.blockWidth_;
        colAlign_ = A.colAlign_;
        rowAlign_ = A.rowAlign_;
        colCut_ = A.colCut_;
        rowCut_ = A.rowCut_;
        colShift_ = A.colShift_;
        rowShift_ = A.rowShift_;
        root_ = A.root_;
        grid_ = A.grid_;
    }
    return *this;
}

template<typename T>
void
AbstractBlockDistMatrix<T>::Empty()
{
    matrix_.Empty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
    blockHeight_ = 0;
    blockWidth_ = 0;
    colAlign_ = 0;
    rowAlign_ = 0;
    colCut_ = 0;
    rowCut_ = 0;
    colConstrained_ = false;
    rowConstrained_ = false;
    rootConstrained_ = false;
}

template<typename T>
void
AbstractBlockDistMatrix<T>::EmptyData()
{
    matrix_.Empty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
}

template<typename T>
void
AbstractBlockDistMatrix<T>::SetGrid( const El::Grid& grid )
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
AbstractBlockDistMatrix<T>::Resize( Int height, Int width )
{
    DEBUG_ONLY(
        CallStackEntry cse("ABDM::Resize");
        AssertNotLocked();
    )
    height_ = height; 
    width_ = width;
    if( Participating() )
        matrix_.Resize_
        ( BlockedLength(height,ColShift(),BlockHeight(),ColCut(),ColStride()),
          BlockedLength(width,RowShift(),BlockWidth(),RowCut(),RowStride()) );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::Resize( Int height, Int width, Int ldim )
{
    DEBUG_ONLY(
        CallStackEntry cse("ABDM::Resize");
        AssertNotLocked();
    )
    height_ = height; 
    width_ = width;
    if( Participating() )
        matrix_.Resize_
        ( BlockedLength(height,ColShift(),BlockHeight(),ColCut(),ColStride()),
          BlockedLength(width,RowShift(),BlockWidth(),RowCut(),RowStride()), 
          ldim );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::MakeConsistent( bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::MakeConsistent"))

    const Int msgLength = 13;
    Int message[msgLength];
    if( CrossRank() == Root() )
    {
        message[ 0] = viewType_;
        message[ 1] = height_;
        message[ 2] = width_;
        message[ 3] = colConstrained_;
        message[ 4] = rowConstrained_;
        message[ 5] = rootConstrained_;
        message[ 6] = blockHeight_;
        message[ 7] = blockWidth_;
        message[ 8] = colAlign_;
        message[ 9] = rowAlign_;
        message[10] = colCut_;
        message[11] = rowCut_;
        message[12] = root_;
    }

    const El::Grid& g = *grid_;
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
    const Int newHeight           = message[ 1]; 
    const Int newWidth            = message[ 2];
    const bool newConstrainedCol  = message[ 3];
    const bool newConstrainedRow  = message[ 4];
    const bool newConstrainedRoot = message[ 5];
    const Int newBlockHeight      = message[ 6];
    const Int newBlockWidth       = message[ 7];
    const Int newColAlign         = message[ 8];
    const Int newRowAlign         = message[ 9];
    const Int newColCut           = message[10];
    const Int newRowCut           = message[11];
    const Int root                = message[12];

    root_            = root;
    viewType_        = newViewType;
    colConstrained_  = newConstrainedCol;
    rowConstrained_  = newConstrainedRow;
    rootConstrained_ = newConstrainedRoot;
    blockHeight_     = newBlockHeight;
    blockWidth_      = newBlockWidth;
    colAlign_        = newColAlign;
    rowAlign_        = newRowAlign;
    colCut_          = newColCut;
    rowCut_          = newRowCut;

    SetShifts();
    Resize( newHeight, newWidth );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::MakeSizeConsistent( bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::MakeSizeConsistent"))

    const Int msgLength = 2;
    Int message[msgLength];
    if( CrossRank() == Root() )
    {
        message[0] = height_;
        message[1] = width_;
    }

    const El::Grid& g = *grid_;
    if( !g.InGrid() && !includingViewers )
        LogicError("Non-participating process called MakeSizeConsistent");
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
    const Int newHeight = message[0]; 
    const Int newWidth  = message[1];
    Resize( newHeight, newWidth );
}

// Realignment
// -----------

template<typename T>
void
AbstractBlockDistMatrix<T>::Align
( Int blockHeight, Int blockWidth, 
  Int colAlign, Int rowAlign, Int colCut, Int rowCut, bool constrain )
{ 
    DEBUG_ONLY(CallStackEntry cse("ABDM::Align"))
    const bool requireChange = 
        blockHeight_ != blockHeight || blockWidth_ != blockWidth ||
        colAlign_    != colAlign    || rowAlign_   != rowAlign   ||
        colCut_      != colCut      || rowCut_     != rowCut;
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
    blockHeight_ = blockHeight;
    blockWidth_ = blockWidth;
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    colCut_ = colCut;
    rowCut_ = rowCut;
    SetShifts();
}

template<typename T>
void
AbstractBlockDistMatrix<T>::AlignCols
( Int blockHeight, Int colAlign, Int colCut, bool constrain )
{ 
    DEBUG_ONLY(CallStackEntry cse("ABDM::AlignCols"))
    const bool requireChange = 
        blockHeight_ != blockHeight || 
        colAlign_    != colAlign    || 
        colCut_      != colCut;
    DEBUG_ONLY(
        if( Viewing() && requireChange )
            LogicError("Tried to realign a view");
    )
    if( requireChange )
        EmptyData();
    if( constrain )
        colConstrained_ = true;
    blockHeight_ = blockHeight;
    colAlign_ = colAlign;
    colCut_ = colCut;
    SetShifts();
}

template<typename T>
void
AbstractBlockDistMatrix<T>::AlignRows
( Int blockWidth, Int rowAlign, Int rowCut, bool constrain )
{ 
    DEBUG_ONLY(CallStackEntry cse("ABDM::AlignRows"))
    const bool requireChange = 
        blockWidth_ != blockWidth || 
        rowAlign_   != rowAlign   || 
        rowCut_     != rowCut;
    DEBUG_ONLY(
        if( Viewing() && requireChange )
            LogicError("Tried to realign a view");
    )
    if( requireChange )
        EmptyData();
    if( constrain )
        rowConstrained_ = true;
    blockWidth_ = blockWidth;
    rowAlign_ = rowAlign;
    rowCut_ = rowCut;
    SetShifts();
}

template<typename T>
void
AbstractBlockDistMatrix<T>::FreeAlignments() 
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
AbstractBlockDistMatrix<T>::SetRoot( Int root, bool constrain )
{
    DEBUG_ONLY(
        CallStackEntry cse("ABDM::SetRoot");
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
AbstractBlockDistMatrix<T>::AlignWith
( const El::BlockDistData& data, bool constrain )
{ 
    DEBUG_ONLY(CallStackEntry cse("ABDM::AlignWith"))
    AlignColsWith( data, constrain );
    AlignRowsWith( data, constrain );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::AlignColsWith
( const El::BlockDistData& data, bool constrain )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ABDM::AlignColsWith");
        if( colAlign_ != 0 )
            LogicError("Alignment should have been zero");
    )
    SetGrid( *data.grid );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::AlignRowsWith
( const El::BlockDistData& data, bool constrain )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ABDM::AlignRowsWith");
        if( rowAlign_ != 0 )
            LogicError("Alignment should have been zero");
    )
    SetGrid( *data.grid );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::AlignAndResize
( Int blockHeight, Int blockWidth, 
  Int colAlign, Int rowAlign, Int colCut, Int rowCut,
  Int height, Int width, bool force, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::AlignAndResize"))
    if( !Viewing() )
    {
        if( force || !ColConstrained() )
        {
            blockHeight_ = blockHeight;
            colAlign_ = colAlign;
            colCut_ = colCut;
            SetColShift(); 
        }
        if( force || !RowConstrained() )
        {
            blockWidth_ = blockWidth;
            rowAlign_ = rowAlign;
            rowCut_ = rowCut;
            SetRowShift();
        }
    }
    if( constrain )
    {
        colConstrained_ = true;
        rowConstrained_ = true;
    }
    if( force && 
        (blockHeight_ != blockHeight || blockWidth_ != blockWidth || 
         colAlign_    != colAlign    || rowAlign_   != rowAlign   ||
         colCut_      != colCut      || rowCut_     != rowCut) )
        LogicError("Could not set alignments and cuts"); 
    Resize( height, width );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::AlignColsAndResize
( Int blockHeight, Int colAlign, Int colCut, Int height, Int width, 
  bool force, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::AlignColsAndResize"))
    if( !Viewing() && (force || !ColConstrained()) )
    {
        blockHeight_ = blockHeight;
        colAlign_ = colAlign;
        colCut_ = colCut;
        SetColShift(); 
    }
    if( constrain )
        colConstrained_ = true;
    if( force && 
        (colAlign_ != colAlign || colCut_ != colCut || 
         blockHeight_ != blockHeight) )
        LogicError("Could not set col alignment and cut");
    Resize( height, width );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::AlignRowsAndResize
( Int blockWidth, Int rowAlign, Int rowCut, Int height, Int width, 
  bool force, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::AlignRowsAndResize"))
    if( !Viewing() && (force || !RowConstrained()) )
    {
        blockWidth_ = blockWidth;
        rowAlign_ = rowAlign;
        rowCut_ = rowCut;
        SetRowShift(); 
    }
    if( constrain )
        rowConstrained_ = true;
    if( force && 
        (rowAlign_ != rowAlign || rowCut_ != rowCut ||
         blockWidth_ != blockWidth) )
        LogicError("Could not set row alignment and cut");
    Resize( height, width );
}

// Buffer attachment
// -----------------

template<typename T>
void
AbstractBlockDistMatrix<T>::Attach
( Int height, Int width, const El::Grid& g, 
  Int blockHeight, Int blockWidth, 
  Int colAlign, Int rowAlign, Int colCut, Int rowCut,
  T* buffer, Int ldim, Int root )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::Attach"))
    Empty();

    grid_ = &g;
    root_ = root;
    height_ = height;
    width_ = width;
    blockHeight_ = blockHeight;
    blockWidth_ = blockWidth;
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    colCut_ = colCut;
    rowCut_ = rowCut;
    colConstrained_ = true;
    rowConstrained_ = true;
    viewType_ = VIEW;
    SetShifts();
    if( Participating() )
    {
        Int localHeight = 
            BlockedLength(height,colShift_,blockHeight,colCut,ColStride());
        Int localWidth = 
            BlockedLength(width,rowShift_,blockWidth,rowCut,RowStride());
        matrix_.Attach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T>
void
AbstractBlockDistMatrix<T>::Attach
( Int height, Int width, const El::Grid& g,
  Int blockHeight, Int blockWidth, 
  Int colAlign, Int rowAlign, Int colCut, Int rowCut, El::Matrix<T>& A, 
  Int root )
{
    // TODO: Assert that the local dimensions are correct
    Attach
    ( height, width, g, blockHeight, blockWidth, 
      colAlign, rowAlign, colCut, rowCut, A.Buffer(), A.LDim(), root );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::LockedAttach
( Int height, Int width, const El::Grid& g, 
  Int blockHeight, Int blockWidth, 
  Int colAlign, Int rowAlign, Int colCut, Int rowCut,
  const T* buffer, Int ldim, Int root )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::LockedAttach"))
    Empty();

    grid_ = &g;
    root_ = root;
    height_ = height;
    width_ = width;
    blockHeight_ = blockHeight;
    blockWidth_ = blockWidth;
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    colCut_ = colCut;
    rowCut_ = rowCut;
    colConstrained_ = true;
    rowConstrained_ = true;
    viewType_ = LOCKED_VIEW;
    SetShifts();
    if( Participating() )
    {
        Int localHeight = 
            BlockedLength(height,colShift_,blockHeight,colCut,ColStride());
        Int localWidth = 
            BlockedLength(width,rowShift_,blockWidth,rowCut,RowStride());
        matrix_.LockedAttach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T>
void
AbstractBlockDistMatrix<T>::LockedAttach
( Int height, Int width, const El::Grid& g,
  Int blockHeight, Int blockWidth, 
  Int colAlign, Int rowAlign, Int colCut, Int rowCut, const El::Matrix<T>& A,
  Int root )
{
    // TODO: Assert that the local dimensions are correct
    LockedAttach
    ( height, width, g, blockHeight, blockWidth, 
      colAlign, rowAlign, colCut, rowCut, A.LockedBuffer(), A.LDim(), root );
}

// Basic queries
// =============

// Global matrix information
// -------------------------

template<typename T>
Int AbstractBlockDistMatrix<T>::Height() const { return height_; }
template<typename T>
Int AbstractBlockDistMatrix<T>::Width() const { return width_; }

template<typename T>
Int AbstractBlockDistMatrix<T>::DiagonalLength( Int offset ) const
{ return El::DiagonalLength(height_,width_,offset); }

template<typename T>
bool AbstractBlockDistMatrix<T>::Viewing() const 
{ return IsViewing( viewType_ ); }
template<typename T>
bool AbstractBlockDistMatrix<T>::Locked() const 
{ return IsLocked( viewType_ ); }

// Local matrix information
// ------------------------

template<typename T>
Int AbstractBlockDistMatrix<T>::LocalHeight() const { return matrix_.Height(); }
template<typename T>
Int AbstractBlockDistMatrix<T>::LocalWidth() const { return matrix_.Width(); }
template<typename T>
Int AbstractBlockDistMatrix<T>::LDim() const { return matrix_.LDim(); }

template<typename T>
El::Matrix<T>& 
AbstractBlockDistMatrix<T>::Matrix() { return matrix_; }
template<typename T>
const El::Matrix<T>& 
AbstractBlockDistMatrix<T>::LockedMatrix() const { return matrix_; }

template<typename T>
size_t
AbstractBlockDistMatrix<T>::AllocatedMemory() const 
{ return matrix_.MemorySize(); }

template<typename T>
T*
AbstractBlockDistMatrix<T>::Buffer() { return matrix_.Buffer(); }

template<typename T>
T*
AbstractBlockDistMatrix<T>::Buffer( Int iLoc, Int jLoc )
{ return matrix_.Buffer(iLoc,jLoc); }

template<typename T>
const T*
AbstractBlockDistMatrix<T>::LockedBuffer() const
{ return matrix_.LockedBuffer(); }

template<typename T>
const T*
AbstractBlockDistMatrix<T>::LockedBuffer( Int iLoc, Int jLoc ) const
{ return matrix_.LockedBuffer(iLoc,jLoc); }

// Distribution information
// ------------------------

template<typename T>
const El::Grid& AbstractBlockDistMatrix<T>::Grid() const { return *grid_; }

template<typename T>
Int AbstractBlockDistMatrix<T>::BlockHeight() const { return blockHeight_; }
template<typename T>
Int AbstractBlockDistMatrix<T>::BlockWidth() const { return blockWidth_; }

template<typename T>
Int AbstractBlockDistMatrix<T>::ColAlign() const { return colAlign_; }
template<typename T>
Int AbstractBlockDistMatrix<T>::RowAlign() const { return rowAlign_; }

template<typename T>
Int AbstractBlockDistMatrix<T>::ColCut() const { return colCut_; }
template<typename T>
Int AbstractBlockDistMatrix<T>::RowCut() const { return rowCut_; }

template<typename T>
Int AbstractBlockDistMatrix<T>::ColShift() const { return colShift_; }
template<typename T>
Int AbstractBlockDistMatrix<T>::RowShift() const { return rowShift_; }

template<typename T>
bool AbstractBlockDistMatrix<T>::ColConstrained() const 
{ return colConstrained_; }
template<typename T>
bool AbstractBlockDistMatrix<T>::RowConstrained() const 
{ return rowConstrained_; }
template<typename T>
bool AbstractBlockDistMatrix<T>::RootConstrained() const
{ return rootConstrained_; }

template<typename T>
bool AbstractBlockDistMatrix<T>::Participating() const
{ return grid_->InGrid() && (CrossRank()==root_); }

template<typename T>
Int AbstractBlockDistMatrix<T>::RowOwner( Int i ) const
{ return (((i+ColCut())/BlockHeight())+ColAlign()) % ColStride(); }
template<typename T>
Int AbstractBlockDistMatrix<T>::ColOwner( Int j ) const
{ return (((j+RowCut())/BlockWidth())+RowAlign()) % RowStride(); }
template<typename T>
Int AbstractBlockDistMatrix<T>::Owner( Int i, Int j ) const
{ return RowOwner(i)+ColOwner(j)*ColStride(); }

template<typename T>
Int AbstractBlockDistMatrix<T>::LocalRow( Int i ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ABDM::LocalRow");
        if( !IsLocalRow(i) )
            LogicError("Requested local index of non-local row");
    )
    return LocalRowOffset(i);
}

template<typename T>
Int AbstractBlockDistMatrix<T>::LocalCol( Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ABDM::LocalCol");
        if( !IsLocalCol(j) )
            LogicError("Requested local index of non-local column");
    )
    return LocalColOffset(j);
}

template<typename T>
Int AbstractBlockDistMatrix<T>::LocalRowOffset( Int i ) const
{ return BlockedLength_
         ( i, ColShift(), BlockHeight(), ColCut(), ColStride() ); }
template<typename T>
Int AbstractBlockDistMatrix<T>::LocalColOffset( Int j ) const
{ return BlockedLength_( j, RowShift(), BlockWidth(), RowCut(), RowStride() ); }

template<typename T>
Int AbstractBlockDistMatrix<T>::GlobalRow( Int iLoc ) const
{ 
    return GlobalBlockedIndex
           (iLoc,ColShift(),BlockHeight(),ColCut(),ColStride()); 
}
template<typename T>
Int AbstractBlockDistMatrix<T>::GlobalCol( Int jLoc ) const
{ 
    return GlobalBlockedIndex
           (jLoc,RowShift(),BlockWidth(),RowCut(),RowStride()); 
}

template<typename T>
bool AbstractBlockDistMatrix<T>::IsLocalRow( Int i ) const
{ return Participating() && RowOwner(i) == ColRank(); }
template<typename T>
bool AbstractBlockDistMatrix<T>::IsLocalCol( Int j ) const
{ return Participating() && ColOwner(j) == RowRank(); }
template<typename T>
bool AbstractBlockDistMatrix<T>::IsLocal( Int i, Int j ) const
{ return IsLocalRow(i) && IsLocalCol(j); }

template<typename T>
mpi::Comm AbstractBlockDistMatrix<T>::PartialColComm() const
{ return ColComm(); }
template<typename T>
mpi::Comm AbstractBlockDistMatrix<T>::PartialRowComm() const
{ return RowComm(); }

template<typename T>
mpi::Comm AbstractBlockDistMatrix<T>::PartialUnionColComm() const
{ return mpi::COMM_SELF; }
template<typename T>
mpi::Comm AbstractBlockDistMatrix<T>::PartialUnionRowComm() const
{ return mpi::COMM_SELF; }

template<typename T>
Int AbstractBlockDistMatrix<T>::PartialColStride() const
{ return ColStride(); }
template<typename T>
Int AbstractBlockDistMatrix<T>::PartialRowStride() const
{ return RowStride(); }

template<typename T>
Int AbstractBlockDistMatrix<T>::PartialUnionColStride() const
{ return 1; }
template<typename T>
Int AbstractBlockDistMatrix<T>::PartialUnionRowStride() const
{ return 1; }

template<typename T>
Int AbstractBlockDistMatrix<T>::ColRank() const { return mpi::Rank(ColComm()); }
template<typename T>
Int AbstractBlockDistMatrix<T>::RowRank() const { return mpi::Rank(RowComm()); }

template<typename T>
Int AbstractBlockDistMatrix<T>::PartialColRank() const
{ return mpi::Rank(PartialColComm()); }
template<typename T>
Int AbstractBlockDistMatrix<T>::PartialRowRank() const
{ return mpi::Rank(PartialRowComm()); }

template<typename T>
Int AbstractBlockDistMatrix<T>::PartialUnionColRank() const
{ return mpi::Rank(PartialUnionColComm()); }
template<typename T>
Int AbstractBlockDistMatrix<T>::PartialUnionRowRank() const
{ return mpi::Rank(PartialUnionRowComm()); }

template<typename T>
Int AbstractBlockDistMatrix<T>::DistRank() const
{ return mpi::Rank(DistComm()); }
template<typename T>
Int AbstractBlockDistMatrix<T>::CrossRank() const
{ return mpi::Rank(CrossComm()); }
template<typename T>
Int AbstractBlockDistMatrix<T>::RedundantRank() const
{ return mpi::Rank(RedundantComm()); }

template<typename T>
Int AbstractBlockDistMatrix<T>::Root() const { return root_; }

// Single-entry manipulation
// =========================

// Global entry manipulation
// -------------------------

template<typename T>
T
AbstractBlockDistMatrix<T>::Get( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ABDM::Get");
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
AbstractBlockDistMatrix<T>::GetRealPart( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ABDM::GetRealPart");
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
AbstractBlockDistMatrix<T>::GetImagPart( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ABDM::GetImagPart");
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
AbstractBlockDistMatrix<T>::Set( Int i, Int j, T value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::Set"))
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
AbstractBlockDistMatrix<T>::SetRealPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::SetRealPart"))
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
AbstractBlockDistMatrix<T>::SetImagPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::SetImagPart"))
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
AbstractBlockDistMatrix<T>::Update( Int i, Int j, T value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::Update"))
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
AbstractBlockDistMatrix<T>::UpdateRealPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::UpdateRealPart"))
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
AbstractBlockDistMatrix<T>::UpdateImagPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::UpdateImagPart"))
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
AbstractBlockDistMatrix<T>::MakeReal( Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::MakeReal"))
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
AbstractBlockDistMatrix<T>::Conjugate( Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::Conjugate"))
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
AbstractBlockDistMatrix<T>::GetLocal( Int i, Int j ) const
{ return matrix_.Get(i,j); }

template<typename T>
Base<T>
AbstractBlockDistMatrix<T>::GetLocalRealPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetRealPart(iLoc,jLoc); }

template<typename T>
Base<T>
AbstractBlockDistMatrix<T>::GetLocalImagPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetImagPart(iLoc,jLoc); }

template<typename T>
void
AbstractBlockDistMatrix<T>::SetLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Set(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractBlockDistMatrix<T>::SetLocalRealPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.SetRealPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractBlockDistMatrix<T>::SetLocalImagPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.SetImagPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractBlockDistMatrix<T>::UpdateLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Update(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractBlockDistMatrix<T>::UpdateLocalRealPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.UpdateRealPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractBlockDistMatrix<T>::UpdateLocalImagPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.UpdateImagPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractBlockDistMatrix<T>::MakeLocalReal( Int iLoc, Int jLoc )
{ matrix_.MakeReal( iLoc, jLoc ); }

template<typename T>
void
AbstractBlockDistMatrix<T>::ConjugateLocal( Int iLoc, Int jLoc )
{ matrix_.Conjugate( iLoc, jLoc ); }

// Diagonal manipulation
// =====================

template<typename T>
void
AbstractBlockDistMatrix<T>::MakeDiagonalReal( Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::MakeDiagonalReal"))
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
AbstractBlockDistMatrix<T>::ConjugateDiagonal( Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::ConjugateDiagonal"))
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
AbstractBlockDistMatrix<T>::GetSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  DistMatrix<T,STAR,STAR>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::GetSubmatrix"))
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
AbstractBlockDistMatrix<T>::GetRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  DistMatrix<Base<T>,STAR,STAR>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::GetRealPartOfSubmatrix"))
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
AbstractBlockDistMatrix<T>::GetImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd, 
  DistMatrix<Base<T>,STAR,STAR>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::GetImagPartOfSubmatrix"))
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
AbstractBlockDistMatrix<T>::GetSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::GetSubmatrix"))
    DistMatrix<T,STAR,STAR> ASub( Grid() );
    GetSubmatrix( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
DistMatrix<Base<T>,STAR,STAR>
AbstractBlockDistMatrix<T>::GetRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::GetRealPartOfSubmatrix"))
    DistMatrix<Base<T>,STAR,STAR> ASub( Grid() );
    GetRealPartOfSubmatrix( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
DistMatrix<Base<T>,STAR,STAR>
AbstractBlockDistMatrix<T>::GetImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::GetImagPartOfSubmatrix"))
    DistMatrix<Base<T>,STAR,STAR> ASub( Grid() );
    GetImagPartOfSubmatrix( rowInd, colInd, ASub );
    return ASub;
}

template<typename T>
void 
AbstractBlockDistMatrix<T>::SetSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const DistMatrix<T,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::SetSubmatrix"))
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
AbstractBlockDistMatrix<T>::SetRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const DistMatrix<Base<T>,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::SetRealPartOfSubmatrix"))
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
AbstractBlockDistMatrix<T>::SetImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  const DistMatrix<Base<T>,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::SetImagPartOfSubmatrix"))
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
AbstractBlockDistMatrix<T>::UpdateSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  T alpha, const DistMatrix<T,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::UpdateSubmatrix"))
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
AbstractBlockDistMatrix<T>::UpdateRealPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::UpdateRealPartOfSubmatrix"))
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
AbstractBlockDistMatrix<T>::UpdateImagPartOfSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
  Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::UpdateImagPartOfSubmatrix"))
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
AbstractBlockDistMatrix<T>::MakeSubmatrixReal
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::MakeSubmatrixReal"))
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
AbstractBlockDistMatrix<T>::ConjugateSubmatrix
( const std::vector<Int>& rowInd, const std::vector<Int>& colInd )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::ConjugateSubmatrix"))
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
AbstractBlockDistMatrix<T>::GetLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
  El::Matrix<T>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::GetLocalSubmatrix"))
    LockedMatrix().GetSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::GetRealPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
  El::Matrix<Base<T>>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::GetRealPartOfLocalSubmatrix"))
    LockedMatrix().GetRealPartOfSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::GetImagPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc, 
  El::Matrix<Base<T>>& ASub ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::GetImagPartOfLocalSubmatrix"))
    LockedMatrix().GetImagPartOfSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void 
AbstractBlockDistMatrix<T>::SetLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  const El::Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::SetLocalSubmatrix"))
    Matrix().SetSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void 
AbstractBlockDistMatrix<T>::SetRealPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  const El::Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::SetRealPartOfLocalSubmatrix"))
    Matrix().SetRealPartOfSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void 
AbstractBlockDistMatrix<T>::SetImagPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  const El::Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::SetImagPartOfLocalSubmatrix"))
    Matrix().SetImagPartOfSubmatrix( rowIndLoc, colIndLoc, ASub );
}

template<typename T>
void 
AbstractBlockDistMatrix<T>::UpdateLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  T alpha, const El::Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::UpdateLocalSubmatrix"))
    Matrix().UpdateSubmatrix( rowIndLoc, colIndLoc, alpha, ASub );
}

template<typename T>
void 
AbstractBlockDistMatrix<T>::UpdateRealPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  Base<T> alpha, const El::Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::UpdateRealPartOfLocalSubmatrix"))
    Matrix().UpdateRealPartOfSubmatrix( rowIndLoc, colIndLoc, alpha, ASub );
}

template<typename T>
void 
AbstractBlockDistMatrix<T>::UpdateImagPartOfLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
  Base<T> alpha, const El::Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::UpdateImagPartOfLocalSubmatrix"))
    Matrix().UpdateImagPartOfSubmatrix( rowIndLoc, colIndLoc, alpha, ASub );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::MakeLocalSubmatrixReal
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::MakeLocalSubmatrixReal"))
    Matrix().MakeSubmatrixReal( rowIndLoc, colIndLoc );
}

template<typename T>
void
AbstractBlockDistMatrix<T>::ConjugateLocalSubmatrix
( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::ConjugateLocalSubmatrix"))
    Matrix().ConjugateSubmatrix( rowIndLoc, colIndLoc );
}

// Sum the local matrix over a particular communicator
// ===================================================
// NOTE: The matrix dimensions *must* be uniform over the communicator.

template<typename T>
void
AbstractBlockDistMatrix<T>::SumOver( mpi::Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::SumOver"))
    if( !Participating() )
        return;

    const Int localHeight = LocalHeight();
    const Int localWidth = LocalWidth();
    const Int localSize = mpi::Pad( localHeight*localWidth );
    T* sumBuf = auxMemory_.Require( localSize );   

    // Pack
    T* buf = Buffer();
    const Int ldim = LDim(); 
    EL_PARALLEL_FOR
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const T* thisCol = &buf[jLoc*ldim];
        T* sumCol = &sumBuf[jLoc*localHeight];
        MemCopy( sumCol, thisCol, localHeight );
    }

    // AllReduce sum
    mpi::AllReduce( sumBuf, localSize, comm );

    // Unpack
    EL_PARALLEL_FOR
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
AbstractBlockDistMatrix<T>::ComplainIfReal() const
{ 
    if( !IsComplex<T>::val )
        LogicError("Called complex-only routine with real data");
}

template<typename T>
void
AbstractBlockDistMatrix<T>::AssertNotLocked() const
{
    if( Locked() )
        LogicError("Assertion that matrix not be a locked view failed");
}

template<typename T>
void
AbstractBlockDistMatrix<T>::AssertNotStoringData() const
{
    if( matrix_.MemorySize() > 0 )
        LogicError("Assertion that matrix not be storing data failed");
}

template<typename T>
void
AbstractBlockDistMatrix<T>::AssertValidEntry( Int i, Int j ) const
{
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
        LogicError
        ("Entry (",i,",",j,") is out of bounds of ",Height(),
         " x ",Width()," matrix");
}

template<typename T>
void
AbstractBlockDistMatrix<T>::AssertValidSubmatrix
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
AbstractBlockDistMatrix<T>::AssertSameGrid( const El::Grid& grid ) const
{
    if( Grid() != grid )
        LogicError("Assertion that grids match failed");
}

template<typename T> 
void
AbstractBlockDistMatrix<T>::AssertSameSize( Int height, Int width ) const
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
AbstractBlockDistMatrix<T>::ShallowSwap( AbstractBlockDistMatrix<T>& A )
{
    matrix_.ShallowSwap( A.matrix_ );
    auxMemory_.ShallowSwap( A.auxMemory_ );
    std::swap( viewType_, A.viewType_ );
    std::swap( height_ , A.height_ );
    std::swap( width_, A.width_ );
    std::swap( colConstrained_, A.colConstrained_ );
    std::swap( rowConstrained_, A.rowConstrained_ );
    std::swap( rootConstrained_, A.rootConstrained_ );
    std::swap( blockHeight_, A.blockHeight_ );
    std::swap( blockWidth_, A.blockWidth_ );
    std::swap( colAlign_, A.colAlign_ );
    std::swap( rowAlign_, A.rowAlign_ );
    std::swap( colCut_, A.colCut_ );
    std::swap( rowCut_, A.rowCut_ );
    std::swap( colShift_, A.colShift_ );
    std::swap( rowShift_, A.rowShift_ );
    std::swap( root_, A.root_ );
    std::swap( grid_, A.grid_ );
}

// Modify the distribution metadata
// ================================

template<typename T>
void
AbstractBlockDistMatrix<T>::SetShifts()
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
AbstractBlockDistMatrix<T>::SetColShift()
{
    if( Participating() )
        colShift_ = Shift(ColRank(),colAlign_,ColStride());
    else
        colShift_ = 0;
}

template<typename T>
void
AbstractBlockDistMatrix<T>::SetRowShift()
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
( const AbstractBlockDistMatrix<T>& AL, const AbstractBlockDistMatrix<T>& AR )
{
    if( AL.Height() != AR.Height() )    
        LogicError
        ("1x2 not conformant:\n",
         DimsString(AL,"Left"),"\n",DimsString(AR,"Right"));
    if( AL.ColAlign() != AR.ColAlign() || AL.ColCut() != AR.ColCut() )
        LogicError("1x2 is misaligned");
}

template<typename T> 
void
AssertConforming2x1
( const AbstractBlockDistMatrix<T>& AT, const AbstractBlockDistMatrix<T>& AB )
{
    if( AT.Width() != AB.Width() )
        LogicError
        ("2x1 is not conformant:\n",
         DimsString(AT,"Top"),"\n",DimsString(AB,"Bottom"));
    if( AT.RowAlign() != AB.RowAlign() || AT.RowCut() != AB.RowCut() )
        LogicError("2x1 is not aligned");
}

template<typename T> 
void
AssertConforming2x2
( const AbstractBlockDistMatrix<T>& ATL, const AbstractBlockDistMatrix<T>& ATR, 
  const AbstractBlockDistMatrix<T>& ABL, const AbstractBlockDistMatrix<T>& ABR )
{
    if( ATL.Width() != ABL.Width() || ATR.Width() != ABR.Width() ||
        ATL.Height() != ATR.Height() || ABL.Height() != ABR.Height() )
        LogicError
        ("2x2 is not conformant:\n",
         DimsString(ATL,"TL"),"\n",DimsString(ATR,"TR"),"\n",
         DimsString(ABL,"BL"),"\n",DimsString(ABR,"BR"));
    if( ATL.ColAlign() != ATR.ColAlign() || ATL.ColCut() != ATR.ColCut() ||
        ABL.ColAlign() != ABR.ColAlign() || ABL.ColCut() != ABR.ColCut() ||
        ATL.RowAlign() != ABL.RowAlign() || ATL.RowCut() != ABL.RowCut() ||
        ATR.RowAlign() != ABR.RowAlign() || ATR.RowCut() != ABR.RowCut() )
        LogicError("2x2 set of matrices must aligned to combine");
}

// Instantiations for {Int,Real,Complex<Real>} for each Real in {float,double}
// ###########################################################################

#ifndef EL_RELEASE
 #define PROTO(T) \
  template class AbstractBlockDistMatrix<T>;\
  template void AssertConforming1x2\
  ( const AbstractBlockDistMatrix<T>& AL,  \
    const AbstractBlockDistMatrix<T>& AR );\
  template void AssertConforming2x1\
  ( const AbstractBlockDistMatrix<T>& AT,  \
    const AbstractBlockDistMatrix<T>& AB );\
  template void AssertConforming2x2\
  ( const AbstractBlockDistMatrix<T>& ATL, \
    const AbstractBlockDistMatrix<T>& ATR, \
    const AbstractBlockDistMatrix<T>& ABL, \
    const AbstractBlockDistMatrix<T>& ABR )
#else
 #define PROTO(T) template class AbstractBlockDistMatrix<T>
#endif
 
PROTO(Int);
#ifndef EL_DISABLE_FLOAT
PROTO(float);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef EL_DISABLE_COMPLEX
#endif // ifndef EL_DISABLE_FLOAT

PROTO(double);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El
