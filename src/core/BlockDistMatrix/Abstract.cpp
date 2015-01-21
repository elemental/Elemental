/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
AbstractBlockDistMatrix<T>::AbstractBlockDistMatrix
( const El::Grid& g, int root )
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
( const El::Grid& g, Int blockHeight, Int blockWidth, int root )
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
    if( Viewing() || A.Viewing() )
    {
        Copy( A, *this );
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
void AbstractBlockDistMatrix<T>::Empty()
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
void AbstractBlockDistMatrix<T>::EmptyData()
{
    matrix_.Empty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
}

template<typename T>
void AbstractBlockDistMatrix<T>::SetGrid( const El::Grid& grid )
{
    if( grid_ != &grid )
    {
        Empty();
        grid_ = &grid; 
        SetShifts();
    }
}

template<typename T>
void AbstractBlockDistMatrix<T>::Resize( Int height, Int width )
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
void AbstractBlockDistMatrix<T>::Resize( Int height, Int width, Int ldim )
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
void AbstractBlockDistMatrix<T>::MakeConsistent( bool includingViewers )
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
    const int root                = message[12];

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
void AbstractBlockDistMatrix<T>::MakeSizeConsistent( bool includingViewers )
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
void AbstractBlockDistMatrix<T>::Align
( Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut, bool constrain )
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
void AbstractBlockDistMatrix<T>::AlignCols
( Int blockHeight, int colAlign, Int colCut, bool constrain )
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
void AbstractBlockDistMatrix<T>::AlignRows
( Int blockWidth, int rowAlign, Int rowCut, bool constrain )
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
void AbstractBlockDistMatrix<T>::FreeAlignments() 
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
void AbstractBlockDistMatrix<T>::SetRoot( int root, bool constrain )
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
void AbstractBlockDistMatrix<T>::AlignWith
( const El::BlockDistData& data, bool constrain )
{ 
    DEBUG_ONLY(CallStackEntry cse("ABDM::AlignWith"))
    AlignColsWith( data, constrain );
    AlignRowsWith( data, constrain );
}

template<typename T>
void AbstractBlockDistMatrix<T>::AlignAndResize
( Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut,
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
void AbstractBlockDistMatrix<T>::AlignColsAndResize
( Int blockHeight, int colAlign, Int colCut, Int height, Int width, 
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
void AbstractBlockDistMatrix<T>::AlignRowsAndResize
( Int blockWidth, int rowAlign, Int rowCut, Int height, Int width, 
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
void AbstractBlockDistMatrix<T>::Attach
( Int height, Int width, const El::Grid& g, 
  Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut,
  T* buffer, Int ldim, int root )
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
void AbstractBlockDistMatrix<T>::Attach
( Int height, Int width, const El::Grid& g,
  Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut, El::Matrix<T>& A, 
  int root )
{
    // TODO: Assert that the local dimensions are correct
    Attach
    ( height, width, g, blockHeight, blockWidth, 
      colAlign, rowAlign, colCut, rowCut, A.Buffer(), A.LDim(), root );
}

template<typename T>
void AbstractBlockDistMatrix<T>::LockedAttach
( Int height, Int width, const El::Grid& g, 
  Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut,
  const T* buffer, Int ldim, int root )
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
void AbstractBlockDistMatrix<T>::LockedAttach
( Int height, Int width, const El::Grid& g,
  Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut, const El::Matrix<T>& A,
  int root )
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
int AbstractBlockDistMatrix<T>::ColAlign() const { return colAlign_; }
template<typename T>
int AbstractBlockDistMatrix<T>::RowAlign() const { return rowAlign_; }

template<typename T>
Int AbstractBlockDistMatrix<T>::ColCut() const { return colCut_; }
template<typename T>
Int AbstractBlockDistMatrix<T>::RowCut() const { return rowCut_; }

template<typename T>
int AbstractBlockDistMatrix<T>::ColShift() const { return colShift_; }
template<typename T>
int AbstractBlockDistMatrix<T>::RowShift() const { return rowShift_; }

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
int AbstractBlockDistMatrix<T>::RowOwner( Int i ) const
{ return int((((i+ColCut())/BlockHeight())+ColAlign()) % ColStride()); }
template<typename T>
int AbstractBlockDistMatrix<T>::ColOwner( Int j ) const
{ return int((((j+RowCut())/BlockWidth())+RowAlign()) % RowStride()); }
template<typename T>
int AbstractBlockDistMatrix<T>::Owner( Int i, Int j ) const
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
int AbstractBlockDistMatrix<T>::PartialColStride() const { return ColStride(); }
template<typename T>
int AbstractBlockDistMatrix<T>::PartialRowStride() const { return RowStride(); }

template<typename T>
int AbstractBlockDistMatrix<T>::PartialUnionColStride() const { return 1; }
template<typename T>
int AbstractBlockDistMatrix<T>::PartialUnionRowStride() const { return 1; }

template<typename T>
int AbstractBlockDistMatrix<T>::ColRank() const { return mpi::Rank(ColComm()); }
template<typename T>
int AbstractBlockDistMatrix<T>::RowRank() const { return mpi::Rank(RowComm()); }

template<typename T>
int AbstractBlockDistMatrix<T>::PartialColRank() const
{ return mpi::Rank(PartialColComm()); }
template<typename T>
int AbstractBlockDistMatrix<T>::PartialRowRank() const
{ return mpi::Rank(PartialRowComm()); }

template<typename T>
int AbstractBlockDistMatrix<T>::PartialUnionColRank() const
{ return mpi::Rank(PartialUnionColComm()); }
template<typename T>
int AbstractBlockDistMatrix<T>::PartialUnionRowRank() const
{ return mpi::Rank(PartialUnionRowComm()); }

template<typename T>
int AbstractBlockDistMatrix<T>::DistRank() const
{ return mpi::Rank(DistComm()); }
template<typename T>
int AbstractBlockDistMatrix<T>::CrossRank() const
{ return mpi::Rank(CrossComm()); }
template<typename T>
int AbstractBlockDistMatrix<T>::RedundantRank() const
{ return mpi::Rank(RedundantComm()); }

template<typename T>
int AbstractBlockDistMatrix<T>::Root() const { return root_; }

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
        const int owner = Owner( i, j );
        if( owner == DistRank() )
            value = GetLocal( LocalRow(i), LocalCol(j) );
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
        const int owner = Owner( i, j );
        if( owner == DistRank() )
            value = GetLocalRealPart( LocalRow(i), LocalCol(j) );
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
            const int owner = Owner( i, j );
            if( owner == DistRank() )
                value = GetLocalRealPart( LocalRow(i), LocalCol(j) );
            mpi::Broadcast( value, owner, DistComm() );
        }
        mpi::Broadcast( value, Root(), CrossComm() );
    }
    else
        value = 0;
    return value;
}

template<typename T>
void AbstractBlockDistMatrix<T>::Set( Int i, Int j, T value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::Set"))
    if( IsLocal(i,j) )
        SetLocal( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void AbstractBlockDistMatrix<T>::SetRealPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::SetRealPart"))
    if( IsLocal(i,j) )
        SetLocalRealPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void AbstractBlockDistMatrix<T>::SetImagPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::SetImagPart"))
    if( IsLocal(i,j) )
        SetLocalImagPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void AbstractBlockDistMatrix<T>::Update( Int i, Int j, T value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::Update"))
    if( IsLocal(i,j) )
        UpdateLocal( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void AbstractBlockDistMatrix<T>::UpdateRealPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::UpdateRealPart"))
    if( IsLocal(i,j) )
        UpdateLocalRealPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void AbstractBlockDistMatrix<T>::UpdateImagPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::UpdateImagPart"))
    if( IsLocal(i,j) )
        UpdateLocalImagPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void AbstractBlockDistMatrix<T>::MakeReal( Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::MakeReal"))
    if( IsLocal(i,j) )
        MakeLocalReal( LocalRow(i), LocalCol(j) );
}

template<typename T>
void AbstractBlockDistMatrix<T>::Conjugate( Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::Conjugate"))
    if( IsLocal(i,j) )
        ConjugateLocal( LocalRow(i), LocalCol(j) );
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
void AbstractBlockDistMatrix<T>::SetLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Set(iLoc,jLoc,alpha); }

template<typename T>
void AbstractBlockDistMatrix<T>::SetLocalRealPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.SetRealPart(iLoc,jLoc,alpha); }

template<typename T>
void AbstractBlockDistMatrix<T>::SetLocalImagPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.SetImagPart(iLoc,jLoc,alpha); }

template<typename T>
void AbstractBlockDistMatrix<T>::UpdateLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Update(iLoc,jLoc,alpha); }

template<typename T>
void AbstractBlockDistMatrix<T>::UpdateLocalRealPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.UpdateRealPart(iLoc,jLoc,alpha); }

template<typename T>
void AbstractBlockDistMatrix<T>::UpdateLocalImagPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.UpdateImagPart(iLoc,jLoc,alpha); }

template<typename T>
void AbstractBlockDistMatrix<T>::MakeLocalReal( Int iLoc, Int jLoc )
{ matrix_.MakeReal( iLoc, jLoc ); }

template<typename T>
void AbstractBlockDistMatrix<T>::ConjugateLocal( Int iLoc, Int jLoc )
{ matrix_.Conjugate( iLoc, jLoc ); }

// Diagonal manipulation
// =====================
template<typename T>
bool AbstractBlockDistMatrix<T>::DiagonalAlignedWith
( const El::BlockDistData& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::DiagonalAlignedWith"))
    // TODO: Ensure blocksize is compatible...the blocksizes needed for a 
    //       diagonal distribution are variable except for special cases.
    LogicError("This routine is not yet written");
    return false;
}

template<typename T>
int AbstractBlockDistMatrix<T>::DiagonalRoot( Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::DiagonalRoot"))
    LogicError("This routine is not yet written");
    return 0;
}

template<typename T>
int AbstractBlockDistMatrix<T>::DiagonalAlign( Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::DiagonalAlign"))
    LogicError("This routine is not yet written");
    return 0;
}

template<typename T>
void AbstractBlockDistMatrix<T>::AlignColsWith
( const El::BlockDistData& data, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::AlignColsWith"))
    SetGrid( *data.grid );
    SetRoot( data.root );
    if( data.colDist == ColDist() || data.colDist == PartialColDist() )
        AlignCols( data.blockHeight, data.colAlign, data.colCut, constrain );
    else if( data.rowDist == ColDist() || data.rowDist == PartialColDist() )
        AlignCols( data.blockWidth, data.rowAlign, data.rowCut, constrain );
    else if( data.colDist == PartialUnionColDist() )
        AlignCols
        ( data.blockHeight, data.colAlign % ColStride(), data.colCut,
          constrain );
    else if( data.rowDist == PartialUnionColDist() )
        AlignCols
        ( data.blockWidth, data.rowAlign % ColStride(), data.rowCut,
          constrain );
    DEBUG_ONLY(
        else if( ColDist()    != CollectedColDist() && 
                 data.colDist != CollectedColDist() && 
                 data.rowDist != CollectedColDist() )
            LogicError("Nonsensical alignment");
    )
}

template<typename T>
void AbstractBlockDistMatrix<T>::AlignRowsWith
( const El::BlockDistData& data, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("ABDM::AlignRowsWith"))
    SetGrid( *data.grid );
    SetRoot( data.root );
    if( data.colDist == RowDist() || data.colDist == PartialRowDist() )
        AlignRows( data.blockHeight, data.colAlign, data.colCut, constrain );
    else if( data.rowDist == RowDist() || data.rowDist == PartialRowDist() )
        AlignRows( data.blockWidth, data.rowAlign, data.rowCut, constrain );
    else if( data.colDist == PartialUnionRowDist() )
        AlignRows
        ( data.blockHeight, data.colAlign % RowStride(), data.colCut,
          constrain );
    else if( data.rowDist == PartialUnionRowDist() )
        AlignRows
        ( data.blockWidth, data.rowAlign % RowStride(), data.rowCut,
          constrain );
    DEBUG_ONLY(
        else if( RowDist()    != CollectedRowDist() && 
                 data.colDist != CollectedRowDist() && 
                 data.rowDist != CollectedRowDist() )
            LogicError("Nonsensical alignment");
    )
}

// Assertions
// ==========

template<typename T>
void AbstractBlockDistMatrix<T>::ComplainIfReal() const
{ 
    if( !IsComplex<T>::val )
        LogicError("Called complex-only routine with real data");
}

template<typename T>
void AbstractBlockDistMatrix<T>::AssertNotLocked() const
{
    if( Locked() )
        LogicError("Assertion that matrix not be a locked view failed");
}

template<typename T>
void AbstractBlockDistMatrix<T>::AssertNotStoringData() const
{
    if( matrix_.MemorySize() > 0 )
        LogicError("Assertion that matrix not be storing data failed");
}

template<typename T>
void AbstractBlockDistMatrix<T>::AssertValidEntry( Int i, Int j ) const
{
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
        LogicError
        ("Entry (",i,",",j,") is out of bounds of ",Height(),
         " x ",Width()," matrix");
}

template<typename T>
void AbstractBlockDistMatrix<T>::AssertValidSubmatrix
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
void AbstractBlockDistMatrix<T>::AssertSameSize( Int height, Int width ) const
{
    if( Height() != height || Width() != width )
        LogicError("Assertion that matrices be the same size failed");
}

// Private section
// ###############

// Exchange metadata with another matrix
// =====================================

template<typename T>
void AbstractBlockDistMatrix<T>::ShallowSwap( AbstractBlockDistMatrix<T>& A )
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
void AbstractBlockDistMatrix<T>::SetShifts()
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
void AbstractBlockDistMatrix<T>::SetColShift()
{
    if( Participating() )
        colShift_ = Shift(ColRank(),colAlign_,ColStride());
    else
        colShift_ = 0;
}

template<typename T>
void AbstractBlockDistMatrix<T>::SetRowShift()
{
    if( Participating() )
        rowShift_ = Shift(RowRank(),rowAlign_,RowStride());
    else
        rowShift_ = 0;
}

// Outside of class
// ----------------

template<typename T> 
void AssertConforming1x2
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
void AssertConforming2x1
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
void AssertConforming2x2
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
    const AbstractBlockDistMatrix<T>& ABR );
#else
 #define PROTO(T) template class AbstractBlockDistMatrix<T>;
#endif

#include "El/macros/Instantiate.h"

} // namespace El
