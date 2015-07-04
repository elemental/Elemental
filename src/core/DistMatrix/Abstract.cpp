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

template<typename Ring>
AbstractDistMatrix<Ring>::AbstractDistMatrix( const El::Grid& grid, int root )
: viewType_(OWNER),
  height_(0), width_(0),
  matrix_(0,0,true),
  colConstrained_(false), rowConstrained_(false), rootConstrained_(false),
  colAlign_(0), rowAlign_(0),
  root_(root), grid_(&grid)
{ }

template<typename Ring>
AbstractDistMatrix<Ring>::AbstractDistMatrix( AbstractDistMatrix<Ring>&& A ) 
EL_NOEXCEPT
: viewType_(A.viewType_),
  height_(A.height_), width_(A.width_), 
  colConstrained_(A.colConstrained_), rowConstrained_(A.rowConstrained_),
  rootConstrained_(A.rootConstrained_),
  colAlign_(A.colAlign_), rowAlign_(A.rowAlign_),
  colShift_(A.colShift_), rowShift_(A.rowShift_), 
  root_(A.root_),
  grid_(A.grid_)
{ matrix_.ShallowSwap( A.matrix_ ); }

// Optional to override
// --------------------

template<typename Ring>
AbstractDistMatrix<Ring>::~AbstractDistMatrix() { }

// Assignment and reconfiguration
// ==============================
template<typename Ring>
void
AbstractDistMatrix<Ring>::Empty()
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::EmptyData()
{
    matrix_.Empty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetGrid( const El::Grid& grid )
{
    if( grid_ != &grid )
    {
        grid_ = &grid; 
        Empty();
    }
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::Resize( Int height, Int width )
{
    DEBUG_ONLY(
      CSE cse("ADM::Resize");
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::Resize( Int height, Int width, Int ldim )
{
    DEBUG_ONLY(
      CSE cse("ADM::Resize");
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::MakeConsistent( bool includingViewers )
{
    DEBUG_ONLY(CSE cse("ADM::MakeConsistent"))

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

    const auto& g = *grid_;
    if( !g.InGrid() && !includingViewers )
        LogicError("Non-participating process called MakeConsistent");
    if( g.InGrid() )
    {
        // TODO: Ensure roots are consistent within each cross communicator
        mpi::Broadcast( message, msgLength, Root(), CrossComm() );
    }
    if( includingViewers )
    {
        const Int vcRoot = g.VCToViewing(0);
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
    const int root                = message[8];

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

template<typename Ring>
void
AbstractDistMatrix<Ring>::MakeSizeConsistent( bool includingViewers )
{
    DEBUG_ONLY(CSE cse("ADM::MakeSizeConsistent"))

    const Int msgSize = 2;
    Int message[msgSize];
    if( CrossRank() == Root() )
    {
        message[0] = height_;
        message[1] = width_;
    }

    const auto& g = *grid_;
    if( !g.InGrid() && !includingViewers )
        LogicError("Non-participating process called MakeSizeConsistent");
    if( g.InGrid() )
        mpi::Broadcast( message, msgSize, Root(), CrossComm() );
    if( includingViewers )
    {
        const Int vcRoot = g.VCToViewing(0);
        mpi::Broadcast( message, msgSize, vcRoot, g.ViewingComm() );
    }
    const Int newHeight = message[0]; 
    const Int newWidth  = message[1];
    Resize( newHeight, newWidth );
}

// Realignment
// -----------

template<typename Ring>
void
AbstractDistMatrix<Ring>::Align( int colAlign, int rowAlign, bool constrain )
{ 
    DEBUG_ONLY(CSE cse("ADM::Align"))
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::AlignCols( int colAlign, bool constrain )
{ 
    DEBUG_ONLY(
      CSE cse("ADM::AlignCols");
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::AlignRows( int rowAlign, bool constrain )
{ 
    DEBUG_ONLY(
      CSE cse("ADM::AlignRows");
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::FreeAlignments() 
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetRoot( int root, bool constrain )
{
    DEBUG_ONLY(
      CSE cse("ADM::SetRoot");
      if( root < 0 || root >= mpi::Size(CrossComm()) )
          LogicError("Invalid root");
    )
    if( root != root_ )
        Empty();
    root_ = root;
    if( constrain )
        rootConstrained_ = true;
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::AlignWith
( const El::DistData& data, bool constrain, bool allowMismatch )
{ 
    DEBUG_ONLY(CSE cse("ADM::AlignWith"))
    AlignColsWith( data, constrain, allowMismatch );
    AlignRowsWith( data, constrain, allowMismatch );
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::AlignColsWith
( const El::DistData& data, bool constrain, bool allowMismatch )
{
    DEBUG_ONLY(CSE cse("ADM::AlignColsWith"))
    SetGrid( *data.grid );
    SetRoot( data.root );
    if(      data.colDist == ColDist() || data.colDist == PartialColDist() )
        AlignCols( data.colAlign, constrain );
    else if( data.rowDist == ColDist() || data.rowDist == PartialColDist() )
        AlignCols( data.rowAlign, constrain );
    else if( data.colDist == PartialUnionColDist() )
        AlignCols( data.colAlign % ColStride(), constrain );
    else if( data.rowDist == PartialUnionColDist() )
        AlignCols( data.rowAlign % ColStride(), constrain );
    else if( ColDist()    != CollectedColDist() && 
             data.colDist != CollectedColDist() && 
             data.rowDist != CollectedColDist() && !allowMismatch )
        LogicError("Nonsensical alignment");
}

template<typename Ring>
void AbstractDistMatrix<Ring>::AlignRowsWith
( const El::DistData& data, bool constrain, bool allowMismatch )
{
    DEBUG_ONLY(CSE cse("ADM::AlignRowsWith"))
    SetGrid( *data.grid );
    SetRoot( data.root );
    if(      data.colDist == RowDist() || data.colDist == PartialRowDist() )
        AlignRows( data.colAlign, constrain );
    else if( data.rowDist == RowDist() || data.rowDist == PartialRowDist() )
        AlignRows( data.rowAlign, constrain );
    else if( data.colDist == PartialUnionRowDist() )
        AlignRows( data.colAlign % RowStride(), constrain );
    else if( data.rowDist == PartialUnionRowDist() )
        AlignRows( data.rowAlign % RowStride(), constrain );
    else if( RowDist()    != CollectedRowDist() && 
             data.colDist != CollectedRowDist() && 
             data.rowDist != CollectedRowDist() && !allowMismatch )
        LogicError("Nonsensical alignment");
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::AlignAndResize
( int colAlign, int rowAlign, Int height, Int width, 
  bool force, bool constrain )
{
    DEBUG_ONLY(CSE cse("ADM::AlignAndResize"))
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::AlignColsAndResize
( int colAlign, Int height, Int width, bool force, bool constrain )
{
    DEBUG_ONLY(CSE cse("ADM::AlignColsAndResize"))
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::AlignRowsAndResize
( int rowAlign, Int height, Int width, bool force, bool constrain )
{
    DEBUG_ONLY(CSE cse("ADM::AlignRowsAndResize"))
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::Attach
( Int height, Int width, const El::Grid& g, 
  int colAlign, int rowAlign, Ring* buffer, Int ldim, int root )
{
    DEBUG_ONLY(CSE cse("ADM::Attach"))
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::Attach
( Int height, Int width, const El::Grid& g,
  int colAlign, int rowAlign, El::Matrix<Ring>& A, int root )
{
    // TODO: Assert that the local dimensions are correct
    Attach( height, width, g, colAlign, rowAlign, A.Buffer(), A.LDim(), root );
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::Attach( const El::Grid& g, El::Matrix<Ring>& A )
{
    DEBUG_ONLY(CSE cse("ADM::Attach"))
    if( g.Size() != 1 )
        LogicError("Assumed a grid size of one");
    Attach( A.Height(), A.Width(), g, 0, 0, A.Buffer(), A.LDim() );
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::LockedAttach
( Int height, Int width, const El::Grid& g, 
  int colAlign, int rowAlign, const Ring* buffer, Int ldim, int root )
{
    DEBUG_ONLY(CSE cse("ADM::LockedAttach"))
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::LockedAttach
( Int height, Int width, const El::Grid& g, 
  int colAlign, int rowAlign, const El::Matrix<Ring>& A, int root )
{
    // TODO: Assert that the local dimensions are correct
    LockedAttach
    ( height, width, g, colAlign, rowAlign, A.LockedBuffer(), A.LDim(), root );
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::LockedAttach( const El::Grid& g, const El::Matrix<Ring>& A )
{
    DEBUG_ONLY(CSE cse("ADM::LockedAttach"))
    if( g.Size() != 1 )
        LogicError("Assumed a grid size of one");
    LockedAttach( A.Height(), A.Width(), g, 0, 0, A.LockedBuffer(), A.LDim() );
}

// Operator overloading
// ====================

// Copy
// ----
template<typename Ring>
const AbstractDistMatrix<Ring>&
AbstractDistMatrix<Ring>::operator=( const AbstractDistMatrix<Ring>& A )
{
    DEBUG_ONLY(CSE cse("ADM::operator=(ADM&)"))
    Copy( A, *this );
    return *this;
}

template<typename Ring>
const AbstractDistMatrix<Ring>&
AbstractDistMatrix<Ring>::operator=( const DistMultiVec<Ring>& A )
{
    DEBUG_ONLY(CSE cse("ADM::operator=(DMV&)"))
    Copy( A, *this );
    return *this;
}

// Move assignment
// ---------------
template<typename Ring>
AbstractDistMatrix<Ring>& 
AbstractDistMatrix<Ring>::operator=( AbstractDistMatrix<Ring>&& A )
{
    DEBUG_ONLY(CSE cse("ADM::operator=(ADM&&)"))
    if( Viewing() || A.Viewing() )
    {
        Copy( A, *this );
    }
    else
    {
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

// Rescaling
// ---------
template<typename Ring>
const AbstractDistMatrix<Ring>&
AbstractDistMatrix<Ring>::operator*=( Ring alpha )
{
    DEBUG_ONLY(CSE cse("ADM::operator*=(T)"))    
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename Ring>
const AbstractDistMatrix<Ring>&
AbstractDistMatrix<Ring>::operator+=( const AbstractDistMatrix<Ring>& A )
{
    DEBUG_ONLY(CSE cse("ADM::operator+="))
    Axpy( Ring(1), A, *this );
    return *this;
}

template<typename Ring>
const AbstractDistMatrix<Ring>&
AbstractDistMatrix<Ring>::operator-=( const AbstractDistMatrix<Ring>& A )
{
    DEBUG_ONLY(CSE cse("ADM::operator-="))
    Axpy( Ring(-1), A, *this );
    return *this;
}

// Basic queries
// =============

// Global matrix information
// -------------------------

template<typename Ring>
Int AbstractDistMatrix<Ring>::Height() const { return height_; }
template<typename Ring>
Int AbstractDistMatrix<Ring>::Width() const { return width_; }

template<typename Ring>
Int AbstractDistMatrix<Ring>::DiagonalLength( Int offset ) const
{ return El::DiagonalLength(height_,width_,offset); }

template<typename Ring>
bool AbstractDistMatrix<Ring>::Viewing() const { return IsViewing( viewType_ ); }
template<typename Ring>
bool AbstractDistMatrix<Ring>::Locked() const { return IsLocked( viewType_ ); }

// Local matrix information
// ------------------------

template<typename Ring>
Int AbstractDistMatrix<Ring>::LocalHeight() const { return matrix_.Height(); }
template<typename Ring>
Int AbstractDistMatrix<Ring>::LocalWidth() const { return matrix_.Width(); }
template<typename Ring>
Int AbstractDistMatrix<Ring>::LDim() const { return matrix_.LDim(); }

template<typename Ring>
El::Matrix<Ring>& 
AbstractDistMatrix<Ring>::Matrix() { return matrix_; }
template<typename Ring>
const El::Matrix<Ring>& 
AbstractDistMatrix<Ring>::LockedMatrix() const { return matrix_; }

template<typename Ring>
size_t
AbstractDistMatrix<Ring>::AllocatedMemory() const { return matrix_.MemorySize(); }

template<typename Ring>
Ring*
AbstractDistMatrix<Ring>::Buffer() { return matrix_.Buffer(); }

template<typename Ring>
Ring*
AbstractDistMatrix<Ring>::Buffer( Int iLoc, Int jLoc )
{ return matrix_.Buffer(iLoc,jLoc); }

template<typename Ring>
const Ring*
AbstractDistMatrix<Ring>::LockedBuffer() const
{ return matrix_.LockedBuffer(); }

template<typename Ring>
const Ring*
AbstractDistMatrix<Ring>::LockedBuffer( Int iLoc, Int jLoc ) const
{ return matrix_.LockedBuffer(iLoc,jLoc); }

// Distribution information
// ------------------------

template<typename Ring>
const El::Grid& AbstractDistMatrix<Ring>::Grid() const { return *grid_; }

template<typename Ring>
int AbstractDistMatrix<Ring>::ColAlign() const { return colAlign_; }
template<typename Ring>
int AbstractDistMatrix<Ring>::RowAlign() const { return rowAlign_; }

template<typename Ring>
int AbstractDistMatrix<Ring>::ColShift() const { return colShift_; }
template<typename Ring>
int AbstractDistMatrix<Ring>::RowShift() const { return rowShift_; }

template<typename Ring>
bool AbstractDistMatrix<Ring>::ColConstrained() const { return colConstrained_; }
template<typename Ring>
bool AbstractDistMatrix<Ring>::RowConstrained() const { return rowConstrained_; }
template<typename Ring>
bool AbstractDistMatrix<Ring>::RootConstrained() const { return rootConstrained_; }

template<typename Ring>
bool AbstractDistMatrix<Ring>::Participating() const
{ return grid_->InGrid() && (CrossRank()==root_); }

template<typename Ring>
int AbstractDistMatrix<Ring>::RowOwner( Int i ) const
{
    if( i == END ) i = height_ - 1;
    return int((i+ColAlign()) % ColStride());
}

template<typename Ring>
int AbstractDistMatrix<Ring>::ColOwner( Int j ) const
{ 
    if( j == END ) j = width_ - 1;
    return int((j+RowAlign()) % RowStride()); 
}

template<typename Ring>
int AbstractDistMatrix<Ring>::Owner( Int i, Int j ) const
{ return RowOwner(i)+ColOwner(j)*ColStride(); }

template<typename Ring>
Int AbstractDistMatrix<Ring>::LocalRow( Int i ) const
{ 
    DEBUG_ONLY(
      CSE cse("ADM::LocalRow");
      if( !IsLocalRow(i) )
          LogicError("Requested local index of non-local row");
    )
    return LocalRowOffset(i);
}

template<typename Ring>
Int AbstractDistMatrix<Ring>::LocalCol( Int j ) const
{
    DEBUG_ONLY(
      CSE cse("ADM::LocalCol");
      if( !IsLocalCol(j) )
          LogicError("Requested local index of non-local column");
    )
    return LocalColOffset(j);
}

template<typename Ring>
Int AbstractDistMatrix<Ring>::LocalRowOffset( Int i ) const
{ 
    if( i == END ) i = height_ - 1;
    return Length_(i,ColShift(),ColStride()); 
}

template<typename Ring>
Int AbstractDistMatrix<Ring>::LocalColOffset( Int j ) const
{ 
    if( j == END ) j = width_ - 1;
    return Length_(j,RowShift(),RowStride()); 
}

template<typename Ring>
Int AbstractDistMatrix<Ring>::GlobalRow( Int iLoc ) const
{ 
    if( iLoc == END ) iLoc = LocalHeight() - 1;
    return ColShift() + iLoc*ColStride(); 
}

template<typename Ring>
Int AbstractDistMatrix<Ring>::GlobalCol( Int jLoc ) const
{ 
    if( jLoc == END ) jLoc = LocalWidth() - 1;
    return RowShift() + jLoc*RowStride(); 
}

template<typename Ring>
bool AbstractDistMatrix<Ring>::IsLocalRow( Int i ) const
{ return Participating() && RowOwner(i) == ColRank(); }
template<typename Ring>
bool AbstractDistMatrix<Ring>::IsLocalCol( Int j ) const
{ return Participating() && ColOwner(j) == RowRank(); }
template<typename Ring>
bool AbstractDistMatrix<Ring>::IsLocal( Int i, Int j ) const
{ return IsLocalRow(i) && IsLocalCol(j); }

template<typename Ring>
mpi::Comm AbstractDistMatrix<Ring>::PartialColComm() const { return ColComm(); }
template<typename Ring>
mpi::Comm AbstractDistMatrix<Ring>::PartialRowComm() const { return RowComm(); }

template<typename Ring>
mpi::Comm AbstractDistMatrix<Ring>::PartialUnionColComm() const
{ return mpi::COMM_SELF; }
template<typename Ring>
mpi::Comm AbstractDistMatrix<Ring>::PartialUnionRowComm() const
{ return mpi::COMM_SELF; }

template<typename Ring>
int AbstractDistMatrix<Ring>::PartialColStride() const { return ColStride(); }
template<typename Ring>
int AbstractDistMatrix<Ring>::PartialRowStride() const { return RowStride(); }

template<typename Ring>
int AbstractDistMatrix<Ring>::PartialUnionColStride() const { return 1; }
template<typename Ring>
int AbstractDistMatrix<Ring>::PartialUnionRowStride() const { return 1; }

template<typename Ring>
int AbstractDistMatrix<Ring>::ColRank() const { return mpi::Rank(ColComm()); }
template<typename Ring>
int AbstractDistMatrix<Ring>::RowRank() const { return mpi::Rank(RowComm()); }

template<typename Ring>
int AbstractDistMatrix<Ring>::PartialColRank() const
{ return mpi::Rank(PartialColComm()); }
template<typename Ring>
int AbstractDistMatrix<Ring>::PartialRowRank() const
{ return mpi::Rank(PartialRowComm()); }

template<typename Ring>
int AbstractDistMatrix<Ring>::PartialUnionColRank() const
{ return mpi::Rank(PartialUnionColComm()); }
template<typename Ring>
int AbstractDistMatrix<Ring>::PartialUnionRowRank() const
{ return mpi::Rank(PartialUnionRowComm()); }

template<typename Ring>
int AbstractDistMatrix<Ring>::DistRank() const
{ return mpi::Rank(DistComm()); }
template<typename Ring>
int AbstractDistMatrix<Ring>::CrossRank() const
{ return mpi::Rank(CrossComm()); }
template<typename Ring>
int AbstractDistMatrix<Ring>::RedundantRank() const
{ return mpi::Rank(RedundantComm()); }

template<typename Ring>
int AbstractDistMatrix<Ring>::Root() const { return root_; }

// Single-entry manipulation
// =========================

// Global entry manipulation
// -------------------------

template<typename Ring>
Ring
AbstractDistMatrix<Ring>::Get( Int i, Int j ) const
{
    DEBUG_ONLY(
      CSE cse("ADM::Get");
      if( !grid_->InGrid() )
          LogicError("Get should only be called in-grid");
    )
    Ring value;
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

template<typename Ring>
Base<Ring>
AbstractDistMatrix<Ring>::GetRealPart( Int i, Int j ) const
{
    DEBUG_ONLY(
      CSE cse("ADM::GetRealPart");
      if( !grid_->InGrid() )
          LogicError("Get should only be called in-grid");
    )
    Base<Ring> value;
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

template<typename Ring>
Base<Ring>
AbstractDistMatrix<Ring>::GetImagPart( Int i, Int j ) const
{
    DEBUG_ONLY(
      CSE cse("ADM::GetImagPart");
      if( !grid_->InGrid() )
          LogicError("Get should only be called in-grid");
    )
    Base<Ring> value;
    if( IsComplex<Ring>::val )
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::Set( Int i, Int j, Ring value )
{
    DEBUG_ONLY(CSE cse("ADM::Set"))
    if( IsLocal(i,j) )
        SetLocal( LocalRow(i), LocalCol(j), value );
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::Set( const Entry<Ring>& entry )
{ Set( entry.i, entry.j, entry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetRealPart( Int i, Int j, Base<Ring> value )
{
    DEBUG_ONLY(CSE cse("ADM::SetRealPart"))
    if( IsLocal(i,j) )
        SetLocalRealPart( LocalRow(i), LocalCol(j), value );
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetRealPart( const Entry<Base<Ring>>& entry )
{ SetRealPart( entry.i, entry.j, entry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetImagPart( Int i, Int j, Base<Ring> value )
{
    DEBUG_ONLY(CSE cse("ADM::SetImagPart"))
    if( IsLocal(i,j) )
        SetLocalImagPart( LocalRow(i), LocalCol(j), value );
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetImagPart( const Entry<Base<Ring>>& entry )
{ SetImagPart( entry.i, entry.j, entry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::Update( Int i, Int j, Ring value )
{
    DEBUG_ONLY(CSE cse("ADM::Update"))
    if( IsLocal(i,j) )
        UpdateLocal( LocalRow(i), LocalCol(j), value );
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::Update( const Entry<Ring>& entry )
{ Update( entry.i, entry.j, entry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::UpdateRealPart( Int i, Int j, Base<Ring> value )
{
    DEBUG_ONLY(CSE cse("ADM::UpdateRealPart"))
    if( IsLocal(i,j) )
        UpdateLocalRealPart( LocalRow(i), LocalCol(j), value );
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::UpdateRealPart( const Entry<Base<Ring>>& entry )
{ UpdateRealPart( entry.i, entry.j, entry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::UpdateImagPart( Int i, Int j, Base<Ring> value )
{
    DEBUG_ONLY(CSE cse("ADM::UpdateImagPart"))
    if( IsLocal(i,j) )
        UpdateLocalImagPart( LocalRow(i), LocalCol(j), value );
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::UpdateImagPart( const Entry<Base<Ring>>& entry )
{ UpdateImagPart( entry.i, entry.j, entry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::MakeReal( Int i, Int j )
{
    DEBUG_ONLY(CSE cse("ADM::MakeReal"))
    if( IsLocal(i,j) )
        MakeLocalReal( LocalRow(i), LocalCol(j) );
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::Conjugate( Int i, Int j )
{
    DEBUG_ONLY(CSE cse("ADM::Conjugate"))
    if( IsLocal(i,j) )
        ConjugateLocal( LocalRow(i), LocalCol(j) );
}

// Batch remote updates
// --------------------
template<typename Ring>
void AbstractDistMatrix<Ring>::Reserve( Int numRemoteUpdates )
{ 
    DEBUG_ONLY(CSE cse("AbstractDistMatrix::Reserve"))
    remoteUpdates_.reserve( numRemoteUpdates ); 
}

template<typename Ring>
void AbstractDistMatrix<Ring>::QueueUpdate( const Entry<Ring>& entry )
{
    DEBUG_ONLY(CSE cse("AbstractDistMatrix::QueueUpdate"))
    if( IsLocal(entry.i,entry.j) )
        Update( entry );
    else
        remoteUpdates_.push_back( entry );
}

template<typename Ring>
void AbstractDistMatrix<Ring>::QueueUpdate( Int i, Int j, Ring value )
{ QueueUpdate( Entry<Ring>{i,j,value} ); }

template<typename Ring>
void AbstractDistMatrix<Ring>::ProcessQueues()
{
    DEBUG_ONLY(CSE cse("AbstractDistMatrix::ProcessQueues"))
    const auto& g = Grid();
    mpi::Comm comm = g.ViewingComm();
    const int commSize = mpi::Size( comm );

    // Compute the metadata
    // ====================
    vector<int> sendCounts(commSize,0);
    for( const auto& entry : remoteUpdates_ )
    {
        const int owner = 
          g.VCToViewing( 
            g.CoordsToVC(ColDist(),RowDist(),Owner(entry.i,entry.j),Root())
          );
        ++sendCounts[owner];
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    vector<Entry<Ring>> sendBuf(totalSend);
    auto offs = sendOffs;
    for( const auto& entry : remoteUpdates_ )
    {
        const int owner = 
          g.VCToViewing( 
            g.CoordsToVC(ColDist(),RowDist(),Owner(entry.i,entry.j),Root())
          );
        sendBuf[offs[owner]++] = entry;
    }

    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    Int recvBufSize = recvBuf.size();
    mpi::Broadcast( recvBufSize, 0, RedundantComm() );
    recvBuf.resize( recvBufSize );
    mpi::Broadcast( recvBuf.data(), recvBufSize, 0, RedundantComm() );
    for( const auto& entry : recvBuf )
        Update( entry );
}

// Local entry manipulation
// ------------------------

template<typename Ring>
Ring AbstractDistMatrix<Ring>::GetLocal( Int iLoc, Int jLoc ) const
{ return matrix_.Get(iLoc,jLoc); }

template<typename Ring>
Base<Ring> AbstractDistMatrix<Ring>::GetLocalRealPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetRealPart(iLoc,jLoc); }

template<typename Ring>
Base<Ring> AbstractDistMatrix<Ring>::GetLocalImagPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetImagPart(iLoc,jLoc); }

template<typename Ring>
void AbstractDistMatrix<Ring>::SetLocal( Int iLoc, Int jLoc, Ring alpha )
{ matrix_.Set(iLoc,jLoc,alpha); }

template<typename Ring>
void AbstractDistMatrix<Ring>::SetLocal( const Entry<Ring>& localEntry )
{ SetLocal( localEntry.i, localEntry.j, localEntry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetLocalRealPart( Int iLoc, Int jLoc, Base<Ring> alpha )
{ matrix_.SetRealPart(iLoc,jLoc,alpha); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetLocalRealPart( const Entry<Base<Ring>>& localEntry )
{ SetLocalRealPart( localEntry.i, localEntry.j, localEntry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetLocalImagPart( Int iLoc, Int jLoc, Base<Ring> alpha )
{ matrix_.SetImagPart(iLoc,jLoc,alpha); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetLocalImagPart( const Entry<Base<Ring>>& localEntry )
{ SetLocalImagPart( localEntry.i, localEntry.j, localEntry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::UpdateLocal( Int iLoc, Int jLoc, Ring alpha )
{ matrix_.Update(iLoc,jLoc,alpha); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::UpdateLocal( const Entry<Ring>& localEntry )
{ UpdateLocal( localEntry.i, localEntry.j, localEntry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::UpdateLocalRealPart
( Int iLoc, Int jLoc, Base<Ring> alpha )
{ matrix_.UpdateRealPart(iLoc,jLoc,alpha); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::UpdateLocalRealPart( const Entry<Base<Ring>>& localEntry )
{ UpdateLocalRealPart( localEntry.i, localEntry.j, localEntry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::UpdateLocalImagPart
( Int iLoc, Int jLoc, Base<Ring> alpha )
{ matrix_.UpdateImagPart(iLoc,jLoc,alpha); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::UpdateLocalImagPart( const Entry<Base<Ring>>& localEntry )
{ UpdateLocalImagPart( localEntry.i, localEntry.j, localEntry.value ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::MakeLocalReal( Int iLoc, Int jLoc )
{ matrix_.MakeReal( iLoc, jLoc ); }

template<typename Ring>
void
AbstractDistMatrix<Ring>::ConjugateLocal( Int iLoc, Int jLoc )
{ matrix_.Conjugate( iLoc, jLoc ); }

// Diagonal manipulation
// =====================
template<typename Ring>
bool AbstractDistMatrix<Ring>::DiagonalAlignedWith
( const El::DistData& d, Int offset ) const
{
    DEBUG_ONLY(CSE cse("ADM::DiagonalAlignedWith"))
    if( Grid() != *d.grid )
        return false;

    const Int diagRoot = DiagonalRoot(offset);
    if( diagRoot != d.root )
        return false;

    const int diagAlign = DiagonalAlign(offset);
    const Dist UDiag = DiagCol( ColDist(), RowDist() ); 
    const Dist VDiag = DiagRow( ColDist(), RowDist() );
    if( d.colDist == UDiag && d.rowDist == VDiag )
        return d.colAlign == diagAlign;
    else if( d.colDist == VDiag && d.rowDist == UDiag )
        return d.rowAlign == diagAlign;
    else
        return false;
}

template<typename Ring>
int AbstractDistMatrix<Ring>::DiagonalRoot( Int offset ) const
{
    DEBUG_ONLY(CSE cse("ADM::DiagonalRoot"))
    const auto& grid = Grid();

    if( ColDist() == MC && RowDist() == MR )
    {
        // Result is an [MD,* ] or [* ,MD]
        int owner;
        if( offset >= 0 )
        {
            const int procRow = ColAlign();
            const int procCol = (RowAlign()+offset) % RowStride();
            owner = procRow + ColStride()*procCol;
        }
        else
        {
            const int procRow = (ColAlign()-offset) % ColStride();
            const int procCol = RowAlign();
            owner = procRow + ColStride()*procCol;
        }
        return grid.Diag(owner);
    }
    else if( ColDist() == MR && RowDist() == MC )
    {
        // Result is an [MD,* ] or [* ,MD]
        int owner;
        if( offset >= 0 )
        {
            const int procCol = ColAlign();
            const int procRow = (RowAlign()+offset) % RowStride();
            owner = procRow + ColStride()*procCol;
        }
        else
        {
            const int procCol = (ColAlign()-offset) % ColStride();
            const int procRow = RowAlign();
            owner = procRow + ColStride()*procCol;
        }
        return grid.Diag(owner);
    }
    else
        return Root();
}

template<typename Ring>
int AbstractDistMatrix<Ring>::DiagonalAlign( Int offset ) const
{
    DEBUG_ONLY(CSE cse("ADM::DiagonalAlign"))
    const auto& grid = Grid();

    if( ColDist() == MC && RowDist() == MR )
    {
        // Result is an [MD,* ] or [* ,MD]
        int owner;
        if( offset >= 0 )
        {
            const int procRow = ColAlign();
            const int procCol = (RowAlign()+offset) % RowStride();
            owner = procRow + ColStride()*procCol;
        }
        else
        {
            const int procRow = (ColAlign()-offset) % ColStride();
            const int procCol = RowAlign();
            owner = procRow + ColStride()*procCol;
        }
        return grid.DiagRank(owner);
    }
    else if( ColDist() == MR && RowDist() == MC )
    {
        // Result is an [MD,* ] or [* ,MD]
        int owner;
        if( offset >= 0 )
        {
            const int procCol = ColAlign();
            const int procRow = (RowAlign()+offset) % RowStride();
            owner = procRow + ColStride()*procCol;
        }
        else
        {
            const int procCol = (ColAlign()-offset) % ColStride();
            const int procRow = RowAlign();
            owner = procRow + ColStride()*procCol;
        }
        return grid.DiagRank(owner);
    }
    else if( ColDist() == STAR )
    {
        // Result is a [V,* ] or [* ,V]
        if( offset >= 0 )
            return (RowAlign()+offset) % RowStride();
        else
            return RowAlign();
    }
    else
    {
        // Result is [U,V] or [V,U], where V is either STAR or CIRC
        if( offset >= 0 )
            return ColAlign();
        else
            return (ColAlign()-offset) % ColStride();
    }
}

// Assertions
// ==========

template<typename Ring>
void
AbstractDistMatrix<Ring>::ComplainIfReal() const
{ 
    if( !IsComplex<Ring>::val )
        LogicError("Called complex-only routine with real data");
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::AssertNotLocked() const
{
    if( Locked() )
        LogicError("Assertion that matrix not be a locked view failed");
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::AssertNotStoringData() const
{
    if( matrix_.MemorySize() > 0 )
        LogicError("Assertion that matrix not be storing data failed");
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::AssertValidEntry( Int i, Int j ) const
{
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
        LogicError
        ("Entry (",i,",",j,") is out of bounds of ",Height(),
         " x ",Width()," matrix");
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::AssertValidSubmatrix
( Int i, Int j, Int height, Int width ) const
{
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
    if( i < 0 || j < 0 )
        LogicError("Indices of submatrix were negative");
    if( height < 0 || width < 0 )
        LogicError("Dimensions of submatrix were negative");
    if( (i+height) > Height() || (j+width) > Width() )
        LogicError
        ("Submatrix is out of bounds: accessing up to (",i+height-1,
         ",",j+width-1,") of ",Height()," x ",Width()," matrix");
}

template<typename Ring> 
void
AbstractDistMatrix<Ring>::AssertSameSize( Int height, Int width ) const
{
    if( Height() != height || Width() != width )
        LogicError("Assertion that matrices be the same size failed");
}

// Private section
// ###############

// Exchange metadata with another matrix
// =====================================

template<typename Ring>
void 
AbstractDistMatrix<Ring>::ShallowSwap( AbstractDistMatrix<Ring>& A )
{
    matrix_.ShallowSwap( A.matrix_ );
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetShifts()
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

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetColShift()
{
    if( Participating() )
        colShift_ = Shift(ColRank(),colAlign_,ColStride());
    else
        colShift_ = 0;
}

template<typename Ring>
void
AbstractDistMatrix<Ring>::SetRowShift()
{
    if( Participating() )
        rowShift_ = Shift(RowRank(),rowAlign_,RowStride());
    else
        rowShift_ = 0;
}

// Outside of class
// ----------------

template<typename Ring> 
void
AssertConforming1x2
( const AbstractDistMatrix<Ring>& AL, const AbstractDistMatrix<Ring>& AR )
{
    if( AL.Height() != AR.Height() )    
        LogicError
        ("1x2 not conformant:\n",
         DimsString(AL,"Left"),"\n",DimsString(AR,"Right"));
    if( AL.ColAlign() != AR.ColAlign() )
        LogicError("1x2 is misaligned");
}

template<typename Ring> 
void
AssertConforming2x1
( const AbstractDistMatrix<Ring>& AT, const AbstractDistMatrix<Ring>& AB )
{
    if( AT.Width() != AB.Width() )
        LogicError
        ("2x1 is not conformant:\n",
         DimsString(AT,"Top"),"\n",DimsString(AB,"Bottom"));
    if( AT.RowAlign() != AB.RowAlign() )
        LogicError("2x1 is not aligned");
}

template<typename Ring> 
void
AssertConforming2x2
( const AbstractDistMatrix<Ring>& ATL, const AbstractDistMatrix<Ring>& ATR, 
  const AbstractDistMatrix<Ring>& ABL, const AbstractDistMatrix<Ring>& ABR ) 
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

#ifndef EL_RELEASE
 #define PROTO(Ring) \
  template class AbstractDistMatrix<Ring>;\
  template void AssertConforming1x2\
  ( const AbstractDistMatrix<Ring>& AL, const AbstractDistMatrix<Ring>& AR );\
  template void AssertConforming2x1\
  ( const AbstractDistMatrix<Ring>& AT, const AbstractDistMatrix<Ring>& AB );\
  template void AssertConforming2x2\
  ( const AbstractDistMatrix<Ring>& ATL, const AbstractDistMatrix<Ring>& ATR,\
    const AbstractDistMatrix<Ring>& ABL, const AbstractDistMatrix<Ring>& ABR );
#else
 #define PROTO(Ring) template class AbstractDistMatrix<Ring>;
#endif

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
