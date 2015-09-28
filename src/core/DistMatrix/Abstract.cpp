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
AbstractDistMatrix<T>::AbstractDistMatrix( const El::Grid& grid, int root )
: root_(root), grid_(&grid)
{ }

template<typename T>
AbstractDistMatrix<T>::AbstractDistMatrix( AbstractDistMatrix<T>&& A ) 
EL_NO_EXCEPT
: viewType_(A.viewType_),
  height_(A.height_),
  width_(A.width_), 
  colConstrained_(A.colConstrained_),
  rowConstrained_(A.rowConstrained_),
  rootConstrained_(A.rootConstrained_),
  colAlign_(A.colAlign_),
  rowAlign_(A.rowAlign_),
  colShift_(A.colShift_),
  rowShift_(A.rowShift_), 
  root_(A.root_),
  grid_(A.grid_)
{ matrix_.ShallowSwap( A.matrix_ ); }

template<typename T>
AbstractDistMatrix<T>::~AbstractDistMatrix() { }

// Assignment and reconfiguration
// ==============================
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
AbstractDistMatrix<T>::SetGrid( const El::Grid& grid )
{
    if( grid_ != &grid )
    {
        grid_ = &grid; 
        Empty();
    }
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
AbstractDistMatrix<T>::MakeSizeConsistent( bool includingViewers )
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

template<typename T>
void
AbstractDistMatrix<T>::SetRoot( int root, bool constrain )
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

// Operator overloading
// ====================

// Move assignment
// ---------------
template<typename T>
AbstractDistMatrix<T>& 
AbstractDistMatrix<T>::operator=( AbstractDistMatrix<T>&& A )
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
template<typename T>
const AbstractDistMatrix<T>&
AbstractDistMatrix<T>::operator*=( T alpha )
{
    DEBUG_ONLY(CSE cse("ADM::operator*=(T)"))    
    Scale( alpha, *this );
    return *this;
}

// Basic queries
// =============

// Global matrix information
// -------------------------

template<typename T>
Int AbstractDistMatrix<T>::Height() const EL_NO_EXCEPT { return height_; }
template<typename T>
Int AbstractDistMatrix<T>::Width() const EL_NO_EXCEPT { return width_; }

template<typename T>
Int AbstractDistMatrix<T>::DiagonalLength( Int offset ) const EL_NO_EXCEPT
{ return El::DiagonalLength(height_,width_,offset); }

template<typename T>
bool AbstractDistMatrix<T>::Viewing() const EL_NO_EXCEPT 
{ return IsViewing( viewType_ ); }
template<typename T>
bool AbstractDistMatrix<T>::Locked() const EL_NO_EXCEPT 
{ return IsLocked( viewType_ ); }

// Local matrix information
// ------------------------

template<typename T>
Int AbstractDistMatrix<T>::LocalHeight() const EL_NO_EXCEPT
{ return matrix_.Height(); }
template<typename T>
Int AbstractDistMatrix<T>::LocalWidth() const EL_NO_EXCEPT
{ return matrix_.Width(); }
template<typename T>
Int AbstractDistMatrix<T>::LDim() const EL_NO_EXCEPT
{ return matrix_.LDim(); }

template<typename T>
El::Matrix<T>& 
AbstractDistMatrix<T>::Matrix() EL_NO_EXCEPT { return matrix_; }
template<typename T>
const El::Matrix<T>& 
AbstractDistMatrix<T>::LockedMatrix() const EL_NO_EXCEPT { return matrix_; }

template<typename T>
size_t
AbstractDistMatrix<T>::AllocatedMemory() const EL_NO_EXCEPT 
{ return matrix_.MemorySize(); }

template<typename T>
T*
AbstractDistMatrix<T>::Buffer() EL_NO_RELEASE_EXCEPT
{ return matrix_.Buffer(); }

template<typename T>
T*
AbstractDistMatrix<T>::Buffer( Int iLoc, Int jLoc ) EL_NO_RELEASE_EXCEPT
{ return matrix_.Buffer(iLoc,jLoc); }

template<typename T>
const T*
AbstractDistMatrix<T>::LockedBuffer() const EL_NO_EXCEPT
{ return matrix_.LockedBuffer(); }

template<typename T>
const T*
AbstractDistMatrix<T>::LockedBuffer( Int iLoc, Int jLoc ) const EL_NO_EXCEPT
{ return matrix_.LockedBuffer(iLoc,jLoc); }

// Distribution information
// ------------------------

template<typename T>
const El::Grid& AbstractDistMatrix<T>::Grid() const EL_NO_EXCEPT
{ return *grid_; }

template<typename T>
int AbstractDistMatrix<T>::ColAlign() const EL_NO_EXCEPT { return colAlign_; }
template<typename T>
int AbstractDistMatrix<T>::RowAlign() const EL_NO_EXCEPT { return rowAlign_; }

template<typename T>
int AbstractDistMatrix<T>::ColShift() const EL_NO_EXCEPT { return colShift_; }
template<typename T>
int AbstractDistMatrix<T>::RowShift() const EL_NO_EXCEPT { return rowShift_; }

template<typename T>
bool AbstractDistMatrix<T>::ColConstrained() const EL_NO_EXCEPT
{ return colConstrained_; }
template<typename T>
bool AbstractDistMatrix<T>::RowConstrained() const EL_NO_EXCEPT
{ return rowConstrained_; }
template<typename T>
bool AbstractDistMatrix<T>::RootConstrained() const EL_NO_EXCEPT
{ return rootConstrained_; }

template<typename T>
bool AbstractDistMatrix<T>::Participating() const EL_NO_RELEASE_EXCEPT
{ return grid_->InGrid() && (CrossRank()==root_); }

template<typename T>
int AbstractDistMatrix<T>::Owner( Int i, Int j ) const EL_NO_EXCEPT
{ return RowOwner(i)+ColOwner(j)*ColStride(); }

template<typename T>
Int AbstractDistMatrix<T>::LocalRow( Int i ) const EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_ONLY(
      CSE cse("ADM::LocalRow");
      if( !IsLocalRow(i) )
          LogicError("Requested local index of non-local row");
    )
    return LocalRowOffset(i);
}

template<typename T>
Int AbstractDistMatrix<T>::LocalCol( Int j ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(
      CSE cse("ADM::LocalCol");
      if( !IsLocalCol(j) )
          LogicError("Requested local index of non-local column");
    )
    return LocalColOffset(j);
}

template<typename T>
bool AbstractDistMatrix<T>::IsLocalRow( Int i ) const EL_NO_RELEASE_EXCEPT
{ return Participating() && RowOwner(i) == ColRank(); }
template<typename T>
bool AbstractDistMatrix<T>::IsLocalCol( Int j ) const EL_NO_RELEASE_EXCEPT
{ return Participating() && ColOwner(j) == RowRank(); }
template<typename T>
bool AbstractDistMatrix<T>::IsLocal( Int i, Int j ) const EL_NO_RELEASE_EXCEPT
{ return IsLocalRow(i) && IsLocalCol(j); }

template<typename T>
mpi::Comm AbstractDistMatrix<T>::PartialColComm() const EL_NO_EXCEPT
{ return ColComm(); }
template<typename T>
mpi::Comm AbstractDistMatrix<T>::PartialRowComm() const EL_NO_EXCEPT
{ return RowComm(); }

template<typename T>
mpi::Comm AbstractDistMatrix<T>::PartialUnionColComm() const EL_NO_EXCEPT
{ return mpi::COMM_SELF; }
template<typename T>
mpi::Comm AbstractDistMatrix<T>::PartialUnionRowComm() const EL_NO_EXCEPT
{ return mpi::COMM_SELF; }

template<typename T>
int AbstractDistMatrix<T>::PartialColStride() const EL_NO_EXCEPT
{ return ColStride(); }
template<typename T>
int AbstractDistMatrix<T>::PartialRowStride() const EL_NO_EXCEPT
{ return RowStride(); }

template<typename T>
int AbstractDistMatrix<T>::PartialUnionColStride() const EL_NO_EXCEPT
{ return 1; }
template<typename T>
int AbstractDistMatrix<T>::PartialUnionRowStride() const EL_NO_EXCEPT
{ return 1; }

template<typename T>
int AbstractDistMatrix<T>::ColRank() const EL_NO_RELEASE_EXCEPT
{ return mpi::Rank(ColComm()); }
template<typename T>
int AbstractDistMatrix<T>::RowRank() const EL_NO_RELEASE_EXCEPT
{ return mpi::Rank(RowComm()); }

template<typename T>
int AbstractDistMatrix<T>::PartialColRank() const EL_NO_RELEASE_EXCEPT
{ return mpi::Rank(PartialColComm()); }
template<typename T>
int AbstractDistMatrix<T>::PartialRowRank() const EL_NO_RELEASE_EXCEPT
{ return mpi::Rank(PartialRowComm()); }

template<typename T>
int AbstractDistMatrix<T>::PartialUnionColRank() const EL_NO_RELEASE_EXCEPT
{ return mpi::Rank(PartialUnionColComm()); }
template<typename T>
int AbstractDistMatrix<T>::PartialUnionRowRank() const EL_NO_RELEASE_EXCEPT
{ return mpi::Rank(PartialUnionRowComm()); }

template<typename T>
int AbstractDistMatrix<T>::DistRank() const EL_NO_RELEASE_EXCEPT
{ return mpi::Rank(DistComm()); }
template<typename T>
int AbstractDistMatrix<T>::CrossRank() const EL_NO_RELEASE_EXCEPT
{ return mpi::Rank(CrossComm()); }
template<typename T>
int AbstractDistMatrix<T>::RedundantRank() const EL_NO_RELEASE_EXCEPT
{ return mpi::Rank(RedundantComm()); }

template<typename T>
int AbstractDistMatrix<T>::Root() const EL_NO_EXCEPT { return root_; }

// Single-entry manipulation
// =========================

// Global entry manipulation
// -------------------------

template<typename T>
T
AbstractDistMatrix<T>::Get( Int i, Int j ) const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(
      CSE cse("ADM::Get");
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
AbstractDistMatrix<T>::GetRealPart( Int i, Int j ) const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(
      CSE cse("ADM::GetRealPart");
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
AbstractDistMatrix<T>::GetImagPart( Int i, Int j ) const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(
      CSE cse("ADM::GetImagPart");
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
void
AbstractDistMatrix<T>::Set( Int i, Int j, T value )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("ADM::Set"))
    if( IsLocal(i,j) )
        SetLocal( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::Set( const Entry<T>& entry )
EL_NO_RELEASE_EXCEPT
{ Set( entry.i, entry.j, entry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::SetRealPart( Int i, Int j, Base<T> value )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("ADM::SetRealPart"))
    if( IsLocal(i,j) )
        SetLocalRealPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::SetRealPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ SetRealPart( entry.i, entry.j, entry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::SetImagPart( Int i, Int j, Base<T> value )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("ADM::SetImagPart"))
    if( IsLocal(i,j) )
        SetLocalImagPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::SetImagPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ SetImagPart( entry.i, entry.j, entry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::Update( Int i, Int j, T value )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("ADM::Update"))
    if( IsLocal(i,j) )
        UpdateLocal( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::Update( const Entry<T>& entry )
EL_NO_RELEASE_EXCEPT
{ Update( entry.i, entry.j, entry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateRealPart( Int i, Int j, Base<T> value )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("ADM::UpdateRealPart"))
    if( IsLocal(i,j) )
        UpdateLocalRealPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::UpdateRealPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ UpdateRealPart( entry.i, entry.j, entry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateImagPart( Int i, Int j, Base<T> value )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("ADM::UpdateImagPart"))
    if( IsLocal(i,j) )
        UpdateLocalImagPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::UpdateImagPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ UpdateImagPart( entry.i, entry.j, entry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::MakeReal( Int i, Int j )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("ADM::MakeReal"))
    if( IsLocal(i,j) )
        MakeLocalReal( LocalRow(i), LocalCol(j) );
}

template<typename T>
void
AbstractDistMatrix<T>::Conjugate( Int i, Int j )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("ADM::Conjugate"))
    if( IsLocal(i,j) )
        ConjugateLocal( LocalRow(i), LocalCol(j) );
}

// Batch remote updates
// --------------------
template<typename T>
void AbstractDistMatrix<T>::Reserve( Int numRemoteUpdates )
{ 
    DEBUG_ONLY(CSE cse("AbstractDistMatrix::Reserve"))
    remoteUpdates_.reserve( numRemoteUpdates ); 
}

template<typename T>
void AbstractDistMatrix<T>::QueueUpdate( const Entry<T>& entry )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("AbstractDistMatrix::QueueUpdate"))
    if( IsLocal(entry.i,entry.j) )
        Update( entry );
    else
        remoteUpdates_.push_back( entry );
}

template<typename T>
void AbstractDistMatrix<T>::QueueUpdate( Int i, Int j, T value )
EL_NO_RELEASE_EXCEPT
{ QueueUpdate( Entry<T>{i,j,value} ); }

template<typename T>
void AbstractDistMatrix<T>::ProcessQueues()
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
    vector<Entry<T>> sendBuf(totalSend);
    auto offs = sendOffs;
    for( const auto& entry : remoteUpdates_ )
    {
        const int owner = 
          g.VCToViewing( 
            g.CoordsToVC(ColDist(),RowDist(),Owner(entry.i,entry.j),Root())
          );
        sendBuf[offs[owner]++] = entry;
    }
    SwapClear( remoteUpdates_ );

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

template<typename T>
T AbstractDistMatrix<T>::GetLocal( Int iLoc, Int jLoc ) const
EL_NO_RELEASE_EXCEPT
{ return matrix_.Get(iLoc,jLoc); }

template<typename T>
Base<T> AbstractDistMatrix<T>::GetLocalRealPart( Int iLoc, Int jLoc ) const
EL_NO_RELEASE_EXCEPT
{ return matrix_.GetRealPart(iLoc,jLoc); }

template<typename T>
Base<T> AbstractDistMatrix<T>::GetLocalImagPart( Int iLoc, Int jLoc ) const
EL_NO_RELEASE_EXCEPT
{ return matrix_.GetImagPart(iLoc,jLoc); }

template<typename T>
void AbstractDistMatrix<T>::SetLocal( Int iLoc, Int jLoc, T alpha )
EL_NO_RELEASE_EXCEPT
{ matrix_.Set(iLoc,jLoc,alpha); }

template<typename T>
void AbstractDistMatrix<T>::SetLocal( const Entry<T>& localEntry )
EL_NO_RELEASE_EXCEPT
{ SetLocal( localEntry.i, localEntry.j, localEntry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocalRealPart( Int iLoc, Int jLoc, Base<T> alpha )
EL_NO_RELEASE_EXCEPT
{ matrix_.SetRealPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocalRealPart( const Entry<Base<T>>& localEntry )
EL_NO_RELEASE_EXCEPT
{ SetLocalRealPart( localEntry.i, localEntry.j, localEntry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocalImagPart( Int iLoc, Int jLoc, Base<T> alpha )
EL_NO_RELEASE_EXCEPT
{ matrix_.SetImagPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocalImagPart( const Entry<Base<T>>& localEntry )
EL_NO_RELEASE_EXCEPT
{ SetLocalImagPart( localEntry.i, localEntry.j, localEntry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocal( Int iLoc, Int jLoc, T alpha )
EL_NO_RELEASE_EXCEPT
{ matrix_.Update(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocal( const Entry<T>& localEntry )
EL_NO_RELEASE_EXCEPT
{ UpdateLocal( localEntry.i, localEntry.j, localEntry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocalRealPart
( Int iLoc, Int jLoc, Base<T> alpha )
EL_NO_RELEASE_EXCEPT
{ matrix_.UpdateRealPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocalRealPart( const Entry<Base<T>>& localEntry )
EL_NO_RELEASE_EXCEPT
{ UpdateLocalRealPart( localEntry.i, localEntry.j, localEntry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocalImagPart
( Int iLoc, Int jLoc, Base<T> alpha )
EL_NO_RELEASE_EXCEPT
{ matrix_.UpdateImagPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocalImagPart( const Entry<Base<T>>& localEntry )
EL_NO_RELEASE_EXCEPT
{ UpdateLocalImagPart( localEntry.i, localEntry.j, localEntry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::MakeLocalReal( Int iLoc, Int jLoc )
EL_NO_RELEASE_EXCEPT
{ matrix_.MakeReal( iLoc, jLoc ); }

template<typename T>
void
AbstractDistMatrix<T>::ConjugateLocal( Int iLoc, Int jLoc )
EL_NO_RELEASE_EXCEPT
{ matrix_.Conjugate( iLoc, jLoc ); }

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
    if( i == END ) i = height_ - 1;
    if( j == END ) j = width_ - 1;
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

// Instantiations for {Int,Real,Complex<Real>} for each Real in {float,double}
// ###########################################################################

#define PROTO(T) template class AbstractDistMatrix<T>;

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
