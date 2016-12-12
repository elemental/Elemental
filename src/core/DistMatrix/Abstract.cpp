/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1/Copy.hpp>
#include <El/blas_like/level1/Scale.hpp>

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
AbstractDistMatrix<T>::Empty( bool freeMemory )
{
    matrix_.Empty_( freeMemory );
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
    colAlign_ = 0;
    rowAlign_ = 0;
    colConstrained_ = false;
    rowConstrained_ = false;
    rootConstrained_ = false;
    SetShifts();

    SwapClear( remoteUpdates );
}

template<typename T>
void
AbstractDistMatrix<T>::EmptyData( bool freeMemory )
{
    matrix_.Empty_( freeMemory );
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
    SwapClear( remoteUpdates );
}

template<typename T>
void
AbstractDistMatrix<T>::SetGrid( const El::Grid& grid )
{
    if( grid_ != &grid )
    {
        grid_ = &grid;
        Empty(false);
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
    EL_DEBUG_CSE

    const Int msgSize = 2;
    Int message[msgSize];
    if( CrossRank() == Root() )
    {
        message[0] = height_;
        message[1] = width_;
    }

    const auto& grid = *grid_;
    if( !grid.InGrid() && !includingViewers )
        LogicError("Non-participating process called MakeSizeConsistent");
    if( grid.InGrid() )
        mpi::Broadcast( message, msgSize, Root(), CrossComm() );
    if( includingViewers )
    {
        const Int vcRoot = grid.VCToViewing(0);
        mpi::Broadcast( message, msgSize, vcRoot, grid.ViewingComm() );
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
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( root < 0 || root >= CrossSize() )
          LogicError("Invalid root");
    )
    if( root != root_ )
        EmptyData(false);
    root_ = root;
    if( constrain )
        rootConstrained_ = true;
    SetShifts();
}

// Operator overloading
// ====================

// Move assignment
// ---------------
template<typename T>
AbstractDistMatrix<T>&
AbstractDistMatrix<T>::operator=( AbstractDistMatrix<T>&& A )
{
    EL_DEBUG_CSE
    if( Viewing() || A.Viewing() )
    {
        El::Copy( A, *this );
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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( !IsLocalRow(i) )
          LogicError
          ("Row ",i," is owned by ",RowOwner(i),", not ",ColRank());
    )
    return LocalRowOffset(i);
}

template<typename T>
Int AbstractDistMatrix<T>::LocalCol( Int j ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( !IsLocalCol(j) )
          LogicError
          ("Column ",j," is owned by ",ColOwner(j),", not ",RowRank());
    )
    return LocalColOffset(j);
}

template<typename T>
Int AbstractDistMatrix<T>::LocalRow( Int i, int rowOwner ) const
EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( RowOwner(i) != rowOwner )
          LogicError
          ("Row ",i,"is owned by ",RowOwner(i)," not ",rowOwner);
    )
    return LocalRowOffset(i,rowOwner);
}

template<typename T>
Int AbstractDistMatrix<T>::LocalCol( Int j, int colOwner ) const
EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( ColOwner(j) != colOwner )
          LogicError
          ("Column ",j,"is owned by ",ColOwner(j),", not ",colOwner);
    )
    return LocalColOffset(j,colOwner);
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
int AbstractDistMatrix<T>::Root() const EL_NO_EXCEPT { return root_; }

template<typename T>
El::DistData AbstractDistMatrix<T>::DistData() const
{ return El::DistData(*this); }

// Single-entry manipulation
// =========================

// Global entry manipulation
// -------------------------

template<typename T>
T
AbstractDistMatrix<T>::Get( Int i, Int j ) const
EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( !grid_->InGrid() )
          LogicError("Get should only be called in-grid");
    )
    Base<T> value;
    if( IsComplex<T>::value )
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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
    if( IsLocal(i,j) )
        SetLocalRealPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::SetRealPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ SetRealPart( entry.i, entry.j, entry.value ); }

template<typename T>
void AbstractDistMatrix<T>::SetImagPart( Int i, Int j, Base<T> value )
EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    if( IsLocal(i,j) )
        SetLocalImagPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void AbstractDistMatrix<T>::SetImagPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ SetImagPart( entry.i, entry.j, entry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::Update( Int i, Int j, T value )
EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
    if( IsLocal(i,j) )
        UpdateLocalRealPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::UpdateRealPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ UpdateRealPart( entry.i, entry.j, entry.value ); }

template<typename T>
void AbstractDistMatrix<T>::UpdateImagPart( Int i, Int j, Base<T> value )
EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    if( IsLocal(i,j) )
        UpdateLocalImagPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void AbstractDistMatrix<T>::UpdateImagPart( const Entry<Base<T>>& entry )
EL_NO_RELEASE_EXCEPT
{ UpdateImagPart( entry.i, entry.j, entry.value ); }

template<typename T>
void
AbstractDistMatrix<T>::MakeReal( Int i, Int j )
EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    if( IsLocal(i,j) )
        MakeLocalReal( LocalRow(i), LocalCol(j) );
}

template<typename T>
void
AbstractDistMatrix<T>::Conjugate( Int i, Int j )
EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    if( IsLocal(i,j) )
        ConjugateLocal( LocalRow(i), LocalCol(j) );
}

// Batch remote updates
// --------------------
template<typename T>
void AbstractDistMatrix<T>::Reserve( Int numRemoteUpdates )
{
    EL_DEBUG_CSE
    const Int currSize = remoteUpdates.size();
    remoteUpdates.reserve( currSize+numRemoteUpdates );
}

template<typename T>
void AbstractDistMatrix<T>::QueueUpdate( const Entry<T>& entry )
EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    // NOTE: We cannot always simply locally update since it can (and has)
    //       lead to the processors in the same redundant communicator having
    //       different results after ProcessQueues()
    if( RedundantSize() == 1 && IsLocal(entry.i,entry.j) )
        UpdateLocal( LocalRow(entry.i), LocalCol(entry.j), entry.value );
    else
        remoteUpdates.push_back( entry );
}

template<typename T>
void AbstractDistMatrix<T>::QueueUpdate( Int i, Int j, T value )
EL_NO_RELEASE_EXCEPT
{ QueueUpdate( Entry<T>{i,j,value} ); }

template<typename T>
void AbstractDistMatrix<T>::ProcessQueues( bool includeViewers )
{
    EL_DEBUG_CSE
    const auto& grid = Grid();
    const Dist colDist = ColDist();
    const Dist rowDist = RowDist();
    const Int totalSend = remoteUpdates.size();

    // We will first push to redundant rank 0
    const int redundantRoot = 0;

    // Compute the metadata
    // ====================
    mpi::Comm comm;
    vector<int> sendCounts, owners(totalSend);
    if( includeViewers )
    {
        comm = grid.ViewingComm();
        const int viewingSize = mpi::Size( grid.ViewingComm() );
        sendCounts.resize(viewingSize,0);
        for( Int k=0; k<totalSend; ++k )
        {
            const Entry<T>& entry = remoteUpdates[k];
            const int distOwner = Owner(entry.i,entry.j);
            const int vcOwner =
              grid.CoordsToVC(colDist,rowDist,distOwner,redundantRoot);
            owners[k] = grid.VCToViewing(vcOwner);
            ++sendCounts[owners[k]];
        }
    }
    else
    {
        if( !Participating() )
            return;
        comm = grid.VCComm();
        const int vcSize = mpi::Size( grid.VCComm() );
        sendCounts.resize(vcSize,0);
        for( Int k=0; k<totalSend; ++k )
        {
            const Entry<T>& entry = remoteUpdates[k];
            const int distOwner = Owner(entry.i,entry.j);
            owners[k] =
              grid.CoordsToVC(colDist,rowDist,distOwner,redundantRoot);
            ++sendCounts[owners[k]];
        }
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    Scan( sendCounts, sendOffs );
    vector<Entry<T>> sendBuf(totalSend);
    auto offs = sendOffs;
    for( Int k=0; k<totalSend; ++k )
        sendBuf[offs[owners[k]]++] = remoteUpdates[k];
    SwapClear( remoteUpdates );

    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    Int recvBufSize = recvBuf.size();
    mpi::Broadcast( recvBufSize, redundantRoot, RedundantComm() );
    recvBuf.resize( recvBufSize );
    mpi::Broadcast
    ( recvBuf.data(), recvBufSize, redundantRoot, RedundantComm() );
    // TODO: Make this loop faster
    for( const auto& entry : recvBuf )
        UpdateLocal( LocalRow(entry.i), LocalCol(entry.j), entry.value );
}

template<typename T>
void AbstractDistMatrix<T>::ReservePulls( Int numPulls ) const
{
    EL_DEBUG_CSE
    remotePulls_.reserve( numPulls );
}

template<typename T>
void AbstractDistMatrix<T>::QueuePull( Int i, Int j ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    remotePulls_.push_back( ValueInt<Int>{i,j} );
}

template<typename T>
void AbstractDistMatrix<T>::ProcessPullQueue( T* pullBuf, bool includeViewers ) const
{
    EL_DEBUG_CSE
    const auto& grid = Grid();
    const Dist colDist = ColDist();
    const Dist rowDist = RowDist();
    const int root = Root();
    const Int totalRecv = remotePulls_.size();

    // Compute the metadata
    // ====================
    mpi::Comm comm;
    int commSize;
    vector<int> recvCounts, owners(totalRecv);
    if( includeViewers )
    {
        comm = grid.ViewingComm();
        commSize = mpi::Size( comm );
        recvCounts.resize(commSize,0);
        for( Int k=0; k<totalRecv; ++k )
        {
            const auto& valueInt = remotePulls_[k];
            const Int i = valueInt.value;
            const Int j = valueInt.index;
            const int distOwner = Owner(i,j);
            const int vcOwner = grid.CoordsToVC(colDist,rowDist,distOwner,root);
            const int owner = grid.VCToViewing(vcOwner);
            owners[k] = owner;
            ++recvCounts[owner];
        }
    }
    else
    {
        if( !Participating() )
            return;
        comm = grid.VCComm();
        commSize = mpi::Size( comm );
        recvCounts.resize(commSize,0);
        for( Int k=0; k<totalRecv; ++k )
        {
            const auto& valueInt = remotePulls_[k];
            const Int i = valueInt.value;
            const Int j = valueInt.index;
            const int distOwner = Owner(i,j);
            const int owner = grid.CoordsToVC(colDist,rowDist,distOwner,root);
            owners[k] = owner;
            ++recvCounts[owner];
        }
    }
    vector<int> recvOffs;
    Scan( recvCounts, recvOffs );
    vector<int> sendCounts(commSize);
    mpi::AllToAll( recvCounts.data(), 1, sendCounts.data(), 1, comm );
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );

    auto offs = recvOffs;
    vector<ValueInt<Int>> recvCoords(totalRecv);
    for( Int k=0; k<totalRecv; ++k )
        recvCoords[offs[owners[k]]++] = remotePulls_[k];
    vector<ValueInt<Int>> sendCoords(totalSend);
    mpi::AllToAll
    ( recvCoords.data(), recvCounts.data(), recvOffs.data(),
      sendCoords.data(), sendCounts.data(), sendOffs.data(), comm );

    // Pack the data
    // =============
    vector<T> sendBuf;
    FastResize( sendBuf, totalSend );
    for( Int k=0; k<totalSend; ++k )
    {
        const Int i = sendCoords[k].value;
        const Int j = sendCoords[k].index;
        sendBuf[k] = GetLocal( LocalRow(i), LocalCol(j) );
    }

    // Exchange and unpack the data
    // ============================
    vector<T> recvBuf;
    FastResize( recvBuf, totalRecv );
    mpi::AllToAll
    ( sendBuf.data(), sendCounts.data(), sendOffs.data(),
      recvBuf.data(), recvCounts.data(), recvOffs.data(), comm );
    offs = recvOffs;
    for( Int k=0; k<totalRecv; ++k )
        pullBuf[k] = recvBuf[offs[owners[k]]++];
    SwapClear( remotePulls_ );
}

template<typename T>
void AbstractDistMatrix<T>::ProcessPullQueue( vector<T>& pullVec, bool includeViewers ) const
{
    EL_DEBUG_CSE
    pullVec.resize( remotePulls_.size() );
    ProcessPullQueue( pullVec.data(), includeViewers );
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
void AbstractDistMatrix<T>::SetLocalImagPart
( Int iLoc, Int jLoc, Base<T> alpha )
EL_NO_RELEASE_EXCEPT
{ matrix_.SetImagPart(iLoc,jLoc,alpha); }

template<typename T>
void AbstractDistMatrix<T>::SetLocalImagPart
( const Entry<Base<T>>& localEntry )
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
void AbstractDistMatrix<T>::UpdateLocalImagPart
( Int iLoc, Int jLoc, Base<T> alpha )
EL_NO_RELEASE_EXCEPT
{ matrix_.UpdateImagPart(iLoc,jLoc,alpha); }

template<typename T>
void AbstractDistMatrix<T>::UpdateLocalImagPart
( const Entry<Base<T>>& localEntry )
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

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
