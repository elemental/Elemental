/*
   Copyright 2009-2011, Jack Poulson.
   All rights reserved.

   Copyright 2011-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright 2013-2014, Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   Copyright 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {

// Constructors and destructors
// ============================

template<typename T>
DistSparseMatrix<T>::DistSparseMatrix( mpi::Comm comm )
: distGraph_(comm)
{ }

template<typename T>
DistSparseMatrix<T>::DistSparseMatrix( Int height, Int width, mpi::Comm comm )
: distGraph_(height,width,comm)
{ }

template<typename T>
DistSparseMatrix<T>::DistSparseMatrix( const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::DistMultiVec"))
    distGraph_.numSources_ = -1;
    distGraph_.numTargets_ = -1;
    distGraph_.comm_ = mpi::COMM_WORLD;
    if( &A != this )
        *this = A;
    DEBUG_ONLY(
      else
          LogicError("Tried to construct DistMultiVec via itself");
    )
}

template<typename T>
DistSparseMatrix<T>::~DistSparseMatrix()
{ }

// Assignment and reconfiguration
// ==============================

// Change the matrix size
// ----------------------
template<typename T>
void DistSparseMatrix<T>::Empty( bool freeMemory )
{
    distGraph_.Empty( freeMemory );
    if( freeMemory )
        SwapClear( vals_ );
    else
        vals_.resize( 0 );
    multMeta.Clear();

    SwapClear( remoteVals_ );
}

template<typename T>
void DistSparseMatrix<T>::Resize( Int height, Int width )
{
    distGraph_.Resize( height, width );
    vals_.resize( 0 );

    SwapClear( remoteVals_ );
}

// Change the distribution
// -----------------------
template<typename T>
void DistSparseMatrix<T>::SetComm( mpi::Comm comm )
{ 
    if( Comm() == comm )
        return;
    distGraph_.SetComm( comm ); 
    vals_.resize( 0 );

    SwapClear( remoteVals_ );
}

// Assembly
// --------
template<typename T>
void DistSparseMatrix<T>::Reserve( Int numLocalEntries, Int numRemoteEntries )
{ 
    const Int currSize = vals_.size();
    const Int currRemoteSize = remoteVals_.size();

    distGraph_.Reserve( numLocalEntries, numRemoteEntries );
    vals_.reserve( currSize+numLocalEntries );
    remoteVals_.reserve( currRemoteSize+numRemoteEntries );
}

template<typename T>
void DistSparseMatrix<T>::FreezeSparsity() EL_NO_EXCEPT
{ distGraph_.frozenSparsity_ = true; }
template<typename T>
void DistSparseMatrix<T>::UnfreezeSparsity() EL_NO_EXCEPT
{ distGraph_.frozenSparsity_ = false; }
template<typename T>
bool DistSparseMatrix<T>::FrozenSparsity() const EL_NO_EXCEPT
{ return distGraph_.frozenSparsity_; }

template<typename T>
void DistSparseMatrix<T>::Update( Int row, Int col, T value )
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::Update"))
    QueueUpdate( row, col, value, true );
    ProcessLocalQueues();
}

template<typename T>
void DistSparseMatrix<T>::Update( const Entry<T>& entry )
{ Update( entry.i, entry.j, entry.value ); }

template<typename T>
void DistSparseMatrix<T>::UpdateLocal( Int localRow, Int col, T value )
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::UpdateLocal"))
    QueueLocalUpdate( localRow, col, value );
    ProcessLocalQueues();
}

template<typename T>
void DistSparseMatrix<T>::UpdateLocal( const Entry<T>& localEntry )
{ UpdateLocal( localEntry.i, localEntry.j, localEntry.value ); }

template<typename T>
void DistSparseMatrix<T>::Zero( Int row, Int col )
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::Zero"))
    QueueZero( row, col, true );
    ProcessLocalQueues();
}

template<typename T>
void DistSparseMatrix<T>::ZeroLocal( Int localRow, Int col )
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::ZeroLocal"))
    QueueLocalZero( localRow, col );
    ProcessLocalQueues();
}

template<typename T>
void DistSparseMatrix<T>::QueueUpdate( Int row, Int col, T value, bool passive )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::QueueUpdate"))
    // TODO: Use FrozenSparsity()
    if( row == END ) row = Height() - 1;
    if( col == END ) col = Width() - 1;
    if( row >= FirstLocalRow() && row < FirstLocalRow()+LocalHeight() )
    {
        QueueLocalUpdate( row-FirstLocalRow(), col, value );
    }
    else if( !passive )
    {
        distGraph_.remoteSources_.push_back( row ); 
        distGraph_.remoteTargets_.push_back( col );
        remoteVals_.push_back( value );
    }
}

template<typename T>
void DistSparseMatrix<T>::QueueUpdate( const Entry<T>& entry, bool passive )
EL_NO_RELEASE_EXCEPT
{ QueueUpdate( entry.i, entry.j, entry.value, passive ); }

template<typename T>
void DistSparseMatrix<T>::QueueLocalUpdate( Int localRow, Int col, T value )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::QueueLocalUpdate"))
    if( FrozenSparsity() )
    {
        const Int offset = distGraph_.Offset( localRow, col );
        vals_[offset] += value;
    }
    else
    {
        distGraph_.QueueLocalConnection( localRow, col );
        vals_.push_back( value );
        multMeta.ready = false;
    }
}

template<typename T>
void DistSparseMatrix<T>::QueueLocalUpdate( const Entry<T>& localEntry )
EL_NO_RELEASE_EXCEPT
{ QueueLocalUpdate( localEntry.i, localEntry.j, localEntry.value ); }

template<typename T>
void DistSparseMatrix<T>::QueueZero( Int row, Int col, bool passive )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::QueueZero"))
    if( row == END ) row = Height() - 1;
    if( col == END ) col = Width() - 1;
    if( row >= FirstLocalRow() && row < FirstLocalRow()+LocalHeight() )
        QueueLocalZero( row-FirstLocalRow(), col );
    else if( !passive )
        distGraph_.remoteRemovals_.push_back( pair<Int,Int>(row,col) );
}

template<typename T>
void DistSparseMatrix<T>::QueueLocalZero( Int localRow, Int col )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::QueueZero"))
    if( FrozenSparsity() )
    {
        const Int offset = distGraph_.Offset( localRow, col );
        vals_[offset] = 0;
    }
    else
    {
        distGraph_.QueueLocalDisconnection( localRow, col );
        multMeta.ready = false;
    }
}

template<typename T>
void DistSparseMatrix<T>::ProcessQueues()
{
    DEBUG_ONLY(
      CSE cse("DistSparseMatrix::ProcessQueues");
      if( distGraph_.sources_.size() != distGraph_.targets_.size() || 
          distGraph_.targets_.size() != vals_.size() )
          LogicError("Inconsistent sparse matrix buffer sizes");
    )

    // Send the remote updates
    // =======================
    const int commSize = distGraph_.commSize_;
    {
        // Compute the send counts
        // -----------------------
        vector<int> sendCounts(commSize,0);
        for( auto s : distGraph_.remoteSources_ )
            ++sendCounts[RowOwner(s)];
        // Pack the send data
        // ------------------
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        auto offs = sendOffs;
        vector<Entry<T>> sendBuf(totalSend);
        for( Int i=0; i<totalSend; ++i )
        {
            const int owner = RowOwner(distGraph_.remoteSources_[i]);
            sendBuf[offs[owner]++] = 
                Entry<T>
                { distGraph_.remoteSources_[i],
                  distGraph_.remoteTargets_[i], remoteVals_[i] };
        }
        SwapClear( distGraph_.remoteSources_ );
        SwapClear( distGraph_.remoteTargets_ );
        SwapClear( remoteVals_ );
        // Exchange and unpack
        // -------------------
        auto recvBuf=
          mpi::AllToAll( sendBuf, sendCounts, sendOffs, distGraph_.comm_ );
        if( !FrozenSparsity() )
            Reserve( NumLocalEntries()+recvBuf.size() );
        for( auto& entry : recvBuf )
            QueueUpdate( entry );
    }

    // Send the remote entry removals
    // ==============================
    {
        // Compute the send counts
        // -----------------------
        vector<int> sendCounts(commSize,0);
        const Int numRemoteRemovals = distGraph_.remoteRemovals_.size();
        for( Int i=0; i<numRemoteRemovals; ++i )
            ++sendCounts[RowOwner(distGraph_.remoteRemovals_[i].first)];
        // Pack the send data
        // ------------------
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        auto offs = sendOffs;
        vector<Int> sendRows(totalSend), sendCols(totalSend);
        for( Int i=0; i<totalSend; ++i )
        {
            const int owner = RowOwner(distGraph_.remoteRemovals_[i].first);
            sendRows[offs[owner]] = distGraph_.remoteRemovals_[i].first;
            sendCols[offs[owner]] = distGraph_.remoteRemovals_[i].second;
            ++offs[owner];
        }
        SwapClear( distGraph_.remoteRemovals_ );
        // Exchange and unpack
        // -------------------
        auto recvRows = 
          mpi::AllToAll(sendRows,sendCounts,sendOffs,distGraph_.comm_);
        auto recvCols = 
          mpi::AllToAll(sendCols,sendCounts,sendOffs,distGraph_.comm_);
        const Int totalRecv = recvRows.size();
        for( Int i=0; i<totalRecv; ++i )
            QueueZero( recvRows[i], recvCols[i] );
    }

    // Ensure that the kept local triplets are sorted and combined
    // ===========================================================
    ProcessLocalQueues();
}

template<typename T>
void DistSparseMatrix<T>::ProcessLocalQueues()
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::ProcessLocalQueues"))
    if( distGraph_.locallyConsistent_ )
        return;

    Int numRemoved = 0;
    const Int numLocalEntries = vals_.size();
    vector<Entry<T>> entries( numLocalEntries );
    if( distGraph_.markedForRemoval_.size() != 0 )
    {
        for( Int s=0; s<numLocalEntries; ++s )
        {
            pair<Int,Int> candidate(distGraph_.sources_[s],
                                    distGraph_.targets_[s]);
            if( distGraph_.markedForRemoval_.find(candidate) ==
                distGraph_.markedForRemoval_.end() )
            {
                entries[s-numRemoved].i = distGraph_.sources_[s];
                entries[s-numRemoved].j = distGraph_.targets_[s];
                entries[s-numRemoved].value = vals_[s];
            }
            else
            {
                ++numRemoved;
            }
        }
        SwapClear( distGraph_.markedForRemoval_ );
        entries.resize( numLocalEntries-numRemoved );
    }
    else
    {
        for( Int s=0; s<numLocalEntries; ++s )
            entries[s] = Entry<T>{distGraph_.sources_[s],
                                  distGraph_.targets_[s],vals_[s]};
    }
    std::sort( entries.begin(), entries.end(), CompareEntries );
    const Int numSorted = entries.size();

    // Combine duplicates
    // ------------------
    Int lastUnique=0;
    for( Int s=1; s<numSorted; ++s )
    {
        if( entries[s].i != entries[lastUnique].i ||
            entries[s].j != entries[lastUnique].j )
            entries[++lastUnique] = entries[s];
        else
            entries[lastUnique].value += entries[s].value;
    }
    const Int numUnique = lastUnique+1;

    entries.resize( numUnique );
    distGraph_.sources_.resize( numUnique );
    distGraph_.targets_.resize( numUnique );
    vals_.resize( numUnique );
    for( Int s=0; s<numUnique; ++s )
    {
        distGraph_.sources_[s] = entries[s].i;
        distGraph_.targets_[s] = entries[s].j;
        vals_[s] = entries[s].value;
    }
    distGraph_.ComputeSourceOffsets();
    distGraph_.locallyConsistent_ = true;
}

// Operator overloading
// ====================

// Make a copy
// -----------
template<typename T>
const DistSparseMatrix<T>& 
DistSparseMatrix<T>::operator=( const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::operator="))
    distGraph_ = A.distGraph_;
    vals_ = A.vals_;
    remoteVals_ = A.remoteVals_;
    multMeta = A.multMeta;
    return *this;
}

// Make a copy of a submatrix
// --------------------------
template<typename T>
DistSparseMatrix<T>
DistSparseMatrix<T>::operator()
( Range<Int> I, Range<Int> J ) const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::operator()"))
    DistSparseMatrix<T> ASub(this->Comm());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
DistSparseMatrix<T>
DistSparseMatrix<T>::operator()
( const vector<Int>& I, Range<Int> J ) const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::operator()"))
    DistSparseMatrix<T> ASub(this->Comm());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
DistSparseMatrix<T>
DistSparseMatrix<T>::operator()
( Range<Int> I, const vector<Int>& J ) const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::operator()"))
    DistSparseMatrix<T> ASub(this->Comm());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
DistSparseMatrix<T>
DistSparseMatrix<T>::operator()
( const vector<Int>& I, const vector<Int>& J ) const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::operator()"))
    DistSparseMatrix<T> ASub(this->Comm());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

// Rescaling
// ---------
template<typename T>
const DistSparseMatrix<T>& DistSparseMatrix<T>::operator*=( T alpha )
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::operator*=( T )"))
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename T>
const DistSparseMatrix<T>&
DistSparseMatrix<T>::operator+=( const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::operator+=( const DSM& )"))
    Axpy( T(1), A, *this );
    return *this;
}

template<typename T>
const DistSparseMatrix<T>&
DistSparseMatrix<T>::operator-=( const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::operator-=( const DSM& )"))
    Axpy( T(-1), A, *this );
    return *this;
}

// Queries
// =======

// High-level information
// ----------------------
template<typename T>
Int DistSparseMatrix<T>::Height() const EL_NO_EXCEPT
{ return distGraph_.NumSources(); }
template<typename T>
Int DistSparseMatrix<T>::Width() const EL_NO_EXCEPT
{ return distGraph_.NumTargets(); }
template<typename T>
Int DistSparseMatrix<T>::NumEntries() const EL_NO_EXCEPT
{ return distGraph_.NumEdges(); }

template<typename T>
El::DistGraph& DistSparseMatrix<T>::DistGraph() EL_NO_EXCEPT
{ return distGraph_; }
template<typename T>
const El::DistGraph& DistSparseMatrix<T>::LockedDistGraph() const EL_NO_EXCEPT
{ return distGraph_; }

template<typename T>
Int DistSparseMatrix<T>::FirstLocalRow() const EL_NO_EXCEPT
{ return distGraph_.FirstLocalSource(); }

template<typename T>
Int DistSparseMatrix<T>::LocalHeight() const EL_NO_EXCEPT
{ return distGraph_.NumLocalSources(); }

template<typename T>
Int DistSparseMatrix<T>::NumLocalEntries() const EL_NO_EXCEPT
{ return distGraph_.NumLocalEdges(); }

template<typename T>
Int DistSparseMatrix<T>::Capacity() const EL_NO_EXCEPT
{ return distGraph_.Capacity(); }

template<typename T>
bool DistSparseMatrix<T>::LocallyConsistent() const EL_NO_EXCEPT
{ return distGraph_.LocallyConsistent(); }

// Distribution information
// ------------------------
template<typename T>
mpi::Comm DistSparseMatrix<T>::Comm() const EL_NO_EXCEPT
{ return distGraph_.Comm(); }
template<typename T>
Int DistSparseMatrix<T>::Blocksize() const EL_NO_EXCEPT
{ return distGraph_.Blocksize(); }

template<typename T>
int DistSparseMatrix<T>::RowOwner( Int i ) const EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_ONLY(CSE cse("DistSparseMatrix::RowOwner"))
    return distGraph_.SourceOwner(i); 
}

template<typename T>
Int DistSparseMatrix<T>::GlobalRow( Int iLoc ) const EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_ONLY(CSE cse("DistSparseMatrix::GlobalRow"))
    return distGraph_.GlobalSource(iLoc); 
}

template<typename T>
Int DistSparseMatrix<T>::LocalRow( Int i ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::LocalRow"))
    return distGraph_.LocalSource(i);
}

// Detailed local information
// --------------------------
template<typename T>
Int DistSparseMatrix<T>::Row( Int localInd ) const
EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_ONLY(CSE cse("DistSparseMatrix::Row"))
    return distGraph_.Source( localInd );
}

template<typename T>
Int DistSparseMatrix<T>::Col( Int localInd ) const
EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_ONLY(CSE cse("DistSparseMatrix::Col"))
    return distGraph_.Target( localInd );
}

template<typename T>
Int DistSparseMatrix<T>::RowOffset( Int localRow ) const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::RowOffset"))
    return distGraph_.SourceOffset( localRow );
}

template<typename T>
Int DistSparseMatrix<T>::Offset( Int localRow, Int col ) const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::Offset"))
    return distGraph_.Offset( localRow, col );
}

template<typename T>
Int DistSparseMatrix<T>::NumConnections( Int localRow ) const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::NumConnections"))
    return distGraph_.NumConnections( localRow );
}

template<typename T>
double DistSparseMatrix<T>::Imbalance() const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::Imbalance"))
    return distGraph_.Imbalance(); 
}

template<typename T>
T DistSparseMatrix<T>::Value( Int localInd ) const
EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_ONLY(
      CSE cse("DistSparseMatrix::Value");
      if( localInd < 0 || localInd >= (Int)vals_.size() )
          LogicError("Entry number out of bounds");
      AssertLocallyConsistent();
    )
    return vals_[localInd];
}

template<typename T>
T DistSparseMatrix<T>::GetLocal( Int localRow, Int col )
const EL_NO_RELEASE_EXCEPT
{
    if( localRow == END ) localRow = LocalHeight() - 1;
    if( col == END ) col = distGraph_.numTargets_ - 1;
    Int index = Offset( localRow, col );
    if( Row(index) != GlobalRow(localRow) || Col(index) != col )
        return T(0); 
    else
        return Value( index ); 
}

template<typename T>
void DistSparseMatrix<T>::Set( Int row, Int col, T val ) EL_NO_RELEASE_EXCEPT
{
    if( row == END ) row = distGraph_.numSources_ - 1;
    if( col == END ) col = distGraph_.numTargets_ - 1;
    // NOTE: This routine currently behaves passively
    const Int firstLocalRow = FirstLocalRow();
    const Int localHeight = LocalHeight();
    if( row >= firstLocalRow && row < firstLocalRow+localHeight )
    {
        const Int localRow = row - firstLocalRow;
        Int index = Offset( localRow, col );
        if( Row(index) == row && Col(index) == col )
        {
            vals_[index] = val;
        }
        else
        {
            QueueLocalUpdate( localRow, col, val );
            ProcessLocalQueues();
        }
        ProcessQueues();
    }
}

template<typename T>
Int* DistSparseMatrix<T>::SourceBuffer() EL_NO_EXCEPT
{ return distGraph_.SourceBuffer(); }
template<typename T>
Int* DistSparseMatrix<T>::TargetBuffer() EL_NO_EXCEPT
{ return distGraph_.TargetBuffer(); }
template<typename T>
Int* DistSparseMatrix<T>::OffsetBuffer() EL_NO_EXCEPT
{ return distGraph_.OffsetBuffer(); }
template<typename T>
T* DistSparseMatrix<T>::ValueBuffer() EL_NO_EXCEPT
{ return vals_.data(); }

template<typename T>
const Int* DistSparseMatrix<T>::LockedSourceBuffer() const EL_NO_EXCEPT
{ return distGraph_.LockedSourceBuffer(); }

template<typename T>
const Int* DistSparseMatrix<T>::LockedTargetBuffer() const EL_NO_EXCEPT
{ return distGraph_.LockedTargetBuffer(); }

template<typename T>
const Int* DistSparseMatrix<T>::LockedOffsetBuffer() const EL_NO_EXCEPT
{ return distGraph_.LockedOffsetBuffer(); }

template<typename T>
const T* DistSparseMatrix<T>::LockedValueBuffer() const EL_NO_EXCEPT
{ return vals_.data(); }

template<typename T>
void DistSparseMatrix<T>::ForceNumLocalEntries( Int numLocalEntries )
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::NumLocalEntries"))
    distGraph_.ForceNumLocalEdges( numLocalEntries );
    vals_.resize( numLocalEntries );
}

template<typename T>
void DistSparseMatrix<T>::ForceConsistency( bool consistent ) EL_NO_EXCEPT
{ distGraph_.ForceConsistency(consistent); }

// Auxiliary routines
// ==================
template<typename T>
void DistSparseMatrix<T>::AssertConsistent() const
{ 
    Int locallyConsistent = ( LocallyConsistent() ? 1 : 0 );
    Int consistent = 
      mpi::AllReduce( locallyConsistent, mpi::BINARY_OR, Comm() );
    if( !consistent )
        LogicError("Distributed sparse matrix must be consistent");
}

template<typename T>
void DistSparseMatrix<T>::AssertLocallyConsistent() const
{ 
    if( !LocallyConsistent() )
        LogicError("Distributed sparse matrix must be consistent");
}

template<typename T>
DistSparseMultMeta DistSparseMatrix<T>::InitializeMultMeta() const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::InitializeMultMeta"))
    if( multMeta.ready )
        return multMeta;
    mpi::Comm comm = Comm();
    const int commSize = distGraph_.commSize_;
    auto& meta = multMeta;
 
    // Compute the set of row indices that we need from X in a normal
    // multiply or update of Y in the adjoint case
    const Int* colBuffer = LockedTargetBuffer();
    const Int numLocalEntries = NumLocalEntries();
    vector<ValueInt<Int>> uniqueCols(numLocalEntries);
    for( Int e=0; e<numLocalEntries; ++e )
        uniqueCols[e] = ValueInt<Int>{colBuffer[e],e};
    std::sort( uniqueCols.begin(), uniqueCols.end(), ValueInt<Int>::Lesser );
    meta.colOffs.resize(numLocalEntries);
    {
        Int uniqueOff=-1, lastUnique=-1;
        for( Int e=0; e<numLocalEntries; ++e )    
        {
            if( lastUnique != uniqueCols[e].value )
            {
                ++uniqueOff;
                lastUnique = uniqueCols[e].value;
                uniqueCols[uniqueOff] = uniqueCols[e];
            }
            meta.colOffs[uniqueCols[e].index] = uniqueOff;
        }
        uniqueCols.resize( uniqueOff+1 );
    }
    const Int numRecvInds = uniqueCols.size();
    meta.numRecvInds = numRecvInds;
    vector<Int> recvInds( numRecvInds );
    meta.recvSizes.clear();
    meta.recvSizes.resize( commSize, 0 );
    meta.recvOffs.resize( commSize );
    Int vecBlocksize = Width() / commSize;
    if( vecBlocksize*commSize < Width() || Width() == 0 ) 
        ++vecBlocksize;

    {
        Int off=0, lastOff=0, qPrev=0;
        for( ; off<numRecvInds; ++off )
        {
            const Int j = uniqueCols[off].value;
            const int q = j / vecBlocksize;
            while( qPrev != q )
            {
                meta.recvSizes[qPrev] = off - lastOff;
                meta.recvOffs[qPrev+1] = off;

                lastOff = off;
                ++qPrev;
            }
            recvInds[off] = j;
        }
        while( qPrev != commSize-1 )
        {
            meta.recvSizes[qPrev] = off - lastOff;
            meta.recvOffs[qPrev+1] = off;
            lastOff = off;
            ++qPrev;
        }
        meta.recvSizes[commSize-1] = off - lastOff;
    }

    // Coordinate
    meta.sendSizes.resize( commSize );
    mpi::AllToAll( meta.recvSizes.data(), 1, meta.sendSizes.data(), 1, comm );
    Int numSendInds=0;
    meta.sendOffs.resize( commSize );
    for( int q=0; q<commSize; ++q )
    {
        meta.sendOffs[q] = numSendInds;
        numSendInds += meta.sendSizes[q];
    }
    meta.sendInds.resize( numSendInds );
    mpi::AllToAll
    ( recvInds.data(),      meta.recvSizes.data(), meta.recvOffs.data(),
      meta.sendInds.data(), meta.sendSizes.data(), meta.sendOffs.data(),
      comm );

    meta.numRecvInds = numRecvInds;
    meta.ready = true;

    return meta;
}

template<typename T>
void DistSparseMatrix<T>::MappedSources
( const DistMap& reordering, vector<Int>& mappedSources ) const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::MappedSources"))
    mpi::Comm comm = Comm();
    const int commRank = mpi::Rank( comm );
    Timer timer;
    const bool time = false; 
    const Int localHeight = LocalHeight();
    if( Int(mappedSources.size()) == localHeight )
        return;

    // Get the reordered indices of our local rows of the sparse matrix
    if( time && commRank == 0 )
        timer.Start();
    mappedSources.resize( localHeight );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        mappedSources[iLoc] = GlobalRow(iLoc);
    reordering.Translate( mappedSources );
    if( time && commRank == 0 )
        Output("Source translation: ",timer.Stop()," secs");
}

template<typename T>
void DistSparseMatrix<T>::MappedTargets
( const DistMap& reordering, 
  vector<Int>& mappedTargets,
  vector<Int>& colOffs ) const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::MappedTargets"))
    if( mappedTargets.size() != 0 && colOffs.size() != 0 ) 
        return;

    mpi::Comm comm = Comm();
    const int commRank = mpi::Rank( comm ); 
    Timer timer;
    const bool time = false;

    // Compute the unique set of column indices that our process interacts with
    if( time && commRank == 0 )
        timer.Start();
    const Int* colBuffer = LockedTargetBuffer();
    const Int numLocalEntries = NumLocalEntries();
    colOffs.resize( numLocalEntries );
    vector<ValueInt<Int>> uniqueCols(numLocalEntries);
    for( Int e=0; e<numLocalEntries; ++e )
        uniqueCols[e] = ValueInt<Int>{colBuffer[e],e};
    std::sort( uniqueCols.begin(), uniqueCols.end(), ValueInt<Int>::Lesser );
    {
        Int uniqueOff=-1, lastUnique=-1;
        for( Int e=0; e<numLocalEntries; ++e )
        {
            if( lastUnique != uniqueCols[e].value )
            {
                ++uniqueOff;
                lastUnique = uniqueCols[e].value;
                uniqueCols[uniqueOff] = uniqueCols[e];
            }
            colOffs[uniqueCols[e].index] = uniqueOff;
        }
        uniqueCols.resize( uniqueOff+1 );
    }
    const Int numUniqueCols = uniqueCols.size();
    if( time && commRank == 0 )
        Output("Unique sort: ",timer.Stop()," secs");

    // Get the reordered indices of the targets of our portion of the 
    // distributed sparse matrix
    if( time && commRank == 0 )
        timer.Start();
    mappedTargets.resize( numUniqueCols );
    for( Int e=0; e<numUniqueCols; ++e )
        mappedTargets[e] = uniqueCols[e].value;
    reordering.Translate( mappedTargets );
    if( time && commRank == 0 )
        Output("Target translation: ",timer.Stop()," secs");
}

template<typename T>
bool DistSparseMatrix<T>::CompareEntries( const Entry<T>& a, const Entry<T>& b )
{ return a.i < b.i || (a.i == b.i && a.j < b.j); }

#ifdef EL_INSTANTIATE_CORE
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) EL_EXTERN template class DistSparseMatrix<T>;
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

#undef EL_EXTERN

} // namespace El
