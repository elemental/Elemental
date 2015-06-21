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
#include "El.hpp"

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
DistSparseMatrix<T>::~DistSparseMatrix()
{ }

// Assignment and reconfiguration
// ==============================

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
DistSparseMatrix<T>::operator()( Range<Int> I, Range<Int> J ) const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::operator()"))
    return GetSubmatrix( *this, I, J );
}   

// Change the matrix size
// ----------------------
template<typename T>
void DistSparseMatrix<T>::Empty( bool clearMemory )
{
    distGraph_.Empty( clearMemory );
    if( clearMemory )
        SwapClear( vals_ );
    else
        vals_.resize( 0 );
    multMeta.Clear();
}

template<typename T>
void DistSparseMatrix<T>::Resize( Int height, Int width )
{
    distGraph_.Resize( height, width );
    vals_.resize( 0 );
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
}

// Assembly
// --------
template<typename T>
void DistSparseMatrix<T>::Reserve( Int numLocalEntries, Int numRemoteEntries )
{ 
    distGraph_.Reserve( numLocalEntries, numRemoteEntries );
    vals_.reserve( numLocalEntries );
    remoteVals_.reserve( numRemoteEntries );
}

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
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::QueueUpdate"))
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
{ QueueUpdate( entry.i, entry.j, entry.value, passive ); }

template<typename T>
void DistSparseMatrix<T>::QueueLocalUpdate( Int localRow, Int col, T value )
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::QueueLocalUpdate"))
    distGraph_.QueueLocalConnection( localRow, col );
    vals_.push_back( value );
    multMeta.ready = false;
}

template<typename T>
void DistSparseMatrix<T>::QueueLocalUpdate( const Entry<T>& localEntry )
{ QueueLocalUpdate( localEntry.i, localEntry.j, localEntry.value ); }

template<typename T>
void DistSparseMatrix<T>::QueueZero( Int row, Int col, bool passive )
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
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::QueueZero"))
    distGraph_.QueueLocalDisconnection( localRow, col );
    multMeta.ready = false;
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
    const int commSize = mpi::Size( distGraph_.comm_ );
    {
        // Compute the send counts
        // -----------------------
        vector<int> sendCounts(commSize);
        for( auto s : distGraph_.remoteSources_ )
            ++sendCounts[RowOwner(s)];
        // Pack the send data
        // ------------------
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        auto offs = sendOffs;
        vector<Entry<T>> sendBuf(totalSend);
        for( Int i=0; i<distGraph_.remoteSources_.size(); ++i )
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
        Reserve( NumLocalEntries()+recvBuf.size() );
        for( auto& entry : recvBuf )
            QueueUpdate( entry );
    }

    // Send the remote entry removals
    // ==============================
    {
        // Compute the send counts
        // -----------------------
        vector<int> sendCounts(commSize);
        for( Int i=0; i<distGraph_.remoteRemovals_.size(); ++i )
            ++sendCounts[RowOwner(distGraph_.remoteRemovals_[i].first)];
        // Pack the send data
        // ------------------
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        auto offs = sendOffs;
        vector<Int> sendRows(totalSend), sendCols(totalSend);
        for( Int i=0; i<distGraph_.remoteRemovals_.size(); ++i )
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
        for( Int i=0; i<recvRows.size(); ++i )
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

    const Int numLocalEntries = vals_.size();
    Int numRemoved = 0;
    vector<Entry<T>> entries( numLocalEntries );
    for( Int s=0; s<numLocalEntries; ++s )
    {
        pair<Int,Int> candidate(distGraph_.sources_[s],distGraph_.targets_[s]);
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
    std::sort( entries.begin(), entries.end(), CompareEntries );
    // Combine duplicates
    // ------------------
    Int lastUnique=0;
    for( Int s=1; s<numLocalEntries; ++s )
    {
        if( entries[s].i != entries[lastUnique].i ||
            entries[s].j != entries[lastUnique].j )
        {
            ++lastUnique;
            entries[lastUnique] = entries[s];
        }
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
    distGraph_.ComputeEdgeOffsets();
    distGraph_.locallyConsistent_ = true;
}

// Queries
// =======

// High-level information
// ----------------------
template<typename T>
Int DistSparseMatrix<T>::Height() const { return distGraph_.NumSources(); }
template<typename T>
Int DistSparseMatrix<T>::Width() const { return distGraph_.NumTargets(); }

template<typename T>
El::DistGraph& DistSparseMatrix<T>::DistGraph() { return distGraph_; }
template<typename T>
const El::DistGraph& DistSparseMatrix<T>::LockedDistGraph() const
{ return distGraph_; }

template<typename T>
Int DistSparseMatrix<T>::FirstLocalRow() const
{ return distGraph_.FirstLocalSource(); }

template<typename T>
Int DistSparseMatrix<T>::LocalHeight() const
{ return distGraph_.NumLocalSources(); }

template<typename T>
Int DistSparseMatrix<T>::NumLocalEntries() const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::NumLocalEntries"))
    return distGraph_.NumLocalEdges();
}

template<typename T>
Int DistSparseMatrix<T>::Capacity() const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::Capacity"))
    return distGraph_.Capacity();
}

template<typename T>
bool DistSparseMatrix<T>::LocallyConsistent() const
{ return distGraph_.LocallyConsistent(); }

// Distribution information
// ------------------------
template<typename T>
mpi::Comm DistSparseMatrix<T>::Comm() const { return distGraph_.Comm(); }
template<typename T>
Int DistSparseMatrix<T>::Blocksize() const { return distGraph_.Blocksize(); }

template<typename T>
int DistSparseMatrix<T>::RowOwner( Int i ) const 
{ 
    if( i == END ) i = Height() - 1;
    return distGraph_.SourceOwner(i); 
}

template<typename T>
Int DistSparseMatrix<T>::GlobalRow( Int iLoc ) const
{ 
    DEBUG_ONLY(CSE cse("DistSparseMatrix::GlobalRow"))
    if( iLoc == END ) iLoc = LocalHeight() - 1;
    return distGraph_.GlobalSource(iLoc); 
}

template<typename T>
Int DistSparseMatrix<T>::LocalRow( Int i ) const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::LocalRow"))
    if( i == END ) i = Height() - 1;
    return distGraph_.LocalSource(i);
}

// Detailed local information
// --------------------------
template<typename T>
Int DistSparseMatrix<T>::Row( Int localInd ) const
{ 
    DEBUG_ONLY(CSE cse("DistSparseMatrix::Row"))
    return distGraph_.Source( localInd );
}

template<typename T>
Int DistSparseMatrix<T>::Col( Int localInd ) const
{ 
    DEBUG_ONLY(CSE cse("DistSparseMatrix::Col"))
    return distGraph_.Target( localInd );
}

template<typename T>
Int DistSparseMatrix<T>::EntryOffset( Int localRow ) const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::EntryOffset"))
    if( localRow == END ) localRow = LocalHeight() - 1;
    return distGraph_.EdgeOffset( localRow );
}

template<typename T>
Int DistSparseMatrix<T>::NumConnections( Int localRow ) const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::NumConnections"))
    if( localRow == END ) localRow = LocalHeight() - 1;
    return distGraph_.NumConnections( localRow );
}

template<typename T>
T DistSparseMatrix<T>::Value( Int localInd ) const
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
Int* DistSparseMatrix<T>::SourceBuffer() { return distGraph_.SourceBuffer(); }
template<typename T>
Int* DistSparseMatrix<T>::TargetBuffer() { return distGraph_.TargetBuffer(); }
template<typename T>
Int* DistSparseMatrix<T>::OffsetBuffer() { return distGraph_.OffsetBuffer(); }
template<typename T>
T* DistSparseMatrix<T>::ValueBuffer() { return vals_.data(); }

template<typename T>
const Int* DistSparseMatrix<T>::LockedSourceBuffer() const
{ return distGraph_.LockedSourceBuffer(); }

template<typename T>
const Int* DistSparseMatrix<T>::LockedTargetBuffer() const
{ return distGraph_.LockedTargetBuffer(); }

template<typename T>
const Int* DistSparseMatrix<T>::LockedOffsetBuffer() const
{ return distGraph_.LockedOffsetBuffer(); }

template<typename T>
const T* DistSparseMatrix<T>::LockedValueBuffer() const
{ return vals_.data(); }

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
    const int commSize = mpi::Size( comm );
    auto& meta = multMeta;
 
    // Compute the set of row indices that we need from X in a normal
    // multiply or update of Y in the adjoint case
    const Int numLocalEntries = NumLocalEntries();
    set<Int> indexSet;
    for( Int e=0; e<numLocalEntries; ++e )
        indexSet.insert( Col(e) );
    const Int numRecvInds = indexSet.size();
    vector<Int> recvInds( numRecvInds );
    meta.recvSizes.clear();
    meta.recvSizes.resize( commSize, 0 );
    meta.recvOffs.resize( commSize );
    const Int vecBlocksize = Width() / commSize;
    {
        Int off=0, lastOff=0, qPrev=0;
        for( auto setIt=indexSet.cbegin(); setIt!=indexSet.cend(); ++setIt )
        {
            const Int j = *setIt;
            const Int q = RowToProcess( j, vecBlocksize, commSize );
            while( qPrev != q )
            {
                meta.recvSizes[qPrev] = off - lastOff;
                meta.recvOffs[qPrev+1] = off;

                lastOff = off;
                ++qPrev;
            }
            recvInds[off++] = j;
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

    meta.colOffs.resize( numLocalEntries );
    for( Int s=0; s<numLocalEntries; ++s )
        meta.colOffs[s] = Find( recvInds, Col(s) );
    meta.numRecvInds = numRecvInds;
    meta.ready = true;

    return meta;
}

template<typename T>
bool DistSparseMatrix<T>::CompareEntries( const Entry<T>& a, const Entry<T>& b )
{ return a.i < b.i || (a.i == b.i && a.j < b.j); }

#define PROTO(T) template class DistSparseMatrix<T>;
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
