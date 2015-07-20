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
// TODO: Always duplicate the communicator and do not treat mpi::COMM_WORLD
//       and mpi::COMM_SELF as special cases.

DistGraph::DistGraph( mpi::Comm comm )
: numSources_(0), numTargets_(0)
{ 
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );
    InitializeLocalData();
}

DistGraph::DistGraph( Int numSources, mpi::Comm comm )
: numSources_(numSources), numTargets_(numSources)
{ 
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );
    InitializeLocalData();
}

DistGraph::DistGraph( Int numSources, Int numTargets, mpi::Comm comm )
: numSources_(numSources), numTargets_(numTargets)
{ 
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );
    InitializeLocalData();
}

DistGraph::DistGraph( const Graph& graph )
: numSources_(-1), numTargets_(-1)
{
    DEBUG_ONLY(CSE cse("DistGraph::DistGraph"))
    *this = graph;
}

DistGraph::DistGraph( const DistGraph& graph )
: numSources_(-1), numTargets_(-1), comm_(mpi::COMM_WORLD)
{
    DEBUG_ONLY(CSE cse("DistGraph::DistGraph"))
    if( &graph != this )
        *this = graph;
    else
        LogicError("Tried to construct DistGraph with itself");
}

DistGraph::~DistGraph()
{ 
    if( !mpi::Finalized() )
        if( comm_ != mpi::COMM_WORLD )
            mpi::Free( comm_ );
} 

// Assignment and reconfiguration
// ==============================

// Make a copy
// -----------
const DistGraph& DistGraph::operator=( const Graph& graph )
{
    DEBUG_ONLY(CSE cse("DistGraph::operator="))
    Copy( graph, *this );
    return *this;
}

const DistGraph& DistGraph::operator=( const DistGraph& graph )
{
    DEBUG_ONLY(CSE cse("DistGraph::operator="))
    Copy( graph, *this );
    return *this;
}

// Make a copy of a contiguous subgraph
// ------------------------------------
DistGraph DistGraph::operator()( Range<Int> I, Range<Int> J ) const
{
    DEBUG_ONLY(CSE cse("DistGraph::operator()"))
    return GetSubgraph( *this, I, J );
}

// Change the graph size
// ---------------------
void DistGraph::Empty( bool clearMemory )
{
    numSources_ = 0;
    numTargets_ = 0;
    firstLocalSource_ = 0;
    numLocalSources_ = 0;
    blocksize_ = 0;
    locallyConsistent_ = true;
    if( clearMemory )
    {
        SwapClear( sources_ );
        SwapClear( targets_ );
        SwapClear( localSourceOffsets_ );
    }
    else
    {
        sources_.resize( 0 );
        targets_.resize( 0 );
    }
    localSourceOffsets_.resize( 1 );
    localSourceOffsets_[0] = 0;
}

void DistGraph::Resize( Int numVertices ) 
{ Resize( numVertices, numVertices ); }

void DistGraph::Resize( Int numSources, Int numTargets )
{
    if( numSources_ == numSources && numTargets == numTargets_ )
        return;

    const int commRank = mpi::Rank( comm_ );
    const int commSize = mpi::Size( comm_ );
    numSources_ = numSources;
    numTargets_ = numTargets;
    blocksize_ = numSources/commSize;
    firstLocalSource_ = commRank*blocksize_;
    if( commRank < commSize-1 )
        numLocalSources_ = blocksize_;
    else
        numLocalSources_ = numSources - (commSize-1)*blocksize_;

    sources_.resize( 0 );
    targets_.resize( 0 );

    localSourceOffsets_.resize( numLocalSources_+1 );
    for( Int e=0; e<=numLocalSources_; ++e )
        localSourceOffsets_[e] = 0;

    locallyConsistent_ = true;
}

// Change the distribution
// -----------------------
void DistGraph::InitializeLocalData()
{
    const int commRank = mpi::Rank( comm_ );
    const int commSize = mpi::Size( comm_ );
    blocksize_ = numSources_/commSize;
    firstLocalSource_ = commRank*blocksize_;
    if( commRank < commSize-1 )
        numLocalSources_ = blocksize_;
    else
        numLocalSources_ = numSources_ - (commSize-1)*blocksize_;

    localSourceOffsets_.resize( numLocalSources_+1 );
    for( Int e=0; e<=numLocalSources_; ++e )
        localSourceOffsets_[e] = 0;

    sources_.resize( 0 );
    targets_.resize( 0 );
    locallyConsistent_ = true;
}

void DistGraph::SetComm( mpi::Comm comm )
{
    if( comm == comm_ )
        return;

    if( comm_ != mpi::COMM_WORLD )
        mpi::Free( comm_ );
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );

    InitializeLocalData();
}

// Assembly
// --------
void DistGraph::Reserve( Int numLocalEdges, Int numRemoteEdges )
{ 
    sources_.reserve( numLocalEdges );
    targets_.reserve( numLocalEdges );
    remoteSources_.reserve( numRemoteEdges );
    remoteTargets_.reserve( numRemoteEdges );
}

void DistGraph::Connect( Int source, Int target )
{
    DEBUG_ONLY(CSE cse("DistGraph::Connect"))
    QueueConnection( source, target, true );
    ProcessLocalQueues();
}

void DistGraph::ConnectLocal( Int localSource, Int target )
{
    DEBUG_ONLY(CSE cse("DistGraph::ConnectLocal"))
    QueueLocalConnection( localSource, target );
    ProcessLocalQueues();
}

void DistGraph::Disconnect( Int source, Int target )
{
    DEBUG_ONLY(CSE cse("DistGraph::Disconnect"))
    QueueDisconnection( source, target, true );
    ProcessLocalQueues();
}

void DistGraph::DisconnectLocal( Int localSource, Int target )
{
    DEBUG_ONLY(CSE cse("DistGraph::DisconnectLocal"))
    QueueLocalDisconnection( localSource, target );
    ProcessLocalQueues();
}

void DistGraph::FreezeSparsity() { frozenSparsity_ = true; }
void DistGraph::UnfreezeSparsity() { frozenSparsity_ = false; }
bool DistGraph::FrozenSparsity() const { return frozenSparsity_; }

void DistGraph::QueueConnection( Int source, Int target, bool passive )
{
    DEBUG_ONLY(CSE cse("DistGraph::QueueConnection"))
    if( source == END ) source = numSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    if( source < firstLocalSource_ || 
        source >= firstLocalSource_+numLocalSources_ )
    {
        QueueLocalConnection( source-firstLocalSource_, target );
    }
    else if( !passive )
    {
        remoteSources_.push_back( source );
        remoteTargets_.push_back( target );
    }
}

void DistGraph::QueueLocalConnection( Int localSource, Int target )
{
    DEBUG_ONLY(CSE cse("DistGraph::QueueLocalConnection"))
    DEBUG_ONLY(
      if( NumLocalEdges() == Capacity() )
          cerr << "WARNING: Pushing back without first reserving space" << endl;
    )
    if( localSource == END ) localSource = numLocalSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    if( localSource < 0 || localSource >= numLocalSources_ )
        LogicError
        ("Local source was out of bounds: ",localSource," is not in [0,",
         numLocalSources_,")");
    if( target < 0 || target >= numTargets_ )
        LogicError
        ("Target was out of bounds: ",target," is not in [0,",numTargets_,")");
    if( !FrozenSparsity() )
    {
        sources_.push_back( firstLocalSource_+localSource );
        targets_.push_back( target );
        locallyConsistent_ = false;
    }
}

void DistGraph::QueueDisconnection( Int source, Int target, bool passive )
{
    DEBUG_ONLY(CSE cse("DistGraph::QueueDisconnection"))
    // TODO: Use FrozenSparsity()
    if( source == END ) source = numSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    if( source < firstLocalSource_ || 
        source >= firstLocalSource_+numLocalSources_ )
    {
        QueueLocalDisconnection( source-firstLocalSource_, target );
    }
    else if( !passive )
    {
        remoteRemovals_.push_back( pair<Int,Int>(source,target) );
    }
}

void DistGraph::QueueLocalDisconnection( Int localSource, Int target )
{
    DEBUG_ONLY(CSE cse("DistGraph::QueueLocalDisconnection"))
    // TODO: Use FrozenSparsity()
    if( localSource == END ) localSource = numLocalSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    if( localSource < 0 || localSource >= numLocalSources_ )
        LogicError
        ("Local source was out of bounds: ",localSource," is not in [0,",
         numLocalSources_,")");
    if( target < 0 || target >= numTargets_ )
        LogicError
        ("Target was out of bounds: ",target," is not in [0,",numTargets_,")");
    if( !FrozenSparsity() )
    {
        markedForRemoval_.insert
        ( pair<Int,Int>(firstLocalSource_+localSource,target) );
        locallyConsistent_ = false;
    }
    // else throw error?
}

void DistGraph::ProcessQueues()
{
    DEBUG_ONLY(
      CSE cse("DistGraph::ProcessQueues");
      if( sources_.size() != targets_.size() )
          LogicError("Inconsistent graph buffer sizes");
    )

    // Send the remote edges
    // =====================
    int commSize = mpi::Size( comm_ );
    {
        // Compute the send counts
        // -----------------------
        vector<int> sendCounts(commSize);
        for( auto s : remoteSources_ )
            ++sendCounts[SourceOwner(s)];
        // Pack the send data
        // ------------------
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        auto offs = sendOffs;
        vector<Int> sendSources(totalSend), sendTargets(totalSend);
        for( Int i=0; i<remoteSources_.size(); ++i ) 
        {
            const int owner = SourceOwner(remoteSources_[i]);
            sendSources[offs[owner]] = remoteSources_[i];
            sendTargets[offs[owner]] = remoteTargets_[i]; 
            ++offs[owner];
        }
        SwapClear( remoteSources_ );
        SwapClear( remoteTargets_ );
        // Exchange and unpack
        // -------------------
        auto recvSources = 
          mpi::AllToAll( sendSources, sendCounts, sendOffs, comm_ );
        auto recvTargets = 
          mpi::AllToAll( sendTargets, sendCounts, sendOffs, comm_ ); 
        if( !FrozenSparsity() )
            Reserve( NumLocalEdges()+recvSources.size() );
        for( Int i=0; i<recvSources.size(); ++i )
            QueueConnection( recvSources[i], recvTargets[i] );
    }

    // Send the remote edge removals
    // =============================
    {
        // Compute the send counts
        // -----------------------
        vector<int> sendCounts(commSize);
        for( Int i=0; i<remoteRemovals_.size(); ++i )
            ++sendCounts[SourceOwner(remoteRemovals_[i].first)];
        // Pack the send data
        // ------------------
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        auto offs = sendOffs;
        vector<Int> sendSources(totalSend), sendTargets(totalSend);
        for( Int i=0; i<remoteRemovals_.size(); ++i ) 
        {
            const int owner = SourceOwner(remoteRemovals_[i].first);
            sendSources[offs[owner]] = remoteRemovals_[i].first;
            sendTargets[offs[owner]] = remoteRemovals_[i].second; 
            ++offs[owner];
        }
        SwapClear( remoteRemovals_ );
        // Exchange and unpack
        // -------------------
        auto recvSources = 
          mpi::AllToAll( sendSources, sendCounts, sendOffs, comm_ );
        auto recvTargets = 
          mpi::AllToAll( sendTargets, sendCounts, sendOffs, comm_ ); 
        for( Int i=0; i<recvSources.size(); ++i )
            QueueDisconnection( recvSources[i], recvTargets[i] );
    }

    // Ensure that the kept local edges are sorted and unique
    // ======================================================
    ProcessLocalQueues();
}

void DistGraph::ProcessLocalQueues()
{
    DEBUG_ONLY(CSE cse("DistGraph::ProcessLocalQueues"))
    if( locallyConsistent_ )
        return;

    const Int numLocalEdges = sources_.size();
    Int numRemoved = 0;
    vector<pair<Int,Int>> pairs( numLocalEdges );
    if( markedForRemoval_.size() != 0 )
    {
        for( Int e=0; e<numLocalEdges; ++e )
        {
            pair<Int,Int> candidate(sources_[e],targets_[e]);
            if( markedForRemoval_.find(candidate) == markedForRemoval_.end() )
                pairs[e-numRemoved] = candidate;
            else
                ++numRemoved;
        }
        SwapClear( markedForRemoval_ );
        pairs.resize( numLocalEdges-numRemoved );
    }
    else
    {
        for( Int e=0; e<numLocalEdges; ++e )
            pairs[e] = pair<Int,Int>{sources_[e],targets_[e]};
    }
    std::sort( pairs.begin(), pairs.end() );
    const Int numSorted = pairs.size();

    // Compress out duplicates
    Int lastUnique=0;
    for( Int e=1; e<numSorted; ++e )
        if( pairs[e] != pairs[lastUnique] )
            pairs[++lastUnique] = pairs[e];
    const Int numUnique = lastUnique+1;
    pairs.resize( numUnique );

    sources_.resize( numUnique );
    targets_.resize( numUnique );
    for( Int e=0; e<numUnique; ++e )
    {
        sources_[e] = pairs[e].first;
        targets_[e] = pairs[e].second;
    }
    ComputeSourceOffsets();
    locallyConsistent_ = true;
}

// Basic queries
// =============

// High-level information
// ----------------------
Int DistGraph::NumSources() const { return numSources_; }
Int DistGraph::NumTargets() const { return numTargets_; }
Int DistGraph::FirstLocalSource() const { return firstLocalSource_; }
Int DistGraph::NumLocalSources() const { return numLocalSources_; }

Int DistGraph::NumLocalEdges() const
{
    DEBUG_ONLY(CSE cse("DistGraph::NumLocalEdges"))
    return sources_.size();
}

Int DistGraph::Capacity() const
{
    DEBUG_ONLY(CSE cse("DistGraph::Capacity"))
    return Min(sources_.capacity(),targets_.capacity());
}

bool DistGraph::LocallyConsistent() const { return locallyConsistent_; }

// Distribution information
// ------------------------
mpi::Comm DistGraph::Comm() const { return comm_; }
Int DistGraph::Blocksize() const { return blocksize_; }

int DistGraph::SourceOwner( Int source ) const
{ 
    if( source == END ) source = numSources_ - 1;
    // TODO: Get rid of RowToProcess...
    return RowToProcess( source, Blocksize(), mpi::Size(Comm()) ); 
}

Int DistGraph::GlobalSource( Int sLoc ) const
{
    DEBUG_ONLY(CSE cse("DistGraph::GlobalSource"))
    if( sLoc == END ) sLoc = numLocalSources_ - 1;
    if( sLoc < 0 || sLoc >= NumLocalSources() )
        LogicError("Invalid local source index");
    return sLoc + FirstLocalSource();
}

Int DistGraph::LocalSource( Int s ) const
{
    DEBUG_ONLY(CSE cse("DistGraph::LocalSource"))
    if( s == END ) s = numSources_ - 1;
    if( s < 0 || s >= NumSources() )
        LogicError("Invalid global source index");
    return s - FirstLocalSource();
}

// Detailed local information
// --------------------------
Int DistGraph::Source( Int localEdge ) const
{
    DEBUG_ONLY(
      CSE cse("DistGraph::Source");
      if( localEdge < 0 || localEdge >= (Int)sources_.size() )
          LogicError("Edge number out of bounds");
    )
    return sources_[localEdge];
}

Int DistGraph::Target( Int localEdge ) const
{
    DEBUG_ONLY(
      CSE cse("DistGraph::Target");
      if( localEdge < 0 || localEdge >= (Int)targets_.size() )
          LogicError("Edge number out of bounds");
    )
    return targets_[localEdge];
}

Int DistGraph::SourceOffset( Int localSource ) const
{
    if( localSource == END ) localSource = numLocalSources_ - 1;
    DEBUG_ONLY(
      CSE cse("DistGraph::SourceOffset");
      if( localSource < 0 || localSource > numLocalSources_ )
          LogicError
          ("Out of bounds localSource: ",localSource,
           " is not in [0,",numLocalSources_,")");
      AssertLocallyConsistent();
    )
    return localSourceOffsets_[localSource];
}

Int DistGraph::Offset( Int localSource, Int target ) const
{
    DEBUG_ONLY(CSE cse("DistGraph::Offset"))
    if( localSource == END ) localSource = numLocalSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    const Int* targetBuf = LockedTargetBuffer();
    const Int thisOff = SourceOffset(localSource);
    const Int nextOff = SourceOffset(localSource+1);
    auto it = std::lower_bound( targetBuf+thisOff, targetBuf+nextOff, target );
    return it-targetBuf;
}

Int DistGraph::NumConnections( Int localSource ) const
{
    DEBUG_ONLY(
      CSE cse("DistGraph::NumConnections");
      AssertLocallyConsistent();
    )
    if( localSource == END ) localSource = numLocalSources_ - 1;
    return SourceOffset(localSource+1) - SourceOffset(localSource);
}

Int* DistGraph::SourceBuffer() { return sources_.data(); }
Int* DistGraph::TargetBuffer() { return targets_.data(); }
Int* DistGraph::OffsetBuffer() { return localSourceOffsets_.data(); }

const Int* DistGraph::LockedSourceBuffer() const { return sources_.data(); }
const Int* DistGraph::LockedTargetBuffer() const { return targets_.data(); }
const Int* DistGraph::LockedOffsetBuffer() const 
{ return localSourceOffsets_.data(); }

// Auxiliary routines
// ==================
void DistGraph::AssertConsistent() const
{
    Int locallyConsistent = ( locallyConsistent_ ? 1 : 0 );
    Int consistent = 
      mpi::AllReduce( locallyConsistent, mpi::BINARY_OR, Comm() );
    if( !consistent )
        LogicError("DistGraph was not consistent");
}

void DistGraph::AssertLocallyConsistent() const
{
    if( !locallyConsistent_ )
        LogicError("DistGraph was not consistent");
}

void DistGraph::ComputeSourceOffsets()
{
    DEBUG_ONLY(CSE cse("DistGraph::ComputeSourceOffsets"))
    Int sourceOffset = 0;
    Int prevSource = firstLocalSource_-1;
    localSourceOffsets_.resize( numLocalSources_+1 );
    const Int numLocalEdges = NumLocalEdges();
    const Int* sourceBuf = LockedSourceBuffer();
    for( Int e=0; e<numLocalEdges; ++e )
    {
        const Int source = sourceBuf[e];
        DEBUG_ONLY(
          if( source < prevSource )
              RuntimeError("sources were not properly sorted");
        )
        for( ; prevSource<source; ++prevSource )
            localSourceOffsets_[sourceOffset++] = e;
    }
    for( ; sourceOffset<=numLocalSources_; ++sourceOffset )
        localSourceOffsets_[sourceOffset] = numLocalEdges;
}

} // namespace El
