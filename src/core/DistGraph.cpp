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
#include <El-lite.hpp>

namespace El {

// Constructors and destructors
// ============================
// TODO: Always duplicate the communicator and do not treat mpi::COMM_WORLD
//       and mpi::COMM_SELF as special cases.

DistGraph::DistGraph( mpi::Comm comm )
: numSources_(0), numTargets_(0), 
  commSize_(mpi::Size(comm)), commRank_(mpi::Rank(comm))
{ 
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );
    InitializeLocalData();
}

DistGraph::DistGraph( Int numSources, mpi::Comm comm )
: numSources_(numSources), numTargets_(numSources),
  commSize_(mpi::Size(comm)), commRank_(mpi::Rank(comm))
{ 
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );
    InitializeLocalData();
}

DistGraph::DistGraph( Int numSources, Int numTargets, mpi::Comm comm )
: numSources_(numSources), numTargets_(numTargets),
  commSize_(mpi::Size(comm)), commRank_(mpi::Rank(comm))
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
    DEBUG_CSE
    *this = graph;
}

DistGraph::DistGraph( const DistGraph& graph )
: numSources_(-1), numTargets_(-1), comm_(mpi::COMM_WORLD)
{
    DEBUG_CSE
    if( &graph != this )
        *this = graph;
    DEBUG_ONLY(
      else
          LogicError("Tried to construct DistGraph with itself");
    )
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
    DEBUG_CSE
    Copy( graph, *this );
    return *this;
}

const DistGraph& DistGraph::operator=( const DistGraph& graph )
{
    DEBUG_CSE
    Copy( graph, *this );
    return *this;
}

// Make a copy of a contiguous subgraph
// ------------------------------------
DistGraph DistGraph::operator()( Range<Int> I, Range<Int> J ) const
{
    DEBUG_CSE
    DistGraph subGraph(this->Comm());
    GetSubgraph( *this, I, J, subGraph );
    return subGraph;
}

DistGraph DistGraph::operator()( Range<Int> I, const vector<Int>& J ) const
{
    DEBUG_CSE
    DistGraph subGraph(this->Comm());
    GetSubgraph( *this, I, J, subGraph );
    return subGraph;
}

DistGraph DistGraph::operator()( const vector<Int>& I, Range<Int> J ) const
{
    DEBUG_CSE
    DistGraph subGraph(this->Comm());
    GetSubgraph( *this, I, J, subGraph );
    return subGraph;
}

DistGraph DistGraph::operator()
( const vector<Int>& I, const vector<Int>& J ) const
{
    DEBUG_CSE
    DistGraph subGraph(this->Comm());
    GetSubgraph( *this, I, J, subGraph );
    return subGraph;
}

// Change the graph size
// ---------------------
// TODO: Replace Empty/SoftEmpty in favor of this approach
void DistGraph::Empty( bool freeMemory )
{
    numSources_ = 0;
    numTargets_ = 0;
    numLocalSources_ = 0;
    blocksize_ = 1;
    locallyConsistent_ = true;
    frozenSparsity_ = false;
    if( freeMemory )
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

    SwapClear( remoteSources_ );
    SwapClear( remoteTargets_ );
}

void DistGraph::Resize( Int numVertices ) 
{ Resize( numVertices, numVertices ); }

void DistGraph::Resize( Int numSources, Int numTargets )
{
    if( numSources_ == numSources && numTargets == numTargets_ )
        return;

    frozenSparsity_ = false;

    numSources_ = numSources;
    numTargets_ = numTargets;

    InitializeLocalData();

    SwapClear( remoteSources_ );
    SwapClear( remoteTargets_ );
}

void DistGraph::InitializeLocalData()
{
    blocksize_ = numSources_ / commSize_;
    if( blocksize_*commSize_ < numSources_ || numSources_ == 0 )
        ++blocksize_;

    numLocalSources_ = Min(blocksize_,Max(numSources_-blocksize_*commRank_,0));

    localSourceOffsets_.resize( numLocalSources_+1 );
    for( Int e=0; e<=numLocalSources_; ++e )
        localSourceOffsets_[e] = 0;

    sources_.resize( 0 );
    targets_.resize( 0 );
    locallyConsistent_ = true;
}

// Change the distribution
// -----------------------
void DistGraph::SetComm( mpi::Comm comm )
{
    commSize_ = mpi::Size( comm );
    commRank_ = mpi::Rank( comm );
    if( comm == comm_ )
        return;

    if( comm_ != mpi::COMM_WORLD )
        mpi::Free( comm_ );
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );

    Resize( 0, 0 );
}

// Assembly
// --------
void DistGraph::Reserve( Int numLocalEdges, Int numRemoteEdges )
{ 
    const Int currSize = sources_.size();
    const Int currRemoteSize = remoteSources_.size();
    sources_.reserve( currSize+numLocalEdges );
    targets_.reserve( currSize+numLocalEdges );
    remoteSources_.reserve( currRemoteSize+numRemoteEdges );
    remoteTargets_.reserve( currRemoteSize+numRemoteEdges );
}

void DistGraph::Connect( Int source, Int target )
{
    DEBUG_CSE
    QueueConnection( source, target, true );
    ProcessLocalQueues();
}

void DistGraph::ConnectLocal( Int localSource, Int target )
{
    DEBUG_CSE
    QueueLocalConnection( localSource, target );
    ProcessLocalQueues();
}

void DistGraph::Disconnect( Int source, Int target )
{
    DEBUG_CSE
    QueueDisconnection( source, target, true );
    ProcessLocalQueues();
}

void DistGraph::DisconnectLocal( Int localSource, Int target )
{
    DEBUG_CSE
    QueueLocalDisconnection( localSource, target );
    ProcessLocalQueues();
}

void DistGraph::FreezeSparsity() EL_NO_EXCEPT
{ frozenSparsity_ = true; }
void DistGraph::UnfreezeSparsity() EL_NO_EXCEPT
{ frozenSparsity_ = false; }
bool DistGraph::FrozenSparsity() const EL_NO_EXCEPT
{ return frozenSparsity_; }

void DistGraph::QueueConnection( Int source, Int target, bool passive )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( source == END ) source = numSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    const Int firstLocalSource = blocksize_*commRank_;
    if( source < firstLocalSource || 
        source >= firstLocalSource+numLocalSources_ )
    {
        QueueLocalConnection( source-firstLocalSource, target );
    }
    else if( !passive )
    {
        remoteSources_.push_back( source );
        remoteTargets_.push_back( target );
    }
}

void DistGraph::QueueLocalConnection( Int localSource, Int target )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( NumLocalEdges() == Capacity() )
          cerr << "WARNING: Pushing back without first reserving space" << endl;
    )
    if( localSource == END ) localSource = numLocalSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    DEBUG_ONLY(
      if( localSource < 0 || localSource >= numLocalSources_ )
          LogicError
          ("Local source out of bounds: ",localSource," is not in [0,",
           numLocalSources_,")");
      if( target < 0 || target >= numTargets_ )
          LogicError
          ("Target out of bounds: ",target," is not in [0,",numTargets_,")");
    )
    if( !FrozenSparsity() )
    {
        const Int firstLocalSource = blocksize_*commRank_;
        sources_.push_back( firstLocalSource+localSource );
        targets_.push_back( target );
        locallyConsistent_ = false;
    }
}

void DistGraph::QueueDisconnection( Int source, Int target, bool passive )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    // TODO: Use FrozenSparsity()
    if( source == END ) source = numSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    const Int firstLocalSource = blocksize_*commRank_;
    if( source < firstLocalSource || 
        source >= firstLocalSource+numLocalSources_ )
    {
        QueueLocalDisconnection( source-firstLocalSource, target );
    }
    else if( !passive )
    {
        remoteRemovals_.push_back( pair<Int,Int>(source,target) );
    }
}

void DistGraph::QueueLocalDisconnection( Int localSource, Int target )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    // TODO: Use FrozenSparsity()
    if( localSource == END ) localSource = numLocalSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    DEBUG_ONLY(
      if( localSource < 0 || localSource >= numLocalSources_ )
          LogicError
          ("Local source out of bounds: ",localSource," is not in [0,",
           numLocalSources_,")");
      if( target < 0 || target >= numTargets_ )
          LogicError
          ("Target out of bounds: ",target," is not in [0,",numTargets_,")");
    )
    if( !FrozenSparsity() )
    {
        const Int firstLocalSource = blocksize_*commRank_;
        markedForRemoval_.insert
        ( pair<Int,Int>(firstLocalSource+localSource,target) );
        locallyConsistent_ = false;
    }
    // else throw error?
}

void DistGraph::ProcessQueues()
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( sources_.size() != targets_.size() )
          LogicError("Inconsistent graph buffer sizes");
    )

    // Send the remote edges
    // =====================
    {
        // Compute the send counts
        // -----------------------
        vector<int> sendCounts(commSize_,0);
        for( auto s : remoteSources_ )
            ++sendCounts[SourceOwner(s)];
        // Pack the send data
        // ------------------
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        auto offs = sendOffs;
        vector<Int> sendSources(totalSend), sendTargets(totalSend);
        for( Int i=0; i<totalSend; ++i ) 
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
        const Int totalRecv = recvSources.size();
        for( Int i=0; i<totalRecv; ++i )
            QueueConnection( recvSources[i], recvTargets[i] );
    }

    // Send the remote edge removals
    // =============================
    {
        // Compute the send counts
        // -----------------------
        vector<int> sendCounts(commSize_,0);
        const Int numRemoteRemovals = remoteRemovals_.size();
        for( Int i=0; i<numRemoteRemovals; ++i )
            ++sendCounts[SourceOwner(remoteRemovals_[i].first)];
        // Pack the send data
        // ------------------
        vector<int> sendOffs;
        const int totalSend = Scan( sendCounts, sendOffs );
        auto offs = sendOffs;
        vector<Int> sendSources(totalSend), sendTargets(totalSend);
        for( Int i=0; i<totalSend; ++i ) 
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
        const Int totalRecv = recvSources.size();
        for( Int i=0; i<totalRecv; ++i )
            QueueDisconnection( recvSources[i], recvTargets[i] );
    }

    // Ensure that the kept local edges are sorted and unique
    // ======================================================
    ProcessLocalQueues();
}

void DistGraph::ProcessLocalQueues()
{
    DEBUG_CSE
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
Int DistGraph::NumSources() const EL_NO_EXCEPT
{ return numSources_; }
Int DistGraph::NumTargets() const EL_NO_EXCEPT
{ return numTargets_; }
Int DistGraph::NumEdges() const EL_NO_EXCEPT
{ return mpi::AllReduce( NumLocalEdges(), Comm() ); }
Int DistGraph::FirstLocalSource() const EL_NO_EXCEPT
{ return blocksize_*commRank_; }
Int DistGraph::NumLocalSources() const EL_NO_EXCEPT
{ return numLocalSources_; }

Int DistGraph::NumLocalEdges() const EL_NO_EXCEPT
{ return sources_.size(); }

Int DistGraph::Capacity() const EL_NO_EXCEPT
{ return Min(sources_.capacity(),targets_.capacity()); }

bool DistGraph::LocallyConsistent() const EL_NO_EXCEPT
{ return locallyConsistent_; }

// Distribution information
// ------------------------
mpi::Comm DistGraph::Comm() const EL_NO_EXCEPT
{ return comm_; }
Int DistGraph::Blocksize() const EL_NO_EXCEPT
{ return blocksize_; }

int DistGraph::SourceOwner( Int source ) const EL_NO_RELEASE_EXCEPT
{ 
    if( source == END ) source = numSources_ - 1;
    return source / blocksize_;
}

Int DistGraph::GlobalSource( Int sLoc ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( sLoc == END ) sLoc = numLocalSources_ - 1;
    DEBUG_ONLY(
      if( sLoc < 0 || sLoc >= NumLocalSources() )
          LogicError("Invalid local source index");
    )
    return sLoc + FirstLocalSource();
}

Int DistGraph::LocalSource( Int s ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( s == END ) s = numSources_ - 1;
    DEBUG_ONLY(
      if( s < 0 || s >= NumSources() )
          LogicError("Invalid global source index");
    )
    return s - FirstLocalSource();
}

// Detailed local information
// --------------------------
Int DistGraph::Source( Int localEdge ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( localEdge < 0 || localEdge >= (Int)sources_.size() )
          LogicError("Edge number out of bounds");
    )
    return sources_[localEdge];
}

Int DistGraph::Target( Int localEdge ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( localEdge < 0 || localEdge >= (Int)targets_.size() )
          LogicError("Edge number out of bounds");
    )
    return targets_[localEdge];
}

Int DistGraph::SourceOffset( Int localSource ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( localSource == END ) localSource = numLocalSources_ - 1;
    DEBUG_ONLY(
      if( localSource < 0 || localSource > numLocalSources_ )
          LogicError
          ("Out of bounds localSource: ",localSource,
           " is not in [0,",numLocalSources_,")");
      AssertLocallyConsistent();
    )
    return localSourceOffsets_[localSource];
}

Int DistGraph::Offset( Int localSource, Int target ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( localSource == END ) localSource = numLocalSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    const Int* targetBuf = LockedTargetBuffer();
    const Int thisOff = SourceOffset(localSource);
    const Int nextOff = SourceOffset(localSource+1);
    auto it = std::lower_bound( targetBuf+thisOff, targetBuf+nextOff, target );
    return it-targetBuf;
}

Int DistGraph::NumConnections( Int localSource ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(AssertLocallyConsistent())
    if( localSource == END ) localSource = numLocalSources_ - 1;
    return SourceOffset(localSource+1) - SourceOffset(localSource);
}

double DistGraph::Imbalance() const
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Int numLocalEdges = NumLocalEdges();
    double numEdges = mpi::AllReduce( numLocalEdges, comm_ );
    double maxLocalEdges = mpi::AllReduce( numLocalEdges, mpi::MAX, comm_ );
    double commSize = commSize_;
    return (Max(maxLocalEdges,1)*commSize)/Max(numEdges,1);
}

Int* DistGraph::SourceBuffer() EL_NO_EXCEPT
{ return sources_.data(); }
Int* DistGraph::TargetBuffer() EL_NO_EXCEPT
{ return targets_.data(); }
Int* DistGraph::OffsetBuffer() EL_NO_EXCEPT
{ return localSourceOffsets_.data(); }

const Int* DistGraph::LockedSourceBuffer() const EL_NO_EXCEPT
{ return sources_.data(); }
const Int* DistGraph::LockedTargetBuffer() const EL_NO_EXCEPT
{ return targets_.data(); }
const Int* DistGraph::LockedOffsetBuffer() const EL_NO_EXCEPT
{ return localSourceOffsets_.data(); }

void DistGraph::ForceNumLocalEdges( Int numLocalEdges )
{
    DEBUG_CSE
    sources_.resize( numLocalEdges );
    targets_.resize( numLocalEdges );
    locallyConsistent_ = false;
}

void DistGraph::ForceConsistency( bool consistent ) EL_NO_EXCEPT
{ locallyConsistent_ = consistent; }

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
DistGraphMultMeta DistGraph::InitializeMultMeta() const
{
    DEBUG_ONLY(CSE cse("DistSparseMatrix::InitializeMultMeta"))
    if( multMeta.ready )
        return multMeta;
    mpi::Comm comm = Comm();
    const int commSize = commSize_;
    auto& meta = multMeta;

    // Compute the set of row indices that we need from X in a normal
    // multiply or update of Y in the adjoint case
    const Int* colBuffer = LockedTargetBuffer();
    const Int numLocalEntries = NumLocalEdges();
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
    Int vecBlocksize = NumTargets() / commSize;
    if( vecBlocksize*commSize < NumTargets() || NumTargets() == 0 )
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

void DistGraph::ComputeSourceOffsets()
{
    DEBUG_CSE
    Int sourceOffset = 0;
    Int prevSource = blocksize_*commRank_-1;
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
