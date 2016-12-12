/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2012 Jack Poulson, Lexing Ying, and
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013 Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2014 Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>

namespace El {

// Constructors and destructors
// ============================

Graph::Graph() : numSources_(0), numTargets_(0) { }

Graph::Graph( Int numVertices )
: numSources_(numVertices), numTargets_(numVertices)
{
    sourceOffsets_.resize( numSources_+1 );
    for( Int e=0; e<=numSources_; ++e )
        sourceOffsets_[e] = 0;
}

Graph::Graph( Int numSources, Int numTargets )
: numSources_(numSources), numTargets_(numTargets)
{
    sourceOffsets_.resize( numSources_+1 );
    for( Int e=0; e<=numSources_; ++e )
        sourceOffsets_[e] = 0;
}

Graph::Graph( const Graph& graph )
: numSources_(-1), numTargets_(-1)
{
    EL_DEBUG_CSE
    if( &graph != this )
        *this = graph;
    else
        LogicError("Tried to construct a graph with itself");
}

Graph::Graph( const DistGraph& graph )
: numSources_(-1), numTargets_(-1)
{
    EL_DEBUG_CSE
    *this = graph;
}

Graph::~Graph() { }

// Assignment and reconfiguration
// ==============================

// Making a copy
// -------------
const Graph& Graph::operator=( const Graph& graph )
{
    EL_DEBUG_CSE
    Copy( graph, *this );
    return *this;
}

const Graph& Graph::operator=( const DistGraph& graph )
{
    EL_DEBUG_CSE
    Copy( graph, *this );
    return *this;
}

// Make a copy of a contiguous subgraph
// ------------------------------------
Graph Graph::operator()( Range<Int> I, Range<Int> J ) const
{
    EL_DEBUG_CSE
    Graph subGraph;
    GetSubgraph( *this, I, J, subGraph );
    return subGraph;
}

Graph Graph::operator()( Range<Int> I, const vector<Int>& J ) const
{
    EL_DEBUG_CSE
    Graph subGraph;
    GetSubgraph( *this, I, J, subGraph );
    return subGraph;
}

Graph Graph::operator()( const vector<Int>& I, Range<Int> J ) const
{
    EL_DEBUG_CSE
    Graph subGraph;
    GetSubgraph( *this, I, J, subGraph );
    return subGraph;
}

Graph Graph::operator()( const vector<Int>& I, const vector<Int>& J ) const
{
    EL_DEBUG_CSE
    Graph subGraph;
    GetSubgraph( *this, I, J, subGraph );
    return subGraph;
}

// Changing the graph size
// -----------------------
void Graph::Empty( bool clearMemory )
{
    numSources_ = 0;
    numTargets_ = 0;
    consistent_ = true;
    frozenSparsity_ = false;
    if( clearMemory )
    {
        SwapClear( sources_ );
        SwapClear( targets_ );
        SwapClear( sourceOffsets_ );
    }
    else
    {
        sources_.resize( 0 );
        targets_.resize( 0 );
    }
    sourceOffsets_.resize( 1 );
    sourceOffsets_[0] = 0;
}

void Graph::Resize( Int numVertices )
{ Resize( numVertices, numVertices ); }

void Graph::Resize( Int numSources, Int numTargets )
{
    EL_DEBUG_CSE
    if( numSources_ == numSources && numTargets_ == numTargets )
        return;

    frozenSparsity_ = false;

    numSources_ = numSources;
    numTargets_ = numTargets;
    sources_.resize( 0 );
    targets_.resize( 0 );
    sourceOffsets_.resize( numSources+1 );
    for( Int e=0; e<=numSources; ++e )
        sourceOffsets_[e] = 0;
    consistent_ = true;
}

// Assembly
// --------
void Graph::Reserve( Int numEdges )
{
    const Int currSize = sources_.size();
    sources_.reserve( currSize+numEdges );
    targets_.reserve( currSize+numEdges );
}

void Graph::Connect( Int source, Int target )
{
    EL_DEBUG_CSE
    QueueConnection( source, target );
    ProcessQueues();
}

void Graph::Disconnect( Int source, Int target )
{
    EL_DEBUG_CSE
    QueueDisconnection( source, target );
    ProcessQueues();
}

void Graph::FreezeSparsity() EL_NO_EXCEPT { frozenSparsity_ = true; }
void Graph::UnfreezeSparsity() EL_NO_EXCEPT { frozenSparsity_ = false; }
bool Graph::FrozenSparsity() const EL_NO_EXCEPT { return frozenSparsity_; }

void Graph::QueueConnection( Int source, Int target )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( NumEdges() == Capacity() )
          cerr << "WARNING: Pushing back without first reserving space" << endl;
    )
    if( source == END ) source = numSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    EL_DEBUG_ONLY(
      if( source < 0 || source >= numSources_ )
          LogicError
          ("Source out of bounds: ",source," not in [0,",numSources_,")");
      if( target < 0 || target >= numTargets_ )
          LogicError
          ("Target out of bounds: ",target," not in [0,",numTargets_,")");
    )
    if( !FrozenSparsity() )
    {
        sources_.push_back( source );
        targets_.push_back( target );
        consistent_ = false;
    }
}

void Graph::QueueDisconnection( Int source, Int target )
{
    EL_DEBUG_CSE
    if( source == END ) source = numSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    if( !FrozenSparsity() )
    {
        markedForRemoval_.insert( pair<Int,Int>(source,target) );
        consistent_ = false;
    }
}

void Graph::ProcessQueues()
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( sources_.size() != targets_.size() )
          LogicError("Inconsistent graph buffer sizes");
    )
    if( consistent_ )
        return;

    Int numRemoved=0;
    const Int numEdges = sources_.size();
    vector<pair<Int,Int>> pairs( numEdges );
    if( markedForRemoval_.size() != 0 )
    {
        for( Int e=0; e<numEdges; ++e )
        {
            pair<Int,Int> candidate(sources_[e],targets_[e]);
            if( markedForRemoval_.find(candidate) == markedForRemoval_.end() )
                pairs[e-numRemoved] = candidate;
            else
                ++numRemoved;
        }
        markedForRemoval_.clear();
        pairs.resize( numEdges-numRemoved );
    }
    else
    {
        for( Int e=0; e<numEdges; ++e )
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
    consistent_ = true;
}

// Queries
// =======

Int Graph::NumSources() const EL_NO_EXCEPT { return numSources_; }
Int Graph::NumTargets() const EL_NO_EXCEPT { return numTargets_; }

Int Graph::NumEdges() const EL_NO_EXCEPT
{
    EL_DEBUG_CSE
    return sources_.size();
}

Int Graph::Capacity() const EL_NO_EXCEPT
{
    EL_DEBUG_CSE
    return Min(sources_.capacity(),targets_.capacity());
}

bool Graph::Consistent() const EL_NO_EXCEPT { return consistent_; }

Int Graph::Source( Int edge ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( edge < 0 || edge >= (Int)sources_.size() )
          LogicError("Edge number out of bounds");
    )
    return sources_[edge];
}

Int Graph::Target( Int edge ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( edge < 0 || edge >= (Int)targets_.size() )
          LogicError("Edge number out of bounds");
    )
    return targets_[edge];
}

Int Graph::SourceOffset( Int source ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    if( source == END ) source = numSources_ - 1;
    EL_DEBUG_ONLY(
      if( source < 0 )
          LogicError("Negative source index");
      if( source > numSources_ )
          LogicError
          ("Source index was too large: ",source," is not in [0,",
           numSources_,"]");
      AssertConsistent();
    )
    return sourceOffsets_[source];
}

Int Graph::Offset( Int source, Int target ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    if( source == END ) source = numSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    const Int* targetBuf = LockedTargetBuffer();
    const Int thisOff = SourceOffset(source);
    const Int nextOff = SourceOffset(source+1);
    auto it = std::lower_bound( targetBuf+thisOff, targetBuf+nextOff, target );
    return it-targetBuf;
}

bool Graph::EdgeExists( Int source, Int target ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    if( source == END ) source = numSources_ - 1;
    if( target == END ) target = numTargets_ - 1;
    Int index = Offset( source, target );
    if( Source(index) != source || Target(index) != target )
        return false;
    else
        return true;
}

Int Graph::NumConnections( Int source ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    if( source == END ) source = numSources_ - 1;
    EL_DEBUG_ONLY(AssertConsistent())
    return SourceOffset(source+1) - SourceOffset(source);
}

Int* Graph::SourceBuffer() EL_NO_EXCEPT { return sources_.data(); }
Int* Graph::TargetBuffer() EL_NO_EXCEPT { return targets_.data(); }
Int* Graph::OffsetBuffer() EL_NO_EXCEPT { return sourceOffsets_.data(); }

void Graph::ForceNumEdges( Int numEdges )
{
    EL_DEBUG_CSE
    sources_.resize( numEdges );
    targets_.resize( numEdges );
    consistent_ = false;
}

void Graph::ForceConsistency( bool consistent ) EL_NO_EXCEPT
{ consistent_ = consistent; }

const Int* Graph::LockedSourceBuffer() const EL_NO_EXCEPT
{ return sources_.data(); }
const Int* Graph::LockedTargetBuffer() const EL_NO_EXCEPT
{ return targets_.data(); }
const Int* Graph::LockedOffsetBuffer() const EL_NO_EXCEPT
{ return sourceOffsets_.data(); }

// Auxiliary functions
// ===================

void Graph::ComputeSourceOffsets()
{
    EL_DEBUG_CSE
    Int sourceOffset = 0;
    Int prevSource = -1;
    sourceOffsets_.resize( numSources_+1 );
    const Int numEdges = NumEdges();
    const Int* sourceBuf = LockedSourceBuffer();
    for( Int e=0; e<numEdges; ++e )
    {
        const Int source = sourceBuf[e];
        EL_DEBUG_ONLY(
          if( source < prevSource )
              RuntimeError("sources were not properly sorted");
        )
        for( ; prevSource<source; ++prevSource )
            sourceOffsets_[sourceOffset++] = e;
    }
    for( ; sourceOffset<=numSources_; ++sourceOffset )
        sourceOffsets_[sourceOffset] = numEdges;
}

void Graph::AssertConsistent() const
{
    if( !consistent_ )
        LogicError("Graph was not consistent; run ProcessQueues()");
}

} // namespace El
