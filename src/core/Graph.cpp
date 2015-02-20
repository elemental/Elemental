/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Constructors and destructors
// ============================

Graph::Graph() : numSources_(0), numTargets_(0), consistent_(true) { }

Graph::Graph( Int numVertices )
: numSources_(numVertices), numTargets_(numVertices), consistent_(true)
{ }

Graph::Graph( Int numSources, Int numTargets )
: numSources_(numSources), numTargets_(numTargets), consistent_(true)
{ }

Graph::Graph( const Graph& graph )
{
    DEBUG_ONLY(CallStackEntry cse("Graph::Graph"))
    if( &graph != this )
        *this = graph;
    else
        LogicError("Tried to construct a graph with itself");
}
    
Graph::Graph( const DistGraph& graph )
{
    DEBUG_ONLY(CallStackEntry cse("Graph::Graph"))
    *this = graph;
}

Graph::~Graph() { }

// Assignment and reconfiguration
// ==============================

// Making a copy
// -------------
const Graph& Graph::operator=( const Graph& graph )
{
    DEBUG_ONLY(CallStackEntry cse("Graph::operator="))
    Copy( graph, *this );
    return *this;
}

const Graph& Graph::operator=( const DistGraph& graph )
{
    DEBUG_ONLY(CallStackEntry cse("Graph::operator="))
    Copy( graph, *this );
    return *this;
}

// Changing the graph size
// -----------------------
void Graph::Empty( bool clearMemory )
{
    numSources_ = 0;
    numTargets_ = 0;
    consistent_ = true;
    if( clearMemory )
    {
        SwapClear( sources_ );
        SwapClear( targets_ );
        SwapClear( edgeOffsets_ );
    }
    else
    {
        sources_.resize( 0 );
        targets_.resize( 0 );
        edgeOffsets_.resize( 0 );
    }
}

void Graph::Resize( Int numVertices )
{ Resize( numVertices, numVertices ); }

void Graph::Resize( Int numSources, Int numTargets )
{
    numSources_ = numSources;
    numTargets_ = numTargets;
    sources_.resize( 0 );
    targets_.resize( 0 );
    edgeOffsets_.resize( 0 );
    consistent_ = true;
}

// Assembly
// --------
void Graph::Reserve( Int numEdges )
{ 
    sources_.reserve( numEdges );
    targets_.reserve( numEdges );
}

void Graph::Connect( Int source, Int target )
{
    DEBUG_ONLY(CallStackEntry cse("Graph::Connect"))
    QueueConnection( source, target );
    MakeConsistent();
}

void Graph::Disconnect( Int source, Int target )
{
    DEBUG_ONLY(CallStackEntry cse("Graph::Disconnect"))
    QueueDisconnection( source, target );
    MakeConsistent();
}

void Graph::QueueConnection( Int source, Int target )
{
    DEBUG_ONLY(
      CallStackEntry cse("Graph::QueueConnection");
      const Int capacity = Capacity();
      const Int numEdges = NumEdges();
      if( numEdges == capacity )
          cerr << "WARNING: Pushing back without first reserving space" << endl;
    )
    if( source < 0 || source >= numSources_ )
        LogicError
        ("Source was out of bounds: ",source," is not in [0,",numSources_,")");
    if( target < 0 || target >= numTargets_ )
        LogicError
        ("Target was out of bounds: ",target," is not in [0,",numTargets_,")");
    sources_.push_back( source );
    targets_.push_back( target );
    consistent_ = false;
}

void Graph::QueueDisconnection( Int source, Int target )
{
    DEBUG_ONLY(CallStackEntry cse("Graph::QueueDisconnection"))
    markedForRemoval_.insert( pair<Int,Int>(source,target) );
    consistent_ = false;
}

void Graph::MakeConsistent()
{
    DEBUG_ONLY(
      CallStackEntry cse("Graph::MakeConsistent");
      if( sources_.size() != targets_.size() )
          LogicError("Inconsistent graph buffer sizes");
    )
    if( !consistent_ )
    {
        const Int numEdges = sources_.size();
        // TODO: Consider switching to using the following by default so that
        //       no extra allocation/memcpy is required
        Int numRemoved=0;
        vector<pair<Int,Int>> pairs( numEdges );
        for( Int e=0; e<numEdges; ++e )
        {
            pair<Int,Int> candidate(sources_[e],targets_[e]);
            if( markedForRemoval_.find(candidate) == markedForRemoval_.end() )
            {
                pairs[e-numRemoved].first = sources_[e];
                pairs[e-numRemoved].second = targets_[e];
            }
            else
            {
                ++numRemoved;
            }
        }
        markedForRemoval_.clear();
        pairs.resize( numEdges-numRemoved );
        std::sort( pairs.begin(), pairs.end(), ComparePairs );

        // Compress out duplicates
        Int lastUnique=0;
        for( Int e=1; e<pairs.size(); ++e )
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

        ComputeEdgeOffsets();

        consistent_ = true;
    }
}

// Queries
// =======

Int Graph::NumSources() const { return numSources_; }
Int Graph::NumTargets() const { return numTargets_; }

Int Graph::NumEdges() const
{
    DEBUG_ONLY(CallStackEntry cse("Graph::NumEdges"))
    return sources_.size();
}

Int Graph::Capacity() const
{
    DEBUG_ONLY(CallStackEntry cse("Graph::Capacity"))
    return Min(sources_.capacity(),targets_.capacity());
}

bool Graph::Consistent() const { return consistent_; }

Int Graph::Source( Int edge ) const
{
    DEBUG_ONLY(
      CallStackEntry cse("Graph::Source");
      if( edge < 0 || edge >= (Int)sources_.size() )
          LogicError("Edge number out of bounds");
    )
    return sources_[edge];
}

Int Graph::Target( Int edge ) const
{
    DEBUG_ONLY(
      CallStackEntry cse("Graph::Target");
      if( edge < 0 || edge >= (Int)targets_.size() )
          LogicError("Edge number out of bounds");
    )
    return targets_[edge];
}

Int Graph::EdgeOffset( Int source ) const
{
    DEBUG_ONLY(
      CallStackEntry cse("Graph::EdgeOffset");
      if( source < 0 )
          LogicError("Negative source index");
      if( source > numSources_ )
          LogicError
          ("Source index was too large: ",source," is not in [0,",
           numSources_,"]");
      AssertConsistent();
    )
    return edgeOffsets_[source];
}

Int Graph::NumConnections( Int source ) const
{
    DEBUG_ONLY(
      CallStackEntry cse("Graph::NumConnections");
      AssertConsistent();
    )
    return EdgeOffset(source+1) - EdgeOffset(source);
}

Int* Graph::SourceBuffer() { return sources_.data(); }
Int* Graph::TargetBuffer() { return targets_.data(); }

const Int* Graph::LockedSourceBuffer() const { return sources_.data(); }
const Int* Graph::LockedTargetBuffer() const { return targets_.data(); }

// Auxiliary functions
// ===================

bool Graph::ComparePairs
( const pair<Int,Int>& a, const pair<Int,Int>& b )
{ return a.first < b.first || (a.first  == b.first && a.second < b.second); }

void Graph::ComputeEdgeOffsets()
{
    DEBUG_ONLY(CallStackEntry cse("Graph::ComputeEdgeOffsets"))
    // Compute the edge offsets
    Int sourceOffset = 0;
    Int prevSource = -1;
    edgeOffsets_.resize( numSources_+1 );
    const Int numEdges = NumEdges();
    for( Int edge=0; edge<numEdges; ++edge )
    {
        const Int source = Source( edge );
        DEBUG_ONLY(
            if( source < prevSource )
                RuntimeError("sources were not properly sorted");
        )
        while( source != prevSource )
        {
            edgeOffsets_[sourceOffset++] = edge;
            ++prevSource;
        }
    }
    edgeOffsets_[numSources_] = numEdges;
}

void Graph::AssertConsistent() const
{ 
    if( !consistent_ )
        LogicError("Graph was not consistent; run MakeConsistent()");
}

} // namespace El
