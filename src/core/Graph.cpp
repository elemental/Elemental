/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

Graph::Graph()
: numSources_(0), numTargets_(0), assembling_(false), sorted_(true)
{ }

Graph::Graph( int numVertices )
: numSources_(numVertices), numTargets_(numVertices), 
  assembling_(false), sorted_(true)
{ }

Graph::Graph( int numSources, int numTargets )
: numSources_(numSources), numTargets_(numTargets),
  assembling_(false), sorted_(true)
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

int Graph::NumSources() const { return numSources_; }
int Graph::NumTargets() const { return numTargets_; }

int Graph::NumEdges() const
{
    DEBUG_ONLY(
        CallStackEntry cse("Graph::NumEdges");
        EnsureConsistentSizes();
    )
    return sources_.size();
}

int Graph::Capacity() const
{
    DEBUG_ONLY(
        CallStackEntry cse("Graph::Capacity");
        EnsureConsistentSizes();
        EnsureConsistentCapacities();
    )
    return sources_.capacity();
}

int Graph::Source( int edge ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("Graph::Source");
        if( edge < 0 || edge >= (int)sources_.size() )
            LogicError("Edge number out of bounds");
    )
    EnsureNotAssembling();
    return sources_[edge];
}

int Graph::Target( int edge ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("Graph::Target");
        if( edge < 0 || edge >= (int)targets_.size() )
            LogicError("Edge number out of bounds");
    )
    EnsureNotAssembling();
    return targets_[edge];
}

int Graph::EdgeOffset( int source ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("Graph::EdgeOffset");
        if( source < 0 )
            LogicError("Negative source index");
        if( source > numSources_ )
            LogicError
            ("Source index was too large: ",source," is not in [0,",
             numSources_,"]");
    )
    EnsureNotAssembling();
    return edgeOffsets_[source];
}

int Graph::NumConnections( int source ) const
{
    DEBUG_ONLY(CallStackEntry cse("Graph::NumConnections"))
    return EdgeOffset(source+1) - EdgeOffset(source);
}

int* Graph::SourceBuffer() { return &sources_[0]; }
int* Graph::TargetBuffer() { return &targets_[0]; }

const int* Graph::LockedSourceBuffer() const { return &sources_[0]; }
const int* Graph::LockedTargetBuffer() const { return &targets_[0]; }

const Graph& Graph::operator=( const Graph& graph )
{
    DEBUG_ONLY(CallStackEntry cse("Graph::operator="))
    numSources_ = graph.numSources_;
    numTargets_ = graph.numTargets_;
    sources_ = graph.sources_; 
    targets_ = graph.targets_;

    sorted_ = graph.sorted_;
    assembling_ = graph.assembling_;
    edgeOffsets_ = graph.edgeOffsets_;
    return *this;
}

const Graph& Graph::operator=( const DistGraph& graph )
{
    DEBUG_ONLY(CallStackEntry cse("Graph::operator="))
    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::Size( comm );
    if( commSize != 1 )
        LogicError
        ("Cannot yet construct sequential graph from distributed graph");

    numSources_ = graph.numSources_;
    numTargets_ = graph.numTargets_;
    sources_ = graph.sources_; 
    targets_ = graph.targets_;

    sorted_ = graph.sorted_;
    assembling_ = graph.assembling_;
    edgeOffsets_ = graph.localEdgeOffsets_;
    return *this;
}

bool Graph::ComparePairs
( const std::pair<int,int>& a, const std::pair<int,int>& b )
{ return a.first < b.first || (a.first  == b.first && a.second < b.second); }

void Graph::StartAssembly()
{
    DEBUG_ONLY(CallStackEntry cse("Graph::StartAssembly"))
    EnsureNotAssembling();
    assembling_ = true;
}

void Graph::StopAssembly()
{
    DEBUG_ONLY(CallStackEntry cse("Graph::StopAssembly"))
    if( !assembling_ )
        LogicError("Cannot stop assembly without starting");
    assembling_ = false;

    // Ensure that the connection pairs are sorted
    if( !sorted_ )
    {
        const int numEdges = sources_.size();
        std::vector<std::pair<int,int>> pairs( numEdges );
        for( int e=0; e<numEdges; ++e )
        {
            pairs[e].first = sources_[e];
            pairs[e].second = targets_[e];
        }
        std::sort( pairs.begin(), pairs.end(), ComparePairs );

        // Compress out duplicates
        int lastUnique=0;
        for( int e=1; e<numEdges; ++e )
            if( pairs[e] != pairs[lastUnique] )
                pairs[++lastUnique] = pairs[e];
        const int numUnique = lastUnique+1;

        sources_.resize( numUnique );
        targets_.resize( numUnique );
        for( int e=0; e<numUnique; ++e )
        {
            sources_[e] = pairs[e].first;
            targets_[e] = pairs[e].second;
        }
    }

    ComputeEdgeOffsets();
}

void Graph::ComputeEdgeOffsets()
{
    DEBUG_ONLY(CallStackEntry cse("Graph::ComputeEdgeOffsets"))
    // Compute the edge offsets
    int sourceOffset = 0;
    int prevSource = -1;
    edgeOffsets_.resize( numSources_+1 );
    const int numEdges = NumEdges();
    for( int edge=0; edge<numEdges; ++edge )
    {
        const int source = Source( edge );
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

void Graph::Reserve( int numEdges )
{ 
    sources_.reserve( numEdges );
    targets_.reserve( numEdges );
}

void Graph::Insert( int source, int target )
{
    DEBUG_ONLY(
        CallStackEntry cse("Graph::Insert");
        EnsureConsistentSizes();
        const int capacity = Capacity();
        const int numEdges = NumEdges();
        if( source < 0 || source >= numSources_ )
            LogicError
            ("Source was out of bounds: ",source," is not in [0,",
             numSources_,")");
        if( numEdges == capacity )
            std::cerr << "WARNING: Pushing back without first reserving space" 
                      << std::endl;
    )
    if( !assembling_ )
        LogicError("Must start assembly before pushing back");
    if( sorted_ && sources_.size() != 0 )
    {
        if( source < sources_.back() )
            sorted_ = false;
        if( source == sources_.back() && target < targets_.back() )
            sorted_ = false;
    }
    sources_.push_back( source );
    targets_.push_back( target );
}

void Graph::Empty()
{
    numSources_ = 0;
    numTargets_ = 0;
    SwapClear( sources_ );
    SwapClear( targets_ );
    sorted_ = true;
    assembling_ = false;
    SwapClear( edgeOffsets_ );
}

void Graph::Resize( int numVertices )
{ Resize( numVertices, numVertices ); }

void Graph::Resize( int numSources, int numTargets )
{
    numSources_ = numSources;
    numTargets_ = numTargets;
    SwapClear( sources_ );
    SwapClear( targets_ );
    sorted_ = true;
    assembling_ = false;
    SwapClear( edgeOffsets_ );
}

void Graph::EnsureNotAssembling() const
{
    if( assembling_ )
        LogicError("Should have finished assembling first");
}

void Graph::EnsureConsistentSizes() const
{ 
    if( sources_.size() != targets_.size() )
        LogicError("Inconsistent graph sizes");
}

void Graph::EnsureConsistentCapacities() const
{ 
    if( sources_.capacity() != targets_.capacity() )
        LogicError("Inconsistent graph capacities");
}

} // namespace El
