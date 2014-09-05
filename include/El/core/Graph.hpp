/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CORE_GRAPH_DECL_HPP
#define EL_CORE_GRAPH_DECL_HPP

namespace El {

// Forward declaration
class DistGraph;

class Graph
{
public:
    // Constructors and destructors
    Graph();
    Graph( int numVertices );
    Graph( int numSources, int numTargets );
    Graph( const Graph& graph );
    // NOTE: This requires the DistGraph to be over a single process
    Graph( const DistGraph& graph );
    ~Graph();

    // High-level information
    int NumSources() const;
    int NumTargets() const;

    // Assembly-related routines
    void StartAssembly();
    void StopAssembly();
    void Reserve( int numEdges );
    void Insert( int source, int target );
    int Capacity() const;

    // Data
    int NumEdges() const;
    int Source( int edge ) const;
    int Target( int edge ) const;
    int EdgeOffset( int source ) const;
    int NumConnections( int source ) const;
    int* SourceBuffer();
    int* TargetBuffer();
    const int* LockedSourceBuffer() const;
    const int* LockedTargetBuffer() const;

    // For resizing the graph
    void Empty();
    void Resize( int numVertices );
    void Resize( int numSources, int numTargets );

    // For copying one graph into another
    const Graph& operator=( const Graph& graph );
    // NOTE: This requires the DistGraph to be over a single process
    const Graph& operator=( const DistGraph& graph );

private:
    int numSources_, numTargets_;
    std::vector<int> sources_, targets_;

    // Helpers for local indexing
    bool assembling_, sorted_;
    std::vector<int> edgeOffsets_;
    void ComputeEdgeOffsets();

    static bool ComparePairs
    ( const std::pair<int,int>& a, const std::pair<int,int>& b );

    void EnsureNotAssembling() const;
    void EnsureConsistentSizes() const;
    void EnsureConsistentCapacities() const;

    friend class DistGraph;
    template<typename F> friend class SparseMatrix;
};

} // namespace El

#endif // ifndef EL_CORE_GRAPH_DECL_HPP
