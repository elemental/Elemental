/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CORE_GRAPH_DECL_HPP
#define EL_CORE_GRAPH_DECL_HPP

namespace El {

// Forward declaration
class DistGraph;
template<typename T>
class SparseMatrix;
template<typename T>
class DistSparseMatrix;

class Graph
{
public:
    // Constructors and destructors
    // ============================
    Graph();
    Graph( Int numVertices );
    Graph( Int numSources, Int numTargets );
    Graph( const Graph& graph );
    // NOTE: This requires the DistGraph to be over a single process
    Graph( const DistGraph& graph );
    // TODO: Move constructor
    ~Graph();

    // Assignment and reconfiguration
    // ==============================

    // Make a copy
    // -----------
    // For copying one graph into another
    const Graph& operator=( const Graph& graph );
    // NOTE: This requires the DistGraph to be over a single process
    const Graph& operator=( const DistGraph& graph );
    // TODO: Move assignment

    // Change the size of the graph
    // ----------------------------
    void Empty( bool clearMemory=true );
    void Resize( Int numVertices );
    void Resize( Int numSources, Int numTargets );

    // Assembly
    // --------
    void Reserve( Int numEdges );

    // Safe (but high overhead) edge insertion/removal procedures
    void Connect( Int source, Int target );
    void Disconnect( Int source, Int target );

    // For appending/removing many edges and then forcing consistency at the end
    void QueueConnection( Int source, Int target );
    void QueueDisconnection( Int source, Int target );
    void MakeConsistent();

    // Queries
    // =======
    Int NumSources() const;
    Int NumTargets() const;
    Int NumEdges() const;
    Int Capacity() const;
    bool Consistent() const;

    Int Source( Int edge ) const;
    Int Target( Int edge ) const;
    Int EdgeOffset( Int source ) const;
    Int NumConnections( Int source ) const;
    Int* SourceBuffer();
    Int* TargetBuffer();
    const Int* LockedSourceBuffer() const;
    const Int* LockedTargetBuffer() const;

private:
    Int numSources_, numTargets_;
    vector<Int> sources_, targets_;
    set<pair<Int,Int>> markedForRemoval_;

    // Helpers for local indexing
    bool consistent_;
    vector<Int> edgeOffsets_;
    void ComputeEdgeOffsets();

    static bool ComparePairs( const pair<Int,Int>& a, const pair<Int,Int>& b );

    void AssertConsistent() const;

    friend class DistGraph;
    template<typename F> friend class SparseMatrix;

    friend void Copy( const Graph& A, Graph& B );
    friend void Copy( const Graph& A, DistGraph& B );
    friend void Copy( const DistGraph& A, Graph& B );

    friend void CopyFromRoot( const DistGraph& GDist, Graph& G );
    template<typename U>
    friend void CopyFromRoot
    ( const DistSparseMatrix<U>& ADist, SparseMatrix<U>& A );
};

} // namespace El

#endif // ifndef EL_CORE_GRAPH_DECL_HPP
