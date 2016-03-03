/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CORE_GRAPH_DECL_HPP
#define EL_CORE_GRAPH_DECL_HPP

#include <set>

namespace El {

using std::set;

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

    // Make a copy of a subgraph
    // -------------------------
    Graph operator()( Range<Int> I, Range<Int> J ) const;
    Graph operator()( Range<Int> I, const vector<Int>& J ) const;
    Graph operator()( const vector<Int>& I, Range<Int> J ) const;
    Graph operator()( const vector<Int>& I, const vector<Int>& J ) const;

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

    void FreezeSparsity() EL_NO_EXCEPT;
    void UnfreezeSparsity() EL_NO_EXCEPT;
    bool FrozenSparsity() const EL_NO_EXCEPT;

    // For appending/removing many edges and then forcing consistency at the end
    void QueueConnection( Int source, Int target ) EL_NO_RELEASE_EXCEPT;
    void QueueDisconnection( Int source, Int target );
    void ProcessQueues();

    // For manually modifying/accessing the buffers
    void ForceNumEdges( Int numEdges );
    void ForceConsistency( bool consistent=true ) EL_NO_EXCEPT;
    Int* SourceBuffer() EL_NO_EXCEPT;
    Int* TargetBuffer() EL_NO_EXCEPT;
    Int* OffsetBuffer() EL_NO_EXCEPT;
    const Int* LockedSourceBuffer() const EL_NO_EXCEPT;
    const Int* LockedTargetBuffer() const EL_NO_EXCEPT;
    const Int* LockedOffsetBuffer() const EL_NO_EXCEPT;
    void ComputeSourceOffsets();

    // Queries
    // =======
    Int NumSources() const EL_NO_EXCEPT;
    Int NumTargets() const EL_NO_EXCEPT;
    Int NumEdges() const EL_NO_EXCEPT;
    Int Capacity() const EL_NO_EXCEPT;
    bool Consistent() const EL_NO_EXCEPT;

    Int Source( Int edge ) const EL_NO_RELEASE_EXCEPT;
    Int Target( Int edge ) const EL_NO_RELEASE_EXCEPT;
    Int SourceOffset( Int source ) const EL_NO_RELEASE_EXCEPT;
    Int Offset( Int source, Int target ) const EL_NO_RELEASE_EXCEPT;
    Int NumConnections( Int source ) const EL_NO_RELEASE_EXCEPT;

    bool EdgeExists( Int source, Int target ) const EL_NO_RELEASE_EXCEPT;

    void AssertConsistent() const;

private:
    Int numSources_, numTargets_;
    bool frozenSparsity_ = false;
    vector<Int> sources_, targets_;
    set<pair<Int,Int>> markedForRemoval_;

    // Helpers for local indexing
    bool consistent_=true;
    vector<Int> sourceOffsets_;

    friend class DistGraph;
    template<typename F> friend class SparseMatrix;

    friend void Copy( const Graph& A, Graph& B );
    friend void Copy( const Graph& A, DistGraph& B );
    friend void Copy( const DistGraph& A, Graph& B );

    friend void CopyFromRoot( const DistGraph& GDist, Graph& G );
    template<typename U>
    friend void CopyFromRoot
    ( const DistSparseMatrix<U>& ADist, SparseMatrix<U>& A );

    template<typename U,typename V>
    friend void EntrywiseMap
    ( const SparseMatrix<U>& A, SparseMatrix<V>& B, function<V(U)> func );
};

} // namespace El

#endif // ifndef EL_CORE_GRAPH_DECL_HPP
