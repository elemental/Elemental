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
#ifndef EL_CORE_DISTGRAPH_DECL_HPP
#define EL_CORE_DISTGRAPH_DECL_HPP

namespace El {

// Use a simple 1d distribution where each process owns a fixed number of 
// sources:
//     if last process,  numSources - (commSize-1)*floor(numSources/commSize)
//     otherwise,        floor(numSources/commSize)
class DistGraph
{
public:
    // Constructors and destructors
    // ============================
    DistGraph();
    DistGraph( mpi::Comm comm );
    DistGraph( Int numVertices, mpi::Comm comm );
    DistGraph( Int numSources, Int numTargets, mpi::Comm comm );
    DistGraph( const Graph& graph );
    // TODO: Move constructor
    DistGraph( const DistGraph& graph );
    ~DistGraph();

    // Assignment and reconfiguration
    // ==============================

    // Making a copy
    // -------------
    const DistGraph& operator=( const Graph& graph );
    const DistGraph& operator=( const DistGraph& graph );
    // TODO: Move assignment

    // Changing the graph size
    // -----------------------
    void Empty( bool clearMemory=true );
    void Resize( Int numVertices );
    void Resize( Int numSources, Int numTargets );

    // Changing the distribution
    // -------------------------
    void SetComm( mpi::Comm comm );

    // Assembly
    // --------
    void Reserve( Int numLocalEdges );

    // Safe edge insertion/removal procedure
    void Connect( Int source, Int target );
    void ConnectLocal( Int localSource, Int target );
    void Disconnect( Int source, Int target ); 
    void DisconnectLocal( Int localSource, Int target );

    // For inserting/removing a sequence of edges and then forcing consistency
    void QueueConnection( Int source, Int target );
    void QueueLocalConnection( Int localSource, Int target ); 
    void QueueDisconnection( Int source, Int target );
    void QueueLocalDisconnection( Int localSource, Int target );
    void MakeConsistent();

    // Queries
    // =======

    // High-level data
    // ---------------
    Int NumSources() const;
    Int NumTargets() const;
    Int FirstLocalSource() const;
    Int NumLocalSources() const;
    Int NumLocalEdges() const;
    Int Capacity() const;
    bool Consistent() const;

    // Distribution information
    // ------------------------
    mpi::Comm Comm() const;
    Int Blocksize() const;

    // Detailed local information
    // --------------------------
    Int Source( Int localEdge ) const;
    Int Target( Int localEdge ) const;
    Int EdgeOffset( Int localSource ) const;
    Int NumConnections( Int localSource ) const;
    Int* SourceBuffer();
    Int* TargetBuffer();
    const Int* LockedSourceBuffer() const;
    const Int* LockedTargetBuffer() const;

private:
    Int numSources_, numTargets_;
    mpi::Comm comm_;

    Int blocksize_;
    Int firstLocalSource_, numLocalSources_;

    vector<Int> sources_, targets_;
    set<pair<Int,Int>> markedForRemoval_;

    // Helpers for local indexing
    bool consistent_;
    vector<Int> localEdgeOffsets_;
    void ComputeEdgeOffsets();

    static bool ComparePairs( const pair<Int,Int>& a, const pair<Int,Int>& b );

    void AssertConsistent() const;

    friend class Graph;
    friend void Copy( const Graph& A, DistGraph& B );
    friend void Copy( const DistGraph& A, Graph& B );
    friend void Copy( const DistGraph& A, DistGraph& B );

    template<typename F> friend class DistSparseMatrix;
    template<typename F> friend struct DistSymmFront;

    template<typename U> friend void Syrk
    ( Orientation orientation,
      U alpha, const DistSparseMatrix<U>& A,
      U beta,        DistSparseMatrix<U>& C, bool conjugate );
};

} // namespace El

#endif // ifndef EL_CORE_DISTGRAPH_DECL_HPP
