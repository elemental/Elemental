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
#ifndef EL_CORE_DISTGRAPH_DECL_HPP
#define EL_CORE_DISTGRAPH_DECL_HPP

#include <set>

namespace El {

struct DistGraphMultMeta
{
    bool ready;
    // NOTE: The 'send' and 'recv' roles reverse for adjoint multiplication
    Int numRecvInds;
    vector<int> sendSizes, sendOffs,
                recvSizes, recvOffs;
    vector<Int> sendInds, colOffs;

    DistGraphMultMeta() : ready(false), numRecvInds(0) { }

    void Clear()
    {
        ready = false;
        numRecvInds = 0;
        SwapClear( sendSizes );
        SwapClear( recvSizes );
        SwapClear( sendOffs );
        SwapClear( recvOffs );
        SwapClear( sendInds );
        SwapClear( colOffs );
    }

    const DistGraphMultMeta& operator=( const DistGraphMultMeta& meta )
    {
        ready = meta.ready;
        numRecvInds = meta.numRecvInds;
        sendSizes = meta.sendSizes;
        sendOffs = meta.sendOffs;
        recvSizes = meta.recvSizes;
        recvOffs = meta.recvOffs;
        sendInds = meta.sendInds;
        colOffs = meta.colOffs;
        return *this;
    }
};



using std::set;

// Forward declare ldl::DistFront
namespace ldl { template<typename F> struct DistFront; }

// Use a simple 1d distribution where each process owns a fixed number of
// sources:
//     if last process,  numSources - (commSize-1)*floor(numSources/commSize)
//     otherwise,        floor(numSources/commSize)
class DistGraph
{
public:
    // Constructors and destructors
    // ============================
    DistGraph( const El::Grid& grid=El::Grid::Default() );
    DistGraph( Int numSources, const El::Grid& grid=El::Grid::Default() );
    DistGraph
    ( Int numSources, Int numTargets,
      const El::Grid& grid=El::Grid::Default() );
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

    // Make a copy of a subgraph
    // -------------------------
    DistGraph operator()( Range<Int> I, Range<Int> J ) const;
    DistGraph operator()( Range<Int> I, const vector<Int>& J ) const;
    DistGraph operator()( const vector<Int>& I, Range<Int> J ) const;
    DistGraph operator()( const vector<Int>& I, const vector<Int>& J ) const;

    // Changing the graph size
    // -----------------------
    void Empty( bool freeMemory=true );
    void Resize( Int numVertices );
    void Resize( Int numSources, Int numTargets );

    // Changing the distribution
    // -------------------------
    void SetGrid( const El::Grid& grid );

    // Assembly
    // --------
    void Reserve( Int numLocalEdges, Int numRemoteEdges=0 );

    // Safe edge insertion/removal procedure
    void Connect( Int source, Int target );
    void ConnectLocal( Int localSource, Int target );
    void Disconnect( Int source, Int target );
    void DisconnectLocal( Int localSource, Int target );

    void FreezeSparsity() EL_NO_EXCEPT;
    void UnfreezeSparsity() EL_NO_EXCEPT;
    bool FrozenSparsity() const EL_NO_EXCEPT;

    // For inserting/removing a sequence of edges and then forcing consistency
    void QueueConnection( Int source, Int target, bool passive=false )
    EL_NO_RELEASE_EXCEPT;
    void QueueLocalConnection( Int localSource, Int target )
    EL_NO_RELEASE_EXCEPT;
    void QueueDisconnection( Int source, Int target, bool passive=false )
    EL_NO_RELEASE_EXCEPT;
    void QueueLocalDisconnection( Int localSource, Int target )
    EL_NO_RELEASE_EXCEPT;
    void ProcessQueues();
    void ProcessLocalQueues();

    // For manually modifying/accessing buffers
    void ForceNumLocalEdges( Int numLocalEdges );
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

    // High-level data
    // ---------------
    Int NumSources() const EL_NO_EXCEPT;
    Int NumTargets() const EL_NO_EXCEPT;
    Int NumEdges() const EL_NO_EXCEPT;
    Int FirstLocalSource() const EL_NO_EXCEPT;
    Int NumLocalSources() const EL_NO_EXCEPT;
    Int NumLocalEdges() const EL_NO_EXCEPT;
    Int Capacity() const EL_NO_EXCEPT;
    bool LocallyConsistent() const EL_NO_EXCEPT;

    // Distribution information
    // ------------------------
    const El::Grid& Grid() const EL_NO_EXCEPT;
    Int Blocksize() const EL_NO_EXCEPT;
    int SourceOwner( Int s ) const EL_NO_RELEASE_EXCEPT;
    Int GlobalSource( Int sLoc ) const EL_NO_RELEASE_EXCEPT;
    Int LocalSource( Int s ) const EL_NO_RELEASE_EXCEPT;

    // Detailed local information
    // --------------------------
    Int Source( Int localEdge ) const EL_NO_RELEASE_EXCEPT;
    Int Target( Int localEdge ) const EL_NO_RELEASE_EXCEPT;
    Int SourceOffset( Int localSource ) const EL_NO_RELEASE_EXCEPT;
    Int Offset( Int localSource, Int target ) const EL_NO_RELEASE_EXCEPT;
    Int NumConnections( Int localSource ) const EL_NO_RELEASE_EXCEPT;

    // Return the ratio of the maximum number of local edges to the
    // total number of edges divided by the number of processes
    double Imbalance() const EL_NO_RELEASE_EXCEPT;

    mutable DistGraphMultMeta multMeta;
    DistGraphMultMeta InitializeMultMeta() const;

    void AssertConsistent() const;
    void AssertLocallyConsistent() const;

private:
    Int numSources_, numTargets_;

    // An observing pointer to a pre-existing Grid.
    const El::Grid* grid_=nullptr;

    Int blocksize_;
    Int numLocalSources_;

    bool frozenSparsity_ = false;
    vector<Int> sources_, targets_;
    set<pair<Int,Int>> markedForRemoval_;

    vector<Int> remoteSources_, remoteTargets_;
    vector<pair<Int,Int>> remoteRemovals_;

    void InitializeLocalData();

    // Helpers for local indexing
    bool locallyConsistent_ = true;
    vector<Int> localSourceOffsets_;

    friend class Graph;
    friend void Copy( const Graph& A, DistGraph& B );
    friend void Copy( const DistGraph& A, Graph& B );
    friend void Copy( const DistGraph& A, DistGraph& B );

    template<typename F> friend class DistSparseMatrix;
};

} // namespace El

#endif // ifndef EL_CORE_DISTGRAPH_HPP
