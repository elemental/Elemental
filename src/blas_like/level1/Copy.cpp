/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>

namespace El {

void Copy( const Graph& A, Graph& B )
{
    EL_DEBUG_CSE
    const Int numSources = A.NumSources();
    const Int numTargets = A.NumTargets();

    B.Resize( numSources, numTargets );
    // Directly assign instead of queueing up the individual edges
    B.sources_ = A.sources_;
    B.targets_ = A.targets_;
    B.consistent_ = A.consistent_;
    B.sourceOffsets_ = A.sourceOffsets_;
    B.ProcessQueues();
}

void Copy( const Graph& A, DistGraph& B )
{
    EL_DEBUG_CSE
    const Int numSources = A.NumSources();
    const Int numTargets = A.NumTargets();

    B.SetGrid( Grid::Trivial() );
    B.Resize( numSources, numTargets );
    // Directly assign instead of queueing up the individual edges
    B.sources_ = A.sources_;
    B.targets_ = A.targets_;
    B.locallyConsistent_ = A.consistent_;
    B.localSourceOffsets_ = A.sourceOffsets_;
    B.ProcessLocalQueues();
}

void Copy( const DistGraph& A, Graph& B )
{
    EL_DEBUG_CSE
    const Int numSources = A.NumSources();
    const Int numTargets = A.NumTargets();
    if( A.Grid().Size() != 1 )
        LogicError("Cannot yet construct sequential graph from distributed");

    B.Resize( numSources, numTargets );
    // Directly assign instead of queueing up the individual edges
    B.sources_ = A.sources_;
    B.targets_ = A.targets_;
    B.consistent_ = A.locallyConsistent_;
    B.sourceOffsets_ = A.localSourceOffsets_;
    B.ProcessQueues();
}

void Copy( const DistGraph& A, DistGraph& B )
{
    EL_DEBUG_CSE
    const Int numSources = A.NumSources();
    const Int numTargets = A.NumTargets();

    B.SetGrid( A.Grid() );
    B.Resize( numSources, numTargets );
    // Directly assign instead of queueing up the individual edges
    B.sources_ = A.sources_;
    B.targets_ = A.targets_;
    B.multMeta = A.multMeta;
    B.locallyConsistent_ = A.locallyConsistent_;
    B.localSourceOffsets_ = A.localSourceOffsets_;
    B.ProcessLocalQueues();
}

void CopyFromRoot( const DistGraph& distGraph, Graph& graph )
{
    EL_DEBUG_CSE
    const Grid& grid = distGraph.Grid();
    const int commSize = grid.Size();
    const int commRank = grid.Rank();

    const int numLocalEdges = distGraph.NumLocalEdges();
    vector<int> edgeSizes(commSize);
    mpi::AllGather( &numLocalEdges, 1, edgeSizes.data(), 1, grid.Comm() );
    vector<int> edgeOffsets;
    const int numEdges = Scan( edgeSizes, edgeOffsets );

    graph.Resize( distGraph.NumSources(), distGraph.NumTargets() );
    graph.Reserve( numEdges );
    graph.sources_.resize( numEdges );
    graph.targets_.resize( numEdges );
    mpi::Gather
    ( distGraph.LockedSourceBuffer(), numLocalEdges,
      graph.SourceBuffer(), edgeSizes.data(), edgeOffsets.data(),
      commRank, grid.Comm() );
    mpi::Gather
    ( distGraph.LockedTargetBuffer(), numLocalEdges,
      graph.TargetBuffer(), edgeSizes.data(), edgeOffsets.data(),
      commRank, grid.Comm() );
    graph.ProcessQueues();
}

void CopyFromNonRoot( const DistGraph& distGraph, int root )
{
    EL_DEBUG_CSE
    const Grid& grid = distGraph.Grid();
    const int commSize = grid.Size();
    const int commRank = grid.Rank();
    if( commRank == root )
        LogicError("Root called CopyFromNonRoot");

    const int numLocalEdges = distGraph.NumLocalEdges();
    vector<int> edgeSizes(commSize);
    mpi::AllGather( &numLocalEdges, 1, edgeSizes.data(), 1, grid.Comm() );
    vector<int> edgeOffsets;
    Scan( edgeSizes, edgeOffsets );

    mpi::Gather
    ( distGraph.LockedSourceBuffer(), numLocalEdges,
      (Int*)0, edgeSizes.data(), edgeOffsets.data(), root, grid.Comm() );
    mpi::Gather
    ( distGraph.LockedTargetBuffer(), numLocalEdges,
      (Int*)0, edgeSizes.data(), edgeOffsets.data(), root, grid.Comm() );
}

} // namespace El
