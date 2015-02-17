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

#ifdef EL_HAVE_METIS
#include "metis.h"
#endif

namespace El {

inline void
NestedDissectionRecursion
( const Graph& graph, 
  const vector<Int>& perm,
        Separator& sep, 
        SymmNode& node,
        Int off, 
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NestedDissectionRecursion"))
    if( graph.NumSources() <= ctrl.cutoff )
    {
        // Fill in this node of the local separator tree
        const Int numSources = graph.NumSources();
        sep.off = off;
        sep.inds = perm;
        // TODO: Replace with better deletion mechanism
        SwapClear( sep.children );

        // Fill in this node of the local elimination tree
        node.size = numSources;
        node.off = off;
        // TODO: Replace with better deletion mechanism
        SwapClear( node.children );
        set<Int> lowerStruct;
        for( Int s=0; s<node.size; ++s )
        {
            const Int numConnections = graph.NumConnections( s );
            const Int edgeOff = graph.EdgeOffset( s );
            for( Int t=0; t<numConnections; ++t )
            {
                const Int target = graph.Target( edgeOff+t );
                if( target >= numSources )
                    lowerStruct.insert( off+target );
            }
        }
        CopySTL( lowerStruct, node.lowerStruct );
    }
    else
    {
        // Partition the graph and construct the inverse map
        Graph leftChild, rightChild;
        vector<Int> map;
        const Int sepSize = Bisect( graph, leftChild, rightChild, map, ctrl );
        const Int numSources = graph.NumSources();
        vector<Int> invMap( numSources );
        for( Int s=0; s<numSources; ++s )
            invMap[map[s]] = s;

        // Mostly compute this node of the local separator tree
        // (we will finish computing the separator indices soon)
        sep.off = off + (numSources-sepSize);
        sep.inds.resize( sepSize );
        for( Int s=0; s<sepSize; ++s )
        {
            const Int mappedSource = s + (numSources-sepSize);
            sep.inds[s] = invMap[mappedSource];
        }
    
        // Fill in this node in the local elimination tree
        node.size = sepSize;
        node.off = sep.off;
        set<Int> lowerStruct;
        for( Int s=0; s<sepSize; ++s )
        {
            const Int source = sep.inds[s];
            const Int numConnections = graph.NumConnections( source );
            const Int edgeOff = graph.EdgeOffset( source );
            for( Int t=0; t<numConnections; ++t )
            {
                const Int target = graph.Target( edgeOff+t );
                if( target >= numSources )
                    lowerStruct.insert( off+target );
            }
        }
        CopySTL( lowerStruct, node.lowerStruct );

        // Finish computing the separator indices
        for( Int s=0; s<sepSize; ++s )
            sep.inds[s] = perm[sep.inds[s]];

        // Construct the inverse maps from the child indices to the original
        // degrees of freedom
        const Int leftChildSize = leftChild.NumSources();
        vector<Int> leftPerm( leftChildSize );
        for( Int s=0; s<leftChildSize; ++s )
            leftPerm[s] = perm[invMap[s]];
        const Int rightChildSize = rightChild.NumSources();
        vector<Int> rightPerm( rightChildSize );
        for( Int s=0; s<rightChildSize; ++s )
            rightPerm[s] = perm[invMap[s+leftChildSize]];

        sep.children.resize( 2 );
        node.children.resize( 2 );
        sep.children[0] = new Separator(&sep);
        sep.children[1] = new Separator(&sep);
        node.children[0] = new SymmNode(&node);
        node.children[1] = new SymmNode(&node);
        NestedDissectionRecursion
        ( leftChild, leftPerm, *sep.children[0], *node.children[0], 
          off, ctrl );
        NestedDissectionRecursion
        ( rightChild, rightPerm, *sep.children[1], *node.children[1], 
          off+leftChildSize, ctrl );
    }
}

inline void
NestedDissectionRecursion
( const DistGraph& graph, 
  const DistMap& perm,
        DistSeparator& sep, 
        DistSymmNode& node,
        Int off, 
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NestedDissectionRecursion"))
    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::Size(comm);

    mpi::Dup( comm, sep.comm );
    mpi::Dup( comm, node.comm );

    if( commSize > 1 )
    {
        // Partition the graph and construct the inverse map
        DistGraph child;
        bool childIsOnLeft;
        DistMap map;
        const Int sepSize = Bisect( graph, child, map, childIsOnLeft, ctrl );
        const Int numSources = graph.NumSources();
        const Int childSize = child.NumSources();
        const Int leftChildSize = 
            ( childIsOnLeft ? childSize : numSources-sepSize-childSize );

        DistMap invMap;
        InvertMap( map, invMap );

        // Mostly fill this node of the DistSeparatorTree
        // (we will finish computing the separator indices at the end)
        sep.off = off + (numSources-sepSize);
        sep.inds.resize( sepSize );
        for( Int s=0; s<sepSize; ++s )
            sep.inds[s] = s + (numSources-sepSize);
        invMap.Translate( sep.inds );

        // Fill in this node of the DistSymmElimTree
        node.size = sepSize;
        node.off = sep.off;
        const Int numLocalSources = graph.NumLocalSources();
        const Int firstLocalSource = graph.FirstLocalSource();
        set<Int> localLowerStruct;
        for( Int s=0; s<sepSize; ++s )
        {
            const Int source = sep.inds[s];
            if( source >= firstLocalSource && 
                source < firstLocalSource+numLocalSources )
            {
                const Int localSource = source - firstLocalSource;
                const Int numConnections = graph.NumConnections( localSource );
                const Int localOff = graph.EdgeOffset( localSource );
                for( Int t=0; t<numConnections; ++t )
                {
                    const Int target = graph.Target( localOff+t );
                    if( target >= numSources )
                        localLowerStruct.insert( off+target );
                }
            }
        }
        const int numLocalConnected = localLowerStruct.size();
        vector<int> localConnectedSizes( commSize );
        mpi::AllGather
        ( &numLocalConnected, 1, localConnectedSizes.data(), 1, comm );
        vector<Int> localConnectedVec;
        CopySTL( localLowerStruct, localConnectedVec );
        vector<int> localConnectedOffs;
        const int sumOfLocalConnectedSizes = 
            Scan( localConnectedSizes, localConnectedOffs );
        vector<Int> localConnections( sumOfLocalConnectedSizes );
        mpi::AllGather
        ( localConnectedVec.data(), numLocalConnected,
          localConnections.data(), 
          localConnectedSizes.data(), localConnectedOffs.data(), comm );
        set<Int> lowerStruct
        ( localConnections.begin(), localConnections.end() );
        CopySTL( lowerStruct, node.lowerStruct );

        // Finish computing the separator indices
        perm.Translate( sep.inds );

        // Construct map from child indices to the original ordering
        DistMap newPerm( child.NumSources(), child.Comm() );
        const Int localChildSize = child.NumLocalSources();
        const Int firstLocalChildSource = child.FirstLocalSource();
        if( childIsOnLeft )
            for( Int s=0; s<localChildSize; ++s )
                newPerm.SetLocal( s, s+firstLocalChildSource );
        else
            for( Int s=0; s<localChildSize; ++s )
                newPerm.SetLocal( s, s+firstLocalChildSource+leftChildSize );
        invMap.Extend( newPerm );
        perm.Extend( newPerm );

        // Recurse
        const Int childOff = ( childIsOnLeft ? off : off+leftChildSize );
        sep.child = new DistSeparator(&sep);
        node.child = new DistSymmNode(&node);
        node.child->onLeft = childIsOnLeft;
        NestedDissectionRecursion
        ( child, newPerm, *sep.child, *node.child, childOff, ctrl );
    }
    else
    {
        Graph seqGraph( graph );

        sep.duplicate = new Separator(&sep);
        node.duplicate = new SymmNode(&node);

        NestedDissectionRecursion
        ( seqGraph, perm.Map(), *sep.duplicate, *node.duplicate, off, ctrl );

        // Pull information up from the duplicates
        sep.off = sep.duplicate->off;
        sep.inds = sep.duplicate->inds;
        node.size = node.duplicate->size;
        node.off = node.duplicate->off;
        node.lowerStruct = node.duplicate->lowerStruct;
    }
}

void NestedDissection
( const Graph& graph, 
        vector<Int>& map,
        Separator& sep, 
        SymmNodeInfo& info,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NestedDissection"))
    // NOTE: There is a potential memory leak here if sep or info is reused

    const Int numSources = graph.NumSources();
    vector<Int> perm(numSources);
    for( Int s=0; s<numSources; ++s )
        perm[s] = s;

    SymmNode node;
    NestedDissectionRecursion( graph, perm, sep, node, 0, ctrl );

    // Construct the distributed reordering    
    BuildMap( sep, map );
    DEBUG_ONLY(EnsurePermutation(map))

    // Run the symbolic analysis
    SymmetricAnalysis( node, info );
}

void NestedDissection
( const DistGraph& graph, 
        DistMap& map,
        DistSeparator& sep, 
        DistSymmNodeInfo& info,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NestedDissection"))
    // NOTE: There is a potential memory leak here if sep or info is reused

    DistMap perm( graph.NumSources(), graph.Comm() );
    const Int firstLocalSource = perm.FirstLocalSource();
    const Int numLocalSources = perm.NumLocalSources();
    for( Int s=0; s<numLocalSources; ++s )
        perm.SetLocal( s, s+firstLocalSource );

    DistSymmNode node;
    NestedDissectionRecursion( graph, perm, sep, node, 0, ctrl );

    // Construct the distributed reordering    
    BuildMap( sep, map );
    DEBUG_ONLY(EnsurePermutation(map))

    // Run the symbolic analysis
    SymmetricAnalysis( node, info, ctrl.storeFactRecvInds );
}

Int Bisect
( const Graph& graph, Graph& leftChild, Graph& rightChild,
  vector<Int>& perm, const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Bisect"))
#ifdef EL_HAVE_METIS
    // METIS assumes that there are no self-connections or connections 
    // outside the sources, so we must manually remove them from our graph
    const Int numSources = graph.NumSources();
    const Int numEdges = graph.NumEdges();
    Int numValidEdges = 0;
    for( Int i=0; i<numEdges; ++i )
        if( graph.Source(i) != graph.Target(i) && graph.Target(i) < numSources )
            ++numValidEdges;

    // Fill our connectivity (ignoring self and too-large connections)
    vector<idx_t> xAdj( numSources+1 );
    vector<idx_t> adjacency( numValidEdges );
    Int validCounter=0;
    Int sourceOff=0;
    Int prevSource=-1;
    for( Int edge=0; edge<numEdges; ++edge )
    {
        const Int source = graph.Source( edge );
        const Int target = graph.Target( edge );
        DEBUG_ONLY(
          if( source < prevSource )
              RuntimeError("sources were not properly sorted");
        )
        while( source != prevSource )
        {
            xAdj[sourceOff++] = validCounter;
            ++prevSource;
        }
        if( source != target && target < numSources )
        {
            adjacency[validCounter] = target;
            ++validCounter;
        }
    }
    DEBUG_ONLY(
      if( sourceOff != numSources )
          LogicError("Mistake in xAdj computation");
    )
    xAdj[numSources] = numValidEdges;

    // Call METIS_ComputeVertexSeparator, which is meant to be used by ParMETIS
    idx_t nvtxs = numSources;
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions( options );
    options[METIS_OPTION_NSEPS] = ctrl.numSeqSeps;
    vector<idx_t> part(numSources);
    idx_t sepSize;
    METIS_ComputeVertexSeparator
    ( &nvtxs, xAdj.data(), adjacency.data(), NULL, options, 
      &sepSize, part.data() );
    
    Int sizes[3] = { 0, 0, 0 };
    for( Int s=0; s<numSources; ++s ) 
        ++sizes[part[s]];
    Int offsets[3];
    offsets[0] = 0;
    offsets[1] = sizes[0];
    offsets[2] = sizes[1] + offsets[1];
    perm.resize( numSources ); 
    for( Int s=0; s<numSources; ++s )
        perm[s] = offsets[part[s]]++;
 
    DEBUG_ONLY(EnsurePermutation( perm ))
    BuildChildrenFromPerm
    ( graph, perm, sizes[0], leftChild, sizes[1], rightChild );
    return sizes[2];
#else
    LogicError("METIS interface was not available");
    return -1;
#endif
}

Int Bisect
( const DistGraph& graph, 
        DistGraph& child, 
        DistMap& perm,
        bool& onLeft, 
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Bisect"))
#ifdef EL_HAVE_METIS
    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );
    if( commSize == 1 )
        LogicError
        ("This routine assumes at least two processes are used, "
         "otherwise one child will be lost");

    // (Par)METIS assumes that there are no self-connections or connections 
    // outside the sources, so we must manually remove them from our graph
    const Int numSources = graph.NumSources();
    const Int numLocalEdges = graph.NumLocalEdges();
    Int numLocalValidEdges = 0;
    for( Int i=0; i<numLocalEdges; ++i )
        if( graph.Source(i) != graph.Target(i) && graph.Target(i) < numSources )
            ++numLocalValidEdges;

    // Fill our local connectivity (ignoring self and too-large connections)
    const Int blocksize = graph.Blocksize();
    const Int numLocalSources = graph.NumLocalSources();
    const Int firstLocalSource = graph.FirstLocalSource();
    vector<idx_t> xAdj( numLocalSources+1 );
    vector<idx_t> adjacency( numLocalValidEdges );
    Int validCounter=0;
    Int sourceOff=0;
    Int prevSource=firstLocalSource-1;
    for( Int localEdge=0; localEdge<numLocalEdges; ++localEdge )
    {
        const Int source = graph.Source( localEdge );
        const Int target = graph.Target( localEdge );
        DEBUG_ONLY(
          if( source < prevSource )
              RuntimeError("sources were not properly sorted");
        )
        while( source != prevSource )
        {
            xAdj[sourceOff++] = validCounter;
            ++prevSource;
        }
        if( source != target && target < numSources )
        {
            adjacency[validCounter] = target;
            ++validCounter;
        }
    }
    DEBUG_ONLY(
      if( sourceOff != numLocalSources )
          LogicError("Mistake in xAdj computation");
    )
    xAdj[numLocalSources] = numLocalValidEdges;

    vector<idx_t> sizes(3);
    if( ctrl.sequential )
    {
        // Gather the number of local valid edges on the root process
        vector<Int> edgeSizes( commSize ), edgeOffs;
        mpi::AllGather( &numLocalValidEdges, 1, edgeSizes.data(), 1, comm );
        Int numEdges;
        if( commRank == 0 )
        {
            edgeOffs.resize( commSize );
            numEdges = Scan( edgeSizes, edgeOffs );
        }

        // Gather the edges on the root process (with padding)
        Int maxLocalValidEdges=0;
        for( int q=0; q<commSize; ++q )
            maxLocalValidEdges = Max( maxLocalValidEdges, edgeSizes[q] );
        adjacency.resize( maxLocalValidEdges );
        vector<idx_t> globalAdj;
        if( commRank == 0 )
            globalAdj.resize( maxLocalValidEdges*commSize, 0 );
        mpi::Gather
        ( adjacency.data(), maxLocalValidEdges, 
          globalAdj.data(), maxLocalValidEdges, 0, comm );

        if( commRank == 0 )
        {
            // Remove the padding
            for( int q=1; q<commSize; ++q )
            {
                const Int edgeOff = q*maxLocalValidEdges;
                for( Int j=0; j<edgeSizes[q]; ++j )
                    globalAdj[edgeOffs[q]+j] = globalAdj[edgeOff+j];
            }
            globalAdj.resize( numEdges );
        }

        // Set up the global xAdj vector
        vector<idx_t> globalXAdj;
        // Set up the first commSize*blocksize entries
        if( commRank == 0 )
            globalXAdj.resize( numSources+1 );
        mpi::Gather
        ( xAdj.data(), blocksize, globalXAdj.data(), blocksize, 0, comm );
        if( commRank == 0 )
            for( int q=1; q<commSize; ++q )
                for( Int j=0; j<blocksize; ++j )
                    globalXAdj[q*blocksize+j] += edgeOffs[q];
        // Fix the remaining entries
        const Int numRemaining = numSources - commSize*blocksize;
        if( commRank == commSize-1 )
            mpi::Send( &xAdj[blocksize], numRemaining, 0, comm );
        if( commRank == 0 )
        {
            mpi::Recv
            ( &globalXAdj[commSize*blocksize], numRemaining, commSize-1, comm );
            for( Int j=0; j<numRemaining; ++j )
                globalXAdj[commSize*blocksize+j] += edgeOffs[commSize-1];
            globalXAdj[numSources] = numEdges;
        }

        vector<Int> seqPerm;
        if( commRank == 0 )
        {
            // Call METIS_ComputeVertexSeparator, which is meant for ParMETIS
            idx_t nvtxs = numSources;
            idx_t options[METIS_NOPTIONS];
            METIS_SetDefaultOptions( options );
            options[METIS_OPTION_NSEPS] = ctrl.numSeqSeps;
            vector<idx_t> part(numSources);
            idx_t sepSize;
            METIS_ComputeVertexSeparator
            ( &nvtxs, globalXAdj.data(), globalAdj.data(), NULL, options,
              &sepSize, part.data() );

            for( Int j=0; j<3; ++j )
                sizes[j] = 0;
            for( Int s=0; s<numSources; ++s )
                ++sizes[part[s]];
            Int offsets[3];
            offsets[0] = 0;
            offsets[1] = sizes[0];
            offsets[2] = sizes[1] + offsets[1];
            seqPerm.resize( numSources );
            for( Int s=0; s<numSources; ++s )
                seqPerm[s] = offsets[part[s]]++;
        }

        // Set up space for the distributed permutation
        perm.SetComm( comm );
        perm.Resize( numSources );

        // Distribute the first commSize*blocksize values of the permutation
        mpi::Scatter
        ( seqPerm.data(), blocksize, perm.Buffer(), blocksize, 0, comm );

        // Make sure the last process gets the straggling entries
        if( commRank == 0 )
            mpi::Send
            ( &seqPerm[commSize*blocksize], numRemaining, commSize-1, comm );
        if( commRank == commSize-1 )
            mpi::Recv
            ( perm.Buffer()+blocksize, numLocalSources-blocksize, 0, comm );

        // Broadcast the sizes information from the root
        mpi::Broadcast( (byte*)sizes.data(), 3*sizeof(idx_t), 0, comm );
    }
    else
    {
#ifdef EL_HAVE_PARMETIS
        // Describe the source distribution
        vector<idx_t> vtxDist( commSize+1 );
        for( int i=0; i<commSize; ++i )
            vtxDist[i] = i*blocksize;
        vtxDist[commSize] = graph.NumSources();

        // Create space for the result
        perm.SetComm( comm );
        perm.Resize( numSources );

        vector<idx_t> perm_idx_t( perm.NumLocalSources() );
        // Use the custom ParMETIS interface
        idx_t nparseps = ctrl.numDistSeps;
        real_t imbalance = 1.1;
        ElParallelBisect
        ( vtxDist.data(), xAdj.data(), adjacency.data(), &nparseps, &nseqseps, 
          &imbalance, NULL, perm.Buffer(), sizes.data(), &comm.comm );

        // Since idx_t might be different than Int
        std::copy( perm_idx_t.begin(), perm_idx_t.end(), perm.Buffer() );
#else
        LogicError("ParMETIS was not available");
#endif
    }
    DEBUG_ONLY(EnsurePermutation( perm ))
    BuildChildFromPerm( graph, perm, sizes[0], sizes[1], onLeft, child );
    return sizes[2];
#else
    LogicError("METIS was not available");
    return -1;
#endif
}

void EnsurePermutation( const vector<Int>& map )
{
    DEBUG_ONLY(CallStackEntry cse("EnsurePermutation"))
    const Int numSources = map.size();
    vector<Int> timesMapped( numSources, 0 );
    for( Int i=0; i<numSources; ++i )
        ++timesMapped[map[i]];
    for( Int i=0; i<numSources; ++i )
        if( timesMapped[i] != 1 )
            LogicError
            (timesMapped[i]," vertices were relabeled as ",i,
             " in sequential map");
}

void EnsurePermutation( const DistMap& map )
{
    DEBUG_ONLY(CallStackEntry cse("EnsurePermutation"))
    mpi::Comm comm = map.Comm();
    const int commRank = mpi::Rank( comm );
    const Int numSources = map.NumSources();
    const Int numLocalSources = map.NumLocalSources();
    vector<Int> timesMapped( numSources, 0 );
    for( Int iLocal=0; iLocal<numLocalSources; ++iLocal )
        ++timesMapped[map.GetLocal(iLocal)];
    mpi::Reduce( timesMapped.data(), numSources, MPI_SUM, 0, comm );
    if( commRank == 0 )
        for( Int i=0; i<numSources; ++i )
            if( timesMapped[i] != 1 )
                LogicError
                (timesMapped[i]," vertices were relabeled as ",i,
                 " in parallel map");
}

void BuildChildrenFromPerm
( const Graph& graph, const vector<Int>& perm, 
  Int leftChildSize, Graph& leftChild,
  Int rightChildSize, Graph& rightChild )
{
    DEBUG_ONLY(
      CallStackEntry cse("BuildChildrenFromPerm");
      const Int sepSize = graph.NumSources() - leftChildSize - rightChildSize;
    )
    const Int numSources = graph.NumSources();
    const Int numTargets = graph.NumTargets();

    // Build the inverse permutation
    vector<Int> invPerm( numSources );
    for( Int i=0; i<numSources; ++i )
        invPerm[perm[i]] = i;

    // Get an upper bound on the number of edges in the child graphs
    Int leftChildUpperBound=0, rightChildUpperBound=0;
    for( Int s=0; s<leftChildSize; ++s )
        leftChildUpperBound += graph.NumConnections( invPerm[s] );
    for( Int s=0; s<rightChildSize; ++s )
        rightChildUpperBound += 
            graph.NumConnections( invPerm[s+leftChildSize] );

    // Build the left child's graph
    leftChild.Resize( leftChildSize, numTargets );
    leftChild.Reserve( leftChildUpperBound );
    for( Int s=0; s<leftChildSize; ++s )
    {
        const Int source = s;
        const Int invSource = invPerm[s];
        const Int off = graph.EdgeOffset( invSource );
        const Int numConnections = graph.NumConnections( invSource );
        for( Int t=0; t<numConnections; ++t )
        {
            const Int invTarget = graph.Target( off+t );
            const Int target = ( invTarget < numSources ? 
                                 perm[invTarget] :
                                 invTarget );
            DEBUG_ONLY(
              if( target >= leftChildSize && target < (numSources-sepSize) )
                  LogicError("Invalid bisection, left set touches right set");
            )
            leftChild.QueueConnection( source, target );
        }
    }
    leftChild.MakeConsistent();

    // Build the right child's graph
    rightChild.Resize( rightChildSize, numTargets-leftChildSize );
    rightChild.Reserve( rightChildUpperBound );
    for( Int s=0; s<rightChildSize; ++s )
    {
        const Int source = s+leftChildSize;
        const Int invSource = invPerm[source];
        const Int off = graph.EdgeOffset( invSource );
        const Int numConnections = graph.NumConnections( invSource );
        for( Int t=0; t<numConnections; ++t )
        {
            const Int invTarget = graph.Target( off+t );
            const Int target = ( invTarget < numSources ?
                                 perm[invTarget] :
                                 invTarget );
            DEBUG_ONLY(
              if( target < leftChildSize )
                  LogicError("Invalid bisection, right set touches left set");
            )
            // The targets that are in parent separators do not need to be
            rightChild.QueueConnection
            ( source-leftChildSize, target-leftChildSize );
        }
    }
    rightChild.MakeConsistent();
}

void BuildChildFromPerm
( const DistGraph& graph, const DistMap& perm,
  Int leftChildSize, Int rightChildSize,
  bool& onLeft, DistGraph& child )
{
    DEBUG_ONLY(CallStackEntry cse("BuildChildFromPerm"))
    const Int numTargets = graph.NumTargets();
    const Int numLocalSources = graph.NumLocalSources();
    DEBUG_ONLY(
      const Int numSources = graph.NumSources();
      const Int sepSize = numSources - leftChildSize - rightChildSize;
    )

    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

    // Build the child graph from the partitioned parent
    const int smallTeamSize = commSize/2;
    const int largeTeamSize = commSize - smallTeamSize;
    const bool inSmallTeam = ( commRank < smallTeamSize );
    const bool smallOnLeft = ( leftChildSize <= rightChildSize );
    const int leftTeamSize = ( smallOnLeft ? smallTeamSize : largeTeamSize );
    const int rightTeamSize = ( smallOnLeft ? largeTeamSize : smallTeamSize );
    const int leftTeamOff = ( smallOnLeft ? 0 : smallTeamSize );
    const int rightTeamOff = ( smallOnLeft ? smallTeamSize : 0 );
    onLeft = ( inSmallTeam == smallOnLeft );

    const Int leftTeamBlocksize = leftChildSize / leftTeamSize;
    const Int rightTeamBlocksize = rightChildSize / rightTeamSize;

    // Count how many rows we must send to each process 
    vector<int> rowSendSizes( commSize, 0 );
    for( Int s=0; s<numLocalSources; ++s )
    {
        const Int i = perm.GetLocal(s);
        if( i < leftChildSize )
        {
            const int q = leftTeamOff + 
                RowToProcess( i, leftTeamBlocksize, leftTeamSize );
            ++rowSendSizes[q];
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q = rightTeamOff +
                RowToProcess
                ( i-leftChildSize, rightTeamBlocksize, rightTeamSize );
            ++rowSendSizes[q];
        }
    }

    // Exchange the number of rows
    vector<int> rowRecvSizes( commSize );
    mpi::AllToAll( rowSendSizes.data(), 1, rowRecvSizes.data(), 1, comm );

    // Prepare for the AllToAll to exchange the row indices and 
    // the number of column indices per row
    vector<int> rowSendOffs, rowRecvOffs;
    const int numSendRows = Scan( rowSendSizes, rowSendOffs );
    const int numRecvRows = Scan( rowRecvSizes, rowRecvOffs );

    // Pack the row indices and how many column entries there will be per row
    vector<Int> rowSendLengths( numSendRows );
    vector<Int> rowSendInds( numSendRows );
    vector<int> offs = rowSendOffs;
    for( Int s=0; s<numLocalSources; ++s )
    {
        const Int i = perm.GetLocal(s);
        if( i < leftChildSize )
        {
            const int q = leftTeamOff + 
                RowToProcess( i, leftTeamBlocksize, leftTeamSize );
            rowSendInds[offs[q]] = i;
            rowSendLengths[offs[q]] = graph.NumConnections( s );
            ++offs[q];
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q = rightTeamOff + 
                RowToProcess
                ( i-leftChildSize, rightTeamBlocksize, rightTeamSize );
            rowSendInds[offs[q]] = i;
            rowSendLengths[offs[q]] = graph.NumConnections( s );
            ++offs[q];
        }
    }

    // Perform the row lengths exchange
    vector<Int> rowRecvLengths( numRecvRows );
    mpi::AllToAll
    ( rowSendLengths.data(), rowSendSizes.data(), rowSendOffs.data(),
      rowRecvLengths.data(), rowRecvSizes.data(), rowRecvOffs.data(), comm );

    // Perform the row indices exchange
    vector<Int> rowRecvInds( numRecvRows );
    mpi::AllToAll
    ( rowSendInds.data(), rowSendSizes.data(), rowSendOffs.data(),
      rowRecvInds.data(), rowRecvSizes.data(), rowRecvOffs.data(), comm );
    SwapClear( rowSendInds );

    // Set up for sending the column indices
    int numSendInds=0;
    vector<int> indSendSizes( commSize, 0 ), indSendOffs(commSize);
    for( int q=0; q<commSize; ++q )
    {
        const int numRows = rowSendSizes[q];
        const int off = rowSendOffs[q];
        for( int s=0; s<numRows; ++s )
            indSendSizes[q] += rowSendLengths[off+s];

        indSendOffs[q] = numSendInds;
        numSendInds += indSendSizes[q];
    }
    SwapClear( rowSendLengths );
    int numRecvInds=0;
    vector<int> indRecvSizes( commSize, 0 ), indRecvOffs(commSize);
    for( int q=0; q<commSize; ++q )
    {
        const int numRows = rowRecvSizes[q];
        const int off = rowRecvOffs[q];
        for( int s=0; s<numRows; ++s )
            indRecvSizes[q] += rowRecvLengths[off+s];

        indRecvOffs[q] = numRecvInds;
        numRecvInds += indRecvSizes[q];
    }
    SwapClear( rowSendSizes );
    SwapClear( rowSendOffs );

    // Pack the indices
    vector<Int> sendInds( numSendInds );
    offs = indSendOffs;
    for( Int s=0; s<numLocalSources; ++s )
    {
        const Int i = perm.GetLocal(s);
        if( i < leftChildSize )
        {
            const int q = leftTeamOff + 
                RowToProcess( i, leftTeamBlocksize, leftTeamSize );

            int& off = offs[q];
            const Int numConnections = graph.NumConnections( s );
            const Int localEdgeOff = graph.EdgeOffset( s );
            for( Int j=0; j<numConnections; ++j )
                sendInds[off++] = graph.Target( localEdgeOff+j );
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q = rightTeamOff + 
                RowToProcess
                ( i-leftChildSize, rightTeamBlocksize, rightTeamSize );
               
            int& off = offs[q];
            const Int numConnections = graph.NumConnections( s );
            const Int localEdgeOff = graph.EdgeOffset( s );
            for( Int j=0; j<numConnections; ++j )
                sendInds[off++] = graph.Target( localEdgeOff+j );
        }
    }

    // Send/recv the column indices
    vector<Int> recvInds( numRecvInds );
    mpi::AllToAll
    ( sendInds.data(), indSendSizes.data(), indSendOffs.data(),
      recvInds.data(), indRecvSizes.data(), indRecvOffs.data(), comm );
    SwapClear( sendInds );
    SwapClear( indSendSizes );
    SwapClear( indSendOffs );

    // Get the indices after reordering
    perm.Translate( recvInds );

    // Put the connections into our new graph
    const int childTeamRank = 
        ( onLeft ? commRank-leftTeamOff : commRank-rightTeamOff );
    mpi::Comm childComm;
    mpi::Split( comm, onLeft, childTeamRank, childComm );
    child.SetComm( childComm );
    mpi::Free( childComm );
    if( onLeft )
        child.Resize( leftChildSize, numTargets );
    else
        child.Resize( rightChildSize, numTargets-leftChildSize );

    child.Reserve( recvInds.size() );
    Int off=0;
    for( Int s=0; s<numRecvRows; ++s )
    {
        const Int source = rowRecvInds[s];
        const Int numConnections = rowRecvLengths[s];
        const Int childFirstLocalSource = child.FirstLocalSource();
        for( Int t=0; t<numConnections; ++t )
        {
            const Int target = recvInds[off++];
            if( onLeft )
            {
                DEBUG_ONLY(
                  if( target >= leftChildSize && 
                      target < (numSources-sepSize) )
                      LogicError
                      ("Invalid bisection, left set touches right:\n",
                       "  ",source," touches ",target," and leftChildSize=",
                       leftChildSize);
                )
                child.QueueLocalConnection
                ( source-childFirstLocalSource, target );
            }
            else
            {
                DEBUG_ONLY(
                  if( target < leftChildSize )
                      LogicError
                      ("Invalid bisection, right set touches left set");
                )
                child.QueueLocalConnection
                ( source-leftChildSize-childFirstLocalSource, 
                  target-leftChildSize );
            }
        }
    }
    child.MakeConsistent();
}

void BuildMap
( const Separator& rootSep, 
        vector<Int>& map )
{
    DEBUG_ONLY(CallStackEntry cse("BuildMap"))
    const Int numSources = rootSep.off + rootSep.inds.size();
    map.resize( numSources );

    function<void(const Separator&)> buildMap = 
      [&]( const Separator& sep )
      {
        for( auto* child : sep.children )  
            buildMap( *child );
        for( Int t=0; t<sep.inds.size(); ++t )
            map[sep.inds[t]] = sep.off + t;
      };
    buildMap( rootSep );
}

void BuildMap
( const DistSeparator& rootSep, 
        DistMap& map )
{
    DEBUG_ONLY(CallStackEntry cse("BuildMap"))

    const Int numSources = rootSep.off + rootSep.inds.size();
    mpi::Comm comm = rootSep.comm;
    map.SetComm( comm );
    map.Resize( numSources );

    const int commSize = mpi::Size( comm );
    vector<int> sendSizes( commSize, 0 );
    function<void(const Separator&)> sendSizeLocalAccumulate =
      [&]( const Separator& sep )
      {
        for( const Separator* childSep : sep.children )
            sendSizeLocalAccumulate( *childSep );
        for( Int i : sep.inds )
            ++sendSizes[ map.RowOwner(i) ];
      };
    function<void(const DistSeparator&)> sendSizeAccumulate = 
      [&]( const DistSeparator& sep )
      {
          if( sep.child == nullptr )
          {
              const Separator& dup = *sep.duplicate;
              for( const Separator* childSep : dup.children )
                  sendSizeLocalAccumulate( *childSep );
          }
          else
              sendSizeAccumulate( *sep.child );

          const Int numInds = sep.inds.size();
          const int teamSize = mpi::Size( sep.comm );
          const int teamRank = mpi::Rank( sep.comm );
          const Int numLocalInds = Length( numInds, teamRank, teamSize );
          for( Int tLocal=0; tLocal<numLocalInds; ++tLocal )
          {
              const Int t = teamRank + tLocal*teamSize;
              ++sendSizes[ map.RowOwner(sep.inds[t]) ];
          }
      };
    sendSizeAccumulate( rootSep );

    // Use a single-entry AllToAll to coordinate how many indices will be 
    // exchanges
    vector<int> recvSizes( commSize );
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, comm );

    // Pack the reordered indices
    vector<int> sendOffs;
    const int numSends = Scan( sendSizes, sendOffs );
    vector<Int> sendInds(numSends), sendOrigInds(numSends);
    auto offs = sendOffs;
    function<void(const Separator&)> packRowsLocal =
      [&]( const Separator& sep )
      {
          for( const Separator* childSep : sep.children )    
              packRowsLocal( *childSep );

          const Int numInds = sep.inds.size(); 
          for( Int t=0; t<numInds; ++t )
          {
              const Int i = sep.inds[t];
              const Int iMap = sep.off + t;              
              const int q = map.RowOwner(i);
              sendOrigInds[offs[q]] = i;
              sendInds[offs[q]] = iMap;
              ++offs[q];
          }
      };
    function<void(const DistSeparator&)> packRows = 
      [&]( const DistSeparator& sep )
      {
          if( sep.child == nullptr )
          {
              const Separator& dup = *sep.duplicate;
              for( const Separator* childSep : dup.children )
                  packRowsLocal( *childSep );
          }
          else
              packRows( *sep.child );

          const Int numInds = sep.inds.size();
          const int teamSize = mpi::Size( sep.comm );
          const int teamRank = mpi::Rank( sep.comm );
          const Int numLocalInds = Length( numInds, teamRank, teamSize );
          for( Int tLocal=0; tLocal<numLocalInds; ++tLocal )
          {
              const Int t = teamRank + tLocal*teamSize;
              const Int i = sep.inds[t];
              const Int iMap = sep.off + t;
              const int q = map.RowOwner(i);
              sendOrigInds[offs[q]] = i;
              sendInds[offs[q]] = iMap;
              ++offs[q];
          }
      };
    packRows( rootSep );

    // Perform an AllToAll to exchange the reordered indices
    vector<int> recvOffs;
    const int numRecvs = Scan( recvSizes, recvOffs );
    DEBUG_ONLY(
      const Int numLocalSources = map.NumLocalSources();
      if( numRecvs != numLocalSources )
          LogicError("incorrect number of recv indices");
    )
    vector<Int> recvInds( numRecvs );
    mpi::AllToAll
    ( sendInds.data(), sendSizes.data(), sendOffs.data(),
      recvInds.data(), recvSizes.data(), recvOffs.data(), comm );

    // Perform an AllToAll to exchange the original indices
    vector<Int> recvOrigInds( numRecvs );
    mpi::AllToAll
    ( sendOrigInds.data(), sendSizes.data(), sendOffs.data(),
      recvOrigInds.data(), recvSizes.data(), recvOffs.data(), comm );

    // Unpack the indices
    const Int firstLocalSource = map.FirstLocalSource();
    for( Int s=0; s<numRecvs; ++s )
        map.SetLocal( recvOrigInds[s]-firstLocalSource, recvInds[s] );
}

} // namespace El
