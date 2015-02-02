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
DistributedDepthRecursion
( unsigned commRank, unsigned commSize, unsigned& distDepth )
{
    if( commSize == 1 )
        return;

    ++distDepth;
    const unsigned smallTeamSize = commSize/2;
    const unsigned largeTeamSize = commSize - smallTeamSize;
    if( commRank < smallTeamSize )
        DistributedDepthRecursion( commRank, smallTeamSize, distDepth );
    else
        DistributedDepthRecursion
        ( commRank-smallTeamSize, largeTeamSize, distDepth );
}

int DistributedDepth( mpi::Comm comm )
{
    unsigned commRank = mpi::Rank( comm );
    unsigned commSize = mpi::Size( comm );
    unsigned distDepth = 0;
    DistributedDepthRecursion( commRank, commSize, distDepth );
    return distDepth;
}

inline void
NestedDissectionRecursion
( const Graph& graph, 
  const std::vector<Int>& perm,
        DistSeparatorTree& sepTree, 
        DistSymmElimTree& eTree,
        Int parent, 
        Int off, 
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NestedDissectionRecursion"))
    if( graph.NumSources() <= ctrl.cutoff )
    {
        // Fill in this node of the local separator tree
        const Int numSources = graph.NumSources();
        sepTree.localSepsAndLeaves.push_back( new SepOrLeaf );
        SepOrLeaf& leaf = *sepTree.localSepsAndLeaves.back();
        leaf.parent = parent;
        leaf.off = off;
        leaf.inds = perm;

        // Fill in this node of the local elimination tree
        eTree.localNodes.push_back( new SymmNode );
        SymmNode& node = *eTree.localNodes.back();
        node.size = numSources;
        node.off = off;
        node.parent = parent;
        SwapClear( node.children );
        std::set<Int> connectedAncestors;
        for( Int s=0; s<node.size; ++s )
        {
            const Int numConnections = graph.NumConnections( s );
            const Int edgeOff = graph.EdgeOffset( s );
            for( Int t=0; t<numConnections; ++t )
            {
                const Int target = graph.Target( edgeOff+t );
                if( target >= numSources )
                    connectedAncestors.insert( off+target );
            }
        }
        node.lowerStruct.resize( connectedAncestors.size() );
        std::copy
        ( connectedAncestors.begin(), connectedAncestors.end(), 
          node.lowerStruct.begin() );
    }
    else
    {
        // Partition the graph and construct the inverse map
        Graph leftChild, rightChild;
        std::vector<Int> map;
        const Int sepSize = Bisect( graph, leftChild, rightChild, map, ctrl );
        const Int numSources = graph.NumSources();
        std::vector<Int> inverseMap( numSources );
        for( Int s=0; s<numSources; ++s )
            inverseMap[map[s]] = s;

        // Mostly compute this node of the local separator tree
        // (we will finish computing the separator indices soon)
        sepTree.localSepsAndLeaves.push_back( new SepOrLeaf );
        SepOrLeaf& sep = *sepTree.localSepsAndLeaves.back();
        sep.parent = parent;
        sep.off = off + (numSources-sepSize);
        sep.inds.resize( sepSize );
        for( Int s=0; s<sepSize; ++s )
        {
            const Int mappedSource = s + (numSources-sepSize);
            sep.inds[s] = inverseMap[mappedSource];
        }
    
        // Fill in this node in the local elimination tree
        eTree.localNodes.push_back( new SymmNode );
        SymmNode& node = *eTree.localNodes.back();
        node.size = sepSize;
        node.off = sep.off;
        node.parent = parent;
        node.children.resize( 2 );
        std::set<Int> connectedAncestors;
        for( Int s=0; s<sepSize; ++s )
        {
            const Int source = sep.inds[s];
            const Int numConnections = graph.NumConnections( source );
            const Int edgeOff = graph.EdgeOffset( source );
            for( Int t=0; t<numConnections; ++t )
            {
                const Int target = graph.Target( edgeOff+t );
                if( target >= numSources )
                    connectedAncestors.insert( off+target );
            }
        }
        node.lowerStruct.resize( connectedAncestors.size() );
        std::copy
        ( connectedAncestors.begin(), connectedAncestors.end(), 
          node.lowerStruct.begin() );

        // Finish computing the separator indices
        for( Int s=0; s<sepSize; ++s )
            sep.inds[s] = perm[sep.inds[s]];

        // Construct the inverse maps from the child indices to the original
        // degrees of freedom
        const Int leftChildSize = leftChild.NumSources();
        std::vector<Int> leftPerm( leftChildSize );
        for( Int s=0; s<leftChildSize; ++s )
            leftPerm[s] = perm[inverseMap[s]];
        const Int rightChildSize = rightChild.NumSources();
        std::vector<Int> rightPerm( rightChildSize );
        for( Int s=0; s<rightChildSize; ++s )
            rightPerm[s] = perm[inverseMap[s+leftChildSize]];

        // Update right then left so that, once we later reverse the order 
        // of the nodes, the left node will be ordered first
        const Int parent = eTree.localNodes.size()-1;
        node.children[1] = eTree.localNodes.size();
        NestedDissectionRecursion
        ( rightChild, rightPerm, sepTree, eTree, parent, off+leftChildSize,
          ctrl );
        node.children[0] = eTree.localNodes.size();
        NestedDissectionRecursion
        ( leftChild, leftPerm, sepTree, eTree, parent, off, ctrl );
    }
}

inline void
NestedDissectionRecursion
( const DistGraph& graph, 
  const DistMap& perm,
        DistSeparatorTree& sepTree, 
        DistSymmElimTree& eTree,
        Int depth, 
        Int off, 
        bool onLeft,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NestedDissectionRecursion"))
    const Int distDepth = sepTree.distSeps.size();
    mpi::Comm comm = graph.Comm();
    if( distDepth - depth > 0 )
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

        DistMap inverseMap;
        map.FormInverse( inverseMap );

        // Mostly fill this node of the DistSeparatorTree
        // (we will finish computing the separator indices at the end)
        DistSeparator& sep = sepTree.distSeps[distDepth-1-depth];
        mpi::Dup( comm, sep.comm );
        sep.off = off + (numSources-sepSize);
        sep.inds.resize( sepSize );
        for( Int s=0; s<sepSize; ++s )
            sep.inds[s] = s + (numSources-sepSize);
        inverseMap.Translate( sep.inds );

        // Fill in this node of the DistSymmElimTree
        DistSymmNode& node = eTree.distNodes[distDepth-depth];
        node.size = sepSize;
        node.off = sep.off;
        node.onLeft = onLeft;
        mpi::Dup( comm, node.comm );
        const Int numLocalSources = graph.NumLocalSources();
        const Int firstLocalSource = graph.FirstLocalSource();
        std::set<Int> localConnectedAncestors;
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
                        localConnectedAncestors.insert( off+target );
                }
            }
        }
        const int numLocalConnected = localConnectedAncestors.size();
        const int commSize = mpi::Size( comm );
        std::vector<int> localConnectedSizes( commSize );
        mpi::AllGather
        ( &numLocalConnected, 1, localConnectedSizes.data(), 1, comm );
        std::vector<Int> localConnectedVec( numLocalConnected );
        std::copy
        ( localConnectedAncestors.begin(), localConnectedAncestors.end(), 
          localConnectedVec.begin() );
        std::vector<int> localConnectedOffs;
        const int sumOfLocalConnectedSizes = 
            Scan( localConnectedSizes, localConnectedOffs );
        std::vector<Int> localConnections( sumOfLocalConnectedSizes );
        mpi::AllGather
        ( localConnectedVec.data(), numLocalConnected,
          localConnections.data(), 
          localConnectedSizes.data(), localConnectedOffs.data(), comm );
        std::set<Int> connectedAncestors
        ( localConnections.begin(), localConnections.end() );
        node.lowerStruct.resize( connectedAncestors.size() );
        std::copy
        ( connectedAncestors.begin(), connectedAncestors.end(), 
          node.lowerStruct.begin() );

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
        inverseMap.Extend( newPerm );
        perm.Extend( newPerm );

        // Recurse
        const Int newOff = ( childIsOnLeft ? off : off+leftChildSize );
        NestedDissectionRecursion
        ( child, newPerm, sepTree, eTree, depth+1, newOff, childIsOnLeft, 
          ctrl );
    }
    else if( graph.NumSources() <= ctrl.cutoff )
    {
        // Convert to a sequential graph
        const Int numSources = graph.NumSources();
        Graph seqGraph( graph );

        // Fill in this node of the local separator tree
        sepTree.localSepsAndLeaves.push_back( new SepOrLeaf );
        SepOrLeaf& leaf = *sepTree.localSepsAndLeaves.back();
        leaf.parent = -1;
        leaf.off = off;
        leaf.inds = perm.Map();

        // Fill in this node of the local and distributed parts of the 
        // elimination tree
        eTree.localNodes.push_back( new SymmNode );
        SymmNode& localNode = *eTree.localNodes.back();
        DistSymmNode& distNode = eTree.distNodes[0];
        mpi::Dup( comm, distNode.comm );
        distNode.onLeft = onLeft;
        distNode.size = localNode.size = numSources;
        distNode.off = localNode.off = off;
        localNode.parent = -1;
        SwapClear( localNode.children );
        std::set<Int> connectedAncestors;
        for( Int s=0; s<numSources; ++s )
        {
            const Int numConnections = seqGraph.NumConnections( s );
            const Int edgeOff = seqGraph.EdgeOffset( s );
            for( Int t=0; t<numConnections; ++t )
            {
                const Int target = seqGraph.Target( edgeOff+t );
                if( target >= numSources )
                    connectedAncestors.insert( off+target );
            }
        }
        localNode.lowerStruct.resize( connectedAncestors.size() );
        std::copy
        ( connectedAncestors.begin(), connectedAncestors.end(), 
          localNode.lowerStruct.begin() );    
        distNode.lowerStruct = localNode.lowerStruct;
    }
    else
    {
        // Convert to a sequential graph
        Graph seqGraph( graph );

        // Partition the graph and construct the inverse map
        Graph leftChild, rightChild;
        std::vector<Int> map;
        const Int sepSize = 
            Bisect( seqGraph, leftChild, rightChild, map, ctrl );
        const Int numSources = graph.NumSources();
        std::vector<Int> inverseMap( numSources );
        for( Int s=0; s<numSources; ++s )
            inverseMap[map[s]] = s;

        // Mostly compute this node of the local separator tree
        // (we will finish computing the separator indices soon)
        sepTree.localSepsAndLeaves.push_back( new SepOrLeaf );
        SepOrLeaf& sep = *sepTree.localSepsAndLeaves.back();
        sep.parent = -1;
        sep.off = off + (numSources-sepSize);
        sep.inds.resize( sepSize );
        for( Int s=0; s<sepSize; ++s )
        {
            const Int mappedSource = s + (numSources-sepSize);
            sep.inds[s] = inverseMap[mappedSource];
        }
        
        // Fill in this node in both the local and distributed parts of 
        // the elimination tree
        eTree.localNodes.push_back( new SymmNode );
        SymmNode& localNode = *eTree.localNodes.back();
        DistSymmNode& distNode = eTree.distNodes[0];
        mpi::Dup( comm, distNode.comm );
        distNode.onLeft = onLeft;
        distNode.size = localNode.size = sepSize;
        distNode.off = localNode.off = sep.off;
        localNode.parent = -1;
        localNode.children.resize( 2 );
        std::set<Int> connectedAncestors;
        for( Int s=0; s<sepSize; ++s )
        {
            const Int source = sep.inds[s];
            const Int numConnections = seqGraph.NumConnections( source );
            const Int edgeOff = seqGraph.EdgeOffset( source );
            for( Int t=0; t<numConnections; ++t )
            {
                const Int target = seqGraph.Target( edgeOff+t );
                if( target >= numSources )
                    connectedAncestors.insert( off+target );
            }
        }
        localNode.lowerStruct.resize( connectedAncestors.size() );
        std::copy
        ( connectedAncestors.begin(), connectedAncestors.end(), 
          localNode.lowerStruct.begin() );
        distNode.lowerStruct = localNode.lowerStruct;

        // Finish computing the separator indices
        // (This is a faster version of the Translate member function)
        for( Int s=0; s<sepSize; ++s )
            sep.inds[s] = perm.GetLocal( sep.inds[s] );

        // Construct the inverse maps from the child indices to the original
        // degrees of freedom
        const Int leftChildSize = leftChild.NumSources();
        std::vector<Int> leftPerm( leftChildSize );
        for( Int s=0; s<leftChildSize; ++s )
            leftPerm[s] = perm.GetLocal( inverseMap[s] );
        const Int rightChildSize = rightChild.NumSources();
        std::vector<Int> rightPerm( rightChildSize );
        for( Int s=0; s<rightChildSize; ++s )
            rightPerm[s] = perm.GetLocal( inverseMap[s+leftChildSize] );

        // Update right then left so that, once we later reverse the order 
        // of the nodes, the left node will be ordered first
        const Int parent=0;
        localNode.children[1] = eTree.localNodes.size();
        NestedDissectionRecursion
        ( rightChild, rightPerm, sepTree, eTree, parent, off+leftChildSize, 
          ctrl );
        localNode.children[0] = eTree.localNodes.size();
        NestedDissectionRecursion
        ( leftChild, leftPerm, sepTree, eTree, parent, off, ctrl );
    }
}

void NestedDissection
( const DistGraph& graph, 
        DistMap& map,
        DistSeparatorTree& sepTree, 
        DistSymmInfo& info,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NestedDissection"))
    // NOTE: There is a potential memory leak here if these data structures 
    //       are reused. Their destructors should call a member function which
    //       we can simply call here to clear the data
    DistSymmElimTree eTree;
    SwapClear( eTree.localNodes );
    SwapClear( sepTree.localSepsAndLeaves );

    mpi::Comm comm = graph.Comm();
    const int distDepth = DistributedDepth( comm );
    eTree.distNodes.resize( distDepth+1 );
    sepTree.distSeps.resize( distDepth );

    DistMap perm( graph.NumSources(), graph.Comm() );
    const Int firstLocalSource = perm.FirstLocalSource();
    const Int numLocalSources = perm.NumLocalSources();
    for( Int s=0; s<numLocalSources; ++s )
        perm.SetLocal( s, s+firstLocalSource );
    NestedDissectionRecursion( graph, perm, sepTree, eTree, 0, 0, false, ctrl );

    ReverseOrder( sepTree, eTree );

    // Construct the distributed reordering    
    BuildMap( graph, sepTree, map );
    DEBUG_ONLY(EnsurePermutation( map ))

    // Run the symbolic analysis
    SymmetricAnalysis( eTree, info, ctrl.storeFactRecvInds );
}

Int Bisect
( const Graph& graph ,Graph& leftChild, Graph& rightChild,
  std::vector<Int>& perm, const BisectCtrl& ctrl )
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
    std::vector<idx_t> xAdj( numSources+1 );
    std::vector<idx_t> adjacency( numValidEdges );
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
    std::vector<idx_t> part(numSources);
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
    std::vector<idx_t> xAdj( numLocalSources+1 );
    std::vector<idx_t> adjacency( numLocalValidEdges );
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

    std::vector<idx_t> sizes(3);
    if( ctrl.sequential )
    {
        // Gather the number of local valid edges on the root process
        std::vector<Int> edgeSizes( commSize ), edgeOffs;
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
            maxLocalValidEdges = 
                std::max( maxLocalValidEdges, edgeSizes[q] );
        adjacency.resize( maxLocalValidEdges );
        std::vector<idx_t> globalAdj;
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
        std::vector<idx_t> globalXAdj;
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

        std::vector<Int> seqPerm;
        if( commRank == 0 )
        {
            // Call METIS_ComputeVertexSeparator, which is meant for ParMETIS
            idx_t nvtxs = numSources;
            idx_t options[METIS_NOPTIONS];
            METIS_SetDefaultOptions( options );
            options[METIS_OPTION_NSEPS] = ctrl.numSeqSeps;
            std::vector<idx_t> part(numSources);
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
        std::vector<idx_t> vtxDist( commSize+1 );
        for( int i=0; i<commSize; ++i )
            vtxDist[i] = i*blocksize;
        vtxDist[commSize] = graph.NumSources();

        // Create space for the result
        perm.SetComm( comm );
        perm.Resize( numSources );

        std::vector<idx_t> perm_idx_t( perm.NumLocalSources() );
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

void EnsurePermutation( const std::vector<Int>& map )
{
    DEBUG_ONLY(CallStackEntry cse("EnsurePermutation"))
    const Int numSources = map.size();
    std::vector<Int> timesMapped( numSources, 0 );
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
    std::vector<Int> timesMapped( numSources, 0 );
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

void ReverseOrder( DistSeparatorTree& sepTree, DistSymmElimTree& eTree )
{
    // Reverse the order of the pointers and indices in the elimination and 
    // separator trees (so that the leaves come first)
    const Int numLocalNodes = eTree.localNodes.size();
    const Int lastInd = numLocalNodes-1;
    if( numLocalNodes != 1 )
    {
        // Switch the pointers for the root and last nodes
        SymmNode* rootNode = eTree.localNodes[0];
        SymmNode* lastNode = eTree.localNodes.back();
        eTree.localNodes[0] = lastNode;
        eTree.localNodes.back() = rootNode;
        SepOrLeaf* rootSep = sepTree.localSepsAndLeaves[0];
        SepOrLeaf* lastLeaf = sepTree.localSepsAndLeaves.back();
        sepTree.localSepsAndLeaves[0] = lastLeaf;
        sepTree.localSepsAndLeaves.back() = rootSep;

        // Update their parent indices to what their final values will be
        // (The root node's parent index does not need to be changed.)
        lastNode->parent = lastInd - lastNode->parent;
        lastLeaf->parent = lastInd - lastLeaf->parent;
        // Update their children's indices
        // (The last node will not have children)
        const Int numRootChildren = rootNode->children.size();
        for( Int c=0; c<numRootChildren; ++c )
            rootNode->children[c] = lastInd - rootNode->children[c];
    }
    // Switch the middle nodes (we will miss the middle node if an odd number)
    for( Int s=1; s<numLocalNodes/2; ++s )
    {
        const Int t = lastInd - s;
        // Switch the pointers for the last and right nodes
        SymmNode* leftNode = eTree.localNodes[s];
        SymmNode* rightNode = eTree.localNodes[t];
        SepOrLeaf* leftSepOrLeaf = sepTree.localSepsAndLeaves[s];
        SepOrLeaf* rightSepOrLeaf = sepTree.localSepsAndLeaves[t];
        eTree.localNodes[s] = rightNode;
        eTree.localNodes[t] = leftNode;
        sepTree.localSepsAndLeaves[s] = rightSepOrLeaf;
        sepTree.localSepsAndLeaves[t] = leftSepOrLeaf;
        // Update their parent indices to what their final values will be
        leftNode->parent = lastInd - leftNode->parent;
        rightNode->parent = lastInd - rightNode->parent;
        leftSepOrLeaf->parent = lastInd - leftSepOrLeaf->parent;
        rightSepOrLeaf->parent = lastInd - rightSepOrLeaf->parent;
        // Update their children's indices
        const Int numLeftChildren = leftNode->children.size();
        for( Int c=0; c<numLeftChildren; ++c )
            leftNode->children[c] = lastInd - leftNode->children[c];
        const Int numRightChildren = rightNode->children.size();
        for( Int c=0; c<numRightChildren; ++c )
            rightNode->children[c] = lastInd - rightNode->children[c];
    }
    // Handle the middle node if it exists
    if( numLocalNodes % 2 != 0 )
    {
        const Int midInd = numLocalNodes/2;
        // Update the parent indices to the final values
        SymmNode* middleNode = eTree.localNodes[midInd];
        SepOrLeaf* middleSepOrLeaf = sepTree.localSepsAndLeaves[midInd];
        middleNode->parent = lastInd - middleNode->parent;
        middleSepOrLeaf->parent = lastInd - middleSepOrLeaf->parent;
        // Update the children's indices
        const Int numChildren = middleNode->children.size();
        for( Int c=0; c<numChildren; ++c )
            middleNode->children[c] = lastInd - middleNode->children[c];
    }
}

void BuildChildrenFromPerm
( const Graph& graph, const std::vector<Int>& perm, 
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
    std::vector<Int> inversePerm( numSources );
    for( Int i=0; i<numSources; ++i )
        inversePerm[perm[i]] = i;

    // Get an upper bound on the number of edges in the child graphs
    Int leftChildUpperBound=0, rightChildUpperBound=0;
    for( Int s=0; s<leftChildSize; ++s )
        leftChildUpperBound += graph.NumConnections( inversePerm[s] );
    for( Int s=0; s<rightChildSize; ++s )
        rightChildUpperBound += 
            graph.NumConnections( inversePerm[s+leftChildSize] );

    // Build the left child's graph
    leftChild.Resize( leftChildSize, numTargets );
    leftChild.Reserve( leftChildUpperBound );
    for( Int s=0; s<leftChildSize; ++s )
    {
        const Int source = s;
        const Int inverseSource = inversePerm[s];
        const Int off = graph.EdgeOffset( inverseSource );
        const Int numConnections = graph.NumConnections( inverseSource );
        for( Int t=0; t<numConnections; ++t )
        {
            const Int inverseTarget = graph.Target( off+t );
            const Int target = ( inverseTarget < numSources ? 
                                 perm[inverseTarget] :
                                 inverseTarget );
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
        const Int inverseSource = inversePerm[source];
        const Int off = graph.EdgeOffset( inverseSource );
        const Int numConnections = graph.NumConnections( inverseSource );
        for( Int t=0; t<numConnections; ++t )
        {
            const Int inverseTarget = graph.Target( off+t );
            const Int target = ( inverseTarget < numSources ?
                                 perm[inverseTarget] :
                                 inverseTarget );
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
    std::vector<int> rowSendSizes( commSize, 0 );
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
    std::vector<int> rowRecvSizes( commSize );
    mpi::AllToAll( rowSendSizes.data(), 1, rowRecvSizes.data(), 1, comm );

    // Prepare for the AllToAll to exchange the row indices and 
    // the number of column indices per row
    std::vector<int> rowSendOffs, rowRecvOffs;
    const int numSendRows = Scan( rowSendSizes, rowSendOffs );
    const int numRecvRows = Scan( rowRecvSizes, rowRecvOffs );

    // Pack the row indices and how many column entries there will be per row
    std::vector<Int> rowSendLengths( numSendRows );
    std::vector<Int> rowSendInds( numSendRows );
    std::vector<int> offs = rowSendOffs;
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
    std::vector<Int> rowRecvLengths( numRecvRows );
    mpi::AllToAll
    ( rowSendLengths.data(), rowSendSizes.data(), rowSendOffs.data(),
      rowRecvLengths.data(), rowRecvSizes.data(), rowRecvOffs.data(), comm );

    // Perform the row indices exchange
    std::vector<Int> rowRecvInds( numRecvRows );
    mpi::AllToAll
    ( rowSendInds.data(), rowSendSizes.data(), rowSendOffs.data(),
      rowRecvInds.data(), rowRecvSizes.data(), rowRecvOffs.data(), comm );
    SwapClear( rowSendInds );

    // Set up for sending the column indices
    int numSendInds=0;
    std::vector<int> indSendSizes( commSize, 0 ), indSendOffs(commSize);
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
    std::vector<int> indRecvSizes( commSize, 0 ), indRecvOffs(commSize);
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
    std::vector<Int> sendInds( numSendInds );
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
    std::vector<Int> recvInds( numRecvInds );
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

inline void
BuildMap
( const DistGraph& graph, 
  const DistSeparatorTree& sepTree, 
        DistMap& map )
{
    DEBUG_ONLY(CallStackEntry cse("BuildMap"))
    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::Size( comm );
    const Int numSources = graph.NumSources();

    map.SetComm( comm );
    map.Resize( numSources );

    const Int numLocal = sepTree.localSepsAndLeaves.size();
    // NOTE: The dist separator tree does not double-count the first 
    //       single-process node, but DistSymmInfo does. Thus their number of
    //       distributed nodes is different by one.
    const Int numDist = sepTree.distSeps.size();

    // Traverse local sepTree to count how many indices we should send the
    // final index for
    std::vector<int> sendSizes( commSize, 0 );
    for( int s=0; s<numLocal; ++s )
    {
        const SepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const Int numInds = sepOrLeaf.inds.size();
        for( Int t=0; t<numInds; ++t )
            ++sendSizes[ map.RowOwner(sepOrLeaf.inds[t]) ];
    }
    for( Int s=0; s<numDist; ++s )
    {
        const DistSeparator& sep = sepTree.distSeps[s];
        const Int numInds = sep.inds.size();
        const int teamSize = mpi::Size( sep.comm );
        const int teamRank = mpi::Rank( sep.comm );
        const Int numLocalInds = Length( numInds, teamRank, teamSize );
        for( Int tLocal=0; tLocal<numLocalInds; ++tLocal )
        {
            const Int t = teamRank + tLocal*teamSize;
            ++sendSizes[ map.RowOwner(sep.inds[t]) ];
        }
    }

    // Use a single-entry AllToAll to coordinate how many indices will be 
    // exchanges
    std::vector<int> recvSizes( commSize );
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, comm );

    // Pack the reordered indices
    std::vector<int> sendOffs;
    const int numSends = Scan( sendSizes, sendOffs );
    std::vector<Int> sendInds(numSends), sendOrigInds(numSends);
    auto offs = sendOffs;
    for( int s=0; s<numLocal; ++s )
    {
        const SepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const Int numInds = sepOrLeaf.inds.size();
        for( Int t=0; t<numInds; ++t )
        {
            const Int i = sepOrLeaf.inds[t];
            const Int iMapped = sepOrLeaf.off + t;
            const int q = map.RowOwner(i);
            sendOrigInds[offs[q]] = i;
            sendInds[offs[q]] = iMapped;
            ++offs[q];
        }
    }
    for( Int s=0; s<numDist; ++s )
    {
        const DistSeparator& sep = sepTree.distSeps[s];
        const Int numInds = sep.inds.size();
        const int teamSize = mpi::Size( sep.comm );
        const int teamRank = mpi::Rank( sep.comm );
        const Int numLocalInds = Length( numInds, teamRank, teamSize );
        for( Int tLocal=0; tLocal<numLocalInds; ++tLocal )
        {
            const Int t = teamRank + tLocal*teamSize;
            const Int i = sep.inds[t];
            const Int iMapped = sep.off + t;
            const int q = map.RowOwner(i);
            sendOrigInds[offs[q]] = i;
            sendInds[offs[q]] = iMapped;
            ++offs[q];
        }
    }

    // Perform an AllToAll to exchange the reordered indices
    std::vector<int> recvOffs;
    const int numRecvs = Scan( recvSizes, recvOffs );
    DEBUG_ONLY(
        const Int numLocalSources = map.NumLocalSources();
        if( numRecvs != numLocalSources )
            LogicError("incorrect number of recv indices");
    )
    std::vector<Int> recvInds( numRecvs );
    mpi::AllToAll
    ( sendInds.data(), sendSizes.data(), sendOffs.data(),
      recvInds.data(), recvSizes.data(), recvOffs.data(), comm );

    // Perform an AllToAll to exchange the original indices
    std::vector<Int> recvOrigInds( numRecvs );
    mpi::AllToAll
    ( sendOrigInds.data(), sendSizes.data(), sendOffs.data(),
      recvOrigInds.data(), recvSizes.data(), recvOffs.data(), comm );

    // Unpack the indices
    const Int firstLocalSource = graph.FirstLocalSource();
    for( Int s=0; s<numRecvs; ++s )
        map.SetLocal( recvOrigInds[s]-firstLocalSource, recvInds[s] );
}

} // namespace El
