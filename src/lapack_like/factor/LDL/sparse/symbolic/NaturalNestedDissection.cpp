/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include <set>

namespace El {
namespace ldl {

inline void
NaturalNestedDissectionRecursion
(       Int nx,
        Int ny,
        Int nz,
  const Graph& graph, 
  const vector<Int>& perm,
        Separator& sep, 
        NodeInfo& node,
        Int off, 
        Int cutoff )
{
    DEBUG_ONLY(CSE cse("NaturalNestedDissectionRecursion"))
    const Int numSources = graph.NumSources();
    const Int* offsetBuf = graph.LockedOffsetBuffer();
    const Int* sourceBuf = graph.LockedSourceBuffer();
    const Int* targetBuf = graph.LockedTargetBuffer();
    if( numSources <= cutoff )
    {
        // Filter out the graph of the diagonal block
        Int numValidEdges = 0;
        const Int numEdges = graph.NumEdges();
        for( Int e=0; e<numEdges; ++e )
            if( targetBuf[e] < numSources )
                ++numValidEdges;
        vector<Int> subOffsets(numSources+1), subTargets(Max(numValidEdges,1));
        Int sourceOff = 0;
        Int validCounter = 0;
        Int prevSource = -1;
        for( Int e=0; e<numEdges; ++e )
        {
            const Int source = sourceBuf[e];
            const Int target = targetBuf[e];
            while( source != prevSource )
            {
                subOffsets[sourceOff++] = validCounter;
                ++prevSource;
            }
            if( target < numSources )
                subTargets[validCounter++] = target;
        }
        while( sourceOff <= numSources)
        { subOffsets[sourceOff++] = validCounter; }

        // Technically, SuiteSparse expects column-major storage, but since
        // the matrix is structurally symmetric, it's okay to pass in the 
        // row-major representation
        vector<Int> amdPerm;
        AMDOrder( subOffsets, subTargets, amdPerm );

        // Compute the symbolic factorization of this leaf node using the
        // reordering just computed
        node.LOffsets.resize( numSources+1 );
        node.LParents.resize( numSources );
        vector<Int> LNnz( numSources ), Flag( numSources ),
                    amdPermInv( numSources );
        suite_sparse::ldl::Symbolic 
        ( numSources, subOffsets.data(), subTargets.data(),
          node.LOffsets.data(), node.LParents.data(), LNnz.data(),
          Flag.data(), amdPerm.data(), amdPermInv.data() );

        // Fill in this node of the local separator tree
        sep.off = off;
        sep.inds.resize( numSources );
        for( Int i=0; i<numSources; ++i )
            sep.inds[i] = perm[amdPerm[i]];

        // Fill in this node of the local elimination tree
        node.size = numSources;
        node.off = off;
        set<Int> lowerStruct;
        for( Int s=0; s<node.size; ++s )
        {
            const Int edgeOff = offsetBuf[s];
            const Int numConn = offsetBuf[s+1] - edgeOff;
            for( Int t=0; t<numConn; ++t )
            {
                const Int target = targetBuf[edgeOff+t];
                if( target >= numSources )
                    lowerStruct.insert( off+target );
            }
        }
        CopySTL( lowerStruct, node.origLowerStruct );
    }
    else
    {
        // Partition the graph and construct the inverse map
        Int nxLeft, nyLeft, nzLeft, nxRight, nyRight, nzRight;
        Graph leftChild, rightChild;
        vector<Int> map;
        const Int sepSize = 
            NaturalBisect
            ( nx, ny, nz, graph, 
              nxLeft, nyLeft, nzLeft, leftChild, 
              nxRight, nyRight, nzRight, rightChild, map );
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
            const Int edgeOff = offsetBuf[source];
            const Int numConn = offsetBuf[source+1] - edgeOff;
            for( Int t=0; t<numConn; ++t )
            {
                const Int target = targetBuf[edgeOff+t];
                if( target >= numSources )
                    lowerStruct.insert( off+target );
            }
        }
        CopySTL( lowerStruct, node.origLowerStruct );

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
        node.children[0] = new NodeInfo(&node);
        node.children[1] = new NodeInfo(&node);
        NaturalNestedDissectionRecursion
        ( nxLeft, nyLeft, nzLeft, leftChild, leftPerm, 
          *sep.children[0], *node.children[0], off, cutoff );
        NaturalNestedDissectionRecursion
        ( nxRight, nyRight, nzRight, rightChild, rightPerm, 
          *sep.children[1], *node.children[1], off+leftChildSize, cutoff );
    }
}

inline void
NaturalNestedDissectionRecursion
(       Int nx,
        Int ny,
        Int nz,
  const DistGraph& graph, 
  const DistMap& perm,
        DistSeparator& sep, 
        DistNodeInfo& node,
        Int off, 
        Int cutoff )
{
    DEBUG_ONLY(CSE cse("NaturalNestedDissectionRecursion"))
    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::Size(comm);

    mpi::Dup( comm, sep.comm );
    mpi::Dup( comm, node.comm );

    if( commSize > 1 )
    {
        const Int numLocalSources = graph.NumLocalSources();
        const Int firstLocalSource = graph.FirstLocalSource();
        const Int* offsetBuf = graph.LockedOffsetBuffer();
        const Int* targetBuf = graph.LockedTargetBuffer();

        // Partition the graph and construct the inverse map
        Int nxChild, nyChild, nzChild;
        DistGraph child;
        bool childIsOnLeft;
        DistMap map;
        const Int sepSize = 
            NaturalBisect
            ( nx, ny, nz, graph, nxChild, nyChild, nzChild, child, 
              map, childIsOnLeft );
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

        // Fill in this node of the DistNode
        node.size = sepSize;
        node.off = sep.off;

        set<Int> localStructSet;
        for( Int s=0; s<sepSize; ++s )
        {
            const Int source = sep.inds[s];
            if( source >= firstLocalSource && 
                source < firstLocalSource+numLocalSources )
            {
                const Int localSource = source - firstLocalSource;
                const Int edgeOff = offsetBuf[localSource];
                const Int numConn = offsetBuf[localSource+1] - edgeOff;
                for( Int t=0; t<numConn; ++t )
                {
                    const Int target = targetBuf[edgeOff+t];
                    if( target >= numSources )
                        localStructSet.insert( off+target );
                }
            }
        }
        const int localStructSize = localStructSet.size();
        const int commSize = mpi::Size( comm );
        vector<int> localStructSizes( commSize );
        mpi::AllGather( &localStructSize, 1, localStructSizes.data(), 1, comm );
        vector<Int> localStruct;
        CopySTL( localStructSet, localStruct );
        vector<int> localStructOffs;
        int nonUniqueStructSize = Scan( localStructSizes, localStructOffs );
        vector<Int> nonUniqueStruct( nonUniqueStructSize );
        mpi::AllGather
        ( localStruct.data(), localStructSize,
          nonUniqueStruct.data(), 
          localStructSizes.data(), localStructOffs.data(), comm );
        set<Int> structSet( nonUniqueStruct.begin(), nonUniqueStruct.end() );
        CopySTL( structSet, node.origLowerStruct );

        // Finish computing the separator indices
        perm.Translate( sep.inds );

        // Construct map from child indices to the original ordering
        DistMap newPerm( child.NumSources(), child.Comm() );
        const Int localChildSize = child.NumLocalSources();
        const Int firstLocalChildSource = child.FirstLocalSource();
        auto& newPermLoc = newPerm.Map();
        if( childIsOnLeft )
            for( Int s=0; s<localChildSize; ++s )
                newPermLoc[s] = s+firstLocalChildSource;
        else
            for( Int s=0; s<localChildSize; ++s )
                newPermLoc[s] = s+firstLocalChildSource+leftChildSize;
        invMap.Extend( newPerm );
        perm.Extend( newPerm );

        // Recurse
        const Int newOff = ( childIsOnLeft ? off : off+leftChildSize );
        sep.child = new DistSeparator(&sep);
        node.child = new DistNodeInfo(&node);
        node.child->onLeft = childIsOnLeft;
        NaturalNestedDissectionRecursion
        ( nxChild, nyChild, nzChild, child, newPerm, 
          *sep.child, *node.child, newOff, cutoff );
    }
    else
    {
        Graph seqGraph( graph );

        sep.duplicate = new Separator(&sep);
        node.duplicate = new NodeInfo(&node);

        NaturalNestedDissectionRecursion
        ( nx, ny, nz, seqGraph, perm.Map(), 
          *sep.duplicate, *node.duplicate, off, cutoff );

        // Pull information up from the duplicates
        sep.off = sep.duplicate->off;
        sep.inds = sep.duplicate->inds;
        node.size = node.duplicate->size;
        node.off = node.duplicate->off;
        node.origLowerStruct = node.duplicate->origLowerStruct;
    }
}

void NaturalNestedDissection
(       Int nx,
        Int ny,
        Int nz,
  const Graph& graph, 
        vector<Int>& map,
        Separator& sep, 
        NodeInfo& node,
        Int cutoff )
{
    DEBUG_ONLY(CSE cse("NaturalNestedDissection"))
    // NOTE: There is a potential memory leak here if sep or info is reused

    const Int numSources = graph.NumSources();
    vector<Int> perm( numSources );
    for( Int s=0; s<numSources; ++s )
        perm[s] = s;

    NaturalNestedDissectionRecursion
    ( nx, ny, nz, graph, perm, sep, node, 0, cutoff );

    // Construct the distributed reordering    
    BuildMap( sep, map );
    DEBUG_ONLY(EnsurePermutation( map ))

    // Run the symbolic analysis
    Analysis( node );
}

void NaturalNestedDissection
(       Int nx,
        Int ny,
        Int nz,
  const DistGraph& graph, 
        DistMap& map,
        DistSeparator& sep, 
        DistNodeInfo& node,
        Int cutoff, 
        bool storeFactRecvInds )
{
    DEBUG_ONLY(CSE cse("NaturalNestedDissection"))
    // NOTE: There is a potential memory leak here if sep or info is reused 

    DistMap perm( graph.NumSources(), graph.Comm() );
    const Int firstLocalSource = perm.FirstLocalSource();
    const Int numLocalSources = perm.NumLocalSources();
    for( Int s=0; s<numLocalSources; ++s )
        perm.SetLocal( s, s+firstLocalSource );

    NaturalNestedDissectionRecursion
    ( nx, ny, nz, graph, perm, sep, node, 0, cutoff );

    // Construct the distributed reordering    
    BuildMap( sep, map );
    DEBUG_ONLY(EnsurePermutation(map))

    // Run the symbolic analysis
    Analysis( node, storeFactRecvInds );
}

} // namespace ldl
} // namespace El
