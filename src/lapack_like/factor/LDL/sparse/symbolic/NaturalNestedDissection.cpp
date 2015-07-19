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
        Node& node,
        Int off, 
        Int cutoff )
{
    DEBUG_ONLY(CSE cse("NaturalNestedDissectionRecursion"))
    if( graph.NumSources() <= cutoff )
    {
        // Fill in this node of the local separator tree
        const Int numSources = graph.NumSources();
        sep.off = off;
        sep.inds = perm;

        // Fill in this node of the local elimination tree
        node.size = numSources;
        node.off = off;
        set<Int> lowerStruct;
        for( Int s=0; s<node.size; ++s )
        {
            const Int numConnections = graph.NumConnections( s );
            const Int edgeOff = graph.SourceOffset( s );
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
        Int nxLeft, nyLeft, nzLeft, nxRight, nyRight, nzRight;
        Graph leftChild, rightChild;
        vector<Int> map;
        const Int sepSize = 
            NaturalBisect
            ( nx, ny, nz, graph, 
              nxLeft, nyLeft, nzLeft, leftChild, 
              nxRight, nyRight, nzRight, rightChild, map );
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
            const Int edgeOff = graph.SourceOffset( source );
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
        node.children[0] = new Node(&node);
        node.children[1] = new Node(&node);
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
        DistNode& node,
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
        const Int numLocalSources = graph.NumLocalSources();
        const Int firstLocalSource = graph.FirstLocalSource();
        set<Int> localStructSet;
        for( Int s=0; s<sepSize; ++s )
        {
            const Int source = sep.inds[s];
            if( source >= firstLocalSource && 
                source < firstLocalSource+numLocalSources )
            {
                const Int localSource = source - firstLocalSource;
                const Int numConnect = graph.NumConnections( localSource );
                const Int localOff = graph.SourceOffset( localSource );
                for( Int t=0; t<numConnect; ++t )
                {
                    const Int target = graph.Target( localOff+t );
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
        CopySTL( structSet, node.lowerStruct );

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
        const Int newOff = ( childIsOnLeft ? off : off+leftChildSize );
        sep.child = new DistSeparator(&sep);
        node.child = new DistNode(&node);
        node.child->onLeft = childIsOnLeft;
        NaturalNestedDissectionRecursion
        ( nxChild, nyChild, nzChild, child, newPerm, 
          *sep.child, *node.child, newOff, cutoff );
    }
    else
    {
        Graph seqGraph( graph );

        sep.duplicate = new Separator(&sep);
        node.duplicate = new Node(&node);

        NaturalNestedDissectionRecursion
        ( nx, ny, nz, seqGraph, perm.Map(), 
          *sep.duplicate, *node.duplicate, off, cutoff );

        // Pull information up from the duplicates
        sep.off = sep.duplicate->off;
        sep.inds = sep.duplicate->inds;
        node.size = node.duplicate->size;
        node.off = node.duplicate->off;
        node.lowerStruct = node.duplicate->lowerStruct;
    }
}

void NaturalNestedDissection
(       Int nx,
        Int ny,
        Int nz,
  const Graph& graph, 
        vector<Int>& map,
        Separator& sep, 
        NodeInfo& info,
        Int cutoff )
{
    DEBUG_ONLY(CSE cse("NaturalNestedDissection"))
    // NOTE: There is a potential memory leak here if sep or info is reused

    const Int numSources = graph.NumSources();
    vector<Int> perm( numSources );
    for( Int s=0; s<numSources; ++s )
        perm[s] = s;

    Node node;
    NaturalNestedDissectionRecursion
    ( nx, ny, nz, graph, perm, sep, node, 0, cutoff );

    // Construct the distributed reordering    
    BuildMap( sep, map );
    DEBUG_ONLY(EnsurePermutation( map ))

    // Run the symbolic analysis
    Analysis( node, info );
}

void NaturalNestedDissection
(       Int nx,
        Int ny,
        Int nz,
  const DistGraph& graph, 
        DistMap& map,
        DistSeparator& sep, 
        DistNodeInfo& info,
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

    DistNode node;
    NaturalNestedDissectionRecursion
    ( nx, ny, nz, graph, perm, sep, node, 0, cutoff );

    // Construct the distributed reordering    
    BuildMap( sep, map );
    DEBUG_ONLY(EnsurePermutation(map))

    // Run the symbolic analysis
    Analysis( node, info, storeFactRecvInds );
}

} // namespace ldl
} // namespace El
