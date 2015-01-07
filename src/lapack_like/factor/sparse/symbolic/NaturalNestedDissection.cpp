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

inline void
NaturalNestedDissectionRecursion
(       Int nx,
        Int ny,
        Int nz,
  const Graph& graph, 
  const std::vector<Int>& perm,
        DistSeparatorTree& sepTree, 
        DistSymmElimTree& eTree,
        Int parent, 
        Int off, 
        Int cutoff )
{
    DEBUG_ONLY(CallStackEntry cse("NaturalNestedDissectionRecursion"))
    if( graph.NumSources() <= cutoff )
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
        Int nxLeft, nyLeft, nzLeft, nxRight, nyRight, nzRight;
        Graph leftChild, rightChild;
        std::vector<Int> map;
        const Int sepSize = 
            NaturalBisect
            ( nx, ny, nz, graph, 
              nxLeft, nyLeft, nzLeft, leftChild, 
              nxRight, nyRight, nzRight, rightChild, map );
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
        NaturalNestedDissectionRecursion
        ( nxRight, nyRight, nzRight, rightChild, rightPerm, sepTree, eTree, 
          parent, off+leftChildSize, cutoff );
        node.children[0] = eTree.localNodes.size();
        NaturalNestedDissectionRecursion
        ( nxLeft, nyLeft, nzLeft, leftChild, leftPerm, sepTree, eTree, 
          parent, off, cutoff );
    }
}

inline void
NaturalNestedDissectionRecursion
(       Int nx,
        Int ny,
        Int nz,
  const DistGraph& graph, 
  const DistMap& perm,
        DistSeparatorTree& sepTree, 
        DistSymmElimTree& eTree,
        Int depth, 
        Int off, 
        bool onLeft,
        Int cutoff )
{
    DEBUG_ONLY(CallStackEntry cse("NaturalNestedDissectionRecursion"))
    const int distDepth = sepTree.distSeps.size();
    mpi::Comm comm = graph.Comm();
    if( distDepth - depth > 0 )
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
        Int sumOfLocalConnectedSizes=0;
        std::vector<int> localConnectedOffs( commSize );
        for( int q=0; q<commSize; ++q )
        {
            localConnectedOffs[q] = sumOfLocalConnectedSizes;
            sumOfLocalConnectedSizes += localConnectedSizes[q];
        }
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
        NaturalNestedDissectionRecursion
        ( nxChild, nyChild, nzChild, child, newPerm, sepTree, eTree, depth+1, 
          newOff, childIsOnLeft, cutoff );
    }
    else if( graph.NumSources() <= cutoff )
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
        Int nxLeft, nyLeft, nzLeft, nxRight, nyRight, nzRight;
        Graph leftChild, rightChild;
        std::vector<Int> map;
        const Int sepSize = 
            NaturalBisect
            ( nx, ny, nz, seqGraph, 
              nxLeft, nyLeft, nzLeft, leftChild, 
              nxRight, nyRight, nzRight, rightChild, map );
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
        NaturalNestedDissectionRecursion
        ( nxRight, nyRight, nzRight, rightChild, rightPerm, sepTree, eTree, 
          parent, off+leftChildSize, cutoff );
        localNode.children[0] = eTree.localNodes.size();
        NaturalNestedDissectionRecursion
        ( nxLeft, nyLeft, nzLeft, leftChild, leftPerm, sepTree, eTree, 
          parent, off, cutoff );
    }
}

void NaturalNestedDissection
(       Int nx,
        Int ny,
        Int nz,
  const DistGraph& graph, 
        DistMap& map,
        DistSeparatorTree& sepTree, 
        DistSymmInfo& info,
        Int cutoff, 
        bool storeFactRecvInds )
{
    DEBUG_ONLY(CallStackEntry cse("NaturalNestedDissection"))
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
    NaturalNestedDissectionRecursion
    ( nx, ny, nz, graph, perm, sepTree, eTree, 0, 0, false, cutoff );

    ReverseOrder( sepTree, eTree );

    // Construct the distributed reordering    
    BuildMap( graph, sepTree, map );
    DEBUG_ONLY(EnsurePermutation( map ))

    // Run the symbolic analysis
    SymmetricAnalysis( eTree, info, storeFactRecvInds );
}

Int NaturalBisect
(       Int nx, 
        Int ny, 
        Int nz,
  const Graph& graph, 
        Int& nxLeft,
        Int& nyLeft,
        Int& nzLeft,
        Graph& leftChild, 
        Int& nxRight,
        Int& nyRight,
        Int& nzRight,
        Graph& rightChild,
        std::vector<Int>& perm )
{
    DEBUG_ONLY(CallStackEntry cse("NaturalBisect"))
    const Int numSources = graph.NumSources();
    if( numSources == 0 )
        LogicError("There is no reason to bisect an empty sequential graph");

    Int leftChildSize, rightChildSize, sepSize;
    perm.resize( numSources );
    if( nx >= ny && nx >= nz )
    {
        nxLeft = (nx-1)/2;
        nyLeft = ny;
        nzLeft = nz;
        leftChildSize = nxLeft*nyLeft*nzLeft;

        nxRight = nx-1-nxLeft;
        nyRight = ny;
        nzRight = nz;
        rightChildSize = nxRight*nyRight*nzRight;

        sepSize = ny*nz;

        // Fill the left side
        Int off=0;
        for( Int z=0; z<nz; ++z )
            for( Int y=0; y<ny; ++y )
                for( Int x=0; x<nxLeft; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the right side
        off = leftChildSize;
        for( Int z=0; z<nz; ++z )
            for( Int y=0; y<ny; ++y )
                for( Int x=nxLeft+1; x<nx; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the separator
        off=leftChildSize+rightChildSize;
        for( Int z=0; z<nz; ++z )
            for( Int y=0; y<ny; ++y )
                perm[nxLeft+y*nx+z*nx*ny] = off++;
    }
    else if( ny >= nx && ny >= nz )
    {
        nxLeft = nx;
        nyLeft = (ny-1)/2;
        nzLeft = nz;
        leftChildSize = nxLeft*nyLeft*nzLeft;

        nxRight = nx;
        nyRight = ny-1-nyLeft;
        nzRight = nz;
        rightChildSize = nxRight*nyRight*nzRight;

        sepSize = nx*nz;

        // Fill the left side
        Int off=0;
        for( Int z=0; z<nz; ++z )
            for( Int y=0; y<nyLeft; ++y )
                for( Int x=0; x<nx; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the right side
        off = leftChildSize;
        for( Int z=0; z<nz; ++z )
            for( Int y=nyLeft+1; y<ny; ++y )
                for( Int x=0; x<nx; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the separator
        off=leftChildSize+rightChildSize;
        for( Int z=0; z<nz; ++z )
            for( Int x=0; x<nx; ++x )
                perm[x+nyLeft*nx+z*nx*ny] = off++;
    }
    else
    {
        nxLeft = nx;
        nyLeft = ny;
        nzLeft = (nz-1)/2;
        leftChildSize = nxLeft*nyLeft*nzLeft;

        nxRight = nx;
        nyRight = ny;
        nzRight = nz-1-nzLeft;
        rightChildSize = nxRight*nyRight*nzRight;

        sepSize = nx*ny;

        // Fill the left side
        Int off=0;
        for( Int z=0; z<nzLeft; ++z )
            for( Int y=0; y<ny; ++y )
                for( Int x=0; x<nx; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the right side
        off = leftChildSize;
        for( Int z=nzLeft+1; z<nz; ++z )
            for( Int y=0; y<ny; ++y )
                for( Int x=0; x<nx; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the separator
        off=leftChildSize+rightChildSize;
        for( Int y=0; y<ny; ++y )
            for( Int x=0; x<nx; ++x )
                perm[x+y*nx+nzLeft*nx*ny] = off++;
    }
    DEBUG_ONLY(EnsurePermutation( perm ))
    BuildChildrenFromPerm
    ( graph, perm, leftChildSize, leftChild, rightChildSize, rightChild );
    return sepSize;
}

Int NaturalBisect
(       Int nx,
        Int ny,
        Int nz,
  const DistGraph& graph, 
        Int& nxChild,
        Int& nyChild,
        Int& nzChild,
        DistGraph& child, 
        DistMap& perm,
        bool& onLeft )
{
    DEBUG_ONLY(CallStackEntry cse("NaturalBisect"))
    const Int numSources = graph.NumSources();
    const Int firstLocalSource = graph.FirstLocalSource();
    const Int numLocalSources = graph.NumLocalSources();
    mpi::Comm comm = graph.Comm();
    const Int commSize = mpi::Size( comm );
    if( commSize == 1 )
        LogicError
        ("This routine assumes at least two processes are used, "
         "otherwise one child will be lost");

    Int leftChildSize, rightChildSize, sepSize;
    Int nxLeft, nyLeft, nzLeft, nxRight, nyRight, nzRight;
    perm.SetComm( comm );
    perm.Resize( numSources );
    if( nx != 0 && ny != 0 && nz != 0 )
    {
        if( nx >= ny && nx >= nz )
        {
            nxLeft = (nx-1)/2;
            nyLeft = ny;
            nzLeft = nz;
            leftChildSize = nxLeft*nyLeft*nzLeft;

            nxRight = nx-1-nxLeft;
            nyRight = ny;
            nzRight = nz;
            rightChildSize = nxRight*nyRight*nzRight;

            sepSize = ny*nz;

            const Int rightOff=leftChildSize, 
                      sepOff=leftChildSize+rightChildSize;
            for( Int iLocal=0; iLocal<numLocalSources; ++iLocal )
            {
                const Int i = iLocal + firstLocalSource;
                const Int x = i % nx;
                const Int y = (i/nx) % ny;
                const Int z = i/(nx*ny);
                if( x < nxLeft )
                {
                    const Int xLeft = x;
                    const Int leftInd = xLeft + y*nxLeft + z*nxLeft*ny;
                    perm.SetLocal( iLocal, leftInd );
                }
                else if( x > nxLeft )
                {
                    const Int xRight = x-(nxLeft+1);
                    const Int rightInd = xRight + y*nxRight + z*nxRight*ny;
                    perm.SetLocal( iLocal, rightOff+rightInd );
                }
                else
                {
                    const Int sepInd = y + z*ny;
                    perm.SetLocal( iLocal, sepOff+sepInd );
                }
            }
        }
        else if( ny >= nx && ny >= nz )
        {
            nxLeft = nx;
            nyLeft = (ny-1)/2;
            nzLeft = nz;
            leftChildSize = nxLeft*nyLeft*nzLeft;

            nxRight = nx;
            nyRight = ny-1-nyLeft;
            nzRight = nz;
            rightChildSize = nxRight*nyRight*nzRight;

            sepSize = nx*nz;

            const Int rightOff=leftChildSize, 
                      sepOff=leftChildSize+rightChildSize;
            for( Int iLocal=0; iLocal<numLocalSources; ++iLocal )
            {
                const Int i = iLocal + firstLocalSource;
                const Int x = i % nx;
                const Int y = (i/nx) % ny;
                const Int z = i/(nx*ny);
                if( y < nyLeft )
                {
                    const Int yLeft = y;
                    const Int leftInd = x + yLeft*nx + z*nx*nyLeft;
                    perm.SetLocal( iLocal, leftInd );
                }
                else if( y > nyLeft )
                {
                    const Int yRight = y - (nyLeft+1);
                    const Int rightInd = x + yRight*nx + z*nx*nyRight;
                    perm.SetLocal( iLocal, rightOff+rightInd );
                }
                else
                {
                    const Int sepInd = x + z*nx;
                    perm.SetLocal( iLocal, sepOff+sepInd );
                }
            }
        }
        else
        {
            nxLeft = nx;
            nyLeft = ny;
            nzLeft = (nz-1)/2;
            leftChildSize = nxLeft*nyLeft*nzLeft;

            nxRight = nx;
            nyRight = ny;
            nzRight = nz-1-nzLeft;
            rightChildSize = nxRight*nyRight*nzRight;

            sepSize = nx*ny;

            const Int rightOff=leftChildSize, 
                      sepOff=leftChildSize+rightChildSize;
            for( Int iLocal=0; iLocal<numLocalSources; ++iLocal )
            {
                const Int i = iLocal + firstLocalSource;
                const Int x = i % nx;
                const Int y = (i/nx) % ny;
                const Int z = i/(nx*ny);
                if( z < nzLeft )
                {
                    const Int zLeft = z;
                    const Int leftInd = x + y*nx + zLeft*nx*ny;
                    perm.SetLocal( iLocal, leftInd );
                }
                else if( z > nzLeft )
                {
                    const Int zRight = z - (nzLeft+1);
                    const Int rightInd = x + y*nx + zRight*nx*ny;
                    perm.SetLocal( iLocal, rightOff+rightInd );
                }
                else
                {
                    const Int sepInd = x + y*nx;
                    perm.SetLocal( iLocal, sepOff+sepInd );
                }
            }
        }
    }
    else
    {
        leftChildSize = rightChildSize = sepSize = 0;
        nxLeft = nx;
        nyLeft = ny;
        nzLeft = nz;
        nxRight = nx;
        nyRight = ny;
        nzRight = nz;
    }
    DEBUG_ONLY(EnsurePermutation( perm ))

    BuildChildFromPerm
    ( graph, perm, leftChildSize, rightChildSize, onLeft, child );

    if( onLeft )
    {
        nxChild = nxLeft;
        nyChild = nyLeft;
        nzChild = nzLeft;
    }
    else
    {
        nxChild = nxRight;
        nyChild = nyRight;
        nzChild = nzRight;
    }
    return sepSize;
}

} // namespace El
