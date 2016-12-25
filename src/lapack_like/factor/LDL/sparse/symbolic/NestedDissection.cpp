/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include <set>

namespace El {
namespace ldl {

void AMDOrder
( const vector<Int>& subOffsets,
  const vector<Int>& subTargets,
        vector<Int>& amdPerm,
        double* control,
        double* info )
{
    EL_DEBUG_CSE
    const Int numSources = subOffsets.size()-1;
    // TODO(poulson): Simplify this after templating ElSuiteSparse's AMD
#ifdef EL_USE_64BIT_INTS
    const Int numEdges = subTargets.size();
    vector<int> subOffsets_int( numSources+1 ),
                subTargets_int( numEdges ),
                amdPerm_int( numSources );
    for( Int j=0; j<numSources+1; ++j )
        subOffsets_int[j] = int(subOffsets[j]);
    for( Int j=0; j<numEdges; ++j )
        subTargets_int[j] = int(subTargets[j]);
    const int amdStatus =
      El_amd_order
      ( int(numSources),
        subOffsets_int.data(), subTargets_int.data(),
        amdPerm_int.data(),
        control, info );
    if( amdStatus != EL_AMD_OK )
        RuntimeError("AMD status was ",amdStatus);
    amdPerm.resize( numSources );
    for( Int j=0; j<numSources; ++j )
        amdPerm[j] = Int(amdPerm_int[j]);
#else
    amdPerm.resize( numSources );
    const int amdStatus =
      El_amd_order
      ( numSources, subOffsets.data(), subTargets.data(), amdPerm.data(),
        control, info );
    if( amdStatus != EL_AMD_OK )
        RuntimeError("AMD status was ",amdStatus);
#endif
}

inline bool IsSymmetric( const Graph& graph )
{
    EL_DEBUG_CSE
    // NOTE: We only check within the numSources x numSources upper-left
    const Int numSources = graph.NumSources();
    const Int numEdges = graph.NumEdges();
    bool isSymmetric = true;
    for( Int e=0; e<numEdges; ++e  )
    {
        const Int source = graph.Source(e);
        const Int target = graph.Target(e);
        if( source < numSources && target < numSources )
        {
            if( !graph.EdgeExists(target,source) )
                isSymmetric = false;
        }
    }
    return isSymmetric;
}

inline void
NestedDissectionRecursion
( const Graph& graph,
  const vector<Int>& perm,
        Separator& sep,
        NodeInfo& info,
        Int off,
  const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int numSources = graph.NumSources();
    const Int* offsetBuf = graph.LockedOffsetBuffer();
    const Int* sourceBuf = graph.LockedSourceBuffer();
    const Int* targetBuf = graph.LockedTargetBuffer();
    if( numSources <= ctrl.cutoff )
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
        while( sourceOff <= numSources )
        { subOffsets[sourceOff++] = validCounter; }

        // Technically, SuiteSparse expects column-major storage, but since
        // the matrix is structurally symmetric, it's okay to pass in the
        // row-major representation
        vector<Int> amdPerm;
        AMDOrder( subOffsets, subTargets, amdPerm );

        // Compute the symbolic factorization of this leaf node using the
        // reordering just computed
        info.LOffsets.resize( numSources+1 );
        info.LParents.resize( numSources );
        vector<Int> LNnz( numSources ), Flag( numSources ),
                    amdPermInv( numSources );
        suite_sparse::ldl::Symbolic
        ( numSources, subOffsets.data(), subTargets.data(),
          info.LOffsets.data(), info.LParents.data(), LNnz.data(),
          Flag.data(), amdPerm.data(), amdPermInv.data() );

        // Fill in this node of the local separator tree
        sep.off = off;
        sep.inds.resize( numSources );
        for( Int i=0; i<numSources; ++i )
            sep.inds[i] = perm[amdPerm[i]];
        // TODO(poulson): Replace with better deletion mechanism
        SwapClear( sep.children );

        // Fill in this node of the local elimination tree
        info.size = numSources;
        info.off = off;
        // TODO(poulson): Replace with better deletion mechanism
        SwapClear( info.children );
        set<Int> lowerStruct;
        for( Int s=0; s<info.size; ++s )
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
        CopySTL( lowerStruct, info.origLowerStruct );
    }
    else
    {
        EL_DEBUG_ONLY(
          if( !IsSymmetric(graph) )
          {
              Print( graph, "graph" );
              LogicError("Graph was not symmetric");
          }
        )

        // Partition the graph and construct the inverse map
        Graph leftChild, rightChild;
        vector<Int> map;
        const Int sepSize = Bisect( graph, leftChild, rightChild, map, ctrl );
        vector<Int> invMap( numSources );
        for( Int s=0; s<numSources; ++s )
            invMap[map[s]] = s;

        EL_DEBUG_ONLY(
          if( !IsSymmetric(leftChild) )
          {
              Print( graph, "graph" );
              Print( leftChild, "leftChild" );
              LogicError("Left child was not symmetric");
          }
        )
        EL_DEBUG_ONLY(
          if( !IsSymmetric(rightChild) )
          {
              Print( graph, "graph" );
              Print( rightChild, "rightChild" );
              LogicError("Right child was not symmetric");
          }
        )

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
        info.size = sepSize;
        info.off = sep.off;
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
        CopySTL( lowerStruct, info.origLowerStruct );

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

        sep.children.reserve( 2 );
        info.children.reserve( 2 );

        sep.children.emplace_back( new Separator(&sep) );
        info.children.emplace_back( new NodeInfo(&info) );
        NestedDissectionRecursion
        ( leftChild, leftPerm,
          *sep.children.back(), *info.children.back(),
          off, ctrl );

        sep.children.emplace_back( new Separator(&sep) );
        info.children.emplace_back( new NodeInfo(&info) );
        NestedDissectionRecursion
        ( rightChild, rightPerm,
          *sep.children.back(), *info.children.back(),
          off+leftChildSize, ctrl );
    }
}

inline void
NestedDissectionRecursion
( const DistGraph& graph,
  const DistMap& perm,
        DistSeparator& sep,
        DistNodeInfo& info,
        Int off,
  const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = graph.Grid();
    mpi::Comm comm = grid.Comm();
    const int commSize = grid.Size();

    if( commSize > 1 )
    {
        const Int numLocalSources = graph.NumLocalSources();
        const Int firstLocalSource = graph.FirstLocalSource();
        const Int* offsetBuf = graph.LockedOffsetBuffer();
        const Int* targetBuf = graph.LockedTargetBuffer();

        // Partition the graph and construct the inverse map
        DistGraph child;
        bool childIsOnLeft;
        DistMap map(grid);
        info.child.reset( new DistNodeInfo(&info) );
        unique_ptr<Grid> childGrid;
        const Int sepSize =
          Bisect( graph, childGrid, child, map, childIsOnLeft, ctrl );
        info.child->AssignGrid( childGrid );
        info.child->onLeft = childIsOnLeft;
        const Int numSources = graph.NumSources();
        const Int childSize = child.NumSources();
        const Int leftChildSize =
          childIsOnLeft ? childSize : numSources-sepSize-childSize;

        DistMap invMap(grid);
        InvertMap( map, invMap );

        // Mostly fill this node of the DistSeparatorTree
        // (we will finish computing the separator indices at the end)
        sep.off = off + (numSources-sepSize);
        sep.inds.resize( sepSize );
        for( Int s=0; s<sepSize; ++s )
            sep.inds[s] = s + (numSources-sepSize);
        invMap.Translate( sep.inds );

        // Fill in this node of the DistNode
        info.size = sepSize;
        info.off = sep.off;

        set<Int> localLowerStruct;
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
        CopySTL( lowerStruct, info.origLowerStruct );

        // Finish computing the separator indices
        perm.Translate( sep.inds );

        // Construct map from child indices to the original ordering
        DistMap newPerm( child.NumSources(), child.Grid() );
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
        const Int childOff = childIsOnLeft ? off : off+leftChildSize;
        sep.child.reset( new DistSeparator(&sep) );
        NestedDissectionRecursion
        ( child, newPerm, *sep.child, *info.child, childOff, ctrl );
    }
    else
    {
        Graph seqGraph( graph );

        sep.duplicate.reset( new Separator(&sep) );
        info.duplicate.reset( new NodeInfo(&info) );
        NestedDissectionRecursion
        ( seqGraph, perm.Map(), *sep.duplicate, *info.duplicate, off, ctrl );

        // Pull information up from the duplicates
        sep.off = sep.duplicate->off;
        sep.inds = sep.duplicate->inds;
        info.size = info.duplicate->size;
        info.off = info.duplicate->off;
        info.origLowerStruct = info.duplicate->origLowerStruct;
    }
}

void NestedDissection
( const Graph& graph,
        vector<Int>& map,
        Separator& sep,
        NodeInfo& info,
  const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE

    const Int numSources = graph.NumSources();
    vector<Int> perm(numSources);
    for( Int s=0; s<numSources; ++s )
        perm[s] = s;

    NestedDissectionRecursion( graph, perm, sep, info, 0, ctrl );

    // Construct the distributed reordering
    sep.BuildMap( map );
    EL_DEBUG_ONLY(EnsurePermutation(map))

    // Run the symbolic analysis
    Analysis( info );
}

void NestedDissection
( const DistGraph& graph,
        DistMap& map,
        DistSeparator& sep,
        DistNodeInfo& info,
  const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE

    DistMap perm( graph.NumSources(), graph.Grid() );
    const Int firstLocalSource = perm.FirstLocalSource();
    const Int numLocalSources = perm.NumLocalSources();
    for( Int s=0; s<numLocalSources; ++s )
        perm.SetLocal( s, s+firstLocalSource );

    info.SetRootGrid( graph.Grid() );
    NestedDissectionRecursion( graph, perm, sep, info, 0, ctrl );

    // Construct the distributed reordering
    sep.BuildMap( info, map );
    EL_DEBUG_ONLY(EnsurePermutation(map))

    // Run the symbolic analysis
    Analysis( info, ctrl.storeFactRecvInds );
}

} // namespace ldl
} // namespace El
