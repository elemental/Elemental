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
NestedDissectionRecursion
( const Graph& graph, 
  const vector<Int>& perm,
        Separator& sep, 
        Node& node,
        Int off, 
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CSE cse("ldl::NestedDissectionRecursion"))
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
        DistNode& node,
        Int off, 
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CSE cse("ldl::NestedDissectionRecursion"))
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

        // Fill in this node of the DistNode
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
                const Int localOff = graph.SourceOffset( localSource );
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
        node.child = new DistNode(&node);
        node.child->onLeft = childIsOnLeft;
        NestedDissectionRecursion
        ( child, newPerm, *sep.child, *node.child, childOff, ctrl );
    }
    else
    {
        Graph seqGraph( graph );

        sep.duplicate = new Separator(&sep);
        node.duplicate = new Node(&node);

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
        NodeInfo& info,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CSE cse("ldl::NestedDissection"))
    // NOTE: There is a potential memory leak here if sep or info is reused

    const Int numSources = graph.NumSources();
    vector<Int> perm(numSources);
    for( Int s=0; s<numSources; ++s )
        perm[s] = s;

    Node node;
    NestedDissectionRecursion( graph, perm, sep, node, 0, ctrl );

    // Construct the distributed reordering    
    BuildMap( sep, map );
    DEBUG_ONLY(EnsurePermutation(map))

    // Run the symbolic analysis
    Analysis( node, info );
}

void NestedDissection
( const DistGraph& graph, 
        DistMap& map,
        DistSeparator& sep, 
        DistNodeInfo& info,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CSE cse("ldl::NestedDissection"))
    // NOTE: There is a potential memory leak here if sep or info is reused

    DistMap perm( graph.NumSources(), graph.Comm() );
    const Int firstLocalSource = perm.FirstLocalSource();
    const Int numLocalSources = perm.NumLocalSources();
    for( Int s=0; s<numLocalSources; ++s )
        perm.SetLocal( s, s+firstLocalSource );

    DistNode node;
    NestedDissectionRecursion( graph, perm, sep, node, 0, ctrl );

    // Construct the distributed reordering    
    BuildMap( sep, map );
    DEBUG_ONLY(EnsurePermutation(map))

    // Run the symbolic analysis
    Analysis( node, info, ctrl.storeFactRecvInds );
}

void BuildMap( const Separator& rootSep, vector<Int>& map )
{
    DEBUG_ONLY(CSE cse("ldl::BuildMap"))
    const Int numSources = rootSep.off + rootSep.inds.size();
    map.resize( numSources );

    function<void(const Separator&)> buildMap = 
      [&]( const Separator& sep )
      {
        for( auto* child : sep.children )  
            buildMap( *child );
        for( size_t t=0; t<sep.inds.size(); ++t )
            map[sep.inds[t]] = sep.off + t;
      };
    buildMap( rootSep );
}

void BuildMap( const DistSeparator& rootSep, DistMap& map )
{
    DEBUG_ONLY(CSE cse("ldl::BuildMap"))

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

} // namespace ldl
} // namespace El
