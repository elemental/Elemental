/*
   Copyright 2009-2011, Jack Poulson.
   All rights reserved.

   Copyright 2011-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright 2013-2014, Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   Copyright 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

inline void PairwiseExchangeLowerStruct
( Int& theirSize, vector<Int>& theirLowerStruct,
  const DistSymmNode& node, const DistSymmNodeInfo& childInfo )
{
    // Determine our partner's rank for this exchange in node's communicator
    const int teamRank = mpi::Rank( node.comm );
    const int teamSize = mpi::Size( node.comm );
    const int childTeamRank = mpi::Rank( childInfo.comm );
    const int myTeamSize = mpi::Size( childInfo.comm );
    const int otherTeamSize = teamSize - myTeamSize;
    const bool inFirstTeam = ( teamRank == childTeamRank );
    const int partner =
        ( inFirstTeam ? teamRank+myTeamSize : teamRank-otherTeamSize );

    // SendRecv the message lengths
    const Int mySize = childInfo.size;
    const Int myLowerStructSize = childInfo.lowerStruct.size();
    const Int initialSends[2] = { mySize, myLowerStructSize };
    Int initialRecvs[2];
    mpi::SendRecv
    ( initialSends, 2, partner,
      initialRecvs, 2, partner, node.comm );
    theirSize = initialRecvs[0];
    const Int theirLowerStructSize = initialRecvs[1];

    // Perform the exchange
    theirLowerStruct.resize( theirLowerStructSize );
    mpi::SendRecv
    ( &childInfo.lowerStruct[0], myLowerStructSize, partner,
      &theirLowerStruct[0], theirLowerStructSize, partner, node.comm );
}

inline void BroadcastLowerStruct
( Int& theirSize, vector<Int>& theirLowerStruct,
  const DistSymmNode& node, const DistSymmNodeInfo& childInfo )
{
    // Determine our partner's rank for this exchange in node's communicator
    const int teamRank = mpi::Rank( node.comm );
    const int teamSize = mpi::Size( node.comm );
    const int childTeamRank = mpi::Rank( childInfo.comm );
    const int myTeamSize = mpi::Size( childInfo.comm );
    const int otherTeamSize = teamSize - myTeamSize;
    const bool inFirstTeam = ( teamRank == childTeamRank );

    if( childTeamRank == 0 )
    {
        const int partner =
            ( inFirstTeam ? teamRank+myTeamSize : teamRank-otherTeamSize );

        // SendRecv the message lengths
        const Int mySize = childInfo.size;
        const Int myLowerStructSize = childInfo.lowerStruct.size();
        const Int initialSends[2] = { mySize, myLowerStructSize };
        Int initialRecvs[2];
        mpi::SendRecv
        ( initialSends, 2, partner, 
          initialRecvs, 2, partner, node.comm );
        theirSize = initialRecvs[0];
        const Int theirLowerStructSize = initialRecvs[1];

        // Perform the exchange
        theirLowerStruct.resize( theirLowerStructSize );
        mpi::SendRecv
        ( &childInfo.lowerStruct[0], myLowerStructSize, partner,
          &theirLowerStruct[0], theirLowerStructSize, partner, node.comm );

        // Broadcast the other team's child's sizes
        mpi::Broadcast( initialRecvs, 2, 0, childInfo.comm );

        // Broadcast the other team's child's lower struct
        mpi::Broadcast
        ( &theirLowerStruct[0], theirLowerStructSize, 0, childInfo.comm );
    } 
    else
    {
        // Receive the other team's child's sizes
        Int initialRecvs[2];
        mpi::Broadcast( initialRecvs, 2, 0, childInfo.comm );
        theirSize = initialRecvs[0];
        const Int theirLowerStructSize = initialRecvs[1];

        // Receive the other team's child's lower struct
        theirLowerStruct.resize( theirLowerStructSize );
        mpi::Broadcast
        ( &theirLowerStruct[0], theirLowerStructSize, 0, childInfo.comm );
    }
}

inline void GetLowerStruct
( Int& theirSize, vector<Int>& theirLowerStruct,
  const DistSymmNode& node,
  const DistSymmNodeInfo& childInfo )
{
    const int teamSize = mpi::Size( node.comm );
    const int childTeamSize = mpi::Size( node.child->comm );
    const int leftTeamSize =
        ( node.child->onLeft ? childTeamSize : teamSize-childTeamSize );
    const int rightTeamSize = teamSize - leftTeamSize;
    if( leftTeamSize == rightTeamSize )
        PairwiseExchangeLowerStruct
        ( theirSize, theirLowerStruct, node, childInfo );
    else
        BroadcastLowerStruct
        ( theirSize, theirLowerStruct, node, childInfo );
}

inline void ComputeStructAndRelInds
( Int theirSize, const vector<Int>& theirLowerStruct,
  const DistSymmNode& node, DistSymmNodeInfo& info )
{
    const auto& myLowerStruct = info.child->lowerStruct;
    DEBUG_ONLY(
        if( !IsStrictlySorted(myLowerStruct) )
        {
            if( IsSorted(myLowerStruct) )
                LogicError("Repeat in my lower struct");
            else
                LogicError("My lower struct not sorted");
        }
        if( !IsStrictlySorted(theirLowerStruct) )
        {
            if( IsSorted(theirLowerStruct) )
                LogicError("Repeat in their lower struct");
            else
                LogicError("Their lower struct not sorted");
        }
        if( !IsStrictlySorted(node.lowerStruct) )
        {
            if( IsSorted(node.lowerStruct) )
                LogicError("Repeat in original struct");
            else
                LogicError("Original struct not sorted");
        }
    )

    // Combine the children's structure
    auto childrenStruct = Union( myLowerStruct, theirLowerStruct );

    // Now add in the original lower structure
    auto partialStruct = Union( childrenStruct, node.lowerStruct );

    // Now the node indices
    vector<Int> nodeInds( node.size );
    for( Int i=0; i<node.size; ++i )
        nodeInds[i] = node.off + i;
    auto fullStruct = Union( nodeInds, partialStruct );

    // Construct the relative indices of the original lower structure
    info.origLowerRelInds = RelativeIndices( node.lowerStruct, fullStruct );

    // Construct the relative indices of the children
    info.childSizes.resize(2);
    info.childRelInds.resize(2);
    if( node.child->onLeft )
    {
        info.childSizes[0] = info.child->size;
        info.childSizes[1] = theirSize;
        info.childRelInds[0] = RelativeIndices( myLowerStruct, fullStruct );
        info.childRelInds[1] = RelativeIndices( theirLowerStruct, fullStruct );
    }
    else
    {
        info.childSizes[0] = theirSize;
        info.childSizes[1] = info.child->size;
        info.childRelInds[0] = RelativeIndices( theirLowerStruct, fullStruct );
        info.childRelInds[1] = RelativeIndices( myLowerStruct, fullStruct );
    }

    // Form lower structure of this node by removing the node indices
    const Int lowerStructSize = fullStruct.size() - node.size;
    info.lowerStruct.resize( lowerStructSize );
    for( Int i=0; i<lowerStructSize; ++i )
        info.lowerStruct[i] = fullStruct[node.size+i];
    DEBUG_ONLY(
        // Ensure that the root process computed a lowerStruct of the same size
        Int rootLowerStructSize;
        if( mpi::Rank(node.comm) == 0 )
            rootLowerStructSize = lowerStructSize;
        mpi::Broadcast( &rootLowerStructSize, 1, 0, node.comm );
        if( rootLowerStructSize != lowerStructSize )
            RuntimeError("Root has different lower struct size");
    )
}

Int SymmetricAnalysis( const SymmNode& node, SymmNodeInfo& info, Int myOff )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricAnalysis"))

    // Recurse on the children
    // NOTE: Cleanup of existing info children should be added
    const Int numChildren = node.children.size();
    info.children.resize( numChildren );
    for( Int c=0; c<numChildren; ++c )
    {
        if( node.children[c] == nullptr )
            LogicError("Node child ",c," was nullptr");
        info.children[c] = new SymmNodeInfo(&info);
        myOff = 
            SymmetricAnalysis( *node.children[c], *info.children[c], myOff );
    }
    
    info.size = node.size;
    info.off= node.off;
    info.myOff= myOff; 
    info.origLowerStruct = node.lowerStruct;

    DEBUG_ONLY(
      if( !IsStrictlySorted(node.lowerStruct) )
      {
          if( IsSorted(node.lowerStruct) )
              LogicError("Repeat in original lower struct");
          else
              LogicError("Original lower struct not sorted");
      }
    )

    if( numChildren > 0 )
    {
        // Union the structures of the children with the original structure
        auto fullStruct = node.lowerStruct;
        for( SymmNodeInfo* child : info.children )
        {
            DEBUG_ONLY(
                if( !IsStrictlySorted(child->lowerStruct) )
                {
                    if( IsSorted(child->lowerStruct) )
                        LogicError("Repeat in child lower struct");
                    else
                        LogicError("Child lower struct not sorted");
                }
            )

            // Combine the structures of the children
            fullStruct = Union( fullStruct, child->lowerStruct );
        }

        // Now add in the node indices
        vector<Int> nodeInds( node.size );
        for( Int i=0; i<node.size; ++i )
            nodeInds[i] = node.off+ i;
        fullStruct = Union( fullStruct, nodeInds );

        // Construct the relative indices of the original lower structure
        info.origLowerRelInds = RelativeIndices( node.lowerStruct, fullStruct );

        // Construct the relative indices of the children
        info.childRelInds.resize( numChildren );
        for( Int c=0; c<numChildren; ++c )
            info.childRelInds[c] = 
                RelativeIndices( info.children[c]->lowerStruct, fullStruct );

        // Form lower struct of this node by removing node indices
        // (which take up the first node.size indices of fullStruct)
        const Int lowerStructSize = fullStruct.size()-node.size;
        info.lowerStruct.resize( lowerStructSize );
        for( Int i=0; i<lowerStructSize; ++i )
            info.lowerStruct[i] = fullStruct[node.size+i];
    }
    else
    {
        info.lowerStruct = node.lowerStruct;

        // Construct the trivial relative indices of the original structure
        const Int numOrigLowerInds = node.lowerStruct.size();
        info.origLowerRelInds.resize( numOrigLowerInds );
        for( Int i=0; i<numOrigLowerInds; ++i )
            info.origLowerRelInds[i] = i + info.size;
    }

    return myOff + info.size;
}

//
// This is the part of the analysis that requires fine-grain parallelism.
// For now, we will assume that the distributed part of the elimination 
// tree is binary.
//

void SymmetricAnalysis
( const DistSymmNode& node, DistSymmNodeInfo& info, bool computeFactRecvInds )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricAnalysis"))

    // Duplicate the communicator from the distributed eTree 
    mpi::Dup( node.comm, info.comm );
    info.grid = new Grid( info.comm );

    info.size = node.size;
    info.off = node.off;
    info.onLeft = node.onLeft;
    info.origLowerStruct = node.lowerStruct;
    if( node.duplicate != nullptr )
    {
        Int myOff = 0;
        info.duplicate = new SymmNodeInfo(&info);
        auto& dupInfo = *info.duplicate;
        SymmetricAnalysis( *node.duplicate, dupInfo, myOff );

        // The bottom node was analyzed locally, so just copy its results over
        info.myOff = dupInfo.myOff;
        info.lowerStruct = dupInfo.lowerStruct;
        info.origLowerRelInds = dupInfo.origLowerRelInds;
        info.childRelInds = dupInfo.childRelInds;

        const Int numChildren = dupInfo.childRelInds.size();
        info.childSizes.resize( numChildren );
        for( Int c=0; c<numChildren; ++c )
            info.childSizes[c] = dupInfo.children[c]->size;

        return;
    }

    if( node.child == nullptr )
        LogicError("Node child was nullptr");
    info.child = new DistSymmNodeInfo(&info);
    auto& childInfo = *info.child;
    SymmetricAnalysis( *node.child, childInfo, computeFactRecvInds );

    info.myOff = childInfo.myOff + childInfo.size;

    // Get the lower struct for the child we do not share
    Int theirSize;
    vector<Int> theirLowerStruct;
    GetLowerStruct( theirSize, theirLowerStruct, node, childInfo );

    // Perform one level of symbolic factorization and then compute
    // a wide variety of relative indices
    ComputeStructAndRelInds( theirSize, theirLowerStruct, node, info );
}

} // namespace El
