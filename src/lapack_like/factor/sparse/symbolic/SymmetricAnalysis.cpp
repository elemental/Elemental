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

inline void ComputeMultiVecCommMeta( DistSymmNodeInfo& node )
{
    DEBUG_ONLY(CallStackEntry cse("ComputeMultiVecCommMeta"))
    if( node.child == nullptr )
    {
        node.multiVecMeta.Empty();
        node.multiVecMeta.localOff = node.duplicate->myOff;
        node.multiVecMeta.localSize = node.duplicate->size;
        return;
    }
    ComputeMultiVecCommMeta( *node.child );

    // This is currently assumed (and will eventually be lifted)
    const Int numChildren = 2;
    vector<int> teamSizes(numChildren), teamOffs(numChildren);

    const int teamSize = mpi::Size( node.comm );
    const int teamRank = mpi::Rank( node.comm );
    const auto& childNode = *node.child;
    const int childTeamSize = mpi::Size( childNode.comm );
    const int childTeamRank = mpi::Rank( childNode.comm );
    const bool inFirstTeam = ( childTeamRank == teamRank );
    const bool leftIsFirst = ( childNode.onLeft==inFirstTeam );
    teamSizes[0] = 
        ( childNode.onLeft ? childTeamSize : teamSize-childTeamSize );
    teamSizes[1] = teamSize - teamSizes[0];
    teamOffs[0] = ( leftIsFirst ? 0 : teamSizes[1] );
    teamOffs[1] = ( leftIsFirst ? teamSizes[0] : 0 );

    const auto& myRelInds = 
        ( childNode.onLeft ? node.childRelInds[0] : node.childRelInds[1] );

    // Fill numChildSendInds
    auto& commMeta = node.multiVecMeta;
    commMeta.Empty();
    commMeta.numChildSendInds.resize( teamSize );
    El::MemZero( &commMeta.numChildSendInds[0], teamSize );
    const Int updateSize = childNode.lowerStruct.size();
    {
        const Int align = childNode.size % childTeamSize;
        const Int shift = Shift( childTeamRank, align, childTeamSize );
        const Int localHeight = Length( updateSize, shift, childTeamSize );
        for( Int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
        {
            const Int iChild = shift + iChildLoc*childTeamSize;
            const int destRank = myRelInds[iChild] % teamSize;
            ++commMeta.numChildSendInds[destRank];
        }
    }

    // Compute the solve recv indices
    commMeta.childRecvInds.resize( teamSize );
    for( Int c=0; c<numChildren; ++c )
    {
        const Int numInds = node.childRelInds[c].size();
        vector<Int> inds; 
        for( Int i=0; i<numInds; ++i )
            if( node.childRelInds[c][i] % teamSize == teamRank )
                inds.push_back( i );

        const Int numSolveInds = inds.size();
        for( Int iPre=0; iPre<numSolveInds; ++iPre )
        {
            const Int iChild = inds[iPre];
            const Int i = node.childRelInds[c][iChild];
            const Int iLoc = (i-teamRank) / teamSize;
            const int childRank = (node.childSizes[c]+iChild) % teamSizes[c];
            const int frontRank = teamOffs[c] + childRank;
            commMeta.childRecvInds[frontRank].push_back(iLoc);
        }
    }
    commMeta.localOff = node.child->multiVecMeta.localOff + 
                        node.child->multiVecMeta.localSize;
    commMeta.localSize = Length(node.size,teamRank,teamSize);
}

inline void ComputeFactorCommMeta
( DistSymmNodeInfo& node, bool computeFactRecvInds )
{
    DEBUG_ONLY(CallStackEntry cse("ComputeFactorCommMeta"))
    if( node.child == nullptr )
    {
        node.factorMeta.Empty();
        return;
    }
    auto& childNode = *node.child;
    ComputeFactorCommMeta( childNode, computeFactRecvInds );

    // Fill factorMeta.numChildSendInds 
    auto& commMeta = node.factorMeta;
    commMeta.Empty();
    const int gridHeight = node.grid->Height();
    const int gridWidth = node.grid->Width();
    const int childGridHeight = childNode.grid->Height();
    const int childGridWidth = childNode.grid->Width();
    const int childGridRow = childNode.grid->Row();
    const int childGridCol = childNode.grid->Col();
    const Int mySize = childNode.size;
    const Int updateSize = childNode.lowerStruct.size();
    const int teamSize = mpi::Size( node.comm );
    commMeta.numChildSendInds.resize( teamSize );
    El::MemZero( &commMeta.numChildSendInds[0], teamSize );
    const auto& myRelInds = 
        ( childNode.onLeft ? node.childRelInds[0] : node.childRelInds[1] );
    {
        const Int colAlign = mySize % childGridHeight;
        const Int rowAlign = mySize % childGridWidth;
        const Int colShift = Shift( childGridRow, colAlign, childGridHeight );
        const Int rowShift = Shift( childGridCol, rowAlign, childGridWidth );
        const Int localHeight = Length( updateSize, colShift, childGridHeight );
        const Int localWidth = Length( updateSize, rowShift, childGridWidth );
        for( Int jChildLoc=0; jChildLoc<localWidth; ++jChildLoc )
        {
            const Int jChild = rowShift + jChildLoc*childGridWidth;
            const int destGridCol = myRelInds[jChild] % gridWidth;

            Int localColShift;
            if( colShift > jChild )
                localColShift = 0;
            else if( (jChild-colShift) % childGridHeight == 0 )
                localColShift = (jChild-colShift)/childGridHeight;
            else
                localColShift = (jChild-colShift)/childGridHeight + 1;
            for( Int iChildLoc=localColShift; 
                     iChildLoc<localHeight; ++iChildLoc )
            {
                const Int iChild = colShift + iChildLoc*childGridHeight;
                if( iChild >= jChild )
                {
                    const int destGridRow = myRelInds[iChild] % gridHeight;
                    const int destRank = destGridRow+destGridCol*gridHeight;
                    ++commMeta.numChildSendInds[destRank];
                }
            }
        }
    }

    // Optionally compute the recv indices for the factorization. 
    // This is optional since it requires a nontrivial amount of storage.
    if( computeFactRecvInds )
        ComputeFactRecvInds( node );
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

    ComputeFactorCommMeta( info, computeFactRecvInds );
    
    // This is thankfully independent of the number of right-hand sides,   
    // unlike the 2d equivalent
    ComputeMultiVecCommMeta( info );
}

} // namespace El
