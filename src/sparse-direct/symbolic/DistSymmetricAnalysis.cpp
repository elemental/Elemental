/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

inline void PairwiseExchangeLowerStruct
( Int& theirSize, std::vector<Int>& theirLowerStruct,
  const DistSymmNode& node, const DistSymmNodeInfo& childNodeInfo )
{
    // Determine our partner's rank for this exchange in node's communicator
    const int teamRank = mpi::Rank( node.comm );
    const int teamSize = mpi::Size( node.comm );
    const int childTeamRank = mpi::Rank( childNodeInfo.comm );
    const int myTeamSize = mpi::Size( childNodeInfo.comm );
    const int otherTeamSize = teamSize - myTeamSize;
    const bool inFirstTeam = ( teamRank == childTeamRank );
    const int partner =
        ( inFirstTeam ? teamRank+myTeamSize : teamRank-otherTeamSize );

    // SendRecv the message lengths
    const Int mySize = childNodeInfo.size;
    const Int myLowerStructSize = childNodeInfo.lowerStruct.size();
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
    ( &childNodeInfo.lowerStruct[0], myLowerStructSize, partner,
      &theirLowerStruct[0], theirLowerStructSize, partner, node.comm );
}

inline void BroadcastLowerStruct
( Int& theirSize, std::vector<Int>& theirLowerStruct,
  const DistSymmNode& node, const DistSymmNodeInfo& childNodeInfo )
{
    // Determine our partner's rank for this exchange in node's communicator
    const int teamRank = mpi::Rank( node.comm );
    const int teamSize = mpi::Size( node.comm );
    const int childTeamRank = mpi::Rank( childNodeInfo.comm );
    const int myTeamSize = mpi::Size( childNodeInfo.comm );
    const int otherTeamSize = teamSize - myTeamSize;
    const bool inFirstTeam = ( teamRank == childTeamRank );

    if( childTeamRank == 0 )
    {
        const int partner =
            ( inFirstTeam ? teamRank+myTeamSize : teamRank-otherTeamSize );

        // SendRecv the message lengths
        const Int mySize = childNodeInfo.size;
        const Int myLowerStructSize = childNodeInfo.lowerStruct.size();
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
        ( &childNodeInfo.lowerStruct[0], myLowerStructSize, partner,
          &theirLowerStruct[0], theirLowerStructSize, partner, node.comm );

        // Broadcast the other team's child's sizes
        mpi::Broadcast( initialRecvs, 2, 0, childNodeInfo.comm );

        // Broadcast the other team's child's lower struct
        mpi::Broadcast
        ( &theirLowerStruct[0], theirLowerStructSize, 0, childNodeInfo.comm );
    } 
    else
    {
        // Receive the other team's child's sizes
        Int initialRecvs[2];
        mpi::Broadcast( initialRecvs, 2, 0, childNodeInfo.comm );
        theirSize = initialRecvs[0];
        const Int theirLowerStructSize = initialRecvs[1];

        // Receive the other team's child's lower struct
        theirLowerStruct.resize( theirLowerStructSize );
        mpi::Broadcast
        ( &theirLowerStruct[0], theirLowerStructSize, 0, childNodeInfo.comm );
    }
}

inline void GetLowerStruct
( Int& theirSize, std::vector<Int>& theirLowerStruct,
  const DistSymmNode& node, const DistSymmNode& childNode, 
  const DistSymmNodeInfo& childNodeInfo )
{
    const int teamSize = mpi::Size( node.comm );
    const int childTeamSize = mpi::Size( childNode.comm );
    const int leftTeamSize =
        ( childNode.onLeft ? childTeamSize : teamSize-childTeamSize );
    const int rightTeamSize = teamSize - leftTeamSize;
    if( leftTeamSize == rightTeamSize )
        PairwiseExchangeLowerStruct
        ( theirSize, theirLowerStruct, node, childNodeInfo );
    else
        BroadcastLowerStruct
        ( theirSize, theirLowerStruct, node, childNodeInfo );
}

inline void ComputeStructAndRelInds
( Int theirSize, const std::vector<Int>& theirLowerStruct,
  const DistSymmNode& node,         const DistSymmNode& childNode, 
        DistSymmNodeInfo& nodeInfo, const DistSymmNodeInfo& childNodeInfo )
{
    const std::vector<Int>& myLowerStruct = childNodeInfo.lowerStruct;
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
    std::vector<Int> nodeInds( node.size );
    for( int i=0; i<node.size; ++i )
        nodeInds[i] = node.off + i;
    auto fullStruct = Union( nodeInds, partialStruct );

    // Construct the relative indices of the original lower structure
    nodeInfo.origLowerRelInds = RelativeIndices( node.lowerStruct, fullStruct );

    // Construct the relative indices of the children
    if( childNode.onLeft )
    {
        nodeInfo.leftSize = childNodeInfo.size;
        nodeInfo.rightSize = theirSize;
        nodeInfo.leftRelInds = RelativeIndices( myLowerStruct, fullStruct );
        nodeInfo.rightRelInds = RelativeIndices( theirLowerStruct, fullStruct );
    }
    else
    {
        nodeInfo.leftSize = theirSize;
        nodeInfo.rightSize = childNodeInfo.size;
        nodeInfo.leftRelInds = RelativeIndices( theirLowerStruct, fullStruct );
        nodeInfo.rightRelInds = RelativeIndices( myLowerStruct, fullStruct );
    }

    // Form lower structure of this node by removing the node indices
    const int lowerStructSize = fullStruct.size() - node.size;
    nodeInfo.lowerStruct.resize( lowerStructSize );
    for( int i=0; i<lowerStructSize; ++i )
        nodeInfo.lowerStruct[i] = fullStruct[node.size+i];
    DEBUG_ONLY(
        // Ensure that the root process computed a lowerStruct of the same size
        int rootLowerStructSize;
        if( mpi::Rank(node.comm) == 0 )
            rootLowerStructSize = lowerStructSize;
        mpi::Broadcast( &rootLowerStructSize, 1, 0, node.comm );
        if( rootLowerStructSize != lowerStructSize )
            RuntimeError("Root has different lower struct size");
    )
}

inline void ComputeMultiVecCommMeta( DistSymmInfo& info )
{
    DEBUG_ONLY(CallStackEntry cse("ComputeMultiVecCommMeta"))
    // Handle the interface node
    Int localOff = info.localNodes.back().myOff;
    info.distNodes[0].multiVecMeta.Empty();
    info.distNodes[0].multiVecMeta.localOff = localOff;
    info.distNodes[0].multiVecMeta.localSize = info.localNodes.back().size;
    localOff += info.distNodes[0].multiVecMeta.localSize;

    // Handle the truly distributed nodes

    const int numDist = info.distNodes.size();
    for( int s=1; s<numDist; ++s )
    {
        DistSymmNodeInfo& node = info.distNodes[s];
        const int teamSize = mpi::Size( node.comm );
        const int teamRank = mpi::Rank( node.comm );

        const DistSymmNodeInfo& childNode = info.distNodes[s-1];
        const int childTeamSize = mpi::Size( childNode.comm );
        const int childTeamRank = mpi::Rank( childNode.comm );
        const bool inFirstTeam = ( childTeamRank == teamRank );
        const bool leftIsFirst = ( childNode.onLeft==inFirstTeam );
        const int leftTeamSize =
            ( childNode.onLeft ? childTeamSize : teamSize-childTeamSize );
        const int rightTeamSize = teamSize - leftTeamSize;
        const int leftTeamOff = ( leftIsFirst ? 0 : rightTeamSize );
        const int rightTeamOff = ( leftIsFirst ? leftTeamSize : 0 );

        const std::vector<Int>& myRelInds = 
            ( childNode.onLeft ? node.leftRelInds : node.rightRelInds );

        // Fill numChildSendInds
        MultiVecCommMeta& commMeta = node.multiVecMeta;
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

        const Int numLeftInds = node.leftRelInds.size();
        const Int numRightInds = node.rightRelInds.size();
        std::vector<Int> leftInds, rightInds; 
        for( Int i=0; i<numLeftInds; ++i )
            if( node.leftRelInds[i] % teamSize == teamRank )
                leftInds.push_back( i );
        for( Int i=0; i<numRightInds; ++i )
            if( node.rightRelInds[i] % teamSize == teamRank )
                rightInds.push_back( i );

        //
        // Compute the solve recv indices
        //
        commMeta.childRecvInds.resize( teamSize );

        // Compute the recv indices for the left child 
        const Int numLeftSolveInds = leftInds.size();
        for( Int iPre=0; iPre<numLeftSolveInds; ++iPre )
        {
            const Int iChild = leftInds[iPre];
            const Int iFront = node.leftRelInds[iChild];
            const Int iFrontLoc = (iFront-teamRank) / teamSize;
            const int childRank = (node.leftSize+iChild) % leftTeamSize;
            const int frontRank = leftTeamOff + childRank;
            commMeta.childRecvInds[frontRank].push_back(iFrontLoc);
        }

        // Compute the recv indices for the right child
        const Int numRightSolveInds = rightInds.size();
        for( Int iPre=0; iPre<numRightSolveInds; ++iPre )
        {
            const Int iChild = rightInds[iPre];
            const Int iFront = node.rightRelInds[iChild];
            const Int iFrontLoc = (iFront-teamRank) / teamSize;
            const int childRank = (node.rightSize+iChild) % rightTeamSize;
            const int frontRank = rightTeamOff + childRank;
            commMeta.childRecvInds[frontRank].push_back(iFrontLoc);
        }

        commMeta.localOff = localOff;
        commMeta.localSize = Length(node.size,teamRank,teamSize);
        localOff += commMeta.localSize;
    }
}

inline void ComputeFactorCommMeta
( DistSymmInfo& info, bool computeFactRecvInds )
{
    DEBUG_ONLY(CallStackEntry cse("ComputeFactorCommMeta"))
    info.distNodes[0].factorMeta.Empty();
    const Int numDist = info.distNodes.size();
    for( Int s=1; s<numDist; ++s )
    {
        DistSymmNodeInfo& node = info.distNodes[s];
        const int teamSize = mpi::Size( node.comm );
        const DistSymmNodeInfo& childNode = info.distNodes[s-1];

        // Fill factorMeta.numChildSendInds 
        FactorCommMeta& commMeta = node.factorMeta;
        commMeta.Empty();
        const int gridHeight = node.grid->Height();
        const int gridWidth = node.grid->Width();
        const int childGridHeight = childNode.grid->Height();
        const int childGridWidth = childNode.grid->Width();
        const int childGridRow = childNode.grid->Row();
        const int childGridCol = childNode.grid->Col();
        const Int mySize = childNode.size;
        const Int updateSize = childNode.lowerStruct.size();
        commMeta.numChildSendInds.resize( teamSize );
        El::MemZero( &commMeta.numChildSendInds[0], teamSize );
        const std::vector<Int>& myRelInds = 
            ( childNode.onLeft ? node.leftRelInds : node.rightRelInds );
        {
            const Int colAlign = mySize % childGridHeight;
            const Int rowAlign = mySize % childGridWidth;
            const Int colShift = 
                Shift( childGridRow, colAlign, childGridHeight );
            const Int rowShift = 
                Shift( childGridCol, rowAlign, childGridWidth );
            const Int localHeight = 
                Length( updateSize, colShift, childGridHeight );
            const Int localWidth = 
                Length( updateSize, rowShift, childGridWidth );
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
            ComputeFactRecvInds( node, childNode );
    }
}

//
// This is the part of the analysis that requires fine-grain parallelism.
// For now, we will assume that the distributed part of the elimination 
// tree is binary.
//
void DistSymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, bool computeFactRecvInds )
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmetricAnalysis"))
    const Unsigned numDist = eTree.distNodes.size();
    info.distNodes.resize( numDist );

    // The bottom node was analyzed locally, so just copy its results over
    const SymmNodeInfo& topLocal = info.localNodes.back();
    DistSymmNodeInfo& bottomDist = info.distNodes[0];
    bottomDist.onLeft = eTree.distNodes[0].onLeft;
    mpi::Dup( eTree.distNodes[0].comm, bottomDist.comm );
    bottomDist.grid = new Grid( bottomDist.comm );
    bottomDist.size = topLocal.size;
    bottomDist.off = topLocal.off;
    bottomDist.myOff = topLocal.myOff;
    bottomDist.lowerStruct = topLocal.lowerStruct;
    bottomDist.origLowerStruct = topLocal.origLowerStruct;
    bottomDist.origLowerRelInds = topLocal.origLowerRelInds;
    bottomDist.leftRelInds = topLocal.leftRelInds;
    bottomDist.rightRelInds = topLocal.rightRelInds;
    bottomDist.leftSize = -1; // not needed, could compute though
    bottomDist.rightSize = -1; // not needed, could compute though

    // Perform the distributed part of the symbolic factorization
    Int myOff = bottomDist.myOff + bottomDist.size;
    for( Unsigned s=1; s<numDist; ++s )
    {
        const DistSymmNode& node = eTree.distNodes[s];
        const DistSymmNode& childNode = eTree.distNodes[s-1];
        const DistSymmNodeInfo& childNodeInfo = info.distNodes[s-1];
        DistSymmNodeInfo& nodeInfo = info.distNodes[s];
        nodeInfo.onLeft = node.onLeft;
        nodeInfo.size = node.size;
        nodeInfo.off = node.off;
        nodeInfo.myOff = myOff;
        nodeInfo.origLowerStruct = node.lowerStruct;

        // Duplicate the communicator from the distributed eTree 
        mpi::Dup( node.comm, nodeInfo.comm );
        nodeInfo.grid = new Grid( nodeInfo.comm );

        // Get the lower struct for the child we do not share
        Int theirSize;
        std::vector<Int> theirLowerStruct;
        GetLowerStruct
        ( theirSize, theirLowerStruct, node, childNode, childNodeInfo );

        // Perform one level of symbolic factorization and then compute
        // a wide variety of relative indices
        ComputeStructAndRelInds
        ( theirSize, theirLowerStruct, node, childNode, 
          nodeInfo, childNodeInfo );

        myOff += nodeInfo.size;
    }

    ComputeFactorCommMeta( info, computeFactRecvInds );
    
    // This is thankfully independent of the number of right-hand sides,   
    // unlike the 2d equivalent
    ComputeMultiVecCommMeta( info );
}

// TODO: Simplify this implementation
void ComputeFactRecvInds
( const DistSymmNodeInfo& node, const DistSymmNodeInfo& childNode )
{
    DEBUG_ONLY(CallStackEntry cse("ComputeFactRecvInds"))
    // Communicate to get the grid sizes
    int childGridDims[4];
    GetChildGridDims( node, childNode, childGridDims );
    const int leftGridHeight = childGridDims[0];
    const int leftGridWidth = childGridDims[1];
    const int rightGridHeight = childGridDims[2];
    const int rightGridWidth = childGridDims[3];

    const int teamSize = mpi::Size( node.comm );
    const int teamRank = mpi::Rank( node.comm );
    const bool onLeft = childNode.onLeft;
    const int childTeamSize = mpi::Size( childNode.comm );
    const int leftTeamSize =
        ( onLeft ? childTeamSize : teamSize-childTeamSize );
    const int rightTeamSize = teamSize - leftTeamSize;
    DEBUG_ONLY(
        if( leftTeamSize != leftGridHeight*leftGridWidth )
            RuntimeError("Computed left grid incorrectly");
        if( rightTeamSize != rightGridHeight*rightGridWidth )
            RuntimeError("Computed right grid incorrectly");
    )

    const FactorCommMeta& commMeta = node.factorMeta;
    const int gridHeight = node.grid->Height();
    const int gridWidth = node.grid->Width();
    const int gridRow = node.grid->Row();
    const int gridCol = node.grid->Col();
    const Int numLeftInds = node.leftRelInds.size();
    const Int numRightInds = node.rightRelInds.size();
    std::vector<Int> leftRowInds, leftColInds, rightRowInds, rightColInds;
    for( Int i=0; i<numLeftInds; ++i )
        if( node.leftRelInds[i] % gridHeight == gridRow )
            leftColInds.push_back( i );
    for( Int i=0; i<numLeftInds; ++i )
        if( node.leftRelInds[i] % gridWidth == gridCol )
            leftRowInds.push_back( i );
    for( Int i=0; i<numRightInds; ++i )
        if( node.rightRelInds[i] % gridHeight == gridRow )
            rightColInds.push_back( i );
    for( Int i=0; i<numRightInds; ++i )
        if( node.rightRelInds[i] % gridWidth == gridCol )
            rightRowInds.push_back( i );

    // Compute the recv indices of the left child from each process 
    const int childTeamRank = mpi::Rank( childNode.comm );
    const bool inFirstTeam = ( childTeamRank == teamRank );
    const bool leftIsFirst = ( onLeft==inFirstTeam );
    const int leftTeamOff = ( leftIsFirst ? 0 : rightTeamSize );
    commMeta.childRecvInds.resize( teamSize );
    std::vector<Int>::const_iterator it;
    const Int numLeftColInds = leftColInds.size();
    const Int numLeftRowInds = leftRowInds.size();
    const Int numRightColInds = rightColInds.size();
    const Int numRightRowInds = rightRowInds.size();
    for( Int jPre=0; jPre<numLeftRowInds; ++jPre )
    {
        const Int jChild = leftRowInds[jPre];
        const Int jFront = node.leftRelInds[jChild];
        DEBUG_ONLY(
            if( (jFront-gridCol) % gridWidth != 0 )
                LogicError("Invalid left jFront");
        )
        const Int jFrontLoc = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.leftSize) % leftGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound( leftColInds.begin(), leftColInds.end(), jChild );
        const Int iPreStart = Int(it-leftColInds.begin());
        for( Int iPre=iPreStart; iPre<numLeftColInds; ++iPre )
        {
            const Int iChild = leftColInds[iPre];
            const Int iFront = node.leftRelInds[iChild];
            DEBUG_ONLY(
                if( iChild < jChild )
                    LogicError("Invalid left iChild");
                if( (iFront-gridRow) % gridHeight != 0 )
                    LogicError("Invalid left iFront");
            )
            const Int iFrontLoc = (iFront-gridRow) / gridHeight;

            const int childRow = (iChild+node.leftSize) % leftGridHeight;
            const int childRank = childRow + childCol*leftGridHeight;

            const int frontRank = leftTeamOff + childRank;
            commMeta.childRecvInds[frontRank].push_back(iFrontLoc);
            commMeta.childRecvInds[frontRank].push_back(jFrontLoc);
        }
    }
    
    // Compute the recv indices of the right child from each process 
    const Int rightTeamOff = ( leftIsFirst ? leftTeamSize : 0 );
    for( Int jPre=0; jPre<numRightRowInds; ++jPre )
    {
        const Int jChild = rightRowInds[jPre];
        const Int jFront = node.rightRelInds[jChild];
        DEBUG_ONLY(
            if( (jFront-gridCol) % gridWidth != 0 )
                LogicError("Invalid right jFront");
        )
        const Int jFrontLoc = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.rightSize) % rightGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound
             ( rightColInds.begin(), rightColInds.end(), jChild );
        const Int iPreStart = Int(it-rightColInds.begin());
        for( Int iPre=iPreStart; iPre<numRightColInds; ++iPre )
        {
            const Int iChild = rightColInds[iPre];
            const Int iFront = node.rightRelInds[iChild];
            DEBUG_ONLY(
                if( iChild < jChild )
                    LogicError("Invalid right iChild");
                if( (iFront-gridRow) % gridHeight != 0 )
                    LogicError("Invalid right iFront");
            )
            const Int iFrontLoc = (iFront-gridRow) / gridHeight;

            const int childRow = (iChild+node.rightSize) % rightGridHeight;
            const int childRank = childRow + childCol*rightGridHeight;

            const int frontRank = rightTeamOff + childRank;
            commMeta.childRecvInds[frontRank].push_back(iFrontLoc);
            commMeta.childRecvInds[frontRank].push_back(jFrontLoc);
        }
    }
}

void GetChildGridDims
( const DistSymmNodeInfo& node, const DistSymmNodeInfo& childNode, 
  int* childGridDims )
{
    const bool onLeft = childNode.onLeft;
    const int childTeamRank = mpi::Rank( childNode.comm );
    El::MemZero( childGridDims, 4 );
    if( onLeft && childTeamRank == 0 )
    {
        childGridDims[0] = childNode.grid->Height();
        childGridDims[1] = childNode.grid->Width();
    }
    else if( !onLeft && childTeamRank == 0 )
    {
        childGridDims[2] = childNode.grid->Height();
        childGridDims[3] = childNode.grid->Width();
    }
    mpi::AllReduce( childGridDims, 4, mpi::SUM, node.comm );
}

} // namespace El
