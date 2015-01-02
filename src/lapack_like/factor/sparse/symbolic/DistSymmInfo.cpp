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

DistSymmInfo::~DistSymmInfo()
{
    const Int numDist = distNodes.size();
    for( Int s=0; s<numDist; ++s )
    {
        delete distNodes[s].grid;
        mpi::Free( distNodes[s].comm );
    }
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
